// @File(label = "Analysis results folder", style = "directory") resultsDir
// @Integer(label = "Targets per site", value = 5, min = 1) targetsPerSite
// @Double(label = "Minimum distance between objects in the overview image (pixels, 0 = no constraint)", value = 340, min = 0) minDistPx

// ------------------------------------------------------------------------------
// Fiji: MD HCS target curation (v1.8.0)
// Created: 2026-04-28 | Updated: 2026-04-29
// Author: thom.dehoog@zmb.uzh.ch | ZMB Center for Microscopy and Image Analysis, UZH
//
// If you publish a paper using this macro, please acknowledge.
//
// Description: Curate MetaXpress TargetData CSVs produced by IN Carta segmentation
//   so a downstream targeted acquisition images N cells per site instead of all
//   detected cells. Picks are random per site with a Euclidean minimum-distance
//   constraint between bounding-box centres so re-imaging tiles stay separated.
//   Distances operate directly in the overview pixel units already present in
//   the TargetData CSVs; no calibration file is required.
//
// Behaviour:
//   - First run renames TargetData/ -> TargetData_original/ (idempotent).
//   - Every run reads from TargetData_original/ and rewrites TargetData/ from
//     scratch. The original IN Carta output is never modified.
//   - Every run also writes TargetData_curated/, a mirror of TargetData/ plus
//     curation_changes.csv for auditing what changed.
//   - If TargetData_original/ already exists, TargetData/ must match the
//     previous TargetData_curated/ mirror before it is overwritten. This avoids
//     accidentally deleting newly regenerated IN Carta TargetData.
//   - SummaryInfo, ObjectData, FieldData, WellData are left untouched.
//   - At the end, a summary dialog reports settings used and any sites that
//     ended up under target (with reason: low_cells or constrained).
//
// Expected input:
//   <resultsDir>/
//     TargetData/                                 (IN Carta output; renamed on first run)
//       *.csv                                     (all CSVs in this folder are considered)
//
// Output:
//   <resultsDir>/
//     TargetData/                                 (curated; MetaXpress reads from here)
//     TargetData_curated/                         (curated mirror plus changes log)
//     TargetData_original/                        (preserved IN Carta output)
//
// Picking algorithm:
//   - Random shuffle of detected cells (Fisher-Yates), walk through, accept
//     each candidate whose centre is at least minDistPx (Euclidean) from every
//     already-picked centre. Stop at N or when no more candidates fit.
//   - minDistPx = 0 disables the constraint (operator gets the first N from the
//     random shuffle). Useful as a baseline run.
//   - To choose minDistPx: the bounding-box columns in the CSVs are in overview
//     pixels, so set minDistPx to the diameter (in overview pixels) below which
//     two re-imaging tiles would overlap, plus a comfort margin. For a typical
//     ZMB MD HCS configuration (10x overview at ~0.665 um/px, target acquisition
//     at 60x + 1.5x lens changer with a ~148 um tile side), the tile diagonal is
//     ~210 um = ~316 px. The default of 340 adds a small comfort gap on top.
//
// Per-site outcome:
//   - full        n_picked == N
//   - low_cells   n_detected < N (ran out of candidates)
//   - constrained n_detected >= N but min-distance forbade more picks
//
// Notes:
//   - Pick order is random but the seed is fixed in source so runs are
//     reproducible: same input + same settings always produces the same
//     picks. Edit the `seed` constant near the top of the script to vary.
//     The seed is global across files in a run, so picks for later files
//     depend on earlier ones; if any earlier CSV changes, picks downstream
//     of it will shift even if those CSVs were untouched.
//   - The picks are NOT a uniform random sample of valid N-cell subsets:
//     the shuffle-then-greedy procedure with a min-distance filter biases
//     against cells that cluster near already-picked ones. Fine for picking
//     sites to image, not appropriate as input for spatial-statistics.
//   - The macro requires the object_id and BoundingBox columns in each
//     TargetData CSV; sites without them are skipped with a warning.
// ------------------------------------------------------------------------------

// --- Helpers ---

// Find a column index by exact header name; returns -1 if absent.
function findColumn(headerCols, name) {
    for (i = 0; i < headerCols.length; i++) {
        c = headerCols[i];
        if (endsWith(c, "\r")) c = substring(c, 0, lengthOf(c)-1);
        c = replace(c, "^[ \t]+|[ \t]+$", "");
        if (c == name) return i;
    }
    return -1;
}

// Strip a trailing \r so split-by-\n on CRLF text yields clean rows.
function stripCR(s) {
    if (endsWith(s, "\r")) return substring(s, 0, lengthOf(s)-1);
    return s;
}

// Strip a UTF-8 / UTF-16 BOM if present at the start of a string. Some
// exporters add a BOM to the first byte of the file, which would otherwise
// glue itself to the first column header and break exact-match column lookup.
function stripBom(s) {
    if (lengthOf(s) >= 1 && substring(s, 0, 1) == fromCharCode(65279))
        return substring(s, 1, lengthOf(s));
    if (lengthOf(s) >= 3 &&
        substring(s, 0, 1) == fromCharCode(239) &&
        substring(s, 1, 2) == fromCharCode(187) &&
        substring(s, 2, 3) == fromCharCode(191))
        return substring(s, 3, lengthOf(s));
    return s;
}

// Two-digit zero-padded decimal for date formatting.
function pad2(n) {
    if (n < 10) return "0" + d2s(n, 0);
    return d2s(n, 0);
}

// Split one CSV row, preserving commas inside quoted fields.
function splitCsvRow(row) {
    fields = newArray(0);
    field = "";
    quote = fromCharCode(34);
    inQuotes = false;

    for (ii = 0; ii < lengthOf(row); ii++) {
        ch = substring(row, ii, ii+1);
        if (ch == quote) {
            if (inQuotes && ii+1 < lengthOf(row)) {
                if (substring(row, ii+1, ii+2) == quote) {
                    field = field + quote;
                    ii++;
                } else {
                    inQuotes = !inQuotes;
                }
            } else {
                inQuotes = !inQuotes;
            }
        } else if (ch == "," && !inQuotes) {
            fields = Array.concat(fields, field);
            field = "";
        } else {
            field = field + ch;
        }
    }

    fields = Array.concat(fields, field);
    return fields;
}

// Classify a site outcome based on detected and picked counts.
// Returns "full", "low_cells", or "constrained".
function siteStatus(nDetected, nPicked, target) {
    if (nPicked == target) return "full";
    if (nDetected < target) return "low_cells";
    return "constrained";
}

// Refuse to delete TargetData/ on rerun unless it is the previous curated output.
function assertCurrentTargetDataIsPreviousCuration(curDir, auditDir) {
    if (!File.isDirectory(curDir)) return;
    if (!File.isDirectory(auditDir))
        exit("TargetData_original/ already exists, but TargetData_curated/ is missing. Refusing to overwrite TargetData/ because it may contain newly regenerated IN Carta output.");
    if (!File.exists(auditDir + "curation_changes.csv"))
        exit("TargetData_original/ already exists, but TargetData_curated/curation_changes.csv is missing. Refusing to overwrite TargetData/ because its provenance cannot be verified.");

    curList = getFileList(curDir);
    auditList = getFileList(auditDir);
    for (ii = 0; ii < curList.length; ii++) {
        name = curList[ii];
        if (File.isDirectory(curDir + name))
            exit("Unexpected subfolder in TargetData/: " + curDir + name);
        if (!File.exists(auditDir + name))
            exit("Current TargetData/ does not match the previous TargetData_curated/ mirror. Refusing to overwrite possible new IN Carta output: " + curDir + name);
        if (File.length(curDir + name) != File.length(auditDir + name))
            exit("Current TargetData/ differs from the previous TargetData_curated/ mirror. Refusing to overwrite possible new IN Carta output: " + curDir + name);
        if (File.openAsString(curDir + name) != File.openAsString(auditDir + name))
            exit("Current TargetData/ differs from the previous TargetData_curated/ mirror. Refusing to overwrite possible new IN Carta output: " + curDir + name);
    }
    for (ii = 0; ii < auditList.length; ii++) {
        name = auditList[ii];
        if (name == "curation_changes.csv") continue;
        if (File.isDirectory(auditDir + name))
            exit("Unexpected subfolder in TargetData_curated/: " + auditDir + name);
        if (!File.exists(curDir + name))
            exit("Current TargetData/ is missing a file from the previous TargetData_curated/ mirror. Refusing to overwrite because the current state is inconsistent: " + name);
    }
}

// --- Resolve paths and set up output folders ---
if (!File.isDirectory(resultsDir))
    exit("Not a directory: " + resultsDir);
if (!endsWith(resultsDir, File.separator)) resultsDir += File.separator;

origDir = resultsDir + "TargetData_original" + File.separator;
curDir  = resultsDir + "TargetData" + File.separator;
auditDir = resultsDir + "TargetData_curated" + File.separator;

// First run: rename IN Carta's TargetData/ -> TargetData_original/.
if (!File.isDirectory(origDir)) {
    if (!File.isDirectory(curDir))
        exit("No TargetData/ folder at: " + curDir);
    if (!File.rename(curDir, origDir))
        exit("Could not rename TargetData/ -> TargetData_original/.");
    print("First run: renamed TargetData/ -> TargetData_original/");
} else {
    assertCurrentTargetDataIsPreviousCuration(curDir, auditDir);
}

// Refresh curated output folders (always overwritten from TargetData_original/).
if (File.isDirectory(curDir)) {
    stale = getFileList(curDir);
    for (i = 0; i < stale.length; i++) {
        if (!File.delete(curDir + stale[i]))
            exit("Could not delete stale file from TargetData/: " + curDir + stale[i]);
    }
} else {
    File.makeDirectory(curDir);
    if (!File.isDirectory(curDir))
        exit("Could not create TargetData/ folder: " + curDir);
}
if (File.isDirectory(auditDir)) {
    stale = getFileList(auditDir);
    for (i = 0; i < stale.length; i++) {
        if (!File.delete(auditDir + stale[i]))
            exit("Could not delete stale file from TargetData_curated/: " + auditDir + stale[i]);
    }
} else {
    File.makeDirectory(auditDir);
    if (!File.isDirectory(auditDir))
        exit("Could not create TargetData_curated/ folder: " + auditDir);
}

// Fixed seed for reproducible curation. Edit to draw a different sample.
seed = 42;
random("seed", seed);

// --- Curate each site ---
list = getFileList(origDir);
Array.sort(list);
totalIn = 0;
totalOut = 0;
nFiles = 0;
totalRejectedRows = 0;
underList = newArray(0);  // human-readable lines for sites that ended up under-target
skippedList = newArray(0); // human-readable lines for files that could not be curated
rejectedList = newArray(0); // human-readable lines for files with rejected malformed rows
changeRows = newArray(0);  // tab-separated rows for TargetData_curated/curation_changes.csv

for (f = 0; f < list.length; f++) {
    name = list[f];
    if (!endsWith(name, ".csv")) continue;

    raw = File.openAsString(origDir + name);
    lines = split(raw, "\n");
    if (lines.length < 1 || lengthOf(stripCR(lines[0])) == 0) {
        print("skip (empty or missing header): " + name);
        skippedList = Array.concat(skippedList, name + ": empty or missing header");
        changeRows = Array.concat(changeRows, name + "\tNA\tNA\t0\tskipped_empty_or_missing_header\t");
        continue;
    }

    header = stripBom(stripCR(lines[0]));
    cols = splitCsvRow(header);
    iBX  = findColumn(cols, "T1$AS_FID_Blob_BoundingBoxX");
    iBY  = findColumn(cols, "T1$AS_FID_Blob_BoundingBoxY");
    iBW  = findColumn(cols, "T1$AS_FID_Blob_BoundingBoxWidth");
    iBH  = findColumn(cols, "T1$AS_FID_Blob_BoundingBoxHeight");
    iOID = findColumn(cols, "object_id");
    if (iBX < 0 || iBY < 0 || iBW < 0 || iBH < 0) {
        print("skip (missing bounding-box columns): " + name);
        skippedList = Array.concat(skippedList, name + ": missing bounding-box columns");
        changeRows = Array.concat(changeRows, name + "\tNA\tNA\t0\tskipped_missing_bounding_box_columns\t");
        continue;
    }
    if (iOID < 0) {
        print("skip (missing object_id column): " + name);
        skippedList = Array.concat(skippedList, name + ": missing object_id column");
        changeRows = Array.concat(changeRows, name + "\tNA\tNA\t0\tskipped_missing_object_id_column\t");
        continue;
    }

    rows = newArray(0);
    cx   = newArray(0);
    cy   = newArray(0);
    oid  = newArray(0);
    rejectedRows = 0;
    for (i = 1; i < lines.length; i++) {
        row = stripCR(lines[i]);
        if (lengthOf(row) == 0) continue;
        flds = splitCsvRow(row);
        if (flds.length < cols.length) {
            rejectedRows++;
            continue;
        }
        x = parseFloat(flds[iBX]);
        y = parseFloat(flds[iBY]);
        w = parseFloat(flds[iBW]);
        h = parseFloat(flds[iBH]);
        if (isNaN(x) || isNaN(y) || isNaN(w) || isNaN(h)) {
            rejectedRows++;
            continue;
        }
        rows = Array.concat(rows, row);
        cx   = Array.concat(cx, x + w/2.0);
        cy   = Array.concat(cy, y + h/2.0);
        oid  = Array.concat(oid, flds[iOID]);
    }

    n = rows.length;
    totalIn += n;
    totalRejectedRows += rejectedRows;
    nFiles++;

    // Fisher-Yates shuffle of [0..n-1].
    idx = newArray(n);
    for (i = 0; i < n; i++) idx[i] = i;
    for (i = n-1; i > 0; i--) {
        j = floor(random() * (i+1));
        tmp = idx[i]; idx[i] = idx[j]; idx[j] = tmp;
    }

    // Greedy pick. If minDistPx == 0, accept the first N from the shuffle.
    picked = newArray(0);
    for (k = 0; k < n && picked.length < targetsPerSite; k++) {
        cand = idx[k];
        ok = true;
        if (minDistPx > 0) {
            for (m = 0; m < picked.length; m++) {
                p = picked[m];
                dx = cx[cand] - cx[p];
                dy = cy[cand] - cy[p];
                if (sqrt(dx*dx + dy*dy) < minDistPx) {
                    ok = false;
                    m = picked.length; // break inner loop
                }
            }
        }
        if (ok) picked = Array.concat(picked, cand);
    }

    status = siteStatus(n, picked.length, targetsPerSite);

    // Write the curated CSV (original header, CRLF line endings, picked rows only).
    out = header + "\r\n";
    pickedIdsStr = "";
    for (k = 0; k < picked.length; k++) {
        out = out + rows[picked[k]] + "\r\n";
        if (k > 0) pickedIdsStr = pickedIdsStr + ";";
        pickedIdsStr = pickedIdsStr + oid[picked[k]];
    }
    File.saveString(out, curDir + name);
    File.saveString(out, auditDir + name);

    siteLabel = substring(name, 0, lengthOf(name) - 4);
    if (status != "full")
        underList = Array.concat(underList,
            siteLabel + ":  " + n + " detected, " + picked.length + " picked  (" + status + ")");
    if (rejectedRows > 0)
        rejectedList = Array.concat(rejectedList,
            siteLabel + ":  " + rejectedRows + " rejected malformed rows");
    changeRows = Array.concat(changeRows,
        name + "\t" + n + "\t" + rejectedRows + "\t" + picked.length + "\t" + status + "\t" + pickedIdsStr);

    print(name + ":  " + n + " -> " + picked.length + "  [" + status + "], rejected rows=" + rejectedRows);
    totalOut += picked.length;
}

getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
runDate = d2s(year, 0) + "-" + pad2(month + 1) + "-" + pad2(dayOfMonth) + " " +
    pad2(hour) + ":" + pad2(minute) + ":" + pad2(second);
changes = "ZMB MD HCS target curation changes\r\n";
changes += "Run\t" + runDate + "\r\n";
changes += "Targets per site\t" + targetsPerSite + "\r\n";
changes += "Minimum distance overview px\t" + minDistPx + "\r\n";
changes += "Seed\t" + seed + "\r\n";
changes += "Originals\t" + origDir + "\r\n";
changes += "Operational TargetData\t" + curDir + "\r\n";
changes += "Curated mirror\t" + auditDir + "\r\n";
changes += "\r\n";
changes += "file\tdetected_valid_rows\trejected_rows\tpicked_rows\tstatus\tpicked_object_ids\r\n";
for (i = 0; i < changeRows.length; i++) changes += changeRows[i] + "\r\n";
File.saveString(changes, auditDir + "curation_changes.csv");

// --- Summary in Log window ---
print("");
print("Sites processed: " + nFiles);
print("Files skipped:   " + skippedList.length);
print("Rows rejected:   " + totalRejectedRows);
print("Targets:         " + totalIn + " -> " + totalOut);
print("Settings:        N=" + targetsPerSite + ", min-dist=" + minDistPx + " overview px");
print("Curated:         " + curDir);
print("Curated mirror:  " + auditDir);
print("Originals:       " + origDir);

// --- Summary dialog ---
report  = "Sites processed: " + nFiles + "\n";
report += "Files skipped:   " + skippedList.length + "\n";
report += "Rows rejected:   " + totalRejectedRows + "\n";
report += "Targets:         " + totalIn + " -> " + totalOut + " (N=" + targetsPerSite + " per site)\n";
report += "Min distance:    " + minDistPx + " overview px\n";
report += "\n";
if (underList.length == 0 && skippedList.length == 0 && rejectedList.length == 0) {
    report += "All sites at full N=" + targetsPerSite + ".\n";
} else {
    if (underList.length > 0) {
        report += "Sites under target (" + underList.length + "):\n";
        for (i = 0; i < underList.length; i++) report += "  " + underList[i] + "\n";
        report += "\n";
    }
    if (skippedList.length > 0) {
        report += "Skipped files (" + skippedList.length + "):\n";
        for (i = 0; i < skippedList.length; i++) report += "  " + skippedList[i] + "\n";
        report += "\n";
    }
    if (rejectedList.length > 0) {
        report += "Files with rejected rows (" + rejectedList.length + "):\n";
        for (i = 0; i < rejectedList.length; i++) report += "  " + rejectedList[i] + "\n";
    }
}
report += "\nCurated:         " + curDir + "\n";
report += "Curated mirror:  " + auditDir + "\n";
report += "Changes:         " + auditDir + "curation_changes.csv\n";
report += "Originals:       " + origDir;
showMessage("Target curation summary", report);
