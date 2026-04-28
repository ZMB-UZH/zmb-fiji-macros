// @File(label = "Analysis results folder", style = "directory") resultsDir
// @Integer(label = "Targets per site", value = 5, min = 1) targetsPerSite
// @Double(label = "Minimum distance between objects in the overview image (pixels, 0 = no constraint)", value = 340, min = 0) minDistPx

// ------------------------------------------------------------------------------
// Fiji: MD HCS target curation (v1.4.0)
// Created: 2026-04-28 | Updated: 2026-04-28
// Author: thom.dehoog@zmb.uzh.ch | ZMB Center for Microscopy and Image Analysis, UZH
//
// If you publish a paper using this macro, please acknowledge.
//
// Description: Curate MetaXpress TargetData CSVs produced by Encarta segmentation
//   so a downstream targeted acquisition images N cells per site instead of all
//   detected cells. Picks are random per site with a Euclidean minimum-distance
//   constraint between bounding-box centres so re-imaging tiles stay separated.
//   Distances operate directly in the overview pixel units already present in
//   the TargetData CSVs; no calibration file is required.
//
// Behaviour:
//   - First run renames TargetData/ -> TargetData_original/ (idempotent).
//   - Every run reads from TargetData_original/ and rewrites TargetData/ from
//     scratch. The original Encarta output is never modified.
//   - SummaryInfo, ObjectData, FieldData, WellData are left untouched.
//   - At the end, a summary dialog reports settings used and any sites that
//     ended up under target (with reason: low_cells or constrained).
//
// Expected input:
//   <resultsDir>/
//     TargetData/                                 (Encarta output; renamed on first run)
//       Cells_singleTargetData_R<r>-C<c>-F<f>-Z<z>-T<t>.csv
//
// Output:
//   <resultsDir>/
//     TargetData/                                 (curated; MetaXpress reads from here)
//     TargetData_original/                        (preserved Encarta output)
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
//   - The macro requires the object_id and BoundingBox columns in each
//     TargetData CSV; sites without them are skipped with a warning.
// ------------------------------------------------------------------------------

// --- Helpers ---

// Find a column index by exact header name; returns -1 if absent.
function findColumn(headerCols, name) {
    for (i = 0; i < headerCols.length; i++) {
        c = headerCols[i];
        if (endsWith(c, "\r")) c = substring(c, 0, lengthOf(c)-1);
        if (c == name) return i;
    }
    return -1;
}

// Strip a trailing \r so split-by-\n on CRLF text yields clean rows.
function stripCR(s) {
    if (endsWith(s, "\r")) return substring(s, 0, lengthOf(s)-1);
    return s;
}

// Classify a site outcome based on detected and picked counts.
// Returns "full", "low_cells", or "constrained".
function siteStatus(nDetected, nPicked, target) {
    if (nPicked == target) return "full";
    if (nDetected < target) return "low_cells";
    return "constrained";
}

// --- Resolve paths and set up output folders ---
if (!File.isDirectory(resultsDir))
    exit("Not a directory: " + resultsDir);
if (!endsWith(resultsDir, File.separator)) resultsDir += File.separator;

origDir = resultsDir + "TargetData_original" + File.separator;
curDir  = resultsDir + "TargetData" + File.separator;

// First run: rename Encarta's TargetData/ -> TargetData_original/.
if (!File.isDirectory(origDir)) {
    if (!File.isDirectory(curDir))
        exit("No TargetData/ folder at: " + curDir);
    if (!File.rename(curDir, origDir))
        exit("Could not rename TargetData/ -> TargetData_original/.");
    print("First run: renamed TargetData/ -> TargetData_original/");
}

// Refresh curated output folder (always overwritten from TargetData_original/).
if (File.isDirectory(curDir)) {
    stale = getFileList(curDir);
    for (i = 0; i < stale.length; i++) File.delete(curDir + stale[i]);
} else {
    File.makeDirectory(curDir);
}

// Fixed seed for reproducible curation. Edit to draw a different sample.
seed = 42;
random("seed", seed);

// --- Curate each site ---
list = getFileList(origDir);
totalIn = 0;
totalOut = 0;
nFiles = 0;
underList = newArray(0);  // human-readable lines for sites that ended up under-target

for (f = 0; f < list.length; f++) {
    name = list[f];
    if (!startsWith(name, "Cells_singleTargetData_")) continue;
    if (!endsWith(name, ".csv")) continue;

    raw = File.openAsString(origDir + name);
    lines = split(raw, "\n");
    if (lines.length < 2) continue;

    header = stripCR(lines[0]);
    cols = split(header, ",");
    iBX  = findColumn(cols, "T1$AS_FID_Blob_BoundingBoxX");
    iBY  = findColumn(cols, "T1$AS_FID_Blob_BoundingBoxY");
    iBW  = findColumn(cols, "T1$AS_FID_Blob_BoundingBoxWidth");
    iBH  = findColumn(cols, "T1$AS_FID_Blob_BoundingBoxHeight");
    iOID = findColumn(cols, "object_id");
    if (iBX < 0 || iBY < 0 || iBW < 0 || iBH < 0) {
        print("skip (missing bounding-box columns): " + name);
        continue;
    }
    if (iOID < 0) {
        print("skip (missing object_id column): " + name);
        continue;
    }

    rows = newArray(0);
    cx   = newArray(0);
    cy   = newArray(0);
    oid  = newArray(0);
    for (i = 1; i < lines.length; i++) {
        row = stripCR(lines[i]);
        if (lengthOf(row) == 0) continue;
        flds = split(row, ",");
        if (flds.length < cols.length) continue;
        x = parseFloat(flds[iBX]);
        y = parseFloat(flds[iBY]);
        w = parseFloat(flds[iBW]);
        h = parseFloat(flds[iBH]);
        if (isNaN(x) || isNaN(y) || isNaN(w) || isNaN(h)) continue;
        rows = Array.concat(rows, row);
        cx   = Array.concat(cx, x + w/2.0);
        cy   = Array.concat(cy, y + h/2.0);
        oid  = Array.concat(oid, flds[iOID]);
    }

    n = rows.length;
    totalIn += n;
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

    siteLabel = substring(name, 0, lengthOf(name) - 4);
    if (status != "full")
        underList = Array.concat(underList,
            siteLabel + ":  " + n + " detected, " + picked.length + " picked  (" + status + ")");

    print(name + ":  " + n + " -> " + picked.length + "  [" + status + "]");
    totalOut += picked.length;
}

// --- Summary in Log window ---
print("");
print("Sites processed: " + nFiles);
print("Targets:         " + totalIn + " -> " + totalOut);
print("Settings:        N=" + targetsPerSite + ", min-dist=" + minDistPx + " overview px");
print("Curated:         " + curDir);
print("Originals:       " + origDir);

// --- Summary dialog ---
report  = "Sites processed: " + nFiles + "\n";
report += "Targets:         " + totalIn + " -> " + totalOut + " (N=" + targetsPerSite + " per site)\n";
report += "Min distance:    " + minDistPx + " overview px\n";
report += "\n";
if (underList.length == 0) {
    report += "All sites at full N=" + targetsPerSite + ".\n";
} else {
    report += "Sites under target (" + underList.length + "):\n";
    for (i = 0; i < underList.length; i++) report += "  " + underList[i] + "\n";
}
report += "\nCurated:    " + curDir + "\n";
report += "Originals:  " + origDir;
showMessage("Target curation summary", report);
