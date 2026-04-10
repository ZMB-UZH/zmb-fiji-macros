// @File(label = "Input folder", style = "directory") inputDir
// @File(label = "Output folder", style = "directory") outputDir
// @Boolean(label = "Hyperstack", value = true) saveHyperstack
// @Boolean(label = "Maximum Intensity Projection", value = false) saveMIP

// ------------------------------------------------------------------------------
// Fiji: MD HCS stack assembler (v1.1.0)
// Created: 2026-04-10 | Updated: 2026-04-10
// Author: thom.dehoog@zmb.uzh.ch | ZMB Center for Microscopy and Image Analysis, UZH
//
// If you publish a paper using this macro, please acknowledge.
//
// Description: Assembles multi-channel stacks from individual TIF tiles exported
//   by Molecular Devices ImageXpress (MD HCS).
//   - Auto-detects flat folders or timepoint subfolder structure
//   - Groups tiles by well prefix from filename tags (_w, _z, _s)
//   - Preserves spatial calibration (XY pixel size, z-step) from tile metadata
//   - Preserves acquisition metadata (objective, NA, exposure, channel names)
//   - Outputs hyperstacks and/or maximum intensity projections
//   - Robust tag detection avoids false matches (e.g. "_s" in "_seeding")
//   - Runs in batch mode to minimise memory overhead
//
// Expected filename patterns:
//   prefix_w0_z0.tif
//   prefix_s0_w0_z0.tif
//
// Input:  Folder containing TIF tiles, or folder with timepoint subfolders
// Output: outputDir/hyperstacks/ and/or outputDir/MIPs/
// ------------------------------------------------------------------------------

// --- Validate ---
if (!saveHyperstack && !saveMIP)
    exit("Please select at least one output option (Hyperstack or Maximum Intensity Projection).");

run("Close All");

if (!endsWith(inputDir, File.separator)) inputDir += File.separator;
if (!endsWith(outputDir, File.separator)) outputDir += File.separator;

// --- Create output directories ---
hsDir = outputDir + "hyperstacks" + File.separator;
mipDir = outputDir + "MIPs" + File.separator;
if (saveHyperstack && !File.isDirectory(hsDir))
    File.makeDirectory(hsDir);
if (saveMIP && !File.isDirectory(mipDir))
    File.makeDirectory(mipDir);

// --- Helpers: filename tag parsing ---

// Find position of a tag (e.g. "_s") that is followed by a digit.
// Skips false matches like "_s" inside "_seeding". Returns -1 if not found.
function findTagPos(filename, tag) {
    from = 0;
    while (from < lengthOf(filename)) {
        idx = indexOf(filename, tag, from);
        if (idx < 0) return -1;
        next = idx + lengthOf(tag);
        if (next < lengthOf(filename) && charCodeAt(filename, next) >= 48 && charCodeAt(filename, next) <= 57)
            return idx;
        from = idx + 1;
    }
    return -1;
}

// Extract a numeric index following a tag like "_w", "_z", "_s".
// Returns -1 if the tag (followed by a digit) is not found.
function extractIndex(filename, tag) {
    pos = findTagPos(filename, tag);
    if (pos < 0) return -1;
    start = pos + lengthOf(tag);
    end = start;
    while (end < lengthOf(filename) && charCodeAt(filename, end) >= 48 && charCodeAt(filename, end) <= 57)
        end++;
    return parseInt(substring(filename, start, end));
}

// Extract the well/group prefix before the first _s or _w tag.
// Result stored in global _groupKey (ImageJ macro cannot return strings).
var _groupKey = "";

function extractGroupKey(filename) {
    sPos = findTagPos(filename, "_s");
    wPos = findTagPos(filename, "_w");
    if (sPos >= 0 && (wPos < 0 || sPos < wPos)) {
        _groupKey = substring(filename, 0, sPos);
        return 1;
    }
    if (wPos >= 0) {
        _groupKey = substring(filename, 0, wPos);
        return 1;
    }
    _groupKey = "";
    return 0;
}

// --- Helpers: MetaSeries XML metadata parsing ---

// Parse value attribute of <prop id="propId" ... value="..."/> from MetaSeries XML.
// Result stored in global _xmlValue. Returns 1 if found, 0 otherwise.
var _xmlValue = "";

function parseXmlProp(info, propId) {
    needle = "id=\"" + propId + "\"";
    idx = indexOf(info, needle);
    if (idx < 0) { _xmlValue = ""; return 0; }
    valNeedle = "value=\"";
    valIdx = indexOf(info, valNeedle, idx);
    if (valIdx < 0 || valIdx > idx + 200) { _xmlValue = ""; return 0; }
    valStart = valIdx + lengthOf(valNeedle);
    valEnd = indexOf(info, "\"", valStart);
    _xmlValue = substring(info, valStart, valEnd);
    return 1;
}

// --- Detect folder structure: flat tiles or timepoint subfolders ---
contents = getFileList(inputDir);
numTP = 0;
tpNames = newArray(1000);
for (i = 0; i < contents.length; i++) {
    if (startsWith(contents[i], "timepoint") && endsWith(contents[i], "/")) {
        tpNames[numTP] = replace(contents[i], "/", "");
        numTP++;
    }
}

if (numTP > 0) {
    tpNames = Array.trim(tpNames, numTP);
    Array.sort(tpNames);
    tpDirs = newArray(numTP);
    for (i = 0; i < numTP; i++)
        tpDirs[i] = inputDir + tpNames[i] + File.separator;
    print("Found " + numTP + " timepoint(s)");
} else {
    numTP = 1;
    tpDirs = newArray(1);
    tpNames = newArray(1);
    tpDirs[0] = inputDir;
    tpNames[0] = "";
}

print("Input:  " + inputDir);
print("Output: " + outputDir);
stackCount = 0;

// --- Process ---
setBatchMode(true);

for (tp = 0; tp < numTP; tp++) {
    dir = tpDirs[tp];
    list = getFileList(dir);

    if (tpNames[tp] != "")
        print("Processing " + tpNames[tp] + " (" + list.length + " files)...");
    else
        print("Processing " + list.length + " files...");

    // Collect unique group keys (one per well)
    maxKeys = 0;
    tempKeys = newArray(1000);
    for (i = 0; i < list.length; i++) {
        if (!endsWith(list[i], ".tif")) continue;
        if (extractIndex(list[i], "_w") < 0) continue;
        if (extractIndex(list[i], "_z") < 0) continue;
        extractGroupKey(list[i]);
        key = _groupKey;
        if (key == "") continue;
        found = false;
        for (g = 0; g < maxKeys; g++) {
            if (tempKeys[g] == key) found = true;
        }
        if (!found) {
            tempKeys[maxKeys] = key;
            maxKeys++;
        }
    }

    if (maxKeys == 0) {
        print("  No matching TIF files, skipping.");
        continue;
    }

    keys = Array.trim(tempKeys, maxKeys);
    Array.sort(keys);

    // Process each well group
    for (k = 0; k < keys.length; k++) {
        key = keys[k];
        maxW = -1;
        maxZ = -1;
        maxS = -1;
        hasSites = false;

        for (i = 0; i < list.length; i++) {
            if (!endsWith(list[i], ".tif")) continue;
            extractGroupKey(list[i]);
            if (_groupKey != key) continue;
            w = extractIndex(list[i], "_w");
            z = extractIndex(list[i], "_z");
            if (w < 0 || z < 0) continue;
            if (w > maxW) maxW = w;
            if (z > maxZ) maxZ = z;
            s = extractIndex(list[i], "_s");
            if (s >= 0) {
                hasSites = true;
                if (s > maxS) maxS = s;
            }
        }

        numC = maxW + 1;
        numZ = maxZ + 1;
        numS = 1;
        if (hasSites) numS = maxS + 1;

        // Process each site
        for (s = 0; s < numS; s++) {

            // --- Open all tiles, reading metadata inline ---
            pixelWidth = 0;
            pixelHeight = 0;
            calUnit = "pixel";
            zPos0 = 0;
            zPos1 = 0;
            objective = "";
            na = "";
            ri = "";
            software = "";
            acqTime = "";
            channelInfo = "";

            for (w = 0; w < numC; w++) {
                for (z = 0; z < numZ; z++) {
                    if (hasSites)
                        filename = key + "_s" + s + "_w" + w + "_z" + z + ".tif";
                    else
                        filename = key + "_w" + w + "_z" + z + ".tif";
                    if (!File.exists(dir + filename))
                        exit("Missing file: " + filename);
                    open(dir + filename);

                    // Read metadata from the active tile before opening the next
                    tileInfo = getMetadata("Info");

                    // First tile (w0, z0): XY calibration + acquisition metadata
                    if (w == 0 && z == 0) {
                        if (parseXmlProp(tileInfo, "spatial-calibration-x")) pixelWidth = parseFloat(_xmlValue);
                        if (parseXmlProp(tileInfo, "spatial-calibration-y")) pixelHeight = parseFloat(_xmlValue);
                        if (parseXmlProp(tileInfo, "spatial-calibration-units")) calUnit = _xmlValue;
                        if (parseXmlProp(tileInfo, "z-position")) zPos0 = parseFloat(_xmlValue);
                        if (parseXmlProp(tileInfo, "HCS.ai Objective")) objective = _xmlValue;
                        if (parseXmlProp(tileInfo, "_MagNA_")) na = _xmlValue;
                        if (parseXmlProp(tileInfo, "_MagRI_")) ri = _xmlValue;
                        if (parseXmlProp(tileInfo, "ApplicationVersion")) software = _xmlValue;
                        if (parseXmlProp(tileInfo, "acquisition-time-local")) acqTime = _xmlValue;
                    }

                    // Second z-slice of first channel: z-step
                    if (w == 0 && z == 1) {
                        if (parseXmlProp(tileInfo, "z-position")) zPos1 = parseFloat(_xmlValue);
                    }

                    // First z-slice of each channel: channel name, wavelength, exposure
                    if (z == 0) {
                        chName = "";
                        chWavelength = "";
                        chExposure = "";
                        parseXmlProp(tileInfo, "_IllumSetting_"); chName = _xmlValue;
                        parseXmlProp(tileInfo, "wavelength"); chWavelength = _xmlValue;
                        parseXmlProp(tileInfo, "Exposure Time"); chExposure = _xmlValue;
                        channelInfo += "Channel " + w + ": " + chName;
                        if (chWavelength != "") channelInfo += " (" + chWavelength + " nm)";
                        if (chExposure != "") channelInfo += ", " + chExposure;
                        channelInfo += "\n";
                    }
                }
            }

            // Compute z-step
            zStep = 0;
            if (numZ > 1) zStep = abs(zPos1 - zPos0);

            // Build metadata summary
            metaSummary = "Instrument: Molecular Devices ImageXpress\n";
            if (software != "") metaSummary += "Software: MetaXpress " + software + "\n";
            if (objective != "") metaSummary += "Objective: " + objective + "\n";
            if (na != "") metaSummary += "NA: " + na + "\n";
            if (ri != "") metaSummary += "Refractive index: " + ri + "\n";
            metaSummary += "Pixel size: " + pixelWidth + " x " + pixelHeight + " " + calUnit + "\n";
            if (numZ > 1) metaSummary += "Z-step: " + zStep + " " + calUnit + "\n";
            if (acqTime != "") metaSummary += "Acquisition time: " + acqTime + "\n";
            metaSummary += "\n" + channelInfo;

            // --- Assemble ---
            run("Images to Stack", "use");
            run("Stack to Hyperstack...",
                "order=xyzct channels=" + numC +
                " slices=" + numZ +
                " frames=1 display=Composite");
            hsID = getImageID();

            // Apply calibration
            if (pixelWidth > 0 && pixelHeight > 0) {
                propString = "unit=" + calUnit + " pixel_width=" + pixelWidth + " pixel_height=" + pixelHeight;
                if (zStep > 0) propString += " voxel_depth=" + zStep;
                run("Properties...", propString);
            }

            // Attach acquisition metadata
            setMetadata("Info", metaSummary);

            // Build output name: preserve original prefix, prepend timepoint if present
            if (hasSites)
                saveName = key + "_s" + s;
            else
                saveName = key;
            if (tpNames[tp] != "")
                saveName = tpNames[tp] + "_" + saveName;

            // --- Save outputs ---
            if (saveHyperstack) {
                selectImage(hsID);
                saveAs("Tiff", hsDir + saveName + ".tif");
                print("  Hyperstack: " + saveName + ".tif");
            }

            if (saveMIP) {
                selectImage(hsID);
                run("Z Project...", "projection=[Max Intensity]");
                setMetadata("Info", metaSummary);
                saveAs("Tiff", mipDir + saveName + ".tif");
                close();
                print("  MIP: " + saveName + ".tif");
            }

            selectImage(hsID);
            close();
            stackCount++;
        }
    }
}

setBatchMode(false);
print("Done — " + stackCount + " stack(s) processed.");
