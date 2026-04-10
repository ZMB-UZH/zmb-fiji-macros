# MD HCS stack assembler

Assembles multi-channel stacks from individual TIF tiles exported by Molecular Devices ImageXpress (MD HCS). Outputs hyperstacks and/or maximum intensity projections.

## Features

- Auto-detects flat folders or `timepoint` subfolder structures
- Groups tiles by well prefix from filename tags (`_w`, `_z`, `_s`)
- Robust tag matching avoids false hits (e.g. `_s` inside `_seeding`)
- Preserves spatial calibration and acquisition metadata (see below)

## Metadata preserved

The following metadata is read from the MetaSeries XML embedded in each tile and carried over to the output files:

**Spatial calibration** (applied via `Image > Properties`):
- Pixel size XY (from `spatial-calibration-x/y`)
- Z-step (computed from `z-position` of consecutive slices)
- Calibration unit

**Acquisition metadata** (stored in image info, visible via `Image > Show Info...`):
- Instrument and software version
- Objective, numerical aperture, refractive index
- Channel names, wavelengths, and exposure times
- Acquisition timestamp

This metadata is preserved in the TIFF files and readable by both Fiji (drag and drop) and Bio-Formats.

## Expected input filename patterns

```
prefix_w0_z0.tif
prefix_s0_w0_z0.tif
```

## Output

- `outputDir/hyperstacks/` — multi-channel, multi-z-slice hyperstacks
- `outputDir/MIPs/` — maximum intensity projections (multi-channel composites)

## Usage

1. Open [Fiji](https://fiji.sc)
2. Drag and drop `ZMB_MD_HCS_stack_assembler.ijm` onto the Fiji toolbar, or open via `Plugins > Macros > Run...`
3. Select the input folder containing TIF tiles (flat or with `timepoint` subfolders)
4. Select an output folder
5. Tick **Hyperstack** and/or **Maximum Intensity Projection**
