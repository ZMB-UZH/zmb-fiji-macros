# zmb-fiji-macros

Fiji macros by the [ZMB Center for Microscopy and Image Analysis](https://www.zmb.uzh.ch), University of Zurich.

## Stack assemblers

### [zmb_MD_HCS_stack_assembler.ijm](stack-assemblers/molecular-devices-HCS/zmb_MD_HCS_stack_assembler.ijm)

Assembles multi-channel stacks from individual TIF tiles exported by Molecular Devices ImageXpress (MD HCS). Outputs hyperstacks and/or maximum intensity projections.

**Features**
- Auto-detects flat folders or `timepoint` subfolder structures
- Groups tiles by well prefix from filename tags (`_w`, `_z`, `_s`)
- Robust tag matching avoids false hits (e.g. `_s` inside `_seeding`)

**Expected input filename patterns**
```
prefix_w0_z0.tif
prefix_s0_w0_z0.tif
```

**Output**
- `outputDir/hyperstacks/` — multi-channel, multi-z-slice hyperstacks
- `outputDir/MIPs/` — maximum intensity projections (multi-channel composites)

## Usage

1. Open [Fiji](https://fiji.sc)
2. Drag and drop the `.ijm` file onto the Fiji toolbar, or open via `Plugins > Macros > Run...`
3. A dialog will prompt for input/output folders and output options
