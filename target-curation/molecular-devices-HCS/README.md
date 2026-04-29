# MD HCS target curation

When you run an IN Carta analysis on a Molecular Devices ImageXpress (MD HCS) overview, it can detect dozens of cells per well — usually more than you actually want to re-image at high magnification. This macro thins each site's list down to **N cells per site** (default 5), picked at random, and skips any cell whose re-imaging tile would overlap one that has already been picked. The original IN Carta output is preserved, so you can re-run the macro as often as you like.

## What it does

- Picks N cells per site at random (default 5).
- Skips cells that would land too close to one already picked, so the high-magnification tiles don't overlap.
- Saves the original IN Carta output once, the first time you run it, and re-curates from that copy on every subsequent run.
- Writes a `TargetData_curated/` mirror of the curated `TargetData/` folder, including a `curation_changes.csv` audit file.
- Refuses to overwrite `TargetData/` on rerun unless it matches the previous `TargetData_curated/` mirror.
- Shows a summary at the end listing any sites that didn't reach N, with the reason.
- Reports malformed/rejected rows in the final summary and audit file.
- Does not touch the `SummaryInfo`, `ObjectData`, `FieldData`, or `WellData` CSVs.

## How the picking works

For each site:

1. Read all the cells IN Carta detected.
2. Look at them in a random order.
3. Accept each cell that is far enough from every cell already picked.
4. Stop once N cells are picked, or when no more cells fit.

If you set the minimum distance to `0`, the spacing check is turned off and the macro simply picks the first N cells in random order.

## Choosing the minimum distance

The minimum distance is given in **pixels of the overview image**, the same units the bounding-box columns in the CSVs already use. Two cells closer than this distance are treated as overlapping — only the first one to be picked is kept.

The default of `340 px` is sized for a typical ZMB setup:

- 10× overview camera at about 0.665 µm/px.
- High-magnification acquisition at 60× with a 1.5× lens changer on a 2048 × 2048 sCMOS with 6.5 µm pixels.
- That gives a re-imaging tile of about 148 µm per side, which is roughly 222 px on the overview, or about 316 px along the diagonal. The default adds a small comfort margin on top.

If your setup is different, work it out as:

```
tile_side_um    = imaging_pixel_size_um * imaging_image_width_px
tile_side_px    = tile_side_um / overview_pixel_size_um
tile_diagonal   = tile_side_px * sqrt(2)
min_distance_px = tile_diagonal * margin   (margin >= 1.0; 1.5 leaves a half-tile gap)
```

## What "under target" means

After the run, each site is reported in the summary dialog as one of:

| Status | Meaning |
|---|---|
| `full` | Got exactly N cells. |
| `low_cells` | IN Carta found fewer than N cells in the first place, so the macro just kept all of them. |
| `constrained` | IN Carta found enough cells, but the minimum-distance rule wouldn't let the macro fit N of them. Lower the minimum distance, or accept the smaller picture. |

## Where the files go

Input — the analysis results folder IN Carta produced, e.g. `experiment/Results/Brightest Nuclei_<timestamp>/` or `.../Largest Cells_<timestamp>/`:

```
<resultsDir>/
  TargetData/                                     (IN Carta's output; renamed on first run)
    Cells_singleTargetData_R<r>-C<c>-F<f>-Z<z>-T<t>.csv
```

Output:

```
<resultsDir>/
  TargetData/                                     (curated list; MetaXpress reads from here)
    Cells_singleTargetData_R<r>-C<c>-F<f>-Z<z>-T<t>.csv     (filtered to N rows)
  TargetData_curated/                             (curated mirror for inspection/auditing)
    Cells_singleTargetData_R<r>-C<c>-F<f>-Z<z>-T<t>.csv     (same curated CSVs as TargetData/)
    curation_changes.csv                          (tab-delimited settings plus per-file picked/skipped/rejected summary)
  TargetData_original/                            (IN Carta's original output, untouched)
    Cells_singleTargetData_R<r>-C<c>-F<f>-Z<z>-T<t>.csv
```

## How to run it

1. Open [Fiji](https://fiji.sc).
2. Drag and drop `ZMB_MD_HCS_target_curation.ijm` onto the Fiji toolbar, or open it from `Plugins > Macros > Run...`.
3. In the dialog, pick the analysis results folder (e.g. `Brightest Nuclei_<timestamp>/`).
4. Set how many cells you want per site (default 5).
5. Set the minimum distance between cells, in overview pixels (default 340, or `0` to turn the spacing check off).
6. Run. A summary dialog at the end tells you how many cells were picked per site, and flags any that didn't reach N.

The picks are produced from a fixed random seed (42) in the script, so re-running with the same settings gives you the same picks. To draw a different random sample, edit the `seed` constant near the top of the macro.

On rerun, the macro only overwrites `TargetData/` if it can verify that `TargetData/` still matches the previous `TargetData_curated/` mirror. If IN Carta has regenerated a new `TargetData/` in the same results folder, the macro aborts instead of deleting it.

## Before trusting it in production

The macro replaces the contents of `TargetData/` with the curated list, but the workflow only works if MetaXpress actually re-reads `TargetData/` at the moment the targeted acquisition is queued — not earlier, when IN Carta finished writing. To check: open one of the curated CSVs, hand-edit it to obviously fewer rows, queue the targeted acquisition, and confirm MetaXpress images only the smaller set. If it does, the workflow is reliable.
