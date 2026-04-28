# MD HCS target curation

Curates the per-site `Cells_singleTargetData_*.csv` files produced by an Encarta segmentation analysis on Molecular Devices ImageXpress (MD HCS), so that a downstream targeted acquisition images **N** cells per site instead of all detected cells. Picks are random per site with a Euclidean minimum-distance constraint between bounding-box centres so re-imaging tiles stay separated.

## Features

- Random N-per-site curation with optional minimum-distance constraint (in overview pixels)
- Idempotent first-run rename of `TargetData/` to `TargetData_original/`; every subsequent run re-curates from the preserved original
- Per-site outcome classification (`full`, `low_cells`, `constrained`)
- End-of-run summary dialog that surfaces any sites that ended up under target
- Original `SummaryInfo`, `ObjectData`, `FieldData`, `WellData` files are never modified

## Picking algorithm

For each site:

1. Load all detected cells from `TargetData_original/Cells_singleTargetData_R<r>-C<c>-F<f>-Z<z>-T<t>.csv`.
2. Compute each cell's centre as the bounding-box centre (`BoundingBoxX + Width/2`, `BoundingBoxY + Height/2`).
3. Fisher-Yates shuffle the cell indices.
4. Walk through the shuffled list and accept each candidate whose centre is at least the minimum distance away (Euclidean) from every already-picked centre.
5. Stop when N picks are made or no more candidates fit.

A minimum distance of `0` disables the constraint entirely (operator gets the first N from the random shuffle). Useful as a baseline run for comparison.

## Choosing the minimum distance

The bounding-box columns in the TargetData CSVs are already in **overview pixels**, so the minimum distance is set in the same units. No calibration is needed.

To choose a value, work backwards from the imaging tile size at the target acquisition magnification:

```
tile_side_um    = imaging_pixel_size_um * imaging_image_width_px
tile_side_px    = tile_side_um / overview_pixel_size_um
tile_diagonal   = tile_side_px * sqrt(2)
min_distance_px = tile_diagonal * margin   (margin >= 1.0; 1.5 = comfortable)
```

Example: 60× objective + 1.5× lens changer + 2048 × 2048 sCMOS with 6.5 µm pixels gives an imaging pixel size of `6.5 / (60 * 1.5) ≈ 0.072 µm/px` and a tile side of `0.072 * 2048 ≈ 148 µm`. On a 10× MD HCS overview at ~0.665 µm/px, that's `148 / 0.665 ≈ 222 px` per tile side, `~316 px` along the diagonal. The default of `340 px` adds a small comfort gap on top.

## Per-site outcome

Each curated site is classified and reported in the end-of-run summary dialog:

| Status | Meaning |
|---|---|
| `full` | Picked exactly N cells |
| `low_cells` | The site had fewer than N cells to begin with |
| `constrained` | The site had ≥ N cells, but the minimum-distance constraint forbade more picks |

## Expected input layout

```
<resultsDir>/
  TargetData/                                     (Encarta output; renamed on first run)
    Cells_singleTargetData_R<r>-C<c>-F<f>-Z<z>-T<t>.csv
```

Where `<resultsDir>` is the per-analysis folder produced by an Encarta analysis (e.g. `experiment/Results/Brightest Nuclei_<timestamp>/` or `.../Largest Cells_<timestamp>/`).

## Output

```
<resultsDir>/
  TargetData/                                     (curated; MetaXpress reads from here)
    Cells_singleTargetData_R<r>-C<c>-F<f>-Z<z>-T<t>.csv     (filtered to N rows)
  TargetData_original/                            (preserved Encarta output)
    Cells_singleTargetData_R<r>-C<c>-F<f>-Z<z>-T<t>.csv     (untouched)
```

## Usage

1. Open [Fiji](https://fiji.sc)
2. Drag and drop `ZMB_MD_HCS_target_curation.ijm` onto the Fiji toolbar, or open via `Plugins > Macros > Run...`
3. In the dialog, select the analysis results folder (e.g. `Brightest Nuclei_<timestamp>/`)
4. Set the targets per site (default 5)
5. Set the minimum distance between picks in overview pixels (default 340, `0` to disable)
6. Run; the summary dialog reports settings, totals, and any under-target sites

The random seed is fixed at 42 in source so runs are reproducible. Edit the `seed` constant near the top of the script to vary.

## Verifying the workflow

Before relying on this in production, confirm that MetaXpress reads the `TargetData/` CSVs at queue time (when the targeted acquisition is started), not at Encarta-write time (in which case it would have cached the original, larger set into its database before the curation ran). To check: edit one CSV by hand, queue the targeted acquisition, and confirm the smaller cell set is imaged.
