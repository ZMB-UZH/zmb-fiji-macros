# MD HCS target curation: progress and handoff

Date: 2026-07-24
Branch: `md-hcs-target-curation-v2`

## Objective

Build one reproducible Fiji/Jython target-curation workflow that works for both:

- Nico: an IN Carta analysis of a stitched whole-well montage.
- Babette: an IN Carta analysis of four overlapping acquisition sites per well.

The biological target in Babette is a nucleus overlapping both Green and Red. The
acquisition target CSV must use the schema of the channel with the highest `T` number.

## Acceptance datasets

Babette:

```text
Z:\transfer\Thom\Babette MD\26.19 NPTX2 ASO IF test_20260722_104349\experiment\Results\TriplePositive_Thom_2026-Jul-22-12-20-15-690 - Copy (2)
```

Nico:

```text
Z:\transfer\Thom\Nico MD\NB26-15_Overview10x_DAPI-mScarlet_20260713_152811\experiment_montage\Results\mScarlet Cells_2026-Jul-13-17-23-33-077 - Copy
```

Babette confirmation acquisition:

```text
Z:\transfer\Thom\Babette MD\confimation_position_placing\fijimacro2_20260724_160301
```

Testing was performed on disposable copies under `C:\tmp`; the acceptance inputs were
not modified by automated tests.

## Implemented

### Highest-T acquisition output

Selection and output roles are now independent:

- The chosen base object is used for gating and sampling.
- The channel with the highest detected `T#` supplies the output filename, header,
  metadata row and bounding-box columns.
- Operational `TargetData` contains only that acquisition channel.

Observed real-data output:

- Babette: 160 `Red_...csv` files using `T3$...BoundingBox`; no Nuclei or Green files.
- Nico: 48 `mScarlet cells_...csv` files using `T2$...BoundingBox`; no DAPI files.

### Faithful plate-overview viewer

The viewer now normalizes both representations into well coordinates:

- A montage is already well-local and remains at origin `(0, 0)`.
- Separate fields use stage positions from the nearest `image_metadata_1.csv`.
- Field offsets are converted from micrometres using the overview pixel calibration.
- Every field and FOV uses one shared scale.
- Real overlap between Babette's fields is displayed.
- Each physical well is rendered on its own clipped processor, so drawings cannot
  spill into neighbouring wells.
- Z/time combinations remain separate overview files.

The viewer has a deterministic grid fallback when acquisition metadata is unavailable,
but physical validation should use the metadata-backed view.

### Other fixes already present on the branch

- Strict parsing of the complete `R#-C#-F#-Z#-T#` site identity.
- Separate Z/T overview facets.
- Explicit `full`, `low_cells` and `constrained` statuses.
- Correct user-facing terminology: current sampling is per site/field.
- Previous SURS and ImageJ Macro string-comparison corrections.

## Evidence

### Dataset geometry

Nico:

- One field (`F0`) per well.
- Image extent approximately `8511 x 8513` pixels.
- This is a stitched whole-well montage.

Babette:

- Four fields (`F0` through `F3`) per well.
- Every field is `2304 x 2304` pixels with coordinates restarting at `(0, 0)`.
- Metadata layout:

```text
F0 upper-left    F3 upper-right
F1 lower-left    F2 lower-right
```

- Adjacent field-centre separation: approximately `2073.5` pixels.
- Adjacent-field overlap: approximately `230.5` pixels (`37.4 um`).

### FOV footprint

The FOV rectangle is expressed in analysis-image coordinates:

```text
overview camera width * overview magnification / target magnification
```

- Babette: `2304 * 40 / 60 = 1536` overview pixels.
- Nico: `2304 * 10 / 60 = 384` overview pixels.

The target camera image may still be 2304 pixels wide; these values describe its
physical footprint in the overview coordinate system.

### Confirmation acquisition

All 99 generated Babette targets were matched to the 99 stage positions recorded by
the confirmation acquisition. Predictions used the original field centre plus the
field-local target displacement.

- Median stage-position difference: `0.69 um`.
- 95th percentile: `0.78 um`.
- Maximum: `1.55 um`.

Conclusion: MetaXpress correctly interprets the CSV `field` plus field-local
coordinates. Do **not** add field offsets to acquisition CSV coordinates; that would
apply the translation twice. Field offsets belong in well-level computation and
visualization, followed by conversion back to local coordinates for output.

### Automated and visual validation

- Embedded `target_curation.py` tests: all pass.
- `test_curation_core.py`: all pass.
- `git diff --check`: clean.
- Babette: 160 highest-T CSVs and 161 readable PNGs.
- Nico: 48 highest-T CSVs and 49 readable PNGs.
- Both final plate overviews were manually inspected.

## Unresolved: cross-field duplicate sampling

The improved Babette viewer exposed a real sampling defect. Sampling and disjoint-FOV
placement currently operate independently per site. Consequently, overlap between
F0-F3 can be counted and acquired more than once.

Row 5, Column 1 is a confirmed example:

- F0: 5 acquired cells in 1 FOV.
- F1: 2 acquired cells in 2 FOVs.
- F2: 1 acquired cell in 1 FOV.
- F3: 5 acquired cells in 2 FOVs.
- An F0 and F3 acquired-cell centre pair is only `5.2` well pixels apart.
- An F0 and F1 pair is only `15.3` well pixels apart.
- Six pairs of cross-field FOV rectangles overlap physically.

The existing disjoint-FOV invariant therefore holds only within a site, not within a
physical well.

## Agreed next design

The next change should enforce two well-level invariants:

1. A physical cell enters the sampling population at most once.
2. No two requested target FOV rectangles overlap physically within a well.

Recommended pipeline:

1. Gate objects in their source-field coordinates.
2. Translate qualifying objects into well coordinates using measured field origins.
3. Partition overlap deterministically between fields (nearest field centre, stable
   field-number tie-break), so the overlap is not counted twice.
4. Add a defensive cross-field duplicate check near ownership boundaries.
5. Pool and sample once per physical well.
6. Place and collision-check all FOVs in well coordinates.
7. Assign each accepted FOV to the source field that contains it most safely.
8. Convert the FOV centre back to that field's local coordinates.
9. Write only the highest-T acquisition CSV using the correct field identity.

For Nico, this pipeline degenerates naturally to one montage field at origin `(0, 0)`.

Tests must cover:

- Babette R5-C1 duplicates collapse to one physical candidate.
- The requested count is a per-well maximum.
- No cross-field FOV pair overlaps in well coordinates.
- Output local-to-stage conversion reproduces the same well position.
- Montage behavior and existing Nico selections remain reproducible.
- Missing or inconsistent acquisition metadata fails clearly when multi-field
  normalization is required; it must not silently pretend fields are independent.

## Important scope note

The highest-T output and faithful viewer are implemented and validated. Well-level
deduplication and cross-field FOV placement are **not implemented yet**. Do not treat
the current branch as solving duplicate sampling in overlapping fields.
