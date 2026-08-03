# Well-level sampling: change plan

Date: 2026-07-31
Branch: `md-hcs-target-curation-v2`
Target: `target_curation.py` (single file, sections 1-13)

## Status - 2026-08-03: implemented

The plan below is carried out. Five commits on the branch, none pushed yet:

| Commit | What |
|---|---|
| `990ab16` | Phase 0 - point the goldens at a path that exists, so the control runs |
| `81ab09d` | Blocker 4 - the SURS clamp that under-sampled the far edge |
| `a40bacf` | Blocker 3 - give up a target only when its whole window is covered |
| `a24fa3e` | Phase 1 - `field_spans`, `owning_field`, `pool_well` + tests, nothing wired |
| `26f9012` | Phase 2 - `curate()` samples the well; callers and reports follow |
| `ef0868b` | Report the recorded layout and its overlap, and warn on implausible values |
| `2a4cc98` | Reach the requested count instead of dropping cells the sample chose |
| `6535fd4` | Remove the frame-mate swap: nothing downstream may reconsider the sample |

Departures from the plan, all deliberate:

- `field_origins` became **`field_spans`** and carries each field's size as well as its
  origin, so `owning_field` can require containment. `field_wh` as a separate injected
  parameter is gone; `curate(..., spans=...)` is the one channel.
- The **jitter guard** (blocker 1) is a second pass inside `pool_well`, not a separate
  function, and is bucketed on a grid coarser than any cell - the every-pair form is
  quadratic in a population of tens of thousands of nuclei per well.
- `fov_home_field` was not needed: `owning_field` on the FOV centre is the same
  question.
- Phase 3 could **not** be run as written. The Babette plate and the Nico plate are
  both gone from `Z:\transfer\Thom\`. See below.

## What replaced Phase 3

The byte-identical gate was run, on the plate that is still on the share
(`10306\FIJI_Target_Curation_Test\...\curation test_2026-Apr-29-15-23-15-269`, six
wells, one stitched montage each): curated with the code from `a24fa3e` and from
`26f9012`, comparing per-site cells, tiles, statuses and the whole generated FOV rows,
over four combinations of FOV size, sample size and seed. Identical throughout. Its
numbers are pinned as `CONTROL_*` and that block now runs rather than skipping.

The four-field case is covered by a **synthetic** well - 2x2 fields at 10% overlap,
cells planted in a field interior, in a seam and in the four-way corner. It bites: on
the same fixture `a24fa3e` selects 10 cells for 6, acquires 10 against a request of 3,
and places 7 overlapping cross-field pairs, against 6, 3 and 0 now.

Still open, and needing a real multi-field plate:

- No real four-field acceptance run. `MULTIFIELD_RESULTS` still names the Babette path
  and auto-skips; point it at the next multi-field plate on the share and the three
  invariants run against real geometry.
- The jitter guard's radius criterion has never met real segmentation disagreement.
- Whether MetaXpress accepts a well whose fields yield header-only CSVs (a pooled
  sample of 5 over 4 fields can leave a field with no FOV) is still unknown.
- `GenericDialog` and the SciJava widgets remain GUI-only, as before.

## The defect

`curate()` loops over **sites** (`R#-C#-F#-Z#-T#`) and treats each as a well. With one
stitched montage per well that is correct: site is well, and every invariant holds. With
four acquisition fields per well at 10% overlap it fails three ways:

- `sample_size` is requested per field, so a well asks for up to `4 x n`.
- A cell in the 230.5 px overlap band is a separate object in each field's CSV, so it can
  be gated, sampled and acquired more than once.
- `_place_disjoint` compares FOV rectangles in field-local coordinates, so two FOVs from
  different fields are never checked against each other. R5-C1 has six such physically
  overlapping pairs.

The well-coordinate machinery to fix this already exists: `read_field_origins()` converts
each field's stage position into a well-local pixel origin. It is called from exactly one
place - `_render_plate_overview` - so the viewer knows the fields overlap and the
selection path does not.

## Design decisions (settled)

1. **Ownership by partition, not by matching.** Candidates for an object are the fields
   whose image contains it; the owner is the candidate whose field centre is nearest, ties
   broken by lowest field number. On a regular grid this is exactly "trim half the overlap
   off every shared edge" - 10% overlap trims 5% per side - but it also survives an
   irregular layout, an unequal overlap, or a gap between tiles, where the fixed trim stops
   being a partition. It never has to decide that two detections *are* the same cell: it
   decides by position, so segmentation disagreement between fields is irrelevant.

2. **The well is the sampling unit; `sample_size` counts target cells.** One SURS grid of
   ~n frames over the well, one random offset, one cell per occupied frame.

3. **Order is forced: partition -> pool -> grid -> draw.** A frame on a seam contains band
   cells that appear twice in the pooled list, so each holds two tickets in that frame's
   draw while an interior cell holds one - the seams get over-sampled, which is the exact
   density bias SURS exists to remove. Deduplicating after sampling is worse: it biases
   *and* returns fewer than n.

4. ~~**Strict FOV disjointness stays.**~~ **SUPERSEDED 2026-08-03 by `2a4cc98` +
   `6535fd4`.** The operator's rule is that the requested count must be reached whenever
   it physically can be and fall short only when the well runs out of cells, and that
   minimum overlap is a requirement but *having the correct objects sampled matters
   more*. So a sampled cell that cannot be placed clear is imaged anyway, at the
   least-overlapping point of its window, rather than dropped. `sample_size` remains a
   hard maximum.

   **The governing rule: optimise the placement, never the sample.** `2a4cc98` also
   swapped to another eligible cell of the same SURS frame before accepting overlap.
   That was wrong and `6535fd4` removed it: whether a cell can be placed clear is a
   function of how close its neighbours are, so choosing between cells on it makes
   inclusion depend on local density - the same defect as oversampling and keeping
   whichever cells placed best, only smaller. Confining the choice to one frame limits
   how far it spreads; it does not make it a sample. `_grow`, which merges cells that
   were *already* sampled when they happen to share a FOV, is post-selection and
   therefore free. Cost of the strict rule, measured: 3 overlapping pairs on the
   four-field fixture instead of 2, with 30/30 cells acquired either way.

   The paragraph below is kept for the reasoning it records:

   *This reverses an earlier suggestion in discussion.*
   The goal is fewest acquisitions and least data; disjointness is that goal in its
   strongest form, and `_grow` already merges nearby targets into one FOV rather than
   emitting two overlapping ones. Relaxing it would add redundant data *and* move the
   single-field goldens, destroying the only regression control we have. `sample_size`
   therefore remains a per-well **maximum**, and shortfall keeps its explicit reason
   (`low_cells` vs `constrained`).

5. **Multi-field without acquisition metadata is a hard failure.** It must never silently
   fall back to treating fields as independent.

## Why this is surgical

`gate()`, `eligible()`, `_surs_sample()`, `_grow()`, `_place_disjoint()`,
`select_and_place()`, `fov_row()`, `_write_csv()` and `write_curated_output()` are
**unchanged**. They already do the right thing - they are simply handed a well instead of
a field, in one coordinate frame.

The single structural decision that keeps the blast radius small: **FOV tiles stay
attached to a site and stay in that site's local coordinates.** Well-level pooling is
internal to `curate()`. Consequences:

- `_render_plate_overview` keeps working untouched - it already adds the site origin
  before drawing.
- `write_curated_output` keeps writing one CSV per site untouched.
- `fov_row` keeps emitting local coordinates, which is what MetaXpress needs (the
  confirmation acquisition proved it applies the field offset itself; median 0.69 um,
  max 1.55 um over 99 positions).

## New pure functions (section 4, beside the existing geometry)

All portable CPython 2/3 + Jython 2.7, no I/O, individually testable.

| Function | Contract |
|---|---|
| `well_key(site)` | `(row, col, z, t)` - the physical well and facet, from `parse_site` |
| `to_well(box, origin)` | field-local box -> well-coordinate box |
| `to_local(x, y, origin)` | well point -> field-local point (the output leg) |
| `owning_field(x, y, candidates, origins, field_wh)` | nearest field centre among the fields containing the point; ties -> lowest field number |
| `partition_by_owner(boxes_by_field, origins, field_wh)` | per field, the boxes it owns, plus the pooled well-coordinate population |
| `fov_home_field(cx, cy, origins, field_wh)` | the field whose image most safely contains a placed FOV; ties -> lowest field number |

## Changed functions

**`curate()`** - the site loop becomes two levels: group `files_by_site` by `well_key`,
read and gate each field as today, translate to well coordinates, partition, pool, call
`select_and_place` **once**, then send each FOV home to a field and convert back to local.

Signature gains two injected, optional parameters:

```python
def curate(results_dir, gate_spec, fov_px=256.0, sample_size=5, seed=42,
           neighbourhood=2.0, stage_margin=0.05, out_csv=None, base=None, progress=None,
           origins=None, field_wh=None):
```

Injected rather than read inside `curate()` so the function stays testable without files,
and so the single-field path is provably identical when they are absent. If a well has
more than one field and `origins` is missing, raise (decision 5).

**Two details that protect the goldens and must be implemented deliberately:**

- **Seed enumeration.** Today the seed is `seed + well_index` over `sorted(files_by_site)`
  - a **string** sort. Wells must be enumerated in the order of their first site in that
  same sorted sequence, not by re-sorting parsed keys, so the single-field index sequence
  is bit-identical. (String and numeric order coincide for the current plates only because
  rows are single-digit and columns are zero-padded - do not rely on that silently.)
- **Sampling extent.** `_surs_sample` spans the extent of `search_boxes` (all base
  objects), deliberately - you do not sample where nothing was imaged. Keep that
  definition and compute it over the **pooled** well population. Switching to the measured
  union of field rectangles would be defensible but is a different extent and would move
  the montage goldens.

**`write_curated_output()`** - one string only: the changes-log line
`Cells per site/field` becomes `Cells per well`.

**`_run_curation()` / `run_macro()`** - pass `origins` and `field_wh` through from
`read_field_origins()` and `overview_scale()`.

## Test plan

Written and green **before** any behaviour change. Tests live in section 12 and run in
both CPython and Fiji Jython.

### Phase 0 - build the net (no production code touched)

The real-data goldens (`4253` positives, acquired `165`, FOVs `164`, captured `313`,
extra `148`) currently sit behind `os.path.dirname(__file__)/Results/...`, which does not
exist in the repo, so **the strongest regression test silently skips**. The overview tests
in the same file already reach the dataset by absolute path. Point the golden block at the
same root so the control actually runs.

Then add, against unmodified code:

- Pin the synthetic single-field dataset's tile **coordinates**, not only its counts.
- A synthetic **four-field well** fixture: 2304 px fields, 2073.5 px pitch, cells planted
  in a field interior, in a two-way seam, and in the four-way corner square. Assert the
  present (wrong) numbers, so the defect is documented as a test and its correction is
  visible as a deliberate flip.

Gate: all green, goldens run rather than skip.

### Phase 1 - pure geometry, not yet wired

- Partition is exact: retained regions pairwise disjoint **and** their union equals the
  union of the fields - no gap, no double cover. Cases: 1x1, 2x2 at 10%, 3x3, unequal
  overlap per axis, abutting tiles, a gap between tiles.
- Every object has exactly one owner (totality).
- Midline equivalence: on a regular grid, nearest-centre ownership agrees cell-for-cell
  with a 5% edge trim.
- Tie-break: a centroid exactly on the midline gets one owner, the lower field number -
  never zero, never two.
- Jitter guard: two detections of one physical cell straddling the midline by less than
  the guard distance collapse to a single candidate.
- Round trip: well -> local -> well reproduces the position exactly; and against the
  confirmation acquisition, local + measured field origin reproduces the recorded stage
  position within the observed 1.55 um.
- Degenerate: one field per well trims nothing and returns the input unchanged.
- Multi-field with `origins=None` raises.

Gate: all green, all Phase 0 goldens still identical (this code is not yet called).

### Phase 2 - wire it into `curate()`

- **Nico is byte-identical.** Not just the five golden numbers - the written CSVs compare
  byte-for-byte against the pre-change run. This is the primary regression gate.
- The four-field fixture's assertions flip to the corrected values:
  - duplicate cells in the seam and the corner collapse to one candidate each;
  - the acquired count is a per-well maximum, never `4 x n`;
  - no two FOVs of one well overlap **in well coordinates**;
  - one SURS grid per well - the same random offset governs all four fields.
- Per-area uniformity holds at well level: the existing clump test, repeated with the
  clump straddling a seam.
- `captured` / `extra` / `status` are well-level and self-consistent
  (`captured == acquired + extra`).

### Phase 3 - real data

- Babette re-baselined and the new goldens pinned. Expect the acquisition count to
  **drop** - four fields each asking for n becomes one well asking for n, and seam
  duplicates collapse. That is the fix working, not a regression; none of the current
  Babette numbers are comparable across the change.
- R5-C1 specifically: zero cross-field FOV overlaps, no duplicate pair within the guard
  distance.
- Plate overview re-rendered and inspected - fields still drawn overlapping (real
  geometry), FOVs no longer duplicated across the seam.
- Full suite green in CPython **and** Fiji Jython 2.7.4; macro path green in the Fiji GUI.

## Out of scope

- The `GenericDialog` and SciJava widgets (GUI-only, untestable headless).
- MetaXpress ingestion of authored rows, and the stale measurement columns those rows
  carry.
- Any change to gating, eligibility, `fov_row`, or the output folder convention.
