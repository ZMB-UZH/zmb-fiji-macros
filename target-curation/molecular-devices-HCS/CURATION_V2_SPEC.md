# MD HCS target curation — v2 redesign spec

**Status: IMPLEMENTED + tested (Jython 2.7.4 in Fiji and CPython).** Flow: gate
positives -> Systematic Uniform Random Sampling (SURS, portable seeded PRNG) ->
fewest FOVs acquiring the sample. Files in this folder: `curation_core.py` (logic), `curation_pipeline.py`
(per-well pipeline), `ZMB_MD_HCS_target_curation_v2.py` (Fiji entry: params +
curated CSV + visual report), `test_curation_core.py` (full suite; runs in CPython
and Fiji Jython, results byte-identical across both).
The algorithm was reviewed over two adversarial subagent rounds (one critical
placement bug found + fixed; feasibility-interval method then proved correct
under 600k fuzz trials). Superseded the v1.9.0 class-dropdown edits. All
quantitative claims below are measured against the real dataset
`NB26-15_Overview10x_DAPI-mScarlet` (48 wells, 2 signals: DAPI/nuclei = T1,
mScarlet cells = T2); the pipeline reproduces 4253 positives, 73 FOVs, 16 short
wells and these are pinned as regression tests.

The optional spread knobs (max FOV overlap %, max cells per FOV, min FOVs per
well) are IMPLEMENTED and exposed in the entry, each disableable via a 0/off value
(defaults reproduce fewest-disjoint-FOVs exactly). Verified end-to-end: off -> 73
FOVs, on (overlap 20% / cap 3 / min-FOV 2) -> 78 FOVs.

**Still to do (bench / GUI):** the adaptive GenericDialog (a checkbox per
discovered signal) is GUI-only and untested — the entry currently takes a text
gate string via SciJava params (works headless + auto-dialog in GUI). Plus the
MetaXpress verify items below (bbox centre-vs-corner, authored-row import,
changer/ROI).

## Purpose

Curate IN Carta HCS output so a downstream high-magnification re-acquisition
images a small, well-chosen set of cells per well instead of everything detected.
v1 thinned the per-cell target list at random with a min-distance rule. v2
reframes it as: **select the positive cells per well, then place the fewest
imaging tiles (FOVs) that capture at least N of them.**

## Guiding principle

**Everything is based on segmented objects.** The only primitive is an object's
bounding box (montage pixels), present in every TargetData CSV. No mask TIFFs, no
pixel data, no IN Carta linkage column, no coordinate conversion. Gating,
counting, FOV placement and spacing are all box operations.

## Workflow context (decided)

The curated file is **loaded into MetaXpress by hand** at re-acquisition setup.
Consequences that shape the whole design:
- No auto-read timing window to miss (the original v1 README caveat is void).
- The macro **does not overwrite `TargetData/`**. It writes a NEW curated file;
  IN Carta originals are never renamed or touched. This removes v1's entire
  rename / in-place-overwrite / assert-previous-curation machinery and its risk
  of destroying fresh IN Carta output.

---

## Dialog — three areas

### Area 1 — Signals (scales with the channel count)

Signals are discovered from the data, never hardcoded:
- signal **name** = text before the first `_` in the filename
- signal **T-index** = the single-`T` prefix (`T<n>$`) on its own columns

For this dataset: `T1 = DAPI`, `T2 = mScarlet cells`. A 3-signal analysis renders
3 rows, 4 signals renders 4 — automatically.

Model is **top-down**: the **first channel is the base object / main
compartment** — the thing you image and centre FOVs on (DAPI nuclei). No separate
base picker; base = channel 1 by position. Every other channel is a gate.

| channel | selector | meaning |
|---------|----------|---------|
| 1 (base) | — | the object you image; the geometry source and the anchor |
| 2..N | Yes | base object touches >=1 object of this channel |
| 2..N | No  | base object touches 0 objects of this channel (negative gate) |
| 2..N | left out | channel not used in the gate |

**Gate is evaluated per base object (within ONE compartment).** A base object is
selected iff, for that same object, every Yes channel touches it and every No
channel does not. Triple positive = the *same* nucleus touches a Ch2 AND a Ch3
object — not merely both signals present somewhere in the well. Number of Yes
gates -> single / double / triple / quadruple positive; a No gate gives a matched
**negative control** set (same objective + geometry, directly comparable).

**Overlap = any-touch (`>0`)** — a channel object counts if its box touches the
base box at all, not required to be inside. Measured vs IN Carta's own linkage:

| rule | positives | agreement | false+ |
|------|-----------|-----------|--------|
| **any-touch (>0)** (chosen) | 4253 | 99.77% | 705 |
| centre-in | 3614 | 99.98% | 66 |
| >=50% of box | 3594 | 99.98% | 46 |
| full-contain | 3559 | 100.00% | 11 |

Any-touch adds ~20% positives vs the instrument (diagonal-neighbour boxes);
accepted as simplest and symmetric. Full-containment reproduces IN Carta exactly
if ever wanted.

### Area 2 — Objective for target acquisition

Dropdown showing **Mag / Name / NA / Immersion** from the lab optics list;
resolves the FOV footprint (a rectangle in montage px).

Read from files, not typed:
- overview pixel size: `experiment(_montage)/*.jdce` ->
  `ObjectiveCalibration.PixelWidth` = **0.6473 um/px**
- ROI (2304x2304), binning (1), overview objective/changer (10x / 1x):
  `*.mxprotocol`

`FOV_montage_px = ROI_px * umPerPx(target) / umPerPx(overview)`. Camera pixel
size **cancels** (same camera + binning), so it is never needed.

TO VERIFY (target acquisition, from a protocol not yet seen):
- **changer** 1x vs 1.5x — a 50% FOV swing (256 px vs 384 px). Runtime control,
  default 1.5x.
- **ROI** — assumed full-frame 2304x2304; may be cropped.

### Area 3 — Target selection

**ENFORCED (hard, always applied):**

| parameter | unit | default |
|-----------|------|---------|
| Min cells | count (coverage target; report short if unreachable) | — |
| Neighbourhood ratio | x object size, clear space each side (min context to count) | 2x |
| Stage-accuracy margin | % of FOV per side | 5%/side (10% total) -> inner 90% usable |

**OPTIONAL (off by default, never forced; only tighten further when on):**

| parameter | unit |
|-----------|------|
| Min FOV | min number of FOVs per well (force spatial spread) |
| Max cells per FOV | count cap per tile (force spatial spread) |

Two separate margins, both enforced together:
- **Stage-accuracy** (% of FOV/side): buffer for stage positioning error; shrinks
  usable FOV to inner 90%. (Physics note: stage error is a fixed distance, so as
  %-of-FOV it scales with objective — fine for a fixed target objective, revisit
  or use um if magnifications vary widely.)
- **Neighbourhood** (x object size): biological context around the cell.

**Hard capture test:** a cell counts iff its box + neighbourhood fits within the
FOV shrunk by the stage-accuracy margin.

---

## Selection = SURS, THEN minimal FOVs (the important flow)

Per well: (1) gate the positives, (2) select `sample_size` of the validly-
observable ones by **Systematic Uniform Random Sampling (SURS)** - the stereology
standard, (3) place the **fewest FOVs that acquire that sample** efficiently. All
quantities are per well.

- **SURS**: lay a grid sized for n points over the well's extent with ONE
  uniform-random offset (the random start), and take the nearest not-yet-chosen
  eligible cell at each grid point. Even spatial coverage (systematic) but
  unbiased because the grid position is randomised. Grid = ncol x nrow with
  ncol = round(sqrt(n*W/H)), nrow = ceil(n/ncol); exactly n cells when eligible
  > n. Measured: for n=4 SURS lands one cell in each quadrant ~100% of the time
  vs ~9% for pure random - the systematic even-coverage property.
- A cell is **eligible** only if its box + neighbourhood margin fits a FOV, so a
  sampled cell is never clipped by the field of view (the margin's purpose).
- Uses a portable SplitMix64 PRNG (integer core, 64-bit masking, top-53-bit
  mantissa) because Python's own `random` is NOT reproducible across CPython and
  Jython (verified: `random()` diverges). Seed + well-index; the whole pipeline is
  byte-identical in CPython 3 and Jython 2.7.4. NOTE round() differs across Py2/Py3
  on .5 - use int(x+0.5).
- **Placement** over the fixed sample: fewest FOVs, each centred on the feasibility-
  interval midpoint of its group so every member is provably inside. FOV count is
  minimised, so overlap only occurs where two sampled cells are geometrically
  forced to share ground (`max_overlap` defaults to 1.0 = image every sampled
  observation; lower it toward 0 to prefer strict-disjoint at the cost of dropping
  close sampled cells). This is the "minimal overlap" the operator asked for:
  minimise FOVs -> minimal overlap, but never drop an observation.
- Cells too large to be validly observed are excluded before sampling; a well with
  fewer than `sample_size` eligible positives is `sample_short` (reported).

Real plate (sample_size=5, seed=42, nbhd=2x, margin=5%): 167 sampled, 165 FOVs,
16 under-sampled wells, 0 placement-short — pinned as regression tests, identical
in CPython and Jython.

## Output — generate FOV rows

Because the macro authors the file and it is hand-loaded, output is not limited
to detected objects. **Write one generated row per FOV**, bbox centred on the
group centroid. The file is a **list of FOVs to acquire**, not a thinned cell
list; covered cells are audit metadata.

Generate a valid row ("conform to the format so MetaXpress can't tell"):
template from a real row in the SAME well, keep every identity column correct
(well/row/column/field/z/t), set the bbox to the centroid, fresh object_id.

TO VERIFY on the scope:
- MetaXpress manual import accepts authored rows (near-certain for a coordinate
  import). Fallback if not: anchor each FOV on the real object nearest the
  centroid — measured cost only +2.7% tiles (74 -> 76), identical tiles/well and
  identical 16 short wells.
- Does MetaXpress use the bbox CENTRE or top-left CORNER as the acquisition
  position? Set the box accordingly (`X = cx - W/2` if centre). Wrong choice
  offsets every tile by half a box — visible immediately; one test settles it.
- Exact format + location the manual import expects.

## Output — visual report (first-class)

Feeds an EXPENSIVE re-acquisition, so the operator must SEE the selection before
committing. Also replaces a heavy interactive "Calculate -> preview" UI: a modal
run + report is enough.
- **Per well:** overlay on the actual montage TIFF — all nuclei faint, gated
  positives coloured (imaged vs not-imaged), FOV rectangles + stage-safe inner
  box + tile centre.
- **Plate summary:** 8x6 well grid heat-mapped by tiles placed; the 16 short
  wells flagged.
- Makes the clustering bias visible (dense well, 301 positives, Min=5 -> one
  256 px FOV images 5, 296 ignored) so the operator can raise Min or enable the
  spread knobs.

Mock-up rendered from real data: `visual_report_mock.png` (this folder).

---

## Data findings that constrain the design

1. **Wells are bimodal.** 12/48 wells have 0 positives; 4 have 1-2; the rest have
   >=113. So Min=5 and Min=10 fail on the identical 16/48 wells. The "cannot reach
   Min" report is essential — a third of the plate always comes back short
   (biology: positives concentrate in row D).
2. **Clumping bias.** One FOV captures >1 positive 32% of the time, up to 8.
   Minimising FOVs preferentially grabs clumps.
3. **Neighbourhood 2x default has a cost.** At 60x+1.5x (256 px FOV) only 84.9% of
   positives fit; the excluded ~15% are the biggest boxes = merged multi-nucleus
   blobs (arguably a good filter, but report it). At 60x+1x (384 px) 97.5% fit —
   so the changer changes what the default does.
4. **"Fully captured" never binds here.** Positive boxes max 122 px vs 256 px FOV.
5. **Coordinates are montage pixels**, and the montage is a genuine registration
   (per-well size varies 8288-8525 px), so nominal well-centre + offset would be
   off by ~65 um. The macro stays entirely in montage-pixel space; MetaXpress
   owns the montage->stage mapping.
6. **Segmentation caveat.** An "object" is only as good as IN Carta's
   segmentation, which occasionally fuses several nuclei into one box.

## Decided
- Object/box-based throughout; overlap = any-touch (>0).
- Channel 1 = base/main compartment (top-down); gates are Yes/No/left-out,
  evaluated per compartment; Yes-count -> double/triple/quadruple positive.
- Objective -> FOV; overview scale read from `.jdce`/`.mxprotocol`.
- Minimise FOVs per well s.t. >= Min cells; disjoint by default (max-overlap-%
  tolerance = 0); FOV centred on group centroid; output = generated FOV rows.
- Enforced: Min cells, neighbourhood ratio (2x), stage margin (5%/side).
  Optional: Min FOV, Max cells per FOV.
- Manual load; no overwrite of IN Carta output.
- minDistPx/seed/shuffle removed (340 default was wrong: assumed a 2048 sensor;
  this camera is 2304, so 340 actually permitted overlapping tiles).

## Open / to verify
1. Target changer (1x/1.5x) and ROI — from the target-acquisition protocol (not
   yet seen).
2. MetaXpress manual import: accepts authored rows? bbox centre vs corner? exact
   format + location?
3. "Min FOV" — confirm it means minimum FOV COUNT per well (force spread), not
   spacing/size.
4. Greedy tiebreak rule (define one, e.g. lowest coordinate, for determinism).
5. UI: simple (OK = calculate + visual report) vs interactive preview panel.
   Recommended: simple for v1.

## Validation done (Python ports vs real data, in scratchpad)
- Signal discovery -> exactly [DAPI, mScarlet cells]; selecting DAPI keeps 48
  files. 13/13 logic checks pass.
- Overlap-rule comparison table (above).
- Greedy: 48 wells in 18.6 ms; median 2 tiles/well, max 3; 74 total tiles at
  Min=5; 16 wells short. Centroid vs anchor: 74 vs 76.
- Visual-report mock-up from real positives.
NOT validated: the IJM itself (no Fiji on the dev machine) — dialogs, SciJava
params, protocol reads, and the ImageJ overlay rendering are unexercised. A real
Fiji run on the sandbox copy is required before live use.
