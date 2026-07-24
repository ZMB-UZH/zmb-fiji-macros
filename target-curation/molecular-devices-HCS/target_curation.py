#@ File (label="InCarta analysis results folder", style="directory") results_dir

"""
MD HCS target curation - single-file Fiji macro (Jython 2.7) and CPython module.

Molecular Devices ImageXpress / IN Carta writes one CSV of segmented objects per
signal per acquisition site. This picks, per site/field, a small representative set of cells to
re-image at high magnification, in three steps:

  1. GATE the objects of interest. One class is the base object (the thing gated and
     re-imaged, e.g. the nuclei); it is kept when it is single / double / triple
     positive - i.e. its box overlaps an object of every 'yes' signal and none of the
     'no' signals. In Fiji every class found in the data gets its own dropdown: mark
     one class as the object (base) and each of the rest positive / negative / ignore
     (base defaults to the first channel).

  2. SAMPLE ~`sample_size` of them by SURS (Systematic Uniform Random Sampling, the
     2D stereological grid): lay a grid of ~n frames over the scanned area with ONE
     random offset and take one random cell from each occupied frame. Sampling is
     uniform over AREA, so a dense cluster is sampled by the area it covers, not its
     cell count, and cannot dominate the sample. A cell is eligible only if the FOV
     can give it a clear radius of `neighbourhood` x its enclosing circle, so it is
     never clipped by the field of view.

  3. PLACE DISJOINT FOVs over the sample. Each FOV frames its cell and any further
     sampled cells that also fit the one field; FOVs never overlap, so no ground is
     imaged twice. A sampled cell too close to an already-placed FOV to be framed
     without overlap is dropped (rare, since the sample is spread).

OUTPUT follows the v1 convention, written in place next to the analysis: the first
run preserves IN Carta's `TargetData/` as `TargetData_original/`; every run rewrites
`TargetData/` from scratch with the curated per-site CSVs (the generated FOV-centre
rows MetaXpress re-images) and mirrors them into `TargetData_curated/` with a
`curation_changes.csv` audit log (positives / eligible / acquired / extra / fovs per
site), a per-site overlay PNG under `report/`, and whole-plate overview PNGs.
A rerun refuses to overwrite `TargetData/` unless it still matches the last curated
mirror, so newly regenerated IN Carta output is never clobbered. `acquired` = sampled target cells (imaged whole);
`extra` = bonus positives that also fall entirely inside a field.

In Fiji you point at the InCarta results folder, then a dialog assembles the run in
sections: one dropdown per marker class found in the data to set positives (positive /
negative / ignore); the target objective; and the parameters (cells per site/field,
neighbourhood as % of the cell radius, stage margin, seed). The FOV size is inferred
from the overview image's own metadata (objective, changer, binning, sensor region);
the overview image is located automatically from the InCarta result metadata
(`result_metadata.csv` names it, reachable by walking up to the acquisition), so there
is normally no image input, and the size adapts when the overview magnification
changes. If the image cannot be found (e.g. the results folder was copied away from
the acquisition), the dialog asks for it.

The file is BOTH the Fiji entry (SciJava params above) and its own test suite: run
`python target_curation.py` in CPython to run the tests; open/run it in Fiji to
curate. All maths is stdlib-only and written to be interpreter-independent (a
portable PRNG and int(x+0.5) rounding, because Python's own `random` and `round`
differ between CPython 3 and Jython 2.7), so results are byte-identical in both.

The code below reads top-to-bottom as the pipeline flows: import the data (1) ->
gate it (2) -> work out the FOV size (3) and what fits in a FOV (4) -> draw random
numbers (5) -> sample (6) -> place fields (7) -> combine per site (8) -> write the
output (9) -> curate every site (10) -> the Fiji front end (11), tests (12), dispatch.
"""
from __future__ import division, print_function
import io, os, re

_BB = "AS_FID_Blob_BoundingBox"          # + X / Y / Width / Height
_OWN = re.compile(r"^T(\d+)\$")          # a signal's own-measurement column prefix
_SITE_KEY = "_singleTargetData_"         # filename split: <signal>_singleTargetData_<site>.csv
_SITE_RE = re.compile(r"^R(\d+)-C(\d+)-F(\d+)-Z(\d+)-T(\d+)$")
_MASK64 = (1 << 64) - 1


# --------------------------------------------------------------------------- #
# 1. Import the data - find the analysis CSVs and format them as bounding boxes#
# --------------------------------------------------------------------------- #
def read_csv(path):
    """(header, rows-as-dicts). Tolerates a UTF-8 BOM; drops ragged rows."""
    import csv
    with io.open(path, "r", encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.reader(fh))
    if not rows:
        return [], []
    head = rows[0]
    return head, [dict(zip(head, r)) for r in rows[1:] if len(r) >= len(head)]


def signal_name(fn):
    """A signal's name is the filename text before the first underscore."""
    return fn.split("_")[0]


def _site_of(fn):
    """The well/site token a CSV belongs to: the text after `_singleTargetData_`."""
    stem = fn[:-4]
    return stem.split(_SITE_KEY, 1)[1] if _SITE_KEY in stem else stem


def parse_site(site):
    """Parse `R#-C#-F#-Z#-T#` into integer row/column/field/z/time values.

    Rendering must never silently collapse an unrecognised identifier onto another
    panel, so malformed site tokens fail loudly instead of being partially parsed.
    """
    m = _SITE_RE.match(site)
    if not m:
        raise ValueError("Malformed site identifier: " + site)
    return tuple(int(v) for v in m.groups())


def site_status(eligible_count, acquired_count, requested_count):
    """Classify why a site did or did not reach the requested acquisition count."""
    if acquired_count >= requested_count:
        return "full"
    if eligible_count < requested_count:
        return "low_cells"
    return "constrained"


def _field_grid(site_keys):
    """Map fields to a deterministic, collision-free near-square report grid."""
    keys = sorted(site_keys, key=lambda k: k[2])
    ncols = max(1, int(len(keys) ** 0.5 + 0.999999))
    nrows = max(1, (len(keys) + ncols - 1) // ncols)
    slots = dict((k, (i // ncols, i % ncols)) for i, k in enumerate(keys))
    return slots, nrows, ncols


def build_overview_layout(wells):
    """Pure layout model used by the Fiji renderer and CPython regression tests.

    Returns {(z, time): {site: layout}}, where layout contains the physical plate
    row/column and a unique field subpanel slot. Z/time variants are separate facets.
    """
    facets = {}
    for well in wells:
        key = parse_site(well["site"])
        facets.setdefault((key[3], key[4]), {}).setdefault((key[0], key[1]), []).append(key)

    result = {}
    for facet, physical_wells in facets.items():
        rows = sorted(set(k[0] for k in physical_wells))
        cols = sorted(set(k[1] for k in physical_wells))
        canonical = dict((key[2], key)
                         for keys in physical_wells.values() for key in keys)
        field_slots, subrows, subcols = _field_grid(canonical.values())
        field_slots = dict((key[2], slot) for key, slot in field_slots.items())
        laid_out = {}
        for (row, col), keys in physical_wells.items():
            for key in keys:
                laid_out["R%d-C%d-F%d-Z%d-T%d" % key] = {
                    "row": row, "col": col,
                    "well_row": rows.index(row), "well_col": cols.index(col),
                    "field_row": field_slots[key[2]][0],
                    "field_col": field_slots[key[2]][1],
                    "field_rows": subrows, "field_cols": subcols,
                }
        result[facet] = laid_out
    return result


def dataset_extent(wells):
    """Shared coordinate extent so every report panel uses the same visual scale."""
    boxes_ = [box for well in wells for box in well["base"]]
    return (max([b[0] + b[2] for b in boxes_] or [1.0]),
            max([b[1] + b[3] for b in boxes_] or [1.0]))


def _has_csv(d):
    return os.path.isdir(d) and any(f.lower().endswith(".csv") for f in os.listdir(d))


def find_target_dir(path):
    """Locate the TargetData folder: `path` if it is one (or the TargetData_original
    that a previous run renamed the source to), else a TargetData child, else one
    analysis-folder level down. Preferred over sibling ObjectData/metadata CSVs, which
    also live in the analysis folder but are not the targets."""
    if os.path.basename(os.path.normpath(path)) in ("TargetData", "TargetData_original") and _has_csv(path):
        return path
    cand = os.path.join(path, "TargetData")
    if _has_csv(cand):
        return cand
    for name in sorted(os.listdir(path)):
        cand = os.path.join(path, name, "TargetData")
        if _has_csv(cand):
            return cand
    raise ValueError("No TargetData folder with CSVs under: " + path)


def discover_signals(target_dir):
    """{name: T-index} + names ordered by (T-index, name) so the base is first. A
    signal's T-index is the T<n>$ prefix its own measurement columns carry - read
    from the data rather than hardcoded, so the tool is not tied to these signals."""
    idx = {}
    for fn in sorted(os.listdir(target_dir)):
        if not fn.lower().endswith(".csv"):
            continue
        name = signal_name(fn)
        if name in idx:
            continue
        head, _ = read_csv(os.path.join(target_dir, fn))
        for c in head:
            m = _OWN.match(c)
            if m:
                idx[name] = int(m.group(1))
                break
    return idx, sorted(idx, key=lambda n: (idx[n], n))


def _bb_cols(t):
    return ["T%d$%s%s" % (t, _BB, s) for s in ("X", "Y", "Width", "Height")]


def boxes(rows, t):
    """[(x, y, w, h)] from target T<t>'s own bounding-box columns; skip rows whose
    bbox is empty/non-numeric (placeholder rows)."""
    cols = _bb_cols(t)
    out = []
    for r in rows:
        try:
            out.append(tuple(float(r[c]) for c in cols))
        except (KeyError, ValueError):
            pass
    return out


# --------------------------------------------------------------------------- #
# 2. Gate - keep the objects of interest (box overlap + the gate spec)        #
# --------------------------------------------------------------------------- #
def touches(a, b):
    """True if two boxes overlap with strictly positive area (edge-only contact,
    zero area, does not count as overlap)."""
    return not (a[0] + a[2] <= b[0] or b[0] + b[2] <= a[0] or
                a[1] + a[3] <= b[1] or b[1] + b[3] <= a[1])


def gate(base, gates, require):
    """Keep base boxes where every 'yes' signal touches and every 'no' does not.
    gates = {signal: [boxes]}; require = {signal: True(yes)/False(no)}; a signal
    absent from `require` is ignored (an empty require keeps all base objects)."""
    return [bb for bb in base
            if all(any(touches(bb, g) for g in gates.get(s, [])) == want
                   for s, want in require.items())]


def parse_gate(spec):
    """'mScarlet cells:yes; foo:no' -> {'mScarlet cells': True, 'foo': False}."""
    req = {}
    for part in spec.split(";"):
        sig, sep, want = part.rpartition(":")
        if sep and sig.strip():
            req[sig.strip()] = want.strip().lower() in ("yes", "y", "true", "1", "+")
    return req


def build_gate_spec(markers, roles):
    """Assemble a gate string from per-marker roles chosen in the dialog. `roles` maps
    a marker name to 'positive' / 'negative' / 'ignore'; positives become '<m>:yes',
    negatives '<m>:no', and ignored markers are left out. Kept in `markers` order so
    the spec reads in channel order. The inverse of parse_gate."""
    parts = []
    for m in markers:
        role = roles.get(m, "ignore")
        if role == "positive":
            parts.append("%s:yes" % m)
        elif role == "negative":
            parts.append("%s:no" % m)
    return "; ".join(parts)


# --------------------------------------------------------------------------- #
# 3. FOV size - inferred, per target objective, from the overview image        #
# The target camera always reads the same pixel count; that field's footprint  #
# in the overview montage scales with the overview-to-target magnification.     #
# --------------------------------------------------------------------------- #
# target objective label -> (magnification, changer); total magnification = mag * changer
_OBJECTIVES = {"2x": (2.0, 1.0), "4x": (4.0, 1.0), "10x": (10.0, 1.0),
               "20x": (20.0, 1.0), "40x": (40.0, 1.0), "60x": (60.0, 1.0),
               "60x + 1.5x changer": (60.0, 1.5)}
_OBJECTIVE_ORDER = ["2x", "4x", "10x", "20x", "40x", "60x", "60x + 1.5x changer"]  # dialog order


def read_image_description(path):
    """The TIFF ImageDescription (tag 270) as text - the MetaXpress MetaData XML. A
    minimal pure-Python TIFF read (header + first IFD only), so it runs in both
    CPython and Jython without an imaging library and touches only a few bytes of
    what may be a very large montage."""
    import struct
    with open(path, "rb") as fh:
        head = fh.read(8)
        order = "<" if head[:2] == b"II" else ">"
        if struct.unpack(order + "H", head[2:4])[0] != 42:      # classic TIFF only
            return ""
        fh.seek(struct.unpack(order + "I", head[4:8])[0])       # first IFD
        for _ in range(struct.unpack(order + "H", fh.read(2))[0]):
            entry = fh.read(12)
            tag, _typ, count = struct.unpack(order + "HHI", entry[:8])
            if tag == 270:                                      # ImageDescription
                if count <= 4:
                    raw = entry[8:8 + count]
                else:
                    here = fh.tell()
                    fh.seek(struct.unpack(order + "I", entry[8:12])[0])
                    raw = fh.read(count)
                    fh.seek(here)
                return raw.rstrip(b"\x00").decode("latin-1")
    return ""


def overview_scale(description):
    """(magnification, changer, binning, region_px) parsed from a MetaXpress
    ImageDescription: the overview objective ('10X ...' -> 10), tube-lens/changer
    ('1X' -> 1.0, '1.5X' -> 1.5), camera binning, and the sensor region in pixels.
    Reading these from the image is what lets the FOV adapt to the overview's own
    objective rather than assuming one."""
    def prop(name, default):
        m = re.search(r'"%s"[^>]*value="([^"]*)"' % re.escape(name), description)
        return m.group(1) if m else default

    def lead(text, default):
        m = re.match(r"\s*([0-9.]+)", text)
        return float(m.group(1)) if m else default

    mag = lead(prop("HCS.ai Objective", ""), 0.0)
    changer = lead(prop("HCS.ai Tube Lens", "1X"), 1.0)
    binning = int(lead(prop("camera-binning-x", "1"), 1.0))
    region = re.search(r"Region:\s*([0-9]+)\s*x", description)
    return mag, changer, binning, int(region.group(1)) if region else 2304


def overview_pixel_size(description):
    """Overview pixel size in micrometres, read from the same image metadata."""
    m = re.search(r'"spatial-calibration-x"[^>]*value="([^"]*)"', description)
    return float(m.group(1)) if m else None


def fov_px_for(target_mag, target_changer, mag_ov, changer_ov, binning_ov, region_px):
    """FOV footprint in overview-montage pixels for a target objective. The camera
    reads `region_px` pixels at any magnification; its montage footprint is that
    count scaled by the overview-to-target total-magnification ratio (overview
    binning coarsens the montage, so it enlarges the footprint)."""
    return region_px * (mag_ov * changer_ov) / (target_mag * target_changer * binning_ov)


def _overview_tiff_path(results_dir):
    """Path to the overview image, found from the InCarta result metadata (so the
    operator needn't select it). `result_metadata.csv` / `channel_metadata.csv` in the
    analysis folder carry a `greyscale_image` name like `timepoint0\\...tif`, relative
    to the acquisition root; we read it and walk up from the results folder until that
    file exists. Returns None if it cannot be located."""
    try:
        analysis_dir = os.path.dirname(os.path.normpath(find_target_dir(results_dir)))
    except ValueError:
        return None
    for meta in ("result_metadata.csv", "channel_metadata.csv"):
        path = os.path.join(analysis_dir, meta)
        if not os.path.isfile(path):
            continue
        head, rows = read_csv(path)
        if "greyscale_image" not in head:
            continue
        rel = ""
        for row in rows:
            rel = row.get("greyscale_image", "").strip()
            if rel:
                break
        if not rel:
            continue
        rel = rel.replace("\\", os.sep).replace("/", os.sep)
        d = analysis_dir
        for _ in range(6):                                # walk up to the acquisition root
            cand = os.path.join(d, rel)
            if os.path.isfile(cand):
                return cand
            d = os.path.dirname(d)
    return None


def overview_from_results(results_dir):
    """The overview image's MetaXpress ImageDescription (for overview_scale), located
    automatically from the InCarta result metadata - no separate overview-image input.
    Empty string if the referenced image cannot be found."""
    path = _overview_tiff_path(results_dir)
    return read_image_description(path) if path else ""


def _acquisition_metadata_path(results_dir):
    """Nearest ancestor's field metadata, which defines physical site origins."""
    d = os.path.dirname(os.path.normpath(find_target_dir(results_dir)))
    for _ in range(6):
        path = os.path.join(d, "image_metadata_1.csv")
        if os.path.isfile(path):
            return path
        parent = os.path.dirname(d)
        if parent == d:
            break
        d = parent
    return None


def field_origins(rows, pixel_um):
    """Map (row, col, field, z, time) to a well-local field origin in pixels."""
    centres = {}
    for row in rows:
        try:
            key = (int(row["Row"]), int(row["Column"]), int(row["Field"]),
                   int(row["ZIndex"]), int(row["Timepoint"]))
            centres.setdefault(key, (float(row["PositionXUm"]), float(row["PositionYUm"])))
        except (KeyError, ValueError):
            pass
    groups = {}
    for key, centre in centres.items():
        groups.setdefault((key[0], key[1], key[3], key[4]), []).append(centre)
    origins = {}
    for key, (x, y) in centres.items():
        group = groups[(key[0], key[1], key[3], key[4])]
        origins[key] = ((x - min(p[0] for p in group)) / pixel_um,
                        (y - min(p[1] for p in group)) / pixel_um)
    return origins


def read_field_origins(results_dir, pixel_um=None):
    """Physical field origins for rendering, or an empty mapping if unavailable."""
    path = _acquisition_metadata_path(results_dir)
    pixel_um = pixel_um or overview_pixel_size(overview_from_results(results_dir))
    if not path or not pixel_um:
        return {}
    return field_origins(read_csv(path)[1], pixel_um)


# --------------------------------------------------------------------------- #
# 4. FOV geometry - what fits inside one field, and which cells are eligible    #
# --------------------------------------------------------------------------- #
def _usable_half(fov, stage_margin):
    """Half-width of the FOV that can be relied on once the stage-accuracy margin is
    trimmed from each side - a cell must fall inside this to be sure it is captured
    despite stage-positioning error."""
    return fov * (1.0 - 2.0 * stage_margin) / 2.0


def cell_radius(b):
    """Radius of the smallest circle enclosing the cell's bounding box."""
    return (b[2] * b[2] + b[3] * b[3]) ** 0.5 / 2.0


def _fully_inside(box, cx, cy, half):
    """True if the whole `box` lies within a FOV's usable square (centre (cx, cy),
    half-width `half`), i.e. that cell is imaged whole - not cut off at the edge."""
    return (box[0] >= cx - half and box[0] + box[2] <= cx + half and
            box[1] >= cy - half and box[1] + box[3] <= cy + half)


def eligible(cells, fov, stage_margin, neighbourhood):
    """Indices of cells that can be validly observed: a disk of radius r*(1+
    neighbourhood) - the cell's enclosing circle (radius r) plus a clear margin of
    `neighbourhood` x r around it - fits inside the FOV shrunk by the stage-accuracy
    margin. So a valid cell always gets its context and is never clipped by the FOV."""
    half = _usable_half(fov, stage_margin)
    return [i for i, b in enumerate(cells)
            if cell_radius(b) * (1.0 + neighbourhood) <= half]


# --------------------------------------------------------------------------- #
# 5. Reproducible random numbers                                              #
# --------------------------------------------------------------------------- #
def _rng(seed):
    """SplitMix64 -> uniform floats in [0, 1). The integer core is masked to 64 bits
    and only the top 53 bits form the mantissa, which yields IDENTICAL sequences in
    CPython 3 and Jython 2.7 - Python's own `random` module does not agree across
    the two interpreters, which would break reproducibility of the sample."""
    state = [seed & _MASK64]

    def draw():
        state[0] = (state[0] + 0x9E3779B97F4A7C15) & _MASK64
        z = state[0]
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & _MASK64
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & _MASK64
        z ^= (z >> 31)
        return (z >> 11) / 9007199254740992.0        # / 2**53, exact

    return draw


# --------------------------------------------------------------------------- #
# 6. Sample - Systematic Uniform Random Sampling (2D grid, per-area uniform)   #
# --------------------------------------------------------------------------- #
def _surs_sample(elig, centre_x, centre_y, n, seed, search_boxes):
    """At most n eligible cell indices by Systematic Uniform Random Sampling - the 2D
    stereological grid (per-AREA uniform).

    Lay a grid of ~n square frames over the SCANNED AREA (the extent of
    `search_boxes` - all nuclei, i.e. the region the overview imaged, not just where
    the positives landed), shift it by ONE uniform-random (x, y) offset, and from each
    frame that holds an eligible cell take ONE cell at random. The single random
    offset is the random start that makes it unbiased; one-cell-per-frame is what
    spreads the sample evenly over space.

    Per-AREA uniform, NOT per-cell: every equal-area frame contributes at most one
    cell, so a dense cluster is sampled in proportion to the AREA it covers, not the
    number of cells in it, and cannot dominate the sample. The sample is thus an
    unbiased picture of the well AREA (even spatial coverage) rather than of the cell
    population (which would over-weight dense regions). Fewer than n come back only
    where the well is sparse - an empty frame contributes nothing. Deterministic for
    a seed."""
    if len(elig) <= n:
        return list(elig)

    xs = [b[0] + b[2] / 2.0 for b in search_boxes]
    ys = [b[1] + b[3] / 2.0 for b in search_boxes]
    x0, y0 = min(xs), min(ys)
    width = (max(xs) - x0) or 1.0
    height = (max(ys) - y0) or 1.0

    # ~n square frames shaped to the well's aspect ratio (n_cols * n_rows ~ n)
    n_cols = max(1, int((n * width / height) ** 0.5 + 0.5))    # int(x+0.5): portable rounding
    n_rows = max(1, int(n / n_cols + 0.5))
    step_x, step_y = width / n_cols, height / n_rows

    draw = _rng(seed)
    off_x, off_y = draw(), draw()                              # the one random start (offset, both axes)

    frames = {}                                                # (col, row) -> cell indices in that frame
    for i in elig:
        col = min(n_cols - 1, int((centre_x[i] - x0) / step_x + off_x))
        row = min(n_rows - 1, int((centre_y[i] - y0) / step_y + off_y))
        frames.setdefault((col, row), []).append(i)

    chosen = []                                                # one random cell per occupied frame
    for key in sorted(frames):
        members = frames[key]
        chosen.append(members[int(draw() * len(members))])

    # Rounding the rectangular grid dimensions can make n_cols * n_rows exceed n
    # (for example, n=5 on a roughly square well becomes a 2x3 grid). Keep the
    # spatially uniform frame sample, but randomly discard the surplus frames rather
    # than returning N+1 cells or truncating in coordinate order.
    for i in range(len(chosen) - 1, 0, -1):
        j = int(draw() * (i + 1))
        chosen[i], chosen[j] = chosen[j], chosen[i]
    if len(chosen) > n:
        chosen = chosen[:n]
    return chosen


# --------------------------------------------------------------------------- #
# 7. Place - the fewest disjoint FOVs that acquire the sample                  #
# --------------------------------------------------------------------------- #
def _grow(anchor, available, centre_x, centre_y, half_window):
    """Grow a group of cells that can share ONE field of view, starting from `anchor`.

    A FOV centred at (X, Y) contains a cell when X is within half_window of the
    cell's centre_x and Y within half_window of its centre_y (half_window already
    accounts for the cell's size and its clear margin). So a set of cells fits one
    FOV exactly when the intersection of their per-axis [centre - half_window,
    centre + half_window] intervals is non-empty. We add nearby cells (nearest
    first) as long as that intersection survives, and place the FOV at its midpoint
    - which therefore contains every member, so a cell is only ever listed in a FOV
    that truly covers it.

    Returns (group, fov_centre_x, fov_centre_y)."""
    x_lo, x_hi = centre_x[anchor] - half_window[anchor], centre_x[anchor] + half_window[anchor]
    y_lo, y_hi = centre_y[anchor] - half_window[anchor], centre_y[anchor] + half_window[anchor]
    group = [anchor]

    nearest_first = sorted(
        (a for a in available if a != anchor),
        key=lambda a: (centre_x[a] - centre_x[anchor]) ** 2 + (centre_y[a] - centre_y[anchor]) ** 2)

    for i in nearest_first:
        new_x_lo, new_x_hi = max(x_lo, centre_x[i] - half_window[i]), min(x_hi, centre_x[i] + half_window[i])
        new_y_lo, new_y_hi = max(y_lo, centre_y[i] - half_window[i]), min(y_hi, centre_y[i] + half_window[i])
        if new_x_lo <= new_x_hi and new_y_lo <= new_y_hi:
            x_lo, x_hi, y_lo, y_hi = new_x_lo, new_x_hi, new_y_lo, new_y_hi
            group.append(i)

    return group, (x_lo + x_hi) / 2.0, (y_lo + y_hi) / 2.0


def _place_disjoint(sampled, centre_x, centre_y, half_window, fov):
    """Place disjoint FOVs over the sampled cells by greedy first-fit.

    Walk the sampled cells; for each not-yet-covered one, build the FOV that frames
    it and any further sampled cells that also fit (several per FOV is fine and free).
    Reject the FOV if it would overlap one already placed, so no ground is imaged
    twice: two equal FOV squares overlap iff their centres are less than `fov` apart
    on BOTH axes. A cell that can only be framed by an overlapping FOV is left
    unimaged (rare, because the sample is spread out). Greedy first-fit, not a
    provably minimal cover, but the disjoint sample makes the two coincide in
    practice. Returns [ {cx, cy, covered: [cell_index, ...]} ]."""
    placed, covered = [], set()
    for cell in sampled:
        if cell in covered:
            continue

        remaining = [i for i in sampled if i not in covered]
        group, fov_x, fov_y = _grow(cell, remaining, centre_x, centre_y, half_window)

        overlaps_placed = any(abs(fov_x - t["cx"]) < fov and abs(fov_y - t["cy"]) < fov
                              for t in placed)
        if overlaps_placed:
            continue

        placed.append({"cx": fov_x, "cy": fov_y, "covered": sorted(group)})
        covered.update(group)

    return placed


# --------------------------------------------------------------------------- #
# 8. Select and place - one well's positive cells -> acquired cells + FOVs     #
# --------------------------------------------------------------------------- #
def select_and_place(cells, fov, stage_margin, neighbourhood, n, seed, search_boxes=None):
    """SURS-sample ~n eligible cells, then acquire them with disjoint FOVs.
    `search_boxes` is the scanned area the sampling grid spans (all nuclei); defaults
    to the `cells` themselves. Deterministic and interpreter-independent for a seed.
    Returns (acquired_cells, tiles, n_eligible)."""
    half = _usable_half(fov, stage_margin)
    elig = eligible(cells, fov, stage_margin, neighbourhood)
    if half <= 0 or not elig:
        return [], [], len(elig)
    centre_x = dict((i, cells[i][0] + cells[i][2] / 2.0) for i in elig)
    centre_y = dict((i, cells[i][1] + cells[i][3] / 2.0) for i in elig)

    sampled = _surs_sample(elig, centre_x, centre_y, n, seed, search_boxes or cells)

    # how far a FOV centre may sit from a cell's centre and still frame it + its margin
    half_window = dict((i, half - cell_radius(cells[i]) * (1.0 + neighbourhood)) for i in sampled)
    tiles = _place_disjoint(sampled, centre_x, centre_y, half_window, fov)

    acquired = sorted(set(i for t in tiles for i in t["covered"]))
    return [cells[i] for i in acquired], tiles, len(elig)


# --------------------------------------------------------------------------- #
# 9. Output - generated FOV rows and the curated CSV writer                    #
# --------------------------------------------------------------------------- #
def fov_row(template, t, cx, cy, object_id, size=1.0):
    """A TargetData row templated on a real one (so MetaXpress ingests it as a valid
    target), with target T<t>'s bbox centred at (cx, cy) and a fresh object_id.
    MetaXpress positions off the bbox; `size` is only a nominal marker footprint."""
    row = dict(template)
    row["object_id"] = str(object_id)
    for col, val in zip(_bb_cols(t), (cx - size / 2.0, cy - size / 2.0, size, size)):
        row[col] = repr(val)
    return row


def _write_csv(path, header, rows):
    """Write a CRLF CSV as unicode with minimal RFC-4180 quoting. Hand-rolled
    because Jython's csv module cannot write to a text stream."""
    def field(v):
        s = u"{0}".format(v)
        return u'"' + s.replace(u'"', u'""') + u'"' if any(c in s for c in u',"\n\r') else s
    with io.open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write(u",".join(field(c) for c in header) + u"\r\n")
        for r in rows:
            fh.write(u",".join(field(r.get(c, u"")) for c in header) + u"\r\n")


def _read_bytes(path):
    with open(path, "rb") as fh:
        return fh.read()


def _assert_previous_curation(cur_dir, audit_dir):
    """Refuse to overwrite TargetData/ on a rerun unless it is still byte-for-byte the
    last TargetData_curated/ mirror - so newly regenerated IN Carta output (which would
    NOT match the mirror) is never silently clobbered."""
    if not os.path.isdir(cur_dir):
        return
    if not (os.path.isdir(audit_dir) and os.path.isfile(audit_dir + "curation_changes.csv")):
        raise ValueError("TargetData_original/ exists but the TargetData_curated/ mirror is "
                         "missing; refusing to overwrite TargetData/ (it may be new IN Carta output).")
    for name in os.listdir(cur_dir):
        src = cur_dir + name
        if not os.path.isfile(src):
            continue
        mirror = audit_dir + name
        if not os.path.isfile(mirror) or _read_bytes(src) != _read_bytes(mirror):
            raise ValueError("TargetData/ differs from the last TargetData_curated/ mirror; "
                             "refusing to overwrite possible new IN Carta output: " + name)


def write_curated_output(results_path, gate_spec, fov_px, sample_size=5, seed=42,
                         neighbourhood=2.0, stage_margin=0.05, run_note="", base=None, progress=None):
    """Curate and write the result IN PLACE, following the v1 folder convention.

    The first run renames IN Carta's `TargetData/` to `TargetData_original/` and never
    touches it again; every run reads from there and rewrites `TargetData/` from scratch
    with the curated per-site CSVs (the generated FOV-centre rows, one file per site,
    using the highest-T channel required by MetaXpress). The same files are mirrored into
    `TargetData_curated/`, with a `curation_changes.csv` audit log. A rerun is refused
    unless `TargetData/` still matches the previous mirror, so regenerated IN Carta
    output is not clobbered. `run_note` is the timestamp string for the log (the caller
    supplies it, so this stays portable). Returns (res, cur_dir, audit_dir)."""
    sep = os.sep
    base_dir = results_path if results_path.endswith(sep) else results_path + sep
    orig_dir = base_dir + "TargetData_original" + sep
    cur_dir = base_dir + "TargetData" + sep
    audit_dir = base_dir + "TargetData_curated" + sep

    if not os.path.isdir(orig_dir):                       # first run: preserve the originals
        if not os.path.isdir(cur_dir):
            raise ValueError("No TargetData folder to curate at: " + cur_dir)
        os.rename(cur_dir.rstrip(sep), orig_dir.rstrip(sep))
    else:                                                 # rerun: guard against clobbering new data
        _assert_previous_curation(cur_dir, audit_dir)

    for d in (cur_dir, audit_dir):                        # rewrite both from scratch (files only)
        if os.path.isdir(d):
            for f in os.listdir(d):
                if os.path.isfile(d + f):
                    os.remove(d + f)
        else:
            os.makedirs(d)

    res = curate(orig_dir, gate_spec, fov_px=fov_px, sample_size=sample_size, seed=seed,
                 neighbourhood=neighbourhood, stage_margin=stage_margin, base=base, progress=progress)
    base, output, header = res["base"], res["output"], res["header"]

    changes = [u"ZMB MD HCS target curation changes",
               u"Run\t" + run_note,
               u"Selected object\t" + base,
               u"Acquisition CSV template\t" + output,
               u"Gate\t" + gate_spec,
               u"Cells per site/field\t%d" % sample_size,
               u"Neighbourhood (%% of radius)\t%g" % (neighbourhood * 100.0),
               u"Stage margin (%% per side)\t%g" % (stage_margin * 100.0),
               u"Seed\t%d" % seed,
               u"FOV size (montage px)\t%g" % fov_px,
               u"Originals\t" + orig_dir,
               u"Operational TargetData\t" + cur_dir,
               u"Curated mirror\t" + audit_dir,
               u"",
               u"file\tpositives\teligible\tacquired\textra\tfovs\tstatus"]
    for well in res["wells"]:
        fname = output + _SITE_KEY + well["site"] + ".csv"
        _write_csv(cur_dir + fname, header, well["fov_rows"])
        _write_csv(audit_dir + fname, header, well["fov_rows"])
        changes.append(u"%s\t%d\t%d\t%d\t%d\t%d\t%s" % (
            fname, len(well["selected"]), well["eligible"], len(well["acquired"]),
            well["extra"], len(well["tiles"]), well["status"]))
    with io.open(audit_dir + "curation_changes.csv", "w", encoding="utf-8", newline="") as fh:
        fh.write(u"\r\n".join(changes) + u"\r\n")
    return res, cur_dir, audit_dir


# --------------------------------------------------------------------------- #
# 10. Curate - run the whole pipeline over every acquisition site/field       #
# --------------------------------------------------------------------------- #
def curate(results_dir, gate_spec, fov_px=256.0, sample_size=5, seed=42,
           neighbourhood=2.0, stage_margin=0.05, out_csv=None, base=None, progress=None):
    """Curate every site/field (gate -> SURS -> disjoint FOVs). `base` is the object
    class to gate (default: the first channel). Acquisition rows use the highest-T
    channel's schema, independently of the selected object. The SURS grid spans the whole scanned
    area (all base objects), so the sample is spread over the site the overview imaged,
    not just where the positives landed. The per-site seed is `seed + well_index`, so
    sites are reproducible yet decorrelated (SplitMix64 makes consecutive seeds
    independent). Returns a dict with per-site results and generated FOV rows; writes
    the curated CSV if out_csv is given."""
    target_dir = find_target_dir(results_dir)
    name_index, order = discover_signals(target_dir)
    base = base if base in name_index else order[0]
    base_t = name_index[base]
    output = max(order, key=lambda name: (name_index[name], name))
    output_t = name_index[output]
    require = parse_gate(gate_spec)
    unknown = [s for s in require if s not in name_index]
    if unknown:
        raise ValueError("Gate signal(s) not found: %s (available: %s)" % (unknown, order))

    files_by_site = {}
    for fn in sorted(os.listdir(target_dir)):
        if fn.lower().endswith(".csv"):
            files_by_site.setdefault(_site_of(fn), {})[signal_name(fn)] = fn

    header, wells, fov_rows = None, [], []
    sites = sorted(files_by_site)
    for well_index, site in enumerate(sites):
        if progress:
            progress(well_index + 1, len(sites))       # for the caller's status/heartbeat
        files = files_by_site[site]
        if base not in files:
            continue
        if output not in files:
            raise ValueError("Highest-T acquisition channel %s is missing for site %s" %
                             (output, site))

        head, base_rows = read_csv(os.path.join(target_dir, files[base]))
        base_boxes = boxes(base_rows, base_t)
        output_head, output_rows = read_csv(os.path.join(target_dir, files[output]))
        header = header or output_head

        gate_boxes = {}
        for sig in require:
            rows = read_csv(os.path.join(target_dir, files[sig]))[1] if sig in files else []
            gate_boxes[sig] = boxes(rows, name_index[sig])
        selected = gate(base_boxes, gate_boxes, require)

        acquired, tiles, n_eligible = select_and_place(
            selected, fov_px, stage_margin, neighbourhood, sample_size,
            seed + well_index, search_boxes=base_boxes)

        # positives imaged whole by the FOVs = the sampled targets PLUS any other
        # positive cell that happens to fall entirely inside a field (free extra
        # observations, especially at low mag where the field is large)
        half = _usable_half(fov_px, stage_margin)
        captured_boxes = [b for b in selected
                          if any(_fully_inside(b, t["cx"], t["cy"], half) for t in tiles)]
        extra_boxes = [b for b in captured_boxes if b not in acquired]   # captured but not sampled

        template = output_rows[0] if output_rows else {}
        well_fov_rows = [fov_row(template, output_t, t["cx"], t["cy"],
                                 "FOV_%d" % (len(fov_rows) + i + 1))
                         for i, t in enumerate(tiles)]
        fov_rows.extend(well_fov_rows)
        status = site_status(n_eligible, len(acquired), sample_size)
        wells.append({"site": site, "base": base_boxes, "selected": selected,
                      "acquired": acquired, "captured": len(captured_boxes),
                      "extra": len(extra_boxes), "extra_boxes": extra_boxes,
                      "eligible": n_eligible, "status": status,
                      "sample_short": status != "full",
                      "tiles": tiles, "fov_rows": well_fov_rows})

    if out_csv and header:
        _write_csv(out_csv, header, fov_rows)
    return {"target_dir": target_dir, "base": base, "output": output, "signals": order,
            "wells": wells, "fov_rows": fov_rows, "header": header}


# --------------------------------------------------------------------------- #
# 11. Fiji front end - class dialog, visual report, and the run entry point    #
# ImageJ imports live inside these functions so the module still imports in     #
# plain CPython (for the tests).                                               #
# --------------------------------------------------------------------------- #
def _render_report(res, report_dir, fov_px):
    """One overlay PNG per site (grey nuclei / orange positives / red sampled / blue
    FOV boxes) written under the curated mirror. Fiji-only (ImageJ ColorProcessor)."""
    from ij import IJ, ImagePlus
    from ij.process import ColorProcessor
    from java.awt import Color
    canvas = 1040
    max_x, max_y = dataset_extent(res["wells"])
    scale = min(canvas / max_x, canvas / max_y)
    grey, orange, green, red, blue = (Color(225, 225, 225), Color(244, 165, 130),
                                      Color(26, 152, 80), Color(202, 0, 32), Color(5, 113, 176))

    def render(well):
        ip = ColorProcessor(canvas, canvas)
        ip.setColor(Color.WHITE); ip.fill()

        def dots(bs, colour, r):
            ip.setColor(colour)
            for b in bs:
                ip.fillOval(int((b[0] + b[2] / 2.0) * scale) - r,
                            int((b[1] + b[3] / 2.0) * scale) - r, 2 * r, 2 * r)
        dots(well["base"], grey, 1)          # all nuclei
        dots(well["selected"], orange, 3)    # positive population
        dots(well["extra_boxes"], green, 3)  # bonus positives also captured in a FOV
        dots(well["acquired"], red, 4)       # the sampled cells we image
        ip.setColor(blue); ip.setLineWidth(2)
        f = int(fov_px * scale)
        for t in well["tiles"]:
            ip.drawRect(int(t["cx"] * scale - f / 2.0), int(t["cy"] * scale - f / 2.0), f, f)
        IJ.saveAs(ImagePlus(well["site"], ip), "PNG",
                  os.path.join(report_dir, "report_%s.png" % well["site"]))

    wells = res["wells"]
    for i, well in enumerate(wells):
        IJ.showStatus("MD HCS curation: rendering report %d/%d" % (i + 1, len(wells)))
        IJ.showProgress(i, len(wells))
        if (i + 1) == len(wells) or (i + 1) % 10 == 0:
            IJ.log("  ...rendered %d/%d reports" % (i + 1, len(wells)))
        render(well)
    IJ.showProgress(1.0)


def _render_plate_overview(res, path, fov_px, overview_desc=""):
    """Whole-plate montage in physical well coordinates, with separate Z/T facets."""
    from ij import IJ, ImagePlus
    from ij.process import ColorProcessor
    from java.awt import Color, Font

    if not res["wells"]:
        return

    layout = build_overview_layout(res["wells"])
    if not layout:
        return

    panel, well_top = 300, 19
    grey, orange, green, red, blue = (Color(210, 210, 210), Color(244, 165, 130),
                                      Color(26, 152, 80), Color(202, 0, 32), Color(5, 113, 176))
    by_site = dict((w["site"], w) for w in res["wells"])
    max_x, max_y = dataset_extent(res["wells"])
    origins = read_field_origins(res["target_dir"], overview_pixel_size(overview_desc))
    saved = []

    for facet in sorted(layout):
        slots = layout[facet]
        n_well_rows = max(v["well_row"] for v in slots.values()) + 1
        n_well_cols = max(v["well_col"] for v in slots.values()) + 1
        ip = ColorProcessor(n_well_cols * panel, n_well_rows * panel)
        ip.setColor(Color.WHITE); ip.fill()

        # Draw each physical well frame and header once.
        physical = {}
        for site, slot in slots.items():
            physical.setdefault((slot["well_row"], slot["well_col"]), slot)
        processors = {}
        for panel_key, slot in physical.items():
            wp = ColorProcessor(panel, panel)
            wp.setColor(Color.WHITE); wp.fill()
            wp.setFont(Font("SansSerif", Font.PLAIN, 12))
            wp.setColor(Color(150, 150, 150)); wp.drawRect(0, 0, panel - 1, panel - 1)
            wp.setColor(Color.BLACK)
            wp.drawString("Row %d Col %d" % (slot["row"], slot["col"]), 5, 14)
            processors[panel_key] = wp

        # Use measured stage offsets when available. The deterministic fallback keeps
        # reports usable for exported datasets that no longer have acquisition metadata.
        site_origins = {}
        for site, slot in slots.items():
            key = parse_site(site)
            measured = origins.get(key)
            if measured is None:
                measured = (slot["field_col"] * max_x, slot["field_row"] * max_y)
            site_origins[site] = measured
        extent_x = max(site_origins[s][0] + max_x for s in slots)
        extent_y = max(site_origins[s][1] + max_y for s in slots)
        scale = min((panel - 8.0) / extent_x, (panel - well_top - 4.0) / extent_y)

        # Draw source-field footprints first. Overlap is intentional and represents
        # the actual acquisition geometry within the well.
        for site, slot in slots.items():
            key = parse_site(site)
            wp = processors[(slot["well_row"], slot["well_col"])]
            fx = 4 + int(site_origins[site][0] * scale)
            fy = well_top + int(site_origins[site][1] * scale)
            wp.setColor(Color(190, 190, 190))
            wp.drawRect(fx, fy, int(max_x * scale), int(max_y * scale))
            wp.setColor(Color.BLACK)
            wp.setFont(Font("SansSerif", Font.PLAIN, 9))
            wp.drawString("F%d" % key[2], fx + 3, fy + 10)

        for site, slot in slots.items():
            w = by_site[site]
            wp = processors[(slot["well_row"], slot["well_col"])]
            origin_x, origin_y = site_origins[site]
            ox, oy = 4, well_top

            def put(bs, colour, rad):
                wp.setColor(colour)
                for b in bs:
                    x = ox + int((origin_x + b[0] + b[2] / 2.0) * scale)
                    y = oy + int((origin_y + b[1] + b[3] / 2.0) * scale)
                    wp.fillOval(x - rad, y - rad, 2 * rad + 1, 2 * rad + 1)

            put(w["base"], grey, 0)
            put(w["selected"], orange, 1)
            put(w["extra_boxes"], green, 2)
            put(w["acquired"], red, 2)
            wp.setColor(blue)
            f = int(fov_px * scale)
            for tile in w["tiles"]:
                wp.drawRect(ox + int((origin_x + tile["cx"]) * scale - f / 2.0),
                            oy + int((origin_y + tile["cy"]) * scale - f / 2.0), f, f)

        for (well_row, well_col), wp in processors.items():
            ip.insert(wp, well_col * panel, well_row * panel)

        if len(layout) == 1:
            out_path = path
        else:
            stem, ext = os.path.splitext(path)
            out_path = "%s_Z%d_T%d%s" % (stem, facet[0], facet[1], ext or ".png")
        IJ.saveAs(ImagePlus("plate_Z%d_T%d" % facet, ip), "PNG", out_path)
        saved.append(out_path)
    return saved


def _run_curation(results_path, gate_spec, base, objective, overview_desc,
                  sample_size, seed, neighbourhood, stage_margin):
    """Infer the FOV size from the objective + overview, write the curated output in
    place (v1 convention), and render the per-well overlays. run_macro collects the
    arguments interactively; this does the work so it can also be driven headlessly."""
    from ij import IJ
    from java.util import Date

    mag_t, changer_t = _OBJECTIVES[objective]
    mag_ov, changer_ov, binning_ov, region_px = overview_scale(overview_desc)
    if mag_ov <= 0:
        IJ.error("MD HCS target curation",
                 "Could not read the overview objective/scale.\n"
                 "Expected the InCarta result_metadata.csv (greyscale_image) and the overview\n"
                 "image reachable by walking up from the results folder - is the acquisition\n"
                 "(experiment_montage/timepoint0) present next to the analysis?")
        return None
    fov_px = fov_px_for(mag_t, changer_t, mag_ov, changer_ov, binning_ov, region_px)
    IJ.log("MD HCS curation: gate='%s'" % gate_spec)
    IJ.log("  overview %gx changer %gx binning %d -> target %s FOV=%.0f montage px"
           % (mag_ov, changer_ov, binning_ov, objective, fov_px))

    def tick(done, total):                             # Log heartbeat + status-bar progress
        IJ.showProgress(done, total)
        if done == total or done % 10 == 0:
            IJ.log("  ...curated %d/%d sites" % (done, total))

    IJ.log("  working - curating sites (reading CSVs, sampling, placing FOVs)...")
    IJ.showStatus("MD HCS curation: curating sites...")
    res, cur_dir, audit_dir = write_curated_output(
        results_path, gate_spec, fov_px, sample_size, seed, neighbourhood, stage_margin,
        run_note=str(Date()), base=base, progress=tick)

    IJ.log("  working - rendering %d per-site reports + plate overview..." % len(res["wells"]))
    report_dir = os.path.join(audit_dir, "report")
    if not os.path.isdir(report_dir):
        os.makedirs(report_dir)
    _render_report(res, report_dir, fov_px)
    IJ.showStatus("MD HCS curation: rendering plate overview...")
    _render_plate_overview(res, os.path.join(audit_dir, "plate_overview.png"),
                           fov_px, overview_desc)
    IJ.showStatus("MD HCS curation: done - see the Log")

    acquired = sum(len(w["acquired"]) for w in res["wells"])
    extra = sum(w["extra"] for w in res["wells"])
    fovs = sum(len(w["tiles"]) for w in res["wells"])
    under = sum(1 for w in res["wells"] if w["status"] != "full")
    constrained = sum(1 for w in res["wells"] if w["status"] == "constrained")
    IJ.log("  sites=%d acquired=%d extra-whole=%d FOVs=%d under-sampled=%d constrained=%d"
           % (len(res["wells"]), acquired, extra, fovs, under, constrained))
    IJ.log("  curated TargetData -> %s" % cur_dir)
    IJ.log("  mirror + changes   -> %s" % audit_dir)
    return res


def run_macro():
    from ij import IJ
    from ij.gui import GenericDialog

    results_path = results_dir.getAbsolutePath()
    IJ.log("MD HCS curation: results folder = " + results_path)
    overview_desc = overview_from_results(results_path)     # found via the InCarta metadata
    need_overview = overview_scale(overview_desc)[0] <= 0   # not reachable next to the results

    # discover the marker classes from the ORIGINAL data (preserved on the first run)
    orig = os.path.join(results_path, "TargetData_original")
    scan_dir = orig if os.path.isdir(orig) else find_target_dir(results_path)
    _, order = discover_signals(scan_dir)          # the classes (groups) found in the data

    gd = GenericDialog("MD HCS target curation")
    gd.addMessage("Classes  -  mark ONE as the object to gate,")                    # section 2
    gd.addMessage("the rest positive / negative / ignore:")
    for sig in order:                              # one pull-down per class, however many
        gd.addChoice(sig, ["object", "positive", "negative", "ignore"],
                     "object" if sig == order[0] else "ignore")
    gd.addMessage("Target acquisition")                                   # section 3
    gd.addChoice("Objective", _OBJECTIVE_ORDER, "60x")
    if need_overview:                                                     # fallback if not auto-found
        gd.addMessage("Overview image not found under the selected results folder:")
        gd.addMessage(results_path)
        gd.addMessage("(Cancel and pick the original results folder, or select the image below.)")
        gd.addFileField("Overview image", "")
    gd.addMessage("Parameters")                                           # section 4
    gd.addNumericField("Cells per site/field (target)", 5, 0)
    gd.addNumericField("Neighbourhood (% of cell radius)", 200, 0)
    gd.addNumericField("Stage margin (% per side)", 5, 0)
    gd.addNumericField("Random seed", 42, 0)
    gd.showDialog()
    if gd.wasCanceled():
        return

    roles = dict((sig, gd.getNextChoice()) for sig in order)
    objective = gd.getNextChoice()
    if need_overview:
        overview_desc = read_image_description(gd.getNextString())
    sample_size = int(gd.getNextNumber())
    neighbourhood = gd.getNextNumber() / 100.0
    stage_margin = gd.getNextNumber() / 100.0
    seed = int(gd.getNextNumber())

    objects = [s for s in order if roles[s] == "object"]      # the class marked 'object'
    base = objects[0] if objects else order[0]
    markers = [s for s in order if s != base]
    gate_spec = build_gate_spec(markers, roles)               # markers -> yes/no/ignore

    _run_curation(results_path, gate_spec, base, objective, overview_desc,
                  sample_size, seed, neighbourhood, stage_margin)


# --------------------------------------------------------------------------- #
# 12. Tests - `python target_curation.py` in CPython; also runs in Fiji Jython #
# --------------------------------------------------------------------------- #
def run_tests():
    fails = []

    def check(name, got, want):
        ok = got == want
        print(("  PASS " if ok else "  FAIL ") + name)
        if not ok:
            print("        got  %r\n        want %r" % (got, want))
            fails.append(name)

    def bx(x, y, w, h):
        return (float(x), float(y), float(w), float(h))

    def disjoint(tiles, fov):
        return all(abs(a["cx"] - b["cx"]) >= fov - 1e-9 or abs(a["cy"] - b["cy"]) >= fov - 1e-9
                   for p, a in enumerate(tiles) for b in tiles[p + 1:])

    # Goldens are the values this code produces; pinned so a change is noticed, and
    # verified identical in CPython and Fiji Jython.
    GOLDEN_TILES = [(800.5, 400.5), (200.5, 1000.5), (1200.5, 600.5),
                    (1000.5, 1400.5), (200.5, 1400.5)]
    GOLDEN_ACQUIRED, GOLDEN_FOVS, GOLDEN_CAPTURED, GOLDEN_EXTRA = 165, 164, 313, 148

    print("gating")
    n1, n2 = bx(0, 0, 10, 10), bx(100, 0, 10, 10)
    ch = [bx(5, 5, 4, 4)]                                   # touches n1 only
    check("edge-only contact is not overlap", touches(bx(0, 0, 10, 10), bx(10, 0, 5, 5)), False)
    check("yes gate", gate([n1, n2], {"m": ch}, {"m": True}), [n1])
    check("no gate", gate([n1, n2], {"m": ch}, {"m": False}), [n2])
    check("gate builder -> spec (channel order, ignore dropped)",
          build_gate_spec(["A", "B", "C"], {"A": "positive", "B": "ignore", "C": "negative"}), "A:yes; C:no")
    check("gate builder all-ignore -> empty", build_gate_spec(["A"], {"A": "ignore"}), "")
    check("gate builder ignores the 'object' role",
          build_gate_spec(["A", "B"], {"A": "object", "B": "positive"}), "B:yes")
    check("gate builder round-trips via parse_gate",
          parse_gate(build_gate_spec(["mScarlet cells"], {"mScarlet cells": "positive"})), {"mScarlet cells": True})

    print("site parsing + overview layout")
    check("full site parsed", parse_site("R3-C2-F1-Z0-T4"), (3, 2, 1, 0, 4))
    malformed = False
    try:
        parse_site("R3-C2")
    except ValueError:
        malformed = True
    check("malformed site rejected", malformed, True)
    check("status full", site_status(5, 5, 5), "full")
    check("status low cells", site_status(4, 4, 5), "low_cells")
    check("status constrained", site_status(6, 4, 5), "constrained")

    layout_wells = [{"site": "R3-C2-F%d-Z0-T0" % f} for f in range(4)]
    laid = build_overview_layout(layout_wells)[(0, 0)]
    field_slots = dict((parse_site(site)[2], (v["field_row"], v["field_col"]))
                       for site, v in laid.items())
    check("deterministic field grid", field_slots,
          {0: (0, 0), 1: (0, 1), 2: (1, 0), 3: (1, 1)})
    check("four fields get unique slots", len(set(field_slots.values())), 4)
    uneven = build_overview_layout(
        layout_wells + [{"site": "R3-C3-F0-Z0-T0"}])[(0, 0)]
    check("all wells share one field grid",
          (uneven["R3-C3-F0-Z0-T0"]["field_rows"],
           uneven["R3-C3-F0-Z0-T0"]["field_cols"]), (2, 2))
    check("dataset extent spans every site",
          dataset_extent([{"base": [(0, 0, 10, 20)]},
                          {"base": [(100, 50, 5, 7)]}]), (105, 57))
    metadata_rows = [
        {"Row": "3", "Column": "2", "Field": "0", "ZIndex": "0", "Timepoint": "0",
         "PositionXUm": "10", "PositionYUm": "20"},
        {"Row": "3", "Column": "2", "Field": "1", "ZIndex": "0", "Timepoint": "0",
         "PositionXUm": "10", "PositionYUm": "30"},
        {"Row": "3", "Column": "2", "Field": "2", "ZIndex": "0", "Timepoint": "0",
         "PositionXUm": "25", "PositionYUm": "30"}]
    check("stage metadata becomes well-local pixel origins",
          field_origins(metadata_rows, 0.5),
          {(3, 2, 0, 0, 0): (0.0, 0.0),
           (3, 2, 1, 0, 0): (0.0, 20.0),
           (3, 2, 2, 0, 0): (30.0, 20.0)})

    faceted = build_overview_layout(
        [{"site": "R3-C2-F0-Z0-T0"}, {"site": "R3-C2-F0-Z1-T0"},
         {"site": "R3-C2-F0-Z0-T1"}])
    check("Z/T variants are separate facets", sorted(faceted), [(0, 0), (0, 1), (1, 0)])

    print("SURS + disjoint placement")
    check("oversized cell not eligible", eligible([bx(0, 0, 500, 500)], 100, 0.0, 0.0), [])
    _, tiles, _ = select_and_place([bx(0, 0, 10, 10), bx(20, 0, 10, 10)], 100, 0.0, 0.0, 5, 1)
    check("close pair -> 1 FOV, 2 cells", (len(tiles), len(tiles[0]["covered"])), (1, 2))
    _, tiles, _ = select_and_place([bx(0, 0, 10, 10), bx(1000, 0, 10, 10)], 100, 0.0, 0.0, 5, 1)
    check("far pair -> 2 disjoint FOVs", (len(tiles), disjoint(tiles, 100)), (2, True))
    # reject-and-replace: sample all three; A and B (70 apart, usable 60) cannot share
    # and their FOVs would overlap, so B is dropped and the two disjoint FOVs frame A
    # and the far C.
    got, tiles, _ = select_and_place([bx(0, 0, 4, 4), bx(70, 0, 4, 4), bx(1000, 0, 4, 4)],
                                     100, 0.2, 0.0, 3, 5)
    check("overlap rejected + another cell taken",
          (len(tiles), disjoint(tiles, 100), set(b[0] for b in got)), (2, True, set([0.0, 1000.0])))
    # a spread field (cells 200 px apart so no two share/overlap a 100 px FOV)
    field = [bx(x * 200, y * 200, 1, 1) for x in range(10) for y in range(10)]
    _, ftiles, _ = select_and_place(field, 100, 0.0, 0.0, 9, 3)
    check("all placed FOVs disjoint", disjoint(ftiles, 100), True)
    # per-area uniform: the sample tracks AREA, not cell density. Here 80% of the
    # cells are packed into one tiny clump (a single grid frame) and 25 sit spread
    # out; a per-CELL sample would draw ~80% of its picks from the clump, but the
    # grid draws only that one frame's worth (~1 of 25 = ~4%).
    bg = [bx(x * 800, y * 800, 1, 1) for x in range(5) for y in range(5)]
    clump = [bx(100 + (k % 10) * 2, 100 + (k // 10) * 2, 1, 1) for k in range(100)]
    mixed = bg + clump
    in_clump = set(range(len(bg), len(mixed)))
    cx = dict((i, mixed[i][0] + 0.5) for i in range(len(mixed)))
    cy = dict((i, mixed[i][1] + 0.5) for i in range(len(mixed)))
    per_trial = [sum(1 for i in _surs_sample(list(range(len(mixed))), cx, cy, 25, sd, mixed) if i in in_clump)
                 for sd in range(100)]
    clump_frac = sum(per_trial) / 100.0 / 25.0
    print("    clump share of sample = %.0f%% (its cells are 80%% of all; per-cell would be ~80%%)"
          % (100 * clump_frac))
    check("per-area uniform is not biased by dense clumps", clump_frac < 0.15, True)
    g7 = [(t["cx"], t["cy"]) for t in select_and_place(field, 100, 0.0, 0.0, 5, 7)[1]]
    print("    tiles(seed 7) = %s" % g7)
    check("requested sample is a hard maximum", len(g7), 5)
    check("reproducible", g7, [(t["cx"], t["cy"]) for t in select_and_place(field, 100, 0.0, 0.0, 5, 7)[1]])
    check("portable golden (seed 7)", g7, GOLDEN_TILES)

    print("fov_row")
    row = fov_row({"object_id": "1", "well_label": "B - 3"}, 1, 100.0, 200.0, "FOV_1", 2.0)
    check("fresh id + identity kept", (row["object_id"], row["well_label"]), ("FOV_1", "B - 3"))
    check("bbox centred at point",
          (float(row["T1$AS_FID_Blob_BoundingBoxX"]), float(row["T1$AS_FID_Blob_BoundingBoxY"])), (99.0, 199.0))

    print("objective / overview inference")
    desc = ('<MetaData><prop id="Description" type="string" value="Binning: 1 x 1 '
            'Region: 2304 x 2304, offset at (0, 0)" />'
            '<prop id="camera-binning-x" type="int" value="1" />'
            '<prop id="HCS.ai Objective" type="string" value="10X Plan Apo Lambda D" />'
            '<prop id="HCS.ai Tube Lens" type="string" value="1X" /></MetaData>')
    scale = overview_scale(desc)
    check("overview scale parsed (10x, 1x, bin 1, 2304)", scale, (10.0, 1.0, 1, 2304))
    fov = dict((name, fov_px_for(m, c, scale[0], scale[1], scale[2], scale[3]))
               for name, (m, c) in _OBJECTIVES.items())
    print("    inferred FOV px (10x overview) = %s"
          % ", ".join("%s:%g" % (k, fov[k]) for k in ("2x", "4x", "10x", "20x", "40x", "60x")))
    check("2x  -> 11520 px", fov["2x"], 11520.0)
    check("20x -> 1152 px", fov["20x"], 1152.0)
    check("40x -> 576 px", fov["40x"], 576.0)
    check("60x -> 384 px", fov["60x"], 384.0)
    check("60x + 1.5x changer -> 256 px", fov["60x + 1.5x changer"], 256.0)
    real_ov = os.path.join(
        r"Z:\transfer\Thom\Nico MD\NB26-15_Overview10x_DAPI-mScarlet_20260713_152811",
        "experiment_montage", "timepoint0",
        "NB26-15_Overview10x_DAPI-mScarlet_t0_C09_s0_w0_z0.tif")
    if os.path.isfile(real_ov):
        check("real overview TIFF scale read", overview_scale(read_image_description(real_ov)),
              (10.0, 1.0, 1, 2304))
    else:
        print("  SKIP real overview TIFF (not present)")

    print("overview auto-discovery from the results folder")
    import tempfile, shutil
    ov_tmp = tempfile.mkdtemp()
    try:
        analysis = os.path.join(ov_tmp, "montage", "Results", "analysis")
        os.makedirs(os.path.join(analysis, "TargetData"))
        with io.open(os.path.join(analysis, "TargetData", "DAPI_singleTargetData_S.csv"),
                     "w", encoding="utf-8", newline="") as fh:
            fh.write(u"object_id,T1$AS_FID_Blob_BoundingBoxX,T1$AS_FID_Blob_BoundingBoxY,"
                     u"T1$AS_FID_Blob_BoundingBoxWidth,T1$AS_FID_Blob_BoundingBoxHeight\r\n")
        with io.open(os.path.join(analysis, "result_metadata.csv"), "w", encoding="utf-8", newline="") as fh:
            fh.write(u"well_label,greyscale_image\r\nB - 3,timepoint0\\ov.tif\r\n")   # names it, 2 up
        os.makedirs(os.path.join(ov_tmp, "montage", "timepoint0"))
        ov = os.path.join(ov_tmp, "montage", "timepoint0", "ov.tif")
        with io.open(ov, "w", encoding="utf-8") as fh:
            fh.write(u"x")
        check("overview image found by walking up from results", _overview_tiff_path(analysis), ov)
    finally:
        shutil.rmtree(ov_tmp, ignore_errors=True)
    real_analysis = os.path.join(
        r"Z:\transfer\Thom\Nico MD\NB26-15_Overview10x_DAPI-mScarlet_20260713_152811",
        "experiment_montage", "Results", "mScarlet Cells_2026-Jul-13-17-23-33-077")
    if os.path.isdir(os.path.join(real_analysis, "TargetData")):
        check("overview auto-found + scaled from real results",
              overview_scale(overview_from_results(real_analysis)), (10.0, 1.0, 1, 2304))
    else:
        print("  SKIP overview auto-discovery on real results (not present)")

    print("v1 output convention (synthetic dataset in a temp folder)")
    import tempfile, shutil
    tmp = tempfile.mkdtemp()
    try:
        td = os.path.join(tmp, "TargetData")
        os.makedirs(td)
        site = "R2-C03-F0-Z0-T0"
        dapi = ["object_id,well_label,T1$AS_FID_Blob_BoundingBoxX,T1$AS_FID_Blob_BoundingBoxY,"
                "T1$AS_FID_Blob_BoundingBoxWidth,T1$AS_FID_Blob_BoundingBoxHeight"]
        dapi += ["%d,B - 3,%d,0,10,10" % (i, i * 100) for i in range(4)]     # 4 nuclei
        ms = ["object_id,T2$AS_FID_Blob_BoundingBoxX,T2$AS_FID_Blob_BoundingBoxY,"
              "T2$AS_FID_Blob_BoundingBoxWidth,T2$AS_FID_Blob_BoundingBoxHeight"]
        ms += ["%d,%d,2,6,6" % (i, i * 100 + 2) for i in range(3)]           # positive on nuclei 0,1,2
        for name, lines in (("DAPI", dapi), ("mScarlet cells", ms)):
            with io.open(os.path.join(td, "%s_singleTargetData_%s.csv" % (name, site)),
                         "w", encoding="utf-8", newline="") as fh:
                fh.write(u"\r\n".join(lines) + u"\r\n")

        res, cur_dir, audit_dir = write_curated_output(tmp, "mScarlet cells:yes", 100.0,
                                                       sample_size=5, seed=1, run_note="t")
        n_fov = len(res["wells"][0]["tiles"])
        curated = os.path.join(cur_dir, "mScarlet cells_singleTargetData_%s.csv" % site)
        check("first run preserves originals", os.path.isdir(os.path.join(tmp, "TargetData_original")), True)
        check("curated per-site file written", os.path.isfile(curated), True)
        check("mirror + changes log written", os.path.isfile(os.path.join(audit_dir, "curation_changes.csv")), True)
        check("highest-T target lands in TargetData/", os.path.isfile(curated), True)
        check("lower-T base target is not emitted",
              os.path.isfile(os.path.join(cur_dir, "DAPI_singleTargetData_%s.csv" % site)), False)
        ch, rr = read_csv(curated)
        check("curated schema is highest T", "T2$AS_FID_Blob_BoundingBoxX" in ch, True)
        check("curated schema excludes base T", "T1$AS_FID_Blob_BoundingBoxX" in ch, False)
        check("3 positives -> 3 disjoint FOV rows", (len(rr), n_fov), (3, 3))
        check("result records independent selection/output roles",
              (res["base"], res["output"]), ("DAPI", "mScarlet cells"))
        check("synthetic site status", res["wells"][0]["status"], "low_cells")
        res2, _, _ = write_curated_output(tmp, "mScarlet cells:yes", 100.0, sample_size=5, seed=1, run_note="t2")
        check("rerun allowed (matches mirror) + reproduces", len(res2["wells"][0]["tiles"]), n_fov)
        rb = curate(os.path.join(tmp, "TargetData_original"), "DAPI:yes", fov_px=100.0, base="mScarlet cells")
        check("curate honours a chosen base object (not just the first channel)", rb["base"], "mScarlet cells")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    print("real data (auto-skip)")
    real = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Results",
                        "mScarlet Cells_2026-Jul-13-17-23-33-077")
    if os.path.isdir(os.path.join(real, "TargetData")):
        res = curate(real, "mScarlet cells:yes")
        pos = sum(len(w["selected"]) for w in res["wells"])
        acquired = sum(len(w["acquired"]) for w in res["wells"])
        captured = sum(w["captured"] for w in res["wells"])
        extra = sum(w["extra"] for w in res["wells"])
        fovs = sum(len(w["tiles"]) for w in res["wells"])
        overlaps = sum(1 for wl in res["wells"] for a in range(len(wl["tiles"]))
                       for b in range(a + 1, len(wl["tiles"]))
                       if abs(wl["tiles"][a]["cx"] - wl["tiles"][b]["cx"]) < 256
                       and abs(wl["tiles"][a]["cy"] - wl["tiles"][b]["cy"]) < 256)
        print("    plate: positives=%d acquired=%d extra=%d captured=%d FOVs=%d overlaps=%d"
              % (pos, acquired, extra, captured, fovs, overlaps))
        check("signals discovered", res["signals"], ["DAPI", "mScarlet cells"])
        check("plate positives (golden)", pos, 4253)
        check("no overlapping FOVs on the plate", overlaps, 0)
        check("plate acquired (golden)", acquired, GOLDEN_ACQUIRED)
        check("captured == acquired + extra", captured, acquired + extra)
        check("plate extra/bonus (golden)", extra, GOLDEN_EXTRA)
        check("extra_boxes count matches extra", sum(len(w["extra_boxes"]) for w in res["wells"]), GOLDEN_EXTRA)
        check("plate captured (golden)", captured, GOLDEN_CAPTURED)
        check("plate FOVs (golden)", fovs, GOLDEN_FOVS)
    else:
        print("  SKIP (dataset not present)")

    print("\n" + ("ALL PASS" if not fails else "%d FAILED: %s" % (len(fails), fails)))
    return 1 if fails else 0


# --------------------------------------------------------------------------- #
# 13. Dispatch - the macro inside Fiji (SciJava binds `results_dir`),          #
# otherwise the tests (CPython `python target_curation.py`).                   #
# --------------------------------------------------------------------------- #
def _in_fiji():
    try:
        results_dir
        return True
    except NameError:
        return False


if _in_fiji():
    run_macro()
elif __name__ == "__main__":
    import sys
    sys.exit(run_tests())
