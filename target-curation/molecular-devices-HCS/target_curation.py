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
    """Classify why a well did or did not reach the requested acquisition count.

    `low_cells` is the ordinary shortfall - the well does not hold enough cells that can
    be imaged whole with their margin. `constrained` means it held enough but they could
    not all be acquired; since placement now accepts overlap rather than dropping a cell
    the sample chose, what remains is a well whose sampling grid had fewer occupied
    frames than the count asked for - the cells are there, but not spread over enough of
    the well for an area-uniform draw to reach that many."""
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


def well_of(site):
    """The physical well and facet a site belongs to: (row, column, z, time). The four
    fields of one well share this key; z and time stay separate because they are
    separate images of that well, not more of the same one."""
    row, col, _field, z, time = parse_site(site)
    return row, col, z, time


def field_spans(rows, pixel_um):
    """Map (row, col, field, z, time) to that field's (x, y, width, height) in
    well-local pixels - where the field sits within its well, and how much of the well
    it images.

    Stage positions are recorded in micrometres and image sizes in pixels, so the
    origin is converted with the overview's own pixel size; the well's frame is
    anchored at its top-left-most field. The size comes from the acquisition's own
    `ImageSizeXPx`/`ImageSizeYPx` rather than the camera's sensor region, which is
    recorded unbinned and says nothing about a stitched montage."""
    placed = {}
    for row in rows:
        try:
            key = (int(row["Row"]), int(row["Column"]), int(row["Field"]),
                   int(row["ZIndex"]), int(row["Timepoint"]))
            placed.setdefault(key, (float(row["PositionXUm"]), float(row["PositionYUm"]),
                                    float(row["ImageSizeXPx"]), float(row["ImageSizeYPx"])))
        except (KeyError, ValueError):
            pass
    anchors = {}
    for key in sorted(placed):
        x, y = placed[key][0], placed[key][1]
        well = (key[0], key[1], key[3], key[4])
        at = anchors.get(well)
        anchors[well] = (x, y) if at is None else (min(at[0], x), min(at[1], y))
    spans = {}
    for key in sorted(placed):
        x, y, width, height = placed[key]
        anchor_x, anchor_y = anchors[(key[0], key[1], key[3], key[4])]
        spans[key] = ((x - anchor_x) / pixel_um, (y - anchor_y) / pixel_um, width, height)
    return spans


def read_field_spans(results_dir, pixel_um=None):
    """Measured field spans for this dataset, or an empty mapping when the acquisition
    metadata is not beside the analysis (a stitched montage has none, and needs none)."""
    path = _acquisition_metadata_path(results_dir)
    pixel_um = pixel_um or overview_pixel_size(overview_from_results(results_dir))
    if not path or not pixel_um:
        return {}
    return field_spans(read_csv(path)[1], pixel_um)


def owning_field(x, y, spans, fields, source):
    """Which field owns the point (x, y) of the well, so that ground imaged by more
    than one overlapping field is counted once.

    The owner is the field whose centre is nearest among those whose image actually
    contains the point, ties going to the lowest field number. Containment stops a
    field owning ground it never imaged; nearest-centre puts the boundary on the
    midline between neighbours - for a regular grid exactly half the overlap off each
    shared edge - so every point is owned by the field that sees it furthest from its
    own edge, where the segmentation is most trustworthy. Ownership is decided by
    position, never by matching two detections, so it is unaffected by the two fields
    disagreeing about a cell's outline. `source` is returned when no field contains the
    point, so a detection can never be orphaned by the field geometry."""
    best = None
    for field in fields:
        span = spans.get(field)
        if span is None:
            continue
        ox, oy, width, height = span
        if not (ox <= x <= ox + width and oy <= y <= oy + height):
            continue
        reach = (x - (ox + width / 2.0)) ** 2 + (y - (oy + height / 2.0)) ** 2
        if best is None or (reach, field) < best:
            best = (reach, field)
    return source if best is None else best[1]


def pool_well(by_field, spans):
    """One deduplicated population for a physical well, in well coordinates.

    `by_field` maps a field number to that field's boxes in its own image coordinates.
    Every box is placed in the well and kept only by the field that owns where its
    centre lies, which removes the copies overlapping fields make of the same ground.
    What ownership cannot remove is a cell the two fields segmented slightly
    differently, whose two centres straddle the boundary and are therefore owned by
    different fields; those are collapsed by dropping the higher-numbered field's copy
    when two survivors from different fields sit closer together than the smaller
    one's own radius - a separation two distinct cells cannot have.

    Returns (boxes, sources): boxes in well coordinates, and sources[i] the (field,
    box) it came from, so every result can be reported in that field's own frame
    without ever undoing the translation and losing bits to rounding."""
    fields = sorted(by_field)
    kept = []
    for field in fields:
        span = spans.get(field)
        ox, oy = (span[0], span[1]) if span else (0.0, 0.0)
        for box in by_field[field]:
            moved = (box[0] + ox, box[1] + oy, box[2], box[3])
            if owning_field(moved[0] + moved[2] / 2.0, moved[1] + moved[3] / 2.0,
                            spans, fields, field) == field:
                kept.append((moved, field, box))

    if len(fields) < 2:                    # a stitched montage: nothing to deduplicate
        return [k[0] for k in kept], [(k[1], k[2]) for k in kept]

    # Twins lie within one cell radius of each other, so each box is compared only
    # against those in its own and the neighbouring buckets of a grid coarser than any
    # cell. Comparing every pair would be quadratic in the well's whole population,
    # which a well of tens of thousands of nuclei cannot afford under Jython.
    radius = max([cell_radius(m) for m, _, _ in kept] or [0.0])
    reach = radius * 2.0 if radius > 0 else 1.0
    centres = [(m[0] + m[2] / 2.0, m[1] + m[3] / 2.0) for m, _, _ in kept]
    buckets = {}
    for index, (x, y) in enumerate(centres):
        buckets.setdefault((int(x // reach), int(y // reach)), []).append(index)

    boxes, sources = [], []
    for index, (moved, field, box) in enumerate(kept):
        x, y = centres[index]
        col, row = int(x // reach), int(y // reach)
        twin = False
        for dcol in (-1, 0, 1):
            for drow in (-1, 0, 1):
                for other in buckets.get((col + dcol, row + drow), []):
                    if kept[other][1] >= field:
                        continue
                    ox, oy = centres[other]
                    if ((x - ox) ** 2 + (y - oy) ** 2 <
                            min(cell_radius(moved), cell_radius(kept[other][0])) ** 2):
                        twin = True
        if not twin:
            boxes.append(moved)
            sources.append((field, box))
    return boxes, sources


def field_overlap(spans):
    """How much neighbouring fields of one well overlap, as a fraction of the field
    along each axis. Either value is None when the well is one field wide on that axis;
    a negative value means the fields do not meet and part of the well was never
    imaged. The pitch is the smallest separation between two field origins, so an
    irregular layout is described by its closest pair rather than assumed regular."""
    def along(origins, sizes):
        size = min(sizes)
        steps = sorted(set(abs(a - b) for a in origins for b in origins))
        steps = [s for s in steps if s > size * 1e-6]
        return (size - steps[0]) / size if steps and size > 0 else None

    placed = list(spans.values())
    if len(placed) < 2:
        return None, None
    return (along([p[0] for p in placed], [p[2] for p in placed]),
            along([p[1] for p in placed], [p[3] for p in placed]))


def field_geometry_report(spans):
    """What the acquisition metadata says about the layout, as lines for the log.

    The overlap reported here is the figure the operator set up in MetaXpress,
    recovered from the geometry actually recorded, so a plate whose stage positions or
    overview pixel size are wrong shows as an implausible number instead of being
    curated silently. Pooling itself never needs this: it asks only which fields
    contain a point, which is right for any overlap, none at all, or a gap. But nothing
    else would ever say the recorded layout is not the one on the microscope."""
    by_well = {}
    for key, span in spans.items():
        by_well.setdefault((key[0], key[1], key[3], key[4]), {})[key[2]] = span
    if not by_well:
        return [u"field geometry: none recorded (a stitched montage needs none)"]

    widest = max(len(w) for w in by_well.values())
    lines = [u"field geometry: wells=%d fields=%d max-per-well=%d"
             % (len(by_well), len(spans), widest)]
    if widest < 2:
        return lines + [u"  one field per well - nothing to pool"]

    for name, axis in ((u"x", 0), (u"y", 1)):
        seen = [field_overlap(w)[axis] for w in by_well.values() if len(w) > 1]
        seen = [v * 100.0 for v in seen if v is not None]
        if not seen:
            continue
        lo, hi = min(seen), max(seen)
        lines.append(u"  overlap %s: %s of the field"
                     % (name, u"%.1f%%" % lo if hi - lo < 0.5
                        else u"%.1f-%.1f%%" % (lo, hi)))
        if lo < 0:
            lines.append(u"  WARNING: fields do not meet on %s - part of the well was "
                         u"never imaged, and cells there cannot be sampled" % name)
        elif hi > 50.0:
            lines.append(u"  WARNING: over half the field overlaps on %s - check the "
                         u"overview pixel size and the recorded stage positions" % name)
    return lines


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
    a seed.

    What comes back is the sample, and nothing downstream may reconsider it. Whether a
    cell is easy to place is a function of how close its neighbours are, so choosing
    between cells on that basis - oversampling and keeping whichever placed best, or
    quietly substituting a neighbour for an awkward one - would make inclusion depend on
    local density, which is exactly the bias this grid exists to remove."""
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

    # The grid is an unbounded lattice of `step_x` by `step_y` frames shifted by the
    # random offset, clipped by the scanned area - so the frames at the far edge are
    # partial, exactly as the frames at the near edge are. Clamping the index instead
    # would fold that trailing partial frame into the last full one, leaving one frame
    # (1 + off) steps wide against a first frame of (1 - off): up to a threefold
    # per-area bias that systematically under-samples the right and bottom of the well.
    frames = {}                                                # (col, row) -> cell indices in that frame
    for i in elig:
        col = int((centre_x[i] - x0) / step_x + off_x)
        row = int((centre_y[i] - y0) / step_y + off_y)
        frames.setdefault((col, row), []).append(i)

    chosen = []                                                # one random cell per occupied frame
    for key in sorted(frames):
        members = frames[key]
        chosen.append(members[int(draw() * len(members))])

    # More frames than n can be occupied: rounding the grid dimensions can already
    # exceed n (n=5 on a roughly square well becomes a 2x3 grid), and the random
    # offset adds a partial frame at each far edge. Keep the spatially uniform frame
    # sample, but randomly discard the surplus rather than returning more than n or
    # truncating in coordinate order - a random subset of a uniform sample is uniform.
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
    first) as long as that intersection survives.

    Returns (group, window) where window is (x_lo, x_hi, y_lo, y_hi). Every point in
    it frames the whole group, so a cell is only ever listed in a FOV that truly
    covers it, and the caller is free to choose where inside it the FOV sits."""
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

    return group, (x_lo, x_hi, y_lo, y_hi)


def _clear_of(placed, x, y, fov):
    """True if a FOV centred at (x, y) images no ground an already placed one covers:
    two equal FOV squares overlap exactly when their centres are closer than `fov` on
    BOTH axes. The tolerance absorbs the rounding of a centre computed to sit exactly
    one FOV width from another."""
    return all(abs(x - t["cx"]) >= fov - 1e-9 or abs(y - t["cy"]) >= fov - 1e-9
               for t in placed)


def _candidates(lo, hi, centres, fov):
    """The points on one axis worth testing for a FOV centre inside [lo, hi]: the ends,
    the midpoint, and the lines one FOV width from an already placed centre. Whether a
    FOV clears the placed ones, and by how much it fails to, both change direction only
    at those lines, so testing them is exact rather than a sampling of the window."""
    out = [lo, hi, (lo + hi) / 2.0]
    for c in centres:
        out += [c - fov, c + fov]
    return sorted(set(v for v in out if lo <= v <= hi))


def _clear_point(window, placed, fov):
    """Where to put a FOV inside `window` so it overlaps none already placed, or None.

    Every point of the window frames the whole group, and each placed FOV forbids an
    open square around its centre, so the free part of the window is bounded by the
    window's own edges and by the lines one FOV width from a placed centre. If a clear
    point exists at all, one exists where those lines cross. Prefer the midpoint - it
    holds the group furthest from the FOV edge - and otherwise take the clear crossing
    nearest to it, so the window is declared covered only when it genuinely is rather
    than merely because the midpoint happened to be taken."""
    x_lo, x_hi, y_lo, y_hi = window
    mid_x, mid_y = (x_lo + x_hi) / 2.0, (y_lo + y_hi) / 2.0
    if _clear_of(placed, mid_x, mid_y, fov):
        return mid_x, mid_y

    best = None
    for x in _candidates(x_lo, x_hi, [t["cx"] for t in placed], fov):
        for y in _candidates(y_lo, y_hi, [t["cy"] for t in placed], fov):
            if _clear_of(placed, x, y, fov):
                reach = (x - mid_x) ** 2 + (y - mid_y) ** 2
                if best is None or reach < best[0]:
                    best = (reach, x, y)
    return (best[1], best[2]) if best else None


def _overlap_area(placed, x, y, fov):
    """How much ground a FOV centred at (x, y) would image a second time, summed over
    the FOVs already placed. Zero exactly when it clears them all."""
    total = 0.0
    for t in placed:
        wide = fov - abs(x - t["cx"])
        high = fov - abs(y - t["cy"])
        if wide > 0.0 and high > 0.0:
            total += wide * high
    return total


def _least_overlap_point(window, placed, fov):
    """Where to put a FOV inside `window` when no point in it clears the ones already
    placed: the point that images the least ground twice, ties going to the one nearest
    the midpoint so the group still sits as far from the FOV edge as it can.

    Re-imaging ground is a real cost - more acquisitions, more data, and the same cells
    photographed twice - so it is only ever paid where the caller has established there
    is no alternative, and then at the smallest amount the window allows."""
    x_lo, x_hi, y_lo, y_hi = window
    mid_x, mid_y = (x_lo + x_hi) / 2.0, (y_lo + y_hi) / 2.0
    best = None
    for x in _candidates(x_lo, x_hi, [t["cx"] for t in placed], fov):
        for y in _candidates(y_lo, y_hi, [t["cy"] for t in placed], fov):
            cost = (_overlap_area(placed, x, y, fov), (x - mid_x) ** 2 + (y - mid_y) ** 2)
            if best is None or cost < best[0]:
                best = (cost, x, y)
    return (best[1], best[2]) if best else (mid_x, mid_y)


def _place_disjoint(sampled, centre_x, centre_y, half_window, fov):
    """Place FOVs over the sampled cells by greedy first-fit, disjoint wherever that is
    possible at all.

    Walk the sampled cells; for each not-yet-covered one, build the FOV that frames it
    and any further sampled cells that also fit (several per FOV is fine and free), then
    put that FOV wherever in its feasible window it images no ground already taken.

    When no such point exists, the FOV is placed anyway, at the point of the window that
    images the least ground twice. Overlap is a real cost - more acquisitions, more
    data, the same cells photographed twice - so it is minimised, but it is the cost
    that gets paid, because the alternative is to lose a cell the sample chose. The
    sample is the scientific object here: which cells it names is fixed before placement
    begins and is never revisited, since preferring cells that happen to place well
    would make inclusion depend on how close a cell's neighbours are, and the count is
    meant to fall short only when the well runs out of cells.

    Greedy first-fit, not a provably minimal cover, but the spread sample makes the two
    coincide in practice. Returns [ {cx, cy, covered: [cell_index, ...]} ]."""
    placed, covered = [], set()
    for cell in sampled:
        if cell in covered:
            continue

        remaining = [i for i in sampled if i not in covered]
        group, window = _grow(cell, remaining, centre_x, centre_y, half_window)
        point = _clear_point(window, placed, fov) or _least_overlap_point(window, placed, fov)

        placed.append({"cx": point[0], "cy": point[1], "covered": sorted(group)})
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
                         neighbourhood=2.0, stage_margin=0.05, run_note="", base=None,
                         progress=None, spans=None):
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
                 neighbourhood=neighbourhood, stage_margin=stage_margin, base=base,
                 progress=progress, spans=spans)
    base, output, header = res["base"], res["output"], res["header"]

    changes = [u"ZMB MD HCS target curation changes",
               u"Run\t" + run_note,
               u"Selected object\t" + base,
               u"Acquisition CSV template\t" + output,
               u"Gate\t" + gate_spec,
               u"Cells per well\t%d" % sample_size,
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
           neighbourhood=2.0, stage_margin=0.05, out_csv=None, base=None, progress=None,
           spans=None):
    """Curate every well (gate -> pool the fields -> SURS -> disjoint FOVs). `base` is
    the object class to gate (default: the first channel). Acquisition rows use the
    highest-T channel's schema, independently of the selected object.

    The WELL is the sampling unit, not the acquisition field. A well imaged as several
    overlapping fields arrives as one CSV per field, so `spans` (from
    `read_field_spans`) is used to place them in a common frame, keep each piece of
    ground for the one field that owns it, and sample the well once. Without it a
    four-field well would ask for `sample_size` cells four times over, count the cells
    in the overlap band twice, and never notice that two of its FOVs image the same
    spot. A well of one stitched montage needs no spans and takes exactly the path it
    always did. The SURS grid spans the whole scanned area (all base objects), so the
    sample is spread over the ground the overview imaged, not just where the positives
    landed. The per-well seed is `seed + well_index`, so wells are reproducible yet
    decorrelated (SplitMix64 makes consecutive seeds independent).

    Results are reported per site, in that field's own coordinates: the caller writes
    one CSV per field and the reports draw one panel per field, and both would
    otherwise have to undo the translation. Returns a dict with those per-site results
    and generated FOV rows; writes the curated CSV if out_csv is given."""
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

    # Wells are enumerated in the order their first site appears in the sorted site
    # list rather than by re-sorting the parsed keys, because the per-well seed is
    # `seed + well_index`: a dataset of one field per well must keep the exact index
    # sequence it had when a site was a well, or every sample it ever produced moves.
    well_order, sites_of_well = [], {}
    for site in sorted(files_by_site):
        key = well_of(site)
        if key not in sites_of_well:
            well_order.append(key)
            sites_of_well[key] = []
        sites_of_well[key].append(site)

    header, wells, fov_rows = None, [], []
    half = _usable_half(fov_px, stage_margin)
    for well_index, key in enumerate(well_order):
        if progress:
            progress(well_index + 1, len(well_order))  # for the caller's status/heartbeat

        base_by_field, selected_by_field, templates = {}, {}, {}
        for site in sites_of_well[key]:
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

            field = parse_site(site)[2]
            base_by_field[field] = base_boxes
            selected_by_field[field] = gate(base_boxes, gate_boxes, require)
            templates[field] = (site, output_rows[0] if output_rows else {})
        if not base_by_field:
            continue

        fields = sorted(base_by_field)
        well_spans = {}
        for field in fields:
            span = (spans or {}).get((key[0], key[1], field, key[2], key[3]))
            if span is not None:
                well_spans[field] = span
        if len(fields) > 1 and len(well_spans) < len(fields):
            raise ValueError(
                "Well R%d-C%d (z%d t%d) was acquired as %d fields but the acquisition "
                "metadata gives no geometry for all of them, so they cannot be placed in "
                "one frame. Curating them as independent wells would sample each field "
                "separately and acquire the overlap twice." % (key + (len(fields),)))

        base_pooled, base_sources = pool_well(base_by_field, well_spans)
        selected, sel_sources = pool_well(selected_by_field, well_spans)

        acquired, tiles, n_eligible = select_and_place(
            selected, fov_px, stage_margin, neighbourhood, sample_size,
            seed + well_index, search_boxes=base_pooled)

        # positives imaged whole by the FOVs = the sampled targets PLUS any other
        # positive cell that happens to fall entirely inside a field (free extra
        # observations, especially at low mag where the field is large)
        acquired_idx = sorted(set(i for t in tiles for i in t["covered"]))
        captured_idx = [i for i, b in enumerate(selected)
                        if any(_fully_inside(b, t["cx"], t["cy"], half) for t in tiles)]
        extra_idx = [i for i in captured_idx if selected[i] not in acquired]

        # Hand every result back to the field it came from, in that field's own
        # coordinates. A FOV goes to the field that owns the ground under its centre.
        site_of_field = dict((f, templates[f][0]) for f in fields)
        parts = dict((templates[f][0], {"base": [], "selected": [], "acquired": [],
                                        "extra_boxes": [], "captured": 0, "eligible": 0,
                                        "tiles": [], "fov_rows": []}) for f in fields)
        for field, local in base_sources:
            parts[site_of_field[field]]["base"].append(local)
        for field, local in sel_sources:
            parts[site_of_field[field]]["selected"].append(local)
        for i in acquired_idx:
            parts[site_of_field[sel_sources[i][0]]]["acquired"].append(sel_sources[i][1])
        for i in captured_idx:
            parts[site_of_field[sel_sources[i][0]]]["captured"] += 1
        for i in extra_idx:
            parts[site_of_field[sel_sources[i][0]]]["extra_boxes"].append(sel_sources[i][1])
        for i in eligible(selected, fov_px, stage_margin, neighbourhood):
            parts[site_of_field[sel_sources[i][0]]]["eligible"] += 1

        start = len(fov_rows)
        for i, tile in enumerate(tiles):
            home = owning_field(tile["cx"], tile["cy"], well_spans, fields, fields[0])
            origin = well_spans.get(home, (0.0, 0.0))
            local = dict(tile)
            local["cx"], local["cy"] = tile["cx"] - origin[0], tile["cy"] - origin[1]
            row = fov_row(templates[home][1], output_t, local["cx"], local["cy"],
                          "FOV_%d" % (start + i + 1))
            part = parts[site_of_field[home]]
            part["tiles"].append(local)
            part["fov_rows"].append(row)
            fov_rows.append(row)

        # The status describes the well - whether IT yielded the cells asked of it - so
        # every field of that well reports it, and none of them re-decides it alone.
        status = site_status(n_eligible, len(acquired), sample_size)
        for f in fields:
            part = parts[site_of_field[f]]
            part.update({"site": site_of_field[f], "well": key, "status": status,
                         "sample_short": status != "full",
                         "extra": len(part["extra_boxes"])})
            wells.append(part)

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
    spans = read_field_spans(res["target_dir"], overview_pixel_size(overview_desc))
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
            measured = spans.get(key)
            site_origins[site] = ((measured[0], measured[1]) if measured is not None else
                                  (slot["field_col"] * max_x, slot["field_row"] * max_y))
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

    # Where each acquisition field sits in its well. A stitched montage has no such
    # metadata and needs none; a well of several fields cannot be curated without it,
    # and curate() refuses rather than sampling each field as though it were a well.
    spans = read_field_spans(results_path, overview_pixel_size(overview_desc))
    for line in field_geometry_report(spans):
        IJ.log("  " + line)

    IJ.log("  working - curating wells (reading CSVs, sampling, placing FOVs)...")
    IJ.showStatus("MD HCS curation: curating wells...")
    res, cur_dir, audit_dir = write_curated_output(
        results_path, gate_spec, fov_px, sample_size, seed, neighbourhood, stage_margin,
        run_note=str(Date()), base=base, progress=tick, spans=spans)

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
    # Status belongs to the well, and every field of that well carries it, so the
    # shortfalls are counted over distinct wells rather than once per field.
    by_well = dict((w["well"], w["status"]) for w in res["wells"])
    under = sum(1 for s in by_well.values() if s != "full")
    constrained = sum(1 for s in by_well.values() if s == "constrained")
    IJ.log("  wells=%d sites=%d acquired=%d extra-whole=%d FOVs=%d under-sampled=%d constrained=%d"
           % (len(by_well), len(res["wells"]), acquired, extra, fovs, under, constrained))
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
# The acceptance datasets live on the facility share, not in the repo, so every block
# that uses them skips when they are absent. Nico is the single-field control (one
# stitched montage per well); Babette is the four-field case (2x2, 10% overlap). Nico
# points at the untouched `- Copy`: the sibling folder without the suffix has already
# been curated in place, so its TargetData/ holds generated FOV rows rather than IN
# Carta's original objects.
# The standing single-field control: six wells, one stitched montage each, one signal.
# Read TargetData_original, never TargetData - this folder has been curated in place,
# so its TargetData/ holds generated FOV rows rather than IN Carta's own objects.
CONTROL_RESULTS = (r"Z:\transfer\Thom\10306\FIJI_Target_Curation_Test"
                   r"\0306_H2BmCherry_OVWF10x_TXRED"
                   r"\10306_H2BmCherry_OVWF10x_TXRED_20260429_151322\experiment\Results"
                   r"\curation test_2026-Apr-29-15-23-15-269\TargetData_original")
# A plate acquired as several fields per well, to run the well-level invariants
# against real geometry. The Babette plate this was written for is no longer on the
# share; point this at the next multi-field plate that is.
MULTIFIELD_RESULTS = (r"Z:\transfer\Thom\Babette MD"
                      r"\26.19 NPTX2 ASO IF test_20260722_104349\experiment"
                      r"\Results\TriplePositive_Thom_2026-Jul-22-12-20-15-690")
MULTIFIELD_GATE, MULTIFIELD_FOV = "Green:yes; Red:yes", 1536.0
# The Nico plate the montage goldens below were measured on, also gone from the share.
NICO_RESULTS = (r"Z:\transfer\Thom\Nico MD"
                r"\NB26-15_Overview10x_DAPI-mScarlet_20260713_152811\experiment_montage"
                r"\Results\mScarlet Cells_2026-Jul-13-17-23-33-077 - Copy")


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
    GOLDEN_TILES = [(1800.5, 1800.5), (1800.5, 1000.5), (0.5, 1800.5),
                    (200.5, 1000.5), (400.5, 400.5)]
    # Nico plate, gate `mScarlet cells:yes`, FOV 384 px (60x target from a 10x overview),
    # 5 cells/well, seed 42 - the parameters of the production run.
    GOLDEN_POSITIVES, GOLDEN_ELIGIBLE = 4253, 4247
    GOLDEN_ACQUIRED, GOLDEN_FOVS, GOLDEN_CAPTURED, GOLDEN_EXTRA = 167, 166, 282, 115
    NICO_FOV = 384.0

    # The single-field control plate, FOV 384 px, 5 cells/well, seed 42. Measured on
    # the code as it stood before well-level pooling and confirmed unmoved by it, so
    # these numbers pin the montage path itself, not merely today's output.
    CONTROL_FOV = 384.0
    CONTROL_POSITIVES, CONTROL_ELIGIBLE = 1637, 1637
    CONTROL_ACQUIRED, CONTROL_FOVS, CONTROL_CAPTURED, CONTROL_EXTRA = 30, 28, 240, 210

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
    def meta_row(field, x, y):
        return {"Row": "3", "Column": "2", "Field": str(field), "ZIndex": "0",
                "Timepoint": "0", "PositionXUm": str(x), "PositionYUm": str(y),
                "ImageSizeXPx": "100", "ImageSizeYPx": "80"}

    metadata_rows = [meta_row(0, 10, 20), meta_row(1, 10, 30), meta_row(2, 25, 30)]
    check("stage metadata becomes well-local pixel spans",
          field_spans(metadata_rows, 0.5),
          {(3, 2, 0, 0, 0): (0.0, 0.0, 100.0, 80.0),
           (3, 2, 1, 0, 0): (0.0, 20.0, 100.0, 80.0),
           (3, 2, 2, 0, 0): (30.0, 20.0, 100.0, 80.0)})
    check("a field without a recorded image size is left out",
          field_spans([{"Row": "3", "Column": "2", "Field": "0", "ZIndex": "0",
                        "Timepoint": "0", "PositionXUm": "10", "PositionYUm": "20"}], 0.5), {})
    check("the well key drops the field but keeps the z/time facet",
          well_of("R2-C3-F1-Z0-T4"), (2, 3, 0, 4))

    faceted = build_overview_layout(
        [{"site": "R3-C2-F0-Z0-T0"}, {"site": "R3-C2-F0-Z1-T0"},
         {"site": "R3-C2-F0-Z0-T1"}])
    check("Z/T variants are separate facets", sorted(faceted), [(0, 0), (0, 1), (1, 0)])

    print("pooling overlapping fields into one well")

    def grid(cols, rows, size, pitch_x, pitch_y=None):
        """Fields of one well on a regular grid, numbered row-major."""
        pitch_y = pitch_x if pitch_y is None else pitch_y
        return dict((r * cols + c, (c * pitch_x * 1.0, r * pitch_y * 1.0,
                                    size * 1.0, size * 1.0))
                    for r in range(rows) for c in range(cols))

    def lattice(spans, step):
        """Every lattice point the fields image, with its owner and the fields that
        actually contain it - the sampled stand-in for 'the whole imaged area'."""
        fields = sorted(spans)
        far_x = max(spans[f][0] + spans[f][2] for f in fields)
        far_y = max(spans[f][1] + spans[f][3] for f in fields)
        out = {}
        for i in range(int(far_x / step) + 1):
            for j in range(int(far_y / step) + 1):
                x, y = i * step * 1.0, j * step * 1.0
                inside = [f for f in fields
                          if spans[f][0] <= x <= spans[f][0] + spans[f][2]
                          and spans[f][1] <= y <= spans[f][1] + spans[f][3]]
                if inside:
                    out[(x, y)] = (owning_field(x, y, spans, fields, -1), inside)
        return out

    def trimmed(x, y, cols, rows, size, pitch_x, pitch_y):
        """The owner under the textbook rule for a regular grid - half of every shared
        overlap trimmed off each side - which the nearest-centre partition must
        reproduce cell for cell."""
        col = sum(1 for c in range(cols - 1) if x > c * pitch_x + (size + pitch_x) / 2.0)
        row = sum(1 for r in range(rows - 1) if y > r * pitch_y + (size + pitch_y) / 2.0)
        return row * cols + col

    for label, cols, rows, size, pitch_x, pitch_y in [
            ("one field", 1, 1, 100, 100, 100),
            ("2x2 at 10% overlap", 2, 2, 100, 90, 90),
            ("3x3 at 10% overlap", 3, 3, 100, 90, 90),
            ("unequal overlap per axis", 2, 2, 100, 90, 80),
            ("abutting fields, no overlap", 2, 2, 100, 100, 100)]:
        spans_case = grid(cols, rows, size, pitch_x, pitch_y)
        points = lattice(spans_case, 5)
        orphans = sum(1 for own, inside in points.values() if own not in inside)
        strays = sum(1 for (x, y), (own, _) in points.items()
                     if own != trimmed(x, y, cols, rows, size, pitch_x, pitch_y))
        check("%s: every imaged point is owned by a field that imaged it" % label, orphans, 0)
        check("%s: ownership is exactly a half-overlap trim" % label, strays, 0)

    g22 = grid(2, 2, 100, 90)
    check("a point exactly on the midline goes to the lower field number",
          owning_field(95.0, 50.0, g22, sorted(g22), -1), 0)
    check("the four-way corner point goes to the lowest field number",
          owning_field(95.0, 95.0, g22, sorted(g22), -1), 0)
    gapped = grid(2, 1, 100, 110)
    check("a point in a gap between fields keeps the field it was detected in",
          owning_field(105.0, 50.0, gapped, sorted(gapped), 7), 7)

    seam = {0: (0.0, 0.0, 100.0, 100.0), 1: (90.0, 0.0, 100.0, 100.0)}
    check("one field per well is passed through untouched",
          pool_well({0: [bx(10, 10, 4, 4), bx(50, 50, 6, 6)]}, {0: (0.0, 0.0, 100.0, 100.0)}),
          ([bx(10, 10, 4, 4), bx(50, 50, 6, 6)], [(0, bx(10, 10, 4, 4)), (0, bx(50, 50, 6, 6))]))
    check("a pooled box carries its well position and its field-local original",
          pool_well({1: [bx(10, 20, 4, 4)]}, {1: (200.0, 300.0, 100.0, 100.0)}),
          ([bx(210, 320, 4, 4)], [(1, bx(10, 20, 4, 4))]))
    check("a cell both fields agree on is pooled once",
          len(pool_well({0: [bx(93, 40, 4, 4)], 1: [bx(3, 40, 4, 4)]}, seam)[0]), 1)
    check("one cell the two fields outline differently collapses to the lower field",
          pool_well({0: [bx(92, 40, 4, 4)], 1: [bx(4, 40, 4, 4)]}, seam),
          ([bx(92, 40, 4, 4)], [(0, bx(92, 40, 4, 4))]))
    check("two distinct cells either side of the seam are both kept",
          len(pool_well({0: [bx(90, 40, 4, 4)], 1: [bx(8, 40, 4, 4)]}, seam)[0]), 2)
    check("distinct cells in different fields are both kept",
          len(pool_well({0: [bx(20, 40, 4, 4)], 1: [bx(60, 40, 4, 4)]}, seam)[0]), 2)

    check("the overlap the acquisition recorded is recovered",
          field_overlap(grid(2, 2, 100, 90)), (0.1, 0.1))
    check("an overlap that differs per axis is reported per axis",
          field_overlap(grid(2, 2, 100, 90, 80)), (0.1, 0.2))
    check("fields that merely abut overlap by nothing",
          field_overlap(grid(2, 2, 100, 100)), (0.0, 0.0))
    check("a gap between fields reads as negative overlap",
          field_overlap(grid(2, 1, 100, 110)), (-0.1, None))
    check("one field has no overlap to report",
          field_overlap(grid(1, 1, 100, 100)), (None, None))

    def as_well(spans_of_fields):
        return dict(((2, 3, f, 0, 0), s) for f, s in spans_of_fields.items())

    check("the report names the layout and the overlap it measured",
          field_geometry_report(as_well(grid(2, 2, 100, 90))),
          [u"field geometry: wells=1 fields=4 max-per-well=4",
           u"  overlap x: 10.0% of the field",
           u"  overlap y: 10.0% of the field"])
    check("a montage needs no geometry and the report says so",
          field_geometry_report({}),
          [u"field geometry: none recorded (a stitched montage needs none)"])
    check("one field per well is called out as nothing to pool",
          field_geometry_report(as_well(grid(1, 1, 100, 100))),
          [u"field geometry: wells=1 fields=1 max-per-well=1",
           u"  one field per well - nothing to pool"])
    check("a gap in the layout is warned about, not just reported",
          [l for l in field_geometry_report(as_well(grid(2, 1, 100, 110)))
           if u"WARNING" in l],
          [u"  WARNING: fields do not meet on x - part of the well was never imaged, "
           u"and cells there cannot be sampled"])
    check("an implausible overlap points at the pixel size and stage positions",
          [l for l in field_geometry_report(as_well(grid(2, 1, 100, 40)))
           if u"WARNING" in l],
          [u"  WARNING: over half the field overlaps on x - check the overview pixel "
           u"size and the recorded stage positions"])

    print("SURS + disjoint placement")
    check("oversized cell not eligible", eligible([bx(0, 0, 500, 500)], 100, 0.0, 0.0), [])
    _, tiles, _ = select_and_place([bx(0, 0, 10, 10), bx(20, 0, 10, 10)], 100, 0.0, 0.0, 5, 1)
    check("close pair -> 1 FOV, 2 cells", (len(tiles), len(tiles[0]["covered"])), (1, 2))
    _, tiles, _ = select_and_place([bx(0, 0, 10, 10), bx(1000, 0, 10, 10)], 100, 0.0, 0.0, 5, 1)
    check("far pair -> 2 disjoint FOVs", (len(tiles), disjoint(tiles, 100)), (2, True))
    # genuinely infeasible, not merely inconvenient: sample all three; A and B (70
    # apart, usable 60) cannot share a FOV, and B's whole feasible window lies within
    # one FOV width of A's, so no placement of B clears A. All three cells are sampled,
    # so all three are acquired: B's FOV is placed anyway, overlapping A's, because the
    # count is only allowed to fall short when the cells run out. Nothing here can be
    # swapped for it - each cell is its own frame - so overlap is the only way.
    got, tiles, _ = select_and_place([bx(0, 0, 4, 4), bx(70, 0, 4, 4), bx(1000, 0, 4, 4)],
                                     100, 0.2, 0.0, 3, 5)
    clashes = sum(1 for i, a in enumerate(tiles) for b in tiles[i + 1:]
                  if abs(a["cx"] - b["cx"]) < 100 and abs(a["cy"] - b["cy"]) < 100)
    check("a cell with no clear placement is imaged anyway, not given up",
          (len(tiles), set(b[0] for b in got)), (3, set([0.0, 70.0, 1000.0])))
    check("and it costs exactly one overlapping pair", clashes, 1)
    framed = [t for t in tiles if abs(t["cx"] - 72.0) <= 27.2 and abs(t["cy"] - 2.0) <= 27.2]
    check("the overlapping FOV still frames the cell it was placed for", len(framed), 1)
    check("the overlap taken is the least the window allowed",
          _overlap_area([{"cx": 2.0, "cy": 2.0}], framed[0]["cx"], framed[0]["cy"], 100)
          <= _overlap_area([{"cx": 2.0, "cy": 2.0}], 72.0, 2.0, 100), True)

    # ...and the awkward cell is never quietly traded for an easier neighbour. Cell 1's
    # window lies wholly within one FOV of cell 0 while cell 2's does not, but cell 2
    # was not sampled, so it is cell 1 that gets imaged - with overlap - and cell 2 is
    # left alone. Substituting it would make inclusion depend on how close a cell's
    # neighbours are, which is the density bias the sampling grid exists to remove.
    near = [bx(-20, -20, 40, 40), bx(40, -20, 40, 40), bx(130, -20, 40, 40)]
    near_hw = dict((i, 50.0 - cell_radius(near[i])) for i in range(3))
    near_x = dict((i, near[i][0] + 20.0) for i in range(3))
    near_y = dict((i, near[i][1] + 20.0) for i in range(3))
    kept = _place_disjoint([0, 1], near_x, near_y, near_hw, 100.0)
    check("the sampled cell is imaged, not an easier neighbour of it",
          (len(kept), disjoint(kept, 100.0),
           sorted(i for t in kept for i in t["covered"])), (2, False, [0, 1]))
    # ...whereas a cell whose window merely straddles a placed FOV is shifted clear
    # rather than dropped: centres 88 apart with a window of +/-27 around each, so the
    # midpoint sits 88 from the placed FOV (too close) but x = 102 is both a full FOV
    # width away and still inside the window. Before, this cell was given up.
    got2, tiles2, _ = select_and_place([bx(0, 0, 4, 4), bx(88, 0, 4, 4)], 100, 0.2, 0.0, 2, 5)
    check("cell shifted clear instead of dropped",
          (len(tiles2), disjoint(tiles2, 100), len(got2)), (2, True, 2))
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
    # Per-area uniformity across the WHOLE extent, not just within a clump: on an even
    # lattice each half of the scanned area must take half the picks. Folding the
    # trailing partial frame into the last full one used to draw ~65% from the near half.
    lattice = [bx(x * 10, y * 10, 1, 1) for x in range(41) for y in range(41)]
    lx = dict((i, lattice[i][0] + 0.5) for i in range(len(lattice)))
    ly = dict((i, lattice[i][1] + 0.5) for i in range(len(lattice)))
    mid_x = (min(lx.values()) + max(lx.values())) / 2.0
    mid_y = (min(ly.values()) + max(ly.values())) / 2.0
    picks = [i for sd in range(400)
             for i in _surs_sample(list(range(len(lattice))), lx, ly, 5, sd, lattice)]
    near_x = sum(1 for i in picks if lx[i] < mid_x)
    near_y = sum(1 for i in picks if ly[i] < mid_y)
    print("    halves took %d / %d (x) and %d / %d (y) of %d picks; equal area -> half"
          % (near_x, len(picks) - near_x, near_y, len(picks) - near_y, len(picks)))
    check("per-area uniform across the whole extent, both axes",
          (abs(2.0 * near_x / len(picks) - 1.0) < 0.10,
           abs(2.0 * near_y / len(picks) - 1.0) < 0.10), (True, True))
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
    # Read a real MetaXpress TIFF header, whichever plate is on the share. The optics
    # are the plate's own, so what is checked is that the fields come back usable - a
    # magnification, a binning and a region the FOV arithmetic can be driven with.
    real_ov = (_overview_tiff_path(MULTIFIELD_RESULTS)
               if os.path.isdir(os.path.join(MULTIFIELD_RESULTS, "TargetData")) else None)
    if real_ov and os.path.isfile(real_ov):
        mag_r, changer_r, binning_r, region_r = overview_scale(read_image_description(real_ov))
        check("real overview TIFF scale read",
              (mag_r > 0, changer_r > 0, binning_r >= 1, region_r > 0), (True,) * 4)
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
    if os.path.isdir(os.path.join(MULTIFIELD_RESULTS, "TargetData")):
        real_desc = overview_from_results(MULTIFIELD_RESULTS)
        check("overview auto-found + scaled from real results",
              overview_scale(real_desc)[0] > 0, True)
        check("overview pixel size read from the same image",
              overview_pixel_size(real_desc) > 0, True)
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

    print("one well, four overlapping fields (synthetic)")
    tmp4 = tempfile.mkdtemp()
    try:
        FIELD, PITCH, FOV4 = 2304.0, 2073.5, 200.0        # 2x2 fields at 10% overlap
        origins4 = {0: (0.0, 0.0), 1: (PITCH, 0.0), 2: (0.0, PITCH), 3: (PITCH, PITCH)}
        spans4 = dict(((2, 3, f, 0, 0), (ox, oy, FIELD, FIELD))
                      for f, (ox, oy) in origins4.items())

        # Six physical cells: one in the interior of each field, one in the vertical
        # seam and one in the four-way corner. The seam cell falls in two fields'
        # images and the corner cell in all four, so the four CSVs hold ten rows.
        cells4 = [(500.0, 500.0), (3000.0, 500.0), (500.0, 3000.0), (3000.0, 3000.0),
                  (2150.0, 500.0), (2150.0, 2150.0)]
        td4 = os.path.join(tmp4, "TargetData")
        os.makedirs(td4)
        rows_written = 0
        for f, (ox, oy) in sorted(origins4.items()):
            local = [(x - ox, y - oy) for x, y in cells4
                     if ox <= x <= ox + FIELD and oy <= y <= oy + FIELD]
            rows_written += len(local)
            dapi4 = ["object_id,well_label,T1$AS_FID_Blob_BoundingBoxX,T1$AS_FID_Blob_BoundingBoxY,"
                     "T1$AS_FID_Blob_BoundingBoxWidth,T1$AS_FID_Blob_BoundingBoxHeight"]
            ms4 = ["object_id,T2$AS_FID_Blob_BoundingBoxX,T2$AS_FID_Blob_BoundingBoxY,"
                   "T2$AS_FID_Blob_BoundingBoxWidth,T2$AS_FID_Blob_BoundingBoxHeight"]
            for i, (x, y) in enumerate(local):
                dapi4.append("%d,B - 3,%g,%g,20,20" % (i, x - 10, y - 10))
                ms4.append("%d,%g,%g,12,12" % (i, x - 6, y - 6))
            for name, lines in (("DAPI", dapi4), ("mScarlet cells", ms4)):
                with io.open(os.path.join(td4, "%s_singleTargetData_R2-C03-F%d-Z0-T0.csv"
                                          % (name, f)), "w", encoding="utf-8", newline="") as fh:
                    fh.write(u"\r\n".join(lines) + u"\r\n")
        check("the four CSVs hold ten rows for six cells", rows_written, 10)

        try:
            curate(tmp4, "mScarlet cells:yes", fov_px=FOV4, sample_size=3)
            refused = False
        except ValueError:
            refused = True
        check("four fields without measured geometry is a hard failure", refused, True)

        res4 = curate(tmp4, "mScarlet cells:yes", fov_px=FOV4, sample_size=3, spans=spans4)
        placed4 = []
        for w in res4["wells"]:
            ox, oy = origins4[parse_site(w["site"])[2]]
            placed4 += [(t["cx"] + ox, t["cy"] + oy) for t in w["tiles"]]
        clashes = sum(1 for i, a in enumerate(placed4) for b in placed4[i + 1:]
                      if abs(a[0] - b[0]) < FOV4 and abs(a[1] - b[1]) < FOV4)
        owns_seam = [w["site"] for w in res4["wells"]
                     for b in w["selected"]
                     if abs(b[0] + b[2] / 2.0 - 2150.0) < 1e-9
                     and abs(b[1] + b[3] / 2.0 - 500.0) < 1e-9]

        check("the well is still reported as its four fields", len(res4["wells"]), 4)
        check("a cell imaged by several fields is curated once",
              sum(len(w["selected"]) for w in res4["wells"]), 6)
        check("the sample is drawn once for the well, not once per field",
              sum(len(w["acquired"]) for w in res4["wells"]), 3)
        check("no two of the well's FOVs overlap in well coordinates", clashes, 0)
        check("a seam cell goes to the field whose centre is nearest",
              owns_seam, ["R2-C03-F0-Z0-T0"])
        check("every field of the well reports the well's own status",
              sorted(set(w["status"] for w in res4["wells"])), ["full"])
        check("captured == acquired + extra across the well",
              (sum(w["captured"] for w in res4["wells"]),
               sum(len(w["acquired"]) + w["extra"] for w in res4["wells"])),
              (3, 3))
    finally:
        shutil.rmtree(tmp4, ignore_errors=True)

    print("real data - single-field montage control (auto-skip)")
    if os.path.isdir(CONTROL_RESULTS):
        res = curate(CONTROL_RESULTS, "Nuclei:yes", fov_px=CONTROL_FOV)
        pos = sum(len(w["selected"]) for w in res["wells"])
        acquired = sum(len(w["acquired"]) for w in res["wells"])
        captured = sum(w["captured"] for w in res["wells"])
        extra = sum(w["extra"] for w in res["wells"])
        fovs = sum(len(w["tiles"]) for w in res["wells"])
        print("    plate: positives=%d acquired=%d extra=%d captured=%d FOVs=%d"
              % (pos, acquired, extra, captured, fovs))
        check("control signals discovered", res["signals"], ["Nuclei"])
        check("control is one field per well",
              sorted(set(parse_site(w["site"])[2] for w in res["wells"])), [0])
        check("control positives (golden)", pos, CONTROL_POSITIVES)
        check("control eligible (golden)",
              sum(w["eligible"] for w in res["wells"]), CONTROL_ELIGIBLE)
        check("control acquired (golden)", acquired, CONTROL_ACQUIRED)
        check("control FOVs (golden)", fovs, CONTROL_FOVS)
        check("control captured (golden)", captured, CONTROL_CAPTURED)
        check("control extra/bonus (golden)", extra, CONTROL_EXTRA)
        check("control captured == acquired + extra", captured, acquired + extra)
    else:
        print("  SKIP (dataset not present)")

    print("real data - the retired Nico plate (auto-skip)")
    if os.path.isdir(os.path.join(NICO_RESULTS, "TargetData")):
        res = curate(NICO_RESULTS, "mScarlet cells:yes", fov_px=NICO_FOV)
        pos = sum(len(w["selected"]) for w in res["wells"])
        elig_n = sum(w["eligible"] for w in res["wells"])
        acquired = sum(len(w["acquired"]) for w in res["wells"])
        captured = sum(w["captured"] for w in res["wells"])
        extra = sum(w["extra"] for w in res["wells"])
        fovs = sum(len(w["tiles"]) for w in res["wells"])
        overlaps = sum(1 for wl in res["wells"] for a in range(len(wl["tiles"]))
                       for b in range(a + 1, len(wl["tiles"]))
                       if abs(wl["tiles"][a]["cx"] - wl["tiles"][b]["cx"]) < NICO_FOV
                       and abs(wl["tiles"][a]["cy"] - wl["tiles"][b]["cy"]) < NICO_FOV)
        print("    plate: positives=%d eligible=%d acquired=%d extra=%d captured=%d FOVs=%d overlaps=%d"
              % (pos, elig_n, acquired, extra, captured, fovs, overlaps))
        check("signals discovered", res["signals"], ["DAPI", "mScarlet cells"])
        check("plate positives (golden)", pos, GOLDEN_POSITIVES)
        check("plate eligible (golden)", elig_n, GOLDEN_ELIGIBLE)
        check("one field per well (the montage control)",
              sorted(set(parse_site(w["site"])[2] for w in res["wells"])), [0])
        check("no overlapping FOVs on the plate", overlaps, 0)
        check("plate acquired (golden)", acquired, GOLDEN_ACQUIRED)
        check("captured == acquired + extra", captured, acquired + extra)
        check("plate extra/bonus (golden)", extra, GOLDEN_EXTRA)
        check("extra_boxes count matches extra", sum(len(w["extra_boxes"]) for w in res["wells"]), GOLDEN_EXTRA)
        check("plate captured (golden)", captured, GOLDEN_CAPTURED)
        check("plate FOVs (golden)", fovs, GOLDEN_FOVS)
    else:
        print("  SKIP (dataset not present)")

    # A real multi-field plate, whenever one is on the share. The counts depend on the
    # plate, so what is asserted here are the invariants the well-level path must hold
    # on any of them - the same three the synthetic four-field well pins exactly.
    print("real data - multi-field plate (auto-skip)")
    if os.path.isdir(os.path.join(MULTIFIELD_RESULTS, "TargetData")):
        spans = read_field_spans(
            MULTIFIELD_RESULTS,
            overview_pixel_size(overview_from_results(MULTIFIELD_RESULTS)))
        res = curate(MULTIFIELD_RESULTS, MULTIFIELD_GATE, fov_px=MULTIFIELD_FOV, spans=spans)

        pooled, fields_per_well = {}, {}
        for w in res["wells"]:
            field = parse_site(w["site"])[2]
            ox, oy = spans[parse_site(w["site"])][:2]
            fields_per_well.setdefault(w["well"], set()).add(field)
            cur = pooled.setdefault(w["well"], {"cells": [], "tiles": []})
            cur["cells"] += [(b[0] + b[2] / 2.0 + ox, b[1] + b[3] / 2.0 + oy, field)
                             for b in w["acquired"]]
            cur["tiles"] += [(t["cx"] + ox, t["cy"] + oy, field) for t in w["tiles"]]

        def cross_field(well, key, near):
            items = well[key]
            return sum(1 for i, a in enumerate(items) for b in items[i + 1:]
                       if a[2] != b[2] and near(a, b))

        duplicates = sum(cross_field(w, "cells", lambda a, b:
                                     (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 < 40.0 ** 2)
                         for w in pooled.values())
        xfield = sum(cross_field(w, "tiles", lambda a, b:
                                 abs(a[0] - b[0]) < MULTIFIELD_FOV
                                 and abs(a[1] - b[1]) < MULTIFIELD_FOV)
                     for w in pooled.values())
        over_target = sum(1 for w in pooled.values() if len(w["cells"]) > 5)
        print("    plate: wells=%d fields/well=%s acquired=%d FOVs=%d"
              % (len(pooled), sorted(set(len(v) for v in fields_per_well.values())),
                 sum(len(w["cells"]) for w in pooled.values()),
                 sum(len(w["tiles"]) for w in pooled.values())))
        short = [w["site"] for w in res["wells"] if w["status"] == "constrained"]
        check("more than one field per well", max(len(v) for v in fields_per_well.values()) > 1, True)
        check("no cell is acquired from two fields of one well", duplicates, 0)
        check("no well exceeds the requested cell count", over_target, 0)
        # FOVs may now overlap where that was the only way to reach the requested count,
        # so what must hold is the reason a well falls short: never the placement.
        check("no well falls short while it still has cells to give", short, [])
        print("    cross-field FOV overlaps (allowed only as a last resort) = %d" % xfield)
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
