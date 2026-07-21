#@ File    (label="Analysis results folder", style="directory") results_dir
#@ Integer (label="Cells per well", value=5, min=1) sample_size
#@ Integer (label="Random seed", value=42) seed
#@ String  (label="Target objective", choices={"2x","4x","10x","20x","40x","60x","60x + 1.5x changer"}, value="60x") target_objective
#@ File    (label="Overview image (reads objective + scale)", style="file") overview_image
#@ Float   (label="Neighbourhood (clear radii around each cell)", value=2.0, min=0) neighbourhood
#@ Float   (label="Stage-accuracy margin (fraction per side)", value=0.05, min=0, max=0.49) stage_margin
#@ File    (label="Output folder", style="directory") out_dir

"""
MD HCS target curation - single-file Fiji macro (Jython 2.7) and CPython module.

Molecular Devices ImageXpress / IN Carta writes one CSV of segmented objects per
signal per well. This picks, per well, a small representative set of cells to
re-image at high magnification, in three steps:

  1. GATE the objects of interest. The base object (the FIRST channel, the nuclei) is
     kept when it is single / double / triple positive - i.e. its box overlaps an
     object of every 'yes' signal and none of the 'no' signals. In Fiji the marker
     classes are discovered from the data and their roles (positive / negative /
     ignore) are chosen in a dialog; the base is always the first channel.

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

OUTPUT: a curated CSV of generated FOV rows (hand-load into MetaXpress), a per-well
summary CSV, and a per-well visual report (PNG: nuclei / positives / sampled / FOVs).
The summary reports, per well, `acquired` (the sampled target cells, always imaged
whole) and `extra` (BONUS positives that also fall entirely inside a field, uncut -
not sampled targets; grows at low mag where fields are large). `captured` = acquired
+ extra = every positive imaged whole.

In Fiji the target FOV size is not typed in: pick the target objective and point at
one overview image. The FOV footprint (in overview-montage pixels) is inferred from
the overview's own metadata - its objective, changer, binning and sensor region -
so it adapts automatically when the overview is taken at a different magnification.

The file is BOTH the Fiji entry (SciJava params above) and its own test suite: run
`python target_curation.py` in CPython to run the tests; open/run it in Fiji to
curate. All maths is stdlib-only and written to be interpreter-independent (a
portable PRNG and int(x+0.5) rounding, because Python's own `random` and `round`
differ between CPython 3 and Jython 2.7), so results are byte-identical in both.

The code below reads top-to-bottom as the pipeline flows: import the data (1) ->
gate it (2) -> work out the FOV size (3) and what fits in a FOV (4) -> draw random
numbers (5) -> sample (6) -> place fields (7) -> combine per well (8) -> write the
output (9) -> curate every well (10) -> the Fiji front end (11), tests (12), dispatch.
"""
from __future__ import division, print_function
import io, os, re

_BB = "AS_FID_Blob_BoundingBox"          # + X / Y / Width / Height
_OWN = re.compile(r"^T(\d+)\$")          # a signal's own-measurement column prefix
_SITE_KEY = "_singleTargetData_"         # filename split: <signal>_singleTargetData_<site>.csv
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


def _has_csv(d):
    return os.path.isdir(d) and any(f.lower().endswith(".csv") for f in os.listdir(d))


def find_target_dir(path):
    """Locate the TargetData folder: `path` if it is one, else a TargetData child,
    else one analysis-folder level down. Preferred over sibling ObjectData/metadata
    CSVs, which also live in the analysis folder but are not the targets."""
    if os.path.basename(os.path.normpath(path)) == "TargetData" and _has_csv(path):
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


def fov_px_for(target_mag, target_changer, mag_ov, changer_ov, binning_ov, region_px):
    """FOV footprint in overview-montage pixels for a target objective. The camera
    reads `region_px` pixels at any magnification; its montage footprint is that
    count scaled by the overview-to-target total-magnification ratio (overview
    binning coarsens the montage, so it enlarges the footprint)."""
    return region_px * (mag_ov * changer_ov) / (target_mag * target_changer * binning_ov)


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
    """~n eligible cell indices by Systematic Uniform Random Sampling - the 2D
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


# --------------------------------------------------------------------------- #
# 10. Curate - run the whole pipeline over every well                         #
# --------------------------------------------------------------------------- #
def curate(results_dir, gate_spec, fov_px=256.0, sample_size=5, seed=42,
           neighbourhood=2.0, stage_margin=0.05, out_csv=None):
    """Curate every well (gate -> SURS -> disjoint FOVs). Base = first signal. The
    SURS grid spans the whole scanned area (all nuclei), so the sample is spread over
    the well the overview imaged, not just where the positives landed. The per-well
    seed is `seed + well_index`, so wells are reproducible yet decorrelated (SplitMix64
    makes consecutive seeds independent). Returns a dict with per-well results and
    generated FOV rows; writes the curated CSV if out_csv is given."""
    target_dir = find_target_dir(results_dir)
    name_index, order = discover_signals(target_dir)
    base = order[0]
    base_t = name_index[base]
    require = parse_gate(gate_spec)
    unknown = [s for s in require if s not in name_index]
    if unknown:
        raise ValueError("Gate signal(s) not found: %s (available: %s)" % (unknown, order))

    files_by_site = {}
    for fn in sorted(os.listdir(target_dir)):
        if fn.lower().endswith(".csv"):
            files_by_site.setdefault(_site_of(fn), {})[signal_name(fn)] = fn

    header, wells, fov_rows = None, [], []
    for well_index, site in enumerate(sorted(files_by_site)):
        files = files_by_site[site]
        if base not in files:
            continue

        head, base_rows = read_csv(os.path.join(target_dir, files[base]))
        header = header or head
        base_boxes = boxes(base_rows, base_t)

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
        captured = sum(1 for b in selected
                       if any(_fully_inside(b, t["cx"], t["cy"], half) for t in tiles))

        template = base_rows[0] if base_rows else {}
        for t in tiles:
            fov_rows.append(fov_row(template, base_t, t["cx"], t["cy"],
                                    "FOV_%d" % (len(fov_rows) + 1)))
        wells.append({"site": site, "base": base_boxes, "selected": selected,
                      "acquired": acquired, "captured": captured,
                      "extra": captured - len(acquired), "eligible": n_eligible,
                      "sample_short": n_eligible < sample_size, "tiles": tiles})

    if out_csv and header:
        _write_csv(out_csv, header, fov_rows)
    return {"target_dir": target_dir, "base": base, "signals": order,
            "wells": wells, "fov_rows": fov_rows, "header": header}


# --------------------------------------------------------------------------- #
# 11. Fiji front end - class dialog, visual report, and the run entry point    #
# ImageJ imports live inside these functions so the module still imports in     #
# plain CPython (for the tests).                                               #
# --------------------------------------------------------------------------- #
def _run_curation(results_path, gate_spec, out, fov_px, sample_size, seed, neighbourhood, stage_margin):
    """Curate, then write the curated CSV, per-well report PNGs and summary CSV. Split
    out of run_macro so it can be driven headlessly - the class dialog cannot be."""
    from ij import IJ, ImagePlus
    from ij.process import ColorProcessor
    from java.awt import Color

    curated = os.path.join(out, "TargetData_curated_FOVs.csv")
    res = curate(results_path, gate_spec, fov_px=fov_px, sample_size=sample_size, seed=seed,
                 neighbourhood=neighbourhood, stage_margin=stage_margin, out_csv=curated)

    report_dir = os.path.join(out, "report")
    if not os.path.isdir(report_dir):
        os.makedirs(report_dir)
    canvas = 1040
    grey, orange, red, blue = (Color(225, 225, 225), Color(244, 165, 130),
                               Color(202, 0, 32), Color(5, 113, 176))

    def render(well, path):
        coords = [c for b in well["base"] for c in (b[0] + b[2], b[1] + b[3])]
        scale = canvas / (max(coords) if coords else 1.0)
        ip = ColorProcessor(canvas, canvas)
        ip.setColor(Color.WHITE); ip.fill()

        def dots(bs, colour, r):
            ip.setColor(colour)
            for b in bs:
                ip.fillOval(int((b[0] + b[2] / 2.0) * scale) - r,
                            int((b[1] + b[3] / 2.0) * scale) - r, 2 * r, 2 * r)
        dots(well["base"], grey, 1)          # all nuclei
        dots(well["selected"], orange, 3)    # positive population
        dots(well["acquired"], red, 4)       # the sampled cells we image
        ip.setColor(blue); ip.setLineWidth(2)
        f = int(fov_px * scale)
        for t in well["tiles"]:
            ip.drawRect(int(t["cx"] * scale - f / 2.0), int(t["cy"] * scale - f / 2.0), f, f)
        IJ.saveAs(ImagePlus(well["site"], ip), "PNG", path)

    summary = [u"well,positives,eligible,acquired,extra,fovs,sample_short"]
    for well in res["wells"]:
        summary.append(u"%s,%d,%d,%d,%d,%d,%s" % (
            well["site"], len(well["selected"]), well["eligible"], len(well["acquired"]),
            well["extra"], len(well["tiles"]), well["sample_short"]))
        render(well, os.path.join(report_dir, "report_%s.png" % well["site"]))
    with io.open(os.path.join(out, "curation_summary.csv"), "w", encoding="utf-8", newline="") as fh:
        fh.write(u"\r\n".join(summary) + u"\r\n")

    acquired = sum(len(w["acquired"]) for w in res["wells"])
    extra = sum(w["extra"] for w in res["wells"])
    fovs = sum(len(w["tiles"]) for w in res["wells"])
    under = sum(1 for w in res["wells"] if w["sample_short"])
    IJ.log("  wells=%d acquired=%d extra-whole=%d FOVs=%d under-sampled=%d"
           % (len(res["wells"]), acquired, extra, fovs, under))
    IJ.log("  curated FOVs -> %s" % curated)
    return res


def run_macro():
    from ij import IJ
    from ij.gui import GenericDialog

    results_path = results_dir.getAbsolutePath()
    _, order = discover_signals(find_target_dir(results_path))
    base, markers = order[0], order[1:]           # base = first channel (the nuclei)

    gd = GenericDialog("MD HCS target curation")
    gd.addMessage("Base object (gated): " + base)
    gd.addMessage("Role of each marker class found in the data:")
    for sig in markers:
        gd.addChoice(sig, ["ignore", "positive", "negative"], "ignore")
    gd.showDialog()
    if gd.wasCanceled():
        return
    roles = dict((sig, gd.getNextChoice()) for sig in markers)
    gate_spec = build_gate_spec(markers, roles)

    mag_t, changer_t = _OBJECTIVES[target_objective]
    mag_ov, changer_ov, binning_ov, region_px = overview_scale(
        read_image_description(overview_image.getAbsolutePath()))
    fov_px = fov_px_for(mag_t, changer_t, mag_ov, changer_ov, binning_ov, region_px)

    IJ.log("MD HCS curation: base=%s gate='%s'" % (base, gate_spec))
    IJ.log("  overview %gx changer %gx binning %d -> target %s FOV=%.0f montage px"
           % (mag_ov, changer_ov, binning_ov, target_objective, fov_px))
    _run_curation(results_path, gate_spec, out_dir.getAbsolutePath(), fov_px,
                  int(sample_size), int(seed), neighbourhood, stage_margin)


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
    GOLDEN_TILES = [(400.5, 400.5), (200.5, 1000.5), (200.5, 1400.5),
                    (800.5, 400.5), (1200.5, 600.5), (1000.5, 1400.5)]
    GOLDEN_ACQUIRED, GOLDEN_FOVS, GOLDEN_CAPTURED, GOLDEN_EXTRA = 189, 188, 341, 152

    print("gating")
    n1, n2 = bx(0, 0, 10, 10), bx(100, 0, 10, 10)
    ch = [bx(5, 5, 4, 4)]                                   # touches n1 only
    check("edge-only contact is not overlap", touches(bx(0, 0, 10, 10), bx(10, 0, 5, 5)), False)
    check("yes gate", gate([n1, n2], {"m": ch}, {"m": True}), [n1])
    check("no gate", gate([n1, n2], {"m": ch}, {"m": False}), [n2])
    check("gate builder -> spec (channel order, ignore dropped)",
          build_gate_spec(["A", "B", "C"], {"A": "positive", "B": "ignore", "C": "negative"}), "A:yes; C:no")
    check("gate builder all-ignore -> empty", build_gate_spec(["A"], {"A": "ignore"}), "")
    check("gate builder round-trips via parse_gate",
          parse_gate(build_gate_spec(["mScarlet cells"], {"mScarlet cells": "positive"})), {"mScarlet cells": True})

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
