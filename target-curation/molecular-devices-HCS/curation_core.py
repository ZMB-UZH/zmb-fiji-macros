"""MD HCS target curation v2 - core logic.

Python 2.7 / 3 compatible: unit-tested in CPython, runs in Fiji's Jython 2.7.4.
stdlib only, geometry in montage pixels. See CURATION_V2_SPEC.md.
"""
from __future__ import division, print_function
import csv, io, os, re

_BB = "AS_FID_Blob_BoundingBox"          # + X / Y / Width / Height
_OWN = re.compile(r"^T(\d+)\$")


def read_csv(path):
    with io.open(path, "r", encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.reader(fh))
    if not rows:
        return [], []
    head = rows[0]
    return head, [dict(zip(head, r)) for r in rows[1:] if len(r) >= len(head)]


def signal_name(fn):
    """A signal's name is the filename text before the first underscore."""
    return fn.split("_")[0]


def discover_signals(target_dir):
    """{name: T-index} + names ordered by (T-index, name) so the base (lowest T)
    is first. T-index = the T<n>$ prefix a signal's columns carry."""
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
    """The four bounding-box column names for target T<t>."""
    return ["T%d$%s%s" % (t, _BB, s) for s in ("X", "Y", "Width", "Height")]


def boxes(rows, t):
    """[(x, y, w, h)] from a target's own bbox columns; skip rows whose bbox is
    empty/non-numeric (placeholder rows)."""
    cols = _bb_cols(t)
    out = []
    for r in rows:
        try:
            out.append(tuple(float(r[c]) for c in cols))
        except (KeyError, ValueError):
            pass
    return out


def touches(a, b):
    """Overlap of two boxes with strictly positive area. Edge-only contact
    (shared boundary, zero area) does not count."""
    return not (a[0] + a[2] <= b[0] or b[0] + b[2] <= a[0] or
                a[1] + a[3] <= b[1] or b[1] + b[3] <= a[1])


def gate(base, gates, require):
    """Keep base boxes where each Yes signal touches and each No does not.
    gates = {signal: [boxes]}; require = {signal: True(Yes)/False(No)};
    signals absent from `require` are ignored (an empty require keeps all)."""
    return [bb for bb in base
            if all(any(touches(bb, g) for g in gates.get(s, [])) == want
                   for s, want in require.items())]


_MASK64 = (1 << 64) - 1


def _rng(seed):
    """SplitMix64 as a closure returning uniform floats in [0, 1). Integer core
    with 64-bit masking yields identical sequences in CPython 3 and Jython 2.7 -
    Python's own `random` module does NOT agree across the two interpreters. The
    top 53 bits give an exact double strictly < 1.0 (the full-64-bit form can
    round up to 1.0 and break the shuffle index)."""
    state = [seed & _MASK64]

    def nxt():
        state[0] = (state[0] + 0x9E3779B97F4A7C15) & _MASK64
        z = state[0]
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & _MASK64
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & _MASK64
        z ^= (z >> 31)
        return (z >> 11) / 9007199254740992.0    # 2**53, exact
    return nxt


def eligible(cells, fov, stage_margin, neighbourhood):
    """Indices of cells that can be validly observed alone: box + neighbourhood
    fits the FOV shrunk by the stage-accuracy margin."""
    usable = fov * (1.0 - 2.0 * stage_margin)
    return [i for i, b in enumerate(cells)
            if usable - b[2] * (1.0 + 2.0 * neighbourhood) >= 0
            and usable - b[3] * (1.0 + 2.0 * neighbourhood) >= 0]


def sample_cells(cells, fov, stage_margin, neighbourhood, n, seed):
    """Systematic Uniform Random Sampling (SURS) of n validly-observable cells -
    the stereology standard. Lay a grid sized for n points over the well's extent
    with ONE uniform-random offset (the random start), and take the nearest
    not-yet-chosen eligible cell at each grid point: even spatial coverage, but
    unbiased because the grid's position is randomised. A cell is eligible only if
    its box + neighbourhood fits a FOV, so a sampled cell is never clipped.
    Deterministic and interpreter-independent for a given seed. Returns
    (sampled_cells, n_eligible); exactly n cells when n_eligible > n."""
    elig = eligible(cells, fov, stage_margin, neighbourhood)
    n_elig = len(elig)
    if n_elig <= n:
        return [cells[i] for i in sorted(elig)], n_elig

    cx = [cells[i][0] + cells[i][2] / 2.0 for i in elig]
    cy = [cells[i][1] + cells[i][3] / 2.0 for i in elig]
    x0, x1, y0, y1 = min(cx), max(cx), min(cy), max(cy)
    w = (x1 - x0) or 1.0
    h = (y1 - y0) or 1.0
    ncol = max(1, int((n * w / h) ** 0.5 + 0.5))     # int(x+0.5): portable rounding
    nrow = max(1, (n + ncol - 1) // ncol)            # ceil(n / ncol) -> >= n points
    tx, ty = w / ncol, h / nrow
    nxt = _rng(seed)
    off_x, off_y = nxt(), nxt()                      # one random offset for the grid
    chosen, taken = [], set()
    for i in range(ncol):
        gx = x0 + (i + off_x) * tx
        for j in range(nrow):
            if len(chosen) >= n:
                break
            gy = y0 + (j + off_y) * ty
            best, best_d = -1, None
            for k in range(n_elig):
                if elig[k] in taken:
                    continue
                d = (cx[k] - gx) ** 2 + (cy[k] - gy) ** 2
                if best_d is None or d < best_d:
                    best_d, best = d, k
            taken.add(elig[best])
            chosen.append(elig[best])
    return [cells[i] for i in sorted(chosen)], n_elig


def place_fovs(cells, fov, min_cells, stage_margin, neighbourhood,
               max_overlap=0.0, max_cells_per_fov=None, min_fovs=0):
    """Fewest FOVs covering >= min_cells, each FOV centred so every cell it
    captures is fully inside (box + neighbourhood, within the FOV shrunk by the
    stage-accuracy margin).

    Returns {tiles: [{cx, cy, covered: [cell_index, ...]}], covered, short,
    fovs_short}. `short` is covered < min_cells; `fovs_short` is
    len(tiles) < min_fovs (the requested spread could not be reached).

    A FOV placed at (cx, cy) captures a cell iff cx/cy lie within that cell's
    per-axis half-window of its centre. A GROUP of cells therefore shares a valid
    centre iff the intersection of their per-axis feasibility intervals is
    non-empty; the tile is placed at that intersection's midpoint, so a cell is
    listed as covered only if it is provably inside the FOV actually placed. Every
    cell is assigned to at most one tile, so `covered` never double-counts.

    Greedy and deterministic: each round grows a candidate group outward from
    every uncovered cell and places the group covering the most cells (ties by
    centre then lowest index). Cells too large to fit any FOV are dropped up
    front; the greedy may report `short` in rare configs where an optimum would
    not.

    Optional knobs, all off by default (identity behaviour = fewest disjoint FOVs):
      max_overlap        fraction (0..1) of FOV area two tiles may share; 0 =
                         disjoint. Raise to let a sparse well reach min_cells when
                         disjoint tiles cannot.
      max_cells_per_fov  cap on cells per tile (None or 0 = no cap); forces spread.
      min_fovs           keep placing until at least this many tiles exist, even
                         once min_cells is met; forces spread across the well. If
                         unreachable, the result's `fovs_short` is True.
    """
    usable = fov * (1.0 - 2.0 * stage_margin)
    if usable <= 0:
        return {"tiles": [], "covered": 0, "short": min_cells > 0}

    # (centre_x, centre_y, half_window_x, half_window_y); drop cells that cannot fit
    info = {}
    for i, b in enumerate(cells):
        room_x = usable - b[2] * (1.0 + 2.0 * neighbourhood)
        room_y = usable - b[3] * (1.0 + 2.0 * neighbourhood)
        if room_x >= 0 and room_y >= 0:
            info[i] = (b[0] + b[2] / 2.0, b[1] + b[3] / 2.0, room_x / 2.0, room_y / 2.0)

    remaining = set(info)
    grid = _Grid(info, usable)
    placed, covered = [], 0
    while remaining and (covered < min_cells or len(placed) < min_fovs):
        best = None
        for seed in remaining:
            group, cx, cy = _grow(seed, remaining, info, grid, max_cells_per_fov)
            if _overlaps(cx, cy, placed, fov, max_overlap):
                continue
            key = (len(group), -cx, -cy, -min(group))
            if best is None or key > best[0]:
                best = (key, cx, cy, group)
        if best is None:
            break
        _, cx, cy, group = best
        placed.append({"cx": cx, "cy": cy, "covered": sorted(group)})
        covered += len(group)
        remaining.difference_update(group)
    return {"tiles": placed, "covered": covered, "short": covered < min_cells,
            "fovs_short": len(placed) < min_fovs}


def _overlaps(cx, cy, placed, fov, max_overlap):
    """True if a tile at (cx, cy) shares more than `max_overlap` of its area with
    any placed tile (equal FOV squares). max_overlap=0 requires full disjointness
    (any positive-area overlap is rejected)."""
    area = fov * fov
    for t in placed:
        ix = fov - abs(cx - t["cx"])
        iy = fov - abs(cy - t["cy"])
        if ix > 0 and iy > 0 and (ix * iy) / area > max_overlap + 1e-9:
            return True
    return False


class _Grid(object):
    """Uniform spatial hash on the usable-FOV cell size, for O(1) lookup of the
    cells near a point (the 3x3 block of buckets around it). Any two cells that
    could share a FOV lie within one usable-FOV of each other, so 3x3 suffices."""

    def __init__(self, info, size):
        self.size = size
        self.buckets = {}
        for i, (x, y, _, _) in info.items():
            self.buckets.setdefault((int(x // size), int(y // size)), []).append(i)

    def near(self, x, y):
        bx, by = int(x // self.size), int(y // self.size)
        return [i for dx in (-1, 0, 1) for dy in (-1, 0, 1)
                for i in self.buckets.get((bx + dx, by + dy), [])]


def _grow(seed, remaining, info, grid, cap):
    """Grow a group outward from `seed` (nearest cells first), keeping the
    per-axis feasibility intervals non-empty, and return (group, centre_x,
    centre_y) where the centre is the intersection midpoint (guaranteed to capture
    every group member). Stops at `cap` cells if given."""
    sx, sy, shx, shy = info[seed]
    xlo, xhi, ylo, yhi = sx - shx, sx + shx, sy - shy, sy + shy
    group = [seed]
    near = [i for i in grid.near(sx, sy) if i in remaining and i != seed]
    near.sort(key=lambda i: ((info[i][0] - sx) ** 2 + (info[i][1] - sy) ** 2, i))
    for i in near:
        if cap and len(group) >= cap:
            break
        x, y, hx, hy = info[i]
        nxlo, nxhi = max(xlo, x - hx), min(xhi, x + hx)
        nylo, nyhi = max(ylo, y - hy), min(yhi, y + hy)
        if nxlo <= nxhi and nylo <= nyhi:
            xlo, xhi, ylo, yhi = nxlo, nxhi, nylo, nyhi
            group.append(i)
    return group, (xlo + xhi) / 2.0, (ylo + yhi) / 2.0


def fov_row(template, t, cx, cy, object_id, size, centre=True):
    """A TargetData row templated on a real one, with target T<t>'s bbox placed
    at (cx, cy) (as box centre if `centre`, else top-left) and a fresh object_id.
    The FOV footprint MetaXpress images is set by the objective, not `size`; the
    box is only a position marker, so callers pass a small nominal `size`."""
    row = dict(template)
    row["object_id"] = str(object_id)
    x, y = (cx - size / 2.0, cy - size / 2.0) if centre else (cx, cy)
    for col, val in zip(_bb_cols(t), (x, y, size, size)):
        row[col] = repr(val)
    return row
