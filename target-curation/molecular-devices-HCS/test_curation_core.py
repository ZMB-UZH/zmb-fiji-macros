"""Tests for curation_core. Runs in CPython and Fiji Jython (stdlib only).

Synthetic fixtures with hand-computed expectations and structural invariants,
plus an auto-skipping test against the real NB26-15 dataset if present.
Run: python test_curation_core.py
"""
from __future__ import division, print_function
import io, os, shutil, tempfile
import curation_core as cc
import curation_pipeline as cp

_fails = []
EPS = 1e-9


def check(name, got, want):
    ok = got == want
    print(("  PASS " if ok else "  FAIL ") + name)
    if not ok:
        print("        got  {0!r}\n        want {1!r}".format(got, want))
        _fails.append(name)


def box(x, y, w, h):
    return (float(x), float(y), float(w), float(h))


def assert_valid(name, cells, res, fov, stage_margin, neighbourhood, disjoint=True):
    """Structural invariants that must ALWAYS hold: every covered cell lies inside
    the FOV actually placed; no cell is counted twice; `covered` equals the number
    of distinct covered cells. Tiles are also pairwise disjoint unless overlap is
    allowed (disjoint=False)."""
    usable = fov * (1.0 - 2.0 * stage_margin)
    inside = True
    for t in res["tiles"]:
        for i in t["covered"]:
            b = cells[i]
            cx, cy = b[0] + b[2] / 2.0, b[1] + b[3] / 2.0
            hx = (usable - b[2] * (1.0 + 2.0 * neighbourhood)) / 2.0
            hy = (usable - b[3] * (1.0 + 2.0 * neighbourhood)) / 2.0
            if abs(cx - t["cx"]) > hx + EPS or abs(cy - t["cy"]) > hy + EPS:
                inside = False
    flat = [i for t in res["tiles"] for i in t["covered"]]
    check(name + ": covered cells inside placed FOV", inside, True)
    check(name + ": no double-count", len(flat) == len(set(flat)) == res["covered"], True)
    if disjoint:
        ok = all(
            abs(res["tiles"][a]["cx"] - res["tiles"][b]["cx"]) >= fov - EPS or
            abs(res["tiles"][a]["cy"] - res["tiles"][b]["cy"]) >= fov - EPS
            for a in range(len(res["tiles"])) for b in range(a + 1, len(res["tiles"])))
        check(name + ": tiles disjoint", ok, True)


# --- touches ---------------------------------------------------------------- #
def test_touches():
    print("touches")
    check("overlap", cc.touches(box(0, 0, 10, 10), box(5, 5, 10, 10)), True)
    check("edge-only contact does not count", cc.touches(box(0, 0, 10, 10), box(10, 0, 5, 5)), False)
    check("clear gap", cc.touches(box(0, 0, 10, 10), box(20, 0, 5, 5)), False)


# --- discover_signals ------------------------------------------------------- #
def test_discover():
    print("discover_signals")
    d = tempfile.mkdtemp()
    try:
        def w(fn, header):
            with io.open(os.path.join(d, fn), "w", encoding="utf-8") as fh:
                fh.write(header + u"\n1\n")
        w("DAPI_singleTargetData_R2-C3.csv", u"object_id,T1$AS_FID_Blob_Size")
        w("DAPI_singleTargetData_R2-C4.csv", u"object_id,T1$AS_FID_Blob_Size")
        w("mScarlet cells_singleTargetData_R2-C3.csv", u"object_id,T2$AS_FID_Blob_Size")
        w("notes.txt", u"ignored")                       # non-CSV ignored
        w("summary_singleTargetData_R2-C3.csv", u"object_id,plate_id")  # no T<n>$ -> skipped
        idx, order = cc.discover_signals(d)
        check("index map", idx, {"DAPI": 1, "mScarlet cells": 2})
        check("ordered base-first", order, ["DAPI", "mScarlet cells"])
    finally:
        shutil.rmtree(d)


# --- read_csv --------------------------------------------------------------- #
def test_read_csv():
    print("read_csv")
    d = tempfile.mkdtemp()
    try:
        p = os.path.join(d, "x.csv")
        with io.open(p, "w", encoding="utf-8-sig") as fh:      # writes a BOM
            fh.write(u"a,b,c\r\n1,2,3\r\nragged,row\r\n")
        head, rows = cc.read_csv(p)
        check("BOM stripped from first header", head[0], "a")
        check("ragged row dropped", len(rows), 1)
        check("row parsed to dict", rows[0], {"a": "1", "b": "2", "c": "3"})
    finally:
        shutil.rmtree(d)


# --- boxes ------------------------------------------------------------------ #
def test_boxes():
    print("boxes")
    rows = [
        {"T1$AS_FID_Blob_BoundingBoxX": "10", "T1$AS_FID_Blob_BoundingBoxY": "20",
         "T1$AS_FID_Blob_BoundingBoxWidth": "5", "T1$AS_FID_Blob_BoundingBoxHeight": "6"},
        {"T1$AS_FID_Blob_BoundingBoxX": "", "T1$AS_FID_Blob_BoundingBoxY": "",
         "T1$AS_FID_Blob_BoundingBoxWidth": "", "T1$AS_FID_Blob_BoundingBoxHeight": ""},
        {"other": "col"},                                # missing bbox cols
    ]
    bs = cc.boxes(rows, 1)
    check("skips empty and missing rows", len(bs), 1)
    check("parsed bbox", bs[0][:4], (10.0, 20.0, 5.0, 6.0))


# --- gate ------------------------------------------------------------------- #
def test_gate():
    print("gate")
    n1, n2, n3 = box(0, 0, 10, 10), box(100, 0, 10, 10), box(200, 0, 10, 10)
    ch2 = [box(5, 5, 4, 4)]                              # touches n1 only
    check("Yes selects touching", cc.gate([n1, n2, n3], {"m": ch2}, {"m": True}), [n1])
    check("No selects non-touching", cc.gate([n1, n2, n3], {"m": ch2}, {"m": False}), [n2, n3])
    check("empty require selects all", cc.gate([n1, n2, n3], {"m": ch2}, {}), [n1, n2, n3])
    check("Yes with signal absent from gates -> none",
          cc.gate([n1, n2, n3], {}, {"m": True}), [])
    ch3 = [box(1, 1, 3, 3)]                              # also touches n1
    check("double positive (AND)",
          cc.gate([n1, n2, n3], {"m": ch2, "k": ch3}, {"m": True, "k": True}), [n1])


# --- place_fovs ------------------------------------------------------------- #
def test_place():
    print("place_fovs")
    # close pair -> one FOV covers both
    cells = [box(0, 0, 10, 10), box(20, 0, 10, 10)]
    r = cc.place_fovs(cells, 100, 2, stage_margin=0.0, neighbourhood=0.0)
    check("close pair -> 1 tile", len(r["tiles"]), 1)
    check("both covered, not short", (r["covered"], r["short"]), (2, False))
    assert_valid("close pair", cells, r, 100, 0.0, 0.0)

    # far pair -> two disjoint FOVs
    cells = [box(0, 0, 10, 10), box(1000, 0, 10, 10)]
    r = cc.place_fovs(cells, 100, 2, stage_margin=0.0, neighbourhood=0.0)
    check("far pair -> 2 tiles", len(r["tiles"]), 2)
    assert_valid("far pair", cells, r, 100, 0.0, 0.0)

    # ASYMMETRIC cluster (the over-claim bug): points at x=0,50,73,74,75 span 75 <=
    # usable 100, so one FOV centred in [25,50] captures all five. The old
    # mean-centroid placed at 54.4 and clipped x=0; feasible-interval must not.
    cells = [box(x, 0, 0, 0) for x in (0, 50, 73, 74, 75)]
    r = cc.place_fovs(cells, 100, 5, stage_margin=0.0, neighbourhood=0.0)
    check("asymmetric span fits one FOV", (len(r["tiles"]), r["covered"], r["short"]), (1, 5, False))
    assert_valid("asymmetric", cells, r, 100, 0.0, 0.0)

    # three separated clusters -> three tiles, invariants hold
    cells = ([box(x, 0, 0, 0) for x in (0, 10, 20)] +
             [box(x, 500, 0, 0) for x in (500, 510)] +
             [box(1000, 1000, 0, 0)])
    r = cc.place_fovs(cells, 100, 6, stage_margin=0.0, neighbourhood=0.0)
    check("three clusters -> 3 tiles, all covered", (len(r["tiles"]), r["covered"]), (3, 6))
    assert_valid("three clusters", cells, r, 100, 0.0, 0.0)

    # neighbourhood too large -> cell cannot fit, well is short
    r = cc.place_fovs([box(0, 0, 30, 30)], 100, 1, stage_margin=0.0, neighbourhood=2.0)
    check("oversized cell excluded", (len(r["tiles"]), r["short"]), (0, True))

    # per-axis window: a wide-short cell fits in y but its width may exclude it
    r = cc.place_fovs([box(0, 0, 90, 10)], 100, 1, stage_margin=0.0, neighbourhood=0.0)
    check("wide cell within width fits", (len(r["tiles"]), r["short"]), (1, False))
    r = cc.place_fovs([box(0, 0, 110, 10)], 100, 1, stage_margin=0.0, neighbourhood=0.0)
    check("too-wide cell excluded on x axis", (len(r["tiles"]), r["short"]), (0, True))

    # stage margin shrinks the usable window enough to exclude a cell
    r = cc.place_fovs([box(0, 0, 95, 10)], 100, 1, stage_margin=0.0, neighbourhood=0.0)
    check("fits at margin 0", r["short"], False)
    r = cc.place_fovs([box(0, 0, 95, 10)], 100, 1, stage_margin=0.05, neighbourhood=0.0)
    check("excluded once stage margin applied", r["short"], True)

    # degenerate: no usable area -> nothing fits, honest short
    r = cc.place_fovs([box(0, 0, 0, 0)], 100, 1, stage_margin=0.5, neighbourhood=0.0)
    check("no usable area -> short, no crash", (len(r["tiles"]), r["short"]), (0, True))

    # min unreachable -> short but returns what it could
    r = cc.place_fovs([box(0, 0, 0, 0)], 100, 5, stage_margin=0.0, neighbourhood=0.0)
    check("short when min unreachable", (r["covered"], r["short"]), (1, True))

    # determinism: identical inputs -> identical tile centres
    cells = [box(x, y, 0, 0) for x in (0, 40, 600) for y in (0, 30)]
    a = cc.place_fovs(cells, 100, 10, 0.05, 1.0)
    b = cc.place_fovs(cells, 100, 10, 0.05, 1.0)
    check("deterministic tile centres",
          [(t["cx"], t["cy"]) for t in a["tiles"]], [(t["cx"], t["cy"]) for t in b["tiles"]])
    check("deterministic covered sets",
          [t["covered"] for t in a["tiles"]], [t["covered"] for t in b["tiles"]])


# --- optional knobs --------------------------------------------------------- #
def test_knobs():
    print("optional knobs")
    # overlap tolerance: a stage margin makes two cells unshareable, yet their
    # tiles sit < FOV apart, so disjoint placement can only image one.
    pair = [box(0, 0, 0, 0), box(80, 0, 0, 0)]
    off = cc.place_fovs(pair, 100, 2, 0.2, 0.0)                     # overlap off -> disjoint
    check("overlap off -> pair short", (off["covered"], off["short"]), (1, True))
    tol = cc.place_fovs(pair, 100, 2, 0.2, 0.0, max_overlap=0.25)
    check("overlap tolerance images both",
          (tol["covered"], len(tol["tiles"]), tol["short"]), (2, 2, False))
    assert_valid("overlap tol", pair, tol, 100, 0.2, 0.0, disjoint=False)

    # max cells per FOV (overlap allowed so a tight cluster can be split up)
    cluster = [box(float(x), 0, 0, 0) for x in (0, 5, 10, 15, 20)]
    capped = cc.place_fovs(cluster, 100, 5, 0.0, 0.0, max_overlap=1.0, max_cells_per_fov=2)
    check("no tile exceeds cap", max(len(t["covered"]) for t in capped["tiles"]), 2)
    check("cap: all covered via more tiles", (capped["covered"], len(capped["tiles"])), (5, 3))
    assert_valid("cap", cluster, capped, 100, 0.0, 0.0, disjoint=False)
    check("cap=1 -> one cell per tile",
          max(len(t["covered"]) for t in
              cc.place_fovs(cluster, 100, 5, 0.0, 0.0, max_overlap=1.0, max_cells_per_fov=1)["tiles"]), 1)
    check("cap=0 means no cap (0=off convention)",
          len(cc.place_fovs(cluster, 100, 5, 0.0, 0.0, max_cells_per_fov=0)["tiles"]), 1)

    # min FOVs forces spread beyond what min_cells alone needs
    spread = [box(0, 0, 0, 0), box(500, 0, 0, 0), box(1000, 0, 0, 0)]
    check("min_fovs off -> stops at min_cells",
          len(cc.place_fovs(spread, 100, 1, 0.0, 0.0)["tiles"]), 1)
    forced = cc.place_fovs(spread, 100, 1, 0.0, 0.0, min_fovs=3)
    check("min_fovs forces >= 3 tiles",
          (len(forced["tiles"]), forced["covered"], forced["fovs_short"]), (3, 3, False))
    unmet = cc.place_fovs(spread, 100, 1, 0.0, 0.0, min_fovs=5)
    check("min_fovs unreachable -> fovs_short flagged",
          (len(unmet["tiles"]), unmet["fovs_short"]), (3, True))

    # every knob is disableable: explicit off-values reproduce the default
    a = cc.place_fovs(cluster, 100, 5, 0.0, 0.0, max_overlap=0.0, max_cells_per_fov=None, min_fovs=0)
    b = cc.place_fovs(cluster, 100, 5, 0.0, 0.0)
    check("off knobs == default",
          [t["covered"] for t in a["tiles"]], [t["covered"] for t in b["tiles"]])


# --- fov_row ---------------------------------------------------------------- #
def test_fov_row():
    print("fov_row")
    tmpl = {"object_id": "7", "well_label": "B - 3", "T1$AS_FID_Blob_BoundingBoxX": "0",
            "T1$AS_FID_Blob_BoundingBoxY": "0", "T1$AS_FID_Blob_BoundingBoxWidth": "0",
            "T1$AS_FID_Blob_BoundingBoxHeight": "0"}
    row = cc.fov_row(tmpl, 1, 100.0, 200.0, "FOV_1", 2.0, centre=True)
    check("identity kept", row["well_label"], "B - 3")
    check("fresh object_id", row["object_id"], "FOV_1")
    check("bbox centred at point",
          (float(row["T1$AS_FID_Blob_BoundingBoxX"]), float(row["T1$AS_FID_Blob_BoundingBoxY"])),
          (99.0, 199.0))
    corner = cc.fov_row(tmpl, 1, 100.0, 200.0, "FOV_2", 2.0, centre=False)
    check("corner mode places top-left",
          (float(corner["T1$AS_FID_Blob_BoundingBoxX"]), float(corner["T1$AS_FID_Blob_BoundingBoxY"])),
          (100.0, 200.0))


# --- sampling (systematic uniform random sampling) -------------------------- #
def test_sampling():
    print("sampling (SURS)")
    grid = [box(float(x * 10), float(y * 10), 2.0, 2.0) for x in range(5) for y in range(5)]
    check("eligible excludes oversized",
          len(cc.eligible(grid + [box(0, 0, 500, 500)], 100, 0.0, 0.0)), 25)
    s, n = cc.sample_cells(grid, 100, 0.0, 0.0, 5, 42)
    check("exactly n sampled + eligible count", (len(s), n), (5, 25))
    allc, n2 = cc.sample_cells(grid, 100, 0.0, 0.0, 50, 42)
    check("take all when population <= n", (len(allc), n2), (25, 25))
    xy = [(b[0], b[1]) for b in cc.sample_cells(grid, 100, 0.0, 0.0, 5, 7)[0]]
    print("    sample(seed 7) = %s" % xy)
    check("reproducible for a seed", xy, [(b[0], b[1]) for b in cc.sample_cells(grid, 100, 0.0, 0.0, 5, 7)[0]])
    check("portable sample golden (seed 7)", xy, GOLDEN_SAMPLE)
    check("distinct (no replacement)", len(xy), len(set(xy)))
    check("different seed differs",
          xy != [(b[0], b[1]) for b in cc.sample_cells(grid, 100, 0.0, 0.0, 5, 8)[0]], True)

    # SURS spatial evenness: a 2x2 grid (n=4) lands one cell in each quadrant
    # essentially always; pure random (SRS) would only ~9% of the time.
    field = [box(float(x * 10), float(y * 10), 1.0, 1.0) for x in range(20) for y in range(20)]
    quad = lambda b: (b[0] >= 95.0, b[1] >= 95.0)
    allfour = sum(1 for sd in range(200)
                  if len(set(quad(b) for b in cc.sample_cells(field, 100, 0.0, 0.0, 4, sd)[0])) == 4)
    print("    SURS all-4-quadrants = %.2f  (SRS would be ~0.09)" % (allfour / 200.0))
    check("SURS covers the space evenly", allfour / 200.0 > 0.9, True)

    # edges
    check("n=0 -> empty, eligible counted", cc.sample_cells(grid, 100, 0.0, 0.0, 0, 1), ([], 25))
    check("empty cells -> empty", cc.sample_cells([], 100, 0.0, 0.0, 5, 1), ([], 0))
    check("all ineligible -> empty",
          cc.sample_cells([box(0, 0, 500, 500)], 100, 0.0, 0.0, 5, 1), ([], 0))
    check("n_eligible == n takes all",
          cc.sample_cells([box(float(x), 0, 2, 2) for x in range(5)], 100, 0.0, 0.0, 5, 3)[1], 5)
    snapshot = list(grid)
    cc.sample_cells(grid, 100, 0.0, 0.0, 5, 1)
    check("caller cells not mutated", grid, snapshot)


GOLDEN_SAMPLE = [(10.0, 0.0), (10.0, 10.0), (10.0, 30.0), (30.0, 0.0), (30.0, 10.0)]  # SURS, portable


# --- pipeline --------------------------------------------------------------- #
def _write_target(path, t, boxes):
    cols = ["T%d$AS_FID_Blob_BoundingBox%s" % (t, s) for s in ("X", "Y", "Width", "Height")]
    with io.open(path, "w", encoding="utf-8") as fh:
        fh.write(u",".join([u"object_id"] + cols) + u"\n")
        for i, b in enumerate(boxes):
            fh.write(u",".join([str(i)] + [str(v) for v in b]) + u"\n")


def test_pipeline():
    print("pipeline")
    check("parse yes", cp.parse_gate("m:yes"), {"m": True})
    check("parse yes/no", cp.parse_gate("a:yes; b:no"), {"a": True, "b": False})
    check("parse name with space", cp.parse_gate("mScarlet cells:yes"), {"mScarlet cells": True})
    check("parse no colon -> empty", cp.parse_gate("m"), {})
    check("site parsed", cp._site_of("DAPI_singleTargetData_R2-C3-F0-Z0-T0.csv"), "R2-C3-F0-Z0-T0")
    check("site without key", cp._site_of("x.csv"), "x")
    check("csv plain", cp._csv_field("ab"), "ab")
    check("csv comma quoted", cp._csv_field("a,b"), '"a,b"')
    check("csv quote doubled", cp._csv_field('a"b'), '"a""b"')

    d = tempfile.mkdtemp()
    try:
        td = os.path.join(d, "TargetData")
        os.mkdir(td)
        _write_target(os.path.join(td, "DAPI_singleTargetData_S1.csv"), 1,
                      [(0, 0, 10, 10), (100, 100, 10, 10)])
        _write_target(os.path.join(td, "mScarlet cells_singleTargetData_S1.csv"), 2,
                      [(5, 5, 4, 4)])                        # overlaps first DAPI only
        with io.open(os.path.join(d, "DAPI_ObjectData.csv"), "w", encoding="utf-8") as fh:
            fh.write(u"object_id\n1\n")                      # sibling CSV must be ignored
        check("find_target_dir prefers TargetData", cp.find_target_dir(d), td)

        raised = False
        try:
            cp.curate(d, "nope:yes", fov_px=100, sample_size=1, neighbourhood=0.0, stage_margin=0.0)
        except ValueError:
            raised = True
        check("unknown gate signal raises", raised, True)

        res = cp.curate(d, "mScarlet cells:yes", fov_px=100, sample_size=1,
                        neighbourhood=0.0, stage_margin=0.0)
        check("gate+sample selected the one positive", len(res["fov_rows"]), 1)
        check("base is DAPI", res["base"], "DAPI")
    finally:
        shutil.rmtree(d)


# --- real data (auto-skip) -------------------------------------------------- #
REAL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Results",
                    "mScarlet Cells_2026-Jul-13-17-23-33-077", "TargetData")


def test_real():
    print("real data")
    if not os.path.isdir(REAL):
        print("  SKIP (dataset not present)")
        return
    idx, order = cc.discover_signals(REAL)
    check("signals discovered", order, ["DAPI", "mScarlet cells"])
    pos = 0
    dense = None
    for fn in sorted(os.listdir(REAL)):
        if not fn.startswith("DAPI"):
            continue
        site = fn.split("_singleTargetData_")[1]
        base = cc.boxes(cc.read_csv(os.path.join(REAL, fn))[1], 1)
        m = cc.boxes(cc.read_csv(os.path.join(REAL, "mScarlet cells_singleTargetData_" + site))[1], 2)
        sel = cc.gate(base, {"m": m}, {"m": True})
        pos += len(sel)
        if site.startswith("R4-C4"):
            dense = (sel, cc.place_fovs(sel, 256.0, 5, stage_margin=0.05, neighbourhood=2.0))
    print("    plate positives (strict overlap): {0}".format(pos))
    print("    R4-C4: {0} positives, {1} tiles, covered {2}".format(
        len(dense[0]), len(dense[1]["tiles"]), dense[1]["covered"]))
    check("plate positives (strict-overlap golden)", pos, 4253)
    check("dense well satisfied, not short",
          bool(dense[1]["tiles"]) and not dense[1]["short"], True)
    assert_valid("R4-C4 real", dense[0], dense[1], 256.0, 0.05, 2.0)

    # full pipeline over the plate (gate -> uniform sample -> place); defaults
    res = cp.curate(os.path.dirname(REAL), "mScarlet cells:yes")
    tiles = sum(len(w["result"]["tiles"]) for w in res["wells"])
    sampled = sum(len(w["sampled"]) for w in res["wells"])
    under = sum(1 for w in res["wells"] if w["sample_short"])
    place_short = sum(1 for w in res["wells"] if w["result"]["short"])
    print("    plate: wells=%d sampled=%d FOVs=%d underSampled=%d placeShort=%d"
          % (len(res["wells"]), sampled, tiles, under, place_short))
    check("plate wells", len(res["wells"]), 48)
    check("plate under-sampled wells (golden)", under, 16)
    check("plate sampled cells (golden)", sampled, GOLDEN_SAMPLED)
    check("plate total FOVs (golden)", tiles, GOLDEN_TILES)
    check("every sampled cell imaged (no placement short)", place_short, 0)
    # reproducible across runs (same seed)
    res2 = cp.curate(os.path.dirname(REAL), "mScarlet cells:yes")
    check("plate reproducible FOV total",
          sum(len(w["result"]["tiles"]) for w in res2["wells"]), tiles)


GOLDEN_SAMPLED = 167   # SURS, seed 42, portable across CPython & Jython
GOLDEN_TILES = 165


def main():
    for t in (test_touches, test_discover, test_read_csv, test_boxes, test_gate,
              test_place, test_knobs, test_sampling, test_fov_row, test_pipeline, test_real):
        t()
    print("\n" + ("ALL PASS" if not _fails else "{0} FAILED: {1}".format(len(_fails), _fails)))
    return 1 if _fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
