#@ File   (label="Analysis results folder", style="directory") results_dir
#@ String (label="Gate (e.g. 'mScarlet cells:yes; marker2:no')", value="mScarlet cells:yes") gate_spec
#@ Float  (label="Target FOV size (montage px)", value=256) fov_px
#@ Integer(label="Cells per well (uniform sample)", value=5, min=1) sample_size
#@ Integer(label="Random seed", value=42) seed
#@ Float  (label="Neighbourhood ratio (x object size)", value=2.0, min=0) neighbourhood
#@ Float  (label="Stage-accuracy margin (fraction per side)", value=0.05, min=0, max=0.49) stage_margin
#@ Float  (label="Max FOV overlap % (100 = image all sampled; lower = stricter disjoint)", value=100, min=0, max=100) max_overlap_pct
#@ Integer(label="Max cells per FOV (0 = off)", value=0, min=0) max_cells_per_fov
#@ Integer(label="Min FOVs per well (0 = off)", value=0, min=0) min_fovs
#@ File   (label="Output folder", style="directory") out_dir

"""Fiji entry for MD HCS target curation v2. Per well: gate positives, draw a
uniform random sample of them, place the fewest FOVs covering that sample; writes
a curated FOV list (to hand-load into MetaXpress) and a per-well visual report.
See CURATION_V2_SPEC.md.

The pure logic lives in curation_core / curation_pipeline (unit-tested in CPython
and this Jython); this script is the ImageJ I/O + dialog + report layer only.
"""
import io, os, sys

_here = os.path.dirname(os.path.realpath(__file__))
if _here not in sys.path:
    sys.path.insert(0, _here)

import curation_pipeline as cp
from ij import IJ, ImagePlus
from ij.process import ColorProcessor
from java.awt import Color

CANVAS = 1040                       # report canvas side, px
GREY = Color(225, 225, 225)         # all nuclei
ORANGE = Color(244, 165, 130)       # positive (gated)
RED = Color(202, 0, 32)             # sampled -> imaged
BLUE = Color(5, 113, 176)           # FOV box


def _extent(well):
    coords = [c for b in well["base"] for c in (b[0] + b[2], b[1] + b[3])]
    return max(coords) if coords else 1.0


def render(well, fov_px, out_png):
    scale = CANVAS / _extent(well)
    ip = ColorProcessor(CANVAS, CANVAS)
    ip.setColor(Color.WHITE)
    ip.fill()

    def dot(b, r):
        ip.fillOval(int((b[0] + b[2] / 2.0) * scale) - r,
                    int((b[1] + b[3] / 2.0) * scale) - r, 2 * r, 2 * r)

    ip.setColor(GREY)
    for b in well["base"]:
        dot(b, 1)
    ip.setColor(ORANGE)                 # positive population
    for b in well["selected"]:
        dot(b, 3)
    ip.setColor(RED)                    # the uniform sample we image
    for b in well["sampled"]:
        dot(b, 4)
    ip.setColor(BLUE)
    ip.setLineWidth(2)
    f = int(fov_px * scale)
    for t in well["result"]["tiles"]:
        ip.drawRect(int(t["cx"] * scale - f / 2.0), int(t["cy"] * scale - f / 2.0), f, f)

    IJ.saveAs(ImagePlus(well["site"], ip), "PNG", out_png)


def run():
    results = results_dir.getAbsolutePath()
    out = out_dir.getAbsolutePath()
    curated = os.path.join(out, "TargetData_curated_FOVs.csv")
    res = cp.curate(results, gate_spec, fov_px=fov_px, sample_size=int(sample_size),
                    seed=int(seed), neighbourhood=neighbourhood, stage_margin=stage_margin,
                    max_overlap=max_overlap_pct / 100.0,
                    max_cells_per_fov=(int(max_cells_per_fov) or None),
                    min_fovs=int(min_fovs), out_csv=curated)

    report_dir = os.path.join(out, "report")
    if not os.path.isdir(report_dir):
        os.makedirs(report_dir)
    summary = [u"well,positives,eligible,sampled,fovs,sample_short"]
    for w in res["wells"]:
        summary.append(u"%s,%d,%d,%d,%d,%s" % (w["site"], len(w["selected"]),
                       w["eligible"], len(w["sampled"]), len(w["result"]["tiles"]),
                       w["sample_short"]))
        render(w, fov_px, os.path.join(report_dir, "report_%s.png" % w["site"]))
    with io.open(os.path.join(out, "curation_summary.csv"), "w",
                 encoding="utf-8", newline="") as fh:
        fh.write(u"\r\n".join(summary) + u"\r\n")

    tiles = sum(len(w["result"]["tiles"]) for w in res["wells"])
    sampled = sum(len(w["sampled"]) for w in res["wells"])
    short = sum(1 for w in res["wells"] if w["sample_short"])
    IJ.log("MD HCS curation v2  base=%s  signals=%s" % (res["base"], res["signals"]))
    IJ.log("  wells=%d  sampled cells=%d  FOVs=%d  under-sampled wells=%d"
           % (len(res["wells"]), sampled, tiles, short))
    IJ.log("  curated FOVs -> %s" % curated)
    IJ.log("  report + summary -> %s" % out)
    print("wells=%d sampled=%d FOVs=%d underSampled=%d"
          % (len(res["wells"]), sampled, tiles, short))


run()
