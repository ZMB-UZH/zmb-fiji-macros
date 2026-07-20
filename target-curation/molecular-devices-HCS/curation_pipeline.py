"""MD HCS target curation v2 - per-well pipeline over a Results folder.

Wraps curation_core: locate TargetData, discover signals, gate each well's base
objects, place FOVs, and emit generated FOV rows + per-well results. Python 2.7/3
compatible, stdlib only (runs in Fiji Jython). No ImageJ dependency here - the
visual report lives in the Fiji entry script.
"""
from __future__ import division, print_function
import io, os
import curation_core as cc

_SITE_KEY = "_singleTargetData_"


def find_target_dir(path):
    """Locate the `TargetData` folder: `path` itself if it is one, else a
    `TargetData` child, else one analysis-folder level down. A `TargetData`
    folder is preferred over sibling CSVs (the analysis folder also holds
    ObjectData/FieldData/metadata CSVs that must not be mistaken for it)."""
    if os.path.basename(os.path.normpath(path)) == "TargetData" and _has_csv(path):
        return path
    cand = os.path.join(path, "TargetData")
    if _has_csv(cand):
        return cand
    for name in sorted(os.listdir(path)):
        cand = os.path.join(path, name, "TargetData")
        if _has_csv(cand):
            return cand
    raise ValueError("No TargetData folder with CSVs found under: " + path)


def _has_csv(d):
    return os.path.isdir(d) and any(f.lower().endswith(".csv") for f in os.listdir(d))


def parse_gate(spec):
    """'mScarlet cells:yes; foo:no' -> {'mScarlet cells': True, 'foo': False}."""
    req = {}
    for part in spec.split(";"):
        sig, sep, want = part.rpartition(":")
        if sep and sig.strip():
            req[sig.strip()] = want.strip().lower() in ("yes", "y", "true", "1", "+")
    return req


def _site_of(fn):
    stem = fn[:-4]
    return stem.split(_SITE_KEY, 1)[1] if _SITE_KEY in stem else stem


def curate(results_dir, gate_spec, fov_px=256.0, sample_size=5, seed=42,
           neighbourhood=2.0, stage_margin=0.05, max_overlap=1.0,
           max_cells_per_fov=None, min_fovs=0, out_csv=None):
    """Curate every well: gate positives, draw a uniform random sample of
    sample_size of them, then place the fewest FOVs imaging that sample. base
    object = first signal (channel 1). Returns a dict with per-well results and
    the generated FOV rows; writes the curated CSV if out_csv is given.

    max_overlap defaults to 1.0 so every sampled observation is imaged; the FOV
    count is still minimised, so overlap only occurs where two sampled cells are
    geometrically forced to share ground. Lower it toward 0 to prefer strict
    disjoint FOVs at the cost of possibly not imaging some close sampled cells.
    See place_fovs for max_cells_per_fov and min_fovs."""
    target_dir = find_target_dir(results_dir)
    name_index, order = cc.discover_signals(target_dir)
    base = order[0]
    base_t = name_index[base]
    require = parse_gate(gate_spec)
    unknown = [s for s in require if s not in name_index]
    if unknown:
        raise ValueError("Gate signal(s) not found: %s (available: %s)" % (unknown, order))

    files_by_site = {}
    for fn in sorted(os.listdir(target_dir)):
        if fn.lower().endswith(".csv"):
            files_by_site.setdefault(_site_of(fn), {})[cc.signal_name(fn)] = fn

    header, wells, fov_rows = None, [], []
    for w, site in enumerate(sorted(files_by_site)):
        files = files_by_site[site]
        if base not in files:
            continue
        head, base_rows = cc.read_csv(os.path.join(target_dir, files[base]))
        header = header or head
        base_boxes = cc.boxes(base_rows, base_t)
        gate_boxes = {}
        for sig in require:
            rows = cc.read_csv(os.path.join(target_dir, files[sig]))[1] if sig in files else []
            gate_boxes[sig] = cc.boxes(rows, name_index.get(sig, base_t))
        selected = cc.gate(base_boxes, gate_boxes, require)
        # uniform sample of the positives (seed + well index -> reproducible,
        # decorrelated across wells), then the fewest FOVs covering that sample
        sampled, n_eligible = cc.sample_cells(selected, fov_px, stage_margin,
                                              neighbourhood, sample_size, seed + w)
        result = cc.place_fovs(sampled, fov_px, len(sampled), stage_margin, neighbourhood,
                               max_overlap, max_cells_per_fov, min_fovs)
        template = base_rows[0] if base_rows else {}
        for tile in result["tiles"]:
            fov_rows.append(cc.fov_row(template, base_t, tile["cx"], tile["cy"],
                                       "FOV_%d" % (len(fov_rows) + 1), 1.0))
        wells.append({"site": site, "base": base_boxes, "selected": selected,
                      "sampled": sampled, "eligible": n_eligible,
                      "sample_short": n_eligible < sample_size, "result": result})

    if out_csv and header:
        _write_csv(out_csv, header, fov_rows)
    return {"target_dir": target_dir, "base": base, "signals": order,
            "require": require, "wells": wells, "fov_rows": fov_rows, "header": header}


def _write_csv(path, header, rows):
    """Write CRLF CSV as unicode (portable across CPython and Jython, where
    csv.writer cannot write to an io text stream)."""
    with io.open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write(_csv_line(header))
        for r in rows:
            fh.write(_csv_line([r.get(c, u"") for c in header]))


def _csv_line(fields):
    return u",".join(_csv_field(f) for f in fields) + u"\r\n"


def _csv_field(value):
    s = u"{0}".format(value)
    if any(c in s for c in u',"\n\r'):
        s = u'"' + s.replace(u'"', u'""') + u'"'
    return s
