"""Resolve the HyCal calibration file of a run from database/runinfo/general.json
(same rules as the C++ RunInfoConfig: largest from_run <= run, chained
inheritance over the earlier entries; path relative to the database dir) and
fail loudly when a calibration did not load."""
import json
import os
import sys

DB_DIR = os.path.abspath(os.environ.get("PRAD2_DATABASE_DIR") or os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "database"))
CALIB_DIR = os.path.join(DB_DIR, "calibration")


def resolve_calib(run):
    """Absolute path of the calibration file runinfo assigns to `run`."""
    cfgs = json.load(open(os.path.join(DB_DIR, "runinfo", "general.json")))["configurations"]
    cfgs = sorted(((c.get("from_run", c.get("run_number")), c) for c in cfgs
                   if c.get("from_run", c.get("run_number")) is not None),
                  key=lambda t: t[0])
    f = ""
    for r0, c in cfgs:
        if r0 > run:
            break
        f = c.get("calibration", {}).get("file", f)
    if not f:
        sys.exit(f"[ERROR] runinfo has no calibration file for run {run}")
    path = f if os.path.isabs(f) else os.path.join(DB_DIR, f)
    if not os.path.exists(path):
        sys.exit(f"[ERROR] calibration file for run {run} not found: {path}")
    return path


def check_calib_json(path, min_n=100):
    """Number of modules with factor > 0 in a calibration JSON; exit if < min_n."""
    n = sum(float(m.get("factor", 0.0)) > 0 for m in json.load(open(path)))
    if n < min_n:
        sys.exit(f"[ERROR] calibration {path}: only {n} modules with factor > 0")
    return n


def check_pipeline(p, calib=""):
    """Number of PbWO4 modules with cal_factor > 0 after the pipeline loaded;
    exit if under half (a missing file leaves all of them at 0)."""
    hy = [p.hycal.module(i) for i in range(p.hycal.module_count())]
    pw = [m for m in hy if m.is_pwo4()]
    n = sum(m.cal_factor > 0 for m in pw)
    if n < 0.5 * len(pw):
        sys.exit(f"[ERROR] calibration did not load ({calib or 'runinfo'}): "
                 f"{n}/{len(pw)} PbWO4 modules calibrated")
    return n
