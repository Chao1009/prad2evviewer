#!/usr/bin/env python3
"""p1b_report.py report --out <dir> <scan dir> [--set ...] <gates>: AMP P1b report driver (diagnostic only; the benchmark is not edited).

The benchmark's report with extra ENERGY_DEF values (the columns of a p1b_scan.py scan) and two environment switches, installed on the module
before the workers fork:
  ENERGY_DEF  'sum3' 'sum4' 'sum5' 'sum7' (h_s3, h_s4, h_s5, h_s7: fixed sums around the peak sample), 'cfd5' 'cfd7' (h_c5, h_c7: fixed sums from the
              CFD sample), plus the benchmark's own 'rawpk' 'parpk'; the scale is the benchmark's global MeV/ADC factor of the column
  P1B_MIX     "above:E" the column is used for hits with integral energy h_e >= E MeV and the integral below; "below:E" the reverse;
              "k1" the column only for single-pulse hits (h_k == 1, not extra) and the integral for the others; "kn" the reverse
  P1B_THR     member threshold E MeV: hits below it are removed before the clusterer, for every arm (also integral)
  P1B_CORR    "<tables.npz>:<kinds>" (p1b_corr_table.py output, kinds "phase", "amp" or "phase+amp"): the column is divided by the module's
              sampling-phase curve and/or the amplitude curve g(A) before the scale, i.e. the hit-level calibration against h_a (diagnostic: how much of
              the excess a table can remove)
"""
import os
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
import energy_resolution_benchmark as bm        # noqa: E402

bm.AMP_COLUMN.update({"sum3": "h_s3", "sum4": "h_s4", "sum5": "h_s5", "sum7": "h_s7", "cfd5": "h_c5", "cfd7": "h_c7"})
MIX, THR = os.environ.get("P1B_MIX", ""), float(os.environ.get("P1B_THR", "0"))
CORR = os.environ.get("P1B_CORR", "")
_amp, _cl = bm.amp_hits, bm.cluster_events
_T = {}


def corr_factor(h, fn, geo, col):
    """Hit-level factor curve[module, phase] * g(A) of the P1B_CORR tables (1 where the column is not finite)."""
    path, _, kinds = CORR.partition(":")
    if "t" not in _T:
        _T["t"] = np.load(path)
        z = np.load(fn)
        _T["t50"] = bm.template_columns(bm.PYDIR / str(z["bank"]), z["mod_name"])[0]
    t = _T["t"]
    m = h["h_mod"].astype(int)
    f = np.ones(m.size)
    if "phase" in kinds:
        x = ((h["h_t"].astype(np.float64) - _T["t50"][m]) % bm.DT) / bm.DT
        f *= t["curve"][m, np.minimum((x * 16).astype(int), 15)]
    if "amp" in kinds:
        ga = np.log(t["g_a"])
        f *= np.interp(np.log(np.maximum(col, 1e-3) / 1.0), ga, t["g"])
    return np.where(np.isfinite(col) & np.isfinite(f) & (f > 0), f, 1.0)


def amp_hits(h, fn, geo):
    if CORR:
        col = h["h_amp"].astype(np.float64)
        h = h | {"h_amp": (col / corr_factor(h, fn, geo, col)).astype(h["h_amp"].dtype)}
    out, s, fb = _amp(h, fn, geo)
    if MIX:
        kind, _, v = MIX.partition(":")
        e0 = h["h_e"]
        single = (h["h_k"] == 1) & ~h["h_extra"]
        integral = {"above": lambda: e0 < float(v), "below": lambda: e0 >= float(v), "k1": lambda: ~single, "kn": lambda: single}[kind]()
        out = out | {"h_e": np.where(integral, e0, out["h_e"]).astype(e0.dtype)}
    return out, s, fb


def cluster_events(h, *a, **k):
    if THR > 0:
        h = h | {"h_e": np.where(h["h_e"] < THR, 0.0, h["h_e"]).astype(h["h_e"].dtype)}
    return _cl(h, *a, **k)


if __name__ == "__main__":
    bm.amp_hits, bm.cluster_events = amp_hits, cluster_events
    sys.argv = [sys.argv[0]] + sys.argv[1:]
    bm.main()
