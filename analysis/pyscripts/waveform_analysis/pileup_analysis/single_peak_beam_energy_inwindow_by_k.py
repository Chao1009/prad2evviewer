#!/usr/bin/env python3
"""single_peak_beam_energy_inwindow_by_k.py -- in-window channels by fitted pulse number k, same scans and 3 sigma windows as single_peak_beam_energy_examples.py.

    source setup_env.sh; python3 waveform_analysis/pileup_analysis/single_peak_beam_energy_inwindow_by_k.py [--emin-frac F] [--stem <name>]
In window: main-pulse E = (largest A) x module gain within |E/E_beam - 1| <= window.  Per run type and k class (1, 2, 3, >=4): count, median chi2/dof,
fraction chi2/dof > 20; k = 2 also median r = A_second / A_main and median |dt| (ns).  -> <stem>.tsv (default single_peak_beam_energy_inwindow_by_k_20261116)
--emin-frac F replaces the window by the one-sided cut E > F x E_beam (no upper bound).
"""
import argparse
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE))
import numpy as np                                      # noqa: E402
import single_peak_beam_energy_examples as sp           # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--emin-frac", type=float, default=None, help="replace the window by the one-sided cut E > F x E_beam (default: 3 sigma window)")
ap.add_argument("--stem", default="single_peak_beam_energy_inwindow_by_k_20261116", help="output file stem")
args = ap.parse_args()
pr, pmo = sp.pr, sp.pmo
sample, full = pmo.PYDIR / "output/pileup_rates/monoedge_sample/old", pmo.PYDIR / "output/pileup_rates/scan"
rows = [("run_type", "window", "k", "n", "median_chi2dof", "frac_chi2gt20", "frac_of_window", "median_r_k2", "median_dt_ns_k2", "frac_k2_r_lt_0.1_dt_lt_40")]
for rt in sp.TYPES:
    half = sp.WIN3[rt]
    wlab = f"{half}" if args.emin_frac is None else f"E>{args.emin_frac:g}Ebeam"
    acc = dict(k=[], chi=[], r=[], dt=[])
    for fn in sorted((full if rt in sp.FULL_FIRST else sample).glob("run_*.npz")):
        z = np.load(fn)
        if str(z["run_type"]) != rt:
            continue
        cfg = pr.run_config(int(z["run"]))
        gain = pr.module_gain(cfg, rt, [str(n) for n in z["mod_name"]])
        k, a, t = z["ch_k"], z["ch_a"], z["ch_t"]
        valid = np.arange(a.shape[1])[None, :] < k[:, None]
        am = np.where(valid, a, -1).max(1)
        E = am * gain[z["ch_midx"]] / cfg["beam"]                                   # main-pulse energy / beam energy
        inside = np.abs(E - 1) <= half if args.emin_frac is None else E > args.emin_frac
        sel = np.nonzero((k >= 1) & inside)[0]
        a, t, k = a[sel], t[sel], k[sel]
        o = np.argsort(-np.where(np.arange(6)[None, :] < k[:, None], a, -1), axis=1)
        am_, tm_ = np.take_along_axis(a, o[:, :1], 1)[:, 0], np.take_along_axis(t, o[:, :1], 1)[:, 0]
        a2, t2 = np.take_along_axis(a, o[:, 1:2], 1)[:, 0], np.take_along_axis(t, o[:, 1:2], 1)[:, 0]
        acc["k"].append(k)
        acc["chi"].append(z["ch_chi2"][sel])
        acc["r"].append(np.where(k == 2, a2 / am_, np.nan))
        acc["dt"].append(np.where(k == 2, np.abs(t2 - tm_), np.nan))
    k, chi, r, dt = (np.concatenate(acc[x]) for x in ("k", "chi", "r", "dt"))
    tot = k.size
    for name, m in (("1", k == 1), ("2", k == 2), ("3", k == 3), (">=4", k >= 4), ("all", k >= 1)):
        n = int(m.sum())
        f = lambda v: f"{v:.4g}"
        small = f(np.mean((r[m] < 0.1) & (dt[m] < 40))) if name == "2" and n else ""
        rows.append((rt, wlab, name, str(n), f(np.median(chi[m])) if n else "", f(np.mean(chi[m] > 20)) if n else "", f(n / tot),
                     f(np.median(r[m])) if name == "2" and n else "", f(np.median(dt[m])) if name == "2" and n else "", small))
out = HERE / f"{args.stem}.tsv"
out.write_text("\n".join("\t".join(x) for x in rows) + "\n")
print(out.read_text())
