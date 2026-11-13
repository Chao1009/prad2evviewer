#!/usr/bin/env python3
"""Histogram t_cluster - t_veto for the four veto counters V1..V4 (2x2 panels, one PNG).

t_cluster : seed_time_ns of the benchmark's result_clusters.tsv (--cluster-times):
            seed-module MF (matched-filter) hit time, RAW (the benchmark ran with --time-calib off,
            so time_calib.json is NOT applied; see tcluster_minus_tveto_scan.py for the calibrated variants).
t_veto    : WaveAnalyzer CFD peak time of the veto FADC channel, from the
            module_timing_dump.py npz of the same split (v_ev, v_ch, v_t, v_int).
            Per event and counter, the pulse with the largest integral is used.
The two files are joined on the EVIO event number.

  python3 tcluster_minus_tveto.py CLUSTERS_TSV VETO_NPZ [--config legacy] [--out PNG] [--bin NS] [--label TEXT]
"""
import argparse
import numpy as np
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
ap.add_argument("clusters_tsv")
ap.add_argument("veto_npz")
ap.add_argument("--config", default="legacy", help="config_name rows to use (legacy|gated)")
ap.add_argument("--out", default="figures/tcluster_minus_tveto.png")
ap.add_argument("--label", default="", help="sample label for title / table rows")
ap.add_argument("--bin", type=float, default=2.0, help="histogram bin width in ns (default 2.0)")
args = ap.parse_args()
bins = np.arange(-300, 100, args.bin)


def gauss_bg(x, a, mu, sig, c):
    return a * np.exp(-0.5 * ((x - mu) / sig) ** 2) + c


def fit_peak(x, h):
    """Gauss + constant within +-10 ns of the histogram mode; returns (a, mu, sigma, c per bin)."""
    mode = x[np.argmax(h)]
    w = abs(x - mode) <= 10
    bg = np.median(h)
    return curve_fit(gauss_bg, x[w], h[w], p0=(h.max() - bg, mode, 4.0, bg),
                     bounds=([0, -np.inf, 0.2, 0], np.inf))

# clusters: event number and seed time of the chosen config
cfg = np.genfromtxt(args.clusters_tsv, delimiter="\t", names=True, dtype=None,
                    encoding=None, usecols=("config_name", "evnum", "seed_time_ns"))
cfg = cfg[cfg["config_name"] == args.config]
ev, t_cl = cfg["evnum"], cfg["seed_time_ns"]

d = np.load(args.veto_npz)
v_evnum = d["ev_num"][d["v_ev"]]          # EVIO event number of each veto pulse

fig, axs = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
for k, ax in zip((1, 2, 3, 4), axs.flat):
    m = d["v_ch"] == k
    # largest-integral pulse per event for this counter
    o = np.lexsort((d["v_int"][m], v_evnum[m]))
    e_s, t_s = v_evnum[m][o], d["v_t"][m][o]
    last = np.r_[e_s[1:] != e_s[:-1], True]
    e_v, t_v = e_s[last], t_s[last]
    # join on event number
    i = np.searchsorted(e_v, ev)
    i[i == len(e_v)] = 0
    ok = e_v[i] == ev
    dt = t_cl[ok] - t_v[i[ok]]
    h, edges = np.histogram(dt, bins=bins)
    x = 0.5 * (edges[1:] + edges[:-1])
    ax.stairs(h, edges, fill=True, color="C0")
    ax.set_xlabel("t_cluster - t_veto  (ns)")
    ax.set_ylabel(f"events / {args.bin:g} ns")
    txt, row = f"N = {len(dt)}\nfit failed", f"| {args.label} | V{k} | {len(dt)} | fit failed | | |"
    try:
        p, _ = fit_peak(x, h)
        xs = np.linspace(p[1] - 10, p[1] + 10, 200)
        ax.plot(xs, gauss_bg(xs, *p), "r-")
        txt = f"N = {len(dt)}\n$\\mu$ = {p[1]:.2f} ns\n$\\sigma$ = {abs(p[2]):.2f} ns"
        row = f"| {args.label} | V{k} | {len(dt)} | {p[1]:.2f} | {abs(p[2]):.2f} | {p[3]:.1f} |"
    except Exception:
        pass
    ax.set_title(f"V{k}")
    ax.text(0.03, 0.95, txt, transform=ax.transAxes, va="top")
    print(row)
fig.suptitle(f"{args.label}  ({args.config}, bin {args.bin:g} ns)")
fig.tight_layout()
fig.savefig(args.out, dpi=130)
