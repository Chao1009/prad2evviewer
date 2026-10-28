#!/usr/bin/env python3
"""plot_module_energy.py — per-module calibrated cluster-energy spectra, pooled per run type.

For every PbWO4 module, one PDF page with four panels (fixed order: 0.73 GeV H2,
2.24 GeV H2, 3.49 GeV H2, X17 Ta): the calibrated cluster-energy spectrum of the
clusters seeded in that module, events of all runs of the type pooled, with the
preferred fit overlaid (Crystal Ball for 3.49 GeV H2 and X17 Ta, iterative
+-2 sigma Gaussian for 2.24 and 0.73 GeV H2) and E'(theta_M) marked.  A cover page
gives the per-module sigma/mu summary of each run type.

Input and pool definition are those of plot_module_resolution.py (imported, not
duplicated): a benchmark scan directory <scan>/scan_<run:06d>/split_<n:05d>/ with
result.tsv, result_window.tsv, result_clusters.tsv (benchmark_hycal_timing_scan.py
--write-cluster-tsv); legacy path, seed r < --inner-radius, seed E >= min_seed,
outside the --exclude-rowcol rectangle, finite cluster_E_cal; energy = cluster_E_cal.
A module is fitted only if N > --min-n (N = pool events of the run type).  A module
gets a page if N > --min-n in at least one run type.

Output, <out>/: module_energy.pdf (cover + one page per module, by module id),
  module_energy.tsv (module, id, x, y; per run type N, E' , mu, sigma, s_over_mu, resid, fit).

Usage:
  python3 plot_module_energy.py [--scan-dir output/scan_20261025/pos_all]
      [--out output/module_energy_20261025] [--min-n 1000] [--jobs 4]
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
import plot_module_resolution as pm  # noqa: E402

bench = pm.bench

# fixed panel order, run list and preferred fit per run type (carbon 25308 and 24578 excluded)
PANELS = ["0.73 GeV H2", "2.24 GeV H2", "3.49 GeV H2", "Ta X17"]
FIT = {"0.73 GeV H2": "gauss", "2.24 GeV H2": "gauss", "3.49 GeV H2": "cb", "Ta X17": "cb"}
RUNS = {t: sorted(r for r, tt in pm.RUN_TYPE.items() if tt == t) for t in PANELS}
SLUG = {"0.73 GeV H2": "0p73_H2", "2.24 GeV H2": "2p24_H2", "3.49 GeV H2": "3p49_H2", "Ta X17": "Ta_X17"}
COLOR = {"0.73 GeV H2": "tab:green", "2.24 GeV H2": "tab:orange", "3.49 GeV H2": "tab:blue",
         "Ta X17": "tab:purple"}
FIT_NAME = {"gauss": "Gaussian (+-2sigma)", "cb": "Crystal Ball"}


def do_type(typ, scan_dir, geo, args):
    """Pool the runs of a type, fit every module with N > min-n, histogram every module.
    Returns only light-weight products (no event arrays)."""
    rs = []
    for run in RUNS[typ]:
        d = pm.load_run(scan_dir, run, args)
        if d is None:
            sys.exit(f"run {run}: no complete split in {scan_dir}")
        rs.append({"d": d})
    d = pm.pool_runs(rs)
    hdr, fit = d["hdr"], FIT[typ]
    df = pm.fit_modules(d["mod"], {"post": d["post"]}, d["win"], geo, hdr, args.min_n)
    df = pm.select_fit(df, fit, ("post",), args.min_n)
    pooled = pm.fit_one(d["post"], hdr)
    f = "g" if fit == "gauss" else "cb"
    lo, hi = pm.hist_range(hdr)
    edges = np.arange(lo, hi + hdr["bw"], hdr["bw"])
    order = np.argsort(d["mod"], kind="stable")
    names, start = np.unique(d["mod"][order], return_index=True)
    counts = {m: np.histogram(d["post"][ix], edges)[0]
              for m, ix in zip(names, np.split(order, start[1:])) if m in df.index}
    return dict(typ=typ, df=df, hdr=hdr, edges=edges, counts=counts, pooled=pooled,
                pooled_s=pm.ratio(pooled[f"{f}_sigma"], pooled[f"{f}_mu"]),
                pooled_mu=pooled[f"{f}_mu"], n_pool=len(d["mod"]), runs=RUNS[typ])


def _worker(a):
    return do_type(*a)


# ---- plotting ----------------------------------------------------------------

def panel(ax, res, m, args):
    typ, hdr, df = res["typ"], res["hdr"], res["df"]
    r = df.loc[m]
    edges = res["edges"]
    cen = 0.5 * (edges[:-1] + edges[1:])
    fit = FIT[typ]
    f = "g" if fit == "gauss" else "cb"
    cnt = res["counts"].get(m, np.zeros(len(cen)))
    xl = (hdr["target_ref"] - 500 * hdr["kE"], hdr["target_ref"] + 300 * hdr["kE"])
    ax.stairs(cnt, edges, color=COLOR[typ], lw=1.0)
    ax.set_title(f"{typ}  ({', '.join(str(x) for x in res['runs'])})", fontsize=8, loc="left")
    ax.tick_params(labelsize=7)
    ax.set_xlabel("calibrated cluster energy (MeV)", fontsize=7)
    ax.set_ylabel(f"clusters / {hdr['bw']:.3g} MeV", fontsize=7)
    ymax = max(cnt.max(), 1) if len(cnt) else 1
    txt = f"N = {int(r.N)}"
    if np.isfinite(r.eprime):
        ax.axvline(r.eprime, color="0.35", lw=0.9, ls="--", label="E'")
        xl = (min(xl[0], r.eprime - 100 * hdr["kE"]), max(xl[1], r.eprime + 100 * hdr["kE"]))
    if r.N <= args.min_n:
        txt = f"N={int(r.N)}, below --min-n ({args.min_n}), no fit"
    elif np.isfinite(r.post_sigma):
        y = (bench.gauss(cen, r.post_g_A, r.post_g_mu, r.post_g_sigma) if f == "g" else
             bench.crystal_ball(cen, r.post_cb_A, r.post_cb_mu, r.post_cb_sigma,
                                r.post_cb_alpha, r.post_cb_n))
        with np.errstate(all="ignore"):
            ax.plot(cen, y, color="tab:red", lw=1.3)
        txt = (f"{FIT_NAME[fit]} fit\nmu = {r.post_mu:.1f} MeV\nsigma = {r.post_sigma:.1f} MeV\n"
               f"sigma/mu = {100 * r.post_s_over_mu:.2f} %\nmu/E' - 1 = {100 * r.post_resid:+.2f} %\n"
               f"N = {int(r.N)}")
        ax.axvline(r.post_mu, color="tab:red", lw=0.6, ls=":")
    else:
        txt = f"N = {int(r.N)}\nfit failed"
    ax.set_xlim(*xl)
    ax.set_ylim(0, ymax * 1.18)
    ax.text(0.02, 0.97, txt, transform=ax.transAxes, fontsize=7, va="top", ha="left",
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85))
    if np.isfinite(r.eprime):
        ax.text(0.98, 0.97, f"E' = {r.eprime:.1f} MeV", transform=ax.transAxes, fontsize=7,
                va="top", ha="right", color="0.35")


def cover(results, args, summary):
    fig = plt.figure(figsize=(11, 8.5))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.1], hspace=0.28,
                          left=0.07, right=0.97, top=0.95, bottom=0.07)
    ax0 = fig.add_subplot(gs[0])
    ax0.axis("off")
    ax0.set_title("Per-module calibrated cluster-energy spectra, pooled per run type", fontsize=12,
                  loc="left", fontweight="bold")
    cfg = (f"Configuration R2: spline templates + ampcorr, per-crystal gains output/percrystal_spl,\n"
           f"position correction output/poscorr_20261025, centred rectangle excluded (rows/cols "
           f"{' '.join(str(x) for x in args.exclude_rowcol)}, 1-based), legacy path.\n"
           f"Pool: clusters seeded in a PbWO4 module with r < {args.inner_radius:g} mm, seed E >= min_seed, outside the\n"
           f"rectangle, finite calibrated energy; runs of a type pooled.  Scan: {args.scan_dir}.\n"
           f"Fit: Crystal Ball for 3.49 GeV H2 and Ta X17, iterative +-2 sigma Gaussian for 2.24 and 0.73 GeV H2;\n"
           f"fitted only if N > --min-n = {args.min_n}.  Carbon 25308 and run 24578 excluded.  "
           f"E' = elastic energy at the module angle; dashed line.")
    ax0.text(0, 0.97, cfg, transform=ax0.transAxes, fontsize=8, va="top", family="monospace")
    cols = ["run type", "runs", "fit", "pool events", "modules fitted", "median s/mu (%)",
            "16-84% (%)", "pooled s/E (%)", "pooled mu (MeV)"]
    rows = [[s["typ"], ",".join(str(x) for x in s["runs"][:3]) + ("\n" + ",".join(str(x) for x in s["runs"][3:]) if len(s["runs"]) > 3 else ""), s["fit"], f"{s['n_pool']:,}",
             str(s["n_fitted"]), f"{s['median']:.2f}", f"{s['p16']:.2f}-{s['p84']:.2f}",
             f"{s['pooled_s']:.2f}", f"{s['pooled_mu']:.1f}"] for s in summary]
    tb = ax0.table(cellText=rows, colLabels=cols, loc="lower center", cellLoc="center",
                   bbox=[0, 0.0, 1, 0.38], colWidths=[0.11, 0.2, 0.06, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
    tb.auto_set_font_size(False)
    tb.set_fontsize(7.5)
    ax = fig.add_subplot(gs[1])
    allv = np.concatenate([100 * r["df"].post_s_over_mu.dropna() for r in results.values()])
    bins = np.linspace(max(0, np.percentile(allv, 0.5) * 0.9), np.percentile(allv, 99.5) * 1.1, 36)
    for s in summary:
        v = 100 * results[s["typ"]]["df"].post_s_over_mu.dropna()
        ax.hist(v, bins, histtype="step", lw=1.6, color=COLOR[s["typ"]],
                label=f"{s['typ']} ({s['fit']}): median {s['median']:.2f}%, n={len(v)}")
        ax.axvline(s["pooled_s"], color=COLOR[s["typ"]], ls="--", lw=1.0)
    ax.set(xlabel="per-module sigma / mu (%)   (dashed: pooled sigma/E of the type)", ylabel="modules")
    ax.legend(fontsize=8)
    return fig


# ---- main --------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scan-dir", type=Path, default=Path("output/scan_20261025/pos_all"))
    ap.add_argument("--out", type=Path, default=Path("output/module_energy_20261025"))
    ap.add_argument("--min-n", type=int, default=1000, help="fit only modules with N > min-n pool events")
    ap.add_argument("--jobs", type=int, default=4, help="run types processed in parallel (1 = serial)")
    args = ap.parse_args(argv)
    # the settings load_run needs, fixed to the R2 benchmark pool definition
    args.path, args.inner_radius, args.exclude_rowcol, args.splits = "legacy", 200.0, (15, 20, 15, 20), None
    args.no_pdf = True
    assert pm.GAUSS_WINDOW == "2sigma", pm.GAUSS_WINDOW
    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    geo = pm.load_geometry(RUNS[PANELS[0]][0])     # before forking: pipeline setup is not fork-safe
    jobs = [(t, args.scan_dir, geo, args) for t in PANELS]
    if args.jobs > 1:
        with mp.get_context("fork").Pool(min(args.jobs, len(jobs))) as pool:
            out = pool.map(_worker, jobs, chunksize=1)
    else:
        out = [_worker(j) for j in jobs]
    results = {r["typ"]: r for r in out}
    print(f"fits done ({time.time() - t0:.0f}s)")

    summary = []
    for t in PANELS:
        r = results[t]
        v = 100 * r["df"].post_s_over_mu.dropna()
        summary.append(dict(typ=t, runs=r["runs"], fit=FIT[t], n_pool=r["n_pool"], n_fitted=len(v),
                            median=v.median(), p16=pm.pct(v, 16), p84=pm.pct(v, 84),
                            pooled_s=100 * r["pooled_s"], pooled_mu=r["pooled_mu"]))
        print(f"{t}: {r['n_pool']} events, {len(v)} modules fitted, median sigma/mu {v.median():.2f}% "
              f"[{pm.pct(v, 16):.2f}, {pm.pct(v, 84):.2f}], pooled sigma/E {100 * r['pooled_s']:.2f}%")

    # modules with a page: N > min-n in at least one run type
    ids = geo.module_id
    keep = np.zeros(len(geo), bool)
    for t in PANELS:
        keep |= (results[t]["df"].N.reindex(geo.index) > args.min_n).to_numpy()
    mods = geo[keep].sort_values("module_id").index

    # table
    tab = geo.loc[mods, ["module_id", "x", "y"]].rename(columns={"module_id": "id"})
    for t in PANELS:
        d = results[t]["df"].loc[mods]
        s = SLUG[t]
        tab[f"{s}_N"] = d.N
        tab[f"{s}_eprime"] = d.eprime
        for k in ("mu", "sigma", "s_over_mu", "resid"):
            tab[f"{s}_{k}"] = d[f"post_{k}"]
        tab[f"{s}_fit"] = np.where(d.N > args.min_n, FIT[t], "")
    tab.reset_index().to_csv(args.out / "module_energy.tsv", sep="\t", index=False, float_format="%.5g")

    with PdfPages(args.out / "module_energy.pdf") as pdf, np.errstate(all="ignore"):
        fig = cover(results, args, summary)
        pdf.savefig(fig)
        plt.close(fig)
        for m in mods:
            g = geo.loc[m]
            fig, axs = plt.subplots(2, 2, figsize=(11, 8.5))
            for ax, t in zip(axs.flat, PANELS):
                panel(ax, results[t], m, args)
            fig.suptitle(f"{m}   (id {int(g.module_id)})   x = {g.x:.1f} mm, y = {g.y:.1f} mm", fontsize=11,
                         fontweight="bold")
            fig.tight_layout(rect=(0, 0, 1, 0.96))
            pdf.savefig(fig)
            plt.close(fig)
    print(f"{len(mods)} module pages + cover -> {args.out / 'module_energy.pdf'}  ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
