#!/usr/bin/env python3
"""plot_module_resolution.py — module-by-module energy resolution after calibration.

For every PbWO4 seed module of HyCal, fit the elastic peak of the cluster
energy (Gaussian and Crystal Ball, same windows as the benchmark) before and
after the seed-indexed calibration, and map sigma_M/mu_M and the peak residual
mu_M/E'(theta_M) - 1 over the detector.  A cross-run summary ranks the worst
modules per run type.

Input: a benchmark scan directory (benchmark_hycal_timing_scan.py with
  --write-cluster-tsv) laid out as <scan>/scan_<run:06d>/split_<n:05d>/ with
  result.tsv, result_window.tsv, result_clusters.tsv.  Only splits with all
  three files are used.  The event pool reproduces the benchmark's calibration
  pool: chosen path, seed radius < --inner-radius, seed energy >= the run's
  min_seed_energy, outside the --exclude-rowcol rectangle, finite
  cluster_E_cal.  pre = cluster_E_raw, post = cluster_E_cal on the same events.

Output, <out>/<run>/ (names <run>_<what>):
  sigma_post.png sigma_pre.png   heat map of sigma_M/mu_M (%), common scale
  resid_post.png resid_pre.png   heat map of mu_M/E' - 1 (%)
  N.png                          heat map of pool events per module
  sigma_hist.png                 sigma_M/mu_M distribution + pooled sigma_E/E
  spectra.pdf                    per-module post spectra and fits (--no-pdf)
  modules.tsv                    all per-module numbers
With --group-runs the primary products are per run type, <out>/<slug>/<slug>_*
  (slug 3p49_H2, 2p24_H2, 0p73_H2, 2p24_C, Ta_X17; same files): the events of all
  runs of the type are pooled before the per-module fits, so N and the --min-n
  cut refer to the type total; the per-run products above are still written.
A module is fitted only if N > --min-n; 0 < N <= --min-n is "low_stat" (hatched).
<out>/summary/: runs_summary.tsv, types_summary.tsv (--group-runs), corr_2p24_H2_r.tsv + scatter grid,
  worst_modules.tsv (10 worst modules, W428/W494 and Moller-excluded modules per
  run type; from the pooled type tables with --group-runs).

Usage:
  python3 plot_module_resolution.py --scan-dir output/scan_20261023_pc_clusters \\
      --out output/module_resolution [--runs 24386 25320] [--path gated] \\
      [--fit {auto,gauss,cb}] [--min-n 1000] [--group-runs] [--splits 10 20] [--no-pdf]
"""

from __future__ import annotations

import argparse
import contextlib
import importlib.util
import io
import math
import os
import re
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Circle

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))

# Fit routines and peak shapes live in the benchmark; import them unchanged.
_spec = importlib.util.spec_from_file_location(
    "bench", _HERE / "benchmark_hycal_timing_test_20261021.py")
bench = importlib.util.module_from_spec(_spec)
sys.modules["bench"] = bench
_spec.loader.exec_module(bench)
import _common as C  # noqa: E402

RUN_TYPE = {**{r: "3.49 GeV H2" for r in (24386, 24503, 24512)},
            **{r: "2.24 GeV H2" for r in (25228, 25257, 25260, 25299, 25320)},
            **{r: "0.73 GeV H2" for r in (24590, 24603)},
            25308: "2.24 GeV C",
            **{r: "Ta X17" for r in (25800, 26138, 26230)}}
TYPES = ["3.49 GeV H2", "2.24 GeV H2", "0.73 GeV H2", "2.24 GeV C", "Ta X17"]
GAUSS_RUNS = {25228, 25257, 25260, 25299, 25320, 24590, 24603}
H2_224 = [25228, 25257, 25260, 25299, 25320]
NEEDED = ("result.tsv", "result_window.tsv", "result_clusters.tsv")
FIT_LABEL = {"g": "Gaussian", "cb": "Crystal Ball"}
SPECIAL = ("W428", "W494")
N_WORST = 10
PER_PAGE = 12


# ---- data loading -----------------------------------------------------------

def parse_header(path):
    txt = "".join(l for l in open(path) if l.startswith("#"))

    def num(pat):
        return float(re.search(pat, txt).group(1))
    return dict(kE=num(r"kE=([\d.eE+-]+)"), kS=num(r"kS=([\d.eE+-]+)"),
                ebeam=num(r"Ebeam=([\d.]+)"),
                min_seed=num(r"min_seed_energy=([\d.]+)"),
                bw=num(r"bin_width=([\d.]+)"),
                target_ref=num(r"target_mean_MeV=([\d.]+)"),
                target=re.search(r"kinematics target=(\w+)", txt).group(1))


def complete_splits(run_dir, wanted=None):
    return [d for d in sorted(run_dir.glob("split_*"))
            if (not wanted or int(d.name.split("_")[1]) in wanted)
            and all((d / f).exists() for f in NEEDED)]


def load_run(scan_dir, run, args):
    """Pool events of all complete splits of a run.  Returns dict or None."""
    splits = complete_splits(scan_dir / f"scan_{run:06d}", args.splits)
    if not splits:
        return None
    cfg = "legacy" if args.path == "legacy" else "gated"
    wpath = "legacy" if args.path == "legacy" else "new"
    r2, (r0, r1, c0, c1) = args.inner_radius ** 2, args.exclude_rowcol
    mods, post, pre, win = [], [], [], []
    for d in splits:
        h = parse_header(d / "result.tsv")
        df = pd.read_csv(d / "result_clusters.tsv", sep="\t")
        df = df[df.config_name == cfg]
        inrect = (df.seed_row + 1).between(r0, r1) & (df.seed_col + 1).between(c0, c1)  # 1-based
        df = df[(df.seed_x_mm ** 2 + df.seed_y_mm ** 2 < r2)
                & (df.seed_energy_MeV >= h["min_seed"]) & ~inrect
                & np.isfinite(df.cluster_E_cal_MeV)]
        mods.append(df.seed_module_id.to_numpy())
        post.append(df.cluster_E_cal_MeV.to_numpy())
        pre.append(df.cluster_E_raw_MeV.to_numpy())
        w = pd.read_csv(d / "result_window.tsv", sep="\t", comment="#")
        w = w[w.path == wpath].set_index("module")
        win.append(pd.DataFrame({"eprime": w.eprime_MeV,
                                 "moller": (w.exclude_reason == "moller_gap").astype(int),
                                 "nsplit": 1}))
    win = pd.concat(win).groupby(level=0).agg(eprime=("eprime", "first"),
                                              moller=("moller", "sum"),
                                              nsplit=("nsplit", "sum"))
    return dict(hdr=parse_header(splits[0] / "result.tsv"), mod=np.concatenate(mods),
                post=np.concatenate(post), pre=np.concatenate(pre), win=win,
                splits=splits)


def load_geometry(run):
    sys.stdout.flush()
    sys.stderr.flush()
    saved = os.dup(1), os.dup(2)    # setup_pipeline logs from Python and C++
    null = os.open(os.devnull, os.O_WRONLY)
    os.dup2(null, 1)
    os.dup2(null, 2)
    try:
        p = C.setup_pipeline(evio_path="", run_num=run)
    finally:
        sys.stdout.flush()
        sys.stderr.flush()
        os.dup2(saved[0], 1)
        os.dup2(saved[1], 2)
        for fd in (*saved, null):
            os.close(fd)
    rows = []
    for i in range(p.hycal.module_count()):
        m = p.hycal.module(i)
        if m.is_pwo4():
            rows.append((m.name, m.id, m.row, m.column, m.x, m.y, m.size_x, m.size_y))
    g = pd.DataFrame(rows, columns=["module", "module_id", "row", "col", "x", "y", "sx", "sy"])
    g["r"] = np.hypot(g.x, g.y)
    return g.set_index("module")


# ---- fitting ----------------------------------------------------------------

def hist_range(hdr):
    return hdr["target_ref"] - 800 * hdr["kE"], hdr["target_ref"] + 800 * hdr["kE"]


def ratio(s, m):
    return s / m if (np.isfinite(s) and np.isfinite(m) and m > 0) else float("nan")


def fit_one(vals, hdr):
    """Gaussian and Crystal-Ball fits as in the benchmark post-cal block."""
    lo, hi = hist_range(hdr)
    with contextlib.redirect_stderr(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gmu, gsig, gA, *_ = bench.fit_peak(vals, lo, hi, hdr["bw"], hdr["kS"])
        cmu, csig, ca, cn, cA, *_ = bench.fit_peak_cb(
            vals, lo, hi, hdr["bw"], mu0_hint=None if math.isnan(gmu) else gmu, escale=hdr["kS"])
    return dict(g_mu=gmu, g_sigma=gsig, g_A=gA, cb_mu=cmu, cb_sigma=csig,
                cb_alpha=ca, cb_n=cn, cb_A=cA)


def fit_modules(mod, energies, win, geo, hdr, min_n):
    """Per-module fits for each tag of `energies` ({tag: array}).  Returns a
    DataFrame over all PbWO4 modules."""
    order = np.argsort(mod, kind="stable")
    names, start = np.unique(mod[order], return_index=True)
    groups = dict(zip(names, np.split(order, start[1:])))
    df = geo.join(win, how="left")
    df["N"] = [len(groups.get(m, ())) for m in df.index]
    df["moller"] = df.moller.fillna(0).astype(int)
    df["nsplit"] = df.nsplit.fillna(0).astype(int)
    df["status"] = np.where(df.N == 0, "no_data", np.where(df.N <= min_n, "low_stat", "fitted"))
    cols = {}
    for tag, e in energies.items():
        res = {m: fit_one(e[groups[m]], hdr) for m in df.index[df.N > min_n]}
        for k in ("g_mu", "g_sigma", "g_A", "cb_mu", "cb_sigma", "cb_alpha", "cb_n", "cb_A"):
            cols[f"{tag}_{k}"] = [res[m][k] if m in res else float("nan") for m in df.index]
        for f in ("g", "cb"):
            cols[f"{tag}_{f}_s_over_mu"] = [ratio(cols[f"{tag}_{f}_sigma"][i], cols[f"{tag}_{f}_mu"][i])
                                            for i in range(len(df))]
    return pd.concat([df, pd.DataFrame(cols, index=df.index)], axis=1)


def select_fit(df, fit, tags, min_n):
    """Fill the preferred-fit columns <tag>_mu/_sigma/_s_over_mu/_resid."""
    f = "g" if fit == "gauss" else "cb"
    for tag in tags:
        df[f"{tag}_mu"] = df[f"{tag}_{f}_mu"]
        df[f"{tag}_sigma"] = df[f"{tag}_{f}_sigma"]
        df[f"{tag}_s_over_mu"] = df[f"{tag}_{f}_s_over_mu"]
        df[f"{tag}_resid"] = df[f"{tag}_mu"] / df.eprime - 1
    df.loc[(df.N > min_n) & df.post_s_over_mu.isna(), "status"] = "fit_fail"
    return df


# ---- plotting ---------------------------------------------------------------

def pct(a, p):
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    return float(np.percentile(a, p)) if a.size else float("nan")


def draw_map(geo, df, values, title, label, cmap, vmin, vmax, path, args):
    """Module-square heat map.  NaN values: hatched if low-stat, else grey outline."""
    fig, ax = plt.subplots(figsize=(7.4, 6.4))
    sq = lambda m: Rectangle((geo.x[m] - geo.sx[m] / 2, geo.y[m] - geo.sy[m] / 2), geo.sx[m], geo.sy[m])
    v = pd.Series(np.asarray(values, float), index=df.index)
    good = v.index[v.notna()]
    low = df.index[(df.status == "low_stat") & v.isna()]
    rest = df.index.difference(good).difference(low)
    ax.add_collection(PatchCollection([sq(m) for m in rest], facecolor="none", edgecolor="0.7", lw=0.4))
    ax.add_collection(PatchCollection([sq(m) for m in low], facecolor="0.92", edgecolor="0.45",
                                      lw=0.4, hatch="////"))
    pc = PatchCollection([sq(m) for m in good], cmap=cmap, edgecolor="0.3", lw=0.3)
    pc.set_array(v[good].to_numpy())
    pc.set_clim(vmin, vmax)
    ax.add_collection(pc)
    ex = geo[(geo.row + 1).between(*args.exclude_rowcol[:2])
             & (geo.col + 1).between(*args.exclude_rowcol[2:])]   # 1-based, as the benchmark
    x0, x1 = (ex.x - ex.sx / 2).min(), (ex.x + ex.sx / 2).max()
    y0, y1 = (ex.y - ex.sy / 2).min(), (ex.y + ex.sy / 2).max()
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, ec="crimson", lw=1.4, ls="--"))
    ax.add_patch(Circle((0, 0), args.inner_radius, fill=False, ec="black", lw=1.2, ls=":"))
    for nm in SPECIAL:
        ax.annotate(nm, (geo.x[nm], geo.y[nm]), fontsize=6, ha="center", va="center")
    ax.set(xlim=(-260, 260), ylim=(-260, 260), aspect="equal", xlabel="x (mm)", ylabel="y (mm)")
    fig.colorbar(pc, ax=ax, shrink=0.85).set_label(label)
    ax.set_title(title, fontsize=9)
    ax.text(0.01, 0.01, f"hatched: N<={args.min_n} | grey: no data / fit failed | dashed: excluded rect | "
            f"dotted: r={args.inner_radius:g} mm", transform=ax.transAxes, fontsize=6, va="bottom")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_hist(df, pooled_s, title, path):
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    post, pre = 100 * df.post_s_over_mu.dropna(), 100 * df.pre_s_over_mu.dropna()
    bins = np.linspace(min(pct(pre, 1), pct(post, 1)) * 0.9,
                       max(pct(pre, 99), pct(post, 99), 100 * pooled_s) * 1.1, 40)
    ax.hist(pre, bins, histtype="step", color="0.5", lw=1.2, label=f"pre (n={len(pre)})")
    ax.hist(post, bins, color="tab:blue", alpha=0.7, label=f"post (n={len(post)})")
    ax.axvline(post.median(), color="k", lw=1.2, label=f"post median {post.median():.2f}%")
    ax.axvspan(pct(post, 16), pct(post, 84), color="k", alpha=0.12,
               label=f"16-84%: {pct(post, 16):.2f}-{pct(post, 84):.2f}%")
    ax.axvline(100 * pooled_s, color="tab:red", ls="--", lw=1.4, label=f"pooled sigma/E {100 * pooled_s:.2f}%")
    ax.set(xlabel="sigma_M / mu_M (%)", ylabel="modules")
    ax.legend(fontsize=7)
    ax.set_title(title, fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_spectra(df, mod, post, hdr, fit, label, path, args):
    f = "g" if fit == "gauss" else "cb"
    order = np.argsort(mod, kind="stable")
    names, start = np.unique(mod[order], return_index=True)
    groups = dict(zip(names, np.split(order, start[1:])))
    ids = df[df.N > 0].sort_values("module_id").index
    lo, hi = hist_range(hdr)
    edges = np.arange(lo, hi + hdr["bw"], hdr["bw"])
    cen = 0.5 * (edges[:-1] + edges[1:])
    xl = (hdr["target_ref"] - 500 * hdr["kE"], hdr["target_ref"] + 300 * hdr["kE"])
    n = PER_PAGE
    with PdfPages(path) as pdf, np.errstate(all="ignore"):
        for p0 in range(0, len(ids), n):
            fig, axs = plt.subplots(3, 4, figsize=(14, 9))
            page = ids[p0:p0 + n]
            for ax, m in zip(axs.flat, page):
                r = df.loc[m]
                ax.stairs(np.histogram(post[groups[m]], edges)[0], edges, color="tab:blue")
                txt = f"{m} (id {r.module_id})  N={int(r.N)}"
                if r.N <= args.min_n:
                    txt += "\nlow statistics, no fit"
                elif np.isfinite(r.post_sigma):
                    y = (bench.gauss(cen, r.post_g_A, r.post_g_mu, r.post_g_sigma) if f == "g" else
                         bench.crystal_ball(cen, r.post_cb_A, r.post_cb_mu, r.post_cb_sigma,
                                            r.post_cb_alpha, r.post_cb_n))
                    ax.plot(cen, y, color="tab:red", lw=1)
                    txt += f"\nmu={r.post_mu:.0f} sig={r.post_sigma:.0f} ({100 * r.post_s_over_mu:.2f}%)"
                else:
                    txt += "\nfit failed"
                ax.set_title(txt, fontsize=7, loc="left")
                ax.axvline(hdr["target_ref"], color="0.6", lw=0.6, ls=":")
                ax.set_xlim(*xl)
                ax.tick_params(labelsize=6)
            for ax in list(axs.flat)[len(page):]:
                ax.axis("off")
            fig.suptitle(f"{label}  post-cal spectra  ({FIT_LABEL[f]} fit)  page {p0 // n + 1}", fontsize=9)
            fig.tight_layout(rect=(0, 0, 1, 0.97))
            pdf.savefig(fig)
            plt.close(fig)


def corr_grid(tabs, out):
    runs = list(tabs)
    n = len(runs)
    r = pd.DataFrame(np.nan, index=runs, columns=runs)
    fig, axs = plt.subplots(n, n, figsize=(2.6 * n, 2.6 * n), squeeze=False)
    allv = np.concatenate([100 * t.post_s_over_mu.dropna() for t in tabs.values()])
    lim = (pct(allv, 1) * 0.95, pct(allv, 99) * 1.05)
    for i, a in enumerate(runs):
        for j, b in enumerate(runs):
            ax = axs[i, j]
            ax.tick_params(labelsize=6)
            if i == j:
                r.loc[a, b] = 1.0
                ax.text(0.5, 0.5, str(a), ha="center", va="center", transform=ax.transAxes)
                ax.set_xticks([])
                ax.set_yticks([])
                continue
            both = pd.concat([tabs[b].post_s_over_mu, tabs[a].post_s_over_mu], axis=1).dropna() * 100
            r.loc[a, b] = np.corrcoef(both.iloc[:, 0], both.iloc[:, 1])[0, 1] if len(both) > 2 else np.nan
            ax.plot(both.iloc[:, 0], both.iloc[:, 1], ".", ms=2, alpha=0.5)
            ax.plot(lim, lim, color="0.6", lw=0.6)
            ax.set(xlim=lim, ylim=lim)
            ax.text(0.04, 0.92, f"r={r.loc[a, b]:.2f}  n={len(both)}", transform=ax.transAxes,
                    fontsize=7, va="top")
            if i == n - 1:
                ax.set_xlabel(f"{b}  sigma/mu (%)", fontsize=7)
            if j == 0:
                ax.set_ylabel(f"{a}  sigma/mu (%)", fontsize=7)
    fig.suptitle("2.24 GeV H2: per-module sigma_M/mu_M, run vs run (post-cal)", fontsize=10)
    fig.tight_layout()
    fig.savefig(out / "corr_2p24_H2_scatter.png", dpi=110)
    plt.close(fig)
    r.round(4).to_csv(out / "corr_2p24_H2_r.tsv", sep="\t", index_label="run")
    return r


# ---- per-run and summary ----------------------------------------------------

SLUG = {"3.49 GeV H2": "3p49_H2", "2.24 GeV H2": "2p24_H2", "0.73 GeV H2": "0p73_H2",
        "2.24 GeV C": "2p24_C", "Ta X17": "Ta_X17"}


def pool_runs(rs):
    """Concatenate the events of several runs (same type) into one data dict;
    target_ref = mean over runs, window table grouped by module."""
    hdr = dict(rs[0]["d"]["hdr"])
    hdr["target_ref"] = np.mean([r["d"]["hdr"]["target_ref"] for r in rs])
    w = pd.concat([r["d"]["win"] for r in rs]).groupby(level=0).agg(
        eprime=("eprime", "mean"), moller=("moller", "sum"), nsplit=("nsplit", "sum"))
    cat = lambda k: np.concatenate([r["d"][k] for r in rs])
    return dict(hdr=hdr, mod=cat("mod"), post=cat("post"), pre=cat("pre"), win=w,
                splits=[s for r in rs for s in r["d"]["splits"]])


def analyze(d, name, label, fit, geo, args):
    """Fits, maps and tables of one data dict `d`; products go to <out>/<name>/<name>_*.
    Returns (module table, summary-statistics dict)."""
    t0 = time.time()
    hdr = d["hdr"]
    f = "g" if fit == "gauss" else "cb"
    # the TSV seed row/col must follow the module objects' row/column
    chk = pd.read_csv(d["splits"][0] / "result_clusters.tsv", sep="\t",
                      usecols=["seed_module_id", "seed_row", "seed_col"]).drop_duplicates()
    chk = chk[chk.seed_module_id.isin(geo.index)]
    bad = int(((geo.row.reindex(chk.seed_module_id).to_numpy() != chk.seed_row.to_numpy())
               | (geo.col.reindex(chk.seed_module_id).to_numpy() != chk.seed_col.to_numpy())).sum())
    df = fit_modules(d["mod"], {"post": d["post"], "pre": d["pre"]}, d["win"], geo, hdr, args.min_n)
    df = select_fit(df, fit, ("post", "pre"), args.min_n)
    pp = {t: fit_one(d[t], hdr) for t in ("post", "pre")}
    pool = {(t, k): ratio(pp[t][f"{k}_sigma"], pp[t][f"{k}_mu"]) for t in pp for k in ("g", "cb")}
    ps = {t: pool[(t, f)] for t in pp}
    print(f"{label}: {len(d['splits'])} splits, {len(d['mod'])} pool events; pooled sigma/E "
          f"post G {100 * pool[('post', 'g')]:.2f}% CB {100 * pool[('post', 'cb')]:.2f}% | "
          f"pre G {100 * pool[('pre', 'g')]:.2f}% CB {100 * pool[('pre', 'cb')]:.2f}%  "
          f"(fit={fit}); row/col mismatches vs geometry: {bad}")
    rdir = args.out / str(name)
    rdir.mkdir(parents=True, exist_ok=True)
    fn = lambda w: str(rdir / f"{name}_{w}")
    tag = f"{label} | {hdr['ebeam'] / 1000:.2f} GeV {hdr['target']} | {args.path} | {FIT_LABEL[f]}"
    s_all = 100 * np.concatenate([df.post_s_over_mu.dropna(), df.pre_s_over_mu.dropna()])
    smin, smax = pct(s_all, 1), pct(s_all, 99)
    rmax = max(pct(100 * df.post_resid.abs(), 98), pct(100 * df.pre_resid.abs(), 98), 0.1)
    if np.isfinite(smax):
        for t in ("post", "pre"):
            draw_map(geo, df, 100 * df[f"{t}_s_over_mu"], f"{tag} | {t}-cal sigma_M/mu_M",
                     "sigma_M / mu_M (%)", "viridis", smin, smax, fn(f"sigma_{t}.png"), args)
            draw_map(geo, df, 100 * df[f"{t}_resid"], f"{tag} | {t}-cal peak residual mu_M/E' - 1",
                     "mu_M / E' - 1 (%)", "RdBu_r", -rmax, rmax, fn(f"resid_{t}.png"), args)
        plot_hist(df, ps["post"], f"{tag} | per-module resolution", fn("sigma_hist.png"))
    nvals = df.N.where(df.N > 0)
    draw_map(geo, df, nvals, f"{tag} | pool events per module", "N events", "magma_r",
             0, pct(nvals, 98), fn("N.png"), args)
    if not args.no_pdf:
        plot_spectra(df, d["mod"], d["post"], hdr, fit, label, fn("spectra.pdf"), args)
    df.reset_index().to_csv(fn("modules.tsv"), sep="\t", index=False, float_format="%.5g")
    post, pre = 100 * df.post_s_over_mu.dropna(), 100 * df.pre_s_over_mu.dropna()
    nd = int((df.N > 0).sum())
    st = dict(n_with_data=nd, n_above_min=int((df.N > args.min_n).sum()), n_fitted=len(post),
              n_low_stat=int((df.status == "low_stat").sum()),
              n_fit_fail=int((df.status == "fit_fail").sum()), frac_fitted=len(post) / max(nd, 1),
              post_median=post.median(), post_p16=pct(post, 16), post_p84=pct(post, 84),
              pre_median=pre.median(), pre_p16=pct(pre, 16), pre_p84=pct(pre, 84),
              pooled_post=100 * ps["post"], pooled_pre=100 * ps["pre"])
    print(f"  {st['n_fitted']} fitted, {st['n_low_stat']} low-stat, {st['n_fit_fail']} failed; "
          f"per-module post median {st['post_median']:.2f}% "
          f"[{st['post_p16']:.2f}, {st['post_p84']:.2f}]  ({time.time() - t0:.0f}s)")
    return df, st


def process_run(scan_dir, run, args):
    d = load_run(scan_dir, run, args)
    if d is None:
        print(f"[skip] run {run}: no complete split")
        return None
    hdr = d["hdr"]
    fit = args.fit if args.fit != "auto" else ("gauss" if run in GAUSS_RUNS else "cb")
    geo = load_geometry(run)
    df, st = analyze(d, run, f"run {run}", fit, geo, args)
    row = dict(run=run, type=RUN_TYPE.get(run, f"{hdr['ebeam'] / 1000:.2f} GeV {hdr['target']}"),
               beam_MeV=hdr["ebeam"], fit=fit, n_splits=len(d["splits"]), n_pool=len(d["mod"]), **st)
    return dict(row=row, df=df, geo=geo, d=d)


def process_type(typ, rs, args):
    """Pool the events of all runs of a type before the per-module fits."""
    d = pool_runs(rs)
    fit = rs[0]["row"]["fit"]
    df, st = analyze(d, SLUG[typ], f"{typ} ({len(rs)} runs)", fit, rs[0]["geo"], args)
    row = dict(type=typ, runs=",".join(str(r["row"]["run"]) for r in rs), fit=fit,
               n_runs=len(rs), n_pool=len(d["mod"]), **st)
    return dict(row=row, df=df)


def worst_table(results, args, out, pooled=None):
    """pooled: {type: module table} already fitted on the pooled runs (--group-runs)."""
    pooled = pooled or {}
    rows = []
    for typ in TYPES:
        rs = [r for r in results if r["row"]["type"] == typ]
        if not rs:
            continue
        fit = rs[0]["row"]["fit"]
        if typ in pooled:
            df = pooled[typ]
        elif len(rs) == 1:
            df = rs[0]["df"]
        else:   # pool the events of all runs of the type and refit
            d = pool_runs(rs)
            df = fit_modules(d["mod"], {"post": d["post"]}, d["win"], rs[0]["geo"], d["hdr"], args.min_n)
            df = select_fit(df, fit, ("post",), args.min_n)
        runs = ",".join(str(r["row"]["run"]) for r in rs)
        common = dict(type=typ, runs=runs)
        ok = df[df.post_s_over_mu.notna()].sort_values("post_s_over_mu", ascending=False)
        sel = [("worst", rank, m, r) for rank, (m, r) in enumerate(ok.head(N_WORST).iterrows(), 1)]
        sel += [("watch", 0, m, df.loc[m]) for m in SPECIAL if m in df.index]
        for kind, rank, m, r in sel:
            rows.append(dict(common, kind=kind, rank=rank, module=m, module_id=int(r.module_id),
                             sigma_over_mu_pct=100 * r.post_s_over_mu, sigma_MeV=r.post_sigma,
                             mu_MeV=r.post_mu, N=int(r.N), x_mm=r.x, y_mm=r.y, r_mm=r.r,
                             resid_pct=100 * r.post_resid, special=m in SPECIAL,
                             moller_splits=int(r.moller), n_splits=int(r.nsplit)))
        for m, r in df[df.moller > 0].sort_values("module_id").iterrows():
            rows.append(dict(common, kind="moller_excluded", rank=0, module=m, module_id=int(r.module_id),
                             N=int(r.N), x_mm=r.x, y_mm=r.y, r_mm=r.r, special=m in SPECIAL,
                             moller_splits=int(r.moller), n_splits=int(r.nsplit)))
    w = pd.DataFrame(rows)
    w.to_csv(out / "worst_modules.tsv", sep="\t", index=False, float_format="%.5g")
    for typ, g in w.groupby("type", sort=False):
        print(f"\n== {typ} (runs {g.runs.iloc[0]}): worst {N_WORST} by sigma_M/mu_M (post), then W428/W494 (rank 0)")
        print("  rank module    id  sig/mu%   sig     mu      N   r_mm   res%  flags")
        for _, r in g[g.kind.isin(("worst", "watch"))].iterrows():
            fl = ("W428/W494 " if r.special else "") + \
                 (f"moller {r.moller_splits}/{r.n_splits}" if r.moller_splits else "")
            print(f"  {r['rank']:>4d} {r.module:<6} {r.module_id:>5d} {r.sigma_over_mu_pct:7.2f} "
                  f"{r.sigma_MeV:5.0f} {r.mu_MeV:7.0f} {r.N:>6d} {r.r_mm:6.1f} {r.resid_pct:+6.2f}  {fl}")
        mo = g[g.kind == "moller_excluded"]
        print("  Moller-excluded: " + (", ".join(f"{r.module}({r.moller_splits}/{r.n_splits})"
                                                for _, r in mo.iterrows()) or "none"))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scan-dir", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--runs", type=int, nargs="+", help="default: every scan_* subdir")
    ap.add_argument("--path", choices=("legacy", "gated"), default="legacy")
    ap.add_argument("--min-n", type=int, default=1000, help="fit only modules with N > min-n pool events")
    ap.add_argument("--group-runs", action="store_true",
                    help="primary output: one map set per run type, events of its runs pooled before the fits")
    ap.add_argument("--fit", choices=("auto", "gauss", "cb"), default="auto",
                    help="auto: CB for 3.49 GeV/C/Ta runs, Gaussian for 2.24 and 0.73 GeV H2")
    ap.add_argument("--no-pdf", action="store_true")
    ap.add_argument("--inner-radius", type=float, default=200.0, help="pool seed radius (mm)")
    ap.add_argument("--exclude-rowcol", type=int, nargs=4, default=(15, 20, 15, 20),
                    metavar=("ROW0", "ROW1", "COL0", "COL1"),
                    help="excluded seed row/col rectangle, 1-based inclusive as the benchmark's "
                         "--exclude-rowcol (default 15 20 15 20: two rings around the hole at 17-18; "
                         "scans up to scan_20261023_pc* used the 0-based reading, i.e. 16 21 16 21 here)")
    ap.add_argument("--splits", type=int, nargs="+", help="restrict to these split numbers")
    args = ap.parse_args(argv)
    if args.path == "gated":   # keep the primary (legacy) outputs intact
        args.out = args.out.with_name(args.out.name + "_gated")

    runs = args.runs or sorted(int(p.name.split("_")[1]) for p in args.scan_dir.glob("scan_*") if p.is_dir())
    args.out.mkdir(parents=True, exist_ok=True)
    results = [r for r in (process_run(args.scan_dir, run, args) for run in runs) if r]
    if not results:
        sys.exit("no complete splits found")
    sdir = args.out / "summary"
    sdir.mkdir(exist_ok=True)
    pd.DataFrame([r["row"] for r in results]).to_csv(sdir / "runs_summary.tsv", sep="\t", index=False,
                                                     float_format="%.4f")
    tabs = {r["row"]["run"]: r["df"] for r in results if r["row"]["run"] in H2_224}
    if len(tabs) >= 2:
        print("\n2.24 GeV H2 per-module sigma/mu correlation (Pearson r):")
        print(corr_grid(tabs, sdir).round(2).to_string())
    pooled = {}
    if args.group_runs:
        trs = [(typ, process_type(typ, rs, args)) for typ in TYPES
               if (rs := [r for r in results if r["row"]["type"] == typ])]
        if trs:
            pd.DataFrame([t["row"] for _, t in trs]).to_csv(sdir / "types_summary.tsv", sep="\t",
                                                           index=False, float_format="%.4f")
            pooled = {typ: t["df"] for typ, t in trs}
    worst_table(results, args, sdir, pooled)


if __name__ == "__main__":
    main()
