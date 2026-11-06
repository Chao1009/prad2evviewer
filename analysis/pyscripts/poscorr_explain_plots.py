#!/usr/bin/env python3
"""poscorr_explain_plots.py - explanatory figures for the position ("bowl") correction, one run type.

  python3 poscorr_explain_plots.py --run-type 3p49GeV_H2 --out notes/figures/poscorr_20261105

Run from analysis/pyscripts after `source setup_env.sh`.  Read-only on database/, runinfo/, output/.
Only 3p49GeV_H2 has the form-study inputs (fig06); other slugs would need --chain/--form-dir edits.

Inputs (defaults relative to analysis/pyscripts):
  --chain output/r5ri_20261104
      poscorr/scan/scan_<run>/split_*/{result.tsv,result_window.tsv,result_clusters.tsv}   uncorrected map-scan (gains, no --pos-corr)
      poscorr/maps/<slug>/poscorr_r5ri_bilinear_<slug>_all.json                             poscorr_v1 bilinear map (12x12)
      poscorr/maps/<slug>/poscorr_r5ri_<slug>_all.json                                      poscorr_hybrid_v1 (poly8 + embedded bilinear)
      poscorr/maps/<slug>/poscorr_r5ri_<slug>_all_modules.tsv                               per-module n, own flag, c0..c7
      scan/splits/scan_024386/split_00010/result.tsv                                        final scan (hybrid applied), same gains/events
  --form-dir output/poscorr_form_20261105
      tables/3p49_heldout_per_module.tsv, tables/3p49_sigma.md
  database/hycal_map.json (via plot_module_resolution.load_geometry) for module geometry.
Imports position_corr.py (loader, grid_estimates, trunc_mean, poly_R) unchanged.

Figures: fig01_definition, fig02_bowl_profile, fig02b_bilinear_map, fig03_module_poly8,
fig04_coeff_layout, fig04b_coeff_hist, fig05_spectrum, fig06_form_comparison (+ facts_data.json).
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
from matplotlib.colors import Normalize

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
import position_corr as P   # noqa: E402

bench, G, R = P.bench, P.G, P.R
DPI = 170
CL = 0.6
plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.25, "axes.spines.top": False,
                     "axes.spines.right": False})
COL = {"unc": "#c0392b", "cor": "#1f5fa8", "bil": "#7f8c8d", "even6": "#e69f00", "even4": "#2a9d8f"}
FACTS = {}


def save(fig, out, name):
    p = out / name
    fig.savefig(p, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("wrote", p)


# ---- data ------------------------------------------------------------------

def load_all(chain, slug):
    args = argparse.Namespace(splits="all", path="legacy", keep_moller=False, inner_radius=200.0,
                              exclude_rowcol=(15, 20, 15, 20))
    scan = chain / "poscorr/scan"
    runs = [r for r in sorted(int(p.name.split("_")[1]) for p in scan.glob("scan_*") if p.is_dir())
            if G.runtype_slug(r) == slug]
    geo = R.load_geometry(runs[0])
    d = P.pool_type([P.load_pool(scan, r, args, geo) for r in runs])
    m = np.isfinite(d["r"]) & (d["r"] >= P.HYB_RMIN) & (d["r"] <= P.HYB_RMAX)
    d = {k: (v[m] if isinstance(v, np.ndarray) else v) for k, v in d.items()}
    return d, geo, runs, scan


def tmean(r):
    t = P.trunc_mean(r)
    if t is None:
        return np.nan, np.nan, 0
    mu, sig, nin, _ = t
    return mu, sig / np.sqrt(nin), len(r)


def hyb_R(hyb, u, v, seed):
    return bench.poscorr_eval(hyb, u, v, seed)[0]


# ---- fig01 -----------------------------------------------------------------

def fig01(out):
    fig, ax = plt.subplots(figsize=(7.6, 5.4))
    a = 20.77
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            ax.add_patch(Rectangle(((i - .5) * a, (j - .5) * a), a, a, fc="#dfe7f0" if (i, j) != (0, 0) else "#bcd0e8",
                                   ec="k", lw=1))
    ax.text(0, -1.45 * a, "grey-blue: neighbour modules (3x3 cluster)", ha="center", fontsize=8, color="#555")
    cx, cy = 0.31 * a, 0.22 * a
    ax.plot(0, 0, "ks", ms=8)
    ax.annotate("seed module centre\n(x_seed, y_seed)", (0, 0), (-1.45 * a, 0.2 * a), fontsize=9,
                arrowprops=dict(arrowstyle="->", color="k"))
    ax.plot(cx, cy, "o", color=COL["unc"], ms=9)
    ax.annotate("cluster centroid\n(x_cluster, y_cluster)\nC++ log-weighted, W0 = 5.2", (cx, cy), (0.9 * a, 1.15 * a),
                fontsize=9, color=COL["unc"], arrowprops=dict(arrowstyle="->", color=COL["unc"]))
    ax.add_patch(FancyArrowPatch((0, -0.06 * a), (cx, -0.06 * a), arrowstyle="<->", mutation_scale=12, color="#1f5fa8"))
    ax.add_patch(FancyArrowPatch((cx, 0), (cx, cy), arrowstyle="<->", mutation_scale=12, color="#2a9d8f"))
    ax.text(cx / 2, -0.20 * a, "u", color="#1f5fa8", ha="center", fontsize=13)
    ax.text(cx + 0.06 * a, cy / 2, "v", color="#2a9d8f", fontsize=13)
    ax.axhline(0, color="#999", lw=.5, ls=":"); ax.axvline(0, color="#999", lw=.5, ls=":")
    ax.add_patch(Rectangle((-.5 * a, -.5 * a), a, a, fc="none", ec="#e69f00", lw=2.2))
    ax.text(-.5 * a, -.5 * a - 0.2 * a, "orange box: u, v = +-0.5 (module edge)", color="#b37800", ha="left", fontsize=8.5)
    ax.set_xlim(-2.4 * a, 4.4 * a); ax.set_ylim(-1.6 * a, 1.7 * a); ax.set_aspect("equal")
    ax.set_xlabel("x (module widths, 1 = 20.77 mm)"); ax.set_ylabel("y")
    ax.set_xticks([-a, 0, a]); ax.set_xticklabels(["-1", "0", "+1"]); ax.set_yticks([-a, 0, a]); ax.set_yticklabels(["-1", "0", "+1"])
    ax.grid(False)
    txt = ("u = (x_cluster - x_seed) / size_x\n"
           "v = (y_cluster - y_seed) / size_y\n"
           "r = E_cluster / E'(seed)\n"
           "E'(seed) = elastic kinematic energy\n"
           "   at the seed-module centre\n"
           "correction: E_cluster -> E_cluster / R(u, v)")
    ax.text(1.75 * a, -1.5 * a, txt, fontsize=9, family="monospace", va="bottom",
             bbox=dict(boxstyle="round", fc="#fffbe6", ec="#999"))
    ax.set_title("Position-correction coordinates (u, v) inside the seed module")
    save(fig, out, "fig01_definition.png")


# ---- fig02 -----------------------------------------------------------------

def fig02(out, d, hyb, bil, slug):
    nb = 24
    ed = np.linspace(-CL, CL, nb + 1)
    cen = 0.5 * (ed[1:] + ed[:-1])
    Rh = hyb_R(hyb, d["u"], d["v"], d["seed"])
    rc = d["r"] / Rh
    res = {}
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.2), sharey=True)
    for ax, (xs, cut, lab, key) in zip(axs, ((d["u"], np.abs(d["v"]) < 0.1, "u", "u"), (d["v"], np.abs(d["u"]) < 0.1, "v", "v"))):
        for rr, name, c in ((d["r"], "uncorrected", COL["unc"]), (rc, "after hybrid correction", COL["cor"])):
            y, e, n = [], [], []
            for lo, hi in zip(ed[:-1], ed[1:]):
                s = cut & (xs >= lo) & (xs < hi)
                m, se, nn = tmean(rr[s]) if s.sum() > 30 else (np.nan, np.nan, int(s.sum()))
                y.append(m); e.append(se); n.append(nn)
            ax.errorbar(cen, y, e, fmt="o-", ms=3.5, lw=1, color=c, label=name, capsize=2)
            res[f"{key}_{name}"] = (list(map(float, y)), list(map(float, e)), list(map(int, n)))
        ax.axvline(-.5, color="#999", ls="--", lw=.8); ax.axvline(.5, color="#999", ls="--", lw=.8)
        ax.set_xlabel(f"{lab} = (cluster - seed) / module width" + ("   (|v| < 0.1)" if lab == "u" else "   (|u| < 0.1)"))
    axs[0].set_ylabel("r = E_cluster / E'   (+-2 sigma truncated mean)")
    axs[0].legend(loc="lower center", fontsize=8.5)
    fig.suptitle(f"The bowl: cluster energy vs impact position, {slug}, all seeds pooled "
                 f"({len(d['r'])} clusters, r in [0.85, 1.12])", fontsize=10.5)
    save(fig, out, "fig02_bowl_profile.png")
    res["n_pool"] = int(len(d["r"]))
    res["bins_u"] = list(map(float, cen))
    FACTS["fig02"] = res

    # 2D bilinear map
    Rm = np.array(bil["R"]); fl = np.array(bil["filled"]).astype(bool); n = np.array(bil["n"]); se = np.array(bil["se"])
    ue = np.array(bil["u_edges"]); ve = np.array(bil["v_edges"])
    dmax = max(abs(Rm - 1).max(), 1e-3)
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.8), constrained_layout=True)
    im = axs[0].pcolormesh(ue, ve, Rm.T, cmap="RdBu_r", norm=Normalize(1 - dmax, 1 + dmax))
    for i in range(12):
        for j in range(12):
            if fl[i, j]:
                axs[0].add_patch(Rectangle((ue[i], ve[j]), ue[1] - ue[0], ve[1] - ve[0], fc="none", ec="k", hatch="////", lw=0))
    axs[0].set_aspect("equal"); axs[0].set_xlabel("u"); axs[0].set_ylabel("v"); axs[0].grid(False)
    axs[0].set_title(f"R(u, v): {int((~fl).sum())} fitted, {int(fl.sum())} filled (hatched)", fontsize=9)
    plt.colorbar(im, ax=axs[0], label="R (scale symmetric about 1)")
    im2 = axs[1].pcolormesh(ue, ve, n.T, cmap="viridis")
    axs[1].set_aspect("equal"); axs[1].set_xlabel("u"); axs[1].grid(False)
    axs[1].set_title("events per bin (n)", fontsize=9)
    plt.colorbar(im2, ax=axs[1], label="n")
    save(fig, out, "fig02b_bilinear_map.png")
    FACTS["fig02b"] = dict(Rmin=float(Rm[~fl].min()), Rmax=float(Rm[~fl].max()), R_center=float(Rm[5:7, 5:7].mean()),
                           n_filled=int(fl.sum()), n_fitted=int((~fl).sum()), norm=float(bil["norm"]),
                           se_median=float(np.median(se[~fl])), n_min_fitted=int(n[~fl].min()), n_max=int(n.max()))


# ---- fig03 -----------------------------------------------------------------

def fig03(out, d, hyb, geo, mods):
    picks = [("W565", "own fit, highest statistics"), ("W489", "own fit, just above the 3000 threshold"),
             ("W601", "fallback (n < 3000): embedded bilinear map")]
    fig, axs = plt.subplots(3, 3, figsize=(12.5, 10))
    ff = {}
    for row, (mod, desc) in enumerate(picks):
        s = d["seed"] == mod
        u, v, r = d["u"][s], d["v"][s], d["r"][s]
        g = P.grid_estimates(u, v, r, P.HYB_GRID, CL, P.HYB_MIN_BIN)
        uc, vc, Rv, se, n = g
        sc = float((n * Rv).sum() / n.sum())
        own = mod in hyb["modules"]
        if own:
            c = np.array(hyb["modules"][mod]["c"])
            f = lambda a, b: P.poly_R_clamped(c, a, b, CL)
            name = "poly8"
        else:
            f = lambda a, b: bench.poscorr_eval(hyb["fb"], a, b)[0]
            name = "bilinear fallback"
        t = np.linspace(-CL, CL, 121)
        for k, (lab, sel) in enumerate((("u", (np.abs(vc - 0.1) < 1e-6, np.abs(vc + 0.1) < 1e-6)),
                                        ("v", (np.abs(uc - 0.1) < 1e-6, np.abs(uc + 0.1) < 1e-6)))):
            ax = axs[row, k]
            for sg, sl, cc in ((+0.1, sel[0], "#1f5fa8"), (-0.1, sel[1], "#c0392b")):
                x = (uc if lab == "u" else vc)[sl]
                ax.errorbar(x, (Rv / sc)[sl], (se / sc)[sl], fmt="o", color=cc, capsize=2,
                            label=f"6x6 grid, {'v' if lab == 'u' else 'u'} = {sg:+.1f}")
                ax.plot(t, f(t, np.full_like(t, sg)) if lab == "u" else f(np.full_like(t, sg), t), "-", color=cc, lw=1.2)
            ax.set_xlabel(lab); ax.set_ylabel("R" if k == 0 else "")
            ax.axvline(-.5, color="#999", ls="--", lw=.7); ax.axvline(.5, color="#999", ls="--", lw=.7)
            if k == 0:
                ax.set_title(f"{mod}: {desc}\nn = {len(r)}, cut along u ({name} curves)", fontsize=9)
            else:
                ax.set_title(f"cut along v ({name} curves)", fontsize=9)
            ax.legend(fontsize=7.5)
        ax = axs[row, 2]
        resid = (Rv / sc - f(uc, vc))
        sc_ = ax.scatter(uc, vc, c=100 * resid, cmap="RdBu_r", vmin=-3, vmax=3, s=170, marker="s", ec="k")
        ax.set_xlabel("u"); ax.set_ylabel("v"); ax.set_aspect("equal"); ax.grid(False)
        chi = float(np.sum((resid / (se / sc)) ** 2) / len(resid))
        ax.set_title(f"residual grid - {name} (%); {len(resid)} bins, chi2/bin = {chi:.2f}", fontsize=9)
        plt.colorbar(sc_, ax=ax, label="residual (%)")
        if own:
            txt = "c0..c7 =\n" + "\n".join(f"c{i} = {c[i]:+.4f}" for i in range(8))
            axs[row, 0].text(0.02, 0.03, txt, transform=axs[row, 0].transAxes, fontsize=7, family="monospace", va="bottom",
                             bbox=dict(fc="w", ec="#999", alpha=.9))
        ff[mod] = dict(n=int(len(r)), own=bool(own), n_grid_bins=int(len(Rv)), grid_scale=sc, chi2_per_bin=chi,
                       x=float(geo.x[mod]), y=float(geo.y[mod]), c=(list(map(float, c)) if own else None),
                       resid_rms_pct=float(100 * resid.std()))
    fig.suptitle("Per-module poly8 vs the 6x6 truncated-mean grid (points = bin means / n-weighted mean; lines = fitted surface)", fontsize=10.5)
    fig.tight_layout()
    save(fig, out, "fig03_module_poly8.png")
    FACTS["fig03"] = ff


# ---- fig04 -----------------------------------------------------------------

def fig04(out, geo, mods):
    t = mods.merge(geo[["x", "y", "sx", "sy"]], left_on="module", right_index=True)
    own = t[t.own == 1]; fb = t[t.own == 0]
    fig, axs = plt.subplots(2, 2, figsize=(11.5, 10))
    FACTS["fig04"] = {}
    for ax, (cn, ttl) in zip(axs.ravel(), (("c1", "c1 (u^2 curvature)"), ("c2", "c2 (v^2 curvature)"),
                                           ("c6", "c6 (linear in u)"), ("c7", "c7 (linear in v)"))):
        for _, m in geo.iterrows():
            ax.add_patch(Rectangle((m.x - m.sx / 2, m.y - m.sy / 2), m.sx, m.sy, fc="#f2f2f2", ec="none"))
        vmax = np.percentile(np.abs(own[cn]), 98) if cn in ("c6", "c7") else None
        kw = dict(cmap="RdBu_r", vmin=-vmax, vmax=vmax) if vmax else dict(cmap="viridis")
        sc = ax.scatter(own.x, own.y, c=own[cn], marker="s", s=62, **kw)
        ax.scatter(fb.x, fb.y, marker="x", s=14, color="#888", lw=.7)
        ax.set_aspect("equal"); ax.set_xlim(-130, 130); ax.set_ylim(-130, 130); ax.grid(False)
        ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)"); ax.set_title(ttl)
        plt.colorbar(sc, ax=ax, shrink=.8)
        FACTS["fig04"][cn] = dict(median=float(own[cn].median()), mean=float(own[cn].mean()), min=float(own[cn].min()),
                                  max=float(own[cn].max()), n_own=int(len(own)))
    fig.suptitle(f"Per-module poly8 coefficients over the HyCal layout: {len(own)} own fits (squares), "
                 f"{len(fb)} fallback modules (grey x), 3.49 GeV H2", fontsize=10.5)
    fig.tight_layout()
    save(fig, out, "fig04_coeff_layout.png")
    fig, axs = plt.subplots(2, 3, figsize=(11.5, 6))
    for ax, k in zip(axs.ravel(), ("c1", "c2", "c4", "c5", "c6", "c7")):
        ax.hist(own[k], bins=20, color=COL["cor"], alpha=.8)
        ax.axvline(own[k].median(), color="k", ls="--", lw=1)
        ax.set_xlabel(k); ax.set_title(f"{k}: median {own[k].median():+.3f}", fontsize=9)
        FACTS["fig04"].setdefault(k, dict(median=float(own[k].median()), mean=float(own[k].mean()),
                                          min=float(own[k].min()), max=float(own[k].max()), n_own=int(len(own))))
    fig.suptitle("Distribution of poly8 coefficients over the own-fit modules", fontsize=10.5)
    fig.tight_layout()
    save(fig, out, "fig04b_coeff_hist.png")


# ---- fig05 -----------------------------------------------------------------

def read_result(path):
    """(centres, counts_legacy_post, cb row dict) from a benchmark result.tsv"""
    lines = Path(path).read_text().splitlines()
    i = next(k for k, l in enumerate(lines) if l.startswith("e_centre\tlegacy_pre"))
    h = []
    for l in lines[i + 1:]:
        if not l.strip():
            break
        h.append(l.split("\t"))
    h = np.array(h, float)
    j = next(k for k, l in enumerate(lines) if l.startswith("path\tn_modules"))
    cols = lines[j].split("\t"); vals = lines[j + 1].split("\t")
    row = dict(zip(cols, vals))
    assert row["path"] == "legacy"
    return h[:, 0], h[:, 2], row


def fig05(out, d, hyb, chain):
    fin = chain / "scan/splits/scan_024386/split_00010/result.tsv"
    unc = chain / "poscorr/scan/scan_024386/split_00010/result.tsv"
    e0, c0, r0 = read_result(unc)
    e1, c1, r1 = read_result(fin)
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.4))
    ax = axs[0]
    for e, c, r, name, col in ((e0, c0, r0, "uncorrected", COL["unc"]), (e1, c1, r1, "hybrid poly8 applied", COL["cor"])):
        ax.step(e, c, where="mid", color=col, lw=1.3,
                label=f"{name}: CB sigma_E/E = {100 * float(r['sigma_over_E_cb']):.2f}% (mu {float(r['mu_cb_MeV']):.1f} MeV, n_inner {r['n_inner_events']})")
    ax.set_xlim(3000, 3900); ax.set_xlabel("E (MeV)  inner-ring spectrum after per-module recalibration"); ax.set_ylabel("events / bin")
    ax.legend(fontsize=7.5, loc="upper left"); ax.set_title("Benchmark elastic peak, run 24386 split 10 (legacy path)", fontsize=9.5)
    # cluster-level r
    s = d["splits"]
    ax = axs[1]
    rr = d["r"]; Rh = hyb_R(hyb, d["u"], d["v"], d["seed"])
    # events of split 24386/10 only: rebuild a small pool
    args = argparse.Namespace(splits="all", path="legacy", keep_moller=False, inner_radius=200.0, exclude_rowcol=(15, 20, 15, 20))
    geo = R.load_geometry(24386)
    sd = chain / "poscorr/scan"
    # load one split directly through load_pool's per-split logic by copying the scan dir view
    import tempfile, os
    td = Path(tempfile.mkdtemp(prefix="sup-pcx-", dir=os.environ.get("TMPDIR", "/private/tmp")))
    (td / "scan_024386").mkdir()
    os.symlink((sd / "scan_024386/split_00010").resolve(), td / "scan_024386/split_00010")
    p1 = P.load_pool(td, 24386, args, geo)
    import shutil; shutil.rmtree(td)
    m = (p1["r"] >= P.HYB_RMIN) & (p1["r"] <= P.HYB_RMAX)
    u, v, r, sd_ = p1["u"][m], p1["v"][m], p1["r"][m], p1["seed"][m]
    rc = r / hyb_R(hyb, u, v, sd_)
    bins = np.linspace(0.90, 1.08, 73)
    st = {}
    for x, name, col in ((r, "uncorrected r", COL["unc"]), (rc, "r / R(u, v)", COL["cor"])):
        mu, sig, nin, _ = P.trunc_mean(x)
        ax.hist(x, bins=bins, histtype="step", color=col, lw=1.4, label=f"{name}: trunc. mean {mu:.4f}, sigma {100 * sig:.2f}%")
        st[name] = (mu, sig, nin)
    ax.set_xlabel("E_cluster / E'(seed), seed in pool (r < 200 mm)"); ax.set_ylabel("clusters / bin")
    ax.legend(fontsize=8, loc="upper left"); ax.set_title("Cluster level, same split, all pooled seeds", fontsize=9.5)
    fig.tight_layout()
    save(fig, out, "fig05_spectrum.png")
    FACTS["fig05"] = dict(unc_file=str(unc), cor_file=str(fin),
                          unc={k: r0[k] for k in ("n_modules", "n_inner_events", "mu_cb_MeV", "sigma_cb_MeV", "sigma_over_E_cb")},
                          cor={k: r1[k] for k in ("n_modules", "n_inner_events", "mu_cb_MeV", "sigma_cb_MeV", "sigma_over_E_cb")},
                          cluster_level={k: [float(a) for a in v_[:2]] + [int(v_[2])] for k, v_ in st.items()},
                          n_cluster_split=int(len(r)))


# ---- fig06 -----------------------------------------------------------------

def fig06(out, form):
    t = pd.read_csv(form / "tables/3p49_heldout_per_module.tsv", sep="\t")
    forms = ["bilinear", "poly8", "even6", "even4"]
    cols = {"bilinear": COL["bil"], "poly8": COL["cor"], "even6": COL["even6"], "even4": COL["even4"]}
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw=dict(width_ratios=[1.5, 1]))
    ax = axs[0]
    chi = {f: t[f"chi2_{f}"] / t["n_heldout_bins"] for f in forms}   # chi2/ndf as in 3p49_heldout.md
    bins = np.logspace(np.log10(min(c.min() for c in chi.values()) * 0.9), np.log10(max(c.max() for c in chi.values()) * 1.1), 30)
    med = {}
    for f in forms:
        x = chi[f]
        ax.hist(x, bins=bins, histtype="step", lw=1.6, color=cols[f], label=f"{f}: median {x.median():.3f}")
        med[f] = (float(x.median()), float(x.mean()))
    ax.set_xscale("log"); ax.set_xlabel("held-out chi2/ndf per module (5-fold, 6x6 bins)"); ax.set_ylabel("modules")
    ax.set_title(f"Held-out chi2/ndf, {len(t)} modules with >= 3000 clusters", fontsize=9.5); ax.legend(fontsize=8.5)
    # sigma bars
    txt = (form / "tables/3p49_sigma.md").read_text()
    leg = txt.split("## legacy")[1].split("## gated")[0]
    sig = {}
    for l in leg.splitlines():
        m = re.match(r"\|\s*(poly8|bilinear|even6)\s*\|\s*([\d.]+) \+- ([\d.]+)\s*\|", l)
        if m:
            sig[m.group(1)] = (float(m.group(2)), float(m.group(3)))
    ax = axs[1]
    order = ["bilinear", "even6", "poly8"]
    ax.bar(order, [sig[k][0] for k in order], yerr=[sig[k][1] for k in order], color=[cols[k] for k in order], capsize=4)
    for i, k in enumerate(order):
        ax.text(i, sig[k][0] + 0.012, f"{sig[k][0]:.4f}\n+-{sig[k][1]:.4f}", ha="center", fontsize=8)
    ax.set_ylim(1.5, 1.78); ax.set_ylabel("sigma_E/E, CB (%)")
    ax.set_title("3.49 GeV H2 legacy, mean +- SE over 27 splits", fontsize=9.5); ax.grid(axis="x", visible=False)
    fig.tight_layout()
    save(fig, out, "fig06_form_comparison.png")
    FACTS["fig06"] = dict(chi2_median_mean=med, sigma_cb_legacy=sig, n_modules=int(len(t)))


# ---- main ------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-type", default="3p49GeV_H2")
    ap.add_argument("--out", required=True)
    ap.add_argument("--chain", default="output/r5ri_20261104")
    ap.add_argument("--form-dir", default="output/poscorr_form_20261105")
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    chain, form, slug = Path(a.chain), Path(a.form_dir), a.run_type
    mp = chain / f"poscorr/maps/{slug}"
    bil = json.load(open(mp / f"poscorr_r5ri_bilinear_{slug}_all.json"))
    hyb = bench.poscorr_load(str(mp / f"poscorr_r5ri_{slug}_all.json"))
    mods = pd.read_csv(mp / f"poscorr_r5ri_{slug}_all_modules.tsv", sep="\t")
    d, geo, runs, scan = load_all(chain, slug)
    print(slug, runs, len(d["r"]), "clusters")
    fig01(out)
    fig02(out, d, hyb, bil, slug)
    fig03(out, d, hyb, geo, mods)
    fig04(out, geo, mods)
    fig05(out, d, hyb, chain)
    fig06(out, form)
    FACTS["meta"] = dict(slug=slug, runs=runs, n_pool=int(len(d["r"])), n_splits=len(d["splits"]))
    json.dump(FACTS, open(out / "facts_data.json", "w"), indent=1)


if __name__ == "__main__":
    main()
