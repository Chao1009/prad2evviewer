#!/usr/bin/env python3
"""exotic_p6_amp_vs_w50.py -- P6 of the exotic-waveform study: amplitude vs w50 (the f2 variable), per run type, and the diagnosis of the high-amplitude f2 group.

    source setup_env.sh
    python3 waveform_analysis/pileup_analysis/exotic_events/exotic_p6_amp_vs_w50.py [--jobs 8] [--plot-only] [--cache <pkl>] [--out <dir>]

Flags: the shared offline definition of exotic_p4_module_rates.flags_of with the uniform P2 thresholds; every category here is the EXCLUSIVE one (priority f4 > f5 > f1 > f2 > f3; f2 excludes clipped and truncated
records).  Run 24937 is excluded (user: not a good run).
Population: PbWO4 kept channels (A_main >= KEEP_ADC; gate and LMS veto as in pileup_rates.scan).  Pages (2 x 2 run types):
  (a) fitted A_max (= A_main) vs w50, 2D log-density, ALL kept channels of the sample npz (exotic_sample/, 2 files per run; the only H2 files with amplitudes), full w50 range
      (x axis to 8, the rest in the last bin); w50 = 3 (the f2 cut) and A = 300 ADC lines; red dots: f2-flagged rows.  Three pages: all, k = 1, k >= 2.
  (b) raw peak (max sample - pedestal) vs w50 for the rows with w50 > 2 (the only rows the full-scan H2 pickles keep; the pickles carry neither amplitude nor event numbers,
      so the raw peak cannot be obtained for the H2 full scan).  Raw peak by an EVIO decode-only re-read (exotic_p5_max_amplitude.raw_peaks): H2 from the sample npz,
      Ta X17 from the 19 full-scan npz.  Saturated rows (nclip >= 1: f4 red triangles, f5 magenta squares) are drawn on top of the density of the unclipped rows; f2 rows are open black circles.
  Table (TSV/pages): f2 rows with A > 300 ADC (fraction per run type; k distribution and top modules of that group), w50-bin profile, saturated fraction of the
      kept channels and of those with E_main > 0.5 E_beam, with the elastic-edge amplitude (beam energy / module gain) for context.
  Examples: 20 waveforms of the high-amplitude f2 group (A > 300 ADC, unclipped) at 2.24 GeV H2 and 0.73 GeV H2 (re-fitted with the default fit from EVIO, as P3); and, for all four run types,
      20 waveforms drawn uniformly over the whole f2 population (w50 > 3, any A), with an overview scatter of the population marking the drawn rows.
Outputs (--out): exotic_P6_amp_vs_w50_20261118.pdf and .tsv.
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys
import time
import pickle
import argparse
import tempfile
import multiprocessing as mp
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import numpy as np                                      # noqa: E402
import exotic_p1_quantities as p1                       # noqa: E402
from exotic_p1_quantities import TYPES, RT_LABEL, write_tsv, EXCL_NOTE     # noqa: E402
from exotic_p3_sample import text_page                  # noqa: E402
from exotic_p4_module_rates import flags_of, UNI, wilson, ORDER_TXT     # noqa: E402
from exotic_p5_max_amplitude import raw_peaks, source_files  # noqa: E402

STEM = "exotic_P6_amp_vs_w50_20261118"
A_HI = 300.0                                            # ADC: the high-amplitude f2 group
WX = np.linspace(0.0, 8.0, 161)                         # w50 bins (overflow goes into the last bin)
KEEP = 10.0                                             # pileup_rates.KEEP_ADC
AY = np.logspace(np.log10(KEEP), np.log10(8000.0), 91)
RY = np.logspace(np.log10(3.0), np.log10(4200.0), 81)   # raw peak bins
RX = np.linspace(2.0, 8.0, 97)
WBINS = ((0.0, 1.5), (1.5, 2.0), (2.0, 3.0), (3.0, 4.0), (4.0, 99.0))
NEX, SPARE = 20, 8
SEED = 20261118


def cdict(z, sel):
    """Candidate-row dict for flags_of from a scan npz (kept rows only)."""
    return dict(k=z["ch_k"].astype(np.int64)[sel], db=z["ch_ex_db"].astype(np.float64)[sel], tail=z["ch_ex_tail"].astype(np.float64)[sel], w50=z["ch_w50"].astype(np.float64)[sel],
                c2s=z["ch_ex_chi2s"].astype(np.float64)[sel], c2k=z["ch_ex_chi2k"].astype(np.float64)[sel], nm=z["ch_ex_nm"].astype(np.int64)[sel],
                nclip=z["ch_nclip"].astype(np.int64)[sel], trunc=z["ch_trunc_main"].astype(bool)[sel], s4=z["ch_ex_f4"].astype(bool)[sel])


def work(task):
    """One npz -> partial: (a) 2D histograms and f2 rows (sample files), (b) rows with w50 > 2 with raw peak (H2 sample / Ta full files), f2 rows, saturation counts."""
    fn, kind, evio_dir = task
    import pileup_rates as pr
    t_start = time.time()
    z = np.load(fn)
    run, rt = int(z["run"]), str(z["run_type"])
    midx_all = z["ch_midx"].astype(np.int64)
    sel = z["mod_pwo4"][midx_all]
    nk = int(sel.sum())
    k = z["ch_k"].astype(np.int64)
    a = z["ch_a"].astype(np.float64)
    valid = np.arange(a.shape[1])[None, :] < k[:, None]
    amax = np.where(valid, a, -np.inf).max(axis=1)
    names = [str(n) for n in z["mod_name"]]
    cfg = pr.run_config(run)
    gain = pr.module_gain(cfg, rt, names)
    F = flags_of(cdict(z, sel), UNI)                    # kept rows only, in order
    F[2] = F["cat"] == 2                                # exclusive f2 (f4 > f5 > f1 > f2 > f3)
    i_sel = np.nonzero(sel)[0]
    w50 = z["ch_w50"].astype(np.float64)[sel]
    ks, am, mi = k[sel], amax[sel], midx_all[sel]
    nclip = z["ch_nclip"].astype(np.int64)[sel]
    emain = am * gain[mi] / cfg["beam"]
    part = dict(run=run, rt=rt, kind=kind, n_kept=nk, names=names, wall=0.0, edge=(cfg["beam"] / gain)[z["mod_pwo4"]])
    # saturation counts: all kept / E_main > 0.5 E_beam
    hiE = np.isfinite(emain) & (emain > 0.5)
    part["sat"] = dict(n=nk, clip=int((nclip >= 1).sum()), f4=int(F[4].sum()), f5=int(F[5].sum()), nE=int(hiE.sum()), clipE=int((hiE & (nclip >= 1)).sum()),
                       f4E=int((hiE & F[4]).sum()), f5E=int((hiE & F[5]).sum()))
    # f2 rows (all kinds): fitted amplitude, k, module, w50, row index in the npz
    f2 = F[2]
    part["f2"] = dict(amax=am[f2].astype(np.float32), k=ks[f2].astype(np.int8), midx=mi[f2].astype(np.int32), w50=w50[f2].astype(np.float32),
                      row=i_sel[f2].astype(np.int64), fid=np.full(int(f2.sum()), run * 100000 + int(z["splits"][0]), np.int64), evnum=z["ev_number"][z["ch_ev"].astype(np.int64)[i_sel[f2]]], chi2=z["ch_chi2"].astype(np.float32)[i_sel[f2]])
    if kind == "sample":
        x = np.clip(w50, 0, WX[-1] - 1e-9)
        ok = np.isfinite(x) & np.isfinite(am) & (am > 0)
        h = {}
        for nm, m in (("all", ok), ("k1", ok & (ks == 1)), ("k2", ok & (ks >= 2))):
            h[nm] = np.histogram2d(x[m], np.clip(am[m], AY[0], AY[-1] - 1e-6), bins=(WX, AY))[0].astype(np.int64)
        part["h"] = h
        part["n_w50_gt8"] = int((w50 > WX[-1]).sum())
        prof = []
        for lo, hi in WBINS:
            m = ok & (w50 > lo) & (w50 <= hi) if lo > 0 else ok & (w50 <= hi)
            v = am[m]
            prof.append(dict(lo=lo, hi=hi, n=int(m.sum()), med=float(np.median(v)) if v.size else np.nan, n_hi=int((v > A_HI).sum()), n_k2=int((ks[m] >= 2).sum()),
                             n_clip=int((nclip[m] >= 1).sum())))
        part["prof"] = prof
    # (b) rows with w50 > 2: H2 from the sample, Ta from the full scan
    if (kind == "sample") != ("Ta" in rt):
        c = np.nonzero(np.isfinite(w50) & (w50 > 2.0))[0]
        raw, smax = raw_peaks(z, run, midx_all, i_sel[c], evio_dir)
        part["b"] = dict(w50=w50[c].astype(np.float32), raw=raw.astype(np.float32), smax=smax.astype(np.float32), amax=am[c].astype(np.float32), k=ks[c].astype(np.int8),
                         nclip=nclip[c].astype(np.int16), f2=F[2][c], f4=F[4][c], f5=F[5][c], midx=mi[c].astype(np.int32))
    if True:                                             # file info for the examples (every sample / full-scan npz)
        part["run_info"] = dict(fid=run * 100000 + int(z["splits"][0]), splits=[int(s) for s in z["splits"]], b=z["mod_b_ref"], s=z["mod_sig_ref"], beam=float(cfg["beam"]), gain=gain)
    part["wall"] = time.time() - t_start
    print(f"[work] {os.path.basename(fn)} {rt} ({kind}): {nk} kept, {int(F[2].sum())} f2, {0 if 'b' not in part else part['b']['w50'].size} rows w50>2, {part['wall']:.0f} s", flush=True)
    return part


def gather(parts):
    D = {rt: dict(h={q: 0 for q in ("all", "k1", "k2")}, f2=[], f2full=[], b=[], sat=dict(n=0, clip=0, f4=0, f5=0, nE=0, clipE=0, f4E=0, f5E=0), sat_full=dict(n=0, clip=0, f4=0, f5=0, nE=0, clipE=0, f4E=0, f5E=0),
                  prof=None, nfiles=0, nfull=0, gt8=0, edge=[], names=None, runs=[]) for rt in TYPES}
    for p in parts:
        d = D[p["rt"]]
        d["names"] = p["names"]
        if p["kind"] == "sample":
            d["nfiles"] += 1
            d["runs"].append(p["run"])
            for q in d["h"]:
                d["h"][q] = d["h"][q] + p["h"][q]
            d["f2"].append(dict(p["f2"], run=np.full(p["f2"]["amax"].size, p["run"])))
            for key in d["sat"]:
                d["sat"][key] += p["sat"][key]
            d["gt8"] += p["n_w50_gt8"]
            if d["prof"] is None:
                d["prof"] = [dict(x) for x in p["prof"]]
            else:
                for x, y in zip(d["prof"], p["prof"]):
                    for key in ("n", "n_hi", "n_k2", "n_clip"):
                        x[key] += y[key]
            d["edge"].append(p["edge"])
        else:
            d["nfull"] += 1
            d["f2full"].append(dict(p["f2"], run=np.full(p["f2"]["amax"].size, p["run"])))
            for key in d["sat_full"]:
                d["sat_full"][key] += p["sat"][key]
        if "b" in p:
            d["b"].append(dict(p["b"], run=np.full(p["b"]["w50"].size, p["run"])))
    for rt in TYPES:
        d = D[rt]
        for key in ("f2", "f2full", "b"):
            lst = d[key]
            d[key] = {q: np.concatenate([x[q] for x in lst]) for q in lst[0]} if lst else {}
        d["edge"] = np.nanmedian(np.stack(d["edge"]), axis=0) if d["edge"] else None
    return D


# ---------------------------------------------------------------------------------------------------------------- tables
def f2_table(D):
    rows = []
    for rt in TYPES:
        for src, f in (("sample (2 files/run)", D[rt]["f2"]), ("full scan npz (19 files)", D[rt]["f2full"])):
            if not f or ("H2" in rt and src.startswith("full")):
                continue
            n = f["amax"].size
            hi = f["amax"] > A_HI
            nh = int(hi.sum())
            lo_, hi_ = wilson(nh, n)
            r = dict(run_type=rt, source=src, n_f2=n, n_f2_A_gt300=nh, frac_A_gt300=nh / max(n, 1), lo68=float(lo_), hi68=float(hi_), A_median=float(np.median(f["amax"])) if n else np.nan, frac_lo68=float(lo_), frac_hi68=float(hi_))
            kk = f["k"][hi]
            kl = f["k"][~hi]
            for lab, v in (("hi", kk), ("lo", kl)):
                r[f"k1_{lab}"], r[f"k2_{lab}"], r[f"k3_{lab}"], r[f"k4p_{lab}"] = (int((v == 1).sum()), int((v == 2).sum()), int((v == 3).sum()), int((v >= 4).sum()))
            names = D[rt]["names"]
            mh = f["midx"][hi]
            u, c = np.unique(mh, return_counts=True)
            o = np.argsort(-c)
            for lab, v in (("hi", kk), ("lo", kl)):               # k >= 4 fraction of each group, Wilson 68 %
                c4, nn = int((v >= 4).sum()), int(v.size)
                l4, h4 = wilson(c4, nn) if nn else (float("nan"), float("nan"))
                r[f"frac_k4p_{lab}"], r[f"frac_k4p_{lab}_lo68"], r[f"frac_k4p_{lab}_hi68"] = (c4 / nn if nn else float("nan")), float(l4), float(h4)
            r["n_modules_hi"] = int(u.size)
            r["top_modules_hi"] = "; ".join(f"{names[u[j]]}:{int(c[j])}" for j in o[:6])
            r["top1_share_hi"] = float(c[o[0]] / max(nh, 1)) if nh else float("nan")
            t1l, t1h = wilson(int(c[o[0]]), nh) if nh else (float("nan"), float("nan"))
            r["top1_share_hi_lo68"], r["top1_share_hi_hi68"] = float(t1l), float(t1h)
            ml = f["midx"][~hi]
            u2, c2 = np.unique(ml, return_counts=True)
            o2 = np.argsort(-c2)
            r["top_modules_lo"] = "; ".join(f"{names[u2[j]]}:{int(c2[j])}" for j in o2[:6])
            r["w50_median_hi"] = float(np.median(f["w50"][hi])) if nh else float("nan")
            r["w50_median_lo"] = float(np.median(f["w50"][~hi])) if (~hi).any() else float("nan")
            rows.append(r)
    return rows


def sat_table(D):
    rows = []
    for rt in TYPES:
        for src, s, nf in (("sample", D[rt]["sat"], D[rt]["nfiles"]), ("full scan npz (Ta)", D[rt]["sat_full"], D[rt]["nfull"])):
            if s["n"] == 0:
                continue
            lo, hi = wilson(s["clip"], s["n"])
            loE, hiE = wilson(s["clipE"], s["nE"])
            ea = D[rt]["edge"]
            rows.append(dict(run_type=rt, source=src, n_kept=s["n"], n_clip=s["clip"], pct_clip=100 * s["clip"] / s["n"], lo68=100 * float(lo), hi68=100 * float(hi), n_f4=s["f4"], n_f5=s["f5"],
                             n_E_gt_0p5=s["nE"], n_clip_E_gt_0p5=s["clipE"], pct_clip_E_gt_0p5=100 * s["clipE"] / max(s["nE"], 1), lo68_E=100 * float(loE), hi68_E=100 * float(hiE),
                             n_f4_E=s["f4E"], n_f5_E=s["f5E"], edge_adc_median=float(np.nanmedian(ea)) if ea is not None else float("nan")))
    return rows


def prof_rows(D):
    rows = []
    for rt in TYPES:
        for x in D[rt]["prof"]:
            la, ha = wilson(x["n_hi"], x["n"])
            lk, hk = wilson(x["n_k2"], x["n"])
            rows.append(dict(run_type=rt, w50_lo=x["lo"], w50_hi=x["hi"], n_kept=x["n"], n_A_gt300=x["n_hi"], frac_A_gt300=x["n_hi"] / max(x["n"], 1), n_k_ge2=x["n_k2"],
                             frac_k_ge2=x["n_k2"] / max(x["n"], 1), n_clipped=x["n_clip"], frac_lo68=float(la), frac_hi68=float(ha), frac_k_ge2_lo68=float(lk), frac_k_ge2_hi68=float(hk)))
    return rows


def b_rows(D):
    rows = []
    for rt in TYPES:
        b = D[rt]["b"]
        if not b:
            continue
        ok = np.isfinite(b["raw"])
        sat = (b["nclip"] >= 1)
        hi = ok & ~sat & (b["w50"] > 3)
        rows.append(dict(run_type=rt, n_rows_w50_gt2=b["w50"].size, n_no_raw=int((~ok).sum()), n_saturated=int(sat.sum()), n_f4=int(b["f4"].sum()), n_f5=int(b["f5"].sum()), n_f2=int(b["f2"].sum()),
                         n_unclipped_w50_gt3=int(hi.sum()), n_unclipped_w50_gt3_raw_gt300=int((hi & (b["raw"] > A_HI)).sum()),
                         n_f2_raw_gt300=int((b["f2"] & ok & (b["raw"] > A_HI)).sum()), n_f2_fitA_gt300=int((b["f2"] & (b["amax"] > A_HI)).sum()),
                         source="sample npz" if "H2" in rt else "full scan npz (19 files)"))
    return rows


# ---------------------------------------------------------------------------------------------------------------- plots
def heat_page(pdf, plt, D, key, title):
    from matplotlib.colors import LogNorm
    fig, axs = plt.subplots(2, 2, figsize=(16, 11.5))
    for ax, rt in zip(axs.flat, TYPES):
        h = D[rt]["h"][key].T.astype(float)
        h[h == 0] = np.nan
        pc = ax.pcolormesh(WX, AY, h, norm=LogNorm(vmin=1, vmax=max(np.nanmax(h) if np.isfinite(h).any() else 1, 2)), cmap="viridis", rasterized=True)
        fig.colorbar(pc, ax=ax, fraction=0.045, label="kept channels / bin")
        f = D[rt]["f2"]
        if f and f["amax"].size:
            m = np.ones(f["amax"].size, bool) if key == "all" else ((f["k"] == 1) if key == "k1" else (f["k"] >= 2))
            idx = np.nonzero(m)[0]
            if idx.size > 3000:
                idx = np.random.default_rng(SEED).choice(idx, 3000, replace=False)
            ax.plot(np.clip(f["w50"][idx], 0, WX[-1] - 0.02), f["amax"][idx], "o", ms=3.2, mfc="none", mec="red", mew=0.6, alpha=0.8, label="f2 rows")
        ax.axvline(UNI["w50"], color="w", lw=1.2, ls="--")
        ax.axvline(UNI["w50"], color="k", lw=0.5, ls="--")
        ax.axhline(A_HI, color="orange", lw=1.0, ls=":")
        ax.set_yscale("log")
        ax.set_xlim(0, WX[-1])
        n = int(D[rt]["h"][key].sum())
        ax.set_title(f"{RT_LABEL[rt]}: N = {n}, {D[rt]['nfiles']} sample files" + (f", {D[rt]['gt8']} with w50 > 8 in the last bin" if key == "all" else ""), fontsize=10.5)
        ax.set_xlabel("w50 / w50(template)")
        ax.set_ylabel("fitted A_max [ADC]")
        ax.legend(loc="upper right", fontsize=8)
    fig.suptitle(title + "   (dashed: f2 cut w50 = 3; dotted: A = 300 ADC; red circles: f2-flagged rows, <= 3000 drawn)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    pdf.savefig(fig)
    plt.close(fig)


def raw_page(pdf, plt, D):
    from matplotlib.colors import LogNorm
    fig, axs = plt.subplots(2, 2, figsize=(16, 11.5))
    for ax, rt in zip(axs.flat, TYPES):
        b = D[rt]["b"]
        if not b:
            ax.text(0.5, 0.5, "no rows", transform=ax.transAxes)
            continue
        ok = np.isfinite(b["raw"])
        sat = b["nclip"] >= 1
        u = ok & ~sat
        h = np.histogram2d(np.clip(b["w50"][u], 2.0, RX[-1] - 1e-6), np.clip(b["raw"][u], RY[0], RY[-1] - 1e-6), bins=(RX, RY))[0].T
        h[h == 0] = np.nan
        if np.isfinite(h).any():
            pc = ax.pcolormesh(RX, RY, h, norm=LogNorm(vmin=1, vmax=max(np.nanmax(h), 2)), cmap="viridis", rasterized=True)
            fig.colorbar(pc, ax=ax, fraction=0.045, label="rows / bin (unclipped)")
        f2 = b["f2"] & ok
        idx = np.nonzero(f2)[0]
        if idx.size > 3000:
            idx = np.random.default_rng(SEED).choice(idx, 3000, replace=False)
        ax.plot(np.clip(b["w50"][idx], 2, RX[-1] - 0.02), b["raw"][idx], "o", ms=3.2, mfc="none", mec="k", mew=0.5, alpha=0.7, label=f"f2 rows ({int(f2.sum())})")
        for flg, mk, col, lab in (("f4", "^", "red", "f4 saturated"), ("f5", "s", "magenta", "f5 saturated single")):
            m = b[flg] & ok
            ax.plot(np.clip(b["w50"][m], 2, RX[-1] - 0.02), np.clip(b["raw"][m], RY[0], RY[-1]), mk, ms=5, mfc=col, mec="k", mew=0.4, alpha=0.9, label=f"{lab} ({int(m.sum())})")
        ax.axvline(UNI["w50"], color="k", lw=0.8, ls="--")
        ax.axhline(A_HI, color="orange", lw=1.0, ls=":")
        ax.set_yscale("log")
        ax.set_xlim(2, RX[-1])
        ax.set_title(f"{RT_LABEL[rt]}: {b['w50'].size} rows with w50 > 2 ({'sample npz' if 'H2' in rt else 'full scan, 19 npz'}); {int((~ok).sum())} without raw peak", fontsize=10.5)
        ax.set_xlabel("w50 / w50(template)")
        ax.set_ylabel("raw peak = max sample - pedestal [ADC]")
        ax.legend(loc="lower right", fontsize=8)
    fig.suptitle("(b) raw peak vs w50, rows with w50 > 2 (density: unclipped rows; saturated f4/f5 rows drawn on top; open circles: f2 rows)   dashed: w50 = 3; dotted: 300 ADC", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    pdf.savefig(fig)
    plt.close(fig)


def ci(lo, hi, nd):
    return f"[{lo:.{nd}f}, {hi:.{nd}f}]"


def table_pages(pdf, plt, D, t2, ts, tp, tb):
    L = [f"f2 rows (uniform P2 thresholds; w50 > 3, not truncated, not clipped) by fitted A_max; high group = A > {A_HI:.0f} ADC.  [68 % Wilson interval]", "",
         f"{'run type':<17}{'source':<25}{'n f2':>8}{'n A>300':>9}{'fraction':>10}{'[lo, hi]':>20}{'A med':>8}   k in hi group (1/2/3/4+)  k in lo group (1/2/3/4+)  hi: n mod, top-1 share [lo, hi], w50 med; k>=4 fraction hi / lo [lo, hi]"]
    for r in t2:
        L.append(f"{RT_LABEL[r['run_type']]:<17}{r['source']:<25}{r['n_f2']:>8d}{r['n_f2_A_gt300']:>9d}{r['frac_A_gt300']:>10.3f}{ci(r['lo68'], r['hi68'], 4):>20}{r['A_median']:>8.0f}"
                 f"   {r['k1_hi']:>5d}/{r['k2_hi']}/{r['k3_hi']}/{r['k4p_hi']}      {r['k1_lo']}/{r['k2_lo']}/{r['k3_lo']}/{r['k4p_lo']}      {r['n_modules_hi']}, {r['top1_share_hi']:.2f} {ci(r['top1_share_hi_lo68'], r['top1_share_hi_hi68'], 2)}, {r['w50_median_hi']:.2f}")
        L.append(f"{'':<42}k>=4 fraction: hi {r['frac_k4p_hi']:.3f} {ci(r['frac_k4p_hi_lo68'], r['frac_k4p_hi_hi68'], 3)} ({r['k4p_hi']} of {r['n_f2_A_gt300']}); lo {r['frac_k4p_lo']:.3f} {ci(r['frac_k4p_lo_lo68'], r['frac_k4p_lo_hi68'], 3)} ({r['k4p_lo']} of {r['n_f2'] - r['n_f2_A_gt300']})")
        L.append(f"{'':<42}top modules (hi): {r['top_modules_hi']}")
        L.append(f"{'':<42}top modules (lo): {r['top_modules_lo']}")
    L = [f"*** {EXCL_NOTE} ***  categories exclusive, priority {ORDER_TXT}", ""] + L
    text_page(pdf, plt, "P6 table 1: the high-amplitude f2 group", L, size=9)
    L = ["Saturated pulses (clipped record nclip >= 1 = f4 + f5) of the kept PbWO4 channels; E_main > 0.5 E_beam subset; [68 % interval].  Full-scan values (103 files) are in P4.", "",
         f"{'run type':<17}{'source':<20}{'n kept':>11}{'n clip':>8}{'%':>9}{'[lo, hi]':>20}{'f4':>6}{'f5':>6}{'n E>0.5':>10}{'clip':>6}{'%':>9}{'[lo, hi]':>18}{'edge ADC':>10}"]
    for r in ts:
        L.append(f"{RT_LABEL[r['run_type']]:<17}{r['source']:<20}{r['n_kept']:>11d}{r['n_clip']:>8d}{r['pct_clip']:>9.4f}{ci(r['lo68'], r['hi68'], 4):>20}{r['n_f4']:>6d}{r['n_f5']:>6d}"
                 f"{r['n_E_gt_0p5']:>10d}{r['n_clip_E_gt_0p5']:>6d}{r['pct_clip_E_gt_0p5']:>9.3f}{ci(r['lo68_E'], r['hi68_E'], 3):>18}{r['edge_adc_median']:>10.0f}")
    L += ["", "edge ADC = beam energy / module gain (median over PbWO4 modules), the elastic-edge amplitude; saturation sets in at ~4095 - pedestal.  Percentages in [lo, hi] are in %.", "",
          "w50 profile of all kept channels of the sample (n, % with A_max > 300 ADC [68 % Wilson], % with k >= 2 [68 % Wilson] per w50 bin):",
          f"{'run type':<17}" + "".join(f"{f'{lo:g}-{hi:g}':>50}" for lo, hi in WBINS)]
    for rt in TYPES:
        rows = [x for x in tp if x["run_type"] == rt]
        L.append(f"{RT_LABEL[rt]:<17}" + "".join(f"{str(x['n_kept']) + ', ' + format(100 * x['frac_A_gt300'], '.2f') + ' [' + format(100 * x['frac_lo68'], '.2f') + ', ' + format(100 * x['frac_hi68'], '.2f') + '], ' + format(100 * x['frac_k_ge2'], '.2f') + ' [' + format(100 * x['frac_k_ge2_lo68'], '.2f') + ', ' + format(100 * x['frac_k_ge2_hi68'], '.2f') + ']':>50}" for x in rows))
    L += ["", "(b) rows with w50 > 2 and their raw peaks:",
          f"{'run type':<17}{'rows':>8}{'no raw':>8}{'clipped':>9}{'f4':>6}{'f5':>6}{'f2':>8}{'f2 raw>300':>12}{'f2 fitA>300':>13}{'unclip w50>3':>14}{'raw>300':>9}"]
    for r in tb:
        L.append(f"{RT_LABEL[r['run_type']]:<17}{r['n_rows_w50_gt2']:>8d}{r['n_no_raw']:>8d}{r['n_saturated']:>9d}{r['n_f4']:>6d}{r['n_f5']:>6d}{r['n_f2']:>8d}{r['n_f2_raw_gt300']:>12d}{r['n_f2_fitA_gt300']:>13d}"
                 f"{r['n_unclipped_w50_gt3']:>14d}{r['n_unclipped_w50_gt3_raw_gt300']:>9d}")
    text_page(pdf, plt, "P6 table 2: saturation, w50 profile, raw-peak rows", L, size=9)


# ---------------------------------------------------------------------------------------------------------------- examples
def f2pop(D, rt):
    """The f2 (exclusive) rows of a run type: H2 from the sample npz, Ta X17 from the full-scan npz."""
    return D[rt]["f2full"] if "Ta" in rt else D[rt]["f2"]


def pick_examples(D, rt):
    """High-amplitude group (A > 300 ADC): NEX + SPARE random rows."""
    f = D[rt]["f2"]
    hi = np.nonzero(f["amax"] > A_HI)[0]
    rng = np.random.default_rng(SEED + len(rt))
    return rng.permutation(hi)[:NEX + SPARE], hi.size


def pick_uniform(D, rt):
    """NEX + SPARE rows drawn uniformly over the whole f2 population of the run type (any A)."""
    f = f2pop(D, rt)
    rng = np.random.default_rng(SEED + 7 + len(rt))
    return rng.permutation(f["amax"].size)[:NEX + SPARE], f["amax"].size


def example_records(D, evio_dir, jobs, picks):
    """picks = {(label, rt): idx into f2pop / D[rt]['f2']}; returns {(run, ev, midx): record} after one re-fit pass per file."""
    import exotic_p3_sample as p3
    need = {}                                            # fid -> {evnum: {midx: stored}}
    for (lab, rt), idx in picks.items():
        f = D[rt]["f2"] if lab == "hi" else f2pop(D, rt)
        for j in idx:
            need.setdefault(int(f["fid"][j]), {}).setdefault(int(f["evnum"][j]), {})[int(f["midx"][j])] = dict(k=int(f["k"][j]), chi2=float(f["chi2"][j]))
    tasks = []
    for fid, pk in need.items():
        ri = next(p for p in D["_runinfo"] if p["fid"] == fid)
        run = fid // 100000
        files = [(str(Path(evio_dir) / f"prad_{run:06d}.evio.{s:05d}"), s, ri["b"][i], ri["s"][i]) for i, s in enumerate(ri["splits"])]
        tasks.append(("rand", run, files, pk))
    with mp.get_context("fork").Pool(min(jobs, max(len(tasks), 1))) as pool:
        res = pool.map(p3.work, tasks, chunksize=1)
    return {(r["run"], r["ev"], r["midx"]): r for rr in res for r in rr}


def collect(D, rt, f, idx, recs):
    ex = []
    for j in idx:
        e = recs.get((int(f["run"][j]), int(f["evnum"][j]), int(f["midx"][j])))
        if e is not None and e.get("ok1") and "Y" in e:
            ex.append(dict(e, w50s=float(f["w50"][j]), nm=D[rt]["names"][int(f["midx"][j])], pos=int(j)))
        if len(ex) == NEX:
            break
    return ex


def draw_pages(pdf, plt, rt, ex, head, sub):
    import exotic_p2_flags as p2
    per = 12
    for p0 in range(0, max(len(ex), 1), per):
        chunk = ex[p0:p0 + per]
        fig = plt.figure(figsize=(24, 17))
        outer = fig.add_gridspec(3, 4, hspace=0.5, wspace=0.2, left=0.03, right=0.995, top=0.9, bottom=0.03)
        fig.text(0.03, 0.975, f"{RT_LABEL[rt]}: {head}, examples {p0 + 1}-{p0 + len(chunk)} of {len(ex)}  ({sub})   [{EXCL_NOTE}]", fontsize=13, fontweight="bold")
        fig.text(0.03, 0.95, "black dots: raw samples; red: total default fit; dashed: fitted pulses; grey: fit baseline; dotted green: stretched single pulse (where run); lower panel: data - fit.  Title: (A_max, w50) is the position in the overview scatter.", fontsize=11)
        for j, e in enumerate(chunk):
            gs = outer[j // 4, j % 4].subgridspec(2, 1, height_ratios=[3, 1], hspace=0.12)
            ax, axr = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
            ttl = (f"#{p0 + j + 1} run {e['run']} ev {e['ev']} {e['nm']}  k={e['k']}  w50={e['w50s']:.2f}  chi2/dof={e['chi2']:.2f}\nA_max={e['amain']:.0f} ADC  E={e['amain'] * e['gain']:.0f} MeV ({100 * e['amain'] * e['gain'] / e['beam']:.0f} % Eb)  "
                   f"nclip={e['nclip']}  f2={'Y' if p2.fl(e, 2) else 'N'}  (A, w50) = ({e['amain']:.0f}, {e['w50s']:.2f})\npulses A@t: " + ", ".join(f"{a:.0f}@{t:.0f}" for a, t in zip(e["a"], e["t"])))
            p2.draw(ax, axr, e["Y"], e["base"], e["model"], e["stretch"], ttl, "0.0")
            ax.title.set_fontsize(8)
        pdf.savefig(fig)
        plt.close(fig)


def overview_page(pdf, plt, D, uex):
    """f2 population per run type in (w50, A_max) with the uniformly drawn examples numbered."""
    fig, axs = plt.subplots(2, 2, figsize=(16, 11.5))
    for ax, rt in zip(axs.flat, TYPES):
        f = f2pop(D, rt)
        n = f["amax"].size
        idx = np.arange(n) if n <= 20000 else np.random.default_rng(SEED).choice(n, 20000, replace=False)
        ax.plot(np.clip(f["w50"][idx], 0, 8), f["amax"][idx], ".", ms=2.5, color="0.55", alpha=0.5, label=f"f2 rows ({n}{'' if n <= 20000 else ', 20000 drawn'})")
        for i, e in enumerate(uex[rt]):
            ax.plot(min(e["w50s"], 8), e["amain"], "o", ms=8, mfc="none", mec="red", mew=1.2)
            ax.annotate(str(i + 1), (min(e["w50s"], 8), e["amain"]), fontsize=8, color="red", xytext=(3, 3), textcoords="offset points")
        ax.axhline(A_HI, color="orange", lw=1.0, ls=":")
        ax.axvline(UNI["w50"], color="k", lw=0.8, ls="--")
        ax.set_yscale("log")
        ax.set_xlim(2.5, 8)
        ax.set_xlabel("w50 / w50(template)  (> 8 in the last bin)")
        ax.set_ylabel("fitted A_max [ADC]")
        ax.set_title(f"{RT_LABEL[rt]}: exclusive f2 population ({n}), {len(uex[rt])} uniform examples (red, numbered as on the example pages)", fontsize=10)
        ax.legend(loc="upper right", fontsize=8)
    fig.suptitle(f"w50 > 3 examples: where the drawn rows sit in the f2 population   ({EXCL_NOTE})", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    pdf.savefig(fig)
    plt.close(fig)


def example_rows_tsv(hex_, uex):
    rows = []
    for sect, dd in (("examples", hex_), ("examples_uniform", uex)):
        for rt in TYPES:
            for n, e in enumerate(dd.get(rt, []), 1):
                rows.append(dict(section=sect, run_type=rt, n=n, run=e["run"], ev=e["ev"], module=e["nm"], k=e["k"], w50=e["w50s"], chi2=e["chi2"], amain=e["amain"],
                                 E_MeV=e["amain"] * e["gain"], E_over_Eb=e["amain"] * e["gain"] / e["beam"], nclip=e["nclip"], pulses="; ".join(f"{a:.0f}@{t:.0f}" for a, t in zip(e["a"], e["t"]))))
    return rows


def report(data, out, evio_dir, jobs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    D = data["D"]
    t2, ts, tp, tb = f2_table(D), sat_table(D), prof_rows(D), b_rows(D)
    picks = {("hi", rt): pick_examples(D, rt)[0] for rt in ("2p24GeV_H2", "0p73GeV_H2")}
    nhi = {rt: pick_examples(D, rt)[1] for rt in ("2p24GeV_H2", "0p73GeV_H2")}
    upicks = {("uni", rt): pick_uniform(D, rt)[0] for rt in TYPES}
    npop = {rt: pick_uniform(D, rt)[1] for rt in TYPES}
    recs = example_records(D, evio_dir, jobs, {**picks, **upicks})
    hex_ = {rt: collect(D, rt, D[rt]["f2"], picks[("hi", rt)], recs) for rt in ("2p24GeV_H2", "0p73GeV_H2")}
    uex = {rt: collect(D, rt, f2pop(D, rt), upicks[("uni", rt)], recs) for rt in TYPES}
    with PdfPages(str(Path(out) / (STEM + ".pdf"))) as pdf:
        heat_page(pdf, plt, D, "all", "(a) fitted A_max vs w50, all kept PbWO4 channels (sample npz)")
        heat_page(pdf, plt, D, "k1", "(a) fitted A_max vs w50, k = 1 channels")
        heat_page(pdf, plt, D, "k2", "(a) fitted A_max vs w50, k >= 2 channels")
        raw_page(pdf, plt, D)
        table_pages(pdf, plt, D, t2, ts, tp, tb)
        for rt in ("2p24GeV_H2", "0p73GeV_H2"):
            draw_pages(pdf, plt, rt, hex_[rt], f"high-amplitude f2 group (A > {A_HI:.0f} ADC, unclipped)", f"random of {nhi[rt]} in the sample")
        overview_page(pdf, plt, D, uex)
        for rt in TYPES:
            draw_pages(pdf, plt, rt, uex[rt], "f2 (w50 > 3) examples drawn uniformly over the f2 population, any A", f"random of {npop[rt]} exclusive f2 rows ({'sample npz' if 'H2' in rt else 'full-scan npz'})")
    ex = example_rows_tsv(hex_, uex)
    stem = str(Path(out) / STEM)
    write_tsv(stem + ".tsv", [dict(section="f2_high_group", **r) for r in t2] + [dict(section="saturation", **r) for r in ts] + [dict(section="w50_profile", **r) for r in tp]
              + [dict(section="raw_rows", **r) for r in tb] + ex,
              ["section", "run_type", "source", "n_f2", "n_f2_A_gt300", "frac_A_gt300", "lo68", "hi68", "A_median", "k1_hi", "k2_hi", "k3_hi", "k4p_hi", "k1_lo", "k2_lo", "k3_lo", "k4p_lo", "n_modules_hi",
               "top_modules_hi", "top1_share_hi", "top_modules_lo", "w50_median_hi", "w50_median_lo", "n_kept", "n_clip", "pct_clip", "n_f4", "n_f5", "n_E_gt_0p5", "n_clip_E_gt_0p5", "pct_clip_E_gt_0p5",
               "lo68_E", "hi68_E", "n_f4_E", "n_f5_E", "edge_adc_median", "w50_lo", "w50_hi", "n_A_gt300", "n_k_ge2", "frac_k_ge2", "n_clipped", "n_rows_w50_gt2", "n_no_raw", "n_saturated", "n_f2_raw_gt300",
               "n_f2_fitA_gt300", "n_unclipped_w50_gt3", "n_unclipped_w50_gt3_raw_gt300", "run", "ev", "module", "k", "w50", "chi2", "amain", "E_MeV", "E_over_Eb", "nclip", "pulses", "n",
               "frac_lo68", "frac_hi68", "frac_k_ge2_lo68", "frac_k_ge2_hi68", "frac_k4p_hi", "frac_k4p_hi_lo68", "frac_k4p_hi_hi68", "frac_k4p_lo", "frac_k4p_lo_lo68", "frac_k4p_lo_hi68",
               "top1_share_hi_lo68", "top1_share_hi_hi68"])
    print(f"[report] {stem}.pdf; examples (high-A) " + ", ".join(f"{rt} {len(hex_[rt])}" for rt in hex_) + "; uniform " + ", ".join(f"{rt} {len(uex[rt])}" for rt in uex))
    for r in t2:
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})
    for r in ts:
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})
    for r in tp:
        print(r)
    for r in tb:
        print(r)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-dir", default=str(p1.PYDIR / "output" / "pileup_rates" / "exotic_sample"))
    ap.add_argument("--full-dir", default=str(p1.PYDIR / "output" / "pileup_rates" / "exotic_full"))
    ap.add_argument("--evio-dir", default=p1.EVIO_DIR)
    ap.add_argument("--out", default=str(HERE))
    ap.add_argument("--cache", default=str(Path(tempfile.gettempdir()) / "exotic_p6_cache.pkl"))
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--plot-only", action="store_true")
    a = ap.parse_args()
    if a.plot_only:
        data = pickle.load(open(a.cache, "rb"))
    else:
        t0 = time.time()
        tasks = [(f, kind, a.evio_dir) for f, kind in source_files(a.sample_dir, a.full_dir)]
        tasks = [t for t in tasks if t[1] == "sample" or "Ta" in str(np.load(t[0])["run_type"])]       # the full-scan npz are Ta X17 only; keep the guard explicit
        with mp.get_context("fork").Pool(a.jobs) as pool:
            parts = pool.map(work, tasks, chunksize=1)
        D = gather(parts)
        D["_runinfo"] = [dict(run=p["run"], **p["run_info"]) for p in parts if "run_info" in p]
        data = dict(D=D, elapsed=time.time() - t0)
        pickle.dump(data, open(a.cache, "wb"))
        print(f"[compute] {data['elapsed']:.0f} s")
    report(data, a.out, a.evio_dir, a.jobs)


if __name__ == "__main__":
    main()
