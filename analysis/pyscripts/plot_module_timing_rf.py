#!/usr/bin/env python3
"""plot_module_timing_rf.py - per-module HyCal timing referenced to the accelerator RF.

PURPOSE
    For every PbWO4 module, the RF phase  phi = t_module - t_RF  folded modulo the RF period
    (4.008 ns), from module_timing_dump.py (module pulses) joined to rf_time_dump.py (RF TDC
    hits) by (run, split, event number).  Four panels per page (0.73 GeV H2, 2.24 GeV H2,
    3.49 GeV H2, X17 Ta), plus the phase-vs-pulse-integral dependence (time walk) and summary
    pages.  The per-board RF phase of the readout is stable inside a run but jumps between runs
    (up to ~P/2), so each run's module peak is aligned onto the reference run's peak before the
    runs of a type are summed (see "per-run alignment").

DEFINITIONS
  Event selection   (trig & --trigger-mask) != 0 (default 0xf00), (trig & plot_module_timing.TRIG_EXCL) == 0
                    (no LMS flasher / alpha pulser), event found in the RF file (join by ev_num within a
                    split) and with >= 1 channel-A RF hit.  Rows: every dominant PbWO4 pulse of a selected
                    event, with m_E >= --min-e if --min-e > 0 (output names then get _e<int>, e.g.
                    module_timing_rf_e100.pdf).  With --runs the names get a _test suffix.
  RF period P       per split: mean of the spacings of consecutive channel-A hits (sorted in time, same
                    event) restricted to 120-136 ns, divided by 32 (the recorded edges are the RF divided
                    by 32; expected 4.00771 ns).  The split's P is used for that split's folding; the
                    entries-weighted per-run P is reported in rf_summary.tsv.
  Reference hit     t_rf = the channel-A hit of the event nearest in time to the row's t_m (same as C++
                    tdc::RfTimeData::nearest_a / prad2::ClusterDeltaRf; t_m in the FADC window frame,
                    t_rf in the TDC window frame, as in the dumps).  Vectorised with a keyed
                    searchsorted over the (event, time)-sorted A hits (equivalent to a nan-padded matrix
                    nanargmin).
  Fold / phase      fold(x) = x - P * round(x / P) in (-P/2, P/2];  phi = fold(t_m - t_rf).  Time
                    calibration (tcal): phi_tcal = fold(t_m - off_m - t_rf), off_m = per-module offset of
                    database/hycal_time_offsets/time_calib.json (plot_module_timing.load_offsets, geometry
                    of run 25320).  off_m is a per-module constant, so the tcal histograms are the raw ones
                    translated by -off_m: all tcal outputs are obtained by translation (no separate fill).
  Pass 1            per (module, run): circular histogram of raw phi, NB = 80 bins over one period (bin
                    = P/80 = 0.050 ns, index from phi / P + 0.5); per run: pooled tcal phi histogram over
                    all PbWO4 rows, join counters, A-hit counts, spacing sums.
  Per-run peak      c_mr for runs with N_mr >= --min-run-n: 5-bin circularly smoothed argmax, then
                    fit_core (bg=False, smooth_ns=0.25) on the histogram unwrapped around the argmax
                    (x = argmax centre + (-P/2 .. P/2)); mu if the fit is ok, else the argmax centre.
                    Reference run of (module, type) = the run of that type with the largest N_mr;
                    shift s_mr = fold(c_mr - c_ref).  Runs with N_mr < --min-run-n are excluded from the
                    module's aligned panel (entries counted in N_excl).
  Pass 2            every row of an included run: phi_a = fold(phi - s_mr), filled as fold(phi_a - c_ref)
                    (= fold(phi - c_mr), computed in one step) in 80 bins over (-P/2, P/2] per (module,
                    type): the "aligned" histogram, centred on the module's own peak.  For page modules
                    also the 2D histogram (integral bin x aligned phase) of the vs-integral analysis.
  Aligned fit       fit_core (bg=False, smooth_ns=0.25) on the aligned histogram unwrapped around its own
                    smoothed argmax -> sigma, sigma_err.  RF offset of (module, type) = fold(c_ref +
                    mu_aligned) (raw), fold(offset_raw - off_m) (tcal), reported in (-P/2, P/2].  The
                    plotted x axis is the aligned phase shifted so that the fitted peak sits at the
                    printed offset (= c_ref + y modulo P).
  Unaligned         sigma_unaligned: the same fit of the plain sum over the included runs of the per-run
                    histograms (no shifts), unwrapped around its own argmax.
  Page modules      N (aligned entries) > --min-n in at least one run type; ordered by module id.  Panels
                    with N >= 200 are fitted, fewer are drawn with a note.
  Vs integral       fit integral eu = m_E / cal_factor (cal json of the dump's calib_file, by module name;
                    fallback rows and factor <= 0 excluded; as plot_module_timing_vs_integral.py).
                    Integral axis per type [0, 1.2 * E_BEAM / median PbWO4 cal_factor], NIB = 150 linear
                    bins.  Per integral bin with >= --min-bin entries: median, p16, p84 of the aligned
                    phase from the 80-bin distribution (linear interpolation of the cumulative),
                    err = 1.2533 (p84 - p16) / 2 / sqrt(n); weighted line fit of the medians vs mean
                    integral / 1000 (plot_module_timing_vs_integral.linfit) -> slope in ns per
                    1000 ADC*samples.
  Caveats           pulses from neighbouring bunches (4.008 ns apart) fold onto the same phase, so
                    pile-up/accidentals are not separated by the fold; the flat pedestal under the peak
                    contains them.  A phase that changes inside a run (split to split) is not corrected
                    (see the split offset span in rf_summary.tsv); where that span is large (> ~0.5 ns)
                    the aligned sigma of the type is not better than sigma_unaligned.

Outputs (--out-dir, default output/module_timing_20261025; <s> = [_e<int>][_test])
  module_timing_rf<s>.pdf, module_timing_rf_tcal<s>.pdf   front pages + one page per page module
  module_timing_rf<s>.tsv         per (module, type): N, N_excl, n_runs, offset_raw, offset_tcal, sigma,
                                  sigma_err, fit_ok, sigma_unaligned, fit_ok_unaligned, max_abs_run_shift
  rf_summary<s>.tsv               per run: join statistics, P, run-level tcal offset, run_R, jumps
  module_timing_rf_vs_integral<s>.pdf / .tsv   slope summary page + per page module 2D maps
  The log of the full run is run_rf.log (stdout, tee).

Usage:
  MPLBACKEND=Agg python3 plot_module_timing_rf.py [--dump-dir DIR] [--rf-dir DIR] [--out-dir DIR]
      [--jobs 12] [--min-n 10000] [--min-run-n 200] [--min-bin 50] [--trigger-mask 0xf00]
      [--min-e 0] [--runs R [R ...]]
"""

from __future__ import annotations

import argparse
import csv
import glob
import hashlib
import json
import multiprocessing as mp
import os
import re
import sys
import textwrap
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import plot_module_timing as PMT  # noqa: E402
from plot_module_timing import (LABEL, NMAX_JOBS, RUN_NAME, TAGS, TYPE_OF_RUN, TRIG_EXCL,  # noqa: E402
                                fit_core, text_page)
from plot_module_timing_vs_integral import E_BEAM, NIB, linfit  # noqa: E402

NB = 80                 # phase bins over one RF period
P_NOM = 4.00771         # ns, nominal RF period (fallback only)
N_SUB = 32              # recorded edges = RF divided by 32
KEY_BIG = 1.0e4         # ns, event key stride of the nearest-hit search (|t| < 1e3 ns)
MIN_FIT = 200           # panels with fewer entries are not fitted
JUMP_NS = 0.5           # |s_mr| above this counts as a run-to-run phase jump
VETO_TSV = "module_timing_veto.tsv"


def fold(x, P):
    """x - P * round(x / P): (-P/2, P/2]."""
    return x - P * np.round(x / P)


# ---- workers ----------------------------------------------------------------

def prep(dp, rp, cfg):
    """Join one dump split to its RF split and build, for every selected row, (pm, t_m, phi, t_rf)."""
    d = np.load(dp)
    names = d["module_names"]
    if hashlib.md5("\n".join(map(str, names)).encode()).hexdigest() != cfg["names_hash"]:
        raise RuntimeError(f"{dp}: module_names differ from the geometry")
    ev, trig = d["ev_num"], d["trig"]
    r = np.load(rp)
    rf_ev, off = r["ev_num"], r["rf_off"].astype(np.int64)
    rf_t, rf_ch = r["rf_t"].astype(np.float64), r["rf_ch"]
    n_rf = len(rf_ev)
    pos = np.minimum(np.searchsorted(rf_ev, ev), n_rf - 1)
    matched = rf_ev[pos] == ev
    ha = np.flatnonzero(rf_ch == 0)
    ev_a = np.repeat(np.arange(n_rf), np.diff(off))[ha]
    ta = rf_t[ha]
    o = np.lexsort((ta, ev_a))
    ev_a, ta = ev_a[o], ta[o]
    nA = np.bincount(ev_a, minlength=n_rf)
    mrf = np.zeros(n_rf, bool)
    mrf[pos[matched]] = True
    same = (ev_a[1:] == ev_a[:-1]) & mrf[ev_a[1:]]
    dd = np.diff(ta)[same]
    dd = dd[(dd > 120.0) & (dd < 136.0)]
    P = float(dd.mean() / N_SUB) if dd.size else P_NOM
    sel = ((trig & np.uint32(cfg["mask"])) != 0) & ((trig & np.uint32(TRIG_EXCL)) == 0) & matched & (nA[pos] >= 1)
    out = dict(run=int(d["run"]), split=int(d["split"]), n_dump=len(ev), n_matched=int(matched.sum()),
               n_sel=int(sel.sum()), a_hits=int(nA[pos[matched]].sum()), sp_sum=float(dd.sum()),
               sp_n=int(dd.size), P=P, trig_equal=bool(np.array_equal(trig[matched], r["trig"][pos[matched]])))
    m_ev, m_idx = d["m_ev"].astype(np.int64), d["m_idx"].astype(np.int64)
    m_t, m_E = d["m_t"].astype(np.float64), d["m_E"].astype(np.float64)
    pm_all = cfg["pmap"][m_idx]
    rk = sel[m_ev] & (pm_all >= 0) & np.isfinite(m_t)
    if cfg["min_e"] > 0:
        rk &= m_E >= cfg["min_e"]
    ri = np.flatnonzero(rk)
    ev_row = pos[m_ev[ri]]
    t = m_t[ri]
    K = ev_a * KEY_BIG + ta
    i = np.searchsorted(K, ev_row * KEY_BIG + t)
    il, ih = np.clip(i - 1, 0, len(K) - 1), np.clip(i, 0, len(K) - 1)
    dl = np.where(ev_a[il] == ev_row, np.abs(ta[il] - t), np.inf)
    dh = np.where(ev_a[ih] == ev_row, np.abs(ta[ih] - t), np.inf)
    t_rf = np.where(dh < dl, ta[ih], ta[il])
    out.update(ri=ri, pm=pm_all[ri], t=t, t_rf=t_rf, phi=fold(t - t_rf, P), d=d, m_idx=m_idx[ri], sel=sel)
    return out


def _bin(phi, P):
    return np.clip(np.floor((phi / P + 0.5) * NB).astype(np.int64), 0, NB - 1)


def pass1(job):
    dp, rp, cfg = job
    o = prep(dp, rp, cfg)
    npw = cfg["npw"]
    pm = o["pm"]
    h = np.bincount(pm * NB + _bin(o["phi"], o["P"]), minlength=npw * NB).reshape(npw, NB)
    phit = fold(o["t"] - cfg["off_p"][pm] - o["t_rf"], o["P"])
    ht = np.bincount(_bin(phit, o["P"]), minlength=NB)
    for k in ("d", "ri", "pm", "t", "t_rf", "phi", "m_idx", "sel"):
        del o[k]
    o.update(h=h.astype(np.int32), ht=ht.astype(np.int64))
    return o


def pass2(job):
    dp, rp, cfg = job
    o = prep(dp, rp, cfg)
    npw, npg = cfg["npw"], cfg["npg"]
    d, pm, P = o["d"], o["pm"], o["P"]
    cm = cfg["c_run"][pm]
    k = np.isfinite(cm)
    y = fold(o["phi"] - np.where(k, cm, 0.0), P)
    by = _bin(y, P)
    h1 = np.bincount(pm[k] * NB + by[k], minlength=npw * NB).reshape(npw, NB)
    out = dict(run=o["run"], h1=h1.astype(np.int32))
    cal = {m["name"]: float(m["factor"]) for m in json.load(open(str(d["calib_file"])))}
    cf = np.array([cal.get(str(n), 0.0) for n in d["module_names"]])
    ri = o["ri"]
    fb = d["m_fb"][ri]
    cfr = cf[o["m_idx"]]
    pg = cfg["pg"][pm]
    eu = d["m_E"][ri].astype(np.float64) / np.where(cfr > 0, cfr, 1.0)
    bi = np.floor(eu / cfg["upper"] * NIB).astype(np.int64)
    k2 = k & (pg >= 0) & ~fb & (cfr > 0) & (bi >= 0) & (bi < NIB) & np.isfinite(eu)
    idx = (pg[k2] * NIB + bi[k2]) * NB + by[k2]
    u, c = np.unique(idx, return_counts=True)
    out.update(idx=u, cnt=c.astype(np.int32),
               eus=np.bincount(pg[k2] * NIB + bi[k2], weights=eu[k2], minlength=npg * NIB))
    return out


# ---- fitting ----------------------------------------------------------------

def unwrap(h, P):
    """Histogram h (NB bins over one period) unwrapped around its 5-bin smoothed argmax:
    returns argmax index, argmax centre, x (NB values, centre at index NB//2), counts."""
    nb = len(h)
    sm = sum(np.roll(h, k) for k in range(-2, 3))
    a = int(np.argmax(sm))
    j = np.arange(-(nb // 2), nb - nb // 2)
    cen = ((a + 0.5) / nb - 0.5) * P
    return a, cen, cen + j * (P / nb), np.asarray(h, float)[(a + j) % nb]


def fit_unwrapped(job):
    h, P = job
    a, cen, x, c = unwrap(h, P)
    r = fit_core(x, c, bg=False, smooth_ns=0.25)
    r["a"], r["cen"] = a, cen
    r["mu_u"] = r["mu"] if r["ok"] else cen
    return r


# ---- statistics helpers -----------------------------------------------------

def circ_median(v, P):
    """Median of values on a circle of period P (median of the offsets from the circular mean)."""
    v = v[np.isfinite(v)]
    if not v.size:
        return np.nan
    m0 = np.angle(np.mean(np.exp(2j * np.pi * v / P))) / (2 * np.pi) * P
    return float(fold(m0 + np.median(fold(v - m0, P)), P))


def circ_mad(v, P):
    v = v[np.isfinite(v)]
    if not v.size:
        return np.nan
    m = circ_median(v, P)
    return float(np.median(np.abs(fold(v - m, P))))


def run_pool(pool, func, jobs, label, chunksize=1):
    t0 = time.time()
    res = []
    for k, r in enumerate(pool.imap_unordered(func, jobs, chunksize=chunksize), 1):
        res.append(r)
        if k % 20 == 0 or k == len(jobs):
            print(f"  [{label}] {k}/{len(jobs)} files, {time.time() - t0:.0f} s", flush=True)
    return res


# ---- plotting ---------------------------------------------------------------

def timing_panel(ax, tag, R, ip, kind, P):
    """One aligned-phase panel of module ip (index among PbWO4 modules)."""
    ax.set_title(LABEL[tag], fontsize=10)
    ax.tick_params(labelsize=8)
    N = int(R["N"][ip])
    if R["nrun_type"] == 0:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.text(0.5, 0.5, "no data", transform=ax.transAxes, ha="center", va="center", fontsize=9, color="0.35")
        return
    if N < MIN_FIT:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.text(0.5, 0.5, f"N = {N} (too few)" + (f"\nN_excl = {int(R['N_excl'][ip])}" if R["N_excl"][ip] else ""),
                transform=ax.transAxes, ha="center", va="center", fontsize=9, color="0.35")
        return
    fit = R["fit"][ip]
    w = P / NB
    off_x = R["off_raw"][ip] if kind == "raw" else R["off_tcal"][ip]
    a, cen, x, c = unwrap(R["h1"][ip], P)
    # aligned histogram: unwrap around the aligned argmax, shift so the fitted peak is at the offset
    dx = off_x - fit["mu_u"]
    xs = x + dx
    ax.stairs(c, np.r_[xs - w / 2, xs[-1] + w / 2], color="C0", lw=0.9)
    txt = []
    if fit["ok"]:
        g = lambda xx: fit["amp"] * np.exp(-0.5 * ((xx - off_x) / fit["sigma"]) ** 2)
        xg = np.linspace(xs[0], xs[-1], 500)
        ax.plot(xg, g(xg), color="crimson", lw=0.8, ls="--", alpha=0.6)
        xw = np.linspace(fit["lo"] + dx, fit["hi"] + dx, 200)
        ax.plot(xw, g(xw), color="crimson", lw=1.8)
        txt = [f"mean (offset) = {off_x:.3f} ns", f"sigma = {fit['sigma']:.3f} +- {fit['sigma_err']:.3f} ns"]
    else:
        txt = [f"mean (offset) = {off_x:.3f} ns (argmax)", "fit failed"]
    su = R["sig_un"][ip]
    txt += [f"N = {N}", f"sigma_unaligned = {su:.3f} ns" if np.isfinite(su) else "sigma_unaligned: no fit",
            f"runs used {int(R['n_runs'][ip])}, N_excl = {int(R['N_excl'][ip])}"]
    ax.text(0.97, 0.96, "\n".join(txt), transform=ax.transAxes, ha="right", va="top", fontsize=7.5,
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85))
    ax.set_xlim(xs[0] - w / 2, xs[-1] + w / 2)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("t_module - t_RF (ns, folded)" + (" - tcal offset" if kind == "tcal" else ""), fontsize=8)
    ax.set_ylabel(f"entries / {w:.3f} ns", fontsize=8)


def render_timing(job):
    """Front pages + module pages of module_timing_rf[_tcal].pdf."""
    from matplotlib.backends.backend_pdf import PdfPages
    plt = PMT._mpl()
    path, kind, D = job
    t0 = time.time()
    with PdfPages(path) as pdf:
        for lines, title in D["text"][kind]:
            text_page(pdf, lines, title, size=7.2)
        front_figure(pdf, plt, D, kind)
        for ip in D["page_mods"]:
            fig, axs = plt.subplots(2, 2, figsize=(11, 8.5))
            g = D["gp"]
            fig.suptitle(f"{g['name'][ip]}   id {g['id'][ip]}   x = {g['x'][ip]:.2f} mm, y = {g['y'][ip]:.2f} mm"
                         + ("   (tcal)" if kind == "tcal" else ""), fontsize=13, weight="bold")
            for ax, tag in zip(axs.ravel(), TAGS):
                timing_panel(ax, tag, D["R"][tag], ip, kind, D["P"][tag])
            fig.tight_layout(rect=(0, 0, 1, 0.96))
            pdf.savefig(fig)
            plt.close(fig)
    return path, time.time() - t0


def front_figure(pdf, plt, D, kind):
    fig = plt.figure(figsize=(11, 8.5))
    gs = fig.add_gridspec(2, 4, hspace=0.38, wspace=0.35, left=0.07, right=0.97, top=0.9, bottom=0.08)
    allsig = np.concatenate([D["R"][t]["sigma"][D["R"][t]["N"] > D["min_n"]] for t in TAGS])
    allsig = allsig[np.isfinite(allsig)]
    hi = max(2.5, float(np.percentile(allsig, 99.5)) if allsig.size else 2.5)
    for j, tag in enumerate(TAGS):
        ax = fig.add_subplot(gs[0, j])
        R = D["R"][tag]
        s = R["sigma"][(R["N"] > D["min_n"])]
        s = s[np.isfinite(s)]
        ax.set_title(f"{LABEL[tag]}", fontsize=9)
        if s.size:
            ax.hist(s, bins=np.linspace(0, hi, 41), color="C0", alpha=0.85)
            ax.axvline(np.median(s), color="crimson", lw=1)
            ax.text(0.97, 0.95, f"{s.size} modules\nmedian {np.median(s):.3f} ns", transform=ax.transAxes,
                    ha="right", va="top", fontsize=7.5)
        ax.set_xlabel("sigma_RF aligned (ns)", fontsize=8)
        ax.tick_params(labelsize=7)
    ax = fig.add_subplot(gs[1, 1:3])
    v = D["veto"]
    if v is not None and len(v[0]):
        ax.plot(v[1], v[0], "o", ms=3, alpha=0.7, color="C0")
        lim = max(float(np.max(v[0])), float(np.max(v[1]))) * 1.08
        ax.plot([0, lim], [0, lim], "-", color="0.4", lw=0.8)
        ax.set_xlim(0, lim)
        ax.set_ylim(0, lim)
        ax.text(0.03, 0.97, f"{len(v[0])} modules (N_0p73 > {D['min_n']}, veto fit ok, RF fit ok)\n"
                f"median sigma_RF = {np.median(v[1]):.3f} ns\nmedian sigma_veto = {np.median(v[0]):.3f} ns",
                transform=ax.transAxes, ha="left", va="top", fontsize=8,
                bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85))
    ax.set_xlabel("sigma_RF, 0.73 GeV H2, aligned (ns)", fontsize=9)
    ax.set_ylabel("sigma_veto, 0.73 GeV H2 (ns)", fontsize=9)
    ax.tick_params(labelsize=8)
    ax.set_aspect("equal", adjustable="box")
    fig.suptitle("Per-module sigma of t_module - t_RF (aligned) and comparison with the veto reference", fontsize=11)
    pdf.savefig(fig)
    plt.close(fig)


def render_vs_integral(job):
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import LogNorm
    plt = PMT._mpl()
    path, D = job
    t0 = time.time()
    g = D["gp"]
    with PdfPages(path) as pdf:
        fig, axs = plt.subplots(4, 1, figsize=(11, 8.5), sharex=True)
        for ax, tag in zip(axs, TAGS):
            S = D["S"][tag]
            pm_ = D["page_mods"]
            s, se = S["slope"][pm_], S["slope_err"][pm_]
            ok = np.isfinite(s) & (S["N"][pm_] > D["min_n"])
            ax.errorbar(g["id"][pm_][ok], s[ok], yerr=se[ok], fmt="o", ms=2, lw=0.6, color="C0")
            ax.axhline(0, color="0.5", lw=0.6)
            med = float(np.median(s[ok])) if ok.any() else np.nan
            ax.set_title(f"{LABEL[tag]}: slope of the median RF phase vs fit integral ({ok.sum()} modules with N > {D['min_n']}), "
                         f"median {med:.3f} ns / 1000", fontsize=8.5)
            ax.set_ylabel("slope (ns / 1000)", fontsize=8)
            ax.tick_params(labelsize=7)
            if ok.any():
                lo, hi = np.percentile(s[ok], [1, 99])
                pad = 0.15 * (hi - lo if hi > lo else 1.0) + np.median(se[ok])
                ax.set_ylim(lo - pad, hi + pad)
        axs[-1].set_xlabel("module id", fontsize=8)
        fig.suptitle("Time walk in the RF frame: slope of median(phi - peak) vs fit integral", fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, 0.96))
        pdf.savefig(fig)
        plt.close(fig)
        for ip in D["page_mods"]:
            fig, axs = plt.subplots(2, 2, figsize=(11, 8.5))
            fig.suptitle(f"{g['name'][ip]}   id {g['id'][ip]}   x = {g['x'][ip]:.2f} mm, y = {g['y'][ip]:.2f} mm",
                         fontsize=13, weight="bold")
            for ax, tag in zip(axs.ravel(), TAGS):
                ax.set_title(LABEL[tag], fontsize=10)
                ax.tick_params(labelsize=8)
                S = D["S"][tag]
                pg = D["pg_of"][ip]
                N = int(S["N"][ip])
                if S["nrun_type"] == 0 or N < MIN_FIT:
                    ax.set_xticks([])
                    ax.set_yticks([])
                    ax.text(0.5, 0.5, "no data" if S["nrun_type"] == 0 else f"N = {N} (too few)",
                            transform=ax.transAxes, ha="center", va="center", fontsize=9, color="0.35")
                    continue
                ed, yb = D["edges"][tag], D["yedges"][tag]
                h = np.ma.masked_less(S["h2"][pg], 1)
                if h.count():
                    ax.pcolormesh(ed, yb, h.T, cmap="viridis", norm=LogNorm(vmin=1, vmax=max(h.max(), 2)),
                                  rasterized=True)
                pr = S["prof"][ip]
                if len(pr):
                    ax.plot(pr[:, 1], pr[:, 2], "o-", color="crimson", ms=2.5, lw=1.2, label="median")
                    ax.plot(pr[:, 1], pr[:, 3], "--", color="crimson", lw=0.8, label="16-84 %")
                    ax.plot(pr[:, 1], pr[:, 4], "--", color="crimson", lw=0.8)
                txt = [f"N = {N}"]
                if np.isfinite(S["slope"][ip]):
                    i0 = np.average(pr[:, 1], weights=1 / pr[:, 5] ** 2)
                    m0 = np.average(pr[:, 2], weights=1 / pr[:, 5] ** 2)
                    xs = np.array([pr[0, 1], pr[-1, 1]])
                    ax.plot(xs, m0 + S["slope"][ip] * (xs - i0) / 1000.0, "--", color="darkorange", lw=1.4,
                            label="lin. fit")
                    txt.append(f"slope = {S['slope'][ip]:+.3f} +- {S['slope_err'][ip]:.3f} ns/1000")
                    txt.append(f"{len(pr)} bins, chi2/ndf = {S['chi2'][ip]:.2f}")
                else:
                    txt.append("< 3 bins with enough entries, no fit")
                ax.text(0.97, 0.97, "\n".join(txt), transform=ax.transAxes, ha="right", va="top", fontsize=7.5,
                        bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85))
                ax.set_xlim(ed[0], ed[-1])
                ax.set_ylim(yb[0], yb[-1])
                ax.set_xlabel("fit integral A * I_int (ADC * samples)", fontsize=8)
                ax.set_ylabel("t_module - t_RF - peak (ns, folded)", fontsize=8)
                ax.legend(fontsize=6.5, loc="lower right", framealpha=0.8)
            fig.tight_layout(rect=(0, 0, 1, 0.96))
            pdf.savefig(fig)
            plt.close(fig)
    return path, time.time() - t0


# ---- analysis ---------------------------------------------------------------

def profile_from_hist(h2, edges_y, eus, min_bin):
    """h2 (NIB, NB) counts, eus (NIB,) sum of the integral per bin.  Returns rows
    (n, mean integral, median, p16, p84, median error) of the bins with >= min_bin entries."""
    n = h2.sum(1)
    ok = np.flatnonzero(n >= min_bin)
    if not ok.size:
        return np.zeros((0, 6))
    cum = np.cumsum(h2[ok], axis=1).astype(np.float64)
    w = edges_y[1] - edges_y[0]
    q = []
    for p in (0.5, 0.16, 0.84):
        tgt = p * n[ok]
        k = np.minimum((cum < tgt[:, None]).sum(1), NB - 1)
        prev = np.where(k > 0, cum[np.arange(len(ok)), np.maximum(k - 1, 0)], 0.0)
        cnt = h2[ok, k].astype(np.float64)
        q.append(edges_y[k] + (tgt - prev) / np.maximum(cnt, 1.0) * w)
    med, p16, p84 = q
    err = 1.2533 * 0.5 * (p84 - p16) / np.sqrt(n[ok])
    return np.column_stack([n[ok], eus[ok] / n[ok], med, p16, p84, err])


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dump-dir", default=str(HERE / "output/module_timing_20261025/dump"))
    ap.add_argument("--rf-dir", default=str(HERE / "output/module_timing_20261025/rf"))
    ap.add_argument("--out-dir", default=str(HERE / "output/module_timing_20261025"))
    ap.add_argument("--jobs", type=int, default=NMAX_JOBS, help=f"worker processes (max {NMAX_JOBS})")
    ap.add_argument("--min-n", type=int, default=10000, help="page / statistics threshold on the aligned N")
    ap.add_argument("--min-run-n", type=int, default=200, help="minimum entries of a (module, run) peak")
    ap.add_argument("--min-bin", type=int, default=50, help="minimum entries of an integral bin")
    ap.add_argument("--trigger-mask", type=lambda s: int(s, 0), default=0xF00)
    ap.add_argument("--min-e", type=float, default=0.0, help="module pulse energy cut (MeV, 0 = none)")
    ap.add_argument("--runs", type=int, nargs="+", default=None, help="subset of runs (outputs get _test)")
    ap.add_argument("--geometry-run", type=int, default=25320)
    ap.add_argument("--time-calib-file", default=str(PMT.TIME_CALIB))
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    t_start = time.time()
    nproc = max(1, min(args.jobs, NMAX_JOBS))
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    suf = (f"_e{int(round(args.min_e))}" if args.min_e > 0 else "") + ("_test" if args.runs else "")

    files = sorted(glob.glob(os.path.join(args.dump_dir, "*.npz")))
    jobs_files = []
    for fp in files:
        mt = re.search(RUN_NAME, fp)
        run = int(mt.group(1)) if mt else None
        if run not in TYPE_OF_RUN or (args.runs and run not in args.runs):
            continue
        rp = os.path.join(args.rf_dir, os.path.basename(fp))
        if not os.path.exists(rp):
            sys.exit(f"missing RF file {rp}")
        jobs_files.append((fp, rp, run))
    if not jobs_files:
        sys.exit(f"no usable npz in {args.dump_dir}")
    runs_of = {t: sorted({r for _, _, r in jobs_files if TYPE_OF_RUN[r] == t}) for t in TAGS}
    all_runs = [r for t in TAGS for r in runs_of[t]]
    print(f"{len(jobs_files)} dump files, runs {all_runs}, min-e {args.min_e}, min-n {args.min_n}, "
          f"min-run-n {args.min_run_n}", flush=True)

    # geometry (parent only), names check, time calibration offsets
    geo = PMT.load_geometry(args.geometry_run)
    d0 = np.load(jobs_files[0][0])
    names = d0["module_names"]
    if [str(n) for n in names] != [str(n) for n in geo["name"]]:
        sys.exit("module_names of the dump differ from the hycal geometry")
    pidx = np.flatnonzero(geo["pwo4"])
    npw = len(pidx)
    pmap = np.full(len(names), -1, np.int64)
    pmap[pidx] = np.arange(npw)
    gp = {k: geo[k][pidx] for k in ("name", "id", "x", "y", "sx", "sy")}
    off, n_res, missing = PMT.load_offsets(args.time_calib_file, geo)
    print(f"time calibration: {n_res}/{npw} PbWO4 modules resolved"
          + (f"; offset 0 for {missing}" if missing else ""), flush=True)
    off_p = off[pidx]
    cfg = dict(mask=args.trigger_mask, min_e=args.min_e, npw=npw, pmap=pmap, off_p=off_p,
               names_hash=hashlib.md5("\n".join(map(str, names)).encode()).hexdigest())

    ctx = mp.get_context("spawn")
    pool = ctx.Pool(nproc)

    # ---- pass 1 ----
    t0 = time.time()
    r1 = run_pool(pool, pass1, [(fp, rp, cfg) for fp, rp, _ in jobs_files], "pass1")
    H = {r: np.zeros((npw, NB), np.int64) for r in all_runs}
    S1 = {r: dict(n_splits=0, n_dump=0, n_matched=0, n_sel=0, a_hits=0, sp_sum=0.0, sp_n=0,
                  ht=np.zeros(NB, np.int64), splits=[], trig_ok=True) for r in all_runs}
    for r in r1:
        s = S1[r["run"]]
        H[r["run"]] += r["h"]
        s["n_splits"] += 1
        for k in ("n_dump", "n_matched", "n_sel", "a_hits", "sp_sum", "sp_n"):
            s[k] += r[k]
        s["ht"] += r["ht"]
        s["splits"].append((r["split"], r["P"], r["ht"]))
        s["trig_ok"] &= r["trig_equal"]
    P_run = {r: S1[r]["sp_sum"] / S1[r]["sp_n"] / N_SUB if S1[r]["sp_n"] else P_NOM for r in all_runs}
    P_type = {t: (float(np.mean([P_run[r] for r in runs_of[t]])) if runs_of[t] else P_NOM) for t in TAGS}
    print(f"pass 1 done in {time.time() - t0:.0f} s; P per run: "
          + " ".join(f"{r}:{P_run[r]:.5f}" for r in all_runs), flush=True)
    if not all(S1[r]["trig_ok"] for r in all_runs):
        print("WARNING: trig differs between dump and rf for some matched events", flush=True)

    # ---- per (module, run) peaks ----
    t0 = time.time()
    Nmr = {r: H[r].sum(1) for r in all_runs}
    fit_jobs, fit_key = [], []
    for r in all_runs:
        for m in np.flatnonzero(Nmr[r] >= args.min_run_n):
            fit_jobs.append((H[r][m], P_run[r]))
            fit_key.append((r, m))
    fres = pool.map(fit_unwrapped, fit_jobs, chunksize=64)
    cmr = {r: np.full(npw, np.nan) for r in all_runs}
    n_fail = 0
    for (r, m), fr in zip(fit_key, fres):
        cmr[r][m] = fold(fr["mu_u"], P_run[r])
        n_fail += not fr["ok"]
    print(f"per-run peaks: {len(fit_jobs)} (module, run) fits, {n_fail} fell back to the argmax, "
          f"{time.time() - t0:.0f} s", flush=True)

    # reference runs, shifts, aligned N, excluded entries
    A = {}
    for t in TAGS:
        rs = runs_of[t]
        a = dict(rs=rs, n_run=len(rs))
        if rs:
            Nm = np.stack([Nmr[r] for r in rs])
            has = np.stack([np.isfinite(cmr[r]) for r in rs])
            ref = np.argmax(np.where(has, Nm, -1), axis=0)
            hasref = has.any(0)
            cm = np.stack([cmr[r] for r in rs])
            cref = np.where(hasref, cm[ref, np.arange(npw)], np.nan)
            shift = fold(cm - cref[None, :], P_type[t])
            a.update(Nm=Nm, has=has, ref=ref, hasref=hasref, cref=cref, shift=shift,
                     N=np.where(has, Nm, 0).sum(0), N_excl=np.where(has, 0, Nm).sum(0),
                     n_runs=has.sum(0),
                     maxshift=np.where(hasref, np.where(has, np.abs(shift), -1.0).max(0), np.nan))
        else:
            z = np.zeros(npw, np.int64)
            a.update(N=z, N_excl=z, n_runs=z, cref=np.full(npw, np.nan), maxshift=np.full(npw, np.nan),
                     hasref=np.zeros(npw, bool))
        A[t] = a
    page_mask = np.zeros(npw, bool)
    for t in TAGS:
        page_mask |= A[t]["N"] > args.min_n
    page_mods = np.flatnonzero(page_mask)
    page_mods = page_mods[np.argsort(gp["id"][page_mods], kind="stable")]
    npg = len(page_mods)
    pg_of = np.full(npw, -1, np.int64)
    pg_of[page_mods] = np.arange(npg)
    print(f"page modules (N > {args.min_n} in >= 1 type): {npg}", flush=True)

    # integral axis per type from the calibration json of the dumps
    upper, calibs = {}, {}
    for t in TAGS:
        if not runs_of[t]:
            continue
        cfs = sorted({str(np.load(fp)["calib_file"]) for fp, _, r in jobs_files if TYPE_OF_RUN[r] == t})
        calibs[t] = cfs
        cal = {m["name"]: float(m["factor"]) for m in json.load(open(cfs[0]))}
        f_ = np.array([cal.get(str(n), 0.0) for n in gp["name"]])
        medf = float(np.median(f_[f_ > 0]))
        upper[t] = 1.2 * E_BEAM[t] / medf
        print(f"  {LABEL[t]}: calib {os.path.basename(cfs[0])}" + (f" (+{len(cfs) - 1} more files)" if len(cfs) > 1 else "")
              + f", median PbWO4 factor {medf:.5f}, integral axis [0, {upper[t]:.0f}]", flush=True)

    # ---- pass 2 ----
    t0 = time.time()
    jobs2 = []
    for fp, rp, run in jobs_files:
        t = TYPE_OF_RUN[run]
        c2 = dict(cfg, c_run=cmr[run].copy(), npg=npg, pg=pg_of, upper=upper[t])
        jobs2.append((fp, rp, c2))
    H1 = {t: np.zeros((npw, NB), np.int64) for t in TAGS}
    H2 = {t: np.zeros(npg * NIB * NB, np.int64) for t in TAGS}
    EUS = {t: np.zeros(npg * NIB) for t in TAGS}
    for r in run_pool(pool, pass2, jobs2, "pass2"):
        t = TYPE_OF_RUN[r["run"]]
        H1[t] += r["h1"]
        H2[t][r["idx"]] += r["cnt"]
        EUS[t] += r["eus"]
    print(f"pass 2 done in {time.time() - t0:.0f} s", flush=True)
    for t in TAGS:
        if runs_of[t]:
            assert np.array_equal(H1[t].sum(1), A[t]["N"]), f"aligned N mismatch for {t}"
        H2[t] = H2[t].reshape(npg, NIB, NB)
        EUS[t] = EUS[t].reshape(npg, NIB)

    # ---- fits ----
    t0 = time.time()
    R, S = {}, {}
    for t in TAGS:
        a = A[t]
        P = P_type[t]
        Rt = dict(N=a["N"], N_excl=a["N_excl"], n_runs=a["n_runs"], nrun_type=a["n_run"], h1=H1[t],
                  maxshift=a["maxshift"], cref=a["cref"])
        Hu = np.zeros((npw, NB), np.int64)
        if a["n_run"]:
            for k, r in enumerate(a["rs"]):
                Hu += np.where(a["has"][k][:, None], H[r], 0)
        todo = np.flatnonzero(a["N"] >= MIN_FIT)
        fa = pool.map(fit_unwrapped, [(H1[t][m], P) for m in todo], chunksize=64)
        fu = pool.map(fit_unwrapped, [(Hu[m], P) for m in todo], chunksize=64)
        blank = dict(ok=False, mu_u=np.nan, sigma=np.nan, sigma_err=np.nan, lo=np.nan, hi=np.nan, amp=np.nan)
        fit = [blank] * npw
        sig_un, ok_un = np.full(npw, np.nan), np.zeros(npw, int)
        for m, x, y in zip(todo, fa, fu):
            fit[m] = x
            if y["ok"]:
                sig_un[m], ok_un[m] = y["sigma"], 1
        mu_u = np.array([f_["mu_u"] for f_ in fit])
        Rt.update(fit=fit, sig_un=sig_un, ok_un=ok_un,
                  sigma=np.array([f_["sigma"] if f_["ok"] else np.nan for f_ in fit]),
                  sigma_err=np.array([f_["sigma_err"] if f_["ok"] else np.nan for f_ in fit]),
                  ok=np.array([int(f_["ok"]) for f_ in fit]))
        Rt["off_raw"] = fold(a["cref"] + mu_u, P)
        Rt["off_tcal"] = fold(Rt["off_raw"] - off_p, P)
        R[t] = Rt
    print(f"aligned / unaligned fits done in {time.time() - t0:.0f} s", flush=True)

    # ---- vs integral ----
    t0 = time.time()
    for t in TAGS:
        P = P_type[t]
        yed = (np.arange(NB + 1) / NB - 0.5) * P
        ed = np.linspace(0.0, upper.get(t, 1.0), NIB + 1)
        St = dict(nrun_type=A[t]["n_run"], N=np.zeros(npw, np.int64), slope=np.full(npw, np.nan),
                  slope_err=np.full(npw, np.nan), chi2=np.full(npw, np.nan), nbins=np.zeros(npw, int),
                  prof={}, h2=H2[t].astype(np.int32))
        for ip in page_mods:
            g = pg_of[ip]
            St["N"][ip] = int(H2[t][g].sum())
            pr = profile_from_hist(H2[t][g], yed, EUS[t][g], args.min_bin)
            St["prof"][ip] = pr
            St["nbins"][ip] = len(pr)
            if len(pr) >= 3:
                s, se, chi = linfit(pr[:, 1] / 1000.0, pr[:, 2], pr[:, 5])
                St["slope"][ip], St["slope_err"][ip], St["chi2"][ip] = s, se, chi
        S[t] = St
    print(f"vs-integral profiles done in {time.time() - t0:.0f} s", flush=True)

    # ---- per-run summary ----
    srows = []
    for r in all_runs:
        t = TYPE_OF_RUN[r]
        s = S1[r]
        P = P_run[r]
        x = (np.arange(NB) + 0.5) / NB * P - P / 2
        fr = fit_unwrapped((s["ht"], P))
        Rr = float(np.abs(np.sum(s["ht"] * np.exp(2j * np.pi * x / P))) / max(s["ht"].sum(), 1))
        k = runs_of[t].index(r)
        a = A[t]
        cmpd = a["has"][k] & (a["ref"] != k)
        njump = int((cmpd & (np.abs(a["shift"][k]) > JUMP_NS)).sum())
        # split-to-split spread of the pooled tcal peak inside the run
        sp = [(sn, ht_) for sn, _, ht_ in sorted(s["splits"], key=lambda z: z[0]) if ht_.sum() >= 2000]
        so = [fold(f_["mu_u"] - fr["mu_u"], P) for f_ in pool.map(fit_unwrapped, [(h_, P) for _, h_ in sp])]
        span = float(max(so) - min(so)) if so else np.nan
        srows.append(dict(run=r, type=t, n_splits=s["n_splits"], n_dump_events=s["n_dump"],
                          n_matched=s["n_matched"], match_rate=s["n_matched"] / max(s["n_dump"], 1),
                          n_sel_events=s["n_sel"], rfA_hits_per_event=s["a_hits"] / max(s["n_matched"], 1),
                          P_ns=P, run_offset_tcal=float(fold(fr["mu_u"], P)), run_R=Rr,
                          n_mod_compared=int(cmpd.sum()), n_mod_jump=njump, split_offset_span_tcal=span,
                          n_splits_in_span=len(so)))
    pool.close()
    pool.join()
    cols = ["run", "type", "n_splits", "n_dump_events", "n_matched", "match_rate", "n_sel_events",
            "rfA_hits_per_event", "P_ns", "run_offset_tcal", "run_R", "n_mod_compared", "n_mod_jump",
            "split_offset_span_tcal", "n_splits_in_span"]
    fmt = lambda v: (f"{v:.6g}" if isinstance(v, float) else str(v))
    with open(out_dir / f"rf_summary{suf}.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for q in srows:
            fh.write("\t".join(fmt(q[c]) for c in cols) + "\n")
    print("\nrun summary:")
    print("  " + "  ".join(f"{c[:11]:>11}" for c in cols))
    for q in srows:
        print("  " + "  ".join(f"{fmt(q[c])[:11]:>11}" for c in cols))

    # ---- TSVs ----
    with open(out_dir / f"module_timing_rf{suf}.tsv", "w") as fh:
        fh.write("module\tid\ttype\tN\tN_excl\tn_runs\toffset_raw\toffset_tcal\tsigma\tsigma_err\tfit_ok\t"
                 "sigma_unaligned\tfit_ok_unaligned\tmax_abs_run_shift\n")
        for m in range(npw):
            for t in TAGS:
                q = R[t]
                fh.write("\t".join([str(gp["name"][m]), str(gp["id"][m]), t, str(int(q["N"][m])),
                                    str(int(q["N_excl"][m])), str(int(q["n_runs"][m])),
                                    fmt(float(q["off_raw"][m])), fmt(float(q["off_tcal"][m])),
                                    fmt(float(q["sigma"][m])), fmt(float(q["sigma_err"][m])), str(int(q["ok"][m])),
                                    fmt(float(q["sig_un"][m])), str(int(q["ok_un"][m])),
                                    fmt(float(q["maxshift"][m]))]) + "\n")
    with open(out_dir / f"module_timing_rf_vs_integral{suf}.tsv", "w") as fh:
        fh.write("module\tid\ttype\tN\tn_bins\tslope_ns_per_1000\tslope_err\tchi2_ndf\n")
        for m in page_mods:
            for t in TAGS:
                q = S[t]
                fh.write("\t".join([str(gp["name"][m]), str(gp["id"][m]), t, str(int(q["N"][m])),
                                    str(int(q["nbins"][m])), fmt(float(q["slope"][m])),
                                    fmt(float(q["slope_err"][m])), fmt(float(q["chi2"][m]))]) + "\n")

    # ---- veto comparison (read only) ----
    veto = None
    vp = out_dir / VETO_TSV
    if vp.exists() and "0p73" in runs_of and runs_of["0p73"]:
        rows = list(csv.DictReader(open(vp), delimiter="\t"))
        byname = {str(n): i for i, n in enumerate(gp["name"])}
        sv, sr = [], []
        for q in rows:
            if int(q["N_0p73"]) > args.min_n and int(q["fit_ok_0p73"]) == 1 and q["module"] in byname:
                m = byname[q["module"]]
                if np.isfinite(R["0p73"]["sigma"][m]):
                    sv.append(float(q["sigma_0p73"]))
                    sr.append(float(R["0p73"]["sigma"][m]))
        veto = (np.array(sv), np.array(sr))
        print(f"veto scatter: {len(sv)} modules; median sigma_RF {np.median(sr):.3f}, "
              f"median sigma_veto {np.median(sv):.3f}", flush=True)

    # ---- front pages text ----
    text = {k: front_text(k, args, A, R, S, gp, off_p, P_type, runs_of, srows, page_mods, suf)
            for k in ("raw", "tcal")}
    D = dict(gp=gp, page_mods=page_mods, R=R, P=P_type, text=text, veto=veto, min_n=args.min_n)
    D2 = dict(gp=gp, page_mods=page_mods, S=S, pg_of=pg_of, min_n=args.min_n,
              edges={t: np.linspace(0, upper.get(t, 1.0), NIB + 1) for t in TAGS},
              yedges={t: (np.arange(NB + 1) / NB - 0.5) * P_type[t] for t in TAGS})

    # ---- PDFs (three in parallel) ----
    t0 = time.time()
    jobs = [("t", (str(out_dir / f"module_timing_rf{suf}.pdf"), "raw", D)),
            ("t", (str(out_dir / f"module_timing_rf_tcal{suf}.pdf"), "tcal", D)),
            ("v", (str(out_dir / f"module_timing_rf_vs_integral{suf}.pdf"), D2))]
    with ctx.Pool(3) as p3:
        res = [p3.apply_async(render_timing if k == "t" else render_vs_integral, (j,)) for k, j in jobs]
        for r_ in res:
            path, dt = r_.get()
            print(f"wrote {path} ({dt:.0f} s)", flush=True)
    print(f"PDFs done in {time.time() - t0:.0f} s; total wall time {time.time() - t_start:.0f} s", flush=True)


def front_text(kind, args, A, R, S, gp, off_p, P_type, runs_of, srows, page_mods, suf):
    """Text pages (list of (lines, title)) of module_timing_rf[_tcal].pdf."""
    key = "off_tcal"      # outliers and the type median use the tcal offsets in both PDFs (raw offsets are
                          # module-dependent by construction: the board phases differ)
    f = PMT.f
    W = 128
    wrap = lambda s: textwrap.wrap(s, W, subsequent_indent="    ") or [""]
    L = []
    for s in [
        "phi = t_module - t_RF folded modulo the RF period P; t_RF = channel-A RF hit nearest in time to t_module "
        "(C++ RfTimeData::nearest_a); fold(x) = x - P round(x/P), phi in (-P/2, P/2].",
        f"Events: (trig & {args.trigger_mask:#x}) != 0, no LMS/alpha (trig & 0x3000000), event matched in the RF file, "
        "at least one channel-A hit.  Rows: every dominant PbWO4 pulse"
        + (f" with m_E >= {args.min_e:g} MeV" if args.min_e > 0 else " (no energy cut)") + ".",
        "P per split = mean spacing of consecutive A hits (120-136 ns) / 32 (expected 4.00771 ns); used for that split's fold.",
        "Per (module, run): 80-bin circular histogram of phi; peak c_mr = fit_core (bg=False, smooth 0.25 ns) around the 5-bin "
        f"smoothed argmax if N_mr >= {args.min_run_n}.",
        "Reference run per (module, type) = run of the type with largest N_mr; every run's peak is shifted onto it (s_mr = "
        "fold(c_mr - c_ref)); runs with N_mr below the threshold are excluded (N_excl).",
        "Aligned histogram = sum over runs of fold(phi - c_mr), centred on the module's own peak; sigma = fit_core Gaussian core "
        "(+-2 sigma, bg=False, 0.05 ns bins).  sigma_unaligned = same fit of the plain sum of the raw per-run histograms.",
        "RF offset = fold(c_ref + mu_aligned) in (-P/2, P/2];  tcal offset = fold(offset - off_m) with off_m from "
        "time_calib.json;  tcal histograms are the raw ones translated by -off_m (no separate fill).",
        f"Page modules: aligned N > {args.min_n} in >= 1 type ({len(page_mods)}).  Panels with N >= {MIN_FIT} are fitted.  "
        f"This PDF shows the {'raw' if kind == 'raw' else 'tcal'} offsets; outliers always use the tcal offsets.",
        "CAVEAT: pulses from neighbouring bunches (4.008 ns apart) fold onto the same phase, so pile-up/accidentals are not "
        "separated by the fold; the flat pedestal under the peak contains them.",
        "CAVEAT: per-run alignment: each run's module peak is shifted onto the reference run's (largest-N run) because per-board "
        "phases jump between runs (up to ~P/2); sigma_unaligned is the plain sum.",
        "CAVEAT: a phase change inside a run (split to split) is not corrected; see split_offset_span_tcal below. Where the "
        "span is large (> ~0.5 ns) the aligned sigma of that type is not better than sigma_unaligned.",
    ]:
        L += wrap(s)
    L.append("")
    L.append(f"{'type':<13}{'modules':>8}{'med sigma':>10}{'med sig_un':>11}{'p16':>7}{'p84':>7}"
             f"{'MAD off(tcal)':>14}{'MAD off(raw)':>13}{'med off(tcal)':>14}")
    outl = {}
    for t in TAGS:
        q = R[t]
        P = P_type[t]
        sel = q["N"] > args.min_n
        n = int(sel.sum())
        if not runs_of[t]:
            L.append(f"{LABEL[t]:<13}{'no data':>8}")
            continue
        sg, su = q["sigma"][sel], q["sig_un"][sel]
        offs = q[key][sel]
        offt = q["off_tcal"][sel]
        medoff = circ_median(offs, P)
        L.append(f"{LABEL[t]:<13}{n:>8}{f(np.nanmedian(sg) if n else np.nan):>10}{f(np.nanmedian(su) if n else np.nan):>11}"
                 f"{f(PMT.pct(sg, 16)):>7}{f(PMT.pct(sg, 84)):>7}{f(circ_mad(offt, P)):>14}"
                 f"{f(circ_mad(q['off_raw'][sel], P)):>13}{f(medoff):>14}")
        ms = np.nanmedian(sg) if n else np.nan
        rows = []
        for m in np.flatnonzero(sel):
            why = []
            if np.isfinite(q["sigma"][m]) and q["sigma"][m] > 2 * ms:
                why.append("s")
            d = float(fold(q[key][m] - medoff, P))
            if abs(d) > 1.0:
                why.append("o")
            if why:
                rows.append((gp["name"][m], q[key][m], q["sigma"][m], "".join(why)))
        outl[t] = rows
    L.append("")
    L.append("Run table (P, pooled tcal offset, run_R = |<exp(2 pi i phi/P)>|, modules jumping > 0.5 ns vs the reference run):")
    L.append(f"{'run':>6} {'type':<5}{'P (ns)':>10}{'off_tcal':>10}{'R':>7}{'jump/compared':>15}{'split span':>11}{'match':>8}")
    for q in srows:
        L.append(f"{q['run']:>6} {q['type']:<5}{q['P_ns']:>10.5f}{q['run_offset_tcal']:>10.3f}{q['run_R']:>7.3f}"
                 f"{str(q['n_mod_jump']) + '/' + str(q['n_mod_compared']):>15}{f(q['split_offset_span_tcal'], 2):>11}"
                 f"{q['match_rate']:>8.4f}")
    pages = [(L, f"Module timing vs accelerator RF ({'raw' if kind == 'raw' else 'time-calibrated'})")]
    L2 = [f"Outliers (aligned N > {args.min_n}): sigma > 2 x type median (s) or circular |off_tcal - type median off_tcal| > 1 ns (o); "
          "listed as name(offset/sigma flags), max 25 per type.", ""]
    for t in TAGS:
        rows = outl.get(t)
        if rows is None:
            continue
        L2.append(f"{LABEL[t]}: {len(rows)} outliers")
        items = [f"{n}({o:+.2f}/{s:.2f} {w})" for n, o, s, w in rows[:25]]
        for i in range(0, len(items), 5):
            L2.append("   " + "  ".join(f"{x:<22}" for x in items[i:i + 5]))
        if len(rows) > 25:
            L2.append(f"   +{len(rows) - 25} more (module_timing_rf{suf}.tsv)")
        L2.append("")
    pages.append((L2, "Outliers"))
    return pages


if __name__ == "__main__":
    main()
