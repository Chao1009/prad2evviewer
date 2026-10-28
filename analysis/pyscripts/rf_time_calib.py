#!/usr/bin/env python3
"""rf_time_calib.py - per-run RF time calibration of the HyCal PbWO4 modules.

PURPOSE
    For every run and every PbWO4 module, the RF timing offset o_{m,r} such that
    fold(t_module - t_RF - o) is centred at zero, from module_timing_dump.py (module pulses) joined
    to rf_time_dump.py (RF TDC hits) exactly as plot_module_timing_rf.py (RF:) does, but per
    (run, module) instead of per (type, module).  One JSON per run in the database/hycal_rf_offsets
    schema (offset_ns per module, extra keys ignored by the C++ loader), two TSVs, a validation text
    and a before/after PDF per run.  The continuous deviation D of each module from the time_calib.json
    prior (which needs the +-P/2 branch) is resolved with the leave-one-out event median.  Plan:
    plans/plan_rf_calib_20261026.md.

DEFINITIONS (fold(x) = x - P*round(x/P) in (-P/2, P/2], P = RF period of the run)
  phi               fold(t_m - t_RF) per row, t_RF = nearest channel-A hit (RF.prep unchanged).  Rows: dominant
                    PbWO4 pulses of selected events, fallback rows (m_fb) dropped.
  Tiers             e100: rows with m_E >= E_MIN (100 MeV); allE: all rows.  Histograms (80 bins, RF._bin) per
                    (run, module, tier, half), half = parity of the split's index in the run's sorted split list.
  Estimator         RF.fit_unwrapped: 80-bin histogram unwrapped around its smoothed argmax, fit_core +-2 sigma
                    Gaussian core (argmax centre on failure, flag argmax; error then P/NB).
  o_{m,r}           c of tier e100 if N_e100 >= N_MIN (300); else c of allE if N_allE >= N_MIN (flag allE); else
                    fold(tau_m + Phi_r) (tier pred, if the module has a prior); else unresolved (not in the JSON).
  tau, rho, Phi     tau_m time_calib.json offset (W1123, W1124: 0, flag noprior).  rho = fold(o - tau) of the
                    measured modules, Phi_r = circ_median(rho), rho^r = fold(rho - Phi_r).
  Board phase       circ_median / circ_mad of rho of the measured modules of a (crate, slot) board, >= 4 modules.
  e (event median)  fit_core mean (bg=False, smooth 1.0) of the 600-bin histogram over [-60, 60) ns of
                    t_m - tau_m - t_loo (PMT.event_refs: leave-one-out median over the event's other pulses with
                    E > 10 MeV, in window, >= 3 mates); measured if N >= N_MIN, fit ok, sigma < 2 ns.  Pass A:
                    mates t_j - tau_j.  Pass B: mates t_j - tau_j - D1_j for modules with a firm pass-A branch.
  G_r, D            G_r = circ_median of fold(o - tau - e) over modules with e; its circ_mad is G_mad.
                    D = fold(o - tau - G_r) + k P, k in {-1, 0, 1} nearest to e; margin = second-best - best
                    distance; firm if margin >= 1 ns; no e -> k = 0, flag nobranch.  D1/G1 pass A, D2/G2 final.
  Bimodal           allE histogram: max over |fold(x - c)| > 0.75 ns of (5-bin mean - Gaussian core) >= ratio * amp
                    flags the module; a run with > 30 % flagged measured modules is unusable (quality bimodal).
  Drift             per-split pooled allE residual fold(phi - o_m) over the measured modules (pass B); span of the
                    split centres > 0.10 ns -> quality drift (usable).
  dhalf             fold(c_even - c_odd) of the tier used, both halves >= 150 entries.

Outputs (--out-dir, default output/rf_calib_20261026)
  <run>.json  rf_calib_modules.tsv  rf_calib_runs.tsv  rf_calib_validation.txt  rf_calib_<run>.pdf
  The log of the full run is run.log (stdout, tee).  The script exits 1 if a validation gate fails.

Usage:
  MPLBACKEND=Agg python3 rf_time_calib.py [--runs R [R ...]] [--jobs 12] [--out-dir DIR] [--no-pdf]
      [--bimodal-ratio 0.35]
"""

from __future__ import annotations

import argparse
import csv
import datetime
import glob
import hashlib
import json
import multiprocessing as mp
import os
import re
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import plot_module_timing as PMT  # noqa: E402
import plot_module_timing_rf as RF  # noqa: E402
from plot_module_timing import (EVT_LO, EVT_NB, NMAX_JOBS, RUN_NAME, TYPE_OF_RUN,  # noqa: E402
                                W_BIN, text_page)

REPO = HERE.parent.parent
DUMP_DIR = HERE / "output/module_timing_20261025/dump"
RF_DIR = HERE / "output/module_timing_20261025/rf"
OUT_DIR = HERE / "output/rf_calib_20261026"
HYCAL_MAP = REPO / "database" / "hycal_map.json"
XCHECK = HERE / "output/module_timing_20261025/rf_summary.tsv"
TIME_CALIB = PMT.TIME_CALIB
GEOMETRY_RUN = 25320

NB = RF.NB
E_MIN = 100.0            # MeV, plateau of the walk
N_MIN = 300              # entries for a measured module / event reference
N_HALF = 150             # entries per half for dhalf
REF_MIN_E, REF_MIN_N = 10.0, 3     # event reference mates: E > 10 MeV, >= 3 mates
SIG_E_MAX = 2.0          # ns, a measured event reference needs sigma_e below
FIRM = 1.0               # ns, branch margin
BIM_RATIO = 0.35         # bimodal residual / amplitude (calibrated, see rf_calib_runs.tsv)
BIM_FAR = 0.75           # ns, residual region |fold(x - c)| >
BIM_RUN_FRAC = 0.30      # run unusable above this fraction of flagged measured modules
DRIFT_SPAN = 0.10        # ns
MIN_BOARD = 4
MIN_SPLIT_N = 2000       # entries of a per-split pooled residual histogram to be fitted

X17_RUNS = (25800, 26138, 26230)
CENTRAL = ("W492 W493 W494 W495 W525 W526 W527 W528 W529 W559 W560 W563 W564 W593 W594 W597 W598 "
           "W627 W628 W629 W630 W631 W663 W664 W665 W697").split()
HALO = "W491 W496 W530 W632 W666 W661 W662".split()
FOUR = "W525 W528 W597 W697".split()
TIERS = ("e100", "allE", "pred", "unresolved")
FLAG_ORDER = ("allE", "pred", "noprior", "bimodal", "nobranch", "branch_weak", "branch_changed", "argmax",
              "run_unusable")

TIME_DEF = ("fold(t_module - t_rf - offset_ns) centred at 0 for pulses >= 100 MeV (tier allE: all pulses); "
            "t_module = matched-filter p_time of module_timing_dump 1.0 (20261025), templates {templates}, "
            "no per-module time offsets applied; t_rf = nearest channel-A RF hit; "
            "fold(x) = x - P*round(x/P), P = period_ns")


# ---- workers ----------------------------------------------------------------

def _evt_ref(d, o, tau_full, corr):
    """Per selected row: e_row = t_own - t_loo, ok mask (row in window, finite), E.  Own time t - tau,
    mates' times t - tau - corr (corr is zero in pass A)."""
    m_ev, m_idx = d["m_ev"].astype(np.int64), d["m_idx"].astype(np.int64)
    m_t, m_E, inw = d["m_t"].astype(np.float64), d["m_E"].astype(np.float64), d["m_inwin"]
    t_own = m_t - tau_full[m_idx]
    Eev = np.where(inw & np.isfinite(m_t), m_E, -1.0)
    sel = o["sel"]
    _, t_loo = PMT.event_refs(m_ev, t_own - corr[m_idx], Eev, sel, len(sel), REF_MIN_E, REF_MIN_N)
    ri = o["ri"]
    e = t_own[ri] - t_loo[ri]
    return e, inw[ri] & np.isfinite(e), m_E[ri]


def _event_hist(pm, e, ok, E, keep, npw):
    b = np.floor((np.where(ok, e, 0.0) - EVT_LO) / W_BIN).astype(np.int64)
    k = ok & keep & (b >= 0) & (b < EVT_NB)
    out = np.zeros((2, npw, EVT_NB), np.int32)
    for ti, tk in enumerate((k & (E >= E_MIN), k)):
        out[ti] = np.bincount(pm[tk] * EVT_NB + b[tk], minlength=npw * EVT_NB).reshape(npw, EVT_NB)
    return out


def pass_a(job):
    dp, rp, cfg, half = job
    o = RF.prep(dp, rp, cfg)
    d, npw, pm, P, phi = o["d"], cfg["npw"], o["pm"], o["P"], o["phi"]
    fb = d["m_fb"][o["ri"]]
    e, ok, E = _evt_ref(d, o, cfg["tau_full"], np.zeros(len(cfg["tau_full"])))
    keep = ~fb
    h = np.zeros((2, npw, NB), np.int32)
    b = RF._bin(phi, P)
    for ti, tk in enumerate((keep & (E >= E_MIN), keep)):
        h[ti] = np.bincount(pm[tk] * NB + b[tk], minlength=npw * NB).reshape(npw, NB)
    phit = RF.fold(o["t"] - cfg["off_p"][pm] - o["t_rf"], P)
    ht = np.bincount(RF._bin(phit, P), minlength=NB)
    return dict(run=o["run"], split=o["split"], half=half, P=P, n_dump=o["n_dump"], n_matched=o["n_matched"],
                n_sel=o["n_sel"], a_hits=o["a_hits"], sp_sum=o["sp_sum"], sp_n=o["sp_n"], trig_equal=o["trig_equal"],
                h=h, he=_event_hist(pm, e, ok, E, keep, npw), ht=ht.astype(np.int64),
                nrow=np.bincount(pm, minlength=npw), nfb=np.bincount(pm[fb], minlength=npw))


def pass_b(job):
    dp, rp, cfg, half, corr, off_m, c_ho = job
    o = RF.prep(dp, rp, cfg)
    d, npw, pm, P, phi = o["d"], cfg["npw"], o["pm"], o["P"], o["phi"]
    keep = ~d["m_fb"][o["ri"]]
    e, ok, E = _evt_ref(d, o, cfg["tau_full"], corr)
    k = keep & np.isfinite(off_m[pm])
    hd = np.bincount(RF._bin(RF.fold(phi[k] - off_m[pm[k]], P), P), minlength=NB)
    k = keep & np.isfinite(c_ho[pm])
    hh = np.bincount(RF._bin(RF.fold(phi[k] - c_ho[pm[k]], P), P), minlength=NB)
    return dict(run=o["run"], split=o["split"], half=half, P=P, he=_event_hist(pm, e, ok, E, keep, npw),
                hd=hd.astype(np.int64), hh=hh.astype(np.int64))


def fit_task(job):
    """One RF.fit_unwrapped; for bim=True also the bimodal residual ratio on the same histogram."""
    h, P, bim = job
    r = RF.fit_unwrapped((h, P))
    ok = bool(r["ok"])
    out = dict(ok=ok, c=float(RF.fold(r["mu_u"], P)), sigma=float(r["sigma"]) if ok else np.nan,
               err=float(r["mu_err"]) if ok else P / NB, ratio=np.nan)
    if bim and ok:
        x = ((np.arange(NB) + 0.5) / NB - 0.5) * P
        sm = sum(np.roll(np.asarray(h, float), k) for k in range(-2, 3)) / 5.0
        dist = RF.fold(x - out["c"], P)
        res = sm - r["amp"] * np.exp(-0.5 * (dist / r["sigma"]) ** 2)
        far = np.abs(dist) > BIM_FAR
        out["ratio"] = float(res[far].max() / r["amp"]) if far.any() and r["amp"] > 0 else np.nan
    return out


def efit_task(h):
    n = int(h.sum())
    if n < N_MIN:
        return (False, np.nan, np.nan, n)
    r = PMT.fit_core(PMT.centres(EVT_LO, W_BIN, EVT_NB), h, bg=False, smooth_ns=1.0)
    return (bool(r["ok"] and r["sigma"] < SIG_E_MAX), r["mu"], r["sigma"], n)


# ---- per-run estimation ------------------------------------------------------

def fold(x, P):
    return RF.fold(x, P)


def fit_run_jobs(r, Rn):
    """Fit tasks of a run: tier e100 / allE (both halves summed) for N >= N_MIN, plus the halves of the tier used."""
    P = Rn["P"]
    Hb = Rn["H"].sum(0)                      # (tier, npw, NB)
    ne, na = Hb[0].sum(1), Hb[1].sum(1)
    tu = np.where(ne >= N_MIN, 0, np.where(na >= N_MIN, 1, 2))
    Rn.update(n_e100=ne, n_allE=na, tier=tu)
    keys, jobs = [], []
    for m in np.flatnonzero(ne >= N_MIN):
        keys.append((r, m, 0, 0)); jobs.append((Hb[0, m], P, False))
    for m in np.flatnonzero(na >= N_MIN):
        keys.append((r, m, 1, 0)); jobs.append((Hb[1, m], P, True))
    for m in np.flatnonzero(tu <= 1):
        for hf in (0, 1):
            if Rn["H"][hf, tu[m], m].sum() >= N_HALF:
                keys.append((r, m, int(tu[m]), 1 + hf)); jobs.append((Rn["H"][hf, tu[m], m], P, False))
    return keys, jobs


def estimate_run(Rn, prior, tau, board, fits):
    """Fill offsets, errors, sigmas, dhalf, rho, Phi, board phases, bimodal ratios, walk bias of a run."""
    P, npw = Rn["P"], len(tau)
    tu = Rn["tier"]
    c = np.full((2, 3, npw), np.nan); er = np.full((2, 3, npw), np.nan); sg = np.full((2, 3, npw), np.nan)
    okf = np.ones((2, 3, npw), bool); ratio = np.full(npw, np.nan)
    for (ti, slot, m), f in fits.items():
        c[ti, slot, m], er[ti, slot, m], sg[ti, slot, m], okf[ti, slot, m] = f["c"], f["err"], f["sigma"], f["ok"]
        if ti == 1 and slot == 0:
            ratio[m] = f["ratio"]
    meas = tu <= 1
    sel = lambda a: np.where(tu == 0, a[0], np.where(tu == 1, a[1], np.nan))
    o = sel(c[:, 0]); err = sel(er[:, 0]); sig = sel(sg[:, 0])
    argmax = meas & ~np.where(tu == 0, okf[0, 0], okf[1, 0])
    n_used = np.where(tu == 0, Rn["n_e100"], np.where(tu == 1, Rn["n_allE"], 0))
    ch = np.stack([sel(c[:, 1]), sel(c[:, 2])])
    eh = np.stack([sel(er[:, 1]), sel(er[:, 2])])
    dhalf = fold(ch[0] - ch[1], P)
    dhalf_err = np.sqrt(eh[0] ** 2 + eh[1] ** 2)
    rho = np.where(meas & prior, fold(o - tau, P), np.nan)
    Phi = RF.circ_median(rho, P)
    rho_run = fold(rho - Phi, P)
    pred = (tu == 2) & prior
    o = np.where(pred, fold(tau + Phi, P), o)
    tier = np.where(pred, 2, np.where(tu == 2, 3, tu))
    bphase, bmad = np.full(npw, np.nan), np.full(npw, np.nan)
    for b in np.unique(board):
        s = (board == b) & np.isfinite(rho)
        if s.sum() >= MIN_BOARD:
            bphase[board == b] = RF.circ_median(rho[s], P)
            bmad[board == b] = RF.circ_mad(rho[s], P)
    both = np.isfinite(c[0, 0]) & np.isfinite(c[1, 0])
    walk = float(np.median(fold(c[1, 0][both] - c[0, 0][both], P))) if both.any() else np.nan
    Rn.update(tier=tier, o=o, err=err, sigma=sig, argmax=argmax, n_used=n_used, c_half=ch, dhalf=dhalf,
              dhalf_err=dhalf_err, rho=rho, Phi=Phi, rho_run=rho_run, bphase=bphase, bmad=bmad, ratio=ratio,
              walk=walk, meas=tier <= 1, c_all=c, n_fits=len(fits), n_argmax_fits=int((~okf[np.isfinite(c)]).sum()))


def branch(Rn, tau, prior, e):
    """G, D, margin for a run given the event reference e (nan = none).  Returns dict."""
    P, o = Rn["P"], Rn["o"]
    meas = Rn["meas"]
    use = meas & prior & np.isfinite(e)
    v = np.where(use, fold(o - tau - np.where(use, e, 0.0), P), np.nan)
    G = RF.circ_median(v, P)
    base = np.where(meas & prior, fold(o - tau - G, P), np.nan)
    dist = np.abs(base[None, :] + P * np.array([-1.0, 0.0, 1.0])[:, None] - np.where(use, e, 0.0)[None, :])
    k = np.argmin(dist, axis=0) - 1
    ds = np.sort(dist, axis=0)
    D = np.where(use, base + P * k, base)
    margin = np.where(use, ds[1] - ds[0], np.nan)
    return dict(G=G, G_mad=RF.circ_mad(v, P), D=D, margin=margin, use=use, e=np.where(use, e, np.nan),
                firm=use & (margin >= FIRM), nobranch=meas & prior & ~use)


def efit_arrays(HE, Rn, pool):
    """e fits of the measured modules (tier used); returns e (nan = not measured), sigma_e."""
    npw = len(Rn["tier"])
    ms = np.flatnonzero(Rn["meas"])
    res = pool.map(efit_task, [HE[Rn["tier"][m], m] for m in ms], chunksize=64)
    e, se = np.full(npw, np.nan), np.full(npw, np.nan)
    for m, (ok, mu, s, n) in zip(ms, res):
        if ok:
            e[m], se[m] = mu, s
    return e, se


def tsv_fmt(v, nd=5):
    if isinstance(v, (bool, np.bool_)):
        return str(int(v))
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, str):
        return v
    return "nan" if v is None or not np.isfinite(v) else f"{v:.{nd}f}"


def flags_of(Rn, m):
    f = []
    t = int(Rn["tier"][m])
    if t == 1:
        f.append("allE")
    if t == 2:
        f.append("pred")
    if not Rn["prior"][m]:
        f.append("noprior")
    if Rn["bimodal"][m]:
        f.append("bimodal")
    if Rn["nobranch"][m]:
        f.append("nobranch")
    if Rn["branch_weak"][m]:
        f.append("branch_weak")
    if Rn["branch_changed"][m]:
        f.append("branch_changed")
    if Rn["argmax"][m]:
        f.append("argmax")
    if not Rn["usable"]:
        f.append("run_unusable")
    return [x for x in FLAG_ORDER if x in f]


def jnum(v, nd):
    return None if v is None or not np.isfinite(v) else round(float(v), nd)


def round_off(x, P):
    r = round(float(fold(x, P)), 4)
    if r <= -P / 2:
        r = round(r + P, 4)
    elif r > P / 2:
        r = round(r - P, 4)
    return r


def write_json(path, r, Rn, gp, created):
    P = Rn["P"]
    mods = []
    if Rn["usable"]:
        for m in np.argsort(gp["id"], kind="stable"):
            t = int(Rn["tier"][m])
            if t > 2:
                continue
            mods.append({"name": str(gp["name"][m]), "offset_ns": round_off(Rn["o"][m], P),
                         "module_id": int(gp["id"][m]),
                         "offset_error_ns": jnum(Rn["err"][m], 4) if t < 2 else None,
                         "sigma_ns": jnum(Rn["sigma"][m], 4) if t < 2 else None,
                         "n": int(Rn["n_used"][m]) if t < 2 else None,
                         "tier": TIERS[t], "flags": ",".join(flags_of(Rn, m))})
    j = {"default": 0.0, "modules": mods, "run": int(r), "run_type": TYPE_OF_RUN[r], "period_ns": round(P, 5),
         "usable": bool(Rn["usable"]), "quality": Rn["quality"], "split_span_ns": jnum(Rn["split_span"], 4),
         "time_definition": TIME_DEF.format(templates=Rn["templates"]), "created": created}
    with open(path, "w") as fh:
        json.dump(j, fh, indent=1, allow_nan=False)
        fh.write("\n")
    return len(mods)


# ---- PDFs --------------------------------------------------------------------

def _hist_axes(ax, h, P, title, fit_text=True):
    a, cen, x, c = RF.unwrap(h, P)
    w = P / NB
    ax.stairs(c, np.r_[x - w / 2, x[-1] + w / 2], color="C0", lw=0.9)
    r = RF.fit_unwrapped((h, P))
    txt = [f"N = {int(np.sum(h))}"]
    if r["ok"] and fit_text:
        g = lambda xx: r["amp"] * np.exp(-0.5 * ((xx - r["mu"]) / r["sigma"]) ** 2)
        xg = np.linspace(x[0], x[-1], 400)
        ax.plot(xg, g(xg), color="crimson", lw=0.8, ls="--", alpha=0.6)
        xw = np.linspace(r["lo"], r["hi"], 100)
        ax.plot(xw, g(xw), color="crimson", lw=1.8)
        txt += [f"centre = {fold(r['mu'], P):.3f} ns", f"sigma = {r['sigma']:.3f} ns"]
    ax.set_title(title, fontsize=9)
    ax.tick_params(labelsize=8)
    ax.set_xlabel("ns (folded, unwrapped at the peak)", fontsize=8)
    ax.set_ylabel(f"entries / {w:.3f} ns", fontsize=8)
    ax.set_ylim(bottom=0)
    ax.text(0.97, 0.96, "\n".join(txt), transform=ax.transAxes, ha="right", va="top", fontsize=7.5,
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85))
    return a, cen, x, r


def render_run(job):
    """rf_calib_<run>.pdf: text page, pooled phase (before / prior / after), map, drift, module panels."""
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Rectangle
    plt = PMT._mpl()
    path, r, D = job
    t0 = time.time()
    P, gp = D["P"], D["gp"]
    npage = 0
    with PdfPages(path) as pdf:
        text_page(pdf, D["text"], f"Run {r}: per-run RF calibration", size=7.6); npage += 1
        fig, axs = plt.subplots(1, 3, figsize=(11, 4.2))
        for ax, h, ttl in zip(axs, (D["h_raw"], D["h_tau"], D["h_after"]),
                              ("raw phi (no offsets)", "fold(phi - tau)  (prior only)",
                               "after: held-out residual (pass B)")):
            _hist_axes(ax, h, P, ttl)
        fig.suptitle(f"Run {r} ({D['type']}): pooled phase, allE rows of all PbWO4 modules", fontsize=11)
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig); npage += 1
        fig, axs = plt.subplots(1, 2, figsize=(11, 5.6))
        for ax, key, lim, ttl in zip(axs, ("rho_run", "dhalf"), (2.0, 0.3),
                                     ("rho^r = fold(rho - Phi_r)  (before, ns)", "dhalf = c_even - c_odd  (after, ns)")):
            v = D[key]
            col = plt.get_cmap("coolwarm")((np.clip(np.where(np.isfinite(v), v, 0.0), -lim, lim) + lim) / (2 * lim))
            col[~np.isfinite(v)] = (0.8, 0.8, 0.8, 1.0)
            pc = PatchCollection([Rectangle((gp["x"][i] - gp["sx"][i] / 2, gp["y"][i] - gp["sy"][i] / 2),
                                            gp["sx"][i], gp["sy"][i]) for i in range(len(v))],
                                 facecolors=col, edgecolors="none")
            ax.add_collection(pc)
            ax.set_xlim(gp["x"].min() - 20, gp["x"].max() + 20); ax.set_ylim(gp["y"].min() - 20, gp["y"].max() + 20)
            ax.set_aspect("equal"); ax.set_title(ttl, fontsize=9); ax.tick_params(labelsize=8)
            sm = plt.cm.ScalarMappable(cmap="coolwarm", norm=plt.Normalize(-lim, lim))
            fig.colorbar(sm, ax=ax, fraction=0.04, pad=0.02)
        fig.suptitle(f"Run {r}: HyCal map (grey: not measured)", fontsize=11)
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig); npage += 1
        if D["drift"] is not None:
            sp, cs = D["drift"]
            fig, ax = plt.subplots(figsize=(11, 5))
            ax.plot(sp, cs, "o-", color="C0")
            ax.axhline(np.median(cs), color="0.5", lw=0.8)
            ax.set_xlabel("split"); ax.set_ylabel("centre of the pooled residual (ns)")
            ax.set_title(f"Run {r}: per-split centre of fold(phi - o_m), span {D['split_span']:.3f} ns")
            pdf.savefig(fig); plt.close(fig); npage += 1
        mods = D["panel_mods"]
        for k0 in range(0, len(mods), 4):
            fig, axs = plt.subplots(2, 2, figsize=(11, 8.5))
            for ax, m in zip(axs.ravel(), mods[k0:k0 + 4]):
                t = int(D["tier"][m])
                tt = min(t, 1) if t < 3 else 1
                h = D["Hb"][tt, m]
                title = f"{gp['name'][m]}  id {gp['id'][m]}   tier {TIERS[t]}"
                a, cen, x, rr = _hist_axes(ax, h, P, title, fit_text=t < 2)
                if t == 2:
                    ax.axvline(cen + fold(D["o"][m] - cen, P), color="crimson", lw=1.2)
                f = lambda v, nd=3: "nan" if not np.isfinite(v) else f"{v:.{nd}f}"
                ax.text(0.03, 0.96, "\n".join([
                    f"offset {f(D['o'][m])}  rho^r {f(D['rho_run'][m])}",
                    f"dhalf {f(D['dhalf'][m])} +- {f(D['dhalf_err'][m])}",
                    f"e {f(D['e'][m])}  D {f(D['D'][m])}  margin {f(D['margin'][m], 2)}",
                    f"flags {D['flags'][m] or '-'}"]), transform=ax.transAxes, ha="left", va="top", fontsize=7.5,
                    family="monospace", bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85))
            for ax in axs.ravel()[len(mods[k0:k0 + 4]):]:
                ax.axis("off")
            fig.tight_layout(); pdf.savefig(fig); plt.close(fig); npage += 1
    return path, time.time() - t0, npage


# ---- main --------------------------------------------------------------------

def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", default=str(OUT_DIR))
    ap.add_argument("--dump-dir", default=str(DUMP_DIR), help="module_timing_dump npz directory")
    ap.add_argument("--rf-dir", default=str(RF_DIR), help="rf_time_dump npz directory (same file names)")
    ap.add_argument("--xcheck", default=str(XCHECK), help="rf_summary.tsv of the same dump (gate vi)")
    ap.add_argument("--jobs", type=int, default=NMAX_JOBS, help=f"worker processes (max {NMAX_JOBS})")
    ap.add_argument("--runs", type=int, nargs="+", default=None, help="subset of runs")
    ap.add_argument("--no-pdf", action="store_true")
    ap.add_argument("--bimodal-ratio", type=float, default=BIM_RATIO)
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    t_start = time.time()
    nproc = max(1, min(args.jobs, NMAX_JOBS))
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    dump_dir, rf_dir, xcheck = Path(args.dump_dir), Path(args.rf_dir), Path(args.xcheck)

    # ---- files, pairing, halves ----
    by_run = {}
    for fp in sorted(glob.glob(os.path.join(str(dump_dir), "*.npz"))):
        mt = re.search(RUN_NAME, fp)
        if not mt:
            continue
        run, split = int(mt.group(1)), int(mt.group(2))
        if run not in TYPE_OF_RUN or (args.runs and run not in args.runs):
            continue
        rp = os.path.join(str(rf_dir), os.path.basename(fp))
        if not os.path.exists(rp):
            sys.exit(f"missing RF file {rp}")
        by_run.setdefault(run, []).append((split, fp, rp))
    if not by_run:
        sys.exit(f"no usable npz in {dump_dir}")
    runs = sorted(by_run)
    for r in runs:
        by_run[r].sort()
    print(f"{sum(len(v) for v in by_run.values())} dump files, runs {runs}", flush=True)

    # ---- geometry, prior, boards ----
    geo = PMT.load_geometry(GEOMETRY_RUN)
    d0 = np.load(by_run[runs[0]][0][1])
    names = d0["module_names"]
    if [str(n) for n in names] != [str(n) for n in geo["name"]]:
        sys.exit("module_names of the dump differ from the hycal geometry")
    pidx = np.flatnonzero(geo["pwo4"])
    npw = len(pidx)
    pmap = np.full(len(names), -1, np.int64)
    pmap[pidx] = np.arange(npw)
    gp = {k: geo[k][pidx] for k in ("name", "id", "x", "y", "sx", "sy")}
    tau_full, n_res, missing = PMT.load_offsets(str(TIME_CALIB), geo)
    tau = tau_full[pidx]
    prior = ~np.isin(gp["name"], list(missing))
    print(f"time calibration: {n_res}/{npw} PbWO4 modules resolved" + (f"; no prior: {missing}" if missing else ""),
          flush=True)
    bd = {m["n"]: (int(m["daq"]["crate"]), int(m["daq"]["slot"])) for m in json.load(open(HYCAL_MAP))}
    lost = [str(n) for n in gp["name"] if str(n) not in bd]
    if lost:
        sys.exit(f"PbWO4 modules missing from hycal_map.json: {lost[:10]}")
    bkeys = [bd[str(n)] for n in gp["name"]]
    board = np.array([c * 100 + s for c, s in bkeys])
    cfg = dict(mask=0xF00, min_e=0.0, npw=npw, pmap=pmap, off_p=tau, tau_full=tau_full,
               names_hash=hashlib.md5("\n".join(map(str, names)).encode()).hexdigest())

    ctx = mp.get_context("spawn")
    pool = ctx.Pool(nproc)

    # ---- pass A ----
    t0 = time.time()
    jobs = []
    for r in runs:
        for i, (sp, fp, rp) in enumerate(by_run[r]):
            jobs.append((fp, rp, cfg, i % 2))
    RN = {r: dict(H=np.zeros((2, 2, npw, NB), np.int64), HE=np.zeros((2, npw, EVT_NB), np.int64),
                  ht=np.zeros(NB, np.int64), sp_sum=0.0, sp_n=0, n_dump=0, n_matched=0, n_sel=0, a_hits=0,
                  nrow=np.zeros(npw, np.int64), nfb=np.zeros(npw, np.int64), splits=[], trig_ok=True,
                  templates=str(np.load(by_run[r][0][1])["templates"])) for r in runs}
    for q in RF.run_pool(pool, pass_a, jobs, "passA"):
        s = RN[q["run"]]
        s["H"][q["half"]] += q["h"]
        s["HE"] += q["he"]
        s["ht"] += q["ht"]
        for k in ("sp_sum", "sp_n", "n_dump", "n_matched", "n_sel", "a_hits"):
            s[k] += q[k]
        s["nrow"] += q["nrow"]; s["nfb"] += q["nfb"]
        s["splits"].append((q["split"], q["half"], q["P"]))
        s["trig_ok"] &= q["trig_equal"]
    t_passA = time.time() - t0
    for r in runs:
        s = RN[r]
        s["P"] = s["sp_sum"] / s["sp_n"] / RF.N_SUB if s["sp_n"] else RF.P_NOM
        s["splits"].sort()
        s["type"] = TYPE_OF_RUN[r]
        assert np.array_equal(s["H"][:, 1].sum((0, 2)), s["nrow"] - s["nfb"]), f"run {r}: H[allE] != row counts"
    print(f"pass A done in {t_passA:.0f} s; P per run: " + " ".join(f"{r}:{RN[r]['P']:.5f}" for r in runs)
          + f"; rows {sum(int(RN[r]['nrow'].sum()) for r in runs)}, fb rows {sum(int(RN[r]['nfb'].sum()) for r in runs)}",
          flush=True)
    if not all(RN[r]["trig_ok"] for r in runs):
        print("WARNING: trig differs between dump and rf for some matched events", flush=True)

    # ---- fits ----
    t0 = time.time()
    keys, fjobs = [], []
    for r in runs:
        k_, j_ = fit_run_jobs(r, RN[r])
        keys += k_; fjobs += j_
    fres = pool.map(fit_task, fjobs, chunksize=64)
    fits = {r: {} for r in runs}
    for (r, m, ti, slot), f in zip(keys, fres):
        fits[r][(ti, slot, m)] = f
    for r in runs:
        s = RN[r]
        s.update(prior=prior, tau=tau, board=board, gp=gp)
        estimate_run(s, prior, tau, board, fits[r])
        s["bimodal"] = s["meas"] & np.isfinite(s["ratio"]) & (s["ratio"] >= args.bimodal_ratio)
        s["frac_bimodal"] = float(s["bimodal"].sum() / max(int(s["meas"].sum()), 1))
        s["usable"] = s["frac_bimodal"] <= BIM_RUN_FRAC
        s["quality"] = "ok" if s["usable"] else "bimodal"
    t_fits = time.time() - t0
    n_fit_all = sum(RN[r]["n_fits"] for r in runs)
    n_arg = sum(RN[r]["n_argmax_fits"] for r in runs)
    print(f"fits: {n_fit_all} histogram fits in {t_fits:.0f} s, {n_arg} argmax fallbacks "
          f"({100 * n_arg / max(n_fit_all, 1):.2f} %)", flush=True)
    for r in runs:
        s = RN[r]
        tc = [int((s["tier"] == k).sum()) for k in range(4)]
        print(f"  {r} {s['type']:5} e100 {tc[0]:4d} allE {tc[1]:4d} pred {tc[2]:4d} unresolved {tc[3]:4d}  "
              f"Phi {s['Phi']:+.3f}  walk {s['walk']:+.3f}  bimodal {100 * s['frac_bimodal']:.1f} % -> "
              f"{'usable' if s['usable'] else 'UNUSABLE'}", flush=True)

    # ---- branch: pass A, pass B ----
    t0 = time.time()
    for r in runs:
        s = RN[r]
        s["e1"], s["sig_e1"] = efit_arrays(s["HE"], s, pool)
        b1 = branch(s, tau, prior, s["e1"])
        s.update(G1=b1["G"], G1_mad=b1["G_mad"], D1=b1["D"], margin1=b1["margin"], firm1=b1["firm"])
    t_e1 = time.time() - t0
    t0 = time.time()
    jobs = []
    for r in runs:
        s = RN[r]
        corr = np.zeros(len(tau_full))
        corr[pidx[s["firm1"]]] = s["D1"][s["firm1"]]
        ch = np.where(np.isfinite(s["c_half"][0]) & np.isfinite(s["c_half"][1]), s["c_half"], np.nan)
        for i, (sp, fp, rp) in enumerate(by_run[r]):
            jobs.append((fp, rp, cfg, i % 2, corr, np.where(s["meas"], s["o"], np.nan), ch[1 - (i % 2)]))
    for r in runs:
        RN[r].update(HE2=np.zeros((2, npw, EVT_NB), np.int64), hh=np.zeros(NB, np.int64), drift=[])
    for q in RF.run_pool(pool, pass_b, jobs, "passB"):
        s = RN[q["run"]]
        s["HE2"] += q["he"]
        s["hh"] += q["hh"]
        s["drift"].append((q["split"], q["P"], q["hd"]))
    t_passB = time.time() - t0
    for r in runs:
        s = RN[r]
        P = s["P"]
        s["e2"], s["sig_e2"] = efit_arrays(s["HE2"], s, pool)
        b2 = branch(s, tau, prior, s["e2"])
        s.update(G2=b2["G"], G2_mad=b2["G_mad"], D2=b2["D"], margin2=b2["margin"], firm2=b2["firm"],
                 nobranch=b2["nobranch"])
        s["branch_weak"] = b2["use"] & ~b2["firm"]
        both = s["firm1"] & np.isfinite(s["D2"])
        s["branch_changed"] = both & (np.abs(s["D2"] - s["D1"]) > P / 2)
        sp = sorted(s["drift"], key=lambda z: z[0])
        sl, cs = [], []
        for spn, Ps, hd in sp:
            if hd.sum() >= MIN_SPLIT_N:
                sl.append(spn); cs.append(float(RF.fold(RF.fit_unwrapped((hd, Ps))["mu_u"], Ps)))
        cs = np.array(cs)
        if cs.size:
            dv = fold(cs - RF.circ_median(cs, P), P)
            s["split_span"] = float(dv.max() - dv.min())
        else:
            s["split_span"] = np.nan
        s["drift_pts"] = (np.array(sl), cs)
        if s["usable"] and s["split_span"] > DRIFT_SPAN:
            s["quality"] = "drift"
        # pooled tcal fit (cross-check against rf_summary.tsv)
        s["pooled_tcal"] = float(fold(RF.fit_unwrapped((s["ht"], P))["mu_u"], P))
        print(f"  {r}: G1 {s['G1']:+.3f} G2 {s['G2']:+.3f} (G2-G1 {fold(s['G2'] - s['G1'], P):+.3f})  "
              f"G_mad {s['G1_mad']:.3f}/{s['G2_mad']:.3f}  split_span {s['split_span']:.3f}  quality {s['quality']}",
              flush=True)
    pool.close()
    pool.join()
    print(f"pass B / branch done: e1 fits {t_e1:.0f} s, pass B {t_passB:.0f} s", flush=True)

    # ---- outputs ----
    created = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
    xc = {}
    if xcheck.exists():
        for q in csv.DictReader(open(xcheck), delimiter="\t"):
            xc[int(q["run"])] = float(q["run_offset_tcal"])
    for r in runs:
        s = RN[r]
        s["flags"] = np.array([",".join(flags_of(s, m)) for m in range(npw)])
        n = write_json(out_dir / f"{r}.json", r, s, gp, created)
        print(f"wrote {r}.json: usable {s['usable']}, quality {s['quality']}, {n} modules", flush=True)
    write_modules_tsv(out_dir / "rf_calib_modules.tsv", runs, RN, gp, board, tau)
    ST = run_stats(runs, RN, xc, args.bimodal_ratio)
    write_runs_tsv(out_dir / "rf_calib_runs.tsv", runs, RN, ST, args.bimodal_ratio, xc)
    ok_all = write_validation(out_dir / "rf_calib_validation.txt", runs, RN, ST, gp, board, args, xc)
    print(open(out_dir / "rf_calib_validation.txt").read(), flush=True)

    # ---- PDFs ----
    t_pdf = 0.0
    if not args.no_pdf:
        t0 = time.time()
        byname = {str(n): i for i, n in enumerate(gp["name"])}
        pj = []
        for r in runs:
            s = RN[r]
            panel = []
            if r in X17_RUNS:
                panel = [byname[n] for n in CENTRAL + HALO]
                panel.sort(key=lambda m: gp["id"][m])
            fl = [m for m in range(npw) if m not in panel and any(
                f not in ("allE", "run_unusable") for f in s["flags"][m].split(",") if f)]
            panel += fl
            Hb = s["H"].sum(0)
            D = dict(P=s["P"], gp=gp, type=s["type"], h_raw=Hb[1].sum(0), h_tau=s["ht"], h_after=s["hh"],
                     rho_run=s["rho_run"], dhalf=s["dhalf"], tier=s["tier"], o=s["o"], dhalf_err=s["dhalf_err"],
                     e=s["e2"], D=s["D2"], margin=s["margin2"], flags=s["flags"], Hb=Hb, panel_mods=panel,
                     drift=(s["drift_pts"] if s["quality"] == "drift" else None), split_span=s["split_span"],
                     text=run_text(r, s, ST[r], len(panel)))
            pj.append((str(out_dir / f"rf_calib_{r}.pdf"), r, D))
        with ctx.Pool(3) as p3:
            res = [p3.apply_async(render_run, (j,)) for j in pj]
            for q in res:
                path, dt, npg = q.get()
                print(f"wrote {path} ({dt:.0f} s, {npg} pages)", flush=True)
        t_pdf = time.time() - t0
        try:
            import pypdf
            for path, r, _ in pj:
                pypdf.PdfReader(path)
            print(f"pypdf: all {len(pj)} PDFs open", flush=True)
        except ImportError:
            pass
    print(f"timing: pass A {t_passA:.0f} s, fits {t_fits:.0f} s, e1 {t_e1:.0f} s, pass B {t_passB:.0f} s, "
          f"PDFs {t_pdf:.0f} s, total wall time {time.time() - t_start:.0f} s", flush=True)
    if not ok_all:
        print("GATE FAILED (see rf_calib_validation.txt)", flush=True)
        sys.exit(1)


# ---- statistics, tables ------------------------------------------------------

def run_stats(runs, RN, xc, bim_ratio):
    ST = {}
    for r in runs:
        s = RN[r]
        P = s["P"]
        meas = s["meas"]
        dh = s["dhalf"][meas]
        dh = dh[np.isfinite(dh)]
        e100 = (s["tier"] == 0) & np.isfinite(s["sigma"])
        rr = s["rho_run"][np.isfinite(s["rho_run"])]
        ST[r] = dict(
            n_meas=int(meas.sum()), n_dh=int(dh.size),
            frac_dh=float(np.mean(np.abs(dh) < 0.30)) if dh.size else np.nan,
            med_dh=float(np.median(np.abs(dh))) if dh.size else np.nan,
            frac_rho=float(np.mean(np.abs(rr) < 0.30)) if rr.size else np.nan,
            med_sig=float(np.median(s["sigma"][e100])) if e100.any() else np.nan,
            frac_sig=float(np.mean(s["sigma"][e100] < 0.30)) if e100.any() else np.nan,
            n_tier=[int((s["tier"] == k).sum()) for k in range(4)],
            pool_diff=abs(s["pooled_tcal"] - xc[r]) if r in xc else np.nan)
    return ST


def write_modules_tsv(path, runs, RN, gp, board, tau):
    cols = ("run type module id crate slot n_e100 n_allE tier offset_ns offset_err_ns sigma_ns tau_ns rho_ns "
            "rho_run_ns board_phase_ns board_mad_ns evt_ns dt_tcal_ns margin_ns dhalf_ns dhalf_err_ns flags").split()
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in runs:
            s = RN[r]
            for m in np.argsort(gp["id"], kind="stable"):
                t = int(s["tier"][m])
                meas = t < 2
                row = [r, s["type"], str(gp["name"][m]), int(gp["id"][m]), int(board[m] // 100), int(board[m] % 100),
                       int(s["n_e100"][m]), int(s["n_allE"][m]), TIERS[t],
                       s["o"][m] if t < 3 else np.nan, s["err"][m] if meas else np.nan,
                       s["sigma"][m] if meas else np.nan, tau[m], s["rho"][m], s["rho_run"][m], s["bphase"][m],
                       s["bmad"][m], s["e2"][m], s["D2"][m] if t < 2 else np.nan, s["margin2"][m],
                       s["dhalf"][m] if meas else np.nan, s["dhalf_err"][m] if meas else np.nan,
                       ",".join(flags_of(s, m))]
                fh.write("\t".join(tsv_fmt(v) for v in row) + "\n")


def write_runs_tsv(path, runs, RN, ST, bim_ratio, xc):
    cols = ("run type n_splits P_ns usable quality n_e100 n_allE n_pred n_unresolved walk_bias_ns split_span_ns "
            "frac_bimodal bimodal_ratio run_phase_ns G1_ns G2_ns G_mad_ns pooled_tcal_offset_ns "
            "rf_summary_run_offset_tcal_ns frac_dhalf_lt_0p3 median_abs_dhalf frac_rhorun_lt_0p3 median_sigma_e100 "
            "frac_sigma_lt_0p3").split()
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in runs:
            s, q = RN[r], ST[r]
            row = [r, s["type"], len(s["splits"]), s["P"], s["usable"], s["quality"], *q["n_tier"][:4],
                   s["walk"], s["split_span"], s["frac_bimodal"], bim_ratio, s["Phi"], s["G1"], s["G2"],
                   s["G2_mad"], s["pooled_tcal"], xc.get(r, np.nan), q["frac_dh"], q["med_dh"], q["frac_rho"],
                   q["med_sig"], q["frac_sig"]]
            fh.write("\t".join(tsv_fmt(v, 6) for v in row) + "\n")


def run_text(r, s, q, npanel):
    f = lambda v, nd=3: "nan" if v is None or not np.isfinite(v) else f"{v:.{nd}f}"
    return [f"run {r}  type {s['type']}  P = {s['P']:.5f} ns  splits {len(s['splits'])}",
            f"quality {s['quality']}  usable {s['usable']}  split span {f(s['split_span'])} ns",
            f"tiers: e100 {q['n_tier'][0]}, allE {q['n_tier'][1]}, pred {q['n_tier'][2]}, unresolved {q['n_tier'][3]}",
            f"walk bias (median c_allE - c_e100) {f(s['walk'])} ns;  run phase Phi_r {f(s['Phi'])} ns",
            f"bimodal modules {100 * s['frac_bimodal']:.1f} % of the measured ones",
            f"G1 {f(s['G1'])}  G2 {f(s['G2'])}  G_mad {f(s['G1_mad'])} / {f(s['G2_mad'])} ns",
            f"dhalf: median |dhalf| {f(q['med_dh'])} ns, fraction < 0.3 ns {f(q['frac_dh'])} ({q['n_dh']} modules)",
            f"|rho^r| < 0.3 ns for {f(q['frac_rho'])} of the measured modules",
            f"sigma (tier e100): median {f(q['med_sig'])} ns, fraction < 0.3 ns {f(q['frac_sig'])}",
            f"pooled tcal offset {f(s['pooled_tcal'], 4)} ns  (rf_summary: diff {f(q['pool_diff'], 4)})",
            f"module panels: {npanel} ({'33 X17 modules + ' if r in X17_RUNS else ''}flagged modules)"]


def write_validation(path, runs, RN, ST, gp, board, args, xc):
    byname = {str(n): i for i, n in enumerate(gp["name"])}
    f = lambda v, nd=3: "nan" if v is None or not np.isfinite(v) else f"{v:.{nd}f}"
    L = []
    use = [r for r in runs if RN[r]["usable"]]
    L.append(f"rf_calib_validation.txt  runs {runs}  usable {use}  bimodal ratio {args.bimodal_ratio}")
    L.append(f"created {datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')}")
    # (i)
    L.append("\n(i) per-run sigma of the tier e100 modules (report; expected median 0.18-0.30 ns)")
    L.append("  run   type   n_e100  median sigma  frac sigma<0.3")
    for r in runs:
        q = ST[r]
        L.append(f"  {r}  {RN[r]['type']:5}  {q['n_tier'][0]:6d}  {f(q['med_sig']):>12}  {f(q['frac_sig']):>14}")
    for t in sorted({RN[r]["type"] for r in runs}):
        v = np.concatenate([RN[r]["sigma"][RN[r]["tier"] == 0] for r in use if RN[r]["type"] == t] or [np.zeros(0)])
        v = v[np.isfinite(v)]
        L.append(f"  type {t}: median sigma {f(np.median(v) if v.size else np.nan)} ns over {v.size} (run, module) fits")
    # (ii)
    ok = {}
    L.append("\n(ii) GATE held-out dhalf over the usable runs: median |dhalf| < 0.05 ns and >= 99 % with |dhalf| < 0.30 ns")
    allv = []
    for r in use:
        s = RN[r]
        v = s["dhalf"][s["meas"]]; v = v[np.isfinite(v)]
        allv.append(v)
        L.append(f"  {r}{' (drift)' if s['quality'] == 'drift' else ''}: n {v.size:5d}/{int(s['meas'].sum())}  "
                 f"median |dhalf| {f(np.median(np.abs(v)) if v.size else np.nan)}  "
                 f"frac < 0.30 {f(np.mean(np.abs(v) < 0.3) if v.size else np.nan)}")
    v = np.concatenate(allv) if allv else np.zeros(0)
    n_meas_all = int(sum(RN[r]["meas"].sum() for r in use))
    if v.size:
        md, fr = float(np.median(np.abs(v))), float(np.mean(np.abs(v) < 0.3))
        ok["ii"] = md < 0.05 and fr >= 0.99
        L.append(f"  all usable runs: n {v.size} of {n_meas_all} measured (dhalf defined when both halves >= {N_HALF}), "
                 f"median |dhalf| {md:.4f}, fraction < 0.30 ns {fr:.4f}  ->  {'PASS' if ok['ii'] else 'FAIL'}")
        L.append(f"  (over all {n_meas_all} measured, counting undefined dhalf as failing: "
                 f"{float(np.sum(np.abs(v) < 0.3)) / n_meas_all:.4f}; not gated)")
        if 24503 in use:
            s = RN[24503]; v3 = s["dhalf"][s["meas"]]; v3 = v3[np.isfinite(v3)]
            L.append(f"  24503 alone: median |dhalf| {f(np.median(np.abs(v3)), 4)}, fraction < 0.30 {f(np.mean(np.abs(v3) < 0.3), 4)}, "
                     f"max |dhalf| {f(np.max(np.abs(v3)))}")
    else:
        ok["ii"] = True
        L.append("  no usable run: not evaluated")
    # (iii)
    L.append("\n(iii) prior consistency: fraction of measured modules with |rho^r| < 0.3 ns (report; expected >~ 95 % "
             "at the H2 runs, 25800, 26138)")
    L.append("  run   type   n_meas  Phi_r   frac |rho^r|<0.3   n(|rho^r|>=0.3)")
    for r in runs:
        s, q = RN[r], ST[r]
        L.append(f"  {r}  {s['type']:5}  {q['n_meas']:6d}  {f(s['Phi']):>6}  {f(q['frac_rho']):>16}   "
                 f"{int(np.sum(np.abs(s['rho_run'][np.isfinite(s['rho_run'])]) >= 0.3))}")
    # (iv)
    L.append("\n(iv) board pattern 25800 vs 26230 (report; expected global jump ~ +2.0 ns, every board within +-0.1 ns "
             "except (3,15), (3,8))")
    if 25800 in RN and 26230 in RN:
        a, b = RN[25800], RN[26230]
        P = a["P"]
        both = np.isfinite(a["rho"]) & np.isfinite(b["rho"])
        g = RF.circ_median(fold(a["rho"][both] - b["rho"][both], P), P)
        L.append(f"  modules measured in both: {int(both.sum())}; global jump (circ median of rho_25800 - rho_26230) {g:+.3f} ns")
        L.append("  board  n_both  Delta_b(25800)  Delta_b(26230)  jump - global  mad_b(25800)  mad_b(26230)")
        rows = []
        for bb in np.unique(board):
            s_ = both & (board == bb)
            if s_.sum() < MIN_BOARD:
                continue
            da = RF.circ_median(a["rho"][s_], P); db = RF.circ_median(b["rho"][s_], P)
            rows.append((abs(fold(da - db - g, P)), bb, int(s_.sum()), da, db, fold(da - db - g, P),
                         RF.circ_mad(a["rho"][s_], P), RF.circ_mad(b["rho"][s_], P)))
        rows.sort(key=lambda z: -z[0])
        for _, bb, n, da, db, dj, ma, mb in rows:
            L.append(f"  ({bb // 100},{bb % 100:2d})  {n:5d}  {da:+14.3f}  {db:+14.3f}  {dj:+13.3f}  {ma:12.3f}  {mb:12.3f}")
        out = [(int(bb // 100), int(bb % 100), round(float(dj), 3)) for _, bb, n, da, db, dj, ma, mb in rows if abs(dj) > 0.1]
        L.append(f"  boards with >= {MIN_BOARD} modules: {len(rows)}; outside +-0.1 ns of the global jump: {out}")
    else:
        L.append("  25800 or 26230 not in this run")
    # (v)
    L.append("\n(v) GATE the 33 X17 central + halo modules (26 + 7): in 25800, 26138, 26230 all measured, firm "
             "(margin >= 1 ns) and |D2 - e2| < 0.5 ns; no pass-A-firm module with |D2 - D1| > P/2")
    mods33 = CENTRAL + HALO
    ok["v"] = True
    show = [r for r in runs if r in X17_RUNS or r == 25320]
    for r in show:
        s = RN[r]
        L.append(f"\n  run {r} ({s['type']}): module board tier  rho^r      e2      D2   margin  dhalf  firm  |D2-e2|")
        nfirm, mx, nm = 0, 0.0, 0
        for n in mods33:
            m = byname[n]
            t = int(s["tier"][m])
            measd = t < 2 and np.isfinite(s["e2"][m])
            firm = bool(s["firm2"][m])
            dd = abs(s["D2"][m] - s["e2"][m]) if measd else np.nan
            nfirm += firm; nm += measd
            if measd:
                mx = max(mx, dd)
            if r in X17_RUNS and not (measd and firm and dd < 0.5):
                ok["v"] = False
            L.append(f"    {n:5} ({board[m] // 100},{board[m] % 100:2d}) {TIERS[t]:5} {f(s['rho_run'][m]):>7} "
                     f"{f(s['e2'][m]):>7} {f(s['D2'][m]):>7} {f(s['margin2'][m], 2):>7} {f(s['dhalf'][m]):>6} "
                     f"{'yes' if firm else 'NO':>4}  {f(dd)}")
        L.append(f"  -> run {r}: measured {nm}/33, firm {nfirm}/33, max |D2 - e2| {mx:.3f} ns")
    nch = 0
    for r in use:
        nch += int(RN[r]["branch_changed"].sum())
    if nch:
        ok["v"] = False
    L.append(f"\n  branch_changed among pass-A-firm modules (usable runs): {nch}")
    L.append("  expectations (report, deviations > 0.3 ns beyond the bound listed):")
    for r in runs:
        s = RN[r]
        if not s["usable"]:
            continue
        idx = [byname[n] for n in mods33]
        D = s["D2"][idx]
        if r == 26230:
            lo, hi = -2.6, -0.4
            dev = [(n, float(d)) for n, d in zip(mods33, D) if np.isfinite(d) and not lo - 0.3 <= d <= hi + 0.3]
            L.append(f"    {r}: D2 in [{np.nanmin(D):+.2f}, {np.nanmax(D):+.2f}] (expected [-2.6, -0.4]); beyond +-0.3: {dev}")
        elif r in X17_RUNS:
            quiet = [(n, float(d)) for n, d in zip(mods33, D) if n not in FOUR]
            nq = sum(abs(d) <= 0.4 for _, d in quiet)
            dev = [(n, round(d, 2)) for n, d in quiet if abs(d) > 0.7]
            four = {n: round(float(s["D2"][byname[n]]), 2) for n in FOUR}
            fdev = [n for n, d in four.items() if not -1.5 - 0.3 <= d <= -0.6 + 0.3]
            L.append(f"    {r}: |D2| <= 0.4 for {nq}/29 non-W525/528/597/697 modules (expected 29); beyond 0.7: {dev}; "
                     f"four {four} (expected [-1.5, -0.6]); beyond +-0.3: {fdev}")
        elif s["type"] == "2p24":
            dev = [(n, round(float(d), 2)) for n, d in zip(mods33, D) if np.isfinite(d) and abs(d) > 0.7]
            L.append(f"    {r}: |D2| <= 0.4 for {int(np.sum(np.abs(D[np.isfinite(D)]) <= 0.4))}/33 (expected 33); beyond 0.7: {dev}")
    L.append(f"  (v) -> {'PASS' if ok['v'] else 'FAIL'}")
    # (vi)
    L.append("\n(vi) GATE G_mad = circ_mad(fold(o - tau - e)) < 0.15 ns per usable run (pass A and B), and the pooled "
             "tcal fit reproduces rf_summary.tsv:run_offset_tcal to < 0.002 ns")
    ok["vi"] = True
    L.append("  run   usable  G_mad(A)  G_mad(B)  G1      G2     G2-G1   pooled tcal   rf_summary   |diff|")
    mxd = 0.0
    for r in runs:
        s, q = RN[r], ST[r]
        P = s["P"]
        if s["usable"] and not (s["G1_mad"] < 0.15 and s["G2_mad"] < 0.15):
            ok["vi"] = False
        if r in xc:
            mxd = max(mxd, q["pool_diff"])
            if not q["pool_diff"] < 0.002:
                ok["vi"] = False
        L.append(f"  {r}  {str(s['usable']):6}  {s['G1_mad']:8.3f}  {s['G2_mad']:8.3f}  {s['G1']:+6.3f}  {s['G2']:+6.3f}  "
                 f"{fold(s['G2'] - s['G1'], P):+6.3f}  {s['pooled_tcal']:+11.5f}  {f(xc.get(r, np.nan), 5):>10}  "
                 f"{f(q['pool_diff'], 5)}")
    L.append(f"  max |pooled - rf_summary| {mxd:.5f} ns over {sum(r in xc for r in runs)} runs;  (vi) -> "
             f"{'PASS' if ok['vi'] else 'FAIL'}")
    L.append("\nGATES: " + "  ".join(f"({k}) {'PASS' if v else 'FAIL'}" for k, v in ok.items()))
    with open(path, "w") as fh:
        fh.write("\n".join(L) + "\n")
    return all(ok.values())


if __name__ == "__main__":
    main()
