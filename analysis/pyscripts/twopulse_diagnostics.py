#!/usr/bin/env python3
"""twopulse_diagnostics.py -- is the U0 two-pulse chi2 fit over-fitting?

Diagnostics for twopulse_fit.py (U0, production settings eps 0.02, kap 0, search 80,
min_sep 1, base_t_min 150).  Imports twopulse_fit; changes nothing in it.

All numbers are on the evaluation half of the real-pulse pool (odd pool index);
thresholds come from the calibration half (even index) only.  Primary doublet arm
`real2` (two real clean pulses of one module and half); `real` (noise-free template
injected) is optimistic and shown only where cheap.  Pages 1-6 (2.24 GeV H2), one page
for 0.73 / 3.49 GeV H2, independence, and the small-ratio (r 0.1, 0.05) study.

    MPLBACKEND=Agg python3 twopulse_diagnostics.py
    MPLBACKEND=Agg python3 twopulse_diagnostics.py --quick        (smoke test)
    MPLBACKEND=Agg python3 twopulse_diagnostics.py --replot       (from the pickle cache)

Outputs: output/ml_pileup_20261027/U0/overfit_diagnostics.pdf and .md
Scratch caches (--cache): pool-fit npz, state pickle.
"""
from __future__ import annotations

import argparse
import contextlib
import copy
import glob
import json
import math
import multiprocessing as mp
import pickle
import sys
import time
from pathlib import Path

import numpy as np

import twopulse_fit as tf
from twopulse_fit import DT, NS, TS

HERE = Path(__file__).resolve().parent
OUT = HERE / "output/ml_pileup_20261027/U0"
SCRATCH = Path("/private/tmp/claude-505/-Users-L00338853-work-PRad-prad2evviewer-analysis-pyscripts/"
               "b1de5a72-6ab3-4e2b-8b9b-77e955030d7b/scratchpad")
SEED = 20261027
MAIN = "2.24 GeV H2"
OTHERS = ("0.73 GeV H2", "3.49 GeV H2")
ALL_RT = (MAIN,) + OTHERS
OPT_PROD = dict(eps=0.02, kap=0.0, search=80.0, grid=0.5, min_sep=1.0)
BASE_T_MIN = 150.0
F_PROD, F_01, F_002 = 0.03, 1e-3, 2e-4
P2_DT, P2_R = (4, 8, 12, 20), (1.0, 0.3)
P3_DT = (4, 6, 8, 10, 12, 16, 20, 24)
ROC_DT, ROC_R = (8, 12), (1.0, 0.5, 0.3, 0.1, 0.05)
SM_DT, SM_R = (4, 8, 12, 16, 20, 24, 32), (0.1, 0.05)
KAP_DT = (-32, -16, -8, 8, 16, 32, 64)     # t2 - t1 offsets (ns) of the extra pulse in the analytic null
SMALL_LIMIT = 0.2          # ratios below this use the local small-pulse event table
COL = dict(single="#0072B2", tpl="#009E73", ana="#7f7f7f", dbl="#D55E00", real="#CC79A7",
           thr="k", p1="#0072B2", p2="#D55E00")


def log(*a):
    print(*a, flush=True)


# =================================================================== helpers ==
@contextlib.contextmanager
def forked(bank, workers):
    """Fork a worker pool whose fit bank is `bank` (twopulse_fit._worker reads _BANK at fork)."""
    prev = tf._BANK
    tf._BANK = bank
    pw = mp.get_context("fork").Pool(workers)
    try:
        yield pw
    finally:
        pw.close()
        pw.join()
        tf._BANK = prev


def shifted_bank(bank, s, ramp_ns=0.0):
    """Copy of `bank` whose templates are delayed by s ns (tpl'(x) = tpl(x - s)).

    s < 0 cuts the head of the template, which leaves a step (about 1.4 % of the peak for 0.5 ns)
    at the new onset; ramp_ns > 0 removes that step linearly over ramp_ns ns (value and slope
    corrected), so the shifted template starts smoothly at its onset value.
    """
    k = int(round(s / tf.LUT_DT))
    b = copy.copy(bank)

    def sh(a):
        a = np.asarray(a)
        out = np.zeros_like(a)
        if k > 0:
            out[:, k:] = a[:, :-k]
        elif k < 0:
            out[:, :k] = a[:, -k:]
            out[:, k:] = a[:, -1:]
        else:
            out = a.copy()
        return out
    b.lut, b.dlut = sh(bank.lut), sh(bank.dlut)
    b.t_pk = np.asarray(bank.t_pk) + s
    if ramp_ns > 0 and k < 0:
        c0 = np.asarray(bank.lut)[:, -k] - np.asarray(bank.lut)[:, 0]     # step at the cut head
        nr = int(round(ramp_ns / tf.LUT_DT))
        w = np.clip(1.0 - np.arange(b.lut.shape[1]) / nr, 0.0, 1.0)
        b.lut = (b.lut - c0[:, None] * w[None, :]).astype(b.lut.dtype)
        dw = np.where(np.arange(b.lut.shape[1]) < nr, 1.0 / ramp_ns, 0.0)
        b.dlut = (b.dlut + c0[:, None] * dw[None, :]).astype(b.dlut.dtype)
    return b


def wil(k, n):
    lo, hi = tf.wilson(int(k), int(n))
    return lo, hi


def pcs(x, d=1):
    return f"{100 * x:.{d}f}"


def eff_str(k, n, d=0):
    lo, hi = wil(k, n)
    return f"{100 * k / max(n, 1):.{d}f} [{100 * lo:.{d}f}-{100 * hi:.{d}f}]"


def match_t(R, T1, T2, tol=4.0):
    """Two-pulse fit pulses match truth onsets (T1 earlier, T2 later) within tol (as pair_match)."""
    return (R[:, 13] == 1) & (np.abs(R[:, 9] - T1) <= tol) & (np.abs(R[:, 11] - T2) <= tol)


def match_amp(R, Tb, Ab, Ts, As, tol=4.0, fac=2.0):
    """Two-pulse fit has a pulse on the big truth pulse (onset within tol) AND a second pulse on
    the small truth pulse (onset within tol, amplitude within a factor fac of the truth)."""
    two = R[:, 13] == 1
    t1, t2, a1, a2 = R[:, 9], R[:, 11], R[:, 8], R[:, 10]
    ok_s = lambda a: (a >= As / fac) & (a <= As * fac)
    x = (np.abs(t1 - Tb) <= tol) & (np.abs(t2 - Ts) <= tol) & ok_s(a2)
    y = (np.abs(t2 - Tb) <= tol) & (np.abs(t1 - Ts) <= tol) & ok_s(a1)
    return two & (x | y)


def overlap_coef(a, b, lo=-1.0, hi=7.0, nb=80):
    """Overlap coefficient sum_bins min(p_a, p_b), log10(dchi2) bins (values clipped at 0.1)."""
    bins = np.linspace(lo, hi, nb + 1)
    la = np.log10(np.clip(a, 0.1, None))
    lb = np.log10(np.clip(b, 0.1, None))
    pa = np.histogram(la, bins)[0] / len(a)
    pb = np.histogram(lb, bins)[0] / len(b)
    return float(np.minimum(pa, pb).sum())


def edge_flags(R, y, bank, m):
    """Per-event: fitted t1 on a search-window edge, or separation on the min_sep / upper edge."""
    o = tf.OPT
    g, srch = o["grid"], o["search"]
    c = DT * np.argmax(y, 1) - np.asarray(bank.t_pk)[m]
    c0 = g * np.round((c - srch) / g)
    K = int(round(2 * srch / g)) + 1
    lo, hi = c0 - 2.0, c0 + g * (K - 1) + 2.0
    sep = R[:, 11] - R[:, 9]
    two = R[:, 13] == 1
    e_t = two & ((R[:, 9] < lo + 0.05) | (R[:, 9] > hi - 0.05))
    e_s = two & ((sep < o["min_sep"] + 0.05) | (sep > 2 * srch + 4 - 0.05))
    return e_t | e_s


# ============================================================ small-ratio table ==
def event_table_small(pool, n, dt_nom, r, jit, rng, partner, half=1):
    """Doublets with a SMALL pulse of amplitude r * A_big (r < 1), local to this script.

    The base real pulse is always the big pulse; the small one is injected 50/50 before or
    after it at dt_nom + U(-jit, jit) ns.  Both orders are therefore treated identically and
    the base amplitude is not biased by the A1 + A2 <= 4000 clip except for A_b > 4000/(1+r).
    partner (arm real2): the small pulse is the real pulse of the same module and half whose
    fitted amplitude is nearest (one of the 3 nearest, no factor-2 limit) to r * A_b, rescaled
    to that amplitude (its noise is scaled by the same factor, i.e. some extra noise remains)
    and moved by a whole number of samples (its own sub-sample phase gives the phase jitter).
    Returns the production event-table dict plus Ab, tb, As, ts and the clip-reject count.
    """
    cand = np.nonzero((pool.parity == half) & (pool.has_partner if partner else True))[0]
    cols, got, nrej, ndr = [], 0, 0, 0
    while got < n:
        k = n - got
        idx = rng.choice(cand, k)
        small_first = rng.random(k) < 0.5
        dtt = np.maximum(dt_nom + rng.uniform(-jit, jit, k), 0.5) if jit > 0 else np.full(k, float(dt_nom))
        Ab, tb = pool.A_b[idx], pool.t_b[idx]
        Ai = r * Ab
        ti = np.where(small_first, tb - dtt, tb + dtt)
        ok = (Ab + Ai) <= 4000.0
        ndr += k
        nrej += int((~ok).sum())
        part = np.full(k, -1, np.int64)
        shift = np.zeros(k, np.int64)
        if partner:
            for j in range(k):
                if not ok[j]:
                    continue
                pm = pool.partners[(half, int(pool.midx[idx[j]]))]
                pm = pm[pm != idx[j]]
                lg = np.abs(np.log(pool.A_b[pm] / Ai[j]))
                part[j] = rng.choice(pm[np.argsort(lg)[:3]])
                shift[j] = int(np.rint((ti[j] - pool.t_b[part[j]]) / DT))
                ti[j] = pool.t_b[part[j]] + DT * shift[j]
            dtt = np.where(small_first, tb - ti, ti - tb)
        cols.append(np.stack([idx, small_first, dtt, Ai, ti, part, shift], 1)[ok])
        got += int(ok.sum())
    rr = np.concatenate(cols)[:n]
    idx = rr[:, 0].astype(np.int64)
    sf = rr[:, 1] > 0.5
    dtt, Ai, ti = rr[:, 2], rr[:, 3], rr[:, 4]
    Ab, tb = pool.A_b[idx], pool.t_b[idx]
    ev = dict(idx=idx, inj_first=sf, dt_true=dtt, A_inj=Ai, t_inj=ti,
              A1=np.where(sf, Ai, Ab), t1=np.where(sf, ti, tb),
              A2=np.where(sf, Ab, Ai), t2=np.where(sf, tb, ti),
              Ab=Ab, tb=tb, As=Ai, ts=ti, n_reject=nrej, n_drawn=ndr)
    if partner:
        ev.update(partner=rr[:, 5].astype(np.int64), shift=rr[:, 6].astype(np.int64))
    return ev


# ================================================================== context ==
class Ctx:
    pass


def build_ctx(rt, workers, cache, quick):
    C = Ctx()
    C.rt = rt
    C.irt = ALL_RT.index(rt)
    slug, split = tf.RUNS[rt]
    C.slug, C.split = slug, split
    C.bank, _ = tf.load_bank(tf.B50, slug)
    C.xbank, _ = tf.load_bank(tf.WAS, slug)
    tf._BANK = C.bank
    pool = tf.Pool(slug, split, C.bank, C.xbank, base_t_min=BASE_T_MIN)
    C.pool = pool
    C.pw = mp.get_context("fork").Pool(workers)
    t0 = time.time()
    y, ped = tf.subtract_ped(C.bank, pool.raw)
    cf = cache / f"poolfit_{slug}.npz"
    R = None
    if cf.exists():
        z = np.load(cf)
        if z["R"].shape[0] == len(pool.raw) and z["opt"].tolist() == [OPT_PROD[k] for k in sorted(OPT_PROD)]:
            R = z["R"]
    if R is None:
        R = tf.fit_batch(C.pw, y, pool.midx, pool.sigp)
        np.savez(cf, R=R, opt=np.array([OPT_PROD[k] for k in sorted(OPT_PROD)]))
    keep = (R[:, 3] >= BASE_T_MIN) & (R[:, 2] >= 100) & (R[:, 2] <= 4000)
    tf.slice_pool(pool, R, ped, keep)
    C.cal = pool.parity == 0
    C.ev = pool.parity == 1
    C.Rb = pool.Rb
    log(f"[{rt}] pool {pool.n} bases (cal {int(C.cal.sum())} / eval {int(C.ev.sum())}); {time.time() - t0:.0f} s")
    rng = np.random.default_rng([SEED, C.irt, 1])
    W = tf.make_waves("tpl_noise", pool, C.bank, C.xbank, np.arange(pool.n), rng)
    C.Rt = tf.fit_waves(C.pw, C.bank, W, np.arange(pool.n), pool)
    # the stored WaveAnalyzer rms (sig_ped) is ~2.3x smaller than the sample std of the pulse baseline; pool.noise() follows sig_ped,
    # so tpl_noise has that smaller noise.  Supplementary null with noise matched to the measured baseline std:
    x30 = pool.raw[:, :tf.PED_N]
    ps = np.sqrt(((x30 - x30.mean(1, keepdims=True)) ** 2).sum(1) / (tf.PED_N - 1))
    C.noise_f = float(np.median(ps / pool.sigp))
    Wx = tf.make_waves("nf_same", pool, C.bank, C.xbank, np.arange(pool.n), rng)
    vx = np.sqrt((C.noise_f ** 2) * np.maximum(pool.sigp ** 2 - 1.0 / 12.0, 1e-4) + 1.0 / 12.0)
    C.Rtx = tf.fit_waves(C.pw, C.bank, np.rint(Wx + pool.noise(rng, vx, pool.n)), np.arange(pool.n), pool)
    C.thr = {f: tf.thr_for(C.Rb[C.cal, 12], f) for f in (F_PROD, F_01, F_002)}
    C.thr_tpl = {f: tf.thr_for(C.Rt[C.cal, 12], f) for f in (F_PROD, F_01, F_002)}
    log(f"[{rt}] thr real (0.03, 1e-3, 2e-4) = {[round(C.thr[f], 1) for f in C.thr]}, "
        f"tpl_noise {[round(C.thr_tpl[f], 1) for f in C.thr_tpl]}; {time.time() - t0:.0f} s")
    return C


def close_ctx(C):
    C.pw.close()
    C.pw.join()


def singles_summary(C, rng):
    """Eval-singles quantities for page 2 / 4 / (e): kappa-scaled analytic null, residual scale."""
    pool, bank = C.pool, C.bank
    ie = np.nonzero(C.ev)[0]
    y = pool.raw[ie] - pool.ped[ie][:, None]
    R = C.Rb[ie]
    kap = np.empty(len(ie))
    Chat = pool.noise_L @ pool.noise_L.T          # unit-variance correlated-noise matrix (pedestal acf)
    rms_body = np.empty(len(ie))
    max_abs = np.empty(len(ie))
    o = tf.OPT
    for j, i in enumerate(ie):
        row = bank.lut[pool.midx[i]].astype(np.float64)
        mod = R[j, 2] * tf.tpl(row, None, TS - R[j, 3], bank.T_end)
        w = tf.weights(y[j], pool.sigp[i], o["eps"], o["kap"])
        sg2 = max(pool.sigp[i] ** 2 - 1.0 / 12.0, 1e-4)        # true noise variance (as pool.noise)
        den = max(float(np.sum(w * mod * mod)), 1e-12)
        v = w * mod
        # E[gain] of one extra weighted template (A2 free, t2 fixed): the extra shape orthogonalised (weighted) against
        # the single-pulse columns {s1, ds1}; noise = correlated (acf, variance sg2) + white rounding noise 1/12.
        # Mean over t2 - t1 in KAP_DT.
        sw_ = np.sqrt(w)
        s1, d1 = tf.tpl(row, row, TS - R[j, 3], bank.T_end, True)
        Gm = np.stack([s1, d1], 1) * sw_[:, None]
        Sg = sg2 * Chat + np.eye(NS) / 12.0
        ks = []
        for dd in KAP_DT:
            s2 = tf.tpl(row, None, TS - R[j, 3] - dd, bank.T_end) * sw_
            if s2 @ s2 < 1e-12:
                continue
            q, _ = np.linalg.qr(Gm)
            r2 = s2 - q @ (q.T @ s2)
            nr = float(r2 @ r2)
            if nr < 1e-12:
                continue
            u = r2 * sw_                                   # (W^1/2 r2) -> W r2 = sw_ * (sw_ r2)
            ks.append(float(u @ Sg @ u) / nr)
        kap[j] = float(np.mean(ks)) if ks else np.nan
        res = y[j] - mod
        body = mod > 0.05 * R[j, 2]
        rms_body[j] = math.sqrt(float(np.mean(res[body] ** 2))) if body.any() else np.nan
        max_abs[j] = float(np.max(np.abs(res)))
    null = (np.nan_to_num(kap, nan=float(np.nanmedian(kap)))[:, None] * rng.chisquare(2, (len(ie), 5))).ravel()
    A = pool.A_b[ie]
    return dict(kappa=kap, null=null, rms_body=rms_body, max_abs=max_abs, A=A, sig=pool.sigp[ie],
                A_med=float(np.median(pool.A_b)), sig_med=float(np.median(pool.sigp)),
                rms_rel_med=float(np.nanmedian(rms_body / A)), max_rel_med=float(np.median(max_abs / A)),
                ratio_rms_med=float(np.nanmedian(0.05 * A / rms_body)),
                kappa_med=float(np.nanmedian(kap)),
                ratio_max_med=float(np.median(0.05 * A / max_abs)),
                ratio_sig_med=float(np.median(0.05 * A / pool.sigp[ie])))


# ===================================================================== cells ==
def snr2_small(C, W, y, ev):
    """Optimal-filter SNR^2 upper bound of the smaller pulse alone (ignores the projection onto
    the big pulse): A_s^2 sum_i w_i tpl_i^2, w from the fit sigma model on the measured y."""
    pool, bank = C.pool, C.bank
    out = np.empty(len(y))
    As = np.minimum(ev["A1"], ev["A2"])
    ts = np.where(ev["A1"] < ev["A2"], ev["t1"], ev["t2"])
    for j in range(len(y)):
        m = pool.midx[ev["idx"][j]]
        w = tf.weights(y[j], pool.sigp[ev["idx"][j]], tf.OPT["eps"], tf.OPT["kap"])
        s = tf.tpl(bank.lut[m].astype(np.float64), None, TS - ts[j], bank.T_end)
        out[j] = As[j] ** 2 * float(np.sum(w * s * s))
    return out


def run_cell(C, arm, dt, r, n, rng, mf=False, keepW=False):
    pool = C.pool
    small = r < SMALL_LIMIT
    if small:
        ev = event_table_small(pool, n, dt, r, 2.0, rng, partner=(arm == "real2"))
    else:
        ev = tf.event_table(pool, n, dt, r, 2.0, rng, partner=(arm == "real2"))
    idx = ev["idx"]
    W = tf.make_waves(arm, pool, C.bank, C.xbank, idx, rng, (ev["A_inj"], ev["t_inj"]), ev)
    y, _ = tf.subtract_ped(C.bank, W)
    R = tf.fit_batch(C.pw, y, pool.midx[idx], pool.sigp[idx])
    cell = dict(arm=arm, dt=dt, r=r, small=small, ev=ev, R=R, edge=edge_flags(R, y, C.bank, pool.midx[idx]),
                snr2=snr2_small(C, W, y, ev))
    if mf:
        cell["mf"] = tf.run_mf(C.bank, W, pool.midx[idx])
    if keepW:
        cell["W"] = W
    return cell


def truth_small(ev):
    """Big / small truth (onset, amplitude); for production tables the larger / smaller of the pair."""
    big1 = ev["A1"] >= ev["A2"]
    return (np.where(big1, ev["t1"], ev["t2"]), np.where(big1, ev["A1"], ev["A2"]),
            np.where(big1, ev["t2"], ev["t1"]), np.where(big1, ev["A2"], ev["A1"]))


def mf_eff(cell):
    """MF efficiency (pair_match +-4 ns) of a cell with cell['mf']."""
    tl, al = cell["mf"]
    ev = cell["ev"]
    ok, _ = tf.pair_match(tl, al, ev["t1"], ev["t2"], 4.0)
    return ok


def cell_keys(main=True):
    ks = []
    if main:
        ks += [("real2", dt, 1.0) for dt in P3_DT]
        ks += [("real2", dt, 0.3) for dt in P3_DT]
        ks += [("real2", dt, 0.5) for dt in ROC_DT]
        ks += [("real", dt, r) for dt in ROC_DT for r in (1.0, 0.5, 0.3)]
        ks += [(a, dt, r) for a in ("real2", "real") for dt in SM_DT for r in SM_R]
    return ks


def compute_main(C, n, rng_seed, quick):
    S = dict(rt=C.rt)
    pool = C.pool
    ev_idx = np.nonzero(C.ev)[0]
    rng = np.random.default_rng([SEED, C.irt, 2])
    S["pool"] = dict(A_b=pool.A_b, t_b=pool.t_b, parity=pool.parity, midx=pool.midx, sigp=pool.sigp,
                     event=pool.event, n=pool.n)
    S["Rb"], S["Rt"], S["Rtx"], S["noise_f"] = C.Rb, C.Rt, C.Rtx, C.noise_f
    S["thr"], S["thr_tpl"] = C.thr, C.thr_tpl
    S["sing"] = singles_summary(C, np.random.default_rng([SEED, C.irt, 3]))
    t0 = time.time()
    # MF on the evaluation singles
    nm = 1500 if quick else len(ev_idx)
    tl, al = tf.run_mf(C.bank, pool.raw[ev_idx[:nm]], pool.midx[ev_idx[:nm]])
    S["mf_single_fake"] = (sum(len(t) >= 2 for t in tl), nm)
    log(f"[{C.rt}] MF on {nm} eval singles: {S['mf_single_fake']} with >=2 pulses; {time.time() - t0:.0f} s")
    cells = {}
    keys = cell_keys()
    mf_keys = {("real2", dt, r) for dt in ROC_DT for r in ROC_R} | \
              {(a, dt, r) for a in ("real2", "real") for dt in SM_DT for r in SM_R}
    for i, k in enumerate(keys):
        arm, dt, r = k
        cells[k] = run_cell(C, arm, dt, r, n, rng, mf=(k in mf_keys),
                            keepW=(k == ("real2", 8, 1.0)))
        if i % 10 == 0 or i == len(keys) - 1:
            log(f"[{C.rt}] cell {i + 1}/{len(keys)} {k}  {time.time() - t0:.0f} s")
    S["cells"] = cells
    # examples for page 1 (arrays copied out of the context)
    S["examples"] = build_examples(C, cells[("real2", 8, 1.0)])
    S["robust"] = robustness(C, cells[("real2", 8, 1.0)], np.random.default_rng([SEED, C.irt, 4]))
    S["indep"] = independence(C, cells[("real2", 8, 1.0)])
    return S


def compute_other(C, n, quick):
    S = dict(rt=C.rt)
    pool = C.pool
    rng = np.random.default_rng([SEED, C.irt, 2])
    S["Rb"], S["Rt"], S["Rtx"], S["noise_f"] = C.Rb, C.Rt, C.Rtx, C.noise_f
    S["thr"], S["thr_tpl"] = C.thr, C.thr_tpl
    S["pool"] = dict(A_b=pool.A_b, t_b=pool.t_b, parity=pool.parity, midx=pool.midx, sigp=pool.sigp, n=pool.n)
    S["sing"] = singles_summary(C, np.random.default_rng([SEED, C.irt, 3]))
    S["cells"] = {("real2", 8, r): run_cell(C, "real2", 8, r, n, rng) for r in (1.0, 0.3)}
    return S


# ============================================================ page-1 examples ==
def build_examples(C, cell):
    pool, bank = C.pool, C.bank
    ev, R, W = cell["ev"], cell["R"], cell["W"]
    y, _ = tf.subtract_ped(bank, W)
    d = R[:, 12]
    order = np.argsort(d)
    pick = [order[int(q * (len(d) - 1))] for q in (0.10, 0.50, 0.90)]
    out = []
    for lab, j in zip(("doublet, dchi2 10th pct", "doublet, dchi2 median", "doublet, dchi2 90th pct"), pick):
        m = pool.midx[ev["idx"][j]]
        out.append(dict(kind="dbl", label=lab, y=y[j], R=R[j], m=int(m), lut=bank.lut[m].astype(np.float64),
                        truth=(ev["A1"][j], ev["A2"][j], ev["dt_true"][j], ev["t1"][j], ev["t2"][j])))
    ie = np.nonzero(C.ev)[0]
    ds = C.Rb[ie, 12]
    thr = C.thr[F_PROD]
    o = np.argsort(ds)
    below = ie[o[np.searchsorted(ds[o], np.quantile(ds, 0.90))]]
    cand_fake = ie[o[int(0.995 * (len(o) - 1))]]
    med = ie[o[len(o) // 2]]
    for lab, i in (("single, median dchi2", med), ("single, 90th pct dchi2", below),
                   ("single, 99.5th pct (fake: dchi2 > thr)" if C.Rb[cand_fake, 12] > thr else "single, 99.5th pct", cand_fake)):
        m = pool.midx[i]
        out.append(dict(kind="sgl", label=lab, y=pool.raw[i] - pool.ped[i], R=C.Rb[i], m=int(m),
                        lut=bank.lut[m].astype(np.float64), truth=(C.Rb[i, 2], np.nan, np.nan, C.Rb[i, 3], np.nan)))
    for e in out:
        e["T_end"] = bank.T_end
    return out


# ===================================================== robustness (panel 6) ==
def robustness(C, cell, rng):
    pool, bank, xbank = C.pool, C.bank, C.xbank
    ev, W = cell["ev"], cell["W"]
    idx = ev["idx"]
    n_all = pool.n
    sh = np.asarray(bank.t50) - np.asarray(xbank.t50)
    conds = [("baseline (best-50)", bank, 0.0, False, np.zeros(len(sh))),
             ("Wassim bank", xbank, 0.0, False, -sh),
             ("noise x2", bank, 0.0, True, np.zeros(len(sh))),
             ("template +0.5 ns", shifted_bank(bank, +0.5), 0.0, False, np.full(len(sh), +0.5)),
             ("template -0.5 ns", shifted_bank(bank, -0.5), 0.0, False, np.full(len(sh), -0.5)),
             ("template -0.5 ns, ramped", shifted_bank(bank, -0.5, 2.0), 0.0, False, np.full(len(sh), -0.5))]
    rows = []
    thr_base = C.thr[F_PROD]
    for name, fb, _, noisy, tcorr in conds:
        raw_s, raw_d, sg_s, sg_d = pool.raw, W, pool.sigp, pool.sigp[idx]
        if noisy:
            def add(raw, sig, ii):
                # pool.noise draws sg^2 = v^2 - 1/12 and the final rint adds 1/12: added variance 3 (f sig)^2, f = measured
                # baseline std / sig_ped (about 2.3), so the total noise sigma (rounding included) is 2x the measured one
                v = np.sqrt(3.0) * C.noise_f * sig
                return np.rint(raw + pool.noise(rng, v, len(raw)))
            raw_s = add(pool.raw, pool.sigp, None)
            raw_d = add(W, pool.sigp[idx], None)
            sg_s, sg_d = 2.0 * pool.sigp, 2.0 * pool.sigp[idx]
        if noisy:
            def bstd(x):
                z = x[:, :tf.PED_N] - x[:, :tf.PED_N].mean(1, keepdims=True)
                return np.sqrt((z ** 2).sum(1) / (tf.PED_N - 1))
            chk = float(np.median(bstd(raw_s[:500]) / bstd(pool.raw[:500])))
            row_noise = chk
            log(f"[{C.rt}] noise x2 check: median baseline-std ratio noisy/original (500 singles) = {chk:.2f}")
        ys, _ = tf.subtract_ped(bank, raw_s)
        yd, _ = tf.subtract_ped(bank, raw_d)
        with forked(fb, 18) as pw:
            Rs = tf.fit_batch(pw, ys, pool.midx, sg_s)
            Rd = tf.fit_batch(pw, yd, pool.midx[idx], sg_d)
        for RR, mm in ((Rs, pool.midx), (Rd, pool.midx[idx])):
            for c in (3, 9, 11):
                RR[:, c] = RR[:, c] + tcorr[mm]
        cal_d = Rs[C.cal, 12]
        ev_d = Rs[C.ev, 12]
        thr_t = tf.thr_for(cal_d, F_PROD)
        m2 = match_t(Rd, ev["t1"], ev["t2"])
        row = dict(name=name, med_single=float(np.median(ev_d)), thr_cal_prod=thr_t, thr_base=thr_base)
        for tag, thr in (("prod", thr_base), ("tuned", thr_t)):
            k = int(((Rd[:, 12] > thr) & m2).sum())
            f = int((ev_d > thr).sum())
            row.update({f"eff_{tag}": k / len(idx), f"eff_{tag}_k": k, f"fake_{tag}": f / len(ev_d), f"fake_{tag}_k": f,
                        f"thr_{tag}": thr})
        row["n_eval"] = int(len(ev_d))
        if noisy:
            row["std_ratio"] = row_noise
        row["n_dbl"] = int(len(idx))
        rows.append(row)
        log(f"[{C.rt}] robust {name}: prod eff {row['eff_prod']:.3f} fake {row['fake_prod']:.4f}; "
            f"tuned thr {thr_t:.1f} eff {row['eff_tuned']:.3f} fake {row['fake_tuned']:.4f}")
    return rows


# ==================================================== independence (panel 7) ==
def independence(C, cell):
    """Is the pool split part of the best-50 template build?  Reconstruct the build pool."""
    slug, split = C.slug, C.split
    zf = tf.B50 / slug / f"templates_{slug}.npz"
    z = np.load(zf, allow_pickle=False)
    meta = json.loads(str(z["meta"]))
    inputs = meta["inputs"]
    keys = {(i["run"], i["split"]) for i in inputs}
    run_n, split_n = [int(x) for x in Path(split.replace("_", "/")).parts]
    info = dict(slug=slug, n_inputs=len(inputs), in_build_inputs=(run_n, split_n) in keys,
                runs=meta["runs"], split=split, argv_n_max=meta["argv"][meta["argv"].index("--n-max") + 1],
                sha1=[i["sha1"] for i in inputs if (i["run"], i["split"]) == (run_n, split_n)])
    # reconstruct holdout_energy and the energy-ordered fit pool over all build inputs
    D = dict(tier=[], b0=[], amp=[], E=[], midx=[], evn=[], run=[], split=[], ev=[], name=[])
    for i in inputs:
        p = np.load(i["path"], allow_pickle=False)
        n = len(p["amp"])
        D["tier"].append(p["tier"])
        D["b0"].append(p["band_b0"])
        D["amp"].append(p["amp"])
        D["E"].append(p["mf_energy"])
        D["midx"].append(p["midx"])
        D["evn"].append(p["event_number"])
        D["ev"].append(p["event"])
        D["run"].append(np.full(n, i["run"]))
        D["split"].append(np.full(n, i["split"]))
    D = {k: np.concatenate(v) for k, v in D.items() if v}
    mi_ref = np.asarray(z["module_index"])      # file midx is the global module index; npz rows are PbWO4 only
    lut_row = {int(m): i for i, m in enumerate(mi_ref)}
    D["midx"] = np.array([lut_row.get(int(m), -1) for m in D["midx"]], dtype=np.int64)
    n_max, frac = int(info["argv_n_max"]), 0.2
    n_hold = int(n_max * frac / (1.0 - frac) + 0.5)
    E = D["E"]
    idx = np.nonzero((D["tier"] >= 1) & D["b0"] & np.isfinite(E) & (D["midx"] >= 0))[0]
    order = idx[np.lexsort((idx, -E[idx], D["midx"][idx]))]
    ms = D["midx"][order]
    rank = np.arange(order.size) - np.searchsorted(ms, ms, side="left")
    hold = np.zeros(D["midx"].size, bool)
    hold[order[(rank >= n_max) & (rank < n_max + n_hold)]] = True
    hk = z["holdout_keys"]
    mine = np.stack([D["run"][hold], D["split"][hold], D["evn"][hold], mi_ref[D["midx"][hold]]], 1)
    sa = {tuple(r) for r in hk.tolist()}
    sb = {tuple(r) for r in mine.tolist()}
    info["holdout_reconstructed"] = (len(sb), len(sa), len(sa & sb))
    # fit pool per module: top n_max (by E, stable) of the chain pool that the module's template used
    fb = np.asarray(z["fallback"])
    amp = D["amp"]
    p0 = (~hold) & (D["tier"] >= 1) & D["b0"]
    p1 = (~hold) & (D["tier"] >= 1) & (amp >= 100.0) & (amp <= 2500.0)
    in_pool = np.zeros(D["midx"].size, bool)
    npool_chk = [0, 0]
    allmod = np.argsort(D["midx"], kind="stable")
    ms_all = D["midx"][allmod]
    starts = np.searchsorted(ms_all, np.unique(ms_all), side="left")
    ends = np.append(starts[1:], ms_all.size)
    npool_ref = np.asarray(z["n_pool"])
    for s0, e0 in zip(starts, ends):
        ii = allmod[s0:e0]
        mi = int(D["midx"][ii[0]])
        if mi < 0 or fb[mi] == 2:
            continue
        mask = (p0 if fb[mi] == 0 else p1)[ii]
        pool_i = ii[mask]
        o = np.argsort(-E[pool_i], kind="stable")
        sel = pool_i[o[:n_max]]
        in_pool[sel] = True
        npool_chk[0] += int(sel.size == npool_ref[mi])
        npool_chk[1] += 1
    info["n_pool_match"] = tuple(npool_chk)
    # keys of the pool-split pulses
    rs = (D["run"] == run_n) & (D["split"] == split_n)
    ks_pool = {(int(e), int(m)) for e, m in zip(D["ev"][rs & in_pool], D["midx"][rs & in_pool])}
    ks_hold = {(int(e), int(m)) for e, m in zip(D["ev"][rs & hold], D["midx"][rs & hold])}
    pool = C.pool
    key = [(int(e), int(m)) for e, m in zip(pool.event, pool.midx)]
    ip = np.array([k in ks_pool for k in key])
    ih = np.array([k in ks_hold for k in key])
    ev_m, cal_m = C.ev, C.cal
    info.update(n_pool_rows=int(pool.n), n_in_build_pool=int(ip.sum()), n_in_holdout=int(ih.sum()),
                n_in_build_pool_cal=int((ip & cal_m).sum()), n_in_build_pool_eval=int((ip & ev_m).sum()),
                n_in_holdout_eval=int((ih & ev_m).sum()),
                n_split_in_build_pool=int((rs & in_pool).sum()),
                A_in_med=float(np.median(pool.A_b[ip])) if ip.any() else np.nan, A_all_med=float(np.median(pool.A_b)))
    # same dchi2 behaviour inside / outside the build pool (eval singles)
    thr = C.thr[F_PROD]
    d = C.Rb[:, 12]
    a = ev_m & ip
    b = ev_m & ~ip & ~ih
    info["sgl_in"] = dict(n=int(a.sum()), med=float(np.median(d[a])) if a.any() else np.nan,
                          fake=float((d[a] > thr).mean()) if a.any() else np.nan, k=int((d[a] > thr).sum()))
    info["sgl_out"] = dict(n=int(b.sum()), med=float(np.median(d[b])), fake=float((d[b] > thr).mean()),
                           k=int((d[b] > thr).sum()))
    ia = ev_m & ip & (pool.A_b >= 1500) & (pool.A_b <= 2500)
    ib = ev_m & ~ip & ~ih & (pool.A_b >= 1500) & (pool.A_b <= 2500)
    for nm_, mk in (("sgl_in_A", ia), ("sgl_out_A", ib)):
        info[nm_] = dict(n=int(mk.sum()), med=float(np.median(d[mk])) if mk.any() else np.nan,
                         k=int((d[mk] > thr).sum()), fake=float((d[mk] > thr).mean()) if mk.any() else np.nan)
    # doublets: base in build pool or not
    ev = cell["ev"]
    base_in = ip[ev["idx"]]
    m2 = match_t(cell["R"], ev["t1"], ev["t2"]) & (cell["R"][:, 12] > thr)
    Ab_ = pool.A_b[ev["idx"]]
    hiA = (Ab_ >= 1500) & (Ab_ <= 2500)
    info["dbl_base_in_A"] = dict(n=int((base_in & hiA).sum()), k=int(m2[base_in & hiA].sum()))
    info["dbl_base_out_A"] = dict(n=int((~base_in & hiA).sum()), k=int(m2[~base_in & hiA].sum()))
    info["dbl_base_in"] = dict(n=int(base_in.sum()), k=int(m2[base_in].sum()), eff=float(m2[base_in].mean()) if base_in.any() else np.nan)
    info["dbl_base_out"] = dict(n=int((~base_in).sum()), k=int(m2[~base_in].sum()), eff=float(m2[~base_in].mean()))
    # halves / partners
    info["halves_disjoint"] = bool(np.all((pool.parity == 0) ^ (pool.parity == 1)))
    # partner check on the real2 cell
    part = ev["partner"]
    info["partner_same_half"] = bool(np.all(pool.parity[part] == pool.parity[ev["idx"]]))
    info["partner_same_module"] = bool(np.all(pool.midx[part] == pool.midx[ev["idx"]]))
    info["partner_not_base"] = bool(np.all(part != ev["idx"]))
    info["bases_in_eval"] = bool(np.all(pool.parity[ev["idx"]] == 1))
    info["partner_in_cal_count"] = int((pool.parity[part] == 0).sum())
    info["dup_rows"] = int(len(key) - len(set(key)))
    return info


# ============================================================ analysis tools ==
def thr_curve(S, cell, thr):
    """(fake_eval, eff, thr) for the thresholds thr (values); fake on the eval real singles."""
    evs = S["Rb"][S["pool"]["parity"] == 1, 12]
    ev = cell["ev"]
    m2 = match_t(cell["R"], ev["t1"], ev["t2"])
    d = cell["R"][:, 12]
    thr = np.asarray(thr, float)
    fk = np.array([(evs > t).mean() for t in thr])
    ef = np.array([((d > t) & m2).mean() for t in thr])
    return fk, ef, thr


def stats_at(S, cell, thr):
    evs = S["Rb"][S["pool"]["parity"] == 1, 12]
    ev = cell["ev"]
    m2 = match_t(cell["R"], ev["t1"], ev["t2"])
    d = cell["R"][:, 12]
    k = int(((d > thr) & m2).sum())
    return dict(k=k, n=len(d), eff=k / len(d), fake_k=int((evs > thr).sum()), fake=float((evs > thr).mean()),
                n_s=len(evs), n2=float((d > thr).mean()))


# ================================================================= plotting ==
def setup_mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "legend.fontsize": 6.5,
                         "axes.grid": True, "grid.alpha": 0.25, "figure.dpi": 100})
    return plt


def text_page(plt, pdf, title, lines, fs=7.0):
    fig = plt.figure(figsize=(11, 8.5))
    fig.text(0.04, 0.965, title, fontsize=10, weight="bold", va="top")
    fig.text(0.04, 0.93, "\n".join(lines), family="monospace", fontsize=fs, va="top")  # lines pre-wrapped by caller
    pdf.savefig(fig)
    plt.close(fig)


def fig_examples(plt, pdf, S):
    ex = S["examples"]
    fig = plt.figure(figsize=(13, 9.5))
    gs = fig.add_gridspec(5, 3, height_ratios=[3, 1, 0.45, 3, 1], hspace=0.08, wspace=0.22, top=0.93, bottom=0.05,
                          left=0.06, right=0.985)
    fine = np.arange(0, 400, 0.25)
    for n, e in enumerate(ex):
        r0 = 0 if e["kind"] == "dbl" else 3
        c0 = n % 3
        ax = fig.add_subplot(gs[r0, c0])
        axr = fig.add_subplot(gs[r0 + 1, c0], sharex=ax)
        R, row, T = e["R"], e["lut"], e["T_end"]

        def curve(A, t):
            return A * tf.tpl(row, None, fine - t, T)
        m1 = curve(R[2], R[3])
        two = R[13] == 1
        m2 = curve(R[8], R[9]) + (curve(R[10], R[11]) if two else 0)
        ax.plot(TS, e["y"], "ko", ms=3, label="data (ped-sub.)")
        ax.plot(fine, m1, color=COL["p1"], lw=1.1, ls="--", label="1-pulse fit")
        ax.plot(fine, m2, color=COL["p2"], lw=1.1, label="2-pulse fit")
        if two:
            ax.plot(fine, curve(R[8], R[9]), color=COL["p2"], lw=0.6, ls=":")
            ax.plot(fine, curve(R[10], R[11]), color=COL["p2"], lw=0.6, ls=":")
        tc = R[9] if two else R[3]
        for tt in (e["truth"][3], e["truth"][4]):
            if np.isfinite(tt):
                ax.axvline(tt, color="0.6", lw=0.5, ls="-.")
        lo, hi = tc - 12, tc + 70
        ax.set_xlim(lo, hi)
        ax.set_ylim(None, 1.45 * float(np.max(e["y"])))
        ax.set_ylabel("ADC")
        ax.set_title(e["label"], fontsize=8)
        tr = (f"truth A1 {e['truth'][0]:.0f}, A2 {e['truth'][1]:.0f}, dt {e['truth'][2]:.1f} ns"
              if e["kind"] == "dbl" else f"truth: single, A {R[2]:.0f} (= base fit), dt none")
        ax.text(0.98, 0.97, f"chi2_1 {R[4]:.0f}  chi2_2 {R[7]:.0f}  dchi2 {R[12]:.0f}\n"
                f"fit A1 {R[8]:.0f}  A2 {R[10]:.0f}  dt {(R[11] - R[9]) if two else float('nan'):.1f} ns\n{tr}",
                transform=ax.transAxes, ha="right", va="top", fontsize=6.5, family="monospace",
                bbox=dict(fc="w", ec="0.8", alpha=0.85))
        if n == 0:
            ax.legend(loc="center right", fontsize=6)
        i_in = (TS >= lo) & (TS <= hi)
        axr.plot(TS, e["y"] - np.interp(TS, fine, m1), "o-", color=COL["p1"], ms=2.5, lw=0.6, label="data - 1p")
        axr.plot(TS, e["y"] - np.interp(TS, fine, m2), "s-", color=COL["p2"], ms=2.5, lw=0.6, label="data - 2p")
        axr.axhline(0, color="k", lw=0.5)
        axr.set_ylabel("resid.")
        axr.set_xlabel("time (ns)")
        plt.setp(ax.get_xticklabels(), visible=False)
        if n == 0:
            axr.legend(loc="upper right", fontsize=6)
    fig.suptitle(f"Page 1: example waveforms, {S['rt']} (rows 1-2: real2 doublets dt 8 ns ratio 1; rows 3-4: eval singles). "
                 "grey dash-dot lines: truth onsets", fontsize=9)
    pdf.savefig(fig)
    plt.close(fig)


def hist_panel(ax, S, cells, dt, r, extra_real=None, title_extra="", ana=True):
    """dchi2 histograms of one case; returns the stats dict."""
    evm = S["pool"]["parity"] == 1
    s = S["Rb"][evm, 12]
    t = S["Rt"][evm, 12]
    cell = cells[("real2", dt, r)]
    d = cell["R"][:, 12]
    thr = S["thr"]
    bins = np.logspace(-1, 7, 61)
    cl = lambda v: np.clip(v, 0.1, None)
    w = lambda v: np.full(len(v), 1.0 / len(v))
    ax.hist(cl(s), bins, weights=w(s), histtype="stepfilled", color=COL["single"], alpha=0.35, label="real singles (eval)")
    ax.hist(cl(s), bins, weights=w(s), histtype="step", color=COL["single"], lw=1)
    ax.hist(cl(t), bins, weights=w(t), histtype="step", color=COL["tpl"], lw=1, label="tpl_noise singles")
    if "Rtx" in S:
        tx = S["Rtx"][evm, 12]
        ax.hist(cl(tx), bins, weights=w(tx), histtype="step", color=COL["tpl"], lw=1, ls=":", label="tpl_noise, noise at measured baseline std")
    if ana and "sing" in S:
        a = S["sing"]["null"]
        ax.hist(cl(a), bins, weights=w(a), histtype="step", color=COL["ana"], lw=1, ls="--", label="analytic kappa_c*chi2(k=2)")
    ax.hist(cl(d), bins, weights=w(d), histtype="stepfilled", color=COL["dbl"], alpha=0.35, label="real2 doublets (eval)")
    ax.hist(cl(d), bins, weights=w(d), histtype="step", color=COL["dbl"], lw=1)
    if extra_real is not None:
        dr = extra_real["R"][:, 12]
        ax.hist(cl(dr), bins, weights=w(dr), histtype="step", color=COL["real"], lw=1, ls=":", label="real (template inj.)")
    for f, ls, lab in ((F_PROD, "-", "fake<=3%"), (F_01, "--", "0.1%"), (F_002, ":", "0.02%")):
        ax.axvline(thr[f], color="k", ls=ls, lw=0.9, label=f"thr {lab} = {thr[f]:.0f}")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(1e-4, 1.6)
    ax.set_xlim(0.1, 1e7)
    ax.set_xlabel("dchi2 = chi2_1 - chi2_2")
    below = float((d <= thr[F_PROD]).mean())
    above = float((s > thr[F_PROD]).mean())
    ovl = overlap_coef(s, d)
    below01 = float((d <= thr[F_01]).mean())
    above01 = float((s > thr[F_01]).mean())
    ax.set_title(f"dt {dt} ns, ratio {r:g}", fontsize=8)
    ax.text(0.02, 0.98, f"thr 3%: dbl below {100 * below:.1f}%, sgl above {100 * above:.1f}%\n"
            f"thr 0.1%: dbl below {100 * below01:.1f}%, sgl above {100 * above01:.2f}%\nOVL {ovl:.4f}",
            transform=ax.transAxes, fontsize=6.0, family="monospace", va="top",
            bbox=dict(fc="w", ec="0.8", alpha=0.85))
    return dict(dt=dt, r=r, below=below, above=above, below01=below01, above01=above01, ovl=ovl,
                below_k=int((d <= thr[F_PROD]).sum()), n_d=len(d), above_k=int((s > thr[F_PROD]).sum()), n_s=len(s),
                below01_k=int((d <= thr[F_01]).sum()),
                med_dbl=float(np.median(d)))


def fig_hist(plt, pdf, S, dts, rs, title, big=True, extra=True):
    fig, axs = plt.subplots(len(rs), len(dts), figsize=(3.8 * len(dts), 3.6 * len(rs)), squeeze=False)
    out = []
    for i, r in enumerate(rs):
        for j, dt in enumerate(dts):
            xr = S["cells"].get(("real", dt, r)) if extra else None
            st = hist_panel(axs[i, j], S, S["cells"], dt, r, extra_real=xr if dt == 8 else None)
            out.append(st)
        axs[i, 0].set_ylabel("fraction per bin")
    h, l = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=4, fontsize=7)
    fig.suptitle(title, fontsize=9)
    fig.tight_layout(rect=(0, 0.05, 1, 0.96))
    pdf.savefig(fig)
    plt.close(fig)
    return out


def null_medians(S):
    evm = S["pool"]["parity"] == 1
    return dict(analytic=float(np.median(S["sing"]["null"])), analytic_kappa=float(np.nanmedian(S["sing"]["kappa"])),
                tpl_zero=float(np.mean(S["Rt"][evm, 12] == 0.0)),
                tpl=float(np.median(S["Rt"][evm, 12])), real=float(np.median(S["Rb"][evm, 12])),
                tpl_q97=float(np.quantile(S["Rt"][evm, 12], 0.97)),
                tplx=float(np.median(S["Rtx"][evm, 12])), tplx_q97=float(np.quantile(S["Rtx"][evm, 12], 0.97)),
                noise_f=float(S["noise_f"]), real_q97=float(np.quantile(S["Rb"][evm, 12], 0.97)),
                ana_q97=float(np.quantile(S["sing"]["null"], 0.97)))


def fig_fitparams(plt, pdf, S):
    cells = S["cells"]
    thr = S["thr"][F_PROD]
    fig, axs = plt.subplots(2, 3, figsize=(13, 7.5))
    cmap = plt.get_cmap("viridis")
    stats = []
    for i, r in enumerate((1.0, 0.3)):
        ax1, ax2, ax3 = axs[i]
        for k, dt in enumerate(P3_DT):
            c = cells[("real2", dt, r)]
            R, ev = c["R"], c["ev"]
            col = cmap(k / (len(P3_DT) - 1))
            det = (R[:, 12] > thr)
            fdt = R[:, 11] - R[:, 9]
            sep = ev["dt_true"]
            degen = det & ((fdt < 2) | (R[:, 10] / np.maximum(R[:, 8], 1e-9) < 0.05) |
                           (R[:, 8] / np.maximum(R[:, 10], 1e-9) < 0.05) | c["edge"])
            ok = det & match_t(R, ev["t1"], ev["t2"])
            sel = np.nonzero(det)[0][:250]
            ax1.plot(sep[sel], fdt[sel], ".", ms=2, color=col, alpha=0.5)
            dd = (fdt - sep)[ok]
            ax2.hist(np.clip(dd, -6, 6), np.linspace(-6, 6, 49), histtype="step", color=col, lw=1,
                     weights=np.full(ok.sum(), 1.0 / max(ok.sum(), 1)), label=f"dt {dt}")
            rt = (R[:, 10] / np.maximum(R[:, 8], 1e-9))[det]
            ax3.hist(np.clip(rt, 0.01, 100), np.logspace(-2, 2, 49), histtype="step", color=col, lw=1,
                     weights=np.full(det.sum(), 1.0 / max(det.sum(), 1)))
            stats.append(dict(r=r, dt=dt, n=len(R), n_det=int(det.sum()), n_ok=int(ok.sum()),
                              n_degen=int(degen.sum()), n_sep2=int((det & (fdt < 2)).sum()),
                              n_tinyA2=int((det & (R[:, 10] / np.maximum(R[:, 8], 1e-9) < 0.05)).sum()),
                              n_tinyA1=int((det & (R[:, 8] / np.maximum(R[:, 10], 1e-9) < 0.05)).sum()),
                              n_edge=int((det & c["edge"]).sum()),
                              ddt_bias=float(dd.mean()) if ok.any() else np.nan, ddt_rms=float(dd.std()) if ok.any() else np.nan,
                              ratio_med=float(np.median(rt)) if det.any() else np.nan,
                              ratio_rms=float(np.std(np.log(np.clip(rt, 1e-3, None)))) if det.any() else np.nan))
            ax3.axvline(r, color=col, lw=0.4, ls=":")
        ax1.plot([0, 30], [0, 30], "k-", lw=0.6)
        ax1.set_xlim(0, 30)
        ax1.set_ylim(0, 40)
        ax1.set_xlabel("true dt (ns)")
        ax1.set_ylabel("fitted t2 - t1 (ns)")
        ax1.set_title(f"ratio {r:g}: fitted vs true dt (detected, 250/dt shown)")
        ax2.set_xlabel("fitted dt - true dt (ns), detected & matched (clipped at +-6)")
        ax2.set_ylabel("fraction")
        ax2.set_title("dt residual")
        ax3.set_xscale("log")
        ax3.set_xlabel("fitted A2/A1 (A2 = later pulse); dotted: true ratio")
        ax3.set_ylabel("fraction")
        ax3.set_title("amplitude ratio")
        if i == 0:
            ax2.legend(fontsize=6, ncol=2)
    fig.suptitle(f"Page 3: real2 doublet fit parameters at the production threshold, {S['rt']}", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    pdf.savefig(fig)
    plt.close(fig)
    return stats


def fake_decomp(S, thr):
    """Secondary (smaller-amplitude fitted pulse) of eval singles that fire, vs main."""
    evm = S["pool"]["parity"] == 1
    R = S["Rb"][evm]
    A = S["pool"]["A_b"][evm]
    t = S["pool"]["t_b"][evm]
    fire = R[:, 12] > thr
    a1, a2, t1, t2 = R[:, 8], R[:, 10], R[:, 9], R[:, 11]
    sec_is_2 = a2 <= a1
    a_main = np.where(sec_is_2, a1, a2)
    a_sec = np.where(sec_is_2, a2, a1)
    t_main = np.where(sec_is_2, t1, t2)
    t_sec = np.where(sec_is_2, t2, t1)
    off = t_sec - t_main
    rat = a_sec / np.maximum(a_main, 1e-9)
    return dict(fire=fire, off=off, rat=rat, sep=t2 - t1, a_main=a_main, a_sec=a_sec, A=A, d=R[:, 12],
                sec_after=off > 0)


def fig_fakes(plt, pdf, S):
    thr = S["thr"][F_PROD]
    t01 = S["thr"][F_01]
    F = fake_decomp(S, thr)
    d = F["d"]
    f1, f2 = d > thr, d > t01
    fig, axs = plt.subplots(2, 3, figsize=(13, 7.5))
    ax = axs[0, 0]
    ax.plot(np.abs(F["off"])[~f1][:3000], F["rat"][~f1][:3000], ".", ms=1.5, color="0.7", label="no fire (3000 shown)")
    ax.plot(np.abs(F["off"])[f1 & ~f2], F["rat"][f1 & ~f2], ".", ms=3, color="#E69F00", label="fire: thr 3% only")
    ax.plot(np.abs(F["off"])[f2], F["rat"][f2], "o", ms=3.5, color="#D55E00", mfc="none", label="fire: > thr 0.1%")
    ax.axvline(2, color="k", lw=0.6, ls=":")
    ax.axhline(0.05, color="k", lw=0.6, ls=":")
    ax.set_yscale("log")
    ax.set_ylim(1e-4, 2)
    ax.set_xlabel("|t_sec - t_main| (ns), 2-pulse fit")
    ax.set_ylabel("A_sec / A_main")
    ax.set_title("eval singles: fitted secondary (dotted: 2 ns, 0.05)")
    ax.legend(fontsize=6)
    ax = axs[0, 1]
    b = np.linspace(-60, 80, 57)
    ax.hist(np.clip(F["off"][f1], -60, 80), b, histtype="step", color="#E69F00", lw=1.2, label=f"fire, thr 3% (n={int(f1.sum())})")
    ax.hist(np.clip(F["off"][f2], -60, 80), b, histtype="step", color="#D55E00", lw=1.2, label=f"fire, thr 0.1% (n={int(f2.sum())})")
    ax.set_xlabel("t_sec - t_main (ns); >0 after the main pulse")
    ax.set_ylabel("events")
    ax.set_title("where in time is the fake secondary")
    ax.legend(fontsize=6)
    ax = axs[0, 2]
    edges = np.array([100, 150, 200, 300, 450, 700, 1000, 1500, 2500, 4001])
    cen = np.sqrt(edges[:-1] * edges[1:])
    for fl, c, lab in ((f1, "#E69F00", "thr 3%"), (f2, "#D55E00", "thr 0.1%")):
        k = np.array([fl[(F["A"] >= lo) & (F["A"] < hi)].sum() for lo, hi in zip(edges[:-1], edges[1:])])
        n = np.array([((F["A"] >= lo) & (F["A"] < hi)).sum() for lo, hi in zip(edges[:-1], edges[1:])])
        ci = np.array([wil(a, b_) for a, b_ in zip(k, n)])
        p = k / np.maximum(n, 1)
        ax.errorbar(cen, p, [np.maximum(p - ci[:, 0], 0), np.maximum(ci[:, 1] - p, 0)], fmt="o-", ms=3, color=c, label=lab, capsize=1.5)
    ax.set_xscale("log")
    ax.set_xticks([100, 200, 500, 1000, 2000, 4000])
    ax.set_xticklabels(["100", "200", "500", "1000", "2000", "4000"])
    ax.set_xlabel("base amplitude A (ADC)")
    ax.set_ylabel("fake rate (eval singles)")
    ax.set_title("fake rate vs amplitude (Wilson 68%)")
    ax.legend(fontsize=6)
    ax = axs[1, 0]
    ax.hist(np.clip(F["sep"][f1], 0, 80), np.linspace(0, 80, 41), color="#E69F00", alpha=0.6, label="thr 3%")
    ax.hist(np.clip(F["sep"][f2], 0, 80), np.linspace(0, 80, 41), color="#D55E00", alpha=0.8, label="thr 0.1%")
    ax.set_xlabel("fitted t2 - t1 (ns)")
    ax.set_ylabel("events")
    ax.set_title("fitted separation of firing singles")
    ax.legend(fontsize=6)
    ax = axs[1, 1]
    bb = np.logspace(-2.5, 0, 41)
    ax.hist(F["rat"][f1], bb, color="#E69F00", alpha=0.6, label="thr 3%")
    ax.hist(F["rat"][f2], bb, color="#D55E00", alpha=0.8, label="thr 0.1%")
    ax.set_xscale("log")
    ax.set_xlabel("A_sec / A_main")
    ax.set_ylabel("events")
    ax.set_title("fitted amplitude ratio of firing singles")
    ax.legend(fontsize=6)
    ax = axs[1, 2]
    ax.plot(F["A"][~f1][:4000], d[~f1][:4000], ".", ms=1.5, color="0.7")
    ax.plot(F["A"][f1], d[f1], ".", ms=3, color="#E69F00")
    ax.plot(F["A"][f2], d[f2], "o", ms=3, mfc="none", color="#D55E00")
    ax.axhline(thr, color="k", lw=0.8)
    ax.axhline(t01, color="k", lw=0.8, ls="--")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(0.1, None)
    ax.set_xlabel("base amplitude A (ADC)")
    ax.set_ylabel("dchi2")
    ax.set_title("dchi2 vs amplitude (single-pulse dchi2 scales ~ A^2 eps^2)")
    fig.suptitle(f"Page 4: what do the single-pulse fakes look like? {S['rt']} eval singles", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    pdf.savefig(fig)
    plt.close(fig)
    out = {}
    for tag, fl in (("p3", f1), ("p01", f2)):
        n = int(fl.sum())
        sep2 = fl & (np.abs(F["off"]) < 2)
        tiny = fl & (F["rat"] < 0.05)
        degen = sep2 | tiny
        gen = fl & ~degen
        out[tag] = dict(n=n, n_sep2=int(sep2.sum()), n_tiny=int(tiny.sum()), n_degen=int(degen.sum()), n_gen=int(gen.sum()),
                        after=int((fl & F["sec_after"]).sum()), before=int((fl & ~F["sec_after"]).sum()),
                        off_med=float(np.median(F["off"][fl])) if n else np.nan,
                        absoff_med=float(np.median(np.abs(F["off"][fl]))) if n else np.nan,
                        rat_med=float(np.median(F["rat"][fl])) if n else np.nan,
                        gen_off_q=[float(x) for x in np.quantile(np.abs(F["off"][gen]), [0.1, 0.5, 0.9])] if gen.any() else None,
                        gen_rat_q=[float(x) for x in np.quantile(F["rat"][gen], [0.1, 0.5, 0.9])] if gen.any() else None,
                        A_med=float(np.median(F["A"][fl])) if n else np.nan, A_med_all=float(np.median(F["A"])))
    return out


def fig_roc(plt, pdf, S):
    cells = S["cells"]
    tgrid = np.logspace(0.5, 5.0, 160)
    cal0 = S["Rb"][S["pool"]["parity"] == 0, 12]
    evs = S["Rb"][S["pool"]["parity"] == 1, 12]
    nS = len(evs)
    fig, axs = plt.subplots(1, 2, figsize=(13, 5.6), sharey=True)
    cols = {1.0: "#0072B2", 0.5: "#009E73", 0.3: "#E69F00", 0.1: "#D55E00", 0.05: "#CC79A7"}
    table = []
    for ax, dt in zip(axs, ROC_DT):
        for r in ROC_R:
            c = cells.get(("real2", dt, r))
            if c is None:
                continue
            fk, ef, _ = thr_curve(S, c, tgrid)
            fk_plot = np.maximum(fk, 0.5 / nS)
            ax.plot(fk_plot, ef, "-", color=cols[r], lw=1.2, label=f"real2 r {r:g}")
            cr = cells.get(("real", dt, r))
            if cr is not None:
                fk2, ef2, _ = thr_curve(S, cr, tgrid)
                ax.plot(np.maximum(fk2, 0.5 / nS), ef2, "--", color=cols[r], lw=0.8, label=f"real r {r:g}" if r in (1.0, 0.3, 0.1) else None)
            thr3 = S["thr"][F_PROD]
            for f_, mk_ in ((F_01, "D"), (F_002, "s")):
                sq = stats_at(S, c, S["thr"][f_])
                ax.plot([max(sq["fake"], 0.5 / nS)], [sq["eff"]], mk_, color=cols[r], ms=4.5, mec="k", mew=0.6)
            st = stats_at(S, c, thr3)
            ax.plot([max(st["fake"], 0.5 / nS)], [st["eff"]], "o", color=cols[r], ms=7, mec="k")
            if "mf" in c:
                e = mf_eff(c).mean()
                mfk = S["mf_single_fake"]
                ax.plot([max(mfk[0] / mfk[1], 0.5 / nS)], [e], "*", color=cols[r], ms=11, mec="k")
            row = dict(dt=dt, r=r)
            for tag, f in (("p3", F_PROD), ("p01", F_01), ("p002", F_002)):
                t = tf.thr_for(S["Rb"][S["pool"]["parity"] == 0, 12], f)
                stt = stats_at(S, c, t)
                row[tag] = stt
            row["mf_eff"] = float(mf_eff(c).mean()) if "mf" in c else np.nan
            row["mf_k"] = int(mf_eff(c).sum()) if "mf" in c else 0
            table.append(row)
        ax.set_xscale("log")
        ax.axvline(0.001, color="0.5", lw=0.6, ls="--")
        ax.axvline(0.0002, color="0.5", lw=0.6, ls=":")
        ax.axvline(1.0 / nS, color="r", lw=0.5, ls="-")
        ax.set_xlabel("eval-half fake rate of real singles (threshold from cal-half quantile)")
        ax.set_title(f"dt {dt} ns (jit 2). circle = cal-quantile thr fake<=3 %, diamond = 0.1 %, square = 0.02 %, star = production MF", fontsize=7)
        ax.set_xlim(0.5 / nS, 0.5)
        ax.set_ylim(-0.02, 1.03)
    axs[0].set_ylabel("efficiency (both onsets +-4 ns)")
    axs[0].legend(fontsize=6, ncol=2, loc="lower right")
    axs[1].text(0.02, 0.05, f"red line: 1 event of {nS} eval singles (fake 0 plotted at 0.5 event); 0.02% = {0.0002 * nS:.1f} events",
                transform=axs[1].transAxes, fontsize=6.5, color="r")
    fig.suptitle(f"Page 5: ROC (efficiency vs eval fake rate), {S['rt']}; curves: threshold scanned in value; markers: thresholds from "
                 "the cal-half quantile; grey verticals 0.1 %, 0.02 %; solid = real2, dashed = real (template)", fontsize=8)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    pdf.savefig(fig)
    plt.close(fig)
    return table


def roc_lines(table, nS, n_cal):
    L = ["ROC read-off, real2 (eff +-4 ns [Wilson 68 %]); fake = eval-half rate at the cal-half quantile threshold;",
         f"statistical limit: eval half has {nS} singles, so 0.02 % = {0.0002 * nS:.1f} events (fake rate granularity "
         f"{100 / nS:.3f} %), and the cal-half quantile at 0.02 % rests on {0.0002 * n_cal:.1f} events.",
         "",
         " dt  ratio | fake<=3%: eff [CI]  fake  | 0.1 %: eff [CI]  fake(k) | 0.02 %: eff [CI]  fake(k) | MF eff"]
    for r in table:
        a, b, c = r["p3"], r["p01"], r["p002"]
        L.append(f"{r['dt']:3d}  {r['r']:5g} | {eff_str(a['k'], a['n']):>14s}  {100 * a['fake']:5.2f}% | "
                 f"{eff_str(b['k'], b['n']):>14s}  {100 * b['fake']:5.3f}% ({b['fake_k']}) | "
                 f"{eff_str(c['k'], c['n']):>14s}  {100 * c['fake']:5.3f}% ({c['fake_k']}) | {100 * r['mf_eff']:.0f}%")
    return L


def fig_robust(plt, pdf, rows, rt):
    L = [f"{rt}: real2 doublets dt 8 ns, ratio 1, jit 2 ns, eval half; singles = all eval real pulses.",
         "'production threshold' = baseline cal-half quantile (fake<=3 %); 'retuned' = cal-half quantile (fake<=3 %) of the same condition.",
         "eff = both onsets +-4 ns [Wilson 68 %].  Onsets of the shifted / Wassim fits are mapped back to the best-50 frame.",
         "Noise x2: independent measured-noise added to doublets and singles (total baseline std 2x the measured sample std, = 2.3 sig_ped x 2), fit sigma_ped doubled.",
         "Template +-0.5 ns: a global shift is absorbed by the free onset, so this is a null check, not a robustness test.",
         "",
         f"{'condition':26s} | {'thr(prod)':>9s} {'eff':>14s} {'fake':>14s} | {'thr(retuned)':>12s} {'eff':>14s} {'fake':>14s} | med dchi2 singles"]
    for r in rows:
        L.append(f"{r['name']:26s} | {r['thr_prod']:9.1f} {eff_str(r['eff_prod_k'], r['n_dbl']):>14s} "
                 f"{100 * r['fake_prod']:6.2f}% ({r['fake_prod_k']:3d}) | {r['thr_tuned']:12.1f} "
                 f"{eff_str(r['eff_tuned_k'], r['n_dbl']):>14s} {100 * r['fake_tuned']:6.2f}% ({r['fake_tuned_k']:3d}) | {r['med_single']:.1f}")
    return L


def fig_small_eff(plt, pdf, S):
    cells = S["cells"]
    thr = S["thr"][F_PROD]
    evs = S["Rb"][S["pool"]["parity"] == 1, 12]
    fake_fit = float((evs > thr).mean())
    mfk = S["mf_single_fake"]
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.4))
    tab = []
    for ax, r in zip(axs[:2], SM_R):
        for arm, ls, lab in (("real2", "-", "real2 (partner)"), ("real", "--", "real (template)")):
            e_fit, lo_fit, hi_fit, e_mf, e_true, e_mim = [], [], [], [], [], []
            for dt in SM_DT:
                c = cells[(arm, dt, r)]
                R, ev = c["R"], c["ev"]
                det = R[:, 12] > thr
                m2 = det & match_t(R, ev["t1"], ev["t2"])
                Tb, Ab, Ts, As = truth_small(ev)
                mt = det & match_amp(R, Tb, Ab, Ts, As)
                k = int(m2.sum())
                lo, hi = wil(k, len(R))
                e_fit.append(k / len(R))
                lo_fit.append(lo)
                hi_fit.append(hi)
                e_mf.append(mf_eff(c).mean())
                e_true.append(mt.mean())
                # singles that mimic the injected configuration (either side of each eval single)
                evm = S["pool"]["parity"] == 1
                Rs = S["Rb"][evm]
                Aq, tq = S["pool"]["A_b"][evm], S["pool"]["t_b"][evm]
                fire = Rs[:, 12] > thr
                mim = 0.5 * sum(float((fire & match_amp(Rs, tq, Aq, tq + sg * dt, r * Aq)).mean()) for sg in (-1, 1))
                e_mim.append(mim)
                tab.append(dict(arm=arm, r=r, dt=dt, eff=k / len(R), k=k, n=len(R), mf=e_mf[-1], true_det=float(mt.mean()),
                                mim=mim, n2=float(det.mean()), snr2=float(np.median(c["snr2"])),
                                snr2_q=[float(x) for x in np.quantile(c["snr2"], [0.1, 0.9])],
                                mismatch=float((det & ~mt).mean()),
                                mf_fake=mfk[0] / mfk[1]))
            e_fit = np.array(e_fit)
            ax.errorbar(SM_DT, e_fit, [e_fit - lo_fit, np.array(hi_fit) - e_fit], color=COL["dbl"] if arm == "real2" else COL["real"],
                        ls=ls, marker="o", ms=3, capsize=1.5, label=f"fit, {lab}")
            ax.plot(SM_DT, e_mf, color="k", ls=ls, marker="x", ms=4, label=f"MF, {lab}")
            ax.plot(SM_DT, e_true, color=COL["tpl"], ls=ls, marker="s", ms=3, label=f"fit, secondary on truth (A within x2), {lab}")
            ax.plot(SM_DT, e_mim, color="0.5", ls=ls, marker="^", ms=3, label=f"singles mimicking the injected pair, {lab}" if arm == "real2" else None)
        ax.set_title(f"ratio {r:g}  (fit fake {100 * fake_fit:.1f} %, MF fake {100 * mfk[0] / mfk[1]:.2f} % on eval singles)")
        ax.set_xlabel("onset separation dt (ns)")
        ax.set_ylabel("efficiency (+-4 ns)")
        ax.set_ylim(-0.02, 1.03)
        if r == SM_R[0]:
            ax.legend(fontsize=5.0, loc="lower right")
    ax = axs[2]
    for arm, ls in (("real2", "-"), ("real", "--")):
        for r, c in zip(SM_R, ("#D55E00", "#CC79A7")):
            sn = [np.median(cells[(arm, dt, r)]["snr2"]) for dt in SM_DT]
            ax.plot(SM_DT, sn, ls=ls, marker="o", ms=3, color=c, label=f"{arm} r {r:g}")
    ax.axhline(thr, color="k", lw=0.8, label=f"thr 3% = {thr:.0f}")
    ax.axhline(S["thr"][F_01], color="k", lw=0.8, ls="--", label=f"thr 0.1% = {S['thr'][F_01]:.0f}")
    ax.set_yscale("log")
    ax.set_xlabel("dt (ns)")
    ax.set_ylabel("median small-pulse SNR^2 (upper bound)")
    ax.set_title("noise-only signal strength of the small pulse vs thresholds")
    ax.legend(fontsize=6)
    fig.suptitle(f"Small ratios: fit vs MF, {S['rt']} (eval half)", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    pdf.savefig(fig)
    plt.close(fig)
    return tab


def small_bias(S):
    thr = S["thr"][F_PROD]
    rows = []
    for arm in ("real2", "real"):
        for r in SM_R:
            for dt in SM_DT:
                c = S["cells"][(arm, dt, r)]
                R, ev = c["R"], c["ev"]
                det = R[:, 12] > thr
                Tb, Ab, Ts, As = truth_small(ev)
                ok = det & match_amp(R, Tb, Ab, Ts, As)
                if ok.sum() < 5:
                    rows.append(dict(arm=arm, r=r, dt=dt, n=int(ok.sum()))) 
                    continue
                t1, t2, a1, a2 = R[:, 9], R[:, 11], R[:, 8], R[:, 10]
                # fitted pulses assigned to truth: the one nearer in time to the small truth onset is "small"
                small_is2 = np.abs(t2 - Ts) <= np.abs(t1 - Ts)
                a_s = np.where(small_is2, a2, a1)
                a_b = np.where(small_is2, a1, a2)
                t_s = np.where(small_is2, t2, t1)
                t_b = np.where(small_is2, t1, t2)
                dt_f = np.abs(t_s - t_b)
                dtt = np.abs(Ts - Tb)
                rr = (a_s / As)[ok]
                ar = (a_s / np.maximum(a_b, 1e-9))[ok]
                rows.append(dict(arm=arm, r=r, dt=dt, n=int(ok.sum()), dt_bias=float((dt_f - dtt)[ok].mean()),
                                 dt_rms=float((dt_f - dtt)[ok].std()), as_bias=float(rr.mean() - 1), as_rms=float(rr.std()),
                                 ratio_bias=float(np.median(ar) / r - 1), ratio_rms=float(ar.std() / r),
                                 ab_bias=float((a_b / Ab)[ok].mean() - 1)))
    return rows


def fig_small_hist(plt, pdf, S):
    fig, axs = plt.subplots(2, 2, figsize=(9.5, 7.5))
    out = []
    for i, r in enumerate(SM_R):
        for j, dt in enumerate((8, 20)):
            ax = axs[i, j]
            xr = S["cells"].get(("real", dt, r))
            out.append(hist_panel(ax, S, S["cells"], dt, r, extra_real=xr))
    h, l = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, fontsize=6.5)
    fig.suptitle(f"Small ratios: dchi2 singles vs real2 doublets, {S['rt']} (eval half)", fontsize=9)
    fig.tight_layout(rect=(0, 0.08, 1, 0.95))
    pdf.savefig(fig)
    plt.close(fig)
    return out


def fig_other(plt, pdf, SO):
    fig, axs = plt.subplots(len(SO), 2, figsize=(10.5, 3.8 * len(SO)), squeeze=False)
    out = {}
    for i, S in enumerate(SO):
        for j, r in enumerate((1.0, 0.3)):
            out[(S["rt"], r)] = hist_panel(axs[i, j], S, S["cells"], 8, r)
            axs[i, j].set_title(f"{S['rt']}: dt 8 ns, ratio {r:g}", fontsize=8)
        axs[i, 0].set_ylabel("fraction per bin")
    h, l = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, fontsize=7)
    fig.suptitle("0.73 and 3.49 GeV H2: panel-2 repeat (eval half, thresholds from the cal half of each type)", fontsize=9)
    fig.tight_layout(rect=(0, 0.06, 1, 0.96))
    pdf.savefig(fig)
    plt.close(fig)
    return out


def indep_lines(S):
    I = S["indep"]
    hk = I["holdout_reconstructed"]
    L = [f"Pool split {I['split']} ({I['slug']}); best-50 bank built from {I['n_inputs']} pulse files, runs {I['runs']}.",
         f"Pool split among the build inputs (templates npz meta.inputs, sha1 {I['sha1'][0][:10]}): {I['in_build_inputs']}.",
         f"Build = top {I['argv_n_max']} by mf_energy per module of the eligible candidates of ALL inputs (--pool-order energy), hold-out = next 13.",
         f"Reconstruction check: hold-out keys reproduced {hk[2]}/{hk[1]} (rebuilt {hk[0]}); module pool size == n_pool for {I['n_pool_match'][0]}/{I['n_pool_match'][1]} modules.",
         f"Pool-split pulses in the build fit-pool (pre-outlier-rejection): {I['n_split_in_build_pool']} in the split file; "
         f"{I['n_in_build_pool']} of the {I['n_pool_rows']} diagnostic pool rows ({100 * I['n_in_build_pool'] / I['n_pool_rows']:.1f} %): "
         f"cal {I['n_in_build_pool_cal']}, eval {I['n_in_build_pool_eval']}.",
         f"Pool rows in the build hold-out list (not used to build): {I['n_in_holdout']} (eval {I['n_in_holdout_eval']}).",
         f"Median fitted amplitude: in-bank {I['A_in_med']:.0f} ADC, all pool {I['A_all_med']:.0f} ADC.",
         f"Amplitude-matched (1500-2500 ADC) eval singles: in-bank n {I['sgl_in_A']['n']}, median dchi2 {I['sgl_in_A']['med']:.1f}, fake {I['sgl_in_A']['k']} ({100 * I['sgl_in_A']['fake']:.1f} %); "
         f"out-of-bank n {I['sgl_out_A']['n']}, median {I['sgl_out_A']['med']:.1f}, fake {I['sgl_out_A']['k']} ({100 * I['sgl_out_A']['fake']:.1f} %).",
         f"Eval singles in the build pool: n {I['sgl_in']['n']}, median dchi2 {I['sgl_in']['med']:.1f}, fake {I['sgl_in']['k']}/{I['sgl_in']['n']} "
         f"({100 * I['sgl_in']['fake']:.1f} %); eval singles outside pool and hold-out: n {I['sgl_out']['n']}, median {I['sgl_out']['med']:.1f}, "
         f"fake {I['sgl_out']['k']}/{I['sgl_out']['n']} ({100 * I['sgl_out']['fake']:.1f} %).",
         f"real2 dt 8 r 1 efficiency: base in build pool {eff_str(I['dbl_base_in']['k'], I['dbl_base_in']['n'], 1)} % (n {I['dbl_base_in']['n']}), "
         f"base not in pool {eff_str(I['dbl_base_out']['k'], I['dbl_base_out']['n'], 1)} % (n {I['dbl_base_out']['n']}); "
         f"amplitude-matched (base 1500-2500 ADC): in {I['dbl_base_in_A']['k']}/{I['dbl_base_in_A']['n']}, out {I['dbl_base_out_A']['k']}/{I['dbl_base_out_A']['n']}.",
         f"Cal / eval halves disjoint pool rows (parity of the row index): {I['halves_disjoint']}; duplicated (event, module) keys in the pool: {I['dup_rows']}.",
         f"real2 bases all eval: {I['bases_in_eval']}; partner same module {I['partner_same_module']}, same half {I['partner_same_half']}, "
         f"partner != base {I['partner_not_base']}, partners from the cal half: {I['partner_in_cal_count']}.",
         "Thresholds use only parity-0 rows (thr_for on the cal half of real singles)."]
    return L


def fig_indep(plt, pdf, S, C_extra=None):
    I = S["indep"]
    lines = indep_lines(S)
    text_page(plt, pdf, f"Page 7: independence of the evaluation pulses from the template build, {S['rt']}", textwrap_lines(lines, 140), fs=8.0)


def textwrap_lines(lines, w):
    import textwrap
    out = []
    for l in lines:
        out += textwrap.wrap(l, w, subsequent_indent="    ") or [""]
    return out


def fig_e(plt, pdf, SS):
    L = ["(e) r = 0.05 second pulse: amplitude of a median first pulse versus noise and template-residual scale (eval singles).",
         "A_med = median fitted amplitude of the pool; sigma = module pedestal rms (sig_ped) median; residual = data - 1-pulse fit on the pulse body",
         "(model > 5 % of peak): RMS and max |residual| in ADC; per-event ratios 0.05*A / quantity, median over events.", ""]
    L.append(f"{'run type':14s} {'A_med':>7s} {'A_s=.05A':>9s} {'sig_ped':>8s} {'A_s/sig':>8s} {'rms res/A':>10s} {'A_s/rms':>8s} "
             f"{'max res/A':>10s} {'A_s/max':>8s}")
    rows = []
    for S in SS:
        g = S["sing"]
        rows.append(dict(rt=S["rt"], **{k: g[k] for k in ("A_med", "sig_med", "rms_rel_med", "max_rel_med", "ratio_rms_med", "ratio_max_med", "ratio_sig_med")}))
        L.append(f"{S['rt']:14s} {g['A_med']:7.0f} {0.05 * g['A_med']:9.1f} {g['sig_med']:8.2f} {g['ratio_sig_med']:8.1f} "
                 f"{100 * g['rms_rel_med']:9.2f}% {g['ratio_rms_med']:8.2f} {100 * g['max_rel_med']:9.2f}% {g['ratio_max_med']:8.2f}")
    text_page(plt, pdf, "Small ratios (e): second-pulse amplitude vs noise and template residual", textwrap_lines(L, 170), fs=7.5)
    return rows


def table_small_lines(tab, bias):
    L = ["Small ratios, real2 (partner = nearest-amplitude same-module real pulse rescaled) and real (template injection); "
         "fit = production threshold fake<=3 %; eff = +-4 ns pair match;",
         "true = fit detects AND secondary within +-4 ns of the injected one with amplitude within x2; mismatch = fires but not true; "
         "mimic = eval singles whose fired fit matches the injected configuration (either side).",
         f"{'arm':5s} {'r':>5s} {'dt':>3s} | {'fit eff [CI]':>14s} {'true':>5s} {'mismatch':>8s} {'mimic':>6s} {'n2 rate':>7s} | {'MF eff':>6s} | {'SNR2 med [10-90]':>20s}"]
    for t in tab:
        L.append(f"{t['arm']:5s} {t['r']:5.2f} {t['dt']:3d} | {eff_str(t['k'], t['n']):>14s} {100 * t['true_det']:4.0f}% "
                 f"{100 * t['mismatch']:7.1f}% {100 * t['mim']:5.2f}% {100 * t['n2']:6.1f}% | {100 * t['mf']:5.0f}% | "
                 f"{t['snr2']:8.1f} [{t['snr2_q'][0]:.1f}-{t['snr2_q'][1]:.0f}]")
    return L


def table_bias_lines(bias):
    L = ["Detected small-ratio doublets (fit fires AND secondary on truth): dt and amplitude bias / RMS.",
         "dt = |t_small - t_big| fitted minus true (ns); A_s bias / RMS of A_s/A_s_true - 1; ratio = fitted A_s/A_b over true r, median bias, RMS/r.",
         f"{'arm':5s} {'r':>5s} {'dt':>3s} | {'n':>5s} {'dt bias':>8s} {'dt rms':>7s} {'A_s bias':>9s} {'A_s rms':>8s} {'ratio bias':>10s} {'ratio rms':>9s}"]
    for b in bias:
        if "dt_bias" not in b:
            L.append(f"{b['arm']:5s} {b['r']:5.2f} {b['dt']:3d} | {b['n']:5d}   (too few)")
        else:
            L.append(f"{b['arm']:5s} {b['r']:5.2f} {b['dt']:3d} | {b['n']:5d} {b['dt_bias']:+8.2f} {b['dt_rms']:7.2f} "
                     f"{100 * b['as_bias']:+8.1f}% {100 * b['as_rms']:7.1f}% {100 * b['ratio_bias']:+9.1f}% {100 * b['ratio_rms']:8.1f}%")
    return L


def table_p2_lines(st, rt, nm):
    L = [f"{rt}: dchi2 separation singles (eval) vs real2 doublets (eval).  thr 3% = production (fake<=0.03 cal quantile), 0.1 %, 0.02 % from the same cal half.",
         f"Medians of single-pulse dchi2: analytic kappa_c*chi2(k=2) {nm['analytic']:.3f} (median kappa_c {nm['analytic_kappa']:.3f}; chi2(2) median 1.39), "
         f"tpl_noise singles {nm['tpl']:.2f}, real singles {nm['real']:.2f}.",
         f"97th percentile: analytic {nm['ana_q97']:.1f}, tpl_noise {nm['tpl_q97']:.1f}, real {nm['real_q97']:.1f}.",
         "Analytic null scaling: expected gain of one extra weighted template fit, E = kappa_c * 2 for k = 2 (A2, t2) dof: kappa_c = sg^2 (W m)^T C (W m) / (m^T W m), m = fitted single-pulse model,",
         "W = diag(1/(sig_ped^2 + (0.02 y)^2)) the fit weights (the 2 % eps term makes the model sigma >> the true noise sg on the pulse), C = measured pedestal noise correlation (alternating sign: acf1 -0.33) + 1/12 ADC^2 rounding noise;",
         "kappa_c averaged over extra-pulse offsets t2 - t1 = -32..+64 ns; the extra shape is weight-orthogonalised against the single-pulse columns {s, ds/dt}.  null = kappa_c * chi2(2).  It ignores A2 >= 0 and the look-elsewhere over t2.",
         f"Empirical tpl_noise singles: {100 * nm['tpl_zero']:.0f} % of them have dchi2 exactly 0 (production search returns a pair only if its 0.5 ns-grid gain beats the 1-pulse grid gain), "
         "so the empirical median lies BELOW chi2(2) while its tail lies above: the relevant null is tpl_noise (no template mismatch), shown in the plots.",
         "", f"{'dt':>3s} {'r':>5s} | {'dbl below 3%':>12s} {'sgl above 3%':>12s} | {'dbl below .1%':>13s} {'sgl above .1%':>13s} | {'OVL':>7s} | med dchi2 dbl"]
    for s in st:
        L.append(f"{s['dt']:3d} {s['r']:5g} | {100 * s['below']:11.1f}% {100 * s['above']:11.2f}% | {100 * s['below01']:12.1f}% "
                 f"{100 * s['above01']:12.3f}% | {s['ovl']:7.4f} | {s['med_dbl']:.3g}")
    return L


# ===================================================================== main ==
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--quick", action="store_true", help="smoke test: few events")
    ap.add_argument("--replot", action="store_true", help="redo pdf / md from the pickle cache")
    ap.add_argument("--cache", default=str(SCRATCH))
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--n", type=int, default=2000, help="doublets per cell")
    ap.add_argument("--workers", type=int, default=18)
    ap.add_argument("--rts", default=",".join(ALL_RT))
    a = ap.parse_args(argv)
    out = Path(a.out)
    cache = Path(a.cache)
    cache.mkdir(parents=True, exist_ok=True)
    n = 120 if a.quick else a.n
    tag = "quick" if a.quick else "full"
    pk = cache / f"diag_state_{tag}.pkl"
    tf.OPT.update(OPT_PROD)
    t_all = time.time()
    if a.replot:
        SS = pickle.load(open(pk, "rb"))
    else:
        SS = {}
        for rt in a.rts.split(","):
            C = build_ctx(rt, a.workers, cache, a.quick)
            SS[rt] = compute_main(C, n, SEED, a.quick) if rt == MAIN else compute_other(C, n, a.quick)
            close_ctx(C)
            tf._BANK = None
            pickle.dump(SS, open(pk, "wb"))
            log(f"[{rt}] done  total {time.time() - t_all:.0f} s")
    SS["wall_compute"] = time.time() - t_all
    write_outputs(SS, out, a, n)
    log(f"[done] wall {time.time() - t_all:.0f} s")


def write_outputs(SS, out, a, n):
    plt = setup_mpl()
    from matplotlib.backends.backend_pdf import PdfPages
    out.mkdir(parents=True, exist_ok=True)
    S = SS[MAIN]
    pages = []
    R = {}
    with PdfPages(out / "overfit_diagnostics.pdf") as pdf:
        def mark(s):
            pages.append(s)
        # page 1
        fig_examples(plt, pdf, S)
        mark("1 example waveforms (3 real2 doublets, 3 eval singles)")
        # page 2
        st = fig_hist(plt, pdf, S, P2_DT, P2_R, f"Page 2: dchi2 of eval singles vs real2 doublets, {MAIN}; solid / dashed / dotted vertical: "
                      "thr fake<=3 % (production), 0.1 %, 0.02 % (cal-half quantiles); OVL = overlap coefficient")
        mark("2 dchi2 histograms, dt 4/8/12/20 x ratio 1/0.3")
        nm = null_medians(S)
        text_page(plt, pdf, "Page 2 (cont.): overlap statistics and null medians", textwrap_lines(table_p2_lines(st, MAIN, nm), 170), fs=7.0)
        mark("2b overlap table")
        R["p2"], R["nm"] = st, nm
        # page 3
        R["p3"] = fig_fitparams(plt, pdf, S)
        mark("3 fitted dt / A2/A1 vs truth (real2)")
        L = [f"Degenerate-solution counts among detected (dchi2 > thr 3 %) real2 doublets, {MAIN}; degenerate = fitted dt < 2 ns OR A2/A1 < 0.05 OR A1/A2 < 0.05 OR "
             "t1 / separation on a search / min_sep edge.",
             f"{'ratio':>5s} {'dt':>3s} | {'n':>5s} {'detected':>8s} {'matched':>8s} | {'dt<2':>5s} {'A2/A1<.05':>9s} {'A1/A2<.05':>9s} {'edge':>5s} {'any degen':>9s} | dt bias / rms (ns)  A2/A1 median"]
        for s in R["p3"]:
            L.append(f"{s['r']:5g} {s['dt']:3d} | {s['n']:5d} {s['n_det']:8d} {s['n_ok']:8d} | {s['n_sep2']:5d} {s['n_tinyA2']:9d} {s['n_tinyA1']:9d} {s['n_edge']:5d} "
                     f"{s['n_degen']:9d} | {s['ddt_bias']:+.2f} / {s['ddt_rms']:.2f}      {s['ratio_med']:.3f}")
        text_page(plt, pdf, "Page 3 (cont.): degenerate-solution counts", textwrap_lines(L, 170), fs=7.0)
        mark("3b degeneracy table")
        # page 4
        R["p4"] = fig_fakes(plt, pdf, S)
        mark("4 fake (single) fit parameters")
        # page 5
        R["roc"] = fig_roc(plt, pdf, S)
        nS = int((S["pool"]["parity"] == 1).sum())
        nC = int((S["pool"]["parity"] == 0).sum())
        mark("5 ROC")
        text_page(plt, pdf, "Page 5 (cont.): ROC read-off", textwrap_lines(roc_lines(R["roc"], nS, nC), 170), fs=7.0)
        mark("5b ROC table")
        # page 6
        L6 = fig_robust(plt, pdf, S["robust"], MAIN)
        text_page(plt, pdf, "Page 6: robustness at dt 8 ns, ratio 1 (real2)", textwrap_lines(L6, 200), fs=7.2)
        mark("6 robustness table")
        # page 7
        fig_indep(plt, pdf, S)
        mark("7 independence")
        # extra types
        SO = [SS[rt] for rt in OTHERS if rt in SS]
        R["other"] = {}
        if SO:
            R["other"] = fig_other(plt, pdf, SO)
            mark("8 0.73 / 3.49 GeV H2 panel-2 repeat (dt 8, r 1 / 0.3)")
            L = []
            for S2 in SO:
                nm2 = null_medians(S2)
                L.append(f"{S2['rt']}: medians analytic {nm2['analytic']:.2f}, tpl_noise {nm2['tpl']:.2f}, real {nm2['real']:.2f}; thr 3% {S2['thr'][F_PROD]:.1f}, "
                         f"0.1% {S2['thr'][F_01]:.1f}, 0.02% {S2['thr'][F_002]:.1f}")
                for r in (1.0, 0.3):
                    s = R["other"][(S2["rt"], r)]
                    L.append(f"   dt 8 r {r:g}: dbl below 3% {100 * s['below']:.1f}%, sgl above 3% {100 * s['above']:.2f}%, dbl below 0.1% {100 * s['below01']:.1f}%, "
                             f"sgl above 0.1% {100 * s['above01']:.3f}%, OVL {s['ovl']:.4f}")
            text_page(plt, pdf, "Page 8 (cont.): 0.73 / 3.49 GeV H2 statistics", textwrap_lines(L, 170), fs=7.5)
            mark("8b statistics")
        # small ratios
        R["small_tab"] = fig_small_eff(plt, pdf, S)
        mark("9 small ratios: efficiency fit vs MF, true vs mismatch, SNR^2")
        for arm in ("real2", "real"):
            text_page(plt, pdf, f"Small ratios (a, d): per-cell table, arm {arm}",
                      textwrap_lines(table_small_lines([t for t in R["small_tab"] if t["arm"] == arm], None), 170), fs=7.0)
            mark(f"9b small-ratio table ({arm})")
        R["small_bias"] = small_bias(S)
        text_page(plt, pdf, "Small ratios (d): bias and RMS of detected pairs", textwrap_lines(table_bias_lines(R["small_bias"]), 170), fs=7.0)
        mark("9c bias table")
        R["small_hist"] = fig_small_hist(plt, pdf, S)
        mark("10 small ratios: dchi2 histograms dt 8 / 20")
        SE = [SS[rt] for rt in ALL_RT if rt in SS]
        R["e"] = fig_e(plt, pdf, SE)
        mark("11 (e) amplitude vs noise and template residual")
    write_md(SS, R, out, n)
    log("PDF pages:")
    for p in pages:
        log("  ", p)


# ======================================================================== md ==
def write_md(SS, R, out, n):
    S = SS[MAIN]
    nm = R["nm"]
    P = {(s["dt"], s["r"]): s for s in R["p2"]}
    thr = S["thr"]
    nS = int((S["pool"]["parity"] == 1).sum())
    nC = int((S["pool"]["parity"] == 0).sum())
    roc = {(r["dt"], r["r"]): r for r in R["roc"]}
    p4 = R["p4"]
    I = S["indep"]
    rob = S["robust"]
    tab = {(t["arm"], t["r"], t["dt"]): t for t in R["small_tab"]}
    p8, p03 = P[(8, 1.0)], P[(8, 0.3)]
    pf, pf1 = p4["p3"], p4["p01"]
    p3_8 = [x for x in R["p3"] if x["dt"] == 8 and x["r"] == 1.0][0]
    p3_8b = [x for x in R["p3"] if x["dt"] == 8 and x["r"] == 0.3][0]
    cal_top = np.sort(S["Rb"][S["pool"]["parity"] == 0, 12])[::-1][:3]
    ev_max = float(np.max(S["Rb"][S["pool"]["parity"] == 1, 12]))
    rb = {r["name"]: r for r in rob}
    b_, p5_, m5_ = rb["baseline (best-50)"], rb["template +0.5 ns"], rb["template -0.5 ns"]
    ret = [r for k_, r in rb.items() if k_ != "baseline (best-50)"]
    L = []
    L.append("# U0 two-pulse fit: over-fitting diagnostics (2026-10-27)")
    L.append(f"`twopulse_diagnostics.py` -> `overfit_diagnostics.pdf` (page list in the script log). Production fit (eps 0.02, search 80, min_sep 1); 2.24 GeV H2 unless stated; "
             f"evaluation half only ({nS} real singles, {n} doublets/cell, jit 2 ns), thresholds from the calibration half ({nC}). Primary doublet arm real2 (two real pulses); `real` (noise-free template injected) is optimistic.")
    L.append(f"**Verdict.** No sign of noise over-fitting. At dt >= 8 ns (r >= 0.3) the doublet dchi2 sits {p8['med_dbl'] / thr[F_PROD]:.0f}x (r 1) and {p03['med_dbl'] / thr[F_PROD]:.1f}x (r 0.3) above the production threshold ({thr[F_PROD]:.0f}), "
             f"and the detected pairs recover the injected pair (dt 8 r 1: fitted dt bias {p3_8['ddt_bias']:+.2f} ns, RMS {p3_8['ddt_rms']:.2f} ns, median fitted A2/A1 {p3_8['ratio_med']:.2f}; r 0.3: {p3_8b['ddt_bias']:+.2f} / {p3_8b['ddt_rms']:.2f} ns, {p3_8b['ratio_med']:.2f}). The fake threshold is set by the template mismatch of real pulses (real-single 97th percentile {nm['real_q97']:.0f} vs {nm['tpl_q97']:.1f} for template + noise at sig_ped, {nm['tplx_q97']:.0f} with noise at the measured baseline std), so it is template-specific and must be re-tuned for any other bank or condition.")
    L.append(f"**Singles vs doublets, dt 8 ns.** r 1: doublets below thr {p8['below_k']}/{p8['n_d']}, singles above {100 * p8['above']:.2f} %, OVL {p8['ovl']:.4f}; r 0.3: doublets below {100 * p03['below']:.1f} %, singles above {100 * p03['above']:.2f} %, OVL {p03['ovl']:.4f}. "
             f"Median dchi2: doublets {p8['med_dbl']:.0f} (r 1), {p03['med_dbl']:.0f} (r 0.3); real singles {nm['real']:.1f}. At dt 4: r 1 {100 * P[(4, 1.0)]['below']:.0f} % of doublets below thr (OVL {P[(4, 1.0)]['ovl']:.2f}), r 0.3 {100 * P[(4, 0.3)]['below']:.0f} % (OVL {P[(4, 0.3)]['ovl']:.2f}); no doublet below the production threshold at dt >= 12.")
    L.append(f"**Null dchi2 (eval singles): median / 97th pct.** analytic kappa_c*chi2(2) {nm['analytic']:.2f} / {nm['ana_q97']:.1f} (kappa_c {nm['analytic_kappa']:.2f}: weighted gain of an extra pulse; the 2 % eps term inflates the model sigma on the pulse; "
             f"correlated + rounding noise); tpl_noise singles (no template mismatch) {nm['tpl']:.3f} / {nm['tpl_q97']:.1f} ({100 * nm['tpl_zero']:.0f} % exactly 0: the search returns a pair only if its grid gain beats the 1-pulse grid gain, so the median is below chi2(2) while the tail is heavier); "
             f"tpl_noise with noise matched to the measured baseline std ({nm['noise_f']:.1f}x sig_ped; the stored sig_ped is that much below the sample std of the pulse baselines) {nm['tplx']:.2f} / {nm['tplx_q97']:.0f}; "
             f"real singles {nm['real']:.1f} / {nm['real_q97']:.0f}. Noise-only dchi2 is {nm['tpl_q97']:.0f}-{nm['tplx_q97']:.0f} at the 97th pct, {nm['real_q97'] / nm['tplx_q97']:.0f}-{nm['real_q97'] / nm['tpl_q97']:.0f}x below the real singles: the real-pulse shape mismatch sets the fake level. "
             f"Same at 0.73 / 3.49 GeV (real median {null_medians(SS[OTHERS[0]])['real']:.1f} / {null_medians(SS[OTHERS[1]])['real']:.1f}, thr {SS[OTHERS[0]]['thr'][F_PROD]:.0f} / {SS[OTHERS[1]]['thr'][F_PROD]:.0f}; dt 8 r 1 doublets below thr {R['other'][(OTHERS[0], 1.0)]['below_k']} / {R['other'][(OTHERS[1], 1.0)]['below_k']} of {R['other'][(OTHERS[0], 1.0)]['n_d']}).")
    L.append(f"**Where the fakes live** (eval singles above thr 3 %: {pf['n']}; above 0.1 %: {pf1['n']}). Fitted secondary |dt| < 2 ns or A_sec/A_main < 0.05: {pf['n_degen']} ({100 * pf['n_degen'] / max(pf['n'], 1):.0f} %; 0.1 % set {pf1['n_degen']}/{pf1['n']}); "
             f"genuine-looking {pf['n_gen']}. {pf['after']} of {pf['n']} secondaries lie after the main pulse (median offset {pf['off_med']:+.0f} ns, median ratio {pf['rat_med']:.3f}): a few-percent ripple on the slow tail, i.e. shape mismatch, not a second pulse. "
             f"The 3 % fakes sit at low amplitude (median {pf['A_med']:.0f} ADC, all singles {pf['A_med_all']:.0f}; flat dchi2 cut, smaller pulses have larger relative mismatch); the rare > 0.1 % ones at high amplitude (median {pf1['A_med']:.0f} ADC).")
    def ef(dt, r, tag):
        x = roc[(dt, r)][tag]
        lo, hi = wil(x["k"], x["n"])
        return f"{100 * x['eff']:.0f} [{100 * lo:.0f}-{100 * hi:.0f}]"
    L.append(f"**Efficiency at tighter fake (real2, % [Wilson 68 %]; thr from cal half; eval fake {roc[(8, 1.0)]['p01']['fake_k']}/{nS} = {100 * roc[(8, 1.0)]['p01']['fake']:.2f} % and {roc[(8, 1.0)]['p002']['fake_k']}/{nS} = {100 * roc[(8, 1.0)]['p002']['fake']:.3f} %).** "
             f"0.1 %: dt 8 r 1 / 0.3 = {ef(8, 1.0, 'p01')} / {ef(8, 0.3, 'p01')}; dt 12 r 1 / 0.3 = {ef(12, 1.0, 'p01')} / {ef(12, 0.3, 'p01')}. "
             f"0.02 %: dt 8 r 1 / 0.3 = {ef(8, 1.0, 'p002')} / {ef(8, 0.3, 'p002')}; dt 12 r 1 / 0.3 = {ef(12, 1.0, 'p002')} / {ef(12, 0.3, 'p002')}. "
             f"At 3 % (production): dt 8 {100 * roc[(8, 1.0)]['p3']['eff']:.0f} / {100 * roc[(8, 0.3)]['p3']['eff']:.0f} %; MF {100 * roc[(8, 1.0)]['mf_eff']:.0f} / {100 * roc[(8, 0.3)]['mf_eff']:.0f} %.")
    L.append(f"Statistical limit: 0.02 % of {nS} eval singles = {0.0002 * nS:.1f} events (cal quantile {0.0002 * nC:.1f} events), so the 0.02 % thr ({thr[F_002]:.0f}) is the cal value ranked #{int(math.floor(F_002 * nC)) + 1} from the top (top three {', '.join(f'{v:.0f}' for v in cal_top)}) and +-1 event moves it a lot; "
             f"the efficiency drop between 0.1 % and 0.02 % at dt 8 is the real-single tail (eval max {ev_max:.2g}) reaching the doublet median ({p8['med_dbl']:.0f}). Read the 0.02 % numbers as indicative only.")
    L.append("**Robustness, dt 8 ns r 1 real2 (eff / eval fake at the production thr; thr retuned on the cal half of the same condition).**")
    L.append("| condition | prod thr: eff | fake | retuned thr | eff | fake |")
    L.append("|---|---|---|---|---|---|")
    for r in rob:
        L.append(f"| {r['name']} | {100 * r['eff_prod']:.1f} % | {100 * r['fake_prod']:.2f} % | {r['thr_tuned']:.0f} | {100 * r['eff_tuned']:.1f} % | {100 * r['fake_tuned']:.2f} % |")
    L.append(f"Noise x2: independent measured noise added to doublets and singles (total baseline std 2x the measured one, checked), fit sigma_ped doubled. +0.5 ns is the null check (fake {100 * p5_['fake_prod']:.1f} vs {100 * b_['fake_prod']:.1f} %). -0.5 ns cuts the first 0.5 ns of the template head and is NOT a pure shift: "
             f"fake {100 * m5_['fake_prod']:.0f} % (ramping the cut step away changes little), so dchi2 is sensitive to the first ns of the leading edge; the Wassim bank likewise needs a higher threshold. "
             f"Retuning gives fake {100 * min(r['fake_tuned'] for r in ret):.1f}-{100 * max(r['fake_tuned'] for r in ret):.1f} % at {100 * min(r['eff_tuned'] for r in ret):.0f}-{100 * max(r['eff_tuned'] for r in ret):.0f} % efficiency.")
    sin, sout = I["sgl_in_A"], I["sgl_out_A"]
    L.append(f"**Independence.** Pool split {I['split']} is one of the {I['n_inputs']} build inputs (templates npz meta; RUNS says only 'not in the production synth'), and all {I['n_inputs'] - 1} other splits are inputs too. "
             f"Rebuilt build pool (top {I['argv_n_max']} mf_energy/module; reproduces the hold-out list {I['holdout_reconstructed'][2]}/{I['holdout_reconstructed'][1]}): {I['n_split_in_build_pool']} pulses of this split, of which {I['n_in_build_pool']} of {I['n_pool_rows']} pool rows ({100 * I['n_in_build_pool'] / I['n_pool_rows']:.1f} %; eval {I['n_in_build_pool_eval']}) "
             f"were template-build inputs (before the build's own outlier rejection; median {I['A_in_med']:.0f} ADC vs {I['A_all_med']:.0f}). Their median dchi2 is lower (amplitude-matched 1500-2500 ADC: {sin['med']:.1f} vs {sout['med']:.1f}); tail fake rates are consistent within statistics (in-build {I['sgl_in']['k']}/{I['sgl_in']['n']} vs rest {I['sgl_out']['k']}/{I['sgl_out']['n']}); "
             f"headline numbers use all eval rows, {100 * I['n_in_build_pool_eval'] / nS:.1f} % of which were template-build inputs. Cal / eval are disjoint rows ({I['halves_disjoint']}); real2 bases are eval and partners are same module + same half + not the base ({I['partner_same_module']} / {I['partner_same_half']} / {I['partner_not_base']}); partners from the cal half: {I['partner_in_cal_count']}; thresholds use cal rows only.")
    fk = 100 * (S["Rb"][S["pool"]["parity"] == 1, 12] > thr[F_PROD]).mean()
    L.append(f"**Small ratios r = 0.1 / 0.05 (second pulse smaller; real2; fit eff % / MF eff %; fit fake {fk:.1f} %, MF fake {100 * S['mf_single_fake'][0] / S['mf_single_fake'][1]:.2f} %).**")
    L.append("| arm, r \\ dt (ns) | " + " | ".join(str(d) for d in SM_DT) + " |")
    L.append("|---|" + "---|" * len(SM_DT))
    for arm in ("real2", "real"):
        for r in SM_R:
            L.append(f"| {arm}, {r:g} | " + " | ".join(f"{100 * tab[(arm, r, dt)]['eff']:.0f} / {100 * tab[(arm, r, dt)]['mf']:.0f}" for dt in SM_DT) + " |")
    def tl(r, dt):
        t = tab[("real2", r, dt)]
        return f"{100 * t['true_det']:.0f} / {100 * t['mismatch']:.0f} / {100 * t['mim']:.2f}"
    pr = {r: float(np.median(S["pool"]["A_b"][S["cells"][("real2", 8, r)]["ev"]["partner"]] / S["cells"][("real2", 8, r)]["ev"]["As"])) for r in SM_R}
    bz = [b for b in R["small_bias"] if b["arm"] == "real2" and b.get("n", 0) >= 50 and b["dt"] >= 12]
    L.append(f"Fit detections that are the injected pulse (+-4 ns, A within x2) / fire but not on it / eval singles mimicking the pair, % (dt 8, 20): r 0.1 {tl(0.1, 8)}, {tl(0.1, 20)}; r 0.05 {tl(0.05, 8)}, {tl(0.05, 20)}. "
             f"So the detections are the injected pulse, not template-mismatch ripples; r 0.05 efficiency is limited (48 % at dt 8, 72 % at dt 20, 96 % at dt 32). Detected pairs, dt >= 12 (n >= 50): |dt bias| <= {max(abs(b['dt_bias']) for b in bz):.2f} ns, dt RMS <= {max(b['dt_rms'] for b in bz):.2f} ns, "
             f"A_small bias {100 * min(b['as_bias'] for b in bz):+.0f} to {100 * max(b['as_bias'] for b in bz):+.0f} %, RMS {100 * min(b['as_rms'] for b in bz):.0f}-{100 * max(b['as_rms'] for b in bz):.0f} %. Small pulse before / after the big one 50 / 50; the rescaled real2 small pulse is nearly noise-free (partner/target amplitude median {pr[0.1]:.0f} at r 0.1, {pr[0.05]:.0f} at r 0.05); "
             f"A_b clip rejects {S['cells'][('real2', 8, 0.1)]['ev']['n_reject']}/{S['cells'][('real2', 8, 0.1)]['ev']['n_drawn']}, base A median {np.median(S['cells'][('real2', 8, 0.1)]['ev']['Ab']):.0f} vs eval pool {np.median(S['pool']['A_b'][S['pool']['parity'] == 1]):.0f} ADC.")
    e = R["e"]
    L.append("**r = 0.05 pulse vs noise and template residual (median first pulse, eval singles).** A_s = 0.05 A_med = " + " / ".join(f"{0.05 * x['A_med']:.0f}" for x in e)
             + f" ADC at {' / '.join(x['rt'].split()[0] for x in e)} GeV: {min(x['ratio_sig_med'] for x in e):.0f}-{max(x['ratio_sig_med'] for x in e):.0f} sigma_ped, {min(x['ratio_rms_med'] for x in e):.1f}-{max(x['ratio_rms_med'] for x in e):.1f}x the 1-pulse residual RMS, "
             f"{min(x['ratio_max_med'] for x in e):.1f}-{max(x['ratio_max_med'] for x in e):.1f}x its max: far above the noise, comparable to the template-shape residual, so shape error limits r <= 0.1 at short dt.")
    (out / "overfit_diagnostics.md").write_text("\n".join(L) + "\n")
    log(f"md lines {len(L)}")


if __name__ == "__main__":
    main()
