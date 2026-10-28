#!/usr/bin/env python3
"""cfd_fraction_scan.py - CFD-fraction scan of the PbWO4 time resolution vs RF, compared with the matched filter.

PURPOSE
    The production WaveAnalyzer pickoff uses cfd_fraction = 0.5 (hard-coded).  Is another fraction better
    for the PbWO4 modules, and does any fraction beat the matched-filter (MF) time?  `dump` re-extracts a
    split with the same MF reconstruction as module_timing_dump.py and, for the dominant pulse of every
    PbWO4 channel with m_E >= 50 MeV, also evaluates the WaveAnalyzer pickoff at CFD fractions
    f = 0.2 .. 0.7 (numpy port of the C++ pickoff incl. the smoothing, vectorised over pulses).  `analyze` joins the dump to the
    RF TDC dump and fits the width sigma_t of the RF-folded phase per (module, method).  `validate`
    writes validation.txt (MF reproduction, f = 0.5 port vs C++, sanity vs the RF calibration).

DEFINITIONS
    MF time            m_t of module_timing_dump.py (mf, best-50 spline templates, nsigma 5, min_amp 5,
                       kmax 4, r_model 0.3; dominant = largest Eu), unchanged.  m_E = mod.energize(Eu).
    matched WA peak    the WaveAnalyzer peak of the same channel with the largest integral (first on ties);
                       wa_t = its C++ peak.time (ns), wa_dt = |wa_t - m_t|; wa_far = wa_dt > 8 ns.
    lin_f              linear CFD on the smoothed waveform (WaveAnalyzer.cpp findPeaks, steps 1-2):
                       fallback t = raw_pos + t_sub - 1.5 (quadratic vertex of raw[pos-1..pos+1]); if
                       raw_height > 5 ped_rms the LAST crossing v[j] < f*raw_height <= v[j+1] with
                       v[j+1]-v[j] > 1.5 ped_rms, j in [int_left-1, raw_pos-1], else raw_pos + t_sub.
    cfd_f              full production emulation (steps 1-3): if raw_height > 10 ped_rms, the 4-parameter
                       log-normal LM fit (A, t0, mu, sigma) of fit_log_normal_cfd is run (float32,
                       vectorised, identical bounds / lambda schedule / Cholesky / 60 iterations); if its
                       f-dependent gate (a rising-edge crossing of f*(buf[pos]-ped) from int_left-2) and
                       acceptance window pass, t = t0 + exp(mu - sigma*sqrt(-2 ln f)) replaces lin_f.
                       The fit does not depend on f, so it is run once per pulse and evaluated for all f.
                       cfd_code per f: 0 peak fallback, 1 linear CFD, 2 log-normal fit.
                       Times in ns = t_sample * 1e3 / 250.  float32 arithmetic; summation order and
                       libm differ from the C++ binary, so fit rows agree only to float32-LM precision.
    sigma_t            per (module, method): RF phase phi = fold(t - t_RF) (P = 4.008 ns, nearest channel-A
                       hit chosen with the MF time, same hit for all methods) in 80 bins over one period,
                       sigma = plot_module_timing_rf.fit_unwrapped width (as rf_time_calib.py).  Fallback
                       rows dropped.  All 13 methods use the SAME rows (all methods finite).  Module
                       offsets are irrelevant (the fit centres itself).
    bands              e100: m_E >= 100 MeV, N >= 300 rows;  e500: m_E >= 500 MeV, N >= 150 rows.
                       A module enters a band's median only if all 13 fits are ok.

OUTPUT  (default output/cfd_scan_20261027/)
    dump/<run:06d>_<split:05d>.npz   per event ev_num, trig; per MF row m_ev, m_idx, m_t, m_E, m_npk, m_inwin,
                                     m_fb (as module_timing_dump, no veto arrays); for rows with m_E >= 50
                                     and not m_fb: wa_t, wa_h, wa_ped_rms, wa_dt, wa_far, t_lin[n,6],
                                     t_cfd[n,6], cfd_code[n,6], fracs (NaN / -1 elsewhere); metadata.
    sigma_t_modules.tsv, sigma_t_medians.tsv / .md, cfd_scan.pdf, validation.txt

USAGE
    python3 cfd_fraction_scan.py dump [FILE ...] [--out DIR] [--max-events N] [--parallel 4]
        [--templates-root output/templates_spline_best50_20261026]     (no FILE: the 4 default splits)
    python3 cfd_fraction_scan.py analyze [--dump-dir DIR] [--rf-dir DIR] [--out-dir DIR] [--jobs 4]
    python3 cfd_fraction_scan.py validate [--out-dir DIR] [--ref-dump-dir DIR] [--ref-rf-calib TSV]
"""
import argparse
import csv
import glob
import hashlib
import json
import multiprocessing as mp
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import module_timing_dump as MTD  # noqa: E402

SCRIPT_VERSION = "cfd_fraction_scan 1.0 (20261027)"
OUT_DIR = "output/cfd_scan_20261027"
TEMPLATES_BEST50 = "output/templates_spline_best50_20261026"
REF_DUMP_DIR = "output/best50_20261026/module_timing/dump"
RF_DIR = "output/module_timing_20261025/rf"
REF_RF_CALIB = "output/best50_20261026/rf_calib/rf_calib_modules.tsv"
DEFAULT_SPLITS = ("prad_024590.evio.00010", "prad_025320.evio.00050",
                  "prad_024512.evio.00010", "prad_025800.evio.00100")
FRACS = (0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
E_CFD_MIN = 50.0            # MeV: CFD columns only for rows above this
FAR_NS = 8.0                # |wa_t - m_t| flag
CLK_MHZ = 250.0
QUEUE_ROWS = 20000          # CFD batch size in the dump
BANDS = (("e100", 100.0, 300), ("e500", 500.0, 150))
GEOMETRY_RUN = 25320
F32 = np.float32


# ----------------------------------------------------------------------------
# numpy port of the WaveAnalyzer time pickoff (float32)
# ----------------------------------------------------------------------------
def smooth_np(raw, order):
    """WaveAnalyzer::smooth (triangular kernel, float32) for raw (n,S) uint16.  The prad2py binding
    wa.smooth() segfaults (request() called with the GIL released), hence the numpy port."""
    raw = np.atleast_2d(raw)
    n, S = raw.shape
    if order <= 1:
        return raw.astype(F32)
    val = raw.astype(F32)
    wsum = np.ones(S, F32)
    for j in range(1, order):
        w = F32(1) - F32(j) / F32(order + 1)
        ok = np.zeros(S, bool)
        ok[j:S - j] = True                       # j <= i and i + j < S
        nb = np.zeros((n, S), F32)
        nb[:, j:S - j] = raw[:, :S - 2 * j].astype(np.int32) + raw[:, 2 * j:].astype(np.int32)
        val = val + np.where(ok, w * nb, F32(0))
        wsum = wsum + np.where(ok, F32(2) * w, F32(0))
    return (val / wsum).astype(F32)


def _crossing(v0, v1, thr, ok_j, dv_min=None):
    """Last j with v0[j] < thr <= v1[j] (and v1-v0 > dv_min) inside ok_j.  Returns (found, j, frac)."""
    c = (v0 < thr[:, None]) & (v1 >= thr[:, None]) & ok_j
    if dv_min is not None:
        c &= (v1 - v0) > dv_min[:, None]
    found = c.any(1)
    n = c.shape[1]
    j = n - 1 - np.argmax(c[:, ::-1], axis=1)
    r = np.arange(len(j))
    d = v1[r, j] - v0[r, j]
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = (thr - v0[r, j]) / d
    return found, j, frac


def _model(t, ped, p):
    """log_normal_pulse_value for t (n,W), ped (n,), p (n,4) float32."""
    A, t0, mu, sg = (p[:, k:k + 1] for k in range(4))
    ok = (sg > 0) & (t > t0)
    with np.errstate(all="ignore"):
        dt = np.where(ok, t - t0, F32(1))
        z = (np.log(dt) - mu) / np.where(sg > 0, sg, F32(1))
        val = ped[:, None] + A * np.exp(F32(-0.5) * z * z)
    return np.where(ok, val, ped[:, None]).astype(F32)


def _chol_solve4(Aug, g):
    """Batched 4x4 Cholesky solve (float32, as the C++ cholesky/chol_solve).  Returns (delta, fail)."""
    n = len(Aug)
    L = np.zeros((n, 4, 4), F32)
    fail = np.zeros(n, bool)
    with np.errstate(all="ignore"):
        for i in range(4):
            for j in range(i + 1):
                s = Aug[:, i, j].copy()
                for k in range(j):
                    s = s - L[:, i, k] * L[:, j, k]
                if i == j:
                    fail |= s <= 0
                    L[:, i, i] = np.sqrt(s)
                else:
                    L[:, i, j] = s / L[:, j, j]
        y = np.zeros((n, 4), F32)
        for i in range(4):
            s = g[:, i].copy()
            for k in range(i):
                s = s - L[:, i, k] * y[:, k]
            y[:, i] = s / L[:, i, i]
        x = np.zeros((n, 4), F32)
        for i in range(3, -1, -1):
            s = y[:, i].copy()
            for k in range(i + 1, 4):
                s = s - L[:, k, i] * x[:, k]
            x[:, i] = s / L[:, i, i]
    return x, fail


def lognormal_fit(buf, left, pos, ped, rms):
    """Vectorised fit_log_normal_cfd parameter search (everything except the f-dependent gate/acceptance).

    buf (n,S) float32 smoothed, left/pos int arrays, ped/rms float32 arrays.  Returns dict with p_best (n,4),
    fit_start, fit_stop, valid (fit ran: nfit >= 8 and chi2 finite)."""
    n, S = buf.shape
    fit_start = np.maximum(0, left - 4)
    fit_stop = np.minimum(S - 1, pos + 4)
    nfit = fit_stop - fit_start + 1
    r_ = np.arange(n)
    W = int(max(nfit.max(), 1)) if n else 1
    ar = np.arange(W)
    valid = ar[None, :] < nfit[:, None]
    idx = np.minimum(fit_start[:, None] + ar[None, :], S - 1)
    y = np.take_along_axis(buf, idx, 1)
    t = idx.astype(F32)
    peak_time = pos.astype(F32)
    A_guess = np.maximum(buf[r_, pos] - ped, F32(1))
    lo = np.stack([np.zeros(n, F32), np.maximum(F32(0), left.astype(F32) - F32(1.5)),
                   np.full(n, 0.5, F32), np.full(n, 0.1, F32)], 1)
    hi = np.stack([F32(3) * A_guess, peak_time - F32(1e-3), np.full(n, 2, F32), np.full(n, 0.9, F32)], 1)
    p = np.stack([A_guess, np.maximum(F32(0), left.astype(F32) - F32(0.5)),
                  np.ones(n, F32), np.full(n, 0.5, F32)], 1)
    p = np.minimum(np.maximum(p, lo), hi)
    sigma_fit = np.maximum(rms, F32(1))
    inv_s2 = (F32(1) / (sigma_fit * sigma_fit)).astype(F32)
    dof = np.maximum(1, nfit - 4).astype(F32)

    def chi2_of(pp, rows):
        bad = ~(pp[:, 0] > 0) | ~(pp[:, 3] > 0) | ~(pp[:, 1] < peak_time[rows])
        r = np.where(valid[rows], y[rows] - _model(t[rows], ped[rows], pp), F32(0))
        c = (r * r).sum(1, dtype=F32) * inv_s2[rows] / dof[rows]
        return np.where(bad, F32(np.inf), c).astype(F32)

    allr = np.arange(n)
    chi2 = chi2_of(p, allr)
    ran = (nfit >= 8) & np.isfinite(chi2)
    p_best = p.copy()
    chi2_best = chi2.copy()
    lam = np.full(n, 1e-3, F32)
    done = ~ran
    LAMBDA_MAX = F32(1e10)
    with np.errstate(all="ignore"):
        for _ in range(60):
            a = np.flatnonzero(~done)
            if not a.size:
                break
            pa, ta, ya, va, pedA = p[a], t[a], y[a], valid[a], ped[a]
            f0 = _model(ta, pedA, pa)
            res = np.where(va, ya - f0, F32(0))
            h = np.stack([
                np.maximum(F32(1e-3) * np.maximum(np.abs(pa[:, 0]), F32(1)), F32(1e-3)),
                np.full(len(a), 1e-3, F32),
                np.maximum(F32(1e-3) * np.maximum(np.abs(pa[:, 2]), F32(1)), F32(1e-4)),
                np.maximum(F32(1e-3) * np.maximum(np.abs(pa[:, 3]), F32(1)), F32(1e-4))], 1)
            J = np.zeros((len(a), 4, ta.shape[1]), F32)
            for j in range(4):
                pf = pa.copy()
                pf[:, j] = np.minimum(hi[a, j], pa[:, j] + h[:, j])
                step = pf[:, j] - pa[:, j]
                d = -(_model(ta, pedA, pf) - f0) / np.where(step > 0, step, F32(1))[:, None]
                J[:, j] = np.where(va & (step > 0)[:, None], d, F32(0))
            isa = inv_s2[a]
            g = np.einsum("njw,nw->nj", J, res).astype(F32) * isa[:, None]
            Am = np.einsum("njw,nkw->njk", J, J).astype(F32) * isa[:, None, None]
            Aug = Am + lam[a][:, None, None] * np.eye(4, dtype=F32)[None]
            delta, fail = _chol_solve4(Aug, g)
            # Cholesky failure: only lambda changes
            lf = lam[a] * F32(10)
            lam[a[fail]] = lf[fail]
            done[a[fail & (lf > LAMBDA_MAX)]] = True
            ok = ~fail
            if not ok.any():
                continue
            b = a[ok]
            pb, dl = pa[ok], delta[ok]
            pn = np.minimum(np.maximum(pb - dl, lo[b]), hi[b])
            scale = np.maximum(np.abs(pb), F32(1e-3))
            dpar = (np.abs(pn - pb) / scale).sum(1, dtype=F32)
            cn = chi2_of(pn, b)
            fin = np.isfinite(cn)
            better = fin & (cn < chi2_best[b])
            p_best[b[better]] = pn[better]
            chi2_best[b[better]] = cn[better]
            acc = fin & (cn < chi2[b])
            ba, br = b[acc], b[~acc]
            p[ba] = pn[acc]
            chi2[ba] = cn[acc]
            lam[ba] = np.maximum(lam[ba] / F32(10), F32(1e-12))
            done[ba[dpar[acc] < F32(1e-4)]] = True
            lam[br] = lam[br] * F32(10)
            done[br[lam[br] > LAMBDA_MAX]] = True
    return dict(p_best=p_best, fit_start=fit_start, fit_stop=fit_stop, valid=ran)


def pickoff_times(raw, buf, pos, left, height, ped, rms, fracs=FRACS, clk=CLK_MHZ):
    """Production time pickoff at several CFD fractions.

    raw (n,S) uint16, buf (n,S) float32 smoothed, pos/left int, height/ped/rms float32.
    Returns t_lin, t_cfd (n,F) ns float32 and code (n,F) int8 (0 fallback, 1 linear, 2 log-normal fit)."""
    n, S = buf.shape
    nf = len(fracs)
    r_ = np.arange(n)
    pos = pos.astype(np.int64)
    left = left.astype(np.int64)
    ped, rms, height = ped.astype(F32), rms.astype(F32), height.astype(F32)
    # step 1: quadratic sub-sample vertex
    rawf = raw.astype(F32)
    pm1 = rawf[r_, np.clip(pos - 1, 0, S - 1)]
    p0 = rawf[r_, pos]
    pp1 = rawf[r_, np.clip(pos + 1, 0, S - 1)]
    denom = pm1 - F32(2) * p0 + pp1
    with np.errstate(all="ignore"):
        dl = F32(0.5) * (pm1 - pp1) / denom
    tsub = np.where((pos > 0) & (pos < S - 1) & (denom < F32(-1e-3)),
                    np.clip(dl, F32(-1), F32(1)), F32(0)).astype(F32)
    t_fb = pos.astype(F32) + tsub - F32(1.5)
    v = buf - ped[:, None]
    v0, v1 = v[:, :-1], v[:, 1:]
    j = np.arange(S - 1)[None, :]
    ok_lin = (j >= np.maximum(left - 1, 0)[:, None]) & (j < pos[:, None])
    ok_gate = (j >= np.maximum(left - 2, 0)[:, None]) & (j < pos[:, None])
    h5 = height > F32(5) * rms
    h10 = height > F32(10) * rms
    hb = buf[r_, pos] - ped                        # height used by the fit gate
    rr = np.flatnonzero(h10)
    sub = lognormal_fit(buf[rr], left[rr], pos[rr], ped[rr], rms[rr]) if rr.size else None
    t_lin = np.zeros((n, nf), F32)
    t_cfd = np.zeros((n, nf), F32)
    code = np.zeros((n, nf), np.int8)
    for k, f in enumerate(fracs):
        ff = F32(f)
        found, jj, frac = _crossing(v0, v1, ff * height, ok_lin, F32(1.5) * rms)
        use = h5 & found
        tl = np.where(h5, pos.astype(F32) + tsub, t_fb)
        tl = np.where(use, jj.astype(F32) + np.clip(frac, F32(0), F32(1)), tl).astype(F32)
        cd = np.where(use, np.int8(1), np.int8(0))
        tc = tl.copy()
        if sub is not None:
            gate_r = _crossing(v0[rr], v1[rr], ff * hb[rr], ok_gate[rr])[0]
            p = sub["p_best"]
            root = np.sqrt(F32(-2) * np.log(ff)).astype(F32)
            with np.errstate(all="ignore"):
                tcs = (p[:, 1] + np.exp(p[:, 2] - p[:, 3] * root)).astype(F32)
                tpk = (p[:, 1] + np.exp(p[:, 2])).astype(F32)
            fs = sub["fit_start"].astype(F32)
            acc = (sub["valid"] & gate_r & np.isfinite(tcs) & np.isfinite(tpk)
                   & (tcs >= fs) & (tcs <= pos[rr].astype(F32))
                   & (tpk >= fs) & (tpk <= sub["fit_stop"].astype(F32) + F32(1)))
            ri = rr[acc]
            tc[ri] = tcs[acc]
            cd[ri] = 2
        t_lin[:, k] = tl * F32(1e3) / F32(clk)
        t_cfd[:, k] = tc * F32(1e3) / F32(clk)
        code[:, k] = cd
    return t_lin.astype(F32), t_cfd.astype(F32), code


# ----------------------------------------------------------------------------
# dump: one split
# ----------------------------------------------------------------------------
def process_split(path, out_dir, max_events=0, templates_root=None):
    import _common as C
    import mf_pulse
    from prad2py import dec

    run, out_path = MTD.out_name(path, out_dir)
    if run not in MTD.RUN_CFG:
        sys.exit(f"[ERROR] run {run} not in RUN_CFG")
    tmpl = os.path.join(templates_root or os.path.join(HERE, TEMPLATES_BEST50), MTD.RUN_CFG[run][0])
    calib = os.path.join(MTD.CALIB_DIR, MTD.RUN_CFG[run][1])
    t_start = time.monotonic()

    p = C.setup_pipeline(evio_path=path, hc_calib_file=calib)
    hy = [p.hycal.module(i) for i in range(p.hycal.module_count())]
    names = np.array([m.name for m in hy])
    pwo4 = [bool(m.is_pwo4()) for m in hy]
    wcfg = p.wave_ana.cfg
    overflow = int(wcfg.overflow)
    assert float(wcfg.clk_mhz) == CLK_MHZ, wcfg.clk_mhz
    bank, info = mf_pulse.build_bank(
        [m.name for m in hy], pwo4, tmpl,
        dict(smooth_order=int(wcfg.smooth_order),
             int_tail_ratio=float(wcfg.int_tail_ratio),
             tail_break_n=int(wcfg.tail_break_n)),
        column=mf_pulse.resolve_column(tmpl, None, prefer="25320"),
        int_tail_ratio=None)
    rec = mf_pulse.MFReconstructor(bank, nsigma=5.0, min_amp=5.0, kmax=4,
                                   r_model=0.3, ped_nsamples=int(wcfg.ped_nsamples))
    print(f"[setup] run {run} templates {tmpl} ({info['kind']}) calib {calib}", flush=True)

    ev_num, trig = [], []
    M = {k: [] for k in ("ev", "idx", "t", "E", "npk", "inwin", "fb")}
    M_DT = dict(ev=np.int32, idx=np.int16, t=np.float32, E=np.float32,
                npk=np.int8, inwin=bool, fb=bool)
    M_chunks = []
    nrow = [0]                       # running row number in M (arrival order)
    # CFD queue: per row WaveAnalyzer quantities
    Q = dict(row=[], raw=[], pos=[], left=[], h=[], ped=[], rms=[], wat=[], mt=[])
    R = dict(row=[], wa_t=[], wa_h=[], wa_rms=[], wa_dt=[], t_lin=[], t_cfd=[], code=[])

    def pack_m():
        if M["ev"]:
            M_chunks.append({k: np.asarray(v, dtype=M_DT[k]) for k, v in M.items()})
            for v in M.values():
                v.clear()

    def put_m(ev, idx, t, E, npk, fb):
        M["ev"].append(ev); M["idx"].append(idx); M["t"].append(t)
        M["E"].append(E); M["npk"].append(npk)
        M["inwin"].append(MTD.PRE_LO < t < MTD.PRE_HI); M["fb"].append(fb)
        nrow[0] += 1

    def put_wa_dominant(ev, mod, peaks):
        if not peaks:
            return
        best = peaks[0]
        for pk in peaks[1:]:
            if pk.integral > best.integral:
                best = pk
        put_m(ev, mod.index, best.time, mod.energize(best.integral), len(peaks), True)

    def queue_row(row, samples, mt):
        ped, rms, pks = p.wave_ana.analyze(samples)
        pks = list(pks)
        raw = np.asarray(samples, dtype=np.uint16)
        R_ = dict(row=row, wa_t=np.nan, wa_h=np.nan, wa_rms=float(rms), wa_dt=np.nan)
        if not pks:
            R["row"].append(row); R["wa_t"].append(np.nan); R["wa_h"].append(np.nan)
            R["wa_rms"].append(float(rms)); R["wa_dt"].append(np.nan)
            R["t_lin"].append(np.full(len(FRACS), np.nan, np.float32))
            R["t_cfd"].append(np.full(len(FRACS), np.nan, np.float32))
            R["code"].append(np.full(len(FRACS), -1, np.int8))
            return
        best = pks[0]
        for pk in pks[1:]:
            if pk.integral > best.integral:
                best = pk
        Q["row"].append(row); Q["raw"].append(raw.copy())
        Q["pos"].append(int(best.pos)); Q["left"].append(int(best.left))
        Q["h"].append(float(best.height)); Q["ped"].append(float(ped)); Q["rms"].append(float(rms))
        Q["wat"].append(float(best.time)); Q["mt"].append(float(mt))

    def run_queue():
        if not Q["row"]:
            return
        raw = np.stack(Q["raw"]); buf = smooth_np(raw, int(wcfg.smooth_order))
        pos = np.asarray(Q["pos"]); left = np.asarray(Q["left"])
        h = np.asarray(Q["h"], F32); ped = np.asarray(Q["ped"], F32); rms = np.asarray(Q["rms"], F32)
        tl, tc, cd = pickoff_times(raw, buf, pos, left, h, ped, rms)
        wat = np.asarray(Q["wat"], np.float32)
        R["row"].extend(Q["row"]); R["wa_t"].extend(wat); R["wa_h"].extend(h)
        R["wa_rms"].extend(rms); R["wa_dt"].extend(np.abs(wat.astype(np.float64) - np.asarray(Q["mt"])))
        R["t_lin"].extend(tl); R["t_cfd"].extend(tc); R["code"].extend(cd)
        for v in Q.values():
            v.clear()

    st = dict(Y=None, mid=None, n=0, ns=None)
    row_ev = []

    def flush():
        n = st["n"]
        if n == 0:
            return
        Y, mid = st["Y"][:n], st["mid"][:n]
        res = rec.run(Y, mid, overflow)
        pc = res.p_ch
        Eu = mf_pulse.pulse_energy_units(bank, mid[pc], res.p_t0, res.p_amp, 1.0, None)
        tl = res.p_time
        starts = np.searchsorted(pc, np.arange(n + 1))
        order = np.lexsort((-Eu, pc))
        first = order[starts[:-1][starts[1:] > starts[:-1]]] if len(pc) else []
        for j in first:
            r = int(pc[j])
            mod = hy[int(mid[r])]
            E = mod.energize(float(Eu[j]))
            row = nrow[0]
            put_m(row_ev[r], mod.index, float(tl[j]), E, int(starts[r + 1] - starts[r]), False)
            if E >= E_CFD_MIN:
                queue_row(row, Y[r], float(tl[j]))
        for r in np.nonzero(res.wa_fb)[0]:
            _, _, pks = p.wave_ana.analyze(Y[r])
            put_wa_dominant(row_ev[r], hy[int(mid[r])], list(pks))
        st["n"] = 0
        row_ev.clear()
        pack_m()
        if len(Q["row"]) >= QUEUE_ROWS:
            run_queue()

    ch = dec.EvChannel()
    ch.set_config(p.cfg)
    n_read = n_phys = n_evt = 0
    if ch.open_auto(path) != dec.Status.success:
        sys.exit(f"[ERROR] cannot open {path}")
    done = False
    while not done and ch.read() == dec.Status.success:
        n_read += 1
        if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
            continue
        for i in range(ch.get_n_events()):
            d = ch.decode_event(i, with_ssp=False)
            if not d["ok"]:
                continue
            e = d["event"]
            n_phys += 1
            ev = n_evt
            n_evt += 1
            ev_num.append(int(e.info.event_number))
            trig.append(int(e.info.trigger_bits))
            for ri in range(e.nrocs):
                roc = e.roc(ri)
                if not roc.present:
                    continue
                crate = p.crate_map.get(roc.tag)
                if crate is None:
                    continue
                for s in roc.present_slots():
                    slot = roc.slot(s)
                    for c in slot.present_channels():
                        mod = p.hycal.module_by_daq(crate, s, c)
                        if mod is None:
                            continue
                        cd = slot.channel(c)
                        if cd.nsamples <= 0:
                            continue
                        if not (mod.is_hycal() and pwo4[mod.index]):
                            continue
                        if st["ns"] is None:
                            st["ns"] = cd.nsamples
                            cap = MTD.BLOCK_EVENTS * 256
                            st["Y"] = np.zeros((cap, st["ns"]), dtype=np.uint16)
                            st["mid"] = np.zeros(cap, dtype=np.int64)
                        if cd.nsamples != st["ns"]:
                            _, _, pks = p.wave_ana.analyze(cd.samples)
                            put_wa_dominant(ev, mod, list(pks))
                            continue
                        r = st["n"]
                        if r >= len(st["mid"]):
                            st["Y"] = np.concatenate([st["Y"], np.zeros_like(st["Y"])])
                            st["mid"] = np.concatenate([st["mid"], np.zeros_like(st["mid"])])
                        st["Y"][r] = cd.samples
                        st["mid"][r] = mod.index
                        row_ev.append(ev)
                        st["n"] = r + 1
            if n_evt % MTD.BLOCK_EVENTS == 0:
                flush()
            if max_events > 0 and n_phys >= max_events:
                done = True
                break
    flush()
    pack_m()
    run_queue()
    elapsed = time.monotonic() - t_start

    Mc = {k: (np.concatenate([c[k] for c in M_chunks]) if M_chunks
              else np.zeros(0, dtype=M_DT[k])) for k in M_DT}
    nr = len(Mc["ev"])
    assert nr == nrow[0]
    o = np.argsort(Mc["ev"], kind="stable")
    nF = len(FRACS)
    wa_t = np.full(nr, np.nan, np.float32); wa_h = wa_t.copy(); wa_rms = wa_t.copy(); wa_dt = wa_t.copy()
    t_lin = np.full((nr, nF), np.nan, np.float32); t_cfd = t_lin.copy()
    code = np.full((nr, nF), -1, np.int8)
    if R["row"]:
        rr = np.asarray(R["row"])
        wa_t[rr] = R["wa_t"]; wa_h[rr] = R["wa_h"]; wa_rms[rr] = R["wa_rms"]; wa_dt[rr] = R["wa_dt"]
        t_lin[rr] = np.stack(R["t_lin"]); t_cfd[rr] = np.stack(R["t_cfd"]); code[rr] = np.stack(R["code"])
    out = dict(
        ev_num=np.asarray(ev_num, dtype=np.int64), trig=np.asarray(trig, dtype=np.uint32),
        m_ev=Mc["ev"][o], m_idx=Mc["idx"][o], m_t=Mc["t"][o], m_E=Mc["E"][o], m_npk=Mc["npk"][o],
        m_inwin=Mc["inwin"][o], m_fb=Mc["fb"][o],
        wa_t=wa_t[o], wa_h=wa_h[o], wa_ped_rms=wa_rms[o], wa_dt=wa_dt[o],
        wa_far=(wa_dt[o] > FAR_NS), t_lin=t_lin[o], t_cfd=t_cfd[o], cfd_code=code[o],
        fracs=np.asarray(FRACS, np.float32),
        run=run, split=MTD.split_of(path), module_names=names,
        module_is_pwo4=np.asarray(pwo4, dtype=bool),
        templates=tmpl, calib_file=calib, n_read=n_read, n_phys=n_phys,
        elapsed_s=elapsed, script_version=SCRIPT_VERSION)
    os.makedirs(out_dir, exist_ok=True)
    tmp = out_path + ".tmp"
    with open(tmp, "wb") as fh:
        np.savez(fh, **out)
    os.replace(tmp, out_path)
    print(f"[done] run {run} split {MTD.split_of(path):05d}: {n_read} read, {n_phys} physics, "
          f"{nr} module rows, {len(R['row'])} CFD rows, {elapsed:.1f} s -> {out_path}", flush=True)


def cmd_dump(a):
    out = os.path.abspath(a.out)
    troot = os.path.abspath(a.templates_root)
    files = list(a.files) or [os.path.join(MTD.EVIO_DIR, f) for f in DEFAULT_SPLITS]
    todo = [f for f in files if not os.path.exists(MTD.out_name(f, out)[1])]
    print(f"[driver] {len(files)} split(s), {len(files) - len(todo)} done, {len(todo)} to run", flush=True)
    if a.parallel <= 1 or len(todo) <= 1:
        for f in todo:
            process_split(f, out, a.max_events, troot)
        return
    log_dir = os.path.join(out, "logs")
    os.makedirs(log_dir, exist_ok=True)
    running, queue, n_fail = [], list(todo), 0
    while queue or running:
        while queue and len(running) < a.parallel:
            f = queue.pop(0)
            op = MTD.out_name(f, out)[1]
            log = open(os.path.join(log_dir, os.path.basename(op)[:-4] + ".log"), "w")
            cmd = [sys.executable, os.path.abspath(__file__), "dump", f, "--out", out,
                   "--max-events", str(a.max_events), "--templates-root", troot, "--parallel", "1"]
            running.append((subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT), f, log))
        time.sleep(2.0)
        for item in list(running):
            pr, f, log = item
            if pr.poll() is None:
                continue
            running.remove(item)
            log.close()
            n_fail += pr.returncode != 0
            with open(log.name) as fh:
                last = [ln for ln in fh.read().splitlines() if ln.startswith("[done]")]
            print(last[-1] if pr.returncode == 0 and last else f"[FAIL rc={pr.returncode}] {f}", flush=True)
    print(f"[driver] finished, {n_fail} failed", flush=True)
    sys.exit(1 if n_fail else 0)


# ----------------------------------------------------------------------------
# analyze: sigma_t per (split, module, method)
# ----------------------------------------------------------------------------
def method_names():
    return ["MF"] + [f"lin_{f}" for f in FRACS] + [f"cfd_{f}" for f in FRACS]


def analyze_split(job):
    import plot_module_timing_rf as RF
    dp, rp, cfg = job
    o = RF.prep(dp, rp, cfg)
    d, ri, pm, P, t_rf = o["d"], o["ri"], o["pm"], o["P"], o["t_rf"]
    E = d["m_E"][ri].astype(np.float64)
    keep = ~d["m_fb"][ri]
    T = [o["t"]] + [d["t_lin"][ri, k].astype(np.float64) for k in range(len(FRACS))] \
        + [d["t_cfd"][ri, k].astype(np.float64) for k in range(len(FRACS))]
    fin = keep.copy()
    for t in T:
        fin &= np.isfinite(t)
    npw = cfg["npw"]
    res = dict(run=o["run"], split=o["split"], P=P, n_rows=len(ri), n_paired=int((fin & (E >= 100)).sum()),
               bands={})
    for band, emin, nmin in BANDS:
        k = fin & (E >= emin)
        N = np.bincount(pm[k], minlength=npw)
        mods = np.flatnonzero(N >= nmin)
        H = []
        for t in T:
            b = RF._bin(RF.fold(t[k] - t_rf[k], P), P)
            H.append(np.bincount(pm[k] * RF.NB + b, minlength=npw * RF.NB).reshape(npw, RF.NB))
        sig = np.full((len(mods), len(T)), np.nan)
        for i, m in enumerate(mods):
            for q in range(len(T)):
                r = RF.fit_unwrapped((H[q][m], P))
                if r["ok"]:
                    sig[i, q] = r["sigma"]
        res["bands"][band] = dict(mods=mods, N=N[mods], sigma=sig)
    return res


def mod_geometry():
    import plot_module_timing as PMT
    geo = PMT.load_geometry(GEOMETRY_RUN)
    pidx = np.flatnonzero(geo["pwo4"])
    pmap = np.full(len(geo["name"]), -1, np.int64)
    pmap[pidx] = np.arange(len(pidx))
    return geo, pidx, pmap


def cmd_analyze(a):
    import plot_module_timing as PMT
    out_dir = Path(a.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    geo, pidx, pmap = mod_geometry()
    names_pw = np.asarray(geo["name"])[pidx]
    files = sorted(glob.glob(os.path.join(a.dump_dir, "*.npz")))
    jobs = []
    for fp in files:
        mt = re.search(PMT.RUN_NAME, fp)
        if not mt or int(mt.group(1)) not in PMT.TYPE_OF_RUN:
            continue
        rp = os.path.join(a.rf_dir, os.path.basename(fp))
        d0 = np.load(fp)
        nm = d0["module_names"]
        if [str(n) for n in nm] != [str(n) for n in geo["name"]]:
            sys.exit("module_names of the dump differ from the geometry")
        cfg = dict(mask=0xF00, min_e=0.0, npw=len(pidx), pmap=pmap,
                   names_hash=hashlib.md5("\n".join(map(str, nm)).encode()).hexdigest())
        jobs.append((fp, rp, cfg))
    if not jobs:
        sys.exit(f"no dumps in {a.dump_dir}")
    t0 = time.time()
    if a.jobs <= 1 or len(jobs) == 1:
        R = [analyze_split(j) for j in jobs]
    else:
        with mp.get_context("spawn").Pool(min(a.jobs, 4, len(jobs))) as pool:
            R = pool.map(analyze_split, jobs, chunksize=1)
    print(f"[analyze] {len(R)} splits fitted in {time.time() - t0:.0f} s", flush=True)
    R.sort(key=lambda r: (PMT.TAGS.index(PMT.TYPE_OF_RUN[r["run"]]), r["run"], r["split"]))
    meth = method_names()
    # modules table
    with open(out_dir / "sigma_t_modules.tsv", "w") as fh:
        fh.write("split\trun_type\tmodule\tband\tN\tsigma_MF\t"
                 + "\t".join(f"sigma_lin_{f}" for f in FRACS) + "\t"
                 + "\t".join(f"sigma_cfd_{f}" for f in FRACS) + "\tall_ok\n")
        for r in R:
            tag = PMT.TYPE_OF_RUN[r["run"]]
            for band, _, _ in BANDS:
                B = r["bands"][band]
                for i, m in enumerate(B["mods"]):
                    s = B["sigma"][i]
                    fh.write(f"{r['run']:06d}_{r['split']:05d}\t{tag}\t{names_pw[m]}\t{band}\t{int(B['N'][i])}\t"
                             + "\t".join("nan" if not np.isfinite(x) else f"{x:.5f}" for x in s)
                             + f"\t{int(np.isfinite(s).all())}\n")
    # medians
    rows = []
    for r in R:
        tag = PMT.TYPE_OF_RUN[r["run"]]
        for band, _, _ in BANDS:
            B = r["bands"][band]
            s = B["sigma"]
            okm = np.isfinite(s).all(1) if len(s) else np.zeros(0, bool)
            s = s[okm]
            n = len(s)
            c05 = meth.index("cfd_0.5")
            for q, mname in enumerate(meth):
                f = float(mname.split("_")[1]) if "_" in mname else np.nan
                if n:
                    x = s[:, q]
                    rows.append(dict(
                        split=f"{r['run']:06d}_{r['split']:05d}", run_type=tag, band=band, n_mod=n,
                        method=mname.split("_")[0], f=f, median=np.median(x),
                        p16=np.percentile(x, 16), p84=np.percentile(x, 84),
                        ratio_cfd05=np.median(x / s[:, c05]), ratio_MF=np.median(x / s[:, 0])))
                else:
                    rows.append(dict(split=f"{r['run']:06d}_{r['split']:05d}", run_type=tag, band=band,
                                     n_mod=0, method=mname.split("_")[0], f=f, median=np.nan, p16=np.nan,
                                     p84=np.nan, ratio_cfd05=np.nan, ratio_MF=np.nan))
    cols = ("split", "run_type", "band", "n_mod", "method", "f", "median", "p16", "p84", "ratio_cfd05", "ratio_MF")
    with open(out_dir / "sigma_t_medians.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for x in rows:
            fh.write("\t".join(x[c] if isinstance(x[c], str) else (str(x[c]) if isinstance(x[c], int)
                                                                   else f"{x[c]:.5f}") for c in cols) + "\n")
    write_medians_md(out_dir / "sigma_t_medians.md", rows, R)
    write_pdf(out_dir / "cfd_scan.pdf", rows, R)
    print(f"[analyze] wrote tables and pdf in {out_dir}", flush=True)


def _get(rows, split, band, method, f=None):
    for x in rows:
        if x["split"] == split and x["band"] == band and x["method"] == method and \
                (f is None or abs(x["f"] - f) < 1e-9):
            return x
    return None


def write_medians_md(path, rows, R):
    L = ["# Median sigma_t per split and band (ns, vs RF)", "",
         "Bands: e100 = m_E >= 100 MeV (N >= 300 rows per module); e500 = m_E >= 500 MeV (N >= 150).  "
         "A module enters only if all 13 methods fit ok (same rows for all methods).  Ratios are the median of "
         "per-module ratios.", ""]
    import plot_module_timing as PMT
    for r in R:
        sp = f"{r['run']:06d}_{r['split']:05d}"
        tag = PMT.TYPE_OF_RUN[r["run"]]
        for band, emin, nmin in BANDS:
            mf = _get(rows, sp, band, "MF")
            L += [f"## {sp} ({PMT.LABEL[tag]}), {band}: E >= {emin:.0f} MeV, N >= {nmin}, {mf['n_mod']} modules, "
                  f"MF median {mf['median']:.4f} ns", "",
                  "| f | lin_f | cfd_f | cfd_f/cfd_0.5 | cfd_f/MF | lin_f/MF |", "|---|---|---|---|---|---|"]
            for f in FRACS:
                c, l = _get(rows, sp, band, "cfd", f), _get(rows, sp, band, "lin", f)
                L.append(f"| {f} | {l['median']:.4f} | {c['median']:.4f} | {c['ratio_cfd05']:.3f} | "
                         f"{c['ratio_MF']:.3f} | {l['ratio_MF']:.3f} |")
            L.append("")
    path.write_text("\n".join(L) + "\n")


def write_pdf(path, rows, R):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    import plot_module_timing as PMT
    fr = np.asarray(FRACS)
    with PdfPages(path) as pdf:
        for r in R:
            sp = f"{r['run']:06d}_{r['split']:05d}"
            fig, axs = plt.subplots(1, 2, figsize=(11, 4.6), sharey=False)
            for ax, (band, emin, nmin) in zip(axs, BANDS):
                mf = _get(rows, sp, band, "MF")
                if mf["n_mod"] == 0:
                    ax.text(0.5, 0.5, f"no module with all fits ok\n(N >= {nmin})", ha="center",
                            transform=ax.transAxes)
                else:
                    ax.axhspan(mf["p16"], mf["p84"], color="0.8", alpha=0.6, lw=0,
                               label="MF 16-84 %")
                    ax.axhline(mf["median"], color="0.3", lw=1.5, label=f"MF median {mf['median']:.3f} ns")
                    for meth, col, mk in (("cfd", "tab:blue", "o"), ("lin", "tab:orange", "s")):
                        c = [_get(rows, sp, band, meth, f) for f in FRACS]
                        med = np.array([x["median"] for x in c])
                        ax.fill_between(fr, [x["p16"] for x in c], [x["p84"] for x in c], color=col,
                                        alpha=0.18, lw=0)
                        ax.plot(fr, med, marker=mk, color=col,
                                label=f"{meth}_f median" + (" (production emulation)" if meth == "cfd" else
                                                            " (linear CFD only)"))
                ax.set_xlabel("CFD fraction f")
                ax.set_ylabel("sigma_t vs RF (ns)")
                ax.set_title(f"{band}: E >= {emin:.0f} MeV, N >= {nmin}, {mf['n_mod']} modules", fontsize=10)
                ax.grid(alpha=0.3)
                ax.legend(fontsize=7)
            fig.suptitle(f"{PMT.LABEL[PMT.TYPE_OF_RUN[r['run']]]} - split {sp}: sigma_t vs CFD fraction "
                         f"(median, 16/84 % over modules)", fontsize=11)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)


# ----------------------------------------------------------------------------
# validate
# ----------------------------------------------------------------------------
def cmd_validate(a):
    out_dir = Path(a.out_dir)
    L = []
    pr = lambda s="": L.append(s)
    # (a) MF reproduction
    pr("(a) MF reproduction: new dump vs existing best-50 module_timing dump (m_t, m_E, m_idx, m_fb)")
    for fp in sorted(glob.glob(str(out_dir / "dump" / "*.npz"))):
        ref = os.path.join(a.ref_dump_dir, os.path.basename(fp))
        d = np.load(fp)
        if not os.path.exists(ref):
            pr(f"  {os.path.basename(fp)}: no reference dump")
            continue
        r = np.load(ref)
        sel = r["m_ev"] < len(d["ev_num"])           # --max-events runs: compare the common prefix
        same_shape = int(sel.sum()) == len(d["m_ev"])
        s = f"  {os.path.basename(fp)}: rows new {len(d['m_ev'])} ref {int(sel.sum())}"
        if same_shape:
            dt = np.abs(d["m_t"].astype(np.float64) - r["m_t"][sel])
            de = np.abs(d["m_E"].astype(np.float64) - r["m_E"][sel])
            s += (f"; max|dm_t| {np.nanmax(dt):.3g} ns, max|dm_E| {np.nanmax(de):.3g} MeV, "
                  f"m_ev/m_idx/m_fb equal {np.array_equal(d['m_ev'], r['m_ev'][sel])}/"
                  f"{np.array_equal(d['m_idx'], r['m_idx'][sel])}/{np.array_equal(d['m_fb'], r['m_fb'][sel])}; "
                  f"ev_num/trig equal {np.array_equal(d['ev_num'], r['ev_num'][:len(d['ev_num'])])}/"
                  f"{np.array_equal(d['trig'], r['trig'][:len(d['trig'])])}")
        else:
            s += "; ROW COUNT DIFFERS"
        pr(s)
    # (b) f = 0.5 port vs C++
    pr("")
    pr("(b) f = 0.5: t_cfd[0.5] - wa_t (C++ peak.time) on 025320_00050, rows m_E >= 100 MeV, not m_fb, "
       "matched WA peak within 8 ns")
    fp = out_dir / "dump" / "025320_00050.npz"
    if fp.exists():
        d = np.load(fp)
        k5 = list(np.round(d["fracs"], 3)).index(0.5)
        E = d["m_E"]
        base = (E >= 100) & ~d["m_fb"] & np.isfinite(d["wa_t"]) & ~d["wa_far"]
        pr(f"  rows E>=100 and not fb: {int(((E >= 100) & ~d['m_fb']).sum())}; with matched WA peak <= 8 ns: "
           f"{int(base.sum())}; wa_far rows: {int((((E >= 100) & ~d['m_fb']) & d['wa_far']).sum())}")
        res = d["t_cfd"][:, k5].astype(np.float64) - d["wa_t"]
        code = d["cfd_code"][:, k5]
        for lab, m in (("all", base), ("code 0 fallback", base & (code == 0)),
                       ("code 1 linear", base & (code == 1)), ("code 2 fit", base & (code == 2))):
            x = res[m]
            if x.size == 0:
                pr(f"  {lab:16s} N 0")
                continue
            pr(f"  {lab:16s} N {x.size:7d}  RMS {np.sqrt(np.mean(x ** 2)):.3e} ns  max|d| {np.abs(x).max():.3e} ns  "
               f"frac|d|<1e-3 ns {np.mean(np.abs(x) < 1e-3):.5f}  frac|d|<1e-2 {np.mean(np.abs(x) < 1e-2):.5f}")
        m = base & (code == 2)
        x = np.abs(res[m])
        if x.size:
            pr(f"  code 2 |d| percentiles 50/90/99/99.9: " + " ".join(f"{np.percentile(x, q):.2e}" for q in (50, 90, 99, 99.9)))
        # also the lin rows and the share of codes
        pr("  code shares at f=0.5 (all rows with CFD): " + ", ".join(
            f"{c}: {np.mean(code[code >= 0] == c):.4f}" for c in (0, 1, 2)))
        pr("  lin_0.5 vs wa_t on rows with code 1 (linear step reached, fit not accepted): "
           f"max|d| {np.abs(d['t_lin'][:, k5].astype(np.float64) - d['wa_t'])[base & (code == 1)].max():.3e} ns"
           if (base & (code == 1)).any() else "  no code-1 rows")
    else:
        pr("  dump 025320_00050 missing")
    # (c) sanity vs rf_calib
    pr("")
    pr("(c) sanity: MF e100 sigma of this split vs rf_calib_modules.tsv (run-level, tier e100, all splits) on the "
       "same modules")
    tsv = out_dir / "sigma_t_modules.tsv"
    if tsv.exists() and os.path.exists(a.ref_rf_calib):
        ref = {}
        with open(a.ref_rf_calib) as fh:
            for rw in csv.DictReader(fh, delimiter="\t"):
                if rw["tier"] == "e100" and rw["sigma_ns"] not in ("nan", ""):
                    ref[(int(rw["run"]), rw["module"])] = float(rw["sigma_ns"])
        mine = {}
        with open(tsv) as fh:
            for rw in csv.DictReader(fh, delimiter="\t"):
                if rw["band"] == "e100" and rw["sigma_MF"] != "nan":
                    mine.setdefault(int(rw["split"][:6]), {})[rw["module"]] = (float(rw["sigma_MF"]), int(rw["N"]))
        for run, mm in sorted(mine.items()):
            common = [m for m in mm if (run, m) in ref]
            if not common:
                pr(f"  run {run}: no common modules")
                continue
            a_ = np.array([mm[m][0] for m in common])
            b_ = np.array([ref[(run, m)] for m in common])
            allrun = np.array([v for (r_, m), v in ref.items() if r_ == run])
            pr(f"  run {run}: split MF median {np.median(a_):.4f} vs rf_calib {np.median(b_):.4f} ns on {len(common)} "
               f"common modules (median ratio {np.median(a_ / b_):.3f}, 16/84 % of ratio "
               f"{np.percentile(a_ / b_, 16):.3f}/{np.percentile(a_ / b_, 84):.3f}); rf_calib all-e100 median "
               f"{np.median(allrun):.4f} ns over {len(allrun)} modules")
    else:
        pr("  sigma_t_modules.tsv not yet written")
    (out_dir / "validation.txt").write_text("\n".join(L) + "\n")
    print("\n".join(L))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sp = ap.add_subparsers(dest="cmd", required=True)
    d = sp.add_parser("dump")
    d.add_argument("files", nargs="*")
    d.add_argument("--out", default=os.path.join(HERE, OUT_DIR, "dump"))
    d.add_argument("--max-events", type=int, default=0)
    d.add_argument("--parallel", type=int, default=4)
    d.add_argument("--templates-root", default=os.path.join(HERE, TEMPLATES_BEST50))
    d.set_defaults(fn=cmd_dump)
    z = sp.add_parser("analyze")
    z.add_argument("--dump-dir", default=os.path.join(HERE, OUT_DIR, "dump"))
    z.add_argument("--rf-dir", default=os.path.join(HERE, RF_DIR))
    z.add_argument("--out-dir", default=os.path.join(HERE, OUT_DIR))
    z.add_argument("--jobs", type=int, default=4)
    z.set_defaults(fn=cmd_analyze)
    v = sp.add_parser("validate")
    v.add_argument("--out-dir", default=os.path.join(HERE, OUT_DIR))
    v.add_argument("--ref-dump-dir", default=os.path.join(HERE, REF_DUMP_DIR))
    v.add_argument("--ref-rf-calib", default=os.path.join(HERE, REF_RF_CALIB))
    v.set_defaults(fn=cmd_validate)
    a = ap.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()
