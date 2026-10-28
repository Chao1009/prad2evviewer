#!/usr/bin/env python3
"""plot_module_timing.py - per-module timing histograms from module_timing_dump.py npz files.

For every PbWO4 module, histograms of  dt = t_CFD(module) - t_ref  in four panels
(fixed order 0.73 GeV H2, 2.24 GeV H2, 3.49 GeV H2, X17 Ta), each with an iterative
+-2 sigma Gaussian core fit and mean / sigma / N annotated, plus summary pages.
Two reference choices, one PDF each:
  --ref veto  ->  <out>/module_timing.pdf          (t_ref = in-time veto pulse time)
  --ref evt   ->  <out>/module_timing_evtref.pdf   (t_ref = leave-one-out event median)
With --time-calib the file names get a _tcal suffix.

DEFINITIONS
  Event selection   trig & --trigger-mask (default 0xf00; physics trigger bits), and not an LMS flasher
                    (bit 24) or alpha pulser (bit 25) event (trig & 0x3000000 == 0).
  Module time       t_m = m_t (50 % leading edge, ns) of the module's dominant mf pulse
                    (dump definition); only pulses with m_E >= --min-e if --min-e > 0.
                    --time-calib subtracts the per-module offset of
                    database/hycal_time_offsets/time_calib.json ("resolved": true entries,
                    matched by name or module_id; others 0) from EVERY PbWO4 time, module
                    and reference pulses alike.
  --ref evt         t_ref = leave-one-out median of t over the OTHER PbWO4 dominant pulses
                    of the same event with m_E > --ref-min-e, needing >= --ref-min-n of
                    them (the module itself is excluded exactly, vectorised).
  --ref veto        per event the V1-V4 pulses with v_int > --veto-min-int and
                    |(v_t - t_evt) - D0| < --veto-gate; t_evt = full-event median of the
                    eligible PbWO4 times (m_E > --ref-min-e, >= --ref-min-n pulses); D0 =
                    coincidence-peak position of v_t - t_evt per run type, found from the
                    data (first pass: histogram, peak, Gaussian+const core fit).  Several
                    in-time pulses: largest integral.  t_ref = that v_t.
  --veto-gate none  --ref veto only: writes module_timing_nogate.pdf (+ TSV) INSTEAD of the gated
                    outputs (module_timing.pdf etc. are not rewritten): dt for EVERY
                    (module dominant pulse) x (veto pulse, integral > --veto-min-int) pair of an
                    event, aligned counters, no gate, no largest-integral choice; per panel the
                    in-time peak, accidental level from sidebands 20-60 ns off the peak, S/B in
                    +-10 ns; same pages as module_timing.pdf.  The 10 ns gate is still used for
                    the run selection (D0) and the gated histograms that fix the page list.
  --veto-align      The counters are not time-aligned in the database.  counter (default):
                    before the pick, v_t of Vk -> v_t - c_k, c_k from the unaligned dt medians
                    of each (module, counter) within +-5 ns of the module peak, fitted as
                    a_module + b_counter (c_k = -b_k, entry-weighted mean of c_k = 0), so module
                    time offsets cancel.  centroid: c_k = centroid of v_t - t_evt of Vk minus D0
                    (biased: t_evt carries the raw offsets of the region each counter covers).
                    type: no shift (old behaviour).  Non-default modes get a _<mode>align
                    suffix on the --ref veto outputs (e.g. _typealign).  The pre-alignment
                    output of the previous code is kept as module_timing_valign_type.pdf.
  No-veto run types --ref veto, run types without a usable veto run (X17 Ta): dt = t_m - t_seed,
                    t_seed = time of the highest-m_E PbWO4 pulse of the event other than the
                    module itself (m_E > --ref-min-e); for the event's seed module that is the
                    second-highest pulse (t_seed - t_2nd).  Pulse-level seed (no clusters in
                    the dump).  Same fit, +-60 ns histogram.
  Gaussian core fit histogram argmax + half-maximum width as start, then least-squares
                    Gaussian in [mu-2s, mu+2s] refitted until mu and s change < 1 % of s
                    or 5 refits (0.2 ns bins).
  Module page       N > --min-n in at least one run type (N = entries of the histogram,
                    |dt| range +-60 ns evt, 140 ns window around -D0 veto); panels with
                    N <= --min-n are drawn but not fitted.

Outputs (--out-dir, default output/module_timing_20261025): module_timing{,_evtref}[_tcal].pdf,
  module_timing_{veto,evtref}[_tcal].tsv, module_timing_{veto,evtref}[_tcal]_summary.txt,
  veto_summary[_tcal].tsv (veto availability, every run; drawn on the first pages of the
  veto PDF).

Usage:
  python3 plot_module_timing.py --ref veto [--dump-dir DIR] [--out-dir DIR] [--min-n 10000]
      [--trigger-mask 0xf00] [--min-e 0] [--ref-min-e 10] [--ref-min-n 3]
      [--veto-min-int 2000] [--veto-gate 10] [--time-calib] [--jobs 12]
"""

from __future__ import annotations

import argparse
import contextlib
import glob
import hashlib
import json
import multiprocessing as mp
import os
import re
import sys
import time
import warnings
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

# (tag, label, runs): fixed panel order
RUN_TYPES = [
    ("0p73", "0.73 GeV H2", (24590, 24603)),
    ("2p24", "2.24 GeV H2", (25228, 25257, 25260, 25299, 25320)),
    ("3p49", "3.49 GeV H2", (24386, 24503, 24512)),
    ("x17", "X17 Ta", (25800, 26138, 26230)),
]
TAGS = [t for t, _, _ in RUN_TYPES]
LABEL = {t: l for t, l, _ in RUN_TYPES}
TYPE_OF_RUN = {r: t for t, _, rs in RUN_TYPES for r in rs}
RUN_NAME = r"(\d{6})_(\d{5})\.npz$"

W_BIN = 0.2                  # ns, module histograms
EVT_LO, EVT_NB = -60.0, 600  # dt = t_m - t_loo, [-60, 60]
VETO_HALF = 70.0   # dt = t_m - t_veto histogram: +-70 ns around -D0 (union over the runs of a type)
D_LO, D_W, D_NB = -250.0, 0.5, 1000   # v_t - t_evt histogram
A_LO, A_W, A_NB = 0.0, 1.0, 400       # absolute v_t histogram
DEFAULT_D0 = -35.0           # only used to place the veto histogram window when D0 is unknown
NMAX_JOBS = 12
TIME_CALIB = HERE.parent.parent / "database" / "hycal_time_offsets" / "time_calib.json"
INT_LEVELS = (500.0, 2000.0)
CNT_WIN = 30.0               # ns, per-counter centroid fit window around D0
NG_LO, NG_W, NG_NB = -450.0, 0.5, 1800     # --veto-gate none: all-pairs dt histogram
NG_PK, NG_SB, NG_SIG = 15.0, (20.0, 60.0), 10.0   # peak search +-, sideband |dt - peak| range, signal +- (ns)
CNT_W, CNT_NB, CNT_HALF = 0.1, 300, 15.0   # per (module, counter) dt histograms, -D0 +- 15 ns
CNT_MIN_MOD, CNT_MIN_MC, CNT_CORE = 200, 50, 5.0   # module / module-counter entries; median window +- ns
TRIG_EXCL = 0x3000000  # LMS flasher (bit 24) / alpha pulser (bit 25) events are never selected


# ---- event references (vectorised) ------------------------------------------

def _med_sorted(s, start, m, k):
    """Median of the m values s[start + j + (j >= k)], j = 0..m-1 (value at rank k
    skipped; k >= m -> nothing skipped).  s is sorted within each event block."""
    j1 = (m - 1) // 2
    j2 = m // 2
    i1 = start + j1 + (j1 >= k)
    i2 = start + j2 + (j2 >= k)
    return 0.5 * (s[i1] + s[i2])


def event_refs(m_ev, t, E, sel, n_ev, ref_min_e, ref_min_n):
    """t_evt[n_ev]: full-event median of eligible PbWO4 times (nan if fewer than
    ref_min_n eligible).  t_loo[n_rows]: leave-one-out median for every row of a
    selected event (eligible rows exclude themselves; others use the full set; nan if
    fewer than ref_min_n remain)."""
    ref_min_n = max(int(ref_min_n), 1)
    m_ev = m_ev.astype(np.int64)
    rowsel = sel[m_ev]
    elig = rowsel & (E > ref_min_e)
    ie = np.flatnonzero(elig)
    ev = m_ev[ie]
    tt = t[ie]
    order = np.lexsort((tt, ev))
    s = tt[order]
    es = ev[order]
    cnt = np.bincount(es, minlength=n_ev)
    start = np.cumsum(cnt) - cnt
    rank_sorted = np.arange(len(s)) - start[es]
    rank = np.empty(len(ie), np.int64)
    rank[order] = rank_sorted
    t_evt = np.full(n_ev, np.nan)
    ok = np.flatnonzero(cnt >= ref_min_n)
    if ok.size:
        t_evt[ok] = _med_sorted(s, start[ok], cnt[ok], cnt[ok])
    k_row = np.full(len(m_ev), np.iinfo(np.int64).max // 4, np.int64)
    k_row[ie] = rank
    m_row = cnt[m_ev] - elig
    valid = rowsel & (m_row >= ref_min_n)
    t_loo = np.full(len(m_ev), np.nan)
    iv = np.flatnonzero(valid)
    if iv.size:
        t_loo[iv] = _med_sorted(s, start[m_ev[iv]], m_row[iv], k_row[iv])
    return t_evt, t_loo


def seed_refs(m_ev, t, E, sel, ref_min_e):
    """Per row of a selected event: (t_ref, is_seed).  t_ref = time of the highest-E eligible
    pulse (E > ref_min_e) of the event other than the row itself (the seed; for the seed row the
    second-highest), nan if none.  Ties: first row in dump order."""
    m_ev = m_ev.astype(np.int64)
    ie = np.flatnonzero(sel[m_ev] & (E > ref_min_e))
    t_ref = np.full(len(m_ev), np.nan)
    is_seed = np.zeros(len(m_ev), bool)
    if not ie.size:
        return t_ref, is_seed
    o = ie[np.lexsort((ie, -E[ie], m_ev[ie]))]
    ev = m_ev[o]
    first = np.r_[True, ev[1:] != ev[:-1]]
    i1 = o[first]                                   # seed row per event
    sec = np.r_[False, first[:-1]] & ~first         # second row of the same event
    i2 = o[sec]
    n_ev = int(m_ev.max()) + 1
    s1 = np.full(n_ev, -1, np.int64)
    s2 = np.full(n_ev, -1, np.int64)
    s1[m_ev[i1]] = i1
    s2[m_ev[i2]] = i2
    rows = np.flatnonzero(sel[m_ev])
    ref = s1[m_ev[rows]]
    me = ref == rows
    ref[me] = s2[m_ev[rows[me]]]
    ok = ref >= 0
    t_ref[rows[ok]] = t[ref[ok]]
    is_seed[rows[me]] = True
    return t_ref, is_seed


def pick_veto_idx(v_ev, v_t, v_int, sel, t_evt, d0, gate, min_int, n_ev):
    """Per event the index of the in-time veto pulse with the largest integral (-1 if none),
    and the in-time pulse mask over all pulses."""
    ok = sel[v_ev] & (v_int > min_int) & np.isfinite(t_evt[v_ev])
    it = np.zeros(len(v_ev), bool)
    idx = np.full(n_ev, -1, np.int64)
    if not np.isfinite(d0) or not ok.any():
        return idx, it
    it = ok & (np.abs(v_t - t_evt[v_ev] - d0) < gate)
    ii = np.flatnonzero(it)
    if ii.size:
        o = np.lexsort((v_int[ii], v_ev[ii]))
        ii = ii[o]
        e = v_ev[ii]
        last = np.r_[e[1:] != e[:-1], True]
        idx[e[last]] = ii[last]
    return idx, it


def pick_veto(v_ev, v_ch, v_t, v_int, sel, t_evt, d0, gate, min_int, n_ev):
    """Per event the in-time veto pulse with the largest integral: returns
    (t_veto[n_ev] nan if none, in-time pulse mask over all pulses)."""
    idx, it = pick_veto_idx(v_ev, v_t, v_int, sel, t_evt, d0, gate, min_int, n_ev)
    t_veto = np.full(n_ev, np.nan)
    k = idx >= 0
    t_veto[k] = v_t[idx[k]]
    return t_veto, it


def pass_cnt(job):
    """Per (PbWO4 module, chosen veto counter) histogram of the UNALIGNED dt = t_m - t_veto,
    [lo, lo + CNT_NB * CNT_W), for the module-referenced counter alignment."""
    path, cfg, d0, lo = job
    d, n_ev, sel, m_ev, m_idx, t, E, t_evt, t_loo, v = _load(path, cfg)
    npw, pmap = cfg["npw"], cfg["pmap"]
    out = dict(run=int(d["run"]), h=np.zeros((npw, 4, CNT_NB), np.int32))
    if not np.isfinite(d0):
        return out
    idx = pick_veto_idx(v["ev"], v["t"], v["i"], sel, t_evt, d0, cfg["gate"], cfg["veto_min_int"], n_ev)[0]
    r = sel[m_ev] & (pmap[m_idx] >= 0) & (idx[m_ev] >= 0)
    if cfg.get("min_e", 0) > 0:
        r &= E >= cfg["min_e"]
    k = idx[m_ev[r]]
    dt = t[r] - v["t"][k]
    b = np.floor((dt - lo) / CNT_W).astype(np.int64)
    ok = (b >= 0) & (b < CNT_NB)
    key = (pmap[m_idx[r]][ok] * 4 + v["ch"][k][ok] - 1) * CNT_NB + b[ok]
    out["h"] = np.bincount(key, minlength=npw * 4 * CNT_NB).reshape(npw, 4, CNT_NB).astype(np.int32)
    return out


# ---- workers ----------------------------------------------------------------

def _load(path, cfg):
    d = np.load(path)
    names = d["module_names"]
    h = hashlib.md5("\n".join(map(str, names)).encode()).hexdigest()
    if h != cfg["names_hash"]:
        raise RuntimeError(f"{path}: module_names differ from the geometry")
    n_ev = len(d["trig"])
    trig = d["trig"]
    sel = ((trig & np.uint32(cfg["mask"])) != 0) & ((trig & np.uint32(TRIG_EXCL)) == 0)
    m_idx = d["m_idx"].astype(np.int64)
    t = d["m_t"].astype(np.float64)
    if cfg["off"] is not None:
        t = t - cfg["off"][m_idx]
    E = d["m_E"].astype(np.float64)
    m_ev = d["m_ev"].astype(np.int64)
    t_evt, t_loo = event_refs(m_ev, t, E, sel, n_ev, cfg["ref_min_e"], cfg["ref_min_n"])
    v = dict(ev=d["v_ev"].astype(np.int64), ch=d["v_ch"].astype(np.int64),
             t=d["v_t"].astype(np.float64), i=d["v_int"].astype(np.float64))
    return d, n_ev, sel, m_ev, m_idx, t, E, t_evt, t_loo, v


def _nuniq(a):
    return int(np.unique(a).size)


def pass1(job):
    """Veto availability counts (Delta0-independent) and v_t - t_evt histograms."""
    path, cfg = job
    d, n_ev, sel, m_ev, m_idx, t, E, t_evt, t_loo, v = _load(path, cfg)
    vs = sel[v["ev"]]
    ve, vc, vt, vi = v["ev"][vs], v["ch"][vs], v["t"][vs], v["i"][vs]
    out = dict(run=int(d["run"]), n_read=int(d["n_read"]), n_dump=n_ev, n_sel=int(sel.sum()),
               n_tevt=int(np.isfinite(t_evt[sel]).sum()), n_any=_nuniq(ve),
               n_500=_nuniq(ve[vi > INT_LEVELS[0]]), n_2000=_nuniq(ve[vi > INT_LEVELS[1]]),
               ch_any=[_nuniq(ve[vc == c]) for c in range(1, 5)],
               ch_2000=[_nuniq(ve[(vc == c) & (vi > INT_LEVELS[1])]) for c in range(1, 5)])
    g = vi > cfg["veto_min_int"]
    ok = g & np.isfinite(t_evt[ve])
    dd = vt[ok] - t_evt[ve[ok]]
    out["h_d"] = np.histogram(dd, bins=D_NB, range=(D_LO, D_LO + D_NB * D_W))[0].astype(np.int64)
    vco = vc[ok]
    out["h_dc"] = np.stack([np.histogram(dd[vco == c], bins=D_NB, range=(D_LO, D_LO + D_NB * D_W))[0]
                            for c in range(1, 5)]).astype(np.int64)
    out["h_abs"] = np.histogram(vt[g], bins=A_NB, range=(A_LO, A_LO + A_NB * A_W))[0].astype(np.int64)
    return out


def _hist2d(pidx, dt, lo, w, nb, npw):
    b = np.floor((dt - lo) / w).astype(np.int64)
    ok = (b >= 0) & (b < nb)
    h = np.bincount(pidx[ok] * nb + b[ok], minlength=npw * nb).reshape(npw, nb)
    n_all = np.bincount(pidx, minlength=npw)
    return h.astype(np.int32), n_all.astype(np.int64)


def pass2(job):
    """Per-module histograms for both references plus in-time veto statistics.  off[ch]:
    per-counter veto time shift (index 1..4); seed_ref: the 'veto' histogram holds t_m - t_seed."""
    path, cfg, d0_stat, d0_hist, vlo, vnb, off, seed_ref = job
    d, n_ev, sel, m_ev, m_idx, t, E, t_evt, t_loo, v = _load(path, cfg)
    npw = cfg["npw"]
    pmap = cfg["pmap"]
    rowok = sel[m_ev] & (pmap[m_idx] >= 0)
    if cfg["min_e"] > 0:
        rowok &= E >= cfg["min_e"]
    vt = v["t"] - off[v["ch"]]
    t_stat, it = pick_veto(v["ev"], v["ch"], vt, v["i"], sel, t_evt, d0_stat,
                           cfg["gate"], cfg["veto_min_int"], n_ev)
    if d0_hist == d0_stat or (np.isnan(d0_hist) and np.isnan(d0_stat)):
        t_veto = t_stat
    else:
        t_veto = pick_veto(v["ev"], v["ch"], vt, v["i"], sel, t_evt, d0_hist,
                           cfg["gate"], cfg["veto_min_int"], n_ev)[0]
    out = dict(run=int(d["run"]))
    r = rowok & np.isfinite(t_loo)
    out["h_evt"], out["n_evt"] = _hist2d(pmap[m_idx[r]], t[r] - t_loo[r], EVT_LO, W_BIN, EVT_NB, npw)
    if seed_ref:
        t_s, is_seed = seed_refs(m_ev, t, E, sel, cfg["ref_min_e"])
        r = rowok & np.isfinite(t_s)
        out["h_veto"], out["n_veto"] = _hist2d(pmap[m_idx[r]], t[r] - t_s[r], vlo, W_BIN, vnb, npw)
        out["n_seed"] = np.bincount(pmap[m_idx[r & is_seed]], minlength=npw).astype(np.int64)
    else:
        r = rowok & np.isfinite(t_veto[m_ev])
        out["h_veto"], out["n_veto"] = _hist2d(pmap[m_idx[r]], t[r] - t_veto[m_ev[r]], vlo, W_BIN, vnb, npw)
    if cfg.get("nogate") and not seed_ref and np.isfinite(d0_hist):
        out["h_ng"] = _all_pairs(m_ev, m_idx, t, rowok, v, vt, sel, cfg, n_ev)
    out["n_intime"] = int(np.isfinite(t_stat).sum())
    out["ch_intime"] = [_nuniq(v["ev"][it & (v["ch"] == c)]) for c in range(1, 5)]
    return out


def _all_pairs(m_ev, m_idx, t, rowok, v, vt, sel, cfg, n_ev):
    """--veto-gate none: per module, histogram of t_m - v_t' over EVERY (module row, veto pulse with
    integral > veto_min_int) pair of the same selected event (rows are sorted by event)."""
    ri = np.flatnonzero(rowok)
    cnt = np.bincount(m_ev[ri], minlength=n_ev)
    start = np.cumsum(cnt) - cnt
    vi = np.flatnonzero(sel[v["ev"]] & (v["i"] > cfg["veto_min_int"]))
    k = cnt[v["ev"][vi]]
    rv = np.repeat(vi, k)
    within = np.arange(k.sum()) - np.repeat(np.cumsum(k) - k, k)
    row = ri[start[v["ev"][rv]] + within]
    return _hist2d(cfg["pmap"][m_idx[row]], t[row] - vt[rv], NG_LO, NG_W, NG_NB, cfg["npw"])[0]


def _run_pool(func, jobs, nproc, label):
    t0 = time.time()
    res = []
    if nproc <= 1:
        for k, j in enumerate(jobs, 1):
            res.append(func(j))
    else:
        ctx = mp.get_context("spawn")
        with ctx.Pool(nproc) as pool:
            for k, r in enumerate(pool.imap_unordered(func, jobs), 1):
                res.append(r)
                if k % 20 == 0 or k == len(jobs):
                    print(f"  [{label}] {k}/{len(jobs)} files, {time.time() - t0:.0f} s", flush=True)
    return res


# ---- fitting ----------------------------------------------------------------

def fit_core(x, c, bg=False, smooth_ns=1.0, max_it=5, tol=0.01):
    """Iterative +-2 sigma Gaussian core fit of a histogram (centres x, counts c).
    Start: smoothed argmax and half-maximum width; refit in [mu-2s, mu+2s] until mu and
    s change by < tol*s or max_it refits.  bg=True adds a constant (for a peak on a flat
    background).  Returns dict(ok, mu, sigma, mu_err, sigma_err, lo, hi, nfit, nit, bg)."""
    from scipy.optimize import curve_fit
    res = dict(ok=False, mu=np.nan, sigma=np.nan, mu_err=np.nan, sigma_err=np.nan,
               lo=np.nan, hi=np.nan, nfit=0, nit=0, amp=np.nan, bg=0.0, nbin=0)
    c = np.asarray(c, float)
    w = x[1] - x[0]
    if c.sum() < 20:
        return res
    k = max(3, int(round(smooth_ns / w)) | 1)
    sm = np.convolve(c, np.ones(k) / k, mode="same")
    i0 = int(np.argmax(sm))
    base = float(np.median(c)) if bg else 0.0
    half = base + 0.5 * (sm[i0] - base)
    lo_i = hi_i = i0
    while lo_i > 0 and sm[lo_i - 1] > half:
        lo_i -= 1
    while hi_i < len(sm) - 1 and sm[hi_i + 1] > half:
        hi_i += 1
    mu = float(x[i0])
    sigma = max((x[hi_i] - x[lo_i] + w) / 2.355, 1.5 * w)
    sigma0 = sigma

    def gauss(xx, a, m, s):
        return a * np.exp(-0.5 * ((xx - m) / s) ** 2)

    def gauss_bg(xx, a, m, s, b):
        return gauss(xx, a, m, s) + b

    popt = None
    for it in range(1, max_it + 1):
        win = (x >= mu - 2 * sigma) & (x <= mu + 2 * sigma)
        if win.sum() < 4 or c[win].sum() < 10:
            return res
        cmax = float(c[win].max())
        p0 = [max(cmax - base, 1e-3), mu, sigma] + ([base] if bg else [])
        lower = [0.0, mu - 2 * sigma, 0.3 * w] + ([0.0] if bg else [])
        upper = [max(2 * cmax, 1.0), mu + 2 * sigma, 4 * sigma0] + ([max(cmax, 1.0)] if bg else [])
        p0 = list(np.clip(p0, np.array(lower) + 1e-9, np.array(upper) - 1e-9))
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                popt, pcov = curve_fit(gauss_bg if bg else gauss, x[win], c[win], p0=p0,
                                       sigma=np.sqrt(np.maximum(c[win], 1.0)), maxfev=4000,
                                       bounds=(lower, upper))
        except Exception:
            return res
        mu_n, s_n = float(popt[1]), abs(float(popt[2]))
        if not np.isfinite(s_n) or s_n < 0.3 * w or s_n > 0.5 * (x[-1] - x[0]):
            return res
        conv = abs(mu_n - mu) < tol * s_n and abs(s_n - sigma) < tol * s_n
        lo_w, hi_w = mu - 2 * sigma, mu + 2 * sigma
        mu, sigma = mu_n, s_n
        err = np.sqrt(np.diag(pcov))
        res.update(ok=True, mu=mu, sigma=sigma, mu_err=float(err[1]), sigma_err=float(err[2]),
                   lo=lo_w, hi=hi_w, nfit=int(c[win].sum()), nbin=int(win.sum()), nit=it, amp=float(popt[0]),
                   bg=float(popt[3]) if bg else 0.0)
        if conv:
            break
    return res


def peak_significant(fit, thr=8.0, smin=0.4, smax=3.0):
    """Gaussian+const fit: excess counts in the fit window over the fitted flat level, in
    units of sqrt(expected) > thr, and width within [smin, smax] ns (the coincidence peak is
    1-2 ns wide; narrow spikes and broad bumps are not it).  Tells a coincidence peak from
    a noise maximum."""
    if fit is None or not fit["ok"] or not smin <= fit["sigma"] <= smax:
        return False
    exp = fit["bg"] * fit["nbin"]
    return (fit["nfit"] - exp) / np.sqrt(exp + 1.0) > thr


def centres(lo, w, nb):
    return lo + (np.arange(nb) + 0.5) * w


# ---- veto coincidence peak and counter alignment -----------------------------

def veto_d0(h_by_run, runs_of, gate):
    """D0 per run type from the per-run v_t - t_evt histograms (all counters).  Returns
    d0[tag], d0fit[tag], d0run[run] (fit or None), d0_stat[run] (gate centre for the in-time
    statistics) and d0_hist[run] (same, nan when the run is not used for the veto reference)."""
    xd = centres(D_LO, D_W, D_NB)
    # coincidence peak: searched in |v_t - t_evt| < 100 ns (excludes the early-window artefact
    # at v_t ~ 2-7 ns); a run has a peak only if it is significant (excess > 8 sqrt(flat level), 0.4 < sigma < 3 ns)
    sm = (xd >= -100.0) & (xd <= 100.0)
    d0_fit = lambda h: fit_core(xd[sm], h[sm], bg=True, smooth_ns=3.0) if h[sm].sum() >= 50 else None
    d0run, has_peak = {}, {}
    for tag in TAGS:
        for run in runs_of[tag]:
            fr = d0_fit(h_by_run[run])
            has_peak[run] = peak_significant(fr)
            d0run[run] = fr if has_peak[run] else None
            print(f"    run {run}: peak D0 = {f(fr['mu'], 2) if fr else 'nan'} ns, sigma = "
                  f"{f(fr['sigma'], 2) if fr else 'nan'} ns, significant: {has_peak[run]}", flush=True)
    d0, d0fit = {}, {}
    for tag in TAGS:
        hd = sum((h_by_run[r] for r in runs_of[tag] if has_peak[r]), np.zeros(D_NB, np.int64))
        fit = d0_fit(hd)
        if not peak_significant(fit):
            fit = dict(ok=False, mu=np.nan, sigma=np.nan)
        d0fit[tag] = fit
        d0[tag] = fit["mu"] if fit["ok"] else np.nan
        if runs_of[tag]:
            print(f"  D0[{LABEL[tag]}] = {f(d0[tag], 2)} ns  sigma = {f(fit['sigma'], 2)} ns "
                  f"(runs with a peak: {[r for r in runs_of[tag] if has_peak[r]]})", flush=True)
    # d0_hist = nan (run not used for the veto-referenced histograms) when the run has no
    # significant peak or its peak is not within the gate of the type's D0
    d0_stat, d0_hist = {}, {}
    for tag in TAGS:
        for run in runs_of[tag]:
            d0_stat[run] = d0[tag]
            ok = has_peak[run] and np.isfinite(d0[tag]) and abs(d0run[run]["mu"] - d0[tag]) < gate
            d0_hist[run] = d0[tag] if ok else np.nan
            if tag != "x17" and not ok:
                print(f"    run {run}: excluded from the veto-referenced histograms (no coincidence peak "
                      f"at the gate centre)", flush=True)
    return d0, d0fit, d0run, d0_stat, d0_hist


def centroid_offsets(hc_by_run, runs_of, d0, d0_hist):
    """Per run type and counter: centroid c_k of the counter's own v_t - t_evt coincidence peak
    (usable runs pooled, +-CNT_WIN around D0; Gaussian+const core fit) and c_k - D0.  Biased as a
    counter offset: without time calibration t_evt carries the raw offsets of the HyCal region
    each counter covers.  Returns (off{tag: array(5)}, info{tag: {Vk: dict}})."""
    xd = centres(D_LO, D_W, D_NB)
    off, info = {}, {}
    for tag in TAGS:
        o = np.zeros(5)
        info[tag] = {}
        use = [r for r in runs_of[tag] if np.isfinite(d0_hist[r])]
        if use:
            w = np.abs(xd - d0[tag]) < CNT_WIN
            for c in range(1, 5):
                h = sum(hc_by_run[r][c - 1] for r in use)
                fc = fit_core(xd[w], h[w], bg=True, smooth_ns=3.0) if h[w].sum() >= 50 else None
                good = peak_significant(fc)
                if good:
                    o[c] = fc["mu"] - d0[tag]
                info[tag][f"V{c}"] = dict(centroid=round(float(fc["mu"]), 4) if good else None,
                                          sigma=round(float(fc["sigma"]), 4) if good else None,
                                          centroid_minus_D0=round(float(o[c]), 4))
        off[tag] = o
    return off, info


def _hist_median(x, h):
    c = np.cumsum(h)
    j = int(np.searchsorted(c, 0.5 * c[-1]))
    lo = c[j - 1] if j else 0
    w = x[1] - x[0]
    return x[j] - w / 2 + w * (0.5 * c[-1] - lo) / max(h[j], 1)


def module_offsets(H, lo):
    """Module-referenced counter offsets from H[module, counter, bin] (unaligned dt): per module
    the dt median of each counter within +-CNT_CORE ns of the module's peak, then the additive
    model median[m, k] = a_m + b_k (weights = entries, alternating least squares).  The counter
    shift is c_k = -b_k with sum_k W_k c_k = 0 (W_k = entries of counter k), so the coincidence
    peak stays at D0 on average.  Returns (off array(5), info)."""
    x = lo + (np.arange(CNT_NB) + 0.5) * CNT_W
    k5 = int(round(0.5 / CNT_W)) | 1
    P, Wt = [], []
    for m in range(H.shape[0]):
        tot = H[m].sum(axis=0)
        if tot.sum() < CNT_MIN_MOD:
            continue
        pk = x[int(np.argmax(np.convolve(tot, np.ones(k5) / k5, mode="same")))]
        win = np.abs(x - pk) < CNT_CORE
        p, w = np.full(4, np.nan), np.zeros(4)
        for c in range(4):
            hc = H[m, c][win]
            if hc.sum() >= CNT_MIN_MC:
                p[c], w[c] = _hist_median(x[win], hc), hc.sum()
        if (w > 0).sum() >= 2:
            P.append(np.nan_to_num(p))
            Wt.append(w)
    off = np.zeros(5)
    info = dict(n_modules=len(P))
    if len(P) < 3:
        return off, info
    P, Wt = np.array(P), np.array(Wt)
    b = np.zeros(4)
    for _ in range(200):
        a = (Wt * (P - b)).sum(axis=1) / Wt.sum(axis=1)
        cw = Wt.sum(axis=0)
        b = np.where(cw > 0, (Wt * (P - a[:, None])).sum(axis=0) / np.maximum(cw, 1), 0.0)
    cw = Wt.sum(axis=0)
    b -= (cw * b).sum() / cw.sum()
    res = P - a[:, None] - b
    rms = np.sqrt((Wt * res ** 2).sum() / Wt.sum())
    off[1:] = -b
    for c in range(4):
        info[f"V{c + 1}"] = dict(shift=round(float(-b[c]), 4), n_modules=int((Wt[:, c] > 0).sum()),
                                 entries=int(cw[c]))
    info["residual_rms_ns"] = round(float(rms), 4)
    return off, info


def measure_counter_offsets(files, cfg, runs_of, type_of_run, d0, d0_hist, p1_hc, align, nproc):
    """Veto counter shifts per run type for --veto-align: 'counter' (module-referenced, default),
    'centroid' (v_t - t_evt centroids) or 'type' (none).  Both methods are always measured and
    returned in info."""
    coff, cinfo = centroid_offsets(p1_hc, runs_of, d0, d0_hist)
    lo = {t: -d0[t] - CNT_HALF for t in TAGS}
    run_of = lambda fp: int(re.search(RUN_NAME, fp).group(1))
    jobs = [(fp, cfg, d0_hist[run_of(fp)], lo[type_of_run[run_of(fp)]]) for fp in files
            if np.isfinite(d0_hist[run_of(fp)])]
    H = {}
    for r in _run_pool(pass_cnt, jobs, nproc, "counter"):
        tg = type_of_run[r["run"]]
        H[tg] = H.get(tg, 0) + r["h"]
    off, info = {}, {}
    for t in TAGS:
        moff, minfo = module_offsets(H[t], lo[t]) if t in H else (np.zeros(5), {})
        off[t] = moff if align == "counter" else coff[t] if align == "centroid" else np.zeros(5)
        info[t] = dict(D0=round(float(d0[t]), 4) if np.isfinite(d0[t]) else None,
                       runs_used=[r for r in runs_of[t] if np.isfinite(d0_hist[r])],
                       applied=[round(float(v), 4) for v in off[t][1:]],
                       module_referenced=minfo, vt_minus_tevt_centroid=cinfo[t])
    return off, info


# ---- geometry / time calibration -------------------------------------------

def load_geometry(run):
    """Return dict of full-length (p.hycal index order) arrays: name, id, x, y, sx, sy, pwo4."""
    import _common as C
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
    n = p.hycal.module_count()
    ms = [p.hycal.module(i) for i in range(n)]
    return dict(name=np.array([m.name for m in ms]), id=np.array([m.id for m in ms]),
                x=np.array([m.x for m in ms]), y=np.array([m.y for m in ms]),
                sx=np.array([m.size_x for m in ms]), sy=np.array([m.size_y for m in ms]),
                pwo4=np.array([m.is_pwo4() for m in ms]))


def load_offsets(path, geo):
    j = json.load(open(path))
    by_name, by_id = {}, {}
    for m in j["modules"]:
        if m.get("resolved"):
            by_name[m["name"]] = float(m["offset_ns"])
            by_id[int(m["module_id"])] = float(m["offset_ns"])
    off = np.zeros(len(geo["name"]))
    nres = 0
    missing = []
    for i, (nm, mid) in enumerate(zip(geo["name"], geo["id"])):
        o = by_name.get(str(nm), by_id.get(int(mid)))
        if o is None:
            if geo["pwo4"][i]:
                missing.append(str(nm))
        else:
            off[i] = o
            nres += geo["pwo4"][i]
    return off, nres, missing


# ---- summary helpers --------------------------------------------------------

def pct(a, p):
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    return float(np.percentile(a, p)) if a.size else float("nan")


def f(x, nd=3):
    return "nan" if x is None or not np.isfinite(x) else f"{x:.{nd}f}"


def outliers(names, mean, sigma, N):
    ok = np.isfinite(mean) & np.isfinite(sigma)
    if ok.sum() < 5:
        return []
    mm = np.median(mean[ok])
    ms = np.median(sigma[ok])
    mad = 1.4826 * np.median(np.abs(mean[ok] - mm))
    rows = []
    for i in np.flatnonzero(ok):
        why = []
        if mad > 0 and abs(mean[i] - mm) > 5 * mad:
            why.append(f"|mean-med|={abs(mean[i] - mm):.2f}>5*{mad:.2f}")
        if sigma[i] > 2 * ms:
            why.append(f"sigma>{2 * ms:.2f}")
        if why:
            rows.append((names[i], mean[i], sigma[i], N[i], "; ".join(why)))
    return rows


# ---- plotting ---------------------------------------------------------------

def _mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def text_page(pdf, lines, title=None, size=9):
    plt = _mpl()
    fig = plt.figure(figsize=(11, 8.5))
    if title:
        fig.text(0.04, 0.95, title, fontsize=15, weight="bold", va="top")
    fig.text(0.04, 0.90, "\n".join(lines), fontsize=size, va="top", family="monospace",
             linespacing=1.35)
    pdf.savefig(fig)
    plt.close(fig)


def draw_panel(ax, tag, h, res, n_all, lo, xlabel, min_n, empty_msg, note=None):
    x = centres(lo, W_BIN, len(h))
    N = int(h.sum())
    ax.set_title(LABEL[tag], fontsize=10)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(f"entries / {W_BIN:g} ns", fontsize=8)
    ax.tick_params(labelsize=8)
    if N == 0:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.text(0.5, 0.5, empty_msg, transform=ax.transAxes, ha="center", va="center",
                fontsize=9, color="0.35")
        return
    if res is not None and res["ok"]:
        half = max(8 * res["sigma"], 10.0)
        x0, x1 = res["mu"] - half, res["mu"] + half
    else:
        k = int(round(1.0 / W_BIN)) | 1
        sm = np.convolve(h, np.ones(k) / k, mode="same")
        pk = x[int(np.argmax(sm))]
        x0, x1 = pk - 15.0, pk + 15.0
    sel = (x >= x0) & (x <= x1)
    edges = np.r_[x[sel] - W_BIN / 2, x[sel][-1] + W_BIN / 2]
    ax.stairs(h[sel], edges, color="C0", lw=0.9)
    txt = []
    if res is not None and res["ok"]:
        xs = np.linspace(x0, x1, 500)
        g = lambda xx: res["amp"] * np.exp(-0.5 * ((xx - res["mu"]) / res["sigma"]) ** 2)
        ax.plot(xs, g(xs), color="crimson", lw=0.8, ls="--", alpha=0.6)
        xw = np.linspace(res["lo"], res["hi"], 200)
        ax.plot(xw, g(xw), color="crimson", lw=1.8)
        txt = [f"mean = {res['mu']:.3f} ns", f"sigma = {res['sigma']:.3f} ns", f"N = {N}",
               f"fit [{res['lo']:.2f}, {res['hi']:.2f}] ns"]
    elif N > min_n:
        txt = [f"N = {N}", "fit failed"]
    else:
        txt = [f"N = {N}, below --min-n", "no fit"]
    nout = int(n_all) - int(h[sel].sum())
    if nout > 0:
        txt.append(f"{nout} outside plot range")
    if note:
        txt.append(note)
    ax.text(0.97, 0.96, "\n".join(txt), transform=ax.transAxes, ha="right", va="top", fontsize=8,
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85))
    ax.set_xlim(x0, x1)
    ax.set_ylim(bottom=0)


def module_page(pdf, ref, geo_row, hists, fits, nall, vlo, min_n, ip, seed_tags=(), n_seed=None):
    plt = _mpl()
    fig, axs = plt.subplots(2, 2, figsize=(11, 8.5))
    nm, mid, x, y = geo_row
    fig.suptitle(f"{nm}   id {mid}   x = {x:.2f} mm, y = {y:.2f} mm", fontsize=13, weight="bold")
    for ax, tag in zip(axs.ravel(), TAGS):
        note = None
        if ref == "veto" and tag in seed_tags:
            xl, lo, msg = "t_module - t_seed (ns), no veto", vlo[tag], "no data"
            ns, nt = int(n_seed[tag][ip]), int(nall[tag][ip])
            if nt:
                note = f"module = seed in {100 * ns / nt:.1f} %\n(those: t_seed - t_2nd)"
        elif ref == "veto":
            xl, lo = "t_module - t_veto (ns)", vlo[tag]
            msg = "no veto readout in X17 runs" if tag == "x17" else "no data"
        else:
            xl, lo = "t_module - t_ref,evt (ns)", EVT_LO
            msg = "no data"
        draw_panel(ax, tag, hists[tag][ip], fits[tag][ip] if fits[tag] is not None else None,
                   nall[tag][ip], lo, xl, min_n, msg, note)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    pdf.savefig(fig)
    plt.close(fig)


def summary_vs_id(pdf, ref_label, geo_p, tab, tags):
    plt = _mpl()
    fig, axs = plt.subplots(2, 1, figsize=(11, 8.5), sharex=True)
    marks = {"0p73": ("o", "C0"), "2p24": ("s", "C3"), "3p49": ("^", "C2"), "x17": ("D", "C1")}
    for ax, key, lab in zip(axs, ("mean", "sigma"), ("mean of dt (ns)", "sigma of dt (ns)")):
        allv = []
        for tag in tags:
            v = tab[tag][key]
            ok = np.isfinite(v)
            m, c = marks[tag]
            ax.plot(geo_p["id"][ok], v[ok], m, ms=3, color=c, alpha=0.75, label=f"{LABEL[tag]} ({ok.sum()})")
            allv.append(v[ok])
        allv = np.concatenate(allv) if allv else np.array([0.0, 1.0])
        if allv.size:
            a, b = pct(allv, 0.5), pct(allv, 99.5)
            pad = 0.1 * (b - a if b > a else 1.0)
            ax.set_ylim(a - pad, b + pad)
        ax.set_ylabel(lab)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=4, frameon=False)
    axs[1].set_xlabel("module id")
    fig.suptitle(f"Per-module timing vs module id - {ref_label}", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    pdf.savefig(fig)
    plt.close(fig)


def summary_maps(pdf, ref_label, geo_p, tab, tags, key):
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Rectangle
    plt = _mpl()
    fig, axs = plt.subplots(2, 2, figsize=(11, 8.5))
    cmap = plt.get_cmap("viridis" if key == "mean" else "plasma")
    sq = [Rectangle((geo_p["x"][i] - geo_p["sx"][i] / 2, geo_p["y"][i] - geo_p["sy"][i] / 2),
                    geo_p["sx"][i], geo_p["sy"][i]) for i in range(len(geo_p["x"]))]
    ext = 1.02 * max(np.abs(geo_p["x"]).max() + 12, np.abs(geo_p["y"]).max() + 12)
    for ax, tag in zip(axs.ravel(), TAGS):
        v = tab[tag][key] if tag in tags else np.full(len(sq), np.nan)
        ok = np.isfinite(v)
        ax.add_collection(PatchCollection([sq[i] for i in np.flatnonzero(~ok)], facecolor="0.85",
                                          edgecolor="0.6", lw=0.3))
        if ok.any():
            vmin, vmax = pct(v[ok], 2), pct(v[ok], 98)
            if not vmax > vmin:
                vmax = vmin + 1e-3
            pc = PatchCollection([sq[i] for i in np.flatnonzero(ok)], cmap=cmap, edgecolor="0.3", lw=0.2)
            pc.set_array(v[ok])
            pc.set_clim(vmin, vmax)
            ax.add_collection(pc)
            cb = fig.colorbar(pc, ax=ax, shrink=0.8)
            cb.set_label(("mean" if key == "mean" else "sigma") + " (ns), colour clipped 2-98 %", fontsize=7)
            cb.ax.tick_params(labelsize=7)
        else:
            ax.text(0.5, 0.5, "no fitted modules", transform=ax.transAxes, ha="center", fontsize=9)
        ax.set(xlim=(-ext, ext), ylim=(-ext, ext), aspect="equal")
        ax.set_title(f"{LABEL[tag]}  ({ok.sum()} fitted, grey = no fit)", fontsize=9)
        ax.tick_params(labelsize=7)
    fig.suptitle(f"Module map of {key} of dt - {ref_label}", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    pdf.savefig(fig)
    plt.close(fig)


# ---- veto summary -----------------------------------------------------------

VETO_COLS = ["run_type", "run", "n_phys", "n_tevt", "f_any", "f_gt500", "f_gt2000", "f_intime",
             *[f"f_any_V{c}" for c in range(1, 5)], *[f"f_gt2000_V{c}" for c in range(1, 5)],
             *[f"f_intime_V{c}" for c in range(1, 5)], "d0_type", "sigma_type", "d0_run",
             "sigma_run", "veto_usable", "n_modules_Ngt_minN"]


def veto_rows(p1, p2, d0, d0fit, d0run, d0_hist, p2_mods, runs_of):
    """p1/p2: per-run dicts of summed stats.  Returns list of dict rows."""
    rows = []
    for tag in TAGS:
        if not runs_of.get(tag):
            continue
        for run in list(runs_of[tag]) + ["ALL"]:
            if run == "ALL":
                a = _sum_dicts([p1[r] for r in runs_of[tag]])
                b = _sum_dicts([p2[r] for r in runs_of[tag]])
                fit = d0fit[tag]
                d0r, sr = fit["mu"], fit["sigma"]
                usable = f"{sum(np.isfinite(d0_hist[r]) for r in runs_of[tag])}/{len(runs_of[tag])}"
                nmod = p2_mods.get(tag, "")
            else:
                a, b = p1[run], p2[run]
                fit = d0run[run] or dict(ok=False, mu=np.nan, sigma=np.nan)
                usable = "yes" if np.isfinite(d0_hist[run]) else "no"
                d0r, sr = (fit["mu"], fit["sigma"]) if fit["ok"] else (np.nan, np.nan)
                nmod = ""
            n = max(a["n_sel"], 1)
            r = dict(run_type=LABEL[tag], run=run, n_phys=a["n_sel"], n_tevt=a["n_tevt"],
                     f_any=a["n_any"] / n, f_gt500=a["n_500"] / n, f_gt2000=a["n_2000"] / n,
                     f_intime=b["n_intime"] / n, d0_type=d0[tag], sigma_type=d0fit[tag]["sigma"],
                     d0_run=d0r, sigma_run=sr, veto_usable=usable, n_modules_Ngt_minN=nmod)
            for c in range(4):
                r[f"f_any_V{c + 1}"] = a["ch_any"][c] / n
                r[f"f_gt2000_V{c + 1}"] = a["ch_2000"][c] / n
                r[f"f_intime_V{c + 1}"] = b["ch_intime"][c] / n
            rows.append(r)
    return rows


def _sum_dicts(ds):
    out = {}
    for k in ds[0]:
        if k == "run":
            continue
        v = [d[k] for d in ds]
        out[k] = np.sum(v, axis=0) if isinstance(v[0], (list, np.ndarray)) else sum(v)
    for k in ("ch_any", "ch_2000", "ch_intime"):
        if k in out:
            out[k] = list(out[k])
    return out


def write_veto_tsv(path, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(VETO_COLS) + "\n")
        for r in rows:
            fh.write("\t".join(
                (f"{r[c]:.5g}" if isinstance(r[c], (float, np.floating)) else str(r[c]))
                for c in VETO_COLS) + "\n")


def veto_pages(pdf, rows, p1, d0, d0fit, runs_of, args):
    plt = _mpl()
    w = lambda v: f"{100 * v:6.2f}"
    h1 = (f"{'run type':<12}{'run':>7}{'n_phys':>9}{'tevt ok':>9}{'any%':>7}{'>500%':>7}"
          f"{'>2000%':>8}{'in-time%':>9}{'D0 run':>8}{'sig run':>8}{'D0 type':>8}{'sig typ':>8}{'usable':>7}{'N>minN':>7}")
    L1 = [h1, "-" * len(h1)]
    L2 = [f"{'run type':<12}{'run':>7}" + "".join(f"{'any V' + str(c) + '%':>9}" for c in range(1, 5))
          + "".join(f"{'>2k V' + str(c) + '%':>9}" for c in range(1, 5))
          + "".join(f"{'it V' + str(c) + '%':>9}" for c in range(1, 5))]
    L2.append("-" * len(L2[0]))
    for r in rows:
        run = r["run"] if r["run"] == "ALL" else str(r["run"])
        L1.append(f"{r['run_type']:<12}{run:>7}{r['n_phys']:>9}{r['n_tevt']:>9}{w(r['f_any']):>7}"
                  f"{w(r['f_gt500']):>7}{w(r['f_gt2000']):>8}{w(r['f_intime']):>9}"
                  f"{f(r['d0_run'], 2):>8}{f(r['sigma_run'], 2):>8}{f(r['d0_type'], 2):>8}"
                  f"{f(r['sigma_type'], 2):>8}{r['veto_usable']:>7}{str(r['n_modules_Ngt_minN']):>7}")
        L2.append(f"{r['run_type']:<12}{run:>7}"
                  + "".join(f"{w(r[f'f_any_V{c}']):>9}" for c in range(1, 5))
                  + "".join(f"{w(r[f'f_gt2000_V{c}']):>9}" for c in range(1, 5))
                  + "".join(f"{w(r[f'f_intime_V{c}']):>9}" for c in range(1, 5)))
        if run == "ALL":
            L1.append("")
            L2.append("")
    note = ["", "Fractions are of physics events (event selection).  'any' = any veto pulse, any integral;",
            f"'>2000' = integral > 2000; 'in-time' = integral > {args.veto_min_int:g} and "
            f"|(v_t - t_evt) - D0| < {args.veto_gate:g} ns.",
            "D0 [ns] = position of the coincidence peak of v_t - t_evt (Gaussian+const core fit), sigma its width;",
            "'run' columns from the run's own histogram, 'type' columns from the pooled type (the gate is centred on D0 of the type);",
            "'usable' = run used for the veto-referenced histograms (significant coincidence peak within the gate of the type D0);",
            "the peak search covers |v_t - t_evt| < 100 ns.",
            "X17 Ta runs: no veto readout (crate 3 slot 3 absent).  N>minN = modules with a veto-referenced",
            f"histogram of more than {args.min_n} entries."]
    text_page(pdf, L1 + note, "Veto availability (1/2): event fractions and coincidence peak", size=7.5)
    text_page(pdf, L2 + note[:3], "Veto availability (2/2): per-channel fractions", size=7.5)
    fig, axs = plt.subplots(2, 4, figsize=(11, 8.5))
    xa = centres(A_LO, A_W, A_NB)
    xd = centres(D_LO, D_W, D_NB)
    for j, tag in enumerate(TAGS):
        ra = [r for r in runs_of.get(tag, [])]
        for i in range(2):
            axs[i, j].set_title(LABEL[tag], fontsize=9)
            axs[i, j].tick_params(labelsize=7)
        if tag == "x17" or not ra:
            for i in range(2):
                axs[i, j].set_xticks([])
                axs[i, j].set_yticks([])
                axs[i, j].text(0.5, 0.5, "no veto readout in X17 runs" if tag == "x17" else "no data",
                               transform=axs[i, j].transAxes, ha="center", fontsize=8, color="0.35")
            continue
        ha = sum(p1[r]["h_abs"] for r in ra)
        hd = sum(p1[r]["h_d"] for r in ra)
        axs[0, j].stairs(ha, np.r_[xa - A_W / 2, xa[-1] + A_W / 2], color="C0", lw=0.8)
        axs[0, j].set_yscale("log")
        axs[0, j].set_xlabel(f"v_t (ns), int > {args.veto_min_int:g}", fontsize=8)
        axs[1, j].stairs(hd, np.r_[xd - D_W / 2, xd[-1] + D_W / 2], color="C0", lw=0.8)
        axs[1, j].set_yscale("log")
        axs[1, j].set_xlabel("v_t - t_evt (ns)", fontsize=8)
        axs[1, j].set_xlim(-200, 200)
        if np.isfinite(d0[tag]):
            axs[1, j].axvspan(d0[tag] - args.veto_gate, d0[tag] + args.veto_gate, color="crimson", alpha=0.15)
            axs[1, j].axvline(d0[tag], color="crimson", lw=0.8)
            axs[1, j].text(0.04, 0.96, f"D0 = {d0[tag]:.2f} ns\nsigma = {d0fit[tag]['sigma']:.2f} ns\n"
                           f"gate +-{args.veto_gate:g} ns", transform=axs[1, j].transAxes, va="top", fontsize=7)
    axs[0, 0].set_ylabel("pulses / 1 ns", fontsize=8)
    axs[1, 0].set_ylabel("pulses / 0.5 ns", fontsize=8)
    fig.suptitle("Veto pulse times (all V1-V4, selected events): absolute and relative to the event time",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    pdf.savefig(fig)
    plt.close(fig)


def seed_summary(tag, q, names, dev=3.0):
    """No-veto run type (t_m - t_seed): median fitted sigma, modules with |mean - median| > dev ns."""
    ok = q["ok"]
    if not ok.sum():
        return [f"  {LABEL[tag]} t_module - t_seed: no fitted module"]
    mm = float(np.median(q["mean"][ok]))
    far = [(str(names[i]), q["mean"][i]) for i in np.flatnonzero(ok) if abs(q["mean"][i] - mm) > dev]
    L = [f"  {LABEL[tag]} t_module - t_seed (no veto): {int(ok.sum())} fitted modules, median sigma "
         f"{np.median(q['sigma'][ok]):.3f} ns, median mean {mm:.3f} ns; |mean - median| > {dev:g} ns: {len(far)}"]
    for k in range(0, len(far), 8):
        L.append("    " + "  ".join(f"{n} {m:+.2f}" for n, m in far[k:k + 8]))
    return L


# ---- main -------------------------------------------------------------------

def nogate_stats(h, d0):
    """All-pairs dt histogram (NG binning) of one module: in-time peak (1 ns smoothed argmax within
    +-NG_PK of -D0), its height (per ns), the accidental level (per ns) from the sidebands
    NG_SB[0] < |dt - peak| < NG_SB[1], and S, B, S/B within +-NG_SIG of the peak."""
    x = centres(NG_LO, NG_W, NG_NB)
    r = dict(N=int(h.sum()), peak=np.nan, height=np.nan, level=np.nan, S=np.nan, B=np.nan, SB=np.nan)
    w = np.abs(x + d0) < NG_PK
    if r["N"] == 0 or not w.any():
        return r
    sm = np.convolve(h, np.ones(3) / 3, mode="same")
    i = np.flatnonzero(w)[int(np.argmax(sm[w]))]
    pk = x[i]
    d = np.abs(x - pk)
    sb = (d > NG_SB[0]) & (d < NG_SB[1])
    lev = h[sb].sum() / (sb.sum() * NG_W)
    sig = d < NG_SIG
    B = lev * sig.sum() * NG_W
    S = h[sig].sum() - B
    r.update(peak=float(pk), height=float(sm[i] / NG_W), level=float(lev), S=float(S), B=float(B),
             SB=float(S / B) if B > 0 else np.inf)
    return r


def nogate_output(args, gp, page_mods, H, T, fits, vlo, seed_tags, d0, d0_hist, runs_of, p1, cinfo,
                  pdf_path, tsv_path):
    """--veto-gate none: module pages of the all-pairs dt histograms (X17: t_m - t_seed as in
    module_timing.pdf), front page with the pooled per-type numbers, per-panel TSV."""
    from matplotlib.backends.backend_pdf import PdfPages
    plt = _mpl()
    vt_tags = [t for t in TAGS if t not in seed_tags and any(np.isfinite(d0_hist[r]) for r in runs_of[t])]
    x = centres(NG_LO, NG_W, NG_NB)
    st, xr, pool = {}, {}, {}
    npw = len(gp["id"])
    for t in vt_tags:
        st[t] = [nogate_stats(H[t]["ng"][i], d0[t]) for i in range(npw)]
        hp = H[t]["ng"].sum(axis=0)
        c = np.cumsum(hp) / max(hp.sum(), 1)
        xr[t] = (x[np.searchsorted(c, 0.005)] - 10.0, x[min(np.searchsorted(c, 0.995), NG_NB - 1)] + 10.0)
        n_ev = sum(p1[r]["n_sel"] for r in runs_of[t] if np.isfinite(d0_hist[r]))
        ok = [q for q in st[t] if q["N"] > 0 and np.isfinite(q["S"])]
        npair = sum(q["N"] for q in ok)
        S, B = sum(q["S"] for q in ok), sum(q["B"] for q in ok)
        pool[t] = dict(n_ev=n_ev, ppe=npair / max(n_ev, 1), lev=sum(q["level"] for q in ok) / max(n_ev, 1),
                       fin=S / max(npair, 1), SB=S / B if B > 0 else np.inf, nmod=len(ok))
    with open(tsv_path, "w") as fh:
        fh.write("module\tid\trun_type\tN_pairs\tpeak_ns\tpeak_height_per_ns\taccidental_per_ns\t"
                 "S_10ns\tB_10ns\tS_over_B_10ns\n")
        for i in page_mods:
            for t in vt_tags:
                q = st[t][i]
                fh.write("\t".join([str(gp["name"][i]), str(gp["id"][i]), LABEL[t], str(q["N"]), f(q["peak"], 2),
                                    f(q["height"], 2), f(q["level"], 3), f(q["S"], 1), f(q["B"], 1),
                                    f(q["SB"], 3)]) + "\n")
    hdr = f"{'run type':<13}{'events':>10}{'pairs/event':>13}{'accid./ns/event':>17}{'in-time frac.':>15}{'S/B (+-10 ns)':>15}{'modules':>9}"
    L = [f"dt = t_module - t_veto' for EVERY pair (dominant PbWO4 pulse of a module) x (V1-V4 pulse with integral > "
         f"{args.veto_min_int:g}) of a selected event:",
         f"no coincidence gate, no largest-integral choice.  Veto times aligned per counter (--veto-align "
         f"{args.veto_align}), shifts V1..V4 [ns]:",
         *[f"    {LABEL[t]}: " + " ".join(f"{v:+.2f}" for v in cinfo[t]["applied"]) for t in vt_tags],
         f"Runs: the veto-usable runs of plot_module_timing.py (" + "; ".join(f"{LABEL[t]} "
         + ",".join(str(r) for r in runs_of[t] if np.isfinite(d0_hist[r])) for t in vt_tags) + ").",
         f"Histogram {NG_W:g} ns bins; panel range = 0.5-99.5 % of the pooled pairs of the run type +- 10 ns.",
         f"Per panel: in-time peak = 1.5 ns smoothed maximum within +-{NG_PK:g} ns of -D0; accidental level = mean "
         f"counts/ns in {NG_SB[0]:g} < |dt - peak| < {NG_SB[1]:g} ns",
         f"(dashed line); S = counts in the bins within +-{NG_SIG:g} ns of the peak minus B = level x their width; S/B in that window.",
         "Pooled numbers: sums over all PbWO4 modules with pairs; per event = per selected event of the usable runs;",
         "in-time frac. = S / N_pairs.  X17 Ta (no veto readout): t_module - t_seed as in module_timing.pdf.",
         f"Pages: the same {len(page_mods)} modules, same order, as module_timing.pdf.", "", hdr, "-" * len(hdr)]
    for t in vt_tags:
        q = pool[t]
        L.append(f"{LABEL[t]:<13}{q['n_ev']:>10}{q['ppe']:>13.3f}{q['lev']:>17.5f}{q['fin']:>15.4f}{q['SB']:>15.3f}{q['nmod']:>9}")
    with PdfPages(pdf_path) as pdf:
        text_page(pdf, L, "Module timing without the veto coincidence gate: all module x veto pairs", size=8.5)
        for ip in page_mods:
            fig, axs = plt.subplots(2, 2, figsize=(11, 8.5))
            fig.suptitle(f"{gp['name'][ip]}   id {gp['id'][ip]}   x = {gp['x'][ip]:.2f} mm, y = {gp['y'][ip]:.2f} mm",
                         fontsize=13, weight="bold")
            for ax, tag in zip(axs.ravel(), TAGS):
                if tag in seed_tags:
                    ns, nt = int(H[tag]["n_seed"][ip]), int(T[tag]["n_all"][ip])
                    draw_panel(ax, tag, H[tag]["veto"][ip], fits[tag][ip] if fits[tag] is not None else None,
                               T[tag]["n_all"][ip], vlo[tag], "t_module - t_seed (ns), no veto", args.min_n, "no data",
                               f"module = seed in {100 * ns / nt:.1f} %\n(those: t_seed - t_2nd)" if nt else None)
                    continue
                ax.set_title(LABEL[tag], fontsize=10)
                ax.tick_params(labelsize=8)
                q = st[tag][ip] if tag in st else None
                if q is None or q["N"] == 0:
                    ax.set_xticks([])
                    ax.set_yticks([])
                    ax.text(0.5, 0.5, "no data", transform=ax.transAxes, ha="center", va="center", color="0.35")
                    continue
                sel = (x >= xr[tag][0]) & (x <= xr[tag][1])
                ax.stairs(H[tag]["ng"][ip][sel], np.r_[x[sel] - NG_W / 2, x[sel][-1] + NG_W / 2], color="C0", lw=0.7)
                txt = [f"N_pairs = {q['N']}"]
                if np.isfinite(q["peak"]):
                    ax.axhline(q["level"] * NG_W, color="crimson", ls="--", lw=1.0)
                    txt += [f"in-time peak at {q['peak']:.1f} ns", f"peak height {q['height']:.1f} /ns",
                            f"accidentals {q['level']:.2f} /ns (sidebands)", f"S/B (+-{NG_SIG:g} ns) = {q['SB']:.2f}"]
                ax.text(0.97, 0.96, "\n".join(txt), transform=ax.transAxes, ha="right", va="top", fontsize=8,
                        bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85))
                ax.set_xlim(*xr[tag])
                ax.set_ylim(bottom=0)
                ax.set_xlabel("t_module - t_veto (ns), all pairs, no gate", fontsize=8)
                ax.set_ylabel(f"pairs / {NG_W:g} ns", fontsize=8)
            fig.tight_layout(rect=(0, 0, 1, 0.96))
            pdf.savefig(fig)
            plt.close(fig)
    print("\n".join(L[-len(vt_tags) - 2:]), flush=True)


def veto_windows(runs_of, d0_hist, seed_tags):
    """Histogram range (lo, number of W_BIN bins) of the 'veto' histogram per run type."""
    vlo, vnb = {}, {}
    for tag in TAGS:
        if tag in seed_tags:
            vlo[tag], vnb[tag] = EVT_LO, EVT_NB
            continue
        cs = [-d0_hist[r] for r in runs_of[tag] if np.isfinite(d0_hist[r])] or [-DEFAULT_D0]
        vlo[tag] = float(np.floor(min(cs))) - VETO_HALF
        vnb[tag] = int(round((np.ceil(max(cs)) + VETO_HALF - vlo[tag]) / W_BIN))
    return vlo, vnb


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", choices=("veto", "evt"), default="veto")
    ap.add_argument("--dump-dir", default=str(HERE / "output/module_timing_20261025/dump"),
                    help="directory with <run>_<split>.npz")
    ap.add_argument("--out-dir", default=str(HERE / "output/module_timing_20261025"))
    ap.add_argument("--trigger-mask", type=lambda s: int(s, 0), default=0xF00)
    ap.add_argument("--min-e", type=float, default=0.0, help="module pulse energy cut (MeV, 0 = none)")
    ap.add_argument("--ref-min-e", type=float, default=10.0)
    ap.add_argument("--ref-min-n", type=int, default=3)
    ap.add_argument("--veto-min-int", type=float, default=2000.0)
    ap.add_argument("--veto-gate", type=lambda x: None if x.lower() == "none" else float(x), default=10.0,
                    help="in-time gate half width (ns); 'none' (with --ref veto): write only "
                         "module_timing_nogate.pdf/.tsv (all module x veto-pulse pairs)")
    ap.add_argument("--min-n", type=int, default=10000)
    ap.add_argument("--time-calib", action="store_true",
                    help="subtract per-module offsets of time_calib.json from all PbWO4 times")
    ap.add_argument("--time-calib-file", default=str(TIME_CALIB))
    ap.add_argument("--veto-align", choices=("counter", "centroid", "type"), default="counter",
                    help="veto counter time alignment before the pick: counter = module-referenced "
                         "shifts (default), centroid = v_t - t_evt centroid of each counter to D0, "
                         "type = none (old behaviour)")
    ap.add_argument("--geometry-run", type=int, default=25320)
    ap.add_argument("--jobs", type=int, default=NMAX_JOBS, help=f"worker processes (max {NMAX_JOBS})")
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    t_start = time.time()
    nproc = max(1, min(args.jobs, NMAX_JOBS))
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    suf = "_tcal" if args.time_calib else ""
    rname = "veto" if args.ref == "veto" else "evtref"
    pdf_path = out_dir / (f"module_timing{suf}.pdf" if args.ref == "veto" else f"module_timing_evtref{suf}.pdf")
    tsv_path = out_dir / f"module_timing_{rname}{suf}.tsv"
    sum_path = out_dir / f"module_timing_{rname}{suf}_summary.txt"
    if args.veto_align != "counter":
        suf += f"_{args.veto_align}align"
        if args.ref == "veto":
            pdf_path = out_dir / f"module_timing{suf}.pdf"
            tsv_path = out_dir / f"module_timing_{rname}{suf}.tsv"
            sum_path = out_dir / f"module_timing_{rname}{suf}_summary.txt"
    veto_tsv = out_dir / f"veto_summary{suf}.tsv"
    nogate = args.veto_gate is None
    if nogate:
        if args.ref != "veto":
            sys.exit("--veto-gate none needs --ref veto")
        args.veto_gate = 10.0      # still the gate of the in-time statistics and of the run selection
        pdf_path = out_dir / f"module_timing_nogate{suf}.pdf"
        tsv_path = out_dir / f"module_timing_nogate{suf}.tsv"

    files = sorted(glob.glob(os.path.join(args.dump_dir, "*.npz")))
    jobs_files, runs_of, splits_of = [], {t: [] for t in TAGS}, {t: 0 for t in TAGS}
    for fp in files:
        mt = re.search(RUN_NAME, fp)
        run = int(mt.group(1)) if mt else None
        if run not in TYPE_OF_RUN:
            print(f"skipping {fp}: run not in the run-type table")
            continue
        jobs_files.append(fp)
        tag = TYPE_OF_RUN[run]
        if run not in runs_of[tag]:
            runs_of[tag].append(run)
        splits_of[tag] += 1
    if not jobs_files:
        sys.exit(f"no usable npz in {args.dump_dir}")
    tags_present = [t for t in TAGS if runs_of[t]]
    for t in TAGS:
        runs_of[t].sort()
    print(f"{len(jobs_files)} dump files, types present: {[LABEL[t] for t in tags_present]}", flush=True)

    # geometry, module mapping, optional time calibration
    geo = load_geometry(args.geometry_run)
    d0f = np.load(jobs_files[0])
    names = d0f["module_names"]
    if [str(n) for n in names] != [str(n) for n in geo["name"]]:
        sys.exit("module_names of the dump differ from the hycal geometry")
    pidx = np.flatnonzero(geo["pwo4"])
    npw = len(pidx)
    pmap = np.full(len(names), -1, np.int64)
    pmap[pidx] = np.arange(npw)
    gp = {k: geo[k][pidx] for k in ("name", "id", "x", "y", "sx", "sy")}
    off, n_res, missing = (None, 0, [])
    if args.time_calib:
        off, n_res, missing = load_offsets(args.time_calib_file, geo)
        print(f"time calibration: {n_res}/{npw} PbWO4 modules resolved"
              + (f"; offset 0 for {missing}" if missing else ""), flush=True)
    cfg = dict(mask=args.trigger_mask, min_e=args.min_e, ref_min_e=args.ref_min_e,
               ref_min_n=args.ref_min_n, veto_min_int=args.veto_min_int, gate=args.veto_gate,
               off=off, npw=npw, pmap=pmap, nogate=nogate,
               names_hash=hashlib.md5("\n".join(map(str, names)).encode()).hexdigest())

    # ---- pass 1: veto availability and D0 ----
    t0 = time.time()
    r1 = _run_pool(pass1, [(fp, cfg) for fp in jobs_files], nproc, "pass1")
    p1 = {}
    for r in r1:
        p1.setdefault(r["run"], []).append(r)
    p1 = {run: _sum_dicts(v) for run, v in p1.items()}
    for run, v in p1.items():
        v["run"] = run
    print(f"pass 1 done in {time.time() - t0:.0f} s", flush=True)
    d0, d0fit, d0run, d0_stat, d0_hist = veto_d0({r: p1[r]["h_d"] for r in p1}, runs_of, args.veto_gate)
    voff, cinfo = measure_counter_offsets(jobs_files, cfg, runs_of, TYPE_OF_RUN, d0, d0_hist,
                                         {r: p1[r]["h_dc"] for r in p1}, args.veto_align, nproc)
    for t in TAGS:
        if cinfo[t]["runs_used"]:
            print(f"  veto counter shifts [{LABEL[t]}] applied ({args.veto_align}): "
                  + " ".join(f"V{c} {v:+.2f}" for c, v in enumerate(cinfo[t]["applied"], 1))
                  + f";  v_t - t_evt centroids: " + " ".join(
                      f"{k} {f(v['centroid'] if v['centroid'] is not None else np.nan, 2)}"
                      for k, v in cinfo[t]["vt_minus_tevt_centroid"].items())
                  + f";  module-referenced residual rms {cinfo[t]['module_referenced'].get('residual_rms_ns')} ns "
                  f"({cinfo[t]['module_referenced'].get('n_modules')} modules)", flush=True)
    seed_tags = [t for t in TAGS if args.ref == "veto" and runs_of[t]
                 and not any(np.isfinite(d0_hist[r]) for r in runs_of[t])]
    vlo, vnb = veto_windows(runs_of, d0_hist, seed_tags)

    # ---- pass 2: per-module histograms ----
    t0 = time.time()
    run_of = lambda fp: int(re.search(RUN_NAME, fp).group(1))
    jobs2 = [(fp, cfg, d0_stat[run_of(fp)], d0_hist[run_of(fp)], vlo[TYPE_OF_RUN[run_of(fp)]],
              vnb[TYPE_OF_RUN[run_of(fp)]], voff[TYPE_OF_RUN[run_of(fp)]], TYPE_OF_RUN[run_of(fp)] in seed_tags)
             for fp in jobs_files]
    H = {t: dict(evt=np.zeros((npw, EVT_NB), np.int64), veto=np.zeros((npw, vnb[t]), np.int64),
                 n_evt=np.zeros(npw, np.int64), n_veto=np.zeros(npw, np.int64),
                 n_seed=np.zeros(npw, np.int64), ng=np.zeros((npw, NG_NB), np.int64)) for t in TAGS}
    p2 = {}
    for r in _run_pool(pass2, jobs2, nproc, "pass2"):
        tag = TYPE_OF_RUN[r["run"]]
        H[tag]["evt"] += r["h_evt"]
        H[tag]["veto"] += r["h_veto"]
        H[tag]["n_evt"] += r["n_evt"]
        H[tag]["n_veto"] += r["n_veto"]
        if "n_seed" in r:
            H[tag]["n_seed"] += r["n_seed"]
        if "h_ng" in r:
            H[tag]["ng"] += r["h_ng"]
        p2.setdefault(r["run"], []).append(dict(n_intime=r["n_intime"], ch_intime=r["ch_intime"]))
    p2 = {run: _sum_dicts(v) for run, v in p2.items()}
    print(f"pass 2 done in {time.time() - t0:.0f} s", flush=True)

    # ---- fits ----
    t0 = time.time()
    xm = {"evt": centres(EVT_LO, W_BIN, EVT_NB)}
    tab = {"evt": {}, "veto": {}}
    fits = {"evt": {}, "veto": {}}
    for ref in ("evt", "veto"):
        for tag in TAGS:
            h = H[tag][ref]
            x = centres(EVT_LO, W_BIN, EVT_NB) if ref == "evt" else centres(vlo[tag], W_BIN, vnb[tag])
            Ns = h.sum(axis=1)
            fl = [None] * npw
            for i in np.flatnonzero(Ns > args.min_n):
                fl[i] = fit_core(x, h[i])
            fits[ref][tag] = fl
            g = lambda key: np.array([r[key] if (r is not None and r["ok"]) else np.nan for r in fl])
            tab[ref][tag] = dict(N=Ns, n_all=H[tag]["n_" + ref], mean=g("mu"), sigma=g("sigma"),
                                 ok=np.array([bool(r is not None and r["ok"]) for r in fl]))
    n_veto_mods = {t: int((tab["veto"][t]["N"] > args.min_n).sum()) for t in TAGS
                   if runs_of[t] and t not in seed_tags}
    print(f"fits done in {time.time() - t0:.0f} s", flush=True)

    if nogate:
        T = tab["veto"]
        pm = np.flatnonzero(np.any([T[t]["N"] > args.min_n for t in TAGS], axis=0))
        nogate_output(args, gp, pm[np.argsort(gp["id"][pm])], H, T, fits["veto"], vlo, seed_tags, d0, d0_hist,
                      runs_of, p1, cinfo, pdf_path, tsv_path)
        print(f"total {time.time() - t_start:.0f} s; outputs: {pdf_path}, {tsv_path}")
        return

    # ---- veto summary ----
    vrows = veto_rows(p1, p2, d0, d0fit, d0run, d0_hist, n_veto_mods, runs_of)
    write_veto_tsv(veto_tsv, vrows)

    # ---- selected reference ----
    ref = args.ref
    T = tab[ref]
    ref_hist = {t: H[t][ref] for t in TAGS}
    ref_nall = {t: T[t]["n_all"] for t in TAGS}
    ref_label = ("t_veto (in-time V1-V4 pulse)" if ref == "veto" else "t_ref,evt (leave-one-out event median)")
    page_mods = np.flatnonzero(np.any([T[t]["N"] > args.min_n for t in TAGS], axis=0))
    page_mods = page_mods[np.argsort(gp["id"][page_mods])]
    npages_mod = len(page_mods)

    # TSV
    with open(tsv_path, "w") as fh:
        cols = ["module", "id", "x", "y"]
        for t in TAGS:
            cols += [f"N_{t}", f"mean_{t}", f"sigma_{t}", f"fit_ok_{t}"]
        fh.write("\t".join(cols) + "\n")
        for i in range(npw):
            row = [str(gp["name"][i]), str(gp["id"][i]), f"{gp['x'][i]:.3f}", f"{gp['y'][i]:.3f}"]
            for t in TAGS:
                q = T[t]
                row += [str(int(q["N"][i])), f"{q['mean'][i]:.4f}", f"{q['sigma'][i]:.4f}",
                        str(int(q["ok"][i]))]
            fh.write("\t".join(row) + "\n")

    # summary text
    S = [f"module timing summary, reference = {ref} ({ref_label}), time-calib = {'on' if args.time_calib else 'off'}",
         f"min-n = {args.min_n}; modules with a page: {npages_mod}", ""]
    for t in tags_present:
        q = T[t]
        okm = q["ok"]
        S.append(f"[{LABEL[t]}]  modules with N > min-n: {int((q['N'] > args.min_n).sum())}; fitted: {int(okm.sum())}")
        if okm.sum():
            sg, mn = q["sigma"][okm], q["mean"][okm]
            S.append(f"  sigma: median {np.median(sg):.3f} ns, 16-84 % range [{pct(sg, 16):.3f}, {pct(sg, 84):.3f}] ns")
            S.append(f"  mean : median {np.median(mn):.3f} ns, 16-84 % range [{pct(mn, 16):.3f}, {pct(mn, 84):.3f}] ns")
            orow = outliers(gp["name"], q["mean"], q["sigma"], q["N"])
            S.append(f"  outlier modules (|mean - median| > 5 x MAD-sigma, or sigma > 2 x median): {len(orow)}")
            for nm, m, s, n, why in orow:
                S.append(f"    {nm:<6} mean {m:8.3f}  sigma {s:7.3f}  N {int(n):8d}  ({why})")
        if t in seed_tags:
            S += seed_summary(t, q, gp["name"])
        S.append("")
    S.append("veto reference, modules with N > min-n per run type: "
             + ", ".join(f"{LABEL[t]} {n_veto_mods[t]}" for t in tags_present if t in n_veto_mods))
    S.append("D0 per run type [ns]: " + ", ".join(f"{LABEL[t]} {f(d0[t], 2)} (sigma {f(d0fit[t]['sigma'], 2)})"
                                                  for t in tags_present))
    open(sum_path, "w").write("\n".join(S) + "\n")
    print("\n".join(S), flush=True)

    # ---- PDF ----
    t0 = time.time()
    from matplotlib.backends.backend_pdf import PdfPages
    n_read = {t: sum(p1[r]["n_read"] for r in runs_of[t]) for t in TAGS}
    n_dump = {t: sum(p1[r]["n_dump"] for r in runs_of[t]) for t in TAGS}
    n_sel = {t: sum(p1[r]["n_sel"] for r in runs_of[t]) for t in TAGS}
    with PdfPages(pdf_path) as pdf:
        tc_txt = (f"ON, per-module offsets of time_calib.json subtracted from all PbWO4 times "
                  f"({n_res}/{npw} modules resolved)" if args.time_calib else "off (raw times)")
        cover = [
            f"Reference: {ref_label}",
            "",
            "Definition:",
        ]
        if ref == "veto":
            cover += [
                "  dt = t_m - t_veto, t_m = CFD (50 % leading-edge, mf) time of the module's dominant pulse.",
                f"  t_veto = time of the in-time V1-V4 pulse: integral > {args.veto_min_int:g}, and",
                f"  |(v_t - t_evt) - D0| < {args.veto_gate:g} ns; several: largest integral.",
                f"  t_evt = full-event median of the eligible PbWO4 times (m_E > {args.ref_min_e:g} MeV, >= {args.ref_min_n} pulses).",
                "  D0 = coincidence-peak position of v_t - t_evt per run type, from the data (Gaussian+const fit):",
                "    " + ", ".join(f"{LABEL[t]} {f(d0[t], 2)} ns" for t in tags_present),
                "  gate centre: D0 of the pooled run type (per-run D0 listed in veto_summary" + ("_tcal" if args.time_calib else "") + ".tsv)",
                "  runs not used for the veto reference (no coincidence peak at the gate centre): "
                + (", ".join(str(r) for t in tags_present if t != "x17" for r in runs_of[t]
                             if not np.isfinite(d0_hist[r])) or "none"),
                {"counter": "  veto counters aligned (module-referenced): v_t' = v_t - c_k, c_k from the additive fit of the "
                            "per-module, per-counter",
                 "centroid": "  veto counters aligned to D0 by their own v_t - t_evt centroids (--veto-align centroid):",
                 "type": "  veto counters NOT aligned (--veto-align type):"}[args.veto_align],
                *(["    unaligned dt medians (median = a_module + b_counter, sum W_k c_k = 0); shifts c_k [ns]:"]
                  if args.veto_align == "counter" else []),
                *[f"    {LABEL[t]}: " + " ".join(f"V{c} {v:+.2f}" for c, v in enumerate(cinfo[t]["applied"], 1))
                  for t in tags_present if cinfo[t]["runs_used"]],
            ]
            if seed_tags:
                cover += [
                    f"  {', '.join(LABEL[t] for t in seed_tags)}: no veto readout (crate 3 slot 3 absent); panel = "
                    "t_m - t_seed, t_seed = time of the",
                    f"    highest-m_E PbWO4 pulse of the event other than the module (m_E > {args.ref_min_e:g} MeV; "
                    "pulse-level, no clusters); for the",
                    "    event's seed module the second-highest pulse is used (t_seed - t_2nd; fraction printed per panel).",
                ]
                for t in seed_tags:
                    sl = seed_summary(t, T[t], gp["name"])
                    cover += sl[:4] + ([f"    ... {len(sl) - 4} more lines in {sum_path.name}"] if len(sl) > 4 else [])
        else:
            cover += [
                "  dt = t_m - t_ref, t_m = CFD (50 % leading-edge, mf) time of the module's dominant pulse.",
                f"  t_ref = leave-one-out median of t over the OTHER PbWO4 dominant pulses of the same event",
                f"  with m_E > {args.ref_min_e:g} MeV, requiring >= {args.ref_min_n} such pulses (module excluded exactly).",
            ]
        cover += [
            "",
            "Selection:",
            f"  trigger bits & 0x{args.trigger_mask:x};  module pulse m_E >= {args.min_e:g} MeV"
            + (" (every reconstructed dominant pulse)" if args.min_e <= 0 else "")
            + f";  ref-min-e {args.ref_min_e:g} MeV;  ref-min-n {args.ref_min_n}",
            f"  excluded events: LMS flasher / alpha pulser (trig & 0x{TRIG_EXCL:x})",
            f"  veto: integral > {args.veto_min_int:g}, gate +-{args.veto_gate:g} ns around D0",
            f"  time-calib: {tc_txt}",
            "  fit: iterative Gaussian in [mu - 2 sigma, mu + 2 sigma] (refit until mu, sigma change < 1 % of sigma, max 5)",
            "",
            "Data (module_timing_dump.py npz):",
        ]
        for t in tags_present:
            cover.append(f"  {LABEL[t]:<12} runs {', '.join(map(str, runs_of[t]))};  {splits_of[t]} splits;  "
                         f"{n_dump[t]} events read, {n_sel[t]} pass the event selection")
        cover += [
            "",
            f"min-n = {args.min_n} (a module gets a page if N > min-n in at least one run type; "
            f"panels with N <= min-n are not fitted)",
            f"Modules with a page: {npages_mod} of {npw} PbWO4 modules;  dump dir {os.path.relpath(args.dump_dir, HERE)}",
        ]
        text_page(pdf, cover, f"Per-module timing: dt = t_module - {'t_veto' if ref == 'veto' else 't_ref (event median, leave-one-out)'}",
                  size=9)
        if ref == "veto":
            veto_pages(pdf, vrows, p1, d0, d0fit, runs_of, args)
        for ip in page_mods:
            module_page(pdf, ref, (gp["name"][ip], gp["id"][ip], gp["x"][ip], gp["y"][ip]),
                        ref_hist, {t: fits[ref][t] for t in TAGS}, ref_nall, vlo, args.min_n, ip,
                        seed_tags, {t: H[t]["n_seed"] for t in TAGS})
        lab = f"{ref} reference, time-calib {'on' if args.time_calib else 'off'}" + (
            f"; {', '.join(LABEL[t] for t in seed_tags)}: t_m - t_seed" if seed_tags else "")
        summary_vs_id(pdf, lab, gp, T, [t for t in tags_present if t not in seed_tags])
        summary_maps(pdf, lab, gp, T, tags_present, "mean")
        summary_maps(pdf, lab, gp, T, tags_present, "sigma")
    print(f"pdf written in {time.time() - t0:.0f} s: {pdf_path}", flush=True)
    print(f"total {time.time() - t_start:.0f} s; outputs: {pdf_path}, {tsv_path}, {sum_path}, {veto_tsv}")


if __name__ == "__main__":
    main()
