#!/usr/bin/env python3
"""twopulse_fit.py -- two-pulse chi2 template fit: pile-up information ceiling (U0).

Plan: plans/plan_ml_pileup_20261027.md (U0).  Per pedestal-subtracted waveform a
one-pulse fit (A, t) and a two-pulse fit (A1, t1, A2, t2) on the best-50 spline
template of the module, evaluated at continuous time (0.05 ns LUT of the spline,
linear interpolation).  Dense onset grid (0.5 ns, +-SEARCH ns around the
leading-edge estimate), amplitudes solved linearly at every t / (t1, t2) pair
with A >= 0, then scipy least_squares refinement (analytic Jacobian).  Pulse 2 is
always the later one: the decision statistic dchi2 = chi2_1 - chi2_2 and the
search are symmetric in time order; there is no minimum-separation floor
beyond MIN_SEP (default 1 ns).  N = 2 iff dchi2 > threshold, the threshold
being the (1 - fake) quantile of dchi2 on the calibration half (even base
index) of the single-pulse events, applied unchanged to the evaluation half
(odd base index).  Sigma per sample: sigma_i^2 = sig_ped(module)^2
+ (eps * y_i^+)^2 + (kap * |dy_i/dt|)^2, y the measured waveform (fixed weights
per waveform, identical for the 1- and 2-pulse hypothesis).

Synthetic arms (one best-50 split per run type, never pooled):
  real      real clean pulse (+ real noise) + noise-free best-50 template, random order
  tpl_noise best-50 templates (A, t of the real base fit) + measured-noise Gaussian
  nf_same   noise-free best-50 doublets, fitted with best-50
  nf_xbank  noise-free doublets of the Wassim bank (50% rise time aligned), fitted with best-50
The production MFReconstructor (benchmark settings) is run on the `real` arm.

    MPLBACKEND=Agg python3 twopulse_fit.py --out output/ml_pileup_20261027/U0

Outputs: cells.tsv, eff_vs_dt.pdf, run_<type>.npz caches, summary.json (the
numbers report.md is written from).  numpy / scipy only.
"""
from __future__ import annotations

import argparse
import json
import math
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

import mf_pulse

HERE = Path(__file__).resolve().parent
DT = 4.0
NS = 100
TS = DT * np.arange(NS)
LUT_DT = 0.05
PED_N = 30
TOL = 4.0

B50 = HERE / "output/templates_spline_best50_20261026"
WAS = HERE / "output/templates_wassim_20261027"
RUNS = {  # run type -> (slug, one best-50 split not in the production synth)
    "0.73 GeV H2": ("0p73GeV_H2", "024603_00100"),
    "2.24 GeV H2": ("2p24GeV_H2", "025228_00100"),
    "3.49 GeV H2": ("3p49GeV_H2", "024386_00040"),
}
DTS = (4, 8, 12, 16, 20, 24, 32)
RATIOS = (0.3, 0.5, 1.0, 2.0)
ARMS = ("real", "tpl_noise", "nf_same", "nf_xbank")

# fit options (set from the CLI before the pool forks)
OPT = dict(eps=0.01, kap=0.0, search=80.0, grid=0.5, min_sep=1.0)
_BANK = None        # best-50 SampledTemplateBank (fit and MF)


# ------------------------------------------------------------------ banks --
def load_bank(root, slug):
    npz = root / slug / f"templates_{slug}.npz"
    z = np.load(npz, allow_pickle=False)
    names = [str(x) for x in z["names"]]
    wc = dict(smooth_order=2, int_tail_ratio=0.1, tail_break_n=2)
    bank, info = mf_pulse.build_bank(names, np.ones(len(names), bool), str(npz), wc)
    return bank, info


def tpl(row, drow, t, T_end, want_slope=False):
    """Unit-peak template (and slope) at t ns after onset (any shape); 0 outside."""
    x = t * (1.0 / LUT_DT)
    n = row.size
    i = np.clip(np.floor(x).astype(np.int64), 0, n - 2)
    f = x - i
    inside = (t > 0) & (t < T_end)
    v = np.where(inside, row[i] * (1 - f) + row[i + 1] * f, 0.0)
    if not want_slope:
        return v
    d = np.where(inside, drow[i] * (1 - f) + drow[i + 1] * f, 0.0)
    return v, d


# ------------------------------------------------------------------- fit --
_TRIU = {}


def _triu(K, k):
    key = (K, k)
    if key not in _TRIU:
        _TRIU[key] = np.triu(np.ones((K, K), bool), k)
    return _TRIU[key]


def weights(y, sig_p, eps, kap):
    d = np.zeros_like(y)
    d[1:-1] = 0.5 * (y[2:] - y[:-2])
    var = sig_p * sig_p + (eps * np.maximum(y, 0.0)) ** 2 + (kap * d / DT) ** 2
    return 1.0 / var


NFIT = 14   # chi0, chi1g, A, t, chi1, chi2g, dchi_g, chi2, A1, t1, A2, t2, dchi, ok2


def fit_wave(y, row, drow, T_end, tpk, sig_p, o=None):
    """One- and two-pulse fits of the pedestal-subtracted waveform y.

    Returns [chi0, chi1g, A, t, chi1, chi2g, dchi_g, chi2, A1, t1, A2, t2, dchi,
    ok2]; *g are grid minima, the rest refined; t, t1, t2 are onsets (ns) on the
    bank axis, pulse 2 the later; A2 = 0 / t2 = nan for a one-pulse optimum.
    """
    o = o or OPT
    w = weights(y, sig_p, o["eps"], o["kap"])
    sw = np.sqrt(w)
    gstep = o["grid"]
    srch = o["search"]
    ipk = int(np.argmax(y))
    c = DT * ipk - tpk
    c0 = gstep * np.round((c - srch) / gstep)
    K = int(round(2 * srch / gstep)) + 1
    g = c0 + gstep * np.arange(K)
    idx = np.rint((TS[None, :] - g[:, None]) / LUT_DT).astype(np.int64)
    ok = (idx > 0) & (idx < row.size)
    S = np.where(ok, row[np.clip(idx, 0, row.size - 1)], 0.0)
    wy = w * y
    chi0 = float(np.dot(wy, y))
    b = S @ wy
    Gm = (S * w) @ S.T
    gd = np.maximum(np.diag(Gm), 1e-12)
    gain1 = np.where(b > 0, b * b / gd, 0.0)
    k1 = int(np.argmax(gain1))
    out = np.full(NFIT, np.nan)
    out[0] = chi0
    out[1] = chi0 - gain1[k1]
    msi = max(1, int(math.ceil(o["min_sep"] / gstep - 1e-9)))
    g11, g22 = gd[:, None], gd[None, :]
    det = g11 * g22 - Gm * Gm
    n1 = g22 * b[:, None] - Gm * b[None, :]
    n2 = g11 * b[None, :] - Gm * b[:, None]
    valid = (n1 > 0) & (n2 > 0) & (det > 1e-8 * g11 * g22) & _triu(K, msi)
    gain2 = np.where(valid, (n1 * b[:, None] + n2 * b[None, :]) / np.where(valid, det, 1.0), -1.0)
    kk = int(np.argmax(gain2))
    a_i, b_i = divmod(kk, K)
    g2 = gain2[a_i, b_i]
    have2 = g2 > gain1[k1]
    out[5] = chi0 - max(g2, gain1[k1])
    out[6] = max(g2 - gain1[k1], 0.0)
    lo_t, hi_t = g[0] - 2.0, g[-1] + 2.0

    def r1(p):
        s = tpl(row, drow, TS - p[1], T_end)
        return (y - p[0] * s) * sw

    def j1(p):
        s, ds = tpl(row, drow, TS - p[1], T_end, True)
        J = np.empty((NS, 2))
        J[:, 0] = -s * sw
        J[:, 1] = p[0] * ds * sw
        return J

    A0 = max(b[k1] / gd[k1], 1e-3)
    s1 = least_squares(r1, [A0, g[k1]], jac=j1, bounds=([0, lo_t], [np.inf, hi_t]),
                       x_scale=[max(A0, 1.0), 1.0], xtol=1e-8, ftol=1e-10, max_nfev=40)
    chi1 = 2 * s1.cost
    out[2], out[3], out[4] = s1.x[0], s1.x[1], chi1
    out[7] = chi1
    out[8], out[9], out[10], out[11] = s1.x[0], s1.x[1], 0.0, np.nan
    out[12] = 0.0
    out[13] = 0.0
    if have2:
        det0 = det[a_i, b_i]
        x0 = [n1[a_i, b_i] / det0, g[a_i], n2[a_i, b_i] / det0, max(g[b_i] - g[a_i], o["min_sep"] + 1e-6)]

        def r2(p):
            s1_ = tpl(row, drow, TS - p[1], T_end)
            s2_ = tpl(row, drow, TS - p[1] - p[3], T_end)
            return (y - p[0] * s1_ - p[2] * s2_) * sw

        def j2(p):
            s1_, d1 = tpl(row, drow, TS - p[1], T_end, True)
            s2_, d2 = tpl(row, drow, TS - p[1] - p[3], T_end, True)
            J = np.empty((NS, 4))
            J[:, 0] = -s1_ * sw
            J[:, 1] = (p[0] * d1 + p[2] * d2) * sw
            J[:, 2] = -s2_ * sw
            J[:, 3] = p[2] * d2 * sw
            return J

        s2 = least_squares(r2, x0, jac=j2,
                           bounds=([0, lo_t, 0, o["min_sep"]], [np.inf, hi_t, np.inf, 2 * srch + 4]),
                           x_scale=[max(x0[0], 1.0), 1.0, max(x0[2], 1.0), 1.0],
                           xtol=1e-8, ftol=1e-10, max_nfev=60)
        chi2 = 2 * s2.cost
        if chi2 < chi1:
            out[7] = chi2
            out[8], out[9] = s2.x[0], s2.x[1]
            out[10], out[11] = s2.x[2], s2.x[1] + s2.x[3]
            out[12] = chi1 - chi2
            out[13] = 1.0
    return out


def _worker(args):
    Y, m, sigp, opt = args
    b = _BANK
    res = np.empty((len(m), NFIT))
    for i in range(len(m)):
        mi = int(m[i])
        res[i] = fit_wave(Y[i], b.lut[mi].astype(np.float64), b.dlut[mi].astype(np.float64),
                          b.T_end, b.t_pk[mi], sigp[i], opt)
    return res


def fit_batch(pool, Y, m, sigp, chunk=25):
    """Pedestal-subtracted waveforms Y[n, 100] -> fit table [n, NFIT]."""
    n = len(m)
    tasks = [(Y[i:i + chunk], m[i:i + chunk], sigp[i:i + chunk], dict(OPT)) for i in range(0, n, chunk)]
    if not tasks:
        return np.zeros((0, NFIT))
    if pool is None:
        return np.concatenate([_worker(t) for t in tasks])
    return np.concatenate(pool.map(_worker, tasks, chunksize=1))


def subtract_ped(bank, raw):
    """MF pedestal (first 30 samples, clipped mean): y = raw - ped, ped."""
    rec = mf_pulse.MFReconstructor(bank, ped_nsamples=PED_N)
    ped, _ = rec._pedestal(np.asarray(raw, np.float32))
    return np.asarray(raw, np.float64) - ped[:, None].astype(np.float64), ped


# ----------------------------------------------------------------- pools --
def wilson(k, n, z=1.0):
    if n <= 0:
        return (0.0, 1.0)
    ph = k / n
    den = 1.0 + z * z / n
    ctr = (ph + z * z / (2 * n)) / den
    hw = z * math.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / den
    return (max(0.0, ctr - hw), min(1.0, ctr + hw))


class Pool:
    """Clean real pulses of one best-50 split plus module noise models."""

    def __init__(self, slug, split, bank, xbank, base_t_min=150.0):
        p = np.load(B50 / slug / "pulses" / f"{split}.npz", allow_pickle=False)
        self.slug, self.split = slug, split
        ix = {n: i for i, n in enumerate(bank.names)}
        nm = np.asarray(p["name"])
        s = p["samples"].astype(np.float64)
        fb = set(bank.fallback)
        inbank = np.array([(n in ix) and (n not in fb) for n in nm])
        ok = ((p["mf_npk"] == 1) & (p["mf_sec_ratio"] < 0.05) & (s.max(1) < 4095)
              & (p["amp"] >= 100) & (p["amp"] <= 4000) & inbank)
        self.n_split, self.n_clean0 = len(nm), int(ok.sum())
        # module pedestal rms (median over all of the module's stored pulses)
        sig = {}
        for n in np.unique(nm[inbank]):
            sig[n] = max(float(np.median(p["rms"][nm == n])), mf_pulse.SIGMA_MIN)
        self.sig_mod = sig
        # pooled normalised noise autocorrelation from the pedestal window
        x = s[ok][:, :PED_N]
        x = x - x.mean(1, keepdims=True)
        den = (x * x).sum(1)
        good = den > 0
        acf = [1.0]
        for k in range(1, PED_N):
            r = (x[good][:, :-k] * x[good][:, k:]).sum(1) / den[good]
            acf.append(float(np.median(r)) * PED_N / (PED_N - k))   # median, unbiased-ish
        self.acf = np.array(acf)
        cov = np.zeros((NS, NS))
        for i in range(NS):
            for j in range(NS):
                k = abs(i - j)
                cov[i, j] = self.acf[k] if k < PED_N else 0.0
        ev, V = np.linalg.eigh(cov)
        self.noise_L = V * np.sqrt(np.clip(ev, 0, None))[None, :]
        self.raw = s[ok]
        self.name = nm[ok]
        self.midx = np.array([ix[n] for n in self.name], dtype=np.int64)
        self.sigp = np.array([sig[n] for n in self.name])
        self.xmidx = None
        self.event = p["event"][ok]
        self.base_t_min = base_t_min

    def noise(self, rng, m_sig, n):
        """n correlated-Gaussian noise rows (module rms, pre-rounding), [n, NS]."""
        z = rng.standard_normal((NS, n))
        sg = np.sqrt(np.maximum(m_sig ** 2 - 1.0 / 12.0, 1e-4))
        return (self.noise_L @ z).T * sg[:, None]


def real_pool_fit(pool, bank, pw, base_t_min):
    """One+two-pulse fits of every clean real pulse; keeps t_b >= base_t_min."""
    y, ped = subtract_ped(bank, pool.raw)
    R = fit_batch(pw, y, pool.midx, pool.sigp)
    keep = (R[:, 3] >= base_t_min) & (R[:, 2] >= 100) & (R[:, 2] <= 4000)
    return R, y, ped, keep


def slice_pool(pool, R, ped, keep):
    for k in ("raw", "name", "midx", "sigp", "event"):
        setattr(pool, k, getattr(pool, k)[keep])
    pool.ped = ped[keep].astype(np.float64)
    pool.Rb = R[keep]
    pool.A_b, pool.t_b = pool.Rb[:, 2].copy(), pool.Rb[:, 3].copy()
    pool.n = int(keep.sum())
    pool.parity = np.arange(pool.n) % 2          # 0 = calibration half, 1 = evaluation half
    pool.partners = {}
    for h in (0, 1):
        for i in np.nonzero(pool.parity == h)[0]:
            pool.partners.setdefault((h, int(pool.midx[i])), []).append(i)
    pool.partners = {k: np.array(v) for k, v in pool.partners.items()}
    pool.has_partner = np.array([len(pool.partners[(int(pool.parity[i]), int(pool.midx[i]))]) > 1
                                 for i in range(pool.n)])


def pulse_wave(bank, m, A, t, shift=None):
    """A * template(4 i - t) for events (arrays m, A, t) -> [n, NS]; noise-free.

    shift: optional per-module onset shift (cross-bank alignment: the template of
    `bank` is placed at onset t + shift[m]).
    """
    out = np.zeros((len(m), NS))
    for j in range(len(m)):
        tt = t[j] + (shift[m[j]] if shift is not None else 0.0)
        out[j] = A[j] * tpl(bank.lut[m[j]].astype(np.float64), None, TS - tt, bank.T_end)
    return out


def event_table(pool, n, dt_nom, ratio, jit, rng, half=1, partner=False):
    """Truth table of n doublets from evaluation-half bases (random order and phase).

    The base is the real/template pulse with fitted (A_b, t_b); the injected pulse
    sits dt_true = dt_nom + U(-jit, jit) ns before or after it (50/50) with
    amplitude A_b*ratio if it is the later one, A_b/ratio if it is the earlier
    (ratio = later / earlier).  Events with A1 + A2 > 4000 ADC (clip) are redrawn.
    partner (arm real2): the injected pulse is another real clean pulse of the module
    and half whose fitted amplitude is within a factor 2 of the target (one of the 3
    nearest, rescaled), moved by a whole number of samples; its own sub-sample phase
    supplies the random offset, so t_inj / dt_true are set from its fitted onset.
    """
    cand = np.nonzero((pool.parity == half) & (pool.has_partner if partner else True))[0]
    cols = []
    got = 0
    while got < n:
        k = n - got
        idx = rng.choice(cand, k)
        inj_first = rng.random(k) < 0.5
        dtt = np.maximum(dt_nom + rng.uniform(-jit, jit, k), 0.5) if jit > 0 else np.full(k, float(dt_nom))
        Ab, tb = pool.A_b[idx], pool.t_b[idx]
        Ai = np.where(inj_first, Ab / ratio, Ab * ratio)
        ti = np.where(inj_first, tb - dtt, tb + dtt)
        ok = (Ab + Ai) <= 4000.0
        part = np.full(k, -1, np.int64)
        shift = np.zeros(k, np.int64)
        if partner:
            for j in range(k):
                if not ok[j]:
                    continue
                pm = pool.partners[(int(pool.parity[idx[j]]), int(pool.midx[idx[j]]))]
                pm = pm[pm != idx[j]]
                lg = np.abs(np.log(pool.A_b[pm] / Ai[j]))
                pm = pm[np.argsort(lg)[:3]]
                pm = pm[np.abs(np.log(pool.A_b[pm] / Ai[j])) <= math.log(2.0)]
                if pm.size == 0:
                    ok[j] = False
                    continue
                part[j] = rng.choice(pm)
                shift[j] = int(np.rint((ti[j] - pool.t_b[part[j]]) / DT))
                ti[j] = pool.t_b[part[j]] + DT * shift[j]
            dtt = np.where(inj_first, tb - ti, ti - tb)
        cols.append(np.stack([idx, inj_first, dtt, Ai, ti, part, shift], 1)[ok])
        got += int(ok.sum())
    r = np.concatenate(cols)[:n]
    idx = r[:, 0].astype(np.int64)
    inj_first = r[:, 1] > 0.5
    dtt, Ai, ti = r[:, 2], r[:, 3], r[:, 4]
    Ab, tb = pool.A_b[idx], pool.t_b[idx]
    ev = dict(idx=idx, inj_first=inj_first, dt_true=dtt, A_inj=Ai, t_inj=ti,
              A1=np.where(inj_first, Ai, Ab), t1=np.where(inj_first, ti, tb),
              A2=np.where(inj_first, Ab, Ai), t2=np.where(inj_first, tb, ti))
    if partner:
        ev.update(partner=r[:, 5].astype(np.int64), shift=r[:, 6].astype(np.int64))
    return ev


def make_waves(arm, pool, bank, xbank, idx, rng, inj=None, ev=None):
    """Raw waveforms [n, NS] (pedestal included) for base rows idx, optional injection
    inj = (A_inj, t_inj).  real: real base + noise-free best-50 injection; template
    arms: the base is the module's template at the base-fit (A_b, t_b)."""
    m = pool.midx[idx]
    if arm == "real":
        W = pool.raw[idx].copy()
        if inj is not None:
            W += pulse_wave(bank, m, inj[0], inj[1])
        return np.rint(W)
    if arm == "real2":
        W = pool.raw[idx].copy()
        if inj is not None:
            for j in range(len(idx)):
                pj = ev["partner"][j]
                y2 = (pool.raw[pj] - pool.ped[pj]) * (inj[0][j] / pool.A_b[pj])
                k = int(ev["shift"][j])
                if k >= 0:
                    W[j, k:] += y2[:NS - k]
                else:
                    W[j, :NS + k] += y2[-k:]
                    # the partner's last |k| samples lie beyond its window: noise-free template tail
                    ii = np.arange(NS + k, NS)
                    W[j, ii] += inj[0][j] * tpl(bank.lut[pool.midx[idx[j]]].astype(np.float64), None,
                                                DT * (ii - k) - pool.t_b[pj], bank.T_end)
        return np.rint(W)
    if arm == "nf_xbank":
        sh = bank.t50 - xbank.t50
        W = pulse_wave(xbank, m, pool.A_b[idx], pool.t_b[idx], sh)
        if inj is not None:
            W += pulse_wave(xbank, m, inj[0], inj[1], sh)
    else:
        W = pulse_wave(bank, m, pool.A_b[idx], pool.t_b[idx])
        if inj is not None:
            W += pulse_wave(bank, m, inj[0], inj[1])
    W = pool.ped[idx][:, None] + W
    if arm == "tpl_noise":
        W = np.rint(W + pool.noise(rng, pool.sigp[idx], len(idx)))
    return W


def fit_waves(pw, bank, W, idx, pool):
    y, _ = subtract_ped(bank, W)
    return fit_batch(pw, y, pool.midx[idx], pool.sigp[idx])


# ---------------------------------------------------------------- metrics --
def pair_match(t_list, a_list, T1, T2, tol):
    """Best ordered pair of pulses (a before b) with |t_a - T1| <= tol and |t_b - T2| <= tol.

    t_list / a_list: per-event arrays of pulse onsets / amplitudes.  Returns
    (matched, tA, tB, aA, aB).
    """
    n = len(t_list)
    ok = np.zeros(n, bool)
    o = np.full((n, 4), np.nan)
    for j in range(n):
        t = np.asarray(t_list[j])
        if t.size < 2:
            continue
        a = np.asarray(a_list[j])
        od = np.argsort(t)
        t, a = t[od], a[od]
        best = None
        for p in range(t.size):
            if abs(t[p] - T1[j]) > tol:
                continue
            for q in range(p + 1, t.size):
                if abs(t[q] - T2[j]) <= tol:
                    e = max(abs(t[p] - T1[j]), abs(t[q] - T2[j]))
                    if best is None or e < best[0]:
                        best = (e, p, q)
        if best is not None:
            ok[j] = True
            o[j] = [t[best[1]], t[best[2]], a[best[1]], a[best[2]]]
    return ok, o


def fit_pulses(R, thr):
    """Per-event pulse lists of the fit at threshold thr."""
    n2 = R[:, 12] > thr
    tl, al = [], []
    for j in range(len(R)):
        if n2[j]:
            tl.append([R[j, 9], R[j, 11]])
            al.append([R[j, 8], R[j, 10]])
        else:
            tl.append([R[j, 3]])
            al.append([R[j, 2]])
    return n2, tl, al


def thr_for(dchi_cal, fake):
    d = np.sort(dchi_cal)
    k = int(math.floor(fake * len(d)))
    return float(d[len(d) - 1 - k]) if k < len(d) else 0.0   # at most k events above


# ------------------------------------------------------------- MF method --
def run_mf(bank, W, midx, chunk=2000):
    """Production MF (benchmark settings) on raw waveforms -> per-event onset/amp lists."""
    rec = mf_pulse.MFReconstructor(bank, nsigma=5.0, min_amp=5.0, kmax=4, r_model=0.3,
                                   ped_nsamples=PED_N)
    tl = [[] for _ in range(len(W))]
    al = [[] for _ in range(len(W))]
    for i in range(0, len(W), chunk):
        r = rec.run(W[i:i + chunk], midx[i:i + chunk])
        for c, t, a in zip(r.p_ch, r.p_t0, r.p_amp):
            tl[i + int(c)].append(float(t))
            al[i + int(c)].append(float(a))
    return tl, al


# ---------------------------------------------------------------- cells --
COLS = ["run", "arm", "method", "jit_ns", "dt_ns", "ratio", "n", "thr03", "thr05",
        "eff4_f03", "eff4_f03_lo", "eff4_f03_hi", "eff2_f03", "eff4_f05", "eff4_f05_lo",
        "eff4_f05_hi", "eff2_f05", "fake_f03", "fake_f05", "n_single", "thr_src",
        "eff4_f03_basefirst", "eff4_f03_injfirst", "a1_bias", "a1_rms", "a2_bias", "a2_rms",
        "ddt_bias", "ddt_rms"]


def cell_row(run, arm, method, jit, dt, ratio, ev, dets, singles_n2, thr, thr_src):
    """One TSV row.  dets[f] = (n2 flags, pulse lists t, a) of the doublets at fake level f;
    singles_n2[f] = N=2 flags of the cell's matching singles."""
    n = len(ev["idx"])
    row = dict(run=run, arm=arm, method=method, jit_ns=jit, dt_ns=dt, ratio=ratio, n=n,
               thr03=thr[0], thr05=thr[1], n_single=len(singles_n2[0]), thr_src=thr_src)
    for tag, f in (("f03", 0), ("f05", 1)):
        n2, tl, al = dets[f]
        ok4, o4 = pair_match(tl, al, ev["t1"], ev["t2"], 4.0)
        ok2, _ = pair_match(tl, al, ev["t1"], ev["t2"], 2.0)
        lo, hi = wilson(int(ok4.sum()), n)
        row[f"eff4_{tag}"] = ok4.mean()
        row[f"eff4_{tag}_lo"], row[f"eff4_{tag}_hi"] = lo, hi
        row[f"eff2_{tag}"] = ok2.mean()
        row[f"fake_{tag}"] = float(np.mean(singles_n2[f]))
        if tag == "f03":
            fi = ev["inj_first"]
            row["eff4_f03_basefirst"] = ok4[~fi].mean() if (~fi).any() else np.nan
            row["eff4_f03_injfirst"] = ok4[fi].mean() if fi.any() else np.nan
            if ok4.any():
                r1 = o4[ok4, 2] / ev["A1"][ok4]
                r2 = o4[ok4, 3] / ev["A2"][ok4]
                dd = (o4[ok4, 1] - o4[ok4, 0]) - ev["dt_true"][ok4]
                row.update(a1_bias=r1.mean() - 1, a1_rms=r1.std(), a2_bias=r2.mean() - 1,
                           a2_rms=r2.std(), ddt_bias=dd.mean(), ddt_rms=dd.std())
            else:
                row.update(a1_bias=np.nan, a1_rms=np.nan, a2_bias=np.nan, a2_rms=np.nan,
                           ddt_bias=np.nan, ddt_rms=np.nan)
    return row


def prod_check(rt, slug, pool, bank, pw, thr):
    """Fit and MF on the production synth (50 events/cell; injected pulse later; dt label =
    injected peak - base CFD time, not the onset separation).  Rows: arm 'prodsynth',
    dt_ns = production label, jit_ns = -1; truth = truth_onset_sampled_ns."""
    sd = HERE / "output/pileup_eff_20261024_fix" / slug / "synth"
    man = json.loads((sd / "manifest.json").read_text())
    ix = {n: i for i, n in enumerate(bank.names)}
    rows = []
    for f in man["files"]:
        d = np.load(sd / f["path"], allow_pickle=False)
        names = [str(x) for x in d["base_channel_names"]]
        if any(n not in ix for n in names):
            continue
        m = np.array([ix[n] for n in names], dtype=np.int64)
        wf = d["waveforms"].astype(np.float64)
        sig = np.array([pool.sig_mod.get(n, 0.4) for n in names])
        to = d["truth_onset_sampled_ns"].astype(np.float64)
        n = len(m)
        ev = dict(idx=np.arange(n), inj_first=np.zeros(n, bool), t1=to.min(1), t2=to.max(1),
                  A1=np.full(n, np.nan), A2=np.full(n, np.nan), dt_true=np.abs(to[:, 0] - to[:, 1]),
                  sep_median=float(np.median(np.abs(to[:, 0] - to[:, 1]))))
        y, _ = subtract_ped(bank, wf)
        Rd = fit_batch(pw, y, m, sig)
        dets = [fit_pulses(Rd, thr[0]), fit_pulses(Rd, thr[1])]
        nan = [np.full(n, np.nan)] * 2
        r = cell_row(rt, "prodsynth", "fit", -1.0, float(f["dt_ns"]), float(f["ratio"]), ev, dets, nan,
                     thr[:2], thr[2])
        r["thr_src"] = f"sep_med={ev['sep_median']:.2f}"
        rows.append(r)
        tl, al = run_mf(bank, wf, m)
        r = cell_row(rt, "prodsynth", "mf", -1.0, float(f["dt_ns"]), float(f["ratio"]), ev, [(None, tl, al)] * 2,
                     nan, (np.nan, np.nan), "-")
        r["thr_src"] = f"sep_med={ev['sep_median']:.2f}"
        rows.append(r)
    print(f"[{rt}] prodsynth check: {len(rows) // 2} cells", flush=True)
    return rows




def process_run(rt, args, rng):
    slug, split = RUNS[rt]
    t0 = time.time()
    bank, _ = load_bank(B50, slug)
    xbank, _ = load_bank(WAS, slug)
    global _BANK
    _BANK = bank
    pool = Pool(slug, split, bank, xbank)
    pw = mp.get_context("fork").Pool(args.workers) if args.workers > 1 else None   # forks after _BANK is set
    R, y, ped, keep = real_pool_fit(pool, bank, pw, args.base_t_min)
    slice_pool(pool, R, ped, keep)
    print(f"[{rt}] split {split}: {pool.n_split} stored, {pool.n_clean0} clean, {pool.n} bases "
          f"(t_b >= {args.base_t_min} ns, 100-4000 ADC); cal {int((pool.parity == 0).sum())} / "
          f"eval {int((pool.parity == 1).sum())}; noise acf1..4 {np.round(pool.acf[1:5], 2)}; "
          f"{time.time() - t0:.0f} s", flush=True)
    arms = [a for a in ARMS_ALL if a in args.arms]
    # singles of every arm (all bases), thresholds
    sing, thr = {}, {}
    for arm in arms:
        if arm in ("real", "real2"):
            sing[arm] = pool.Rb
        else:
            W = make_waves(arm, pool, bank, xbank, np.arange(pool.n), rng)
            sing[arm] = fit_waves(pw, bank, W, np.arange(pool.n), pool)
    cal = {a: sing[a][pool.parity == 0, 12] for a in arms}
    for arm in arms:
        src = arm
        if arm == "nf_same":
            src = "tpl_noise" if "tpl_noise" in arms else arm
        thr[arm] = (thr_for(cal[src], 0.03), thr_for(cal[src], 0.05), src)
    for arm in arms:
        ev_ = sing[arm][pool.parity == 1, 12]
        print(f"[{rt}] {arm:9s} thr(0.03)={thr[arm][0]:9.2f} thr(0.05)={thr[arm][1]:9.2f} (from {thr[arm][2]}); "
              f"odd singles fake {np.mean(ev_ > thr[arm][0]):.4f} / {np.mean(ev_ > thr[arm][1]):.4f}; "
              f"chi2_1/ndf median {np.median(sing[arm][:, 4]) / 98:.2f}", flush=True)
    # MF on singles (real arm bases)
    mf_single = None
    if "real" in arms or "real2" in arms:
        mf_single = run_mf(bank, pool.raw, pool.midx)
    cells = [(dt, r, args.jit) for dt in DTS for r in RATIOS]
    cells += [(dt, 1.0, 0.0) for dt in args.fixed_dts]
    rows, store = [], {}
    for ci, (dt, ratio, jit) in enumerate(cells):
        for arm in arms:
            ev = event_table(pool, args.n_eval, dt, ratio, jit, rng, partner=(arm == "real2"))
            idx = ev["idx"]
            inj = (ev["A_inj"], ev["t_inj"])
            W = make_waves(arm, pool, bank, xbank, idx, rng, inj, ev)
            Rd = fit_waves(pw, bank, W, idx, pool)
            dets, s2 = [], []
            for f in (0, 1):
                dets.append(fit_pulses(Rd, thr[arm][f]))
                s2.append(sing[arm][idx, 12] > thr[arm][f])
            rows.append(cell_row(rt, arm, "fit", jit, dt, ratio, ev, dets, s2,
                                 thr[arm][:2], thr[arm][2]))
            if arm in ("real", "real2"):
                tl, al = run_mf(bank, W, pool.midx[idx])
                mfd = [(None, tl, al)] * 2
                ms = [np.array([len(mf_single[0][i]) >= 2 for i in idx])] * 2
                rows.append(cell_row(rt, arm, "mf", jit, dt, ratio, ev, mfd, ms,
                                     (np.nan, np.nan), "-"))
            if arm == "real" and dt == 8 and ratio == 1.0 and jit > 0:
                store["ex8"] = dict(ev=ev, Rd=Rd)
        print(f"[{rt}] cell {ci + 1}/{len(cells)} dt {dt} r {ratio} jit {jit}  {time.time() - t0:.0f} s",
              flush=True)
    if args.prod_check and "real" in arms:
        rows += prod_check(rt, slug, pool, bank, pw, thr["real"])
    if pw is not None:
        pw.close()
    return rows, dict(thr=thr, pool_n=pool.n, n_split=pool.n_split, n_clean=pool.n_clean0,
                      acf=pool.acf[:6].tolist(), wall=time.time() - t0)




ARMS_ALL = ("real", "tpl_noise", "nf_same", "nf_xbank", "real2")


def write_tsv(path, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in rows:
            vals = []
            for c in COLS:
                v = r.get(c, "")
                vals.append(f"{v:.5g}" if isinstance(v, (float, np.floating)) else str(v))
            fh.write("\t".join(vals) + "\n")


def read_tsv(path):
    rows = []
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for ln in fh:
            v = ln.rstrip("\n").split("\t")
            r = {}
            for c, x in zip(hdr, v):
                try:
                    r[c] = float(x)
                except ValueError:
                    r[c] = x
            rows.append(r)
    return rows


def make_report(out):
    """eff_vs_dt.pdf and tables.md from cells.tsv."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    rows = read_tsv(Path(out) / "cells.tsv")
    runs = list(dict.fromkeys(r["run"] for r in rows))
    col = {r: c for r, c in zip(runs, ("#0072B2", "#D55E00", "#009E73"))}

    def sel(**kw):
        return [r for r in rows if all(r[k] == v for k, v in kw.items())]

    def curve(run, arm, method, jit=None, ratio=None):
        rr = [r for r in rows if r["run"] == run and r["arm"] == arm and r["method"] == method
              and (jit is None or r["jit_ns"] == jit) and (ratio is None or r["ratio"] == ratio)]
        rr.sort(key=lambda r: r["dt_ns"])
        return rr

    def page(pdf, arms, title, jit=2.0, f="f03"):
        fig, axs = plt.subplots(1, len(RATIOS), figsize=(13, 3.6), sharey=True)
        for ax, rt in zip(axs, RATIOS):
            for run in runs:
                for arm, ls, mk in arms:
                    for method in (("fit", "mf") if arm in ("real", "real2") else ("fit",)):
                        rr = curve(run, arm, method, jit, rt)
                        if not rr:
                            continue
                        x = [r["dt_ns"] for r in rr]
                        y = np.array([r[f"eff4_{f}"] for r in rr])
                        lo = np.array([r[f"eff4_{f}_lo"] for r in rr])
                        hi = np.array([r[f"eff4_{f}_hi"] for r in rr])
                        st = ls if method == "fit" else "--"
                        ax.errorbar(x, y, [y - lo, hi - y], color=col[run], ls=st, marker=mk if method == "fit" else "x",
                                    ms=4, lw=1.2, capsize=1.5,
                                    label=f"{run} {arm} {method}" if rt == RATIOS[0] else None)
            ax.set_title(f"ratio {rt:g}")
            ax.set_xlabel("dt (ns, nominal)")
            ax.axhline(0.8, color="0.7", lw=0.7)
            ax.set_xlim(2, 34)
            ax.set_ylim(-0.02, 1.03)
        axs[0].set_ylabel("efficiency (both onsets +-4 ns)")
        axs[0].legend(fontsize=5, ncol=1)
        fig.suptitle(title, fontsize=9)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

    with PdfPages(Path(out) / "eff_vs_dt.pdf") as pdf:
        page(pdf, [("real", "-", "o")], "real arm: two-pulse chi2 fit (solid) vs production MF (dashed); fake <= 0.03/event, "
             "dt_true = dt + U(-2,2) ns; error bars Wilson 68 %")
        for arm in ("tpl_noise", "nf_same", "nf_xbank", "real2"):
            if sel(arm=arm):
                page(pdf, [(arm, "-", "o")], f"arm {arm}: fit" + (" (+ MF)" if arm == "real2" else "") + ", fake <= 0.03")
        page(pdf, [("real", "-", "o")], "real arm, fake <= 0.05", f="f05")
        # fixed-dt page (jit 0, ratio 1)
        fig, ax = plt.subplots(1, 2, figsize=(9, 3.6), sharey=True)
        for a_, f in zip(ax, ("f03", "f05")):
            for run in runs:
                for method, ls in (("fit", "-"), ("mf", "--")):
                    rr = curve(run, "real", method, 0.0, 1.0)
                    if rr:
                        a_.plot([r["dt_ns"] for r in rr], [r[f"eff4_{f}"] for r in rr], ls, color=col[run],
                                marker="o", ms=4, label=f"{run} {method}")
            a_.set_title(f"exact dt (jit 0), ratio 1, fake <= {f[1:]}")
            a_.set_xlabel("dt (ns)")
        ax[0].set_ylabel("efficiency +-4 ns")
        ax[0].legend(fontsize=6)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

    L = []
    pc = lambda v: f"{100 * v:.0f}"
    L.append("## real arm, jit 2, fake<=0.03: fit % / MF % (min-max over the 3 run types)\n")
    L.append("| dt | " + " | ".join(f"r {r:g}" for r in RATIOS) + " |")
    L.append("|---|" + "---|" * len(RATIOS))
    for dt in DTS:
        cells = []
        for rt in RATIOS:
            fv = [sel(run=run, arm="real", method="fit", jit_ns=2.0, dt_ns=float(dt), ratio=rt)[0]["eff4_f03"] for run in runs]
            mv = [sel(run=run, arm="real", method="mf", jit_ns=2.0, dt_ns=float(dt), ratio=rt)[0]["eff4_f03"] for run in runs]
            cells.append(f"{pc(min(fv))}-{pc(max(fv))} / {pc(min(mv))}-{pc(max(mv))}")
        L.append(f"| {dt} | " + " | ".join(cells) + " |")
    L.append("\n## dt 8 ns, ratio 1, jit 2 (fake<=0.03 | <=0.05), eff +-4 / +-2 ns, per arm and run\n")
    L.append("| run | arm | method | eff4 f03 [68% CI] | eff2 f03 | eff4 f05 | fake f03 | fake f05 | thr03 |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    for run in runs:
        for arm in ARMS_ALL:
            for method in ("fit", "mf"):
                rr = sel(run=run, arm=arm, method=method, jit_ns=2.0, dt_ns=8.0, ratio=1.0)
                if rr:
                    r = rr[0]
                    L.append(f"| {run} | {arm} | {method} | {pc(r['eff4_f03'])} [{pc(r['eff4_f03_lo'])}-{pc(r['eff4_f03_hi'])}] | "
                             f"{pc(r['eff2_f03'])} | {pc(r['eff4_f05'])} | {r['fake_f03']:.3f} | {r['fake_f05']:.3f} | {r['thr03']:.3g} |")
    L.append("\n## exact dt (jit 0), ratio 1, real arm: fit / MF eff4 at fake<=0.03\n")
    L.append("| dt | " + " | ".join(runs) + " |")
    L.append("|---|" + "---|" * len(runs))
    for dt in sorted({r["dt_ns"] for r in rows if r["jit_ns"] == 0.0}):
        cs = []
        for run in runs:
            f_ = sel(run=run, arm="real", method="fit", jit_ns=0.0, dt_ns=dt, ratio=1.0)
            m_ = sel(run=run, arm="real", method="mf", jit_ns=0.0, dt_ns=dt, ratio=1.0)
            cs.append(f"{pc(f_[0]['eff4_f03'])} / {pc(m_[0]['eff4_f03'])}" if f_ and m_ else "-")
        L.append(f"| {dt:g} | " + " | ".join(cs) + " |")
    L.append("\n## production synth (dt label = injected peak - base CFD time; median onset separation in col 2), "
             "eff4 at fake<=0.03: fit / MF, ratio 0.5 | 1 | 2, 50 events per cell\n")
    L.append("| label dt | onset sep | " + " | ".join(runs) + " |")
    L.append("|---|---|" + "---|" * len(runs))
    for dt in (8.0, 12.0, 16.0, 20.0, 24.0, 32.0):
        cs, sep = [], ""
        for run in runs:
            parts = []
            for rt in (0.5, 1.0, 2.0):
                f_ = sel(run=run, arm="prodsynth", method="fit", dt_ns=dt, ratio=rt)
                m_ = sel(run=run, arm="prodsynth", method="mf", dt_ns=dt, ratio=rt)
                if f_ and m_:
                    parts.append(f"r{rt:g} {pc(f_[0]['eff4_f03'])}/{pc(m_[0]['eff4_f03'])}")
                    sep = str(f_[0]["thr_src"]).replace("sep_med=", "")
            cs.append(", ".join(parts) if parts else "-")
        L.append(f"| {dt:g} | {sep} | " + " | ".join(cs) + " |")
    (Path(out) / "tables.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


def tune_eps(args, out):
    """eps scan (kap 0) on the calibration half only (even base index) of 3.49 GeV H2:
    threshold from the even singles, efficiency on even-base doublets, ratio 1."""
    rt = "3.49 GeV H2"
    slug, split = RUNS[rt]
    bank, _ = load_bank(B50, slug)
    global _BANK
    _BANK = bank
    pool0 = Pool(slug, split, bank, None)
    pw = mp.get_context("fork").Pool(args.workers)
    lines = ["eps\tthr03\tthr05\tdt_nom\teff4_f03\teff4_f05"]
    for eps in (0.0, 0.005, 0.01, 0.02, 0.03, 0.05, 0.08):
        OPT.update(eps=eps, kap=0.0)
        R, y, ped, keep = real_pool_fit(pool0, bank, pw, args.base_t_min)
        pl = Pool.__new__(Pool)
        pl.__dict__.update(pool0.__dict__)
        slice_pool(pl, R, ped, keep)
        rng = np.random.default_rng(1)
        cal = pl.Rb[pl.parity == 0, 12]
        t3, t5 = thr_for(cal, 0.03), thr_for(cal, 0.05)
        for dt in (2, 4, 6, 8):
            ev = event_table(pl, 300, dt, 1.0, args.jit, rng, half=0)
            W = make_waves("real", pl, bank, None, ev["idx"], rng, (ev["A_inj"], ev["t_inj"]))
            Rd = fit_waves(pw, bank, W, ev["idx"], pl)
            e = []
            for t in (t3, t5):
                n2, tl, al = fit_pulses(Rd, t)
                e.append(pair_match(tl, al, ev["t1"], ev["t2"], TOL)[0].mean())
            lines.append(f"{eps}\t{t3:.1f}\t{t5:.1f}\t{dt}\t{e[0]:.3f}\t{e[1]:.3f}")
        print(lines[-4:], flush=True)
    pw.close()
    (out / "tuning_eps.tsv").write_text("\n".join(lines) + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="output/ml_pileup_20261027/U0")
    ap.add_argument("--runs", default="all", help="comma list of run types or 'all' (%s)" % ", ".join(RUNS))
    ap.add_argument("--arms", default=",".join(ARMS_ALL))
    ap.add_argument("--n-eval", type=int, default=1000, help="doublets per cell (evaluation half)")
    ap.add_argument("--jit", type=float, default=2.0, help="dt_true = dt + U(-jit, jit) ns")
    ap.add_argument("--fixed-dts", default="6,8,10,16,20,24", help="extra ratio-1 cells with jit 0 (16,20,24: production-synth phase convention)")
    ap.add_argument("--eps", type=float, default=0.02, help="fractional template-error sigma term")
    ap.add_argument("--kap", type=float, default=0.0)
    ap.add_argument("--search", type=float, default=80.0, help="onset search half-range (ns)")
    ap.add_argument("--min-sep", type=float, default=1.0, help="minimum t2 - t1 (ns)")
    ap.add_argument("--base-t-min", type=float, default=150.0)
    ap.add_argument("--seed", type=int, default=20261027)
    ap.add_argument("--workers", type=int, default=18)
    ap.add_argument("--report-only", action="store_true", help="redo eff_vs_dt.pdf and tables.md from cells.tsv")
    ap.add_argument("--no-prod-check", dest="prod_check", action="store_false",
                    help="skip the production-synth cross-check (arm prodsynth)")
    ap.add_argument("--tune-eps", action="store_true", help="eps scan on the calibration half (3.49 GeV H2), then exit")
    ap.add_argument("--quick", action="store_true", help="few cells, 100 events (smoke test)")
    a = ap.parse_args(argv)
    a.arms = a.arms.split(",")
    a.fixed_dts = [int(x) for x in a.fixed_dts.split(",") if x]
    OPT.update(eps=a.eps, kap=a.kap, search=a.search, min_sep=a.min_sep)
    runs = list(RUNS) if a.runs == "all" else a.runs.split(",")
    global DTS, RATIOS
    if a.quick:
        DTS, RATIOS, a.n_eval, a.fixed_dts = (4, 8, 16), (0.5, 1.0), 100, []
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    if a.report_only:
        return make_report(out)
    if a.tune_eps:
        return tune_eps(a, out)
    t_all = time.time()
    rows, meta = [], {}
    for i, rt in enumerate(runs):
        rng = np.random.default_rng([a.seed, i])
        r, m = process_run(rt, a, rng)
        rows += r
        meta[rt] = m
    write_tsv(out / "cells.tsv", rows)
    json.dump(dict(args=vars(a), opt=OPT, meta=meta, wall_s=time.time() - t_all),
              open(out / "summary.json", "w"), indent=1, default=str)
    make_report(out)
    print(f"[done] {len(rows)} rows -> {out / 'cells.tsv'}  wall {time.time() - t_all:.0f} s")


if __name__ == "__main__":
    main()
