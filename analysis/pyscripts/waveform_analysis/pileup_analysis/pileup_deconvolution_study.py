#!/usr/bin/env python3
"""pileup_deconvolution_study.py -- matched-filter deconvolution efficiency on stage-1 pile-up waveforms.

    python3 waveform_analysis/pileup_analysis/pileup_deconvolution_study.py <stage1.npz>

Standalone (numpy, scipy, matplotlib).  Templates: the generator's own bank (meta "bank",
relative to analysis/pyscripts, sha1-checked), evaluated as pileup_generator_stage1.py does: per-module cubic B-spline
(spline_t / spline_c), tabulated on a 0.05 ns grid (value and slope, float32), linear
interpolation, s = 0 outside (0, T_end), t = ns after the onset t0.

Reconstruction (noise model "mad"): pedestal = clipped mean of the first PED_NSAMPLES samples, sigma = 1.4826 MAD;
a waveform is fitted if its 3-sample moving average exceeds max(3 sigma, 3 ADC).  Samples at 4095 get weight 0; a
clipped waveform with < 6 unclipped samples in a pulse's fit window or < 2 on its leading edge loses its pulses.

Algorithm (class MF; matched filter over onsets t0 = 4 j - phi, phi = 0..3 ns, 24-sample kernels K, A = C / |K|^2):
  A  find_candidates: the matched-filter amplitude vs onset (1 ns grid) is smoothed by a Gaussian of PEAK_SMOOTH_NS;
     its local maxima with prominence >= max(NSIGMA sigma_A, MIN_AMP), height >= PEAK_FRAC x the highest and
     >= MIN_SEP_NS apart are the candidates (the KMAX highest), each moved to the best raw score within +-2 ns
     (TSpectrum's job without ROOT: Gaussian smoothing + prominence instead of Markov smoothing).  For the candidates
     after the highest one, sigma_A^2 grows by (A_h shape_sigma_m(t - t_h) sqrt(a_band_m / A_h))^2, the expected residual scatter of the
     highest pulse (amplitude A_h, onset t_h, both from the highest smoothed peak) at the candidate's position
     (0 outside the template), so a few-% template deviation on a large pulse's tail seeds no second candidate.
     Onsets run to the end of the window (kernel overlap-energy normalisation); a pulse near the end needs LATE_MIN_TAIL
     samples after its peak inside the window.
  B  fit_all: all candidates at once: linear amplitudes (non-positive ones dropped) and a constant baseline offset (clipped to
     +-BASE_CLIP sigma, then the amplitudes are solved again with it fixed), 4 joint Gauss-Newton steps on all (A_k, t0_k), pulses
     closer than MERGE_NS merged (the smaller dropped) and refitted; then every pulse must have A >= NSIGMA sigma_A (sigma_A from
     the sigma(t, A) variance) and raw max - baseline > RAW_NSIGMA sigma within its own window, the weakest failing pulse is
     dropped and the fit redone.  The reported pedestal is the 30-sample baseline plus the offset.
  C  add_shoulders: the best matched-filter onset of the residual, >= MIN_SEP_NS from every pulse, with
     A >= max(NSIGMA sigma_A, MIN_AMP) and >= R_MODEL x every pulse within 96 ns, is added and all pulses refitted
     (B); it stays only if chi2_old - chi2_new > BIC_PENALTY 2 ln N_w (BIC; N_w = unclipped samples in the fit windows).
     Repeated until a candidate fails or KMAX pulses.
  sigma(t): per sample the variance is sigma^2 + (sum_k A_k shape_sigma_m(t - t0_k) sqrt(a_band_m / max(A_k, A_MIN_SCALE)))^2,
     shape_sigma = the bank's pulse-to-pulse scatter around module m's template (fraction of the peak, >= SHAPE_SIGMA_MIN)
     measured at the amplitude a_band_m, scaled so that the variance is proportional to A; used in C's BIC test, in A's
     prominence term above and in the chi2/dof below, not in B.

Per event: main = fitted pulse nearest t1; secondary = the other pulse nearest t2 (either order: dt = t2 - t1 may be negative).  Found
(timing) if |t0 - t2| <= T_TOL_NS, found (timing + amplitude) if also |A / a2 - 1| <= A_TOL.
Efficiency per (module, A1, r, dt) cell = found / n +- sqrt(e (1 - e) / n) (binomial).
r = 0 controls: false-positive rate = fraction of events with >= 2 pulses.  chi2/dof =
sum_i w_i (y_i - ped - model_i)^2 / var_i over the union of the pulses'
32-sample fit windows (w = 0 on clipped samples), dof = n_used - 2 K - 1 (the baseline offset is a fit parameter); var_i = the sigma(t) variance with
sigma^2 = noise_sigma^2 + 1/12 (generator noise + rounding).

Writes efficiency.tsv, efficiency_tables.pdf (per module r x dt table pooled over
A1, then split by A1, then the pull page) and pileup_examples.pdf (EX_MODULE, EX_A1: one page per r with the
first FOUND event of each dt cell (dt negative and positive), then the r = 0 control);
all written to this script's directory; NPROC worker processes.  efficiency.tsv also has the main-pulse amplitude
pull (A1_fit - A1) / A1 per cell (median and 1.4826 MAD; found events, for r = 0 the one-pulse events), plotted
vs dt in the last page of efficiency_tables.pdf, and ped_bias = median (pedestal_fit - pedestal_true) / noise_sigma
per cell (all events; an early secondary, dt <= -32, lies in the pedestal window).
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")                     # one BLAS thread per worker

import hashlib
import json
import multiprocessing as mp
import sys

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks

OUT_DIR = os.path.dirname(os.path.abspath(__file__))   # outputs go next to this script
PYDIR = os.path.dirname(os.path.dirname(OUT_DIR))       # analysis/pyscripts (the bank path in meta is relative to it)
NPROC = 8             # worker processes
T_TOL_NS = 4.0        # secondary timing window |t0 - t2| (1 sample)
T_TOL2_NS = 2.0       # tighter timing window, recorded in the TSV
A_TOL = 0.30          # secondary amplitude window |A / a2 - 1|
PED_NSAMPLES = 30     # baseline samples (0-116 ns; the main onset is at >= 152 ns)
# matched filter (same values in pileup_rates.py)
KMAX = 6              # at most this many pulses per waveform
NSIGMA, MIN_AMP = 5.0, 5.0   # a pulse needs A >= max(NSIGMA sigma_A, MIN_AMP) (stage A: prominence >= that)
MIN_SEP_NS = 8.0      # smallest onset separation of two pulses (= smallest generated dt)
MERGE_NS = 6.5        # post-fit merge below this; 8 ns - 3 sigma_dt, sigma_dt <= 0.5 ns measured on pseudodata
R_MODEL = 0.02        # stage C: a shoulder needs A >= R_MODEL x every pulse within one kernel span (96 ns)
PEAK_SMOOTH_NS = 8.0  # stage A: Gaussian sigma (ns) of the amplitude smoothing before the peak search (2 samples)
PEAK_FRAC = 0.02      # stage A: a candidate's smoothed height >= PEAK_FRAC x the waveform's highest candidate
SHAPE_SIGMA_MIN = 0.005   # floor of shape_sigma (fraction of the peak) in the sigma(t) variance
A_MIN_SCALE = 20.0        # ADC; floor so the sqrt(a_band / A) scaling of shape_sigma does not blow up for tiny candidates
BASE_CLIP = 1.0           # free-baseline offset limited to ±1 sigma: the 30-sample pedestal mean has a 0.2 sigma error, larger offsets mean a pulse in the pedestal window; at 3 sigma the offset absorbed smeared-template residuals of large pulses and created 4-5 ADC false pulses ~105 ns before the main one (control false positives 0.42 % -> 0.08 % at 1 sigma, pseudodata)
RAW_NSIGMA = 5.0          # an accepted pulse needs raw max - baseline > RAW_NSIGMA sigma_ped within its own window
LATE_MIN_TAIL = 6         # a pulse near the window end needs this many samples after its peak inside the window
BIC_PENALTY = 1.0         # stage C: chi2 improvement must exceed BIC_PENALTY x 2 ln N_w; 1, 1.5, 2, 3 give the same pseudodata efficiencies (+8/0.1 72-73 %, +12/0.1 97.3-97.4 %, control false positives max 0.42 %); the criterion (control false positives <= 0.1 %) is met by none, so the smallest is kept
DT, NPH, L, LW, LEAD = 4.0, 4, 24, 32, 3              # ns/sample, phases, kernel / fit window, lead
LUT_DT, SIGMA_MIN, CAND_NSIGMA = 0.05, 0.3, 3.0
CHUNK = 50000         # events per worker task
# pileup_examples.pdf: every (r, dt) cell of EX_MODULE, EX_A1 (one page per r), plus the r = 0 control
EX_MODULE, EX_A1 = "W493", 1000.0
D = {}                # event arrays, set in the parent before the fork
F32 = np.float32


class Bank:
    """Generator templates: B-spline rows on the 0.05 ns grid, linear interpolation.  slut: module m's shape_sigma
    on the same grid (linear in the bank's t_grid, ends held), at least SHAPE_SIGMA_MIN; aband: its a_band."""

    def __init__(self, meta, names):
        from scipy.interpolate import BSpline
        path = os.path.join(PYDIR, meta["bank"])
        with open(path, "rb") as f:
            if hashlib.sha1(f.read()).hexdigest() != meta["bank_sha1"]:
                raise SystemExit(f"[ERROR] {path} changed since the waveforms were generated")
        z = np.load(path, allow_pickle=False)
        k, self.T_end = int(z["spline_k"]), float(z["T_end"])
        self.n = int(round(self.T_end / LUT_DT)) + 1
        tg = LUT_DT * np.arange(self.n)
        row = {str(s): i for i, s in enumerate(z["names"])}
        tb = np.asarray(z["t_grid"], dtype=np.float64)
        if "a_band" not in z.files:
            raise SystemExit(f"[ERROR] {path} has no a_band: rebuild the bank with the current template_extract.py")
        lut, dlut, slut, aband = [], [], [], []
        for nm in names:
            r = row.get(nm)
            if r is None or z["valid"][r] <= 0 or z["fallback"][r] != 0:
                raise SystemExit(f"[ERROR] {nm} has no own template in {path}")
            t = z["spline_t"][r, :z["spline_nt"][r]]
            s = BSpline(t, z["spline_c"][r, :z["spline_nc"][r]], k)
            x = np.minimum(tg, t[-k - 1])                # hold the end value beyond the base interval
            lut.append(s(x))
            dlut.append(s.derivative()(x))
            ss = z["shape_sigma"][r]
            ok = np.isfinite(ss)
            slut.append(np.fmax(np.interp(tg, tb[ok], ss[ok]), SHAPE_SIGMA_MIN))   # floor
            aband.append(z["a_band"][r])
        self.lut, self.dlut, self.slut = (np.asarray(v, dtype=np.float64).astype(F32) for v in (lut, dlut, slut))
        self.aband = np.asarray(aband, dtype=F32)          # amplitude (ADC) at which shape_sigma was measured
        self.t_pk = np.argmax(self.lut, axis=1) * LUT_DT
        self.lut, self.dlut, self.slut = self.lut.reshape(-1), self.dlut.reshape(-1), self.slut.reshape(-1)
        mi = np.arange(len(names))
        K = self.shape(mi[:, None, None], DT * np.arange(L)[None, None, :]
                       + np.arange(NPH)[None, :, None], np.float64)   # K[m, phi, l] = s(4 l + phi)
        self.K, self.K2, self.Kn2 = K.astype(F32), (K * K).astype(F32), (K * K).sum(2).astype(F32)
        self.Kinv = (1.0 / np.sqrt(np.maximum(self.Kn2, 1e-12))).astype(F32)

    def shape(self, m, t, dtype=F32, slope=False, sigma=False):
        """Unit-peak shape (and slope; sigma=True: shape_sigma) at module m, t ns after onset; 0 outside (0, T_end)."""
        t = np.asarray(t, dtype=dtype)
        x = t * dtype(1.0 / LUT_DT)
        i = np.clip(np.floor(x), 0, self.n - 2)
        f = x - i
        b = np.asarray(m, dtype=np.int64) * self.n + i.astype(np.int64)
        ins = (t > 0) & (t < self.T_end)
        out = [np.where(ins, tab[b].astype(dtype) * (1 - f) + tab[b + 1].astype(dtype) * f, dtype(0))
               for tab in ((self.slut,) if sigma else (self.lut, self.dlut) if slope else (self.lut,))]
        return out if slope else out[0]


class MF:
    """Matched filter + joint template fit, in three stages (see pileup_deconvolution_study.py, Algorithm):
    A find_candidates, B fit_all, C add_shoulders; find_pulses runs them in order."""

    def __init__(self, bank):
        self.b = bank
        # K2c[m, phi, l]: sum of K^2 over the first l kernel samples (the overlap energy of a kernel cut by the window end)
        self.K2c = np.concatenate([np.zeros_like(bank.K2[:, :, :1]), np.cumsum(bank.K2, axis=2)], axis=2)

    def basis(self, m, t0, ns):
        return self.b.shape(m[:, None], DT * np.arange(ns, dtype=F32)[None, :] - t0[:, None]).astype(F32)

    @staticmethod
    def wmask(t0, ns):
        """Fit window of a pulse with onset t0: LW samples starting LEAD samples before the onset."""
        w0 = (np.ceil(t0 / DT).astype(np.int64) - LEAD)[:, None]
        ii = np.arange(ns)[None, :]
        return (ii >= w0) & (ii < w0 + LW)

    @staticmethod
    def windows(t0, valid, ns):
        """Union of the fit windows of the valid pulses [n, ns]."""
        out = np.zeros((t0.shape[0], ns), dtype=bool)
        for k in range(KMAX):
            out |= valid[:, k, None] & MF.wmask(t0[:, k], ns)
        return out

    def variance(self, m, t0, a, valid, sig, ns):
        """sigma(t) variance [n, ns]: sig^2 + (sum_k a_k shape_sigma_m(t - t0_k) sqrt(a_band_m / max(a_k, A_MIN_SCALE)))^2
        over the valid pulses (shape_sigma is measured at a_band_m; the variance grows in proportion to a_k)."""
        xt = DT * np.arange(ns, dtype=F32)[None, :]
        ab = self.b.aband[m][:, None]
        ss = sum(np.where(valid[:, k, None], a[:, k, None] * self.b.shape(m[:, None], xt - t0[:, k, None], sigma=True)
                          * np.sqrt(ab / np.maximum(a[:, k, None], A_MIN_SCALE)), 0.0) for k in range(t0.shape[1]))
        return sig[:, None] ** 2 + ss ** 2

    def corr(self, r, Kb):
        return np.matmul(sliding_window_view(r, L, axis=1)[:, :r.shape[1] - L + 1], Kb.transpose(0, 2, 1))

    def score(self, r, w, m):
        """Matched filter of rows r on the 1 ns onset grid (index q = t0 + 3 ns), three arrays [n, 4 ns]: score C / sqrt(D),
        amplitude C / D and sigma_A / sigma = 1 / sqrt(D), with C = sum r K and D = sum w K^2 (= |K|^2 if w is None).  Onsets
        up to the end of the window are tried (kernel samples past it are missing: D is the overlap energy, as for clipped
        samples).  Onsets whose kernel keeps < 5 % of |K|^2, or whose peak has < LATE_MIN_TAIL samples after it inside the
        window, get score -inf and amplitude 0."""
        b, n, ns = self.b, r.shape[0], r.shape[1]
        C = self.corr(np.pad(r, ((0, 0), (0, L - 1))), b.K[m])      # [n, j, phi]: onset t0 = 4 j - phi, j = 0 .. ns - 1
        if w is None:
            D = self.K2c[m][:, :, np.minimum(L, ns - np.arange(ns))].transpose(0, 2, 1)
        else:
            D = self.corr(np.pad(w, ((0, 0), (0, L - 1))), b.K2[m])
        good = D >= 0.05 * b.Kn2[m][:, None, :]
        D = np.maximum(D, 1e-12)
        on_grid = lambda x: np.ascontiguousarray(x[:, :, ::-1]).reshape(n, -1)    # noqa: E731  q = 4 j + 3 - phi
        Sg, Ag, Ig = on_grid(np.where(good, C / np.sqrt(D), -np.inf)), on_grid(np.where(good, C / D, 0.0)), on_grid(1.0 / np.sqrt(D))
        late = (np.arange(4 * ns, dtype=F32) - 3.0)[None, :] + b.t_pk[m][:, None] > DT * (ns - 1 - LATE_MIN_TAIL)
        return np.where(late, -np.inf, Sg), np.where(late, 0.0, Ag), Ig

    @staticmethod
    def onset(Sg, row, q):
        """Onset t0 (ns) of the score maximum at grid index q of row `row`: vertex of the parabola through q - 1, q, q + 1."""
        nq = Sg.shape[1]
        s0, sm, sp = Sg[row, q], Sg[row, np.clip(q - 1, 0, nq - 1)], Sg[row, np.clip(q + 1, 0, nq - 1)]
        with np.errstate(invalid="ignore"):
            den = sm - 2.0 * s0 + sp
            okp = (q > 0) & (q < nq - 1) & (den < 0) & np.isfinite(den)
            delta = np.where(okp, 0.5 * (sm - sp) / np.where(okp, den, 1.0), 0.0)
        return (q - 3).astype(F32) + np.clip(delta, -1.0, 1.0).astype(F32)

    def find_candidates(self, y, w, m, sig):
        """Stage A, one pass over the waveform: local maxima of the matched-filter amplitude smoothed by a Gaussian of
        PEAK_SMOOTH_NS, with prominence >= max(NSIGMA sigma_A, MIN_AMP) (for the candidates after the highest one sigma_A
        grows by the highest pulse's shape scatter at their lag, see the module docstring), height >= PEAK_FRAC x the
        row's highest one and >= MIN_SEP_NS apart; the KMAX highest, each moved to the best raw score within +-2 ns
        and kept if its raw amplitude there is >= max(NSIGMA sigma_A, MIN_AMP).  Returns (t0, valid)."""
        n = y.shape[0]
        Sg, Ag, Ig = self.score(y, w, m)
        nq = Sg.shape[1]
        Asm = gaussian_filter1d(Ag, PEAK_SMOOTH_NS, axis=1, mode="nearest")      # grid step 1 ns
        thr = np.maximum(NSIGMA * sig[:, None] * Ig, MIN_AMP)                   # max(NSIGMA sigma_A, MIN_AMP) per onset
        # a local maximum can pass the prominence cut only if it lies >= thr above the waveform's lowest value.  Most
        # waveforms have one such maximum p, which is also their highest inner point; its prominence is its height minus
        # the higher of the lowest values left and right of it.  scipy's find_peaks (slow, one waveform per call) for
        # the others.
        inner = Asm[:, 1:-1]
        is_max = (inner > Asm[:, :-2]) & (inner >= Asm[:, 2:]) & (inner >= Asm.min(1)[:, None] + thr[:, 1:-1])
        n_max = is_max.sum(1)
        simple = (n_max == 1) & (np.argmax(is_max, axis=1) == np.argmax(inner, axis=1))
        one = np.nonzero(simple)[0]
        p = np.argmax(is_max[one], axis=1) + 1
        lowest_left = np.minimum.accumulate(Asm[one], axis=1)[np.arange(one.size), p]
        lowest_right = np.minimum.accumulate(Asm[one, ::-1], axis=1)[np.arange(one.size), Asm.shape[1] - 1 - p]
        prom = Asm[one, p] - np.maximum(lowest_left, lowest_right)
        ok = prom >= thr[one, p]
        row, q, rank = [one[ok]], [p[ok]], [np.zeros(ok.sum(), dtype=np.int64)]
        for i in np.nonzero((n_max >= 1) & ~simple)[0]:
            qi, pr = find_peaks(Asm[i], distance=int(MIN_SEP_NS), prominence=thr[i])
            pr = pr["prominences"]
            height = Asm[i, qi]
            keep = height >= PEAK_FRAC * height.max(initial=0.0)
            qi, pr = qi[keep], pr[keep]
            o = np.argsort(-Asm[i, qi])                            # highest first
            qi, pr = qi[o], pr[o]
            if qi.size > 1:                                        # later candidates: sigma_A plus the highest pulse's scatter
                sh = self.b.shape(m[i], (qi[1:] - qi[0]).astype(F32), sigma=True) * Ag[i, qi[0]] \
                    * np.sqrt(self.b.aband[m[i]] / max(Ag[i, qi[0]], A_MIN_SCALE))
                good = pr[1:] >= np.maximum(NSIGMA * np.hypot(sig[i] * Ig[i, qi[1:]], sh), MIN_AMP)
                qi = np.concatenate([qi[:1], qi[1:][good]])
            qi = qi[:KMAX]
            row.append(np.full(qi.size, i))
            q.append(qi)
            rank.append(np.arange(qi.size))
        row, q, rank = np.concatenate(row), np.concatenate(q), np.concatenate(rank)
        near = np.clip(q[:, None] + np.arange(-2, 3), 0, nq - 1)   # best raw score within +-2 ns (all 4 phases)
        q = near[np.arange(q.size), np.argmax(Sg[row[:, None], near], axis=1)]
        t0, valid = np.zeros((n, KMAX), F32), np.zeros((n, KMAX), bool)
        t0[row, rank] = self.onset(Sg, row, q)
        valid[row, rank] = Ag[row, q] >= np.maximum(NSIGMA * sig[row] * Ig[row, q], MIN_AMP)
        return t0, valid

    def solve(self, y, w, m, t0, valid, ns, sig):
        """Joint linear fit on the union of fit windows: amplitudes of the valid pulses (non-positive ones dropped) and a
        constant baseline offset db (one more linear parameter, a column of ones over the windows).  If |db| > BASE_CLIP
        sigma, db is clipped to +-BASE_CLIP sigma and the amplitudes are solved again with db fixed.  Returns
        (a, valid, S, wm, db)."""
        c, km = t0.shape
        valid = valid.copy()
        wm = np.zeros((c, ns), dtype=F32)
        S = np.zeros((c, km, ns), dtype=F32)
        for k in range(km):
            rows = np.nonzero(valid[:, k])[0]
            if rows.size:
                S[rows, k] = self.basis(m[rows], t0[rows, k], ns)
                wm[rows] = np.maximum(wm[rows], self.wmask(t0[rows, k], ns).astype(F32))
        if w is not None:
            wm = wm * w
        a, valid, db = self.lin_solve(y, S, wm, valid, np.ones(c, dtype=bool), np.zeros(c, dtype=F32))
        big = np.abs(db) > BASE_CLIP * sig
        if big.any():
            r = np.nonzero(big)[0]
            db[r] = np.clip(db[r], -BASE_CLIP * sig[r], BASE_CLIP * sig[r])
            a[r], valid[r], _ = self.lin_solve(y[r], S[r], wm[r], valid[r], np.zeros(r.size, dtype=bool), db[r])
        return a, valid, S, wm, db

    @staticmethod
    def lin_solve(y, S, wm, valid, free, db):
        """Linear amplitudes (and, for rows with free = True, the baseline offset) by weighted least squares on y - db
        (fixed offset for the other rows); non-positive amplitudes are dropped one at a time.  Returns (a, valid, db)."""
        c, km, ns = S.shape
        valid = valid.copy()
        yf = y - np.where(free, 0.0, db)[:, None]
        eye = np.eye(km + 1, dtype=F32)[None]
        ones = free[:, None, None] * np.ones((c, 1, ns), dtype=F32)
        for _ in range(km):
            Sv = np.concatenate([S * valid[:, :, None], ones], axis=1)
            Sw = Sv * wm[:, None, :]
            G = np.matmul(Sw, Sv.transpose(0, 2, 1))
            ve = np.concatenate([valid, (free & (wm.sum(1) > 0))[:, None]], axis=1)
            G = G + (~ve)[:, :, None] * eye + 1e-6 * eye * G.max(axis=(1, 2), keepdims=True)
            sol = np.linalg.solve(G, np.matmul(Sw, yf[:, :, None]))[:, :, 0].astype(F32)
            anew = np.where(valid, sol[:, :km], 0.0)
            neg = valid & (anew <= 0.0)
            if not neg.any():
                break
            has = neg.any(axis=1)
            valid[np.nonzero(has)[0], np.argmin(np.where(valid, anew, np.inf), axis=1)[has]] = False
        return anew, valid, np.where(free & ve[:, km], sol[:, km], db).astype(F32)

    def fit_all(self, y, w, m, sig, t0, valid, ns):
        """Stage B: joint fit of all pulses of each row (fit_group), rows grouped by their number of pulses K (the cost
        grows with K^2); then every row's weakest pulse failing `insignificant` is dropped and the row refitted, until
        none fails.  Returns (t0, a, valid, model = baseline offset + sum of the fitted pulses [n, ns], db)."""
        t0, valid = t0.copy(), valid.copy()
        a, model, db = np.zeros(t0.shape, F32), np.zeros(y.shape, F32), np.zeros(y.shape[0], F32)
        rows = np.arange(y.shape[0])
        for _ in range(KMAX + 1):
            t0[rows], a[rows], valid[rows], model[rows], db[rows] = self.fit_rows(
                y[rows], None if w is None else w[rows], m[rows], sig[rows], t0[rows], valid[rows], ns)
            bad = self.insignificant(y[rows], None if w is None else w[rows], m[rows], sig[rows], t0[rows], a[rows],
                                     valid[rows], db[rows], ns)
            if not bad.any():
                break
            rows = rows[bad.any(axis=1)]
            valid[rows] &= ~bad[bad.any(axis=1)]
        return t0, a, valid, model, db

    def fit_rows(self, y, w, m, sig, t0, valid, ns):
        """fit_group over the rows grouped by their number of valid pulses (valid pulses first); rows without one: model 0."""
        order = np.argsort(~valid, axis=1, kind="stable")
        t0, valid = np.take_along_axis(t0, order, axis=1), np.take_along_axis(valid, order, axis=1)
        a, model, db = np.zeros(t0.shape, F32), np.zeros(y.shape, F32), np.zeros(y.shape[0], F32)
        for K in range(1, KMAX + 1):
            g = np.nonzero(valid.sum(1) == K)[0]
            if g.size:
                t0[g, :K], a[g, :K], valid[g, :K], model[g], db[g] = self.fit_group(
                    y[g], None if w is None else w[g], m[g], sig[g], t0[g, :K], valid[g, :K], ns)
        t0[~valid], a[~valid] = 0.0, 0.0
        return t0, a, valid, model, db

    def insignificant(self, y, w, m, sig, t0, a, valid, db, ns):
        """[n, KMAX] bool, at most one pulse per row (the weakest): a pulse fails if A < NSIGMA sigma_A (sigma_A = 1 / sqrt(sum over
        its fit window of s^2 / var), var the sigma(t, A) variance) or if its raw maximum above the fitted baseline within its
        own window is <= RAW_NSIGMA sigma (clipped samples count as large)."""
        var = self.variance(m, t0, a, valid, sig, ns)
        xt = DT * np.arange(ns, dtype=F32)[None, :]
        yr = y - db[:, None] if w is None else np.where(w > 0, y - db[:, None], 1e9)
        snr, fail = np.full(t0.shape, np.inf), np.zeros(t0.shape, dtype=bool)
        for k in range(KMAX):
            rows = np.nonzero(valid[:, k])[0]
            if rows.size == 0:
                continue
            win = self.wmask(t0[rows, k], ns) * (1.0 if w is None else w[rows])
            s = self.b.shape(m[rows, None], xt - t0[rows, k, None])
            sa = 1.0 / np.sqrt(np.maximum((win * s * s / var[rows]).sum(1), 1e-12))
            snr[rows, k] = a[rows, k] / sa
            raw = np.where(self.wmask(t0[rows, k], ns), yr[rows], -np.inf).max(axis=1)
            fail[rows, k] = (snr[rows, k] < NSIGMA) | (raw <= RAW_NSIGMA * sig[rows])
        worst = np.argmin(np.where(fail, snr, np.inf), axis=1)
        out = np.zeros(t0.shape, dtype=bool)
        out[np.arange(t0.shape[0]), worst] = fail.any(axis=1)
        return out

    def fit_group(self, y, w, m, sig, t0, valid, ns, nit=4):
        """Joint fit of K pulses per row: linear amplitudes and baseline offset (solve, non-positive amplitudes dropped), nit
        joint Gauss-Newton steps on all (a_k, t0_k) at the current offset, then pulses closer than MERGE_NS merged (the smaller
        dropped) and the Gauss-Newton redone for those rows; final linear solve.  Returns (t0, a, valid, model, db)."""
        c, km = t0.shape
        t0 = t0.astype(F32).copy()
        a, valid, _, wm, db = self.solve(y, w, m, t0, valid, ns, sig)
        xt = DT * np.arange(ns, dtype=F32)[None, :]
        eye = np.eye(2 * km, dtype=F32)[None]
        redo = np.ones(c, dtype=bool)
        for _ in range(2):
            upd = valid & redo[:, None]
            for _ in range(nit):
                Sk, Gk = np.zeros((2, c, km, ns), dtype=F32)
                for k in range(km):
                    if valid[:, k].any():
                        s_, ds_ = self.b.shape(m[:, None], xt - t0[:, k, None], slope=True)
                        Sk[:, k], Gk[:, k] = s_ * valid[:, k, None], -a[:, k, None] * ds_ * valid[:, k, None]
                Jm = np.concatenate([Sk, Gk], axis=1)
                Jw = Jm * wm[:, None, :]
                r = (y - db[:, None] - np.matmul(a[:, None, :], Sk)[:, 0, :]) * wm
                N = np.matmul(Jw, Jm.transpose(0, 2, 1))
                dg = np.diagonal(N, axis1=1, axis2=2)
                inv = np.concatenate([~valid, ~valid], axis=1) | (dg <= 1e-12)
                N = N + eye * (1e-3 * dg[:, :, None] + inv[:, :, None])
                sol = np.linalg.solve(N, np.matmul(Jw, r[:, :, None]))[:, :, 0]
                a = np.where(valid, np.maximum(a + np.where(upd, sol[:, :km], 0.0), 0.2 * a), 0.0).astype(F32)
                t0 = (t0 + np.where(upd, np.clip(sol[:, km:], -2.0, 2.0), 0.0)).astype(F32)
            both = valid[:, :, None] & valid[:, None, :] & ~np.eye(km, dtype=bool)[None]
            close = both & (np.abs(t0[:, :, None] - t0[:, None, :]) < MERGE_NS)
            if not close.any():
                break
            redo = close.any(axis=(1, 2))
            for ci in np.nonzero(redo)[0]:
                while True:
                    pairs = np.argwhere(close[ci] & valid[ci][:, None] & valid[ci][None, :])
                    if pairs.size == 0:
                        break
                    i, j = pairs[0]
                    valid[ci, i if a[ci, i] < a[ci, j] else j] = False
            a = np.where(valid, a, 0.0).astype(F32)
        a, valid, S, _, db = self.solve(y, w, m, t0, valid, ns, sig)
        return t0, a, valid, np.matmul(a[:, None, :], S)[:, 0, :] + db[:, None], db

    def shoulder_candidate(self, r, w, m, sig, valid, t0, a):
        """Best matched-filter onset of the residual r at >= MIN_SEP_NS from every pulse, and whether its amplitude
        passes A >= max(NSIGMA sigma_A, MIN_AMP) and A >= R_MODEL x every pulse within one kernel span (DT L) of it."""
        Sg, Ag, Ig = self.score(r, w, m)
        onsets = np.arange(Sg.shape[1], dtype=F32) - 3.0
        for k in range(KMAX):
            too_close = valid[:, k, None] & (np.abs(onsets[None, :] - t0[:, k, None]) < MIN_SEP_NS)
            Sg = np.where(too_close, -np.inf, Sg)
        row = np.arange(r.shape[0])
        q = np.argmax(Sg, axis=1)
        tn, A = self.onset(Sg, row, q), Ag[row, q]
        ok = np.isfinite(Sg[row, q]) & (A >= np.maximum(NSIGMA * sig * Ig[row, q], MIN_AMP))
        ok &= ~(valid & (np.abs(tn[:, None] - t0) < MIN_SEP_NS)).any(axis=1)
        ok &= A >= R_MODEL * np.where(valid & (np.abs(tn[:, None] - t0) < DT * L), a, 0.0).max(axis=1)
        return tn, ok

    def add_shoulders(self, y, w, m, sig, t0, a, valid, model, db):
        """Stage C, one pulse at a time: the shoulder_candidate of the residual is added, all pulses are refitted
        (fit_all), and the new pulse is kept only if the BIC improves, chi2_old - chi2_new > BIC_PENALTY 2 ln N_w (sigma(t, A)
        variance of the old pulses, N_w = unclipped samples in the union of both fits' windows).  A row stops
        at its first rejected candidate or at KMAX pulses."""
        ns = y.shape[1]
        wt = np.ones_like(y) if w is None else w
        active = valid.sum(1) < KMAX
        while active.any():
            g = np.nonzero(active)[0]
            old = model[g]
            tn, ok = self.shoulder_candidate((y[g] - old) * wt[g], None if w is None else w[g], m[g], sig[g],
                                             valid[g], t0[g], a[g])
            g, tn, old = g[ok], tn[ok], old[ok]
            wg = None if w is None else w[g]
            t1, v1 = t0[g].copy(), valid[g].copy()
            free = np.argmin(v1, axis=1)                               # first empty slot
            t1[np.arange(g.size), free], v1[np.arange(g.size), free] = tn, True
            t1, a1, v1, new, db1 = self.fit_all(y[g], wg, m[g], sig[g], t1, v1, ns)
            win = (self.windows(t0[g], valid[g], ns) | self.windows(t1, v1, ns)) * wt[g]
            var = self.variance(m[g], t0[g], a[g], valid[g], sig[g], ns)
            dchi2 = (win * ((y[g] - old) ** 2 - (y[g] - new) ** 2) / var).sum(1)
            better = (v1.sum(1) > valid[g].sum(1)) & (dchi2 > BIC_PENALTY * 2.0 * np.log(np.maximum(win.sum(1), 1.0)))
            acc = g[better]
            t0[acc], a[acc], valid[acc], model[acc], db[acc] = t1[better], a1[better], v1[better], new[better], db1[better]
            active[:] = False
            active[acc] = valid[acc].sum(1) < KMAX
        return t0, a, valid, db

    def find_pulses(self, y, w, m, sig):
        """Stages A, B, C on pedestal-subtracted rows y [n, ns] (clipped samples zeroed and weight w = 0 there; w None
        without clipping): (t0, a, valid, db) [n, KMAX], db the fitted baseline offset [n]."""
        t0, valid = self.find_candidates(y, w, m, sig)
        t0, a, valid, model, db = self.fit_all(y, w, m, sig, t0, valid, y.shape[1])
        return self.add_shoulders(y, w, m, sig, t0, a, valid, model, db)

    def bad_clipped(self, v, t0, w, m, ns):
        """Clipped rows with < 6 unclipped window samples or < 2 on a leading edge, or no pulse."""
        bad = np.zeros(v.shape[0], dtype=bool)
        tt = DT * np.arange(ns)[None, :]
        for k in range(KMAX):
            rows = np.nonzero(v[:, k])[0]
            if rows.size:
                dtk = tt - t0[rows, k][:, None]
                lead = (dtk > -DT) & (dtk <= self.b.t_pk[m[rows]][:, None])
                bad[rows] |= ((self.wmask(t0[rows, k], ns) * w[rows]).sum(1) < 6) | ((lead * w[rows]).sum(1) < 2)
        return bad | ~v.any(axis=1)

    def run(self, Y, midx, overflow):
        """Waveforms Y [n, ns] -> ped, sigma [n], t0, a, valid [n, KMAX], sample weights [n, ns], dped [n], truncated [n, KMAX].
        ped = the 30-sample baseline + dped, dped the fitted baseline offset (|dped| <= BASE_CLIP sigma); sigma is the 30-sample
        noise; truncated: the pulse's LW-sample fit window is cut by the end of the waveform."""
        Yf = Y.astype(F32)
        n, ns = Yf.shape
        x = Yf[:, :PED_NSAMPLES]
        med = np.median(x, axis=1)
        sig = np.maximum(1.4826 * np.median(np.abs(x - med[:, None]), axis=1), SIGMA_MIN).astype(F32)
        ok = np.abs(x - med[:, None]) <= np.maximum(2.0 * sig, 1.0)[:, None]
        ped = ((x * ok).sum(axis=1) / np.maximum(ok.sum(axis=1), 1)).astype(F32)
        y = Yf - ped[:, None]
        clipped = (Yf >= overflow).any(axis=1)
        cand = ((y[:, :-2] + y[:, 1:-1] + y[:, 2:]) * (1.0 / 3.0)).max(axis=1) > np.maximum(CAND_NSIGMA * sig, 3.0)
        W = (Yf < overflow).astype(F32)
        T0, A, V = np.zeros((n, KMAX), F32), np.zeros((n, KMAX), F32), np.zeros((n, KMAX), bool)
        DB = np.zeros(n, F32)
        for masked in (False, True):
            sel = np.nonzero(cand & (clipped == masked))[0]
            if sel.size == 0:
                continue
            ys, ms, w = y[sel], midx[sel], W[sel] if masked else None
            if masked:
                ys = ys * w
            t0, a, v, db = self.find_pulses(ys, w, ms, sig[sel])
            if masked:
                v &= ~self.bad_clipped(v, t0, w, ms, ns)[:, None]
            T0[sel], A[sel], V[sel], DB[sel] = t0, np.where(v, a, 0.0), v, np.where(v.any(axis=1), db, 0.0)
        trunc = V & (np.ceil(T0 / DT).astype(np.int64) - LEAD + LW > ns)
        return ped + DB, sig, T0.astype(np.float64), A.astype(np.float64), V, W, DB, trunc


def fit(idx):
    """Reconstruct rows idx of D: ped, t0, a, valid, model [n, ns] (+ped), chi2/dof, dof."""
    Y = D["waveforms"][idx]
    ped, _, t0, a, v, W, _, _ = D["mf"].run(Y, D["midx"][idx], D["overflow"])
    ns = Y.shape[1]
    S = np.stack([D["mf"].basis(D["midx"][idx], t0[:, k].astype(F32), ns) for k in range(KMAX)], axis=1)
    model = ped[:, None] + np.einsum("nk,nks->ns", np.where(v, a, 0.0), S)
    mask = MF.windows(t0, v, ns) & (W > 0)
    noise = np.sqrt(D["noise_sigma"][idx].astype(np.float64) ** 2 + 1.0 / 12.0)   # + rounding variance
    var = D["mf"].variance(D["midx"][idx], t0.astype(F32), a, v, noise, ns)
    dof = mask.sum(1) - 2 * v.sum(1) - 1                  # 2 per pulse + the baseline offset
    with np.errstate(invalid="ignore", divide="ignore"):
        chi2 = np.where(v.any(1) & (dof > 0), (mask * (Y - model) ** 2 / var).sum(1) / dof, np.nan)
    return dict(ped=ped, t0=t0, a=a, v=v, S=S, model=model, chi2=chi2, dof=dof)


def match(idx, f=None):
    """Per event of D rows idx: npk, timing / timing+amp flags (4 ns, 2 ns), t0 - t2, A / a2 - 1, chi2."""
    f = fit(idx) if f is None else f
    t0, a, v, ar = f["t0"], f["a"], f["v"], np.arange(idx.size)
    t1, t2, a2 = (D[k][idx].astype(np.float64) for k in ("t1_ns", "t2_ns", "a2"))
    im = np.argmin(np.where(v, np.abs(t0 - t1[:, None]), np.inf), axis=1)       # main pulse
    d2 = np.where(v & np.isfinite(t2)[:, None], np.abs(t0 - t2[:, None]), np.inf)
    d2[ar, im] = np.inf
    js = np.argmin(d2, axis=1)                                                  # secondary
    has = np.isfinite(d2[ar, js])
    with np.errstate(invalid="ignore", divide="ignore"):
        dt = np.where(has, t0[ar, js] - t2, np.nan)
        da = np.where(has, a[ar, js] / a2 - 1.0, np.nan)
        ft, ft2, fa = np.abs(dt) <= T_TOL_NS, np.abs(dt) <= T_TOL2_NS, np.abs(da) <= A_TOL
    with np.errstate(invalid="ignore", divide="ignore"):
        pull = np.where(v.any(1), a[ar, im] / D["a1"][idx] - 1.0, np.nan)                 # (A1_fit - A1) / A1
    pb = (f["ped"] - D["pedestal"][idx]) / D["noise_sigma"][idx]                     # pedestal bias in sigma
    return dict(npk=v.sum(1), ped_bias=pb, ft=ft, fta=ft & fa, ft2=ft2, fta2=ft2 & fa, dt=dt, da=da, pull=pull, chi2=f["chi2"])


def med_mad(x):
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan, np.nan
    m = np.median(x)
    return m, 1.4826 * np.median(np.abs(x - m))


def tabulate(d, r):
    """List of TSV row dicts per (module, a1, ratio, dt)."""
    rows = []
    cells, inv = np.unique(np.stack([d["midx"], d["a1"], d["ratio"], d["dt_ns"]], axis=1), axis=0, return_inverse=True)
    order = np.argsort(inv.ravel(), kind="stable")
    bounds = np.searchsorted(inv.ravel()[order], np.arange(len(cells) + 1))
    for c, (m, a1, ratio, dtn) in enumerate(cells):
        s = order[bounds[c]:bounds[c + 1]]
        ch = r["chi2"][s]
        row = dict(module=d["names"][int(m)], a1=a1, ratio=ratio, dt_ns=dtn, n=s.size,
                   chi2=np.median(ch[np.isfinite(ch)]) if np.isfinite(ch).any() else np.nan)
        if ratio == 0:
            row["nfp"] = int((r["npk"][s] >= 2).sum())
            row["pull"] = r["pull"][s[r["npk"][s] == 1]]                    # control: events with one pulse
        else:
            for k in ("ft", "fta", "ft2", "fta2"):
                row[k] = int(r[k][s].sum())
            sm = s[r["ft"][s]]
            row["dt_med"], row["dt_mad"] = med_mad(r["dt"][sm])
            row["da_med"], row["da_mad"] = med_mad(r["da"][sm])
            row["pull"] = r["pull"][s[r["fta"][s]]]                         # found events
        row["pull_med"], row["pull_mad"] = med_mad(row["pull"])
        row["ped_bias"] = np.median(r["ped_bias"][s])
        rows.append(row)
    return rows


def write_tsv(rows, path):
    cols = ["module", "a1", "ratio", "dt_ns", "n", "eff_t", "err_t", "eff_ta", "err_ta", "eff_t_2ns",
            "eff_ta_2ns", "fp_rate", "dt_med_ns", "dt_mad_ns", "da_med", "da_mad", "chi2dof_med", "pull_med", "pull_mad", "ped_bias"]
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for w in rows:
            e = {k: w.get(k, np.nan) / w["n"] for k in ("ft", "fta", "ft2", "fta2")}
            err = {k: np.sqrt(e[k] * (1 - e[k]) / w["n"]) for k in ("ft", "fta")}
            v = [w["module"], f"{w['a1']:g}", f"{w['ratio']:g}", f"{w['dt_ns']:g}", w["n"],
                 e["ft"], err["ft"], e["fta"], err["fta"], e["ft2"], e["fta2"], w.get("nfp", np.nan) / w["n"]]
            v += [w.get(k, np.nan) for k in ("dt_med", "dt_mad", "da_med", "da_mad", "chi2", "pull_med", "pull_mad", "ped_bias")]
            f.write("\t".join(x if isinstance(x, str) else f"{x:.6g}" for x in v) + "\n")


def grid(rows, module, ratios, dts, a1=None):
    """(eff %, binomial err %, n) [ratio, dt] of the timing + amplitude efficiency, pooled over A1 unless a1 given."""
    k, n = np.zeros((len(ratios), len(dts))), np.zeros((len(ratios), len(dts)))
    for w in rows:
        if w["module"] != module or w["ratio"] == 0 or (a1 is not None and w["a1"] != a1):
            continue
        i, j = ratios.index(w["ratio"]), dts.index(w["dt_ns"])
        k[i, j] += w["fta"]
        n[i, j] += w["n"]
    with np.errstate(invalid="ignore", divide="ignore"):
        e = k / n
        return 100.0 * e, 100.0 * np.sqrt(e * (1 - e) / n), n


def draw(ax, eff, err, ratios, dts, title, fs):
    ax.imshow(eff, cmap="Greens", vmin=0, vmax=100, aspect="auto")
    for (i, j), v in np.ndenumerate(eff):
        ax.text(j, i, "-" if np.isnan(v) else f"{v:.1f} ± {err[i, j]:.1f}", ha="center", va="center",
                fontsize=fs, color="white" if v > 70 else "black")
    ax.set_xticks(range(len(dts)), [f"{x:g}" for x in dts])
    ax.set_yticks(range(len(ratios)), [f"{x:g}" for x in ratios])
    ax.set_xlabel("time gap dt = t2 - t1 (ns)")
    ax.set_ylabel("amplitude ratio r = A2 / A1")
    ax.set_title(title)


def draw_pull(rows, ratios, dts, a1s, plt):
    """Main-pulse amplitude pull (A1_fit - A1) / A1 of the found events, pooled over modules: median with a bar of
    +-1.4826 MAD vs dt, one line per r, one panel per A1."""
    fig, axs = plt.subplots(len(a1s), 1, figsize=(17, 20), sharex=True)
    for ax, a in zip(axs, a1s):
        for i, rt in enumerate(ratios):
            st = [med_mad(np.concatenate([np.empty(0)] + [w["pull"] for w in rows if w["ratio"] == rt and w["a1"] == a and w["dt_ns"] == d]))
                  for d in dts]
            ax.errorbar(np.arange(len(dts)) + 0.08 * (i - len(ratios) / 2), [100 * x[0] for x in st], [100 * x[1] for x in st],
                        fmt="o-", ms=4, lw=1, capsize=2, label=f"r = {rt:g}")
        ax.axhline(0, color="k", lw=0.5)
        ax.set_xticks(range(len(dts)), [f"{x:g}" for x in dts])
        ax.set_ylabel("main-pulse pull (%)")
        ax.set_title(f"A1 = {a:g} ADC")
        ax.legend(ncol=7)
    axs[-1].set_xlabel("time gap dt = t2 - t1 (ns)")
    fig.suptitle("main-pulse amplitude pull (A1_fit - A1) / A1 of the found events, all modules:\nmedian ± 1.4826 MAD")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return fig


def write_pdf(rows, meta, path, plt, PdfPages):
    f32 = lambda v: [float(np.float32(x)) for x in v]     # noqa: E731  (npz grids are float32)
    ratios, dts, a1s = f32(meta["ratios"]), f32(meta["dt_ns"]), f32(meta["a1_adc"])
    defn = (f"secondary found: |t0 - t2| <= {T_TOL_NS:g} ns and |A/A2 - 1| <= {A_TOL * 100:g} %;\n"
            f"cell = eff ± binomial error (%)  (matched filter + template fit, r_model = {R_MODEL:g})")
    with PdfPages(path) as pdf:
        for by_a1 in (False, True):
            for m in [x["name"] for x in meta["modules"]]:
                ctl = [w for w in rows if w["module"] == m and w["ratio"] == 0]
                nfp, nctl = sum(w["nfp"] for w in ctl), sum(w["n"] for w in ctl)
                if not by_a1:
                    fig, ax = plt.subplots(figsize=(18, 7.5))
                    eff, err, n = grid(rows, m, ratios, dts)
                    draw(ax, eff, err, ratios, dts, f"{m}: deconvolution efficiency (%), pooled over "
                         f"A1 = {', '.join(f'{a:g}' for a in a1s)} ADC  ({int(n.max())} events per cell)", 7)
                else:
                    fig, axs = plt.subplots(3, 1, figsize=(18, 19))
                    for ax, a in zip(axs, a1s):
                        eff, err, n = grid(rows, m, ratios, dts, a)
                        draw(ax, eff, err, ratios, dts, f"by A1: {m}, A1 = {a:g} ADC ({int(n.max())}/cell)", 7)
                    fig.suptitle(f"{m}: deconvolution efficiency (%) by A1")
                fig.text(0.01, 0.005, defn + f"\nr = 0 control false-positive rate (>= 2 pulses): "
                         f"{nfp}/{nctl} = {100 * nfp / max(nctl, 1):.3f} %\n(" + ", ".join(
                             f"A1 {w['a1']:g}: {w['nfp']}/{w['n']}" for w in ctl) + ")", fontsize=14)
                fig.tight_layout(rect=(0, 0.13 if not by_a1 else 0.07, 1, 0.97))
                pdf.savefig(fig)
                plt.close(fig)
        fig = draw_pull(rows, ratios, dts, a1s, plt)
        pdf.savefig(fig)
        plt.close(fig)


def write_examples(sel, res, path, plt, PdfPages):
    """One page per ratio r (4 x 5 panels, one per gap dt): the first FOUND event of each (r, dt) cell, or its
    first event labelled MISSED when the cell has no resolved event.  Final page: the r = 0 control."""
    rs, ds = D["ratio"][sel], D["dt_ns"][sel]
    ratios, dts = np.unique(rs[rs > 0]), np.unique(ds[rs > 0])
    cell = (D["module"][sel] == EX_MODULE) & (D["a1"][sel] == F32(EX_A1))

    def pick(r, d):                                       # (event index, found) or None
        m = np.nonzero(cell & (rs == r) & (ds == d))[0]   # cell members, positions in sel
        hit = m[res["fta"][m]] if r > 0 else m
        q = hit if hit.size else m
        return (int(sel[q[0]]), hit.size > 0) if q.size else None

    pages = [(r, [(d, pick(r, d)) for d in dts]) for r in ratios] + [(0.0, [(0.0, pick(0.0, 0.0))])]
    got = [p[0] for _, row in pages for _, p in row if p is not None]
    if not got:
        print(f"[examples] no events for {EX_MODULE}, A1 = {EX_A1:g}; {path} not written")
        return
    idx = np.array(got)
    f = fit(idx)
    mt = match(idx, f)
    pos = {e: p for p, e in enumerate(got)}
    tt = DT * np.arange(D["waveforms"].shape[1])
    nc = 5
    with PdfPages(path) as pdf:
        for r, row in pages:
            ctl = r == 0
            nr, nk = (1, 1) if ctl else (-(-len(row) // nc), nc)
            fig, axs = plt.subplots(nr, nk, figsize=(14, 8.5) if ctl else (6.6 * nk, 6.4 * nr), squeeze=False)
            for ax in axs.flat:
                ax.axis("off")
            for ax, (d, hit_) in zip(axs.flat, row):
                ax.axis("on")
                if hit_ is None:
                    ax.set_title(f"{EX_MODULE}  A1 = {EX_A1:g}  r = {r:g}  dt = {d:g} ns\nno events in cell", fontsize=11)
                    continue
                e, found = hit_
                p = pos[e]
                ax.plot(tt, D["waveforms"][e], "k.", ms=3, label="data")
                ax.plot(tt, f["model"][p], "r-", lw=1, label="fit")
                ks = np.nonzero(f["v"][p])[0][np.argsort(f["t0"][p][f["v"][p]])]        # pulses in time order
                for c, j in zip(("tab:blue", "tab:orange", "tab:green", "tab:purple", "tab:brown", "tab:pink"), ks):
                    ax.plot(tt, f["ped"][p] + f["a"][p, j] * f["S"][p, j], "--", color=c, lw=0.8)
                ax.axvline(D["t1_ns"][e], color="tab:blue", ls=":", lw=1)
                if not ctl:
                    ax.axvline(D["t2_ns"][e], color="tab:orange", ls=":", lw=1)
                ax.set_xlim(0, 4 * D["waveforms"].shape[1])
                pk = "; ".join(f"{f['t0'][p, j]:.1f}/{f['a'][p, j]:.0f}" for j in ks)
                if ctl:
                    res_txt = f"control: {mt['npk'][p]} pulse(s)"
                elif not found:
                    res_txt = "MISSED (no resolved event in cell)"
                else:
                    res_txt = f"secondary {'FOUND' if mt['fta'][p] else 'MISSED'}"
                ax.set_title(f"{D['names'][D['midx'][e]]}  A1 = {D['a1'][e]:g}  r = {r:g}  dt = {d:g} ns\n"
                             f"truth t1/t2 = {D['t1_ns'][e]:.1f}/{'%.1f' % D['t2_ns'][e] if not ctl else '-'} ns, A2 = {D['a2'][e]:.0f}\n"
                             f"fit t0/A: {pk}\nchi2/dof = {f['chi2'][p]:.2f} (dof {f['dof'][p]})  {res_txt}", fontsize=11)
                ax.set_xlabel("t (ns)")
            fig.suptitle(f"{EX_MODULE}, A1 = {EX_A1:g} ADC, r = {r:g}" + (" (control)" if ctl else ""))
            fig.text(0.01, 0.005, "black: data; red: fitted model; dashed: fitted pulses in time order (blue, orange, ...);\n"
                     "dotted: truth t1 (blue), t2 (orange).  Fit t0/A in time order.\n"
                     "Found = |t0 - t2| <= 4 ns and |A/A2 - 1| <= 30 %.  Each panel: first FOUND event of its (r, dt) cell;\n"
                     "a cell with no resolved event shows its first event, labelled MISSED.",
                     fontsize=14)
            fig.tight_layout(rect=(0, 0.11, 1, 0.96))
            pdf.savefig(fig)
            plt.close(fig)


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    z = np.load(sys.argv[1], allow_pickle=False)
    meta = json.loads(str(z["meta"]))
    names = [m["name"] for m in meta["modules"]]
    for k in ("waveforms", "module", "a1", "a2", "ratio", "dt_ns", "t1_ns", "t2_ns", "noise_sigma", "pedestal"):
        D[k] = z[k]
    D["midx"] = np.zeros(D["module"].size, dtype=np.int64)
    for k, nm in enumerate(names):
        D["midx"][D["module"] == nm] = k
    D["names"], D["overflow"], D["mf"] = names, int(meta["overflow"]), MF(Bank(meta, names))
    sel = np.arange(D["a1"].size)
    chunks = [sel[i:i + CHUNK] for i in range(0, sel.size, CHUNK)]
    print(f"[run] {sel.size} events, {len(chunks)} chunks, {NPROC} processes", flush=True)
    with mp.get_context("fork").Pool(NPROC) as pool:
        parts = []
        for i, p in enumerate(pool.imap(match, chunks)):
            parts.append(p)
            print(f"[run] chunk {i + 1}/{len(chunks)}", flush=True)
    r = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    rows = tabulate({k: D[k][sel] for k in ("midx", "a1", "ratio", "dt_ns")} | {"names": names}, r)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    # plot style: large bold axis labels, large legends
    plt.rcParams.update({
        "axes.labelsize": 16, "axes.labelweight": "bold", "axes.titlesize": 16, "axes.titleweight": "bold",
        "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 14, "axes.linewidth": 1.8,
        "xtick.major.width": 1.5, "ytick.major.width": 1.5, "xtick.major.size": 6, "ytick.major.size": 6,
        "lines.linewidth": 2.0, "figure.titlesize": 18, "figure.titleweight": "bold"})
    from matplotlib.backends.backend_pdf import PdfPages
    out = {k: os.path.join(OUT_DIR, k) for k in ("efficiency.tsv", "efficiency_tables.pdf", "pileup_examples.pdf")}
    write_tsv(rows, out["efficiency.tsv"])
    write_pdf(rows, meta, out["efficiency_tables.pdf"], plt, PdfPages)
    write_examples(sel, r, out["pileup_examples.pdf"], plt, PdfPages)
    print(f"[done] {', '.join(out.values())}")


if __name__ == "__main__":
    main()
