#!/usr/bin/env python3
"""crosstalk_study.py -- matched-filter chi-squared per degree of freedom of the seed and neighbour modules of HyCal
clusters in slices of cluster energy, and the constant-fraction-discriminator time of every fitted channel.

    python3 crosstalk_study/crosstalk_study.py scan <out.npz> <evio file> [<evio file> ...]
    python3 crosstalk_study/crosstalk_study.py report <scan dir> <out dir>

scan    all splits of one run: matched filter on every PbWO4 channel above threshold of the gated events; stores per
        channel the fitted pulses, chi-squared per degree of freedom, pedestal and constant-fraction-discriminator time.
report  pools <scan dir>/run_*.npz per run type, forms the clusters of every event and writes the figures and tables.
Chain, definitions and outputs: README.md.  Standalone: numpy, scipy, matplotlib and the compiled prad2py decoder;
source setup_env.sh first ($PRAD2_DATABASE_DIR).
Bad waveforms (any sample equal to 0 or UNDERSHOOT_ADC below the baseline) are dropped before the matched filter and channels whose
main pulse has its fit window cut by the end of the waveform are dropped from every amplitude quantity; both are counted per module
(mod_n_bad, mod_n_trunc) and printed at the end of the scan.
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")                     # one BLAS thread per process

import json
import re
import sys
import time
import warnings
from pathlib import Path

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks

HERE = Path(__file__).resolve().parent
PYDIR = HERE.parent                                     # analysis/pyscripts
BANK_ROOT = PYDIR / "waveform_analysis" / "pulse_templates"
DB_DIR = Path(os.environ.get("PRAD2_DATABASE_DIR", "")) # database: always $PRAD2_DATABASE_DIR (setup_env.sh exports it)
HYCAL_MAP = DB_DIR / "hycal_map.json"
RUNINFO = DB_DIR / "runinfo" / "general.json"          # per-run beam energy, calibration file, dead modules
RUN_TYPE_OF = {
    24590: "0p73GeV_H2", 24603: "0p73GeV_H2",
    25228: "2p24GeV_H2", 25257: "2p24GeV_H2", 25260: "2p24GeV_H2", 25299: "2p24GeV_H2", 25320: "2p24GeV_H2",
    24698: "3p49GeV_H2", 24911: "3p49GeV_H2", 24937: "3p49GeV_H2",
    25800: "2p24GeV_Ta_X17", 26138: "2p24GeV_Ta_X17", 26230: "2p24GeV_Ta_X17",
}
RUN_TYPES = ["0p73GeV_H2", "2p24GeV_H2", "3p49GeV_H2", "2p24GeV_Ta_X17"]
H2_LABELS = {"0p73GeV_H2": "0.73 GeV H2", "2p24GeV_H2": "2.24 GeV H2", "3p49GeV_H2": "3.49 GeV H2"}   # run types of the figures
PHYSICS_ANY = {"2p24GeV_Ta_X17"}       # trigger: (bits & 0xf00) != 0; others bits == 0x100
PHYSICS_BITS_EXACT, PHYSICS_MASK_ANY = 0x100, 0xf00
LMS_BIT = 0x1000000                    # light-monitoring trigger bit: these events are vetoed
NS_EXPECT, OVERFLOW, BLOCK = 100, 4095, 256   # samples per waveform, ADC overflow, events per fit block
MAX_EVENTS = int(os.environ.get("CROSSTALK_MAX_EVENTS", "0"))   # test aid only: stop the scan after this many decoded events (0 = all)
KEEP_ADC = 10.0   # kept channel: fitted main amplitude (baseline-subtracted) >= this; ≈ 9.5 MeV at the median gain; Møller partners and shower neighbours
PRECUT_ADC = 8.0   # matched filter only on channels with raw max - baseline >= this
UNDERSHOOT_ADC = 20.0   # a sample this far below the baseline, or any sample equal to 0, marks a bad waveform
REC = dict(midx=np.int16, ev=np.int32, k=np.int8, a=np.float32, t=np.float32, chi2=np.float32,
           ped=np.float32, sig=np.float32, tcfd=np.float32, nclip=np.int8)   # per kept channel
CFD_FRAC = 0.5         # constant-fraction-discriminator fraction of the peak sample
CFD_SEARCH_NS = 12.0   # peak sample searched within this of the main pulse's template peak (3 samples)
PRE_LO, PRE_HI, CL_RADIUS = 130.0, 200.0, 1.5   # in-window onset range (ns, benchmark default); cluster radius / mean module size
E_SLICES = [10, 50, 100, 200, 500, 1000, 2000, 3000, 4000]   # cluster-energy slice edges (MeV)
EX_CHI2, EX_KEEP, EX_PER_SLICE = 5.0, 12, 6   # examples: seeds with chi-squared per degree of freedom > 5; kept per slice and run; drawn per slice
EX_RUN_TYPE = "2p24GeV_H2"             # run type of the examples PDF
CHI2_FINE = np.r_[0.0, np.logspace(-2, 4, 1201), np.inf]   # fine log binning (200 per decade) for histograms and percentiles
CHI2_LO, CHI2_HI = 0.1, 100.0          # plotted chi-squared per degree of freedom range (log)
DT_EDGES = np.arange(-10.0, 30.01, 0.5)   # bins of t_cfd - t_mf (ns) for the profile
MIN_ENTRIES = 20       # summary curves and profile bins: points with fewer entries are not drawn

# --- matched filter: class MF and its constants copied verbatim from waveform_analysis/pileup_analysis/pileup_deconvolution_study.py;
# class Bank copied from waveform_analysis/pileup_analysis/pileup_rates.py (the data version of the study's Bank: run-type bank file,
# average spline for modules without their own template)
PED_NSAMPLES = 30     # baseline samples (0-116 ns; the main onset is at >= 152 ns)
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
F32 = np.float32


class Bank:
    """Spline bank of one run type on the 0.05 ns grid, linear interpolation.  LUT row 0 is the
    bank's average spline (lead glass, PbWO4 without a valid row or dead); row[midx] maps modules to rows.  slut: the
    row's shape_sigma on the same grid (linear in the bank's t_grid, ends held), at least SHAPE_SIGMA_MIN; aband: the
    row's a_band; row 0: the median over the bank's valid, not dead rows."""

    def __init__(self, path, names, pwo4):
        from scipy.interpolate import BSpline
        z = np.load(path, allow_pickle=False)
        k, self.T_end = int(z["spline_k"]), float(z["T_end"])
        self.n = int(round(self.T_end / LUT_DT)) + 1
        tg = LUT_DT * np.arange(self.n)
        at, ac = z["avg_spline_t"], z["avg_spline_c"]
        knots = [(at[np.isfinite(at)], ac[np.isfinite(ac)])]
        tb, ss = np.asarray(z["t_grid"], dtype=np.float64), z["shape_sigma"]
        if "a_band" not in z.files:
            raise SystemExit(f"[ERROR] {path} has no a_band: rebuild the bank with the current template_extract.py")
        ab = z["a_band"]
        dead = z["dead"] > 0 if "dead" in z.files else np.zeros(len(z["names"]), dtype=bool)
        use = (z["valid"] > 0) & ~dead
        with warnings.catch_warnings():                  # all-NaN columns (outside the template) stay NaN
            warnings.simplefilter("ignore", RuntimeWarning)
            sig = [np.nanmedian(ss[use], axis=0)]
        aband = [float(np.median(ab[use]))]
        row = {str(s): i for i, s in enumerate(z["names"])}
        self.row = np.zeros(len(names), dtype=np.int64)
        for i, (nm, p) in enumerate(zip(names, pwo4)):
            r = row.get(str(nm))
            if p and r is not None and use[r]:
                self.row[i] = len(knots)
                knots.append((z["spline_t"][r, :z["spline_nt"][r]], z["spline_c"][r, :z["spline_nc"][r]]))
                sig.append(ss[r])
                aband.append(ab[r])
        lut, dlut, slut = [], [], []
        for (t, c), sr in zip(knots, sig):
            s = BSpline(t, c, k)
            x = np.minimum(tg, t[-k - 1])                # hold the end value beyond the base interval
            lut.append(s(x))
            dlut.append(s.derivative()(x))
            ok = np.isfinite(sr)
            slut.append(np.fmax(np.interp(tg, tb[ok], sr[ok]), SHAPE_SIGMA_MIN))   # floor
        self.lut, self.dlut, self.slut = (np.asarray(v, dtype=np.float64).astype(F32) for v in (lut, dlut, slut))
        self.aband = np.asarray(aband, dtype=F32)          # amplitude (ADC) at which shape_sigma was measured
        self.t_pk = np.argmax(self.lut, axis=1) * LUT_DT
        self.lut, self.dlut, self.slut = self.lut.reshape(-1), self.dlut.reshape(-1), self.slut.reshape(-1)
        mi = np.arange(len(knots))
        K = self.shape(mi[:, None, None], DT * np.arange(L)[None, None, :]
                       + np.arange(NPH)[None, :, None], np.float64)   # K[m, phi, l] = s(4 l + phi)
        self.K, self.K2, self.Kn2 = K.astype(F32), (K * K).astype(F32), (K * K).sum(2).astype(F32)
        self.Kinv = (1.0 / np.sqrt(np.maximum(self.Kn2, 1e-12))).astype(F32)

    def shape(self, m, t, dtype=F32, slope=False, sigma=False):
        """Unit-peak shape (and slope; sigma=True: shape_sigma) at LUT row m, t ns after onset; 0 outside (0, T_end)."""
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


# noise reference ==============================================================
N_REF_EVENTS = 20000      # first decoded events of an EVIO file that set its per-module noise reference; sparse readout: ~6 waveforms per module per 2000 events
MIN_REF_WAVEFORMS = 20    # modules with fewer reference waveforms take the median sigma_ref of the other modules
START_NSIGMA = 3.0        # a 30-sample baseline farther than this many sigma_ref from b_ref starts the baseline fit at b_ref


def pedestal(Y):
    """Clipped mean of the first PED_NSAMPLES samples of each waveform of Y [n, ns]: samples within 2 sigma (at least 1 ADC) of the
    median, sigma = 1.4826 MAD of those samples."""
    x = Y[:, :PED_NSAMPLES].astype(F32)
    med = np.median(x, axis=1)
    sig = np.maximum(1.4826 * np.median(np.abs(x - med[:, None]), axis=1), SIGMA_MIN)
    ok = np.abs(x - med[:, None]) <= np.maximum(2.0 * sig, 1.0)[:, None]
    return ((x * ok).sum(axis=1) / np.maximum(ok.sum(axis=1), 1)).astype(F32)


def noise_reference(X, mi, n_mod):
    """Per-module noise reference from the pedestal samples X [n, PED_NSAMPLES] of n waveforms of modules mi: n_ref waveforms,
    b_ref = median of their baselines (pedestal), sigma_ref = 1.4826 MAD of the pooled residuals (sample - its waveform's
    baseline).  Modules with n_ref < MIN_REF_WAVEFORMS get b_ref NaN (no baseline override) and the median sigma_ref of the others."""
    o = np.argsort(mi, kind="stable")
    X, mi = X[o], mi[o]
    n_ref = np.bincount(mi, minlength=n_mod)
    end = np.cumsum(n_ref)
    b_ref, sig_ref = np.full(n_mod, np.nan, F32), np.full(n_mod, np.nan, F32)
    for m in np.nonzero(n_ref >= MIN_REF_WAVEFORMS)[0]:
        x = X[end[m] - n_ref[m]:end[m]]
        ped = pedestal(x)
        r = x - ped[:, None]
        b_ref[m] = np.median(ped)
        sig_ref[m] = max(1.4826 * np.median(np.abs(r - np.median(r))), SIGMA_MIN)
    sig_ref[n_ref < MIN_REF_WAVEFORMS] = np.nanmedian(sig_ref) if (n_ref >= MIN_REF_WAVEFORMS).any() else SIGMA_MIN
    return n_ref, b_ref, sig_ref


def run_ref(mf, Y, midx, overflow, sig, b_ref):
    """MF.run with the noise sigma given per waveform (sig [n], the module's sigma_ref) instead of the 30-sample MAD.  The baseline
    fit starts at the 30-sample baseline, or at b_ref [n] when that is farther than START_NSIGMA sigma from it.  Returns the MF.run
    tuple and far [n], the waveforms whose start was replaced by b_ref."""
    Yf = Y.astype(F32)
    n, ns = Yf.shape
    ped = pedestal(Yf)
    far = np.abs(ped - b_ref) > START_NSIGMA * sig
    ped = np.where(far, b_ref, ped).astype(F32)
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
        t0, a, v, db = mf.find_pulses(ys, w, ms, sig[sel])
        if masked:
            v &= ~mf.bad_clipped(v, t0, w, ms, ns)[:, None]
        T0[sel], A[sel], V[sel], DB[sel] = t0, np.where(v, a, 0.0), v, np.where(v.any(axis=1), db, 0.0)
    trunc = V & (np.ceil(T0 / DT).astype(np.int64) - LEAD + LW > ns)
    return ped + DB, sig, T0.astype(np.float64), A.astype(np.float64), V, W, DB, trunc, far


def run_config(run):
    """Run's settings from database/runinfo/general.json (as template_extract's resolve_calib): beam energy (MeV) and calibration file
    from the "configurations" entries with from_run <= run, later ones overriding (inheritance); HyCal dead module names from the
    "hycal_dead_modules" entry with the largest from_run <= run (no chaining), dead_src = that entry's from_run."""
    j = json.load(open(RUNINFO))
    beam, calib = None, ""
    for c in sorted((c for c in j["configurations"] if "from_run" in c), key=lambda c: c["from_run"]):
        if c["from_run"] <= run:
            beam, calib = c.get("beam_energy", beam), c.get("calibration", {}).get("file", calib)
    if beam is None or not calib:
        sys.exit(f"[ERROR] runinfo has no beam energy or calibration file for run {run}")
    d = max((e for e in j["hycal_dead_modules"] if e["from_run"] <= run), key=lambda e: e["from_run"])
    return dict(beam=float(beam), calib=calib, dead=list(d["dead_modules"]), dead_src=f"runinfo/general.json hycal_dead_modules from_run {d['from_run']}")


def bad_waveform(Y):
    """[n] bool: any sample == 0, or any sample < median(first PED_NSAMPLES) - UNDERSHOOT_ADC."""
    base = np.median(Y[:, :PED_NSAMPLES], axis=1)
    return (Y == 0).any(axis=1) | (Y < base[:, None] - UNDERSHOOT_ADC).any(axis=1)


def clipped(Y):
    """Samples at OVERFLOW per channel Y [n, ns]; a channel is clipped when this is > 0."""
    return (Y >= OVERFLOW).sum(axis=1)


def neighbours(names):
    """Boolean matrix [seed, module]: centre within CL_RADIUS x the mean of the two module sizes."""
    geo = {m["n"]: m["geo"] for m in json.load(open(HYCAL_MAP))}
    x, y, sz = (np.array([geo[n][q] for n in names]) for q in ("x", "y", "sx"))
    dist = np.hypot(x[:, None] - x, y[:, None] - y)
    return dist <= CL_RADIUS * 0.5 * (sz[:, None] + sz)


def module_gain(cfg, rt, names):
    """Per module E / A (MeV per ADC) of a run: calibration factor (database file cfg["calib"] of the run, key "factor" by module
    "name") x I_int ("iint_wa" of the run type's bank, median for modules without a valid row), as template_extract's
    E = A I_int factor.  nan for modules without a factor."""
    fac = {m["name"]: float(m.get("factor", 0.0)) for m in json.load(open(DB_DIR / cfg["calib"]))}
    z = np.load(BANK_ROOT / rt / f"templates_{rt}.npz", allow_pickle=False)
    iint = {str(n): v for n, v, ok in zip(z["names"], z["iint_wa"], z["valid"]) if ok > 0}
    med = float(np.median(list(iint.values())))
    g = np.array([fac.get(str(n), 0.0) * iint.get(str(n), med) for n in names])
    return np.where(g > 0, g, np.nan)



def fit(mf, bank, Y, m, sig, b_ref):
    """Matched filter on channels Y (template rows m), noise sigma = sig [n] (the module's sigma_ref), baseline fit start
    b_ref [n] (run_ref); adds far (start replaced by b_ref).  Pulses in time
    order: t (NaN pad), a (0 pad); model = pedestal + fit [n, ns]; tcfd = constant-fraction-discriminator time of the main pulse."""
    ped, sig, t0, a, v, W, _, trunc, far = run_ref(mf, Y, m, OVERFLOW, sig, b_ref)
    n, ns = Y.shape
    S = np.stack([mf.basis(m, t0[:, k].astype(F32), ns) for k in range(KMAX)], axis=1)
    model = ped[:, None] + np.einsum("nk,nks->ns", np.where(v, a, 0.0), S)
    mask = MF.windows(t0, v, ns) & (W > 0)
    noise = np.sqrt(sig.astype(np.float64) ** 2 + 1.0 / 12.0)          # + rounding variance
    var = mf.variance(m, t0.astype(F32), a, v, noise, ns)
    dof = mask.sum(1) - 2 * v.sum(1) - 1                  # 2 per pulse + the baseline offset
    with np.errstate(invalid="ignore", divide="ignore"):
        chi2 = np.where(v.any(1) & (dof > 0), (mask * (Y - model) ** 2 / var).sum(1) / dof, np.nan)
    trunc_main = trunc[np.arange(n), a.argmax(1)]          # main pulse's fit window is cut by the end of the waveform
    o = np.argsort(np.where(v, t0, np.inf), axis=1)
    vs = np.take_along_axis(v, o, axis=1)
    t = np.where(vs, np.take_along_axis(t0, o, axis=1), np.nan)
    a = np.where(vs, np.take_along_axis(a, o, axis=1), 0.0)
    # constant-fraction discriminator on the main (largest) pulse
    r = np.arange(n)
    y = Y - ped[:, None]
    main = a.argmax(1)
    t_peak = np.nan_to_num(t[r, main]) + bank.t_pk[m]
    near = np.abs(DT * np.arange(ns)[None, :] - t_peak[:, None]) <= CFD_SEARCH_NS
    p = np.where(near, y, -np.inf).argmax(1)
    half = CFD_FRAC * y[r, p]
    below = (y < half[:, None]) & (np.arange(ns)[None, :] < p[:, None])
    j = np.where(below, np.arange(ns)[None, :], -1).max(1)               # last sample below half before the peak
    jc = np.maximum(j, 0)
    nclip = clipped(Y)
    with np.errstate(invalid="ignore", divide="ignore"):
        tcfd = DT * (jc + (half - y[r, jc]) / (y[r, jc + 1] - y[r, jc]))
    tcfd = np.where((vs.sum(1) > 0) & (j >= 0) & (nclip == 0), tcfd, np.nan)
    return dict(ped=ped, sig=sig, k=vs.sum(1), t=t, a=a, chi2=chi2, tcfd=tcfd, model=model, nclip=nclip, trunc_main=trunc_main, far=far)


def neighbour_lists(nbr):
    """Per module the indices of its neighbourhood (itself included), padded with -1 [modules, max size]."""
    nl = np.full((nbr.shape[0], int(nbr.sum(1).max())), -1, dtype=np.int64)
    for m, row in enumerate(nbr):
        idx = np.nonzero(row)[0]
        nl[m, :idx.size] = idx
    return nl


def find_clusters(ev, mi, k, a, t, nbr, nl, gain):
    """Every cluster of each event (see the module docstring).  Returns seed (channel index per cluster), energy (MeV per
    cluster, NaN if the seed has no gain), mem (member channel indices, seed included) and cl (cluster index of each member)."""
    valid = np.arange(KMAX) < k[:, None]
    win = valid & (t > PRE_LO) & (t < PRE_HI)
    a_used = np.where(win, a, -1.0).max(1)                     # largest in-window pulse
    el = np.nonzero(win.any(1))[0]
    if el.size == 0:
        return (np.zeros(0, np.int64),) * 4
    key = ev[el].astype(np.int64) * nl.shape[0] + mi[el]     # (event, module) of the in-window channels
    o = np.argsort(key)
    ks, ka = key[o], a_used[el][o]
    peak = np.ones(el.size, dtype=bool)                       # local maximum: no larger in-window channel in the neighbourhood
    for c in range(nl.shape[1]):
        nb = nl[mi[el], c]
        q = ev[el].astype(np.int64) * nl.shape[0] + nb
        p = np.minimum(np.searchsorted(ks, q), ks.size - 1)
        peak &= ~((nb >= 0) & (ks[p] == q) & (ka[p] > a_used[el]))
    free = np.zeros(ev.size, dtype=bool)                      # in-window channels not yet in a cluster
    free[el] = True
    cand = el[peak]                                           # possible seeds
    cl_of_ev = np.full(int(ev.max()) + 1, -1)
    seed, mem, cl = [], [], []
    n_cl = 0
    while cand.size:                                          # one cluster per event and round, largest free seed first
        o = cand[np.lexsort((-a_used[cand], ev[cand]))]
        s = o[np.r_[True, ev[o][1:] != ev[o][:-1]]]
        cl_of_ev[ev[s]] = n_cl + np.arange(s.size)
        f = np.nonzero(free)[0]
        f = f[cl_of_ev[ev[f]] >= n_cl]                        # free channels of the events with a new seed
        m = f[nbr[mi[s[cl_of_ev[ev[f]] - n_cl]], mi[f]]]       # within CL_RADIUS of the seed (seed included)
        free[m] = False
        seed.append(s), mem.append(m), cl.append(cl_of_ev[ev[m]])
        n_cl += s.size
        cand = cand[free[cand]]
    seed, mem, cl = (np.concatenate(x) for x in (seed, mem, cl))
    e = np.bincount(cl, weights=np.nan_to_num(gain[mi[mem]] * a_used[mem]), minlength=seed.size)
    return seed, np.where(np.isfinite(gain[mi[seed]]), e, np.nan), mem, cl


def keep_random(old, new, per):
    """Uniform random sample of up to `per` entries per slice: the entries with the smallest random keys of old + new."""
    if old is not None:
        new = {q: np.concatenate([old[q], new[q]]) for q in new}
    o = np.lexsort((new["key"], new["slice"]))
    s = new["slice"][o]
    i = np.arange(s.size)
    rank = i - np.maximum.accumulate(np.where(np.r_[True, s[1:] != s[:-1]], i, 0))
    return {q: v[o[rank < per]] for q, v in new.items()}


# scan =========================================================================
def scan(out, evios):
    from prad2py import dec, det                         # compiled decoder binding (PYTHONPATH from setup_env.sh)
    t_start = time.monotonic()
    parsed = [re.search(r"prad_0*(\d+)\.evio\.(\d+)", os.path.basename(f)) for f in evios]
    if not all(parsed):
        sys.exit("[ERROR] file names must be prad_0<run>.evio.<split>")
    runs = {int(p.group(1)) for p in parsed}
    if len(runs) != 1:
        sys.exit(f"[ERROR] files of several runs: {sorted(runs)}")
    run = runs.pop()
    rt = RUN_TYPE_OF.get(run)
    if rt is None:
        print(f"[skip] run {run} is not one of {RUN_TYPES}")
        return
    cfg = run_config(run)
    dead = set(cfg["dead"])
    gate = (lambda b: (b & PHYSICS_MASK_ANY) != 0) if rt in PHYSICS_ANY else (lambda b: b == PHYSICS_BITS_EXACT)
    b = det.PipelineBuilder()                            # DAQ config + HyCal map of the run from $PRAD2_DATABASE_DIR
    b.set_run_number(run)
    core = b.build()
    hycal = core.hycal
    crate_map = {int(r.tag): int(r.crate) for r in core.daq_cfg.roc_tags}
    mods = [hycal.module(i) for i in range(hycal.module_count())]
    names = np.array([md.name for md in mods])
    pwo4 = np.array([bool(md.is_pwo4()) for md in mods])
    if not dead <= set(names):
        sys.exit(f"[ERROR] dead modules not in the HyCal map: {sorted(dead - set(names))}")
    bank = Bank(BANK_ROOT / rt / f"templates_{rt}.npz", names, pwo4)
    mf = MF(bank)
    lut = bank.lut.reshape(-1, bank.n)
    i50 = np.argmax(lut >= 0.5, axis=1)                  # template 50 % rise time after onset, per template row
    t50 = LUT_DT * (i50 - 1 + (0.5 - lut[np.arange(len(lut)), i50 - 1]) / (lut[np.arange(len(lut)), i50] - lut[np.arange(len(lut)), i50 - 1]))
    gain = module_gain(cfg, rt, names)                   # MeV per ADC; NaN: no calibration factor
    nbr = neighbours(names)
    nl = neighbour_lists(nbr)
    rng = np.random.default_rng(run)
    print(f"[setup] run {run} {len(evios)} splits {rt}; beam {cfg['beam']:g} MeV; dead {sorted(dead)}", flush=True)

    rec = {q: [] for q in REC}
    ex = None                                            # example reservoir: seeds with chi-squared per degree of freedom > EX_CHI2
    evnum, rows, evs, mis = [], [], [], []
    n_bad, n_trunc = np.zeros(len(names), np.int64), np.zeros(len(names), np.int64)   # per module: bad waveforms, truncated main pulses

    sig_ref = b_ref = None                                   # noise reference of the file being read
    ref_n, ref_b, ref_sig = [], [], []                       # per file: reference waveforms, b_ref, sigma_ref per module
    n_base_ref = np.zeros(len(names), np.int64)                 # per module: waveforms whose baseline fit started at b_ref

    def flush():
        nonlocal ex
        if not rows:
            return
        Y, ev, mi = np.stack(rows), np.array(evs), np.array(mis)
        rows.clear(), evs.clear(), mis.clear()
        bad = bad_waveform(Y)
        n_bad[:] += np.bincount(mi[bad], minlength=len(names))
        s = np.nonzero(~bad & (Y.max(axis=1) - np.median(Y[:, :PED_NSAMPLES], axis=1) >= PRECUT_ADC))[0]
        f = fit(mf, bank, Y[s], bank.row[mi[s]], sig_ref[mi[s]], b_ref[mi[s]])
        n_base_ref[:] += np.bincount(mi[s][f.pop("far")], minlength=len(names))
        big = f["a"].max(axis=1) >= KEEP_ADC
        n_trunc[:] += np.bincount(mi[s][big & f["trunc_main"]], minlength=len(names))
        keep = big & ~f["trunc_main"]
        s = s[keep]
        f = {q: v[keep] for q, v in f.items()}
        for q, v in dict(midx=mi[s], ev=ev[s], **{q: f[q] for q in REC if q in f}).items():
            rec[q].append(v.astype(REC[q]))
        seed, energy, _, _ = find_clusters(ev[s] - ev[s].min(), mi[s], f["k"], f["a"], f["t"], nbr, nl, gain)   # the block holds whole events
        sl = np.digitize(energy, E_SLICES) - 1
        sel = (f["chi2"][seed] > EX_CHI2) & (sl >= 0) & (sl < len(E_SLICES) - 1)
        pick = seed[sel]
        if pick.size:
            new = dict(samples=Y[s[pick]], model=f["model"][pick].astype(np.float32), midx=mi[s[pick]],
                       evnum=np.array([evnum[e] for e in ev[s[pick]]]), energy=energy[sel], slice=sl[sel], key=rng.random(pick.size),
                       **{q: f[q][pick] for q in ("k", "a", "t", "chi2", "tcfd")})
            ex = keep_random(ex, new, EX_KEEP)

    cache = {}
    def module_of(tag, crate, sl, c):
        """HyCal module index of a DAQ channel, -1 for channels that are not scanned (cached)."""
        mi = cache.get((tag, sl, c))
        if mi is None:
            m = hycal.module_by_daq(crate, sl, c)
            mi = cache[(tag, sl, c)] = m.index if (m is not None and m.is_pwo4() and m.name not in dead) else -1
        return mi

    def read_reference(evio, use):
        """noise_reference of the first N_REF_EVENTS decoded events of an EVIO file, over the modules where use is True."""
        ch = dec.EvChannel()
        ch.set_config(core.daq_cfg)
        if ch.open_auto(str(evio)) != dec.Status.success:
            sys.exit(f"[ERROR] cannot open {evio}")
        X, M, rows, mis = [], [], [], []

        def take():
            if rows:
                Y, mi = np.stack(rows), np.array(mis)
                rows.clear(), mis.clear()
                good = ~bad_waveform(Y)
                X.append(Y[good, :PED_NSAMPLES])
                M.append(mi[good])

        n = 0
        while n < N_REF_EVENTS and ch.read() == dec.Status.success:
            if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
                continue
            for i in range(ch.get_n_events()):
                d = ch.decode_event(i, with_ssp=False)
                if not d["ok"]:
                    continue
                n += 1
                fe = d["event"]
                for ri in range(fe.nrocs):
                    roc = fe.roc(ri)
                    crate = crate_map.get(roc.tag) if roc.present else None
                    if crate is None:
                        continue
                    for sl in roc.present_slots():
                        slot = roc.slot(sl)
                        for c in slot.present_channels():
                            mi = module_of(roc.tag, crate, sl, c)
                            cd = slot.channel(c)
                            if mi >= 0 and use[mi] and cd.nsamples == NS_EXPECT:
                                rows.append(np.asarray(cd.samples, dtype=np.uint16))
                                mis.append(mi)
                if n % 200 == 0:
                    take()
                if n >= N_REF_EVENTS:
                    break
        take()
        if not X:                                                  # no usable waveform in the file: every module on the floor sigma
            return np.zeros(len(names), np.int64), np.full(len(names), np.nan, F32), np.full(len(names), SIGMA_MIN, F32)
        return noise_reference(np.concatenate(X), np.concatenate(M), len(names))

    n_seen = n_trig = n_lms = 0
    done = False
    for evio in evios:
        if done:
            break
        flush()                                                # blocks never span two files: each has its own noise reference
        n_ref, b_ref, sig_ref = read_reference(evio, pwo4)
        ref_n.append(n_ref), ref_b.append(b_ref), ref_sig.append(sig_ref)
        ch = dec.EvChannel()
        ch.set_config(core.daq_cfg)
        if ch.open_auto(str(evio)) != dec.Status.success:
            sys.exit(f"[ERROR] cannot open {evio}")
        while not done and ch.read() == dec.Status.success:
            if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
                continue
            for i in range(ch.get_n_events()):
                d = ch.decode_event(i, with_ssp=False)
                if not d["ok"]:
                    continue
                n_seen += 1
                done = MAX_EVENTS > 0 and n_seen >= MAX_EVENTS
                fe = d["event"]
                bits = int(fe.info.trigger_bits)
                if not gate(bits):
                    if done:
                        break
                    continue
                n_trig += 1
                if bits & LMS_BIT:
                    n_lms += 1
                    if done:
                        break
                    continue
                ev = len(evnum)
                evnum.append(int(fe.info.event_number))
                for ri in range(fe.nrocs):
                    roc = fe.roc(ri)
                    crate = crate_map.get(roc.tag) if roc.present else None
                    if crate is None:
                        continue
                    for sl in roc.present_slots():
                        slot = roc.slot(sl)
                        for c in slot.present_channels():
                            mi = module_of(roc.tag, crate, sl, c)
                            cd = slot.channel(c)
                            if mi >= 0 and cd.nsamples == NS_EXPECT:
                                rows.append(np.asarray(cd.samples, dtype=np.uint16))
                                evs.append(ev)
                                mis.append(mi)
                if len(evnum) % BLOCK == 0:
                    flush()
                if n_seen % 20000 == 0:
                    print(f"  [progress] {n_seen} events, {time.monotonic() - t_start:.0f} s", flush=True)
                if done:
                    break
    flush()
    r = {}
    for q in REC:                                        # one field at a time: peak memory = lists + one field
        v = rec.pop(q)
        r[q] = np.concatenate(v) if v else np.zeros((0, KMAX) if q in ("a", "t") else 0, REC[q])
    wall = time.monotonic() - t_start
    sr = np.array(ref_sig)[np.array(ref_n) >= MIN_REF_WAVEFORMS]
    sr = sr if sr.size else np.array([SIGMA_MIN])
    print(f"[noise reference] sigma_ref median {np.median(sr):.3f} ADC (16-84 %: {np.percentile(sr, 16):.3f}-{np.percentile(sr, 84):.3f}); "
          f"{int((np.array(ref_n) < MIN_REF_WAVEFORMS).sum())} module-files on the median sigma_ref; {int(n_base_ref.sum())} baseline starts from b_ref", flush=True)
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    np.savez_compressed(out, run=run, run_type=rt, splits=np.array([int(p.group(2)) for p in parsed]), n_seen=n_seen,
                        n_trig=n_trig, n_lms=n_lms, n_pass=len(evnum), wall_s=wall, beam_energy_mev=cfg["beam"],
                        mod_n_ref=np.array(ref_n), mod_b_ref=np.array(ref_b), mod_sig_ref=np.array(ref_sig), mod_n_base_ref=n_base_ref, dead_modules=np.array(sorted(dead), dtype=str), mod_name=names, mod_pwo4=pwo4,
                        mod_n_bad=n_bad, mod_n_trunc=n_trunc, mod_t50=t50[bank.row].astype(np.float32), ev_number=np.array(evnum, dtype=np.int64),
                        **{f"ch_{q}": v for q, v in r.items()}, **({f"ex_{q}": v for q, v in ex.items()} if ex else {}))
    print(f"[done] {out}: {n_seen} events, {len(evnum)} after trigger gate and light-monitoring veto, "
          f"{r['k'].size} channels kept; {wall:.0f} s", flush=True)
    top = lambda c: ", ".join(f"{names[i]} {c[i]}" for i in np.argsort(-c)[:3])
    print(f"[bad waveforms] {n_bad.sum()} (most: {top(n_bad)}); [truncated main pulses] {n_trunc.sum()} (most: {top(n_trunc)})", flush=True)


# report =======================================================================
def quantiles(h, qs):
    """Quantiles of a histogram on CHI2_FINE (log interpolation inside a bin); NaN if empty."""
    n = h.sum()
    if n == 0:
        return [np.nan] * len(qs)
    c = np.r_[0.0, np.cumsum(h)] / n
    lo, hi = np.log(np.r_[1e-3, CHI2_FINE[1:-1]]), np.log(np.r_[CHI2_FINE[1:-1], 1e5])   # open end bins: 1e-3 and 1e5
    out = []
    for q in qs:
        i = min(np.searchsorted(c, q, side="left") - 1, h.size - 1)
        i = max(i, 0)
        f = (q - c[i]) / max(c[i + 1] - c[i], 1e-300)
        out.append(float(np.exp(lo[i] + f * (hi[i] - lo[i]))))
    return out


def hist(x):
    return np.histogram(x[np.isfinite(x)], CHI2_FINE)[0]


def load(scan_dir):
    """Per run type: fine chi-squared histograms of seeds and neighbours per slice, the seed-neighbour 2D histogram, the
    (t_cfd - t_mf, chi-squared) 2D histograms of seeds and neighbours, counts, and the example reservoirs."""
    files = sorted(Path(scan_dir).glob("run_*.npz"))
    if not files:
        sys.exit(f"[ERROR] no run_*.npz in {scan_dir}")
    ns, nb, nf = len(E_SLICES) - 1, len(DT_EDGES) - 1, CHI2_FINE.size - 1
    lb = np.log10(np.r_[0.0, np.logspace(np.log10(CHI2_LO), np.log10(CHI2_HI), 61), np.inf] + 1e-300)   # 2D seed-neighbour bins
    D, nbr, nl = {}, None, None
    for fn in files:
        z = np.load(fn, allow_pickle=False)
        rt, run = str(z["run_type"]), int(z["run"])
        d = D.setdefault(rt, dict(runs=[], splits=0, wall=0.0, n_pass=0, n_cl=0, n_nogain=0, n_out=0,
                                  seed=np.zeros((ns, nf)), nbr=np.zeros((ns, nf)), sn=np.zeros((lb.size - 1,) * 2),
                                  prof=np.zeros((2, nb, nf)), n_cfd=np.zeros(2), n_all=np.zeros(2), ex=[]))
        d["runs"].append(run)
        d["splits"] += len(z["splits"])
        d["wall"] += float(z["wall_s"])
        d["n_pass"] += int(z["n_pass"])
        names = z["mod_name"]
        if nbr is None:
            nbr = neighbours(names)
            nl = neighbour_lists(nbr)
        gain = module_gain(run_config(run), rt, names)
        ev, mi, k, a, t, chi2 = (z[f"ch_{q}"] for q in ("ev", "midx", "k", "a", "t", "chi2"))
        seed, energy, mem, cl = find_clusters(ev, mi, k, a.astype(np.float64), t, nbr, nl, gain)
        sl = np.where(np.isfinite(energy), np.digitize(np.nan_to_num(energy), E_SLICES) - 1, -1)
        d["n_cl"] += seed.size
        d["n_nogain"] += int((~np.isfinite(energy)).sum())
        d["n_out"] += int((np.isfinite(energy) & ((sl < 0) | (sl >= ns))).sum())
        nb_mem = mem[mem != seed[cl]]                          # neighbours: non-seed members
        nb_sl = sl[cl[mem != seed[cl]]]
        for i in range(ns):
            d["seed"][i] += hist(chi2[seed[sl == i]])
            d["nbr"][i] += hist(chi2[nb_mem[nb_sl == i]])
        ins = (nb_sl >= 0) & (nb_sl < ns)
        xs, ys = chi2[seed[cl[mem != seed[cl]]]][ins], chi2[nb_mem][ins]
        ok = np.isfinite(xs) & np.isfinite(ys)
        d["sn"] += np.histogram2d(np.log10(xs[ok]), np.log10(ys[ok]), [lb, lb])[0]
        main = a.argmax(1)
        dt = z["ch_tcfd"] - (t[np.arange(k.size), main] + z["mod_t50"][mi])   # t_cfd - t_mf (matched-filter 50 % time)
        for j, ch in enumerate((seed[(sl >= 0) & (sl < ns)], nb_mem[ins])):
            ok = np.isfinite(dt[ch]) & np.isfinite(chi2[ch])
            d["n_all"][j] += ch.size
            d["n_cfd"][j] += ok.sum()
            d["prof"][j] += np.histogram2d(dt[ch][ok], chi2[ch][ok], [DT_EDGES, CHI2_FINE])[0]
        if "ex_key" in z.files:
            d["ex"].append({q[3:]: z[q] for q in z.files if q.startswith("ex_")} | {"run": np.full(z["ex_key"].size, run)})
    return D


def slice_label(i):
    return f"{E_SLICES[i]}-{E_SLICES[i + 1]} MeV"


def report(scan_dir, out):
    import matplotlib.pyplot as plt
    # plot style: large bold axis labels, large legends
    plt.rcParams.update({
        "axes.labelsize": 16, "axes.labelweight": "bold", "axes.titlesize": 16, "axes.titleweight": "bold",
        "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 14, "axes.linewidth": 1.8,
        "xtick.major.width": 1.5, "ytick.major.width": 1.5, "xtick.major.size": 6, "ytick.major.size": 6,
        "lines.linewidth": 2.0, "figure.titlesize": 18, "figure.titleweight": "bold"})
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import LogNorm
    D = load(scan_dir)
    os.makedirs(out, exist_ok=True)
    ns = len(E_SLICES) - 1
    h2 = [rt for rt in H2_LABELS if rt in D]
    plot_edges = np.logspace(np.log10(CHI2_LO), np.log10(CHI2_HI), 61)
    fine_c = np.sqrt(CHI2_FINE[1:-2] * CHI2_FINE[2:-1])      # centres of the finite fine bins
    x_label = "chi-squared per degree of freedom"
    lines = ["run_type\tmodules\tslice_MeV\tentries\tmedian\tp90"]
    summary = []
    stats = {}
    for rt, d in D.items():
        summary.append(f"{rt}: runs {sorted(d['runs'])}, {d['splits']} splits, scan wall {d['wall']:.0f} s, {d['n_pass']} events after "
                       f"gate and veto, {d['n_cl']} clusters, {d['n_nogain']} with a seed without calibration factor, "
                       f"{d['n_out']} outside 10-4000 MeV; constant-fraction-discriminator time for "
                       f"{d['n_cfd'][0]:.0f} of {d['n_all'][0]:.0f} seeds and {d['n_cfd'][1]:.0f} of {d['n_all'][1]:.0f} neighbours")
        for role in ("seed", "nbr"):
            for i in range(ns):
                h = d[role][i]
                med, p90 = quantiles(h, [0.5, 0.9])
                stats[rt, role, i] = (int(h.sum()), med, p90)
                lines.append(f"{rt}\t{'seed' if role == 'seed' else 'neighbour'}\t{E_SLICES[i]}-{E_SLICES[i + 1]}\t{int(h.sum())}\t{med:.4g}\t{p90:.4g}")
    with open(Path(out) / "chi2_slices.tsv", "w") as fh:
        fh.write("\n".join(lines) + "\n")

    # (a) the requested figure: seed chi-squared per degree of freedom, one row per H2 run type, one panel per slice
    fig, axs = plt.subplots(len(h2), ns, figsize=(4.4 * ns, 4.2 * len(h2) + 1.5), squeeze=False, sharex=True)
    for r, rt in enumerate(h2):
        for i in range(ns):
            ax = axs[r, i]
            h = D[rt]["seed"][i]
            cnt = np.histogram(fine_c, plot_edges, weights=h[1:-1])[0]
            ax.stairs(cnt, plot_edges, fill=True, alpha=0.6)
            n, med, p90 = stats[rt, "seed", i]
            for v, ls in ((med, "-"), (p90, "--")):
                if np.isfinite(v):
                    ax.axvline(v, color="k", ls=ls, lw=0.8)
            ax.set_xscale("log")
            ax.set_xlim(CHI2_LO, CHI2_HI)
            ax.text(0.97, 0.95, f"entries {n}\nmedian {med:.2f}\n90th percentile {p90:.2f}", transform=ax.transAxes,
                    ha="right", va="top", fontsize=14)
            if r == 0:
                ax.set_title(slice_label(i))
            if i == 0:
                ax.set_ylabel(f"{H2_LABELS[rt]}\nseeds per bin")
            if r == len(h2) - 1:
                ax.set_xlabel("chi-squared per dof")
    fig.suptitle("Seed module: matched-filter chi-squared per degree of freedom in slices of cluster energy\n"
                 "(solid: median, dashed: 90th percentile)")
    fig.tight_layout()
    fig.savefig(Path(out) / "seed_chi2_vs_cluster_energy.png", dpi=120)
    plt.close(fig)

    # companion: median and 90th percentile vs slice, seeds and neighbours
    fig, axs = plt.subplots(1, 2, figsize=(16, 7), sharey=True)
    x = np.arange(ns)
    for ax, role, title in ((axs[0], "seed", "seed module"), (axs[1], "nbr", "neighbour modules (non-seed members)")):
        for c, rt in enumerate(h2):
            med = [stats[rt, role, i][1] if stats[rt, role, i][0] >= MIN_ENTRIES else np.nan for i in range(ns)]
            p90 = [stats[rt, role, i][2] if stats[rt, role, i][0] >= MIN_ENTRIES else np.nan for i in range(ns)]
            ax.plot(x, med, "o-", color=f"C{c}", label=f"{H2_LABELS[rt]} median")
            ax.plot(x, p90, "s--", color=f"C{c}", label=f"{H2_LABELS[rt]} 90th percentile")
        ax.set_xticks(x, [f"{E_SLICES[i]}-{E_SLICES[i + 1]}" for i in x], rotation=30)
        ax.set_xlabel("cluster energy (MeV)")
        ax.set_yscale("log")
        ax.set_title(title)
        ax.grid(alpha=0.3)
    axs[0].set_ylabel(x_label)
    axs[1].legend()
    fig.suptitle("Matched-filter chi-squared per degree of freedom vs cluster energy")
    fig.tight_layout()
    fig.savefig(Path(out) / "chi2_summary_vs_cluster_energy.png", dpi=120)
    plt.close(fig)

    # (c1) chi-squared vs t_cfd - t_mf: median and 16-84 % band per 0.5 ns bin, seeds and neighbours
    fig, axs = plt.subplots(1, len(h2), figsize=(8 * len(h2), 7), squeeze=False, sharey=True)
    ctr = 0.5 * (DT_EDGES[1:] + DT_EDGES[:-1])
    prof_lines = ["run_type\tmodules\tdt_centre_ns\tentries\tmedian\tp16\tp84"]
    for ax, rt in zip(axs[0], h2):
        for j, (role, col) in enumerate((("seed", "C3"), ("neighbour", "C0"))):
            P = D[rt]["prof"][j]
            q = np.array([quantiles(P[b], [0.5, 0.16, 0.84]) if P[b].sum() >= MIN_ENTRIES else [np.nan] * 3 for b in range(ctr.size)])
            ax.plot(ctr, q[:, 0], color=col, label=f"{role} median")
            ax.fill_between(ctr, q[:, 1], q[:, 2], color=col, alpha=0.25, label=f"{role} 16-84 %")
            prof_lines += [f"{rt}\t{role}\t{c:.2f}\t{int(P[b].sum())}\t{q[b, 0]:.4g}\t{q[b, 1]:.4g}\t{q[b, 2]:.4g}" for b, c in enumerate(ctr)]
        ax.set_yscale("log")
        ax.set_title(H2_LABELS[rt])
        ax.set_xlabel("CFD time - matched-filter time (ns)")
        ax.grid(alpha=0.3)
    axs[0, 0].set_ylabel(x_label)
    axs[0, 0].legend()
    fig.suptitle("Matched-filter chi-squared per degree of freedom vs timing difference\n(cluster members, 10-4000 MeV clusters)")
    fig.tight_layout()
    fig.savefig(Path(out) / "chi2_vs_dt_cfd_mf.png", dpi=120)
    plt.close(fig)
    with open(Path(out) / "chi2_vs_dt_cfd_mf.tsv", "w") as fh:
        fh.write("\n".join(prof_lines) + "\n")

    # (c2) seed vs neighbour chi-squared
    lb = np.r_[np.logspace(np.log10(CHI2_LO), np.log10(CHI2_HI), 61)]
    fig, axs = plt.subplots(1, len(h2), figsize=(8 * len(h2), 7), squeeze=False)
    for ax, rt in zip(axs[0], h2):
        H = D[rt]["sn"][1:-1, 1:-1]
        m = ax.pcolormesh(lb, lb, H.T, norm=LogNorm(vmin=1, vmax=max(H.max(), 1)), cmap="viridis")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(f"{H2_LABELS[rt]} ({int(D[rt]['sn'].sum())} neighbours)")
        ax.set_xlabel(f"seed {x_label}")
        ax.set_ylabel(f"neighbour {x_label}")
        cb = fig.colorbar(m, ax=ax)
        cb.set_label("neighbours per bin", fontsize=16, fontweight="bold")
        cb.ax.tick_params(labelsize=14)
    fig.suptitle("Seed vs neighbour matched-filter chi-squared per degree of freedom\n(10-4000 MeV clusters)")
    fig.tight_layout()
    fig.savefig(Path(out) / "seed_vs_neighbour_chi2.png", dpi=120)
    plt.close(fig)

    # (c3) examples: seeds with chi-squared per degree of freedom > EX_CHI2, EX_PER_SLICE per slice, run type EX_RUN_TYPE
    pdf_path = Path(out) / f"examples_seed_chi2_gt{EX_CHI2:g}_{EX_RUN_TYPE}.pdf"
    exs = D.get(EX_RUN_TYPE, {}).get("ex", [])
    E = {q: np.concatenate([e[q] for e in exs]) for q in exs[0]} if exs else None
    names = np.load(sorted(Path(scan_dir).glob("run_*.npz"))[0])["mod_name"]
    rng = np.random.default_rng(0)
    with PdfPages(pdf_path) as pdf:
        for i in range(ns):
            fig, axs = plt.subplots(2, 3, figsize=(24, 14))
            idx = np.nonzero(E["slice"] == i)[0] if E is not None else np.zeros(0, int)
            idx = rng.permutation(idx)[:EX_PER_SLICE]
            for ax, j in zip(axs.flat, idx):
                ts = DT * np.arange(E["samples"].shape[1])
                ax.plot(ts, E["samples"][j], ".", color="k", ms=3, label="raw samples")
                ax.plot(ts, E["model"][j], "-", color="C3", lw=1, label="pedestal + matched-filter fit")
                for tt in E["t"][j][np.isfinite(E["t"][j])]:
                    ax.axvline(tt, color="C0", lw=0.8, ls="--")
                if np.isfinite(E["tcfd"][j]):
                    ax.axvline(E["tcfd"][j], color="C2", lw=0.8, label="constant-fraction-discriminator time")
                ax.axvspan(PRE_LO, PRE_HI, color="0.9", zorder=0)
                ax.set_title(f"run {E['run'][j]} event {E['evnum'][j]} {names[E['midx'][j]]}\ncluster {E['energy'][j]:.0f} MeV, "
                             f"A {E['a'][j].max():.0f} ADC, {E['k'][j]} pulse(s)\nchi-squared per dof {E['chi2'][j]:.1f}",
                             fontsize=12)
                ax.set_xlabel("time (ns)")
                ax.set_ylabel("ADC")
            for ax in list(axs.flat)[len(idx):]:
                ax.axis("off")
            if len(idx):
                axs.flat[0].legend(loc="best")
            else:
                axs.flat[0].text(0.5, 0.5, f"no seed with chi-squared per degree of freedom > {EX_CHI2:g}", ha="center")
            fig.suptitle(f"{EX_RUN_TYPE}: seeds with matched-filter chi-squared per degree of freedom > {EX_CHI2:g},\n"
                         f"cluster {slice_label(i)} (dashed: matched-filter onsets; grey: in-window onset range)")
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
    with open(Path(out) / "summary.txt", "w") as fh:
        fh.write("\n".join(summary) + "\n")
    print("\n".join(summary))
    print("\n".join(lines))
    print(f"[report] wrote {out}/: seed_chi2_vs_cluster_energy.png chi2_summary_vs_cluster_energy.png chi2_slices.tsv "
          f"chi2_vs_dt_cfd_mf.png chi2_vs_dt_cfd_mf.tsv seed_vs_neighbour_chi2.png {pdf_path.name} summary.txt")


def main():
    a = sys.argv[1:]
    if not DB_DIR.name:
        sys.exit("[ERROR] PRAD2_DATABASE_DIR is not set (source setup_env.sh)")
    if len(a) >= 3 and a[0] == "scan":
        scan(a[1], a[2:])
    elif len(a) == 3 and a[0] == "report":
        report(a[1], a[2])
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main()
