#!/usr/bin/env python3
"""energy_resolution_benchmark.py -- PbWO4 energy resolution: legacy versus time-gated island clustering.

    python3 waveform_analysis/energy_resolution_benchmark.py scan <out.npz> <evio file> [<evio file> ...]
    python3 waveform_analysis/energy_resolution_benchmark.py report [--out <dir>] <scan dir> <gate_ns> [<gate_ns> ...]
    python3 waveform_analysis/energy_resolution_benchmark.py summary <out dir>     # <out dir>/summary.pdf from the five reports

Algorithm
  scan (one run, no clustering): physics events -> every PbWO4 channel through the matched filter MF of
  pileup_analysis/pileup_deconvolution_study.py (copied) with the run type's template bank: hit time = onset + t50, energy =
  calibration x amplitude x template integral; lead glass and clipped channels without a pulse use WaveAnalyzer.  Hits with
  PRE_LO < t < PRE_HI are stored.
  report (per run type and split): HyCalCluster on the hits, legacy (one hit per module, no time coincidence) and one path per
  time gate; events with one cluster above SIGNAL_MIN; pool = seed cuts; module set = elastic-window method of
  benchmark_hycal_timing_test_20261021.py.  Corrections come from the legacy pool: per-crystal gains (percrystal_calib.py),
  position correction (position_corr.py), kinematic Moller rejection (--moller-reject kin) applied to the event selection and
  the gain fit; the clusterer's non-linearity factor is switched on or off.  Arms: none, gains + position, + Moller, +
  non-linearity, all on.  Fits of the used-module spectrum: two-peak Crystal Ball (headline), Gaussian within 2 sigma, interquartile
  range; default arm for legacy and every gate, the other arms for legacy and gate ARM_GATE.
  Moller rejection and non-linearity default off (switches).  OUT_OF_SAMPLE: the default arm is also
  evaluated (legacy and every gate) with gains and position correction fitted on the other half of the splits (even/odd ordinal within a
  run), columns *_oos.  Position-correction setting per run type: PC_SETTING.
Standalone: numpy, scipy, matplotlib, prad2py; database from $PRAD2_DATABASE_DIR (source setup_env.sh).  Recipe and outputs:
waveform_analysis/energy_resolution/README.md.
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")                     # one BLAS thread per process

import csv
import hashlib
import json
import math
import multiprocessing as mp
import re
import sys
import time
import warnings
from collections import Counter
from pathlib import Path

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from prad2py import dec, det
from scipy.interpolate import BSpline
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

# --- paths and runs (run lists as run_energy_resolution.sh)
HERE = Path(__file__).resolve().parent
PYDIR = HERE.parent                                     # analysis/pyscripts
BANK_ROOT = HERE / "pulse_templates"                    # <run type>/templates_<run type>.npz
OUT_ROOT = HERE / "energy_resolution"                   # report output: <run type>/
DB_DIR = Path(os.environ.get("PRAD2_DATABASE_DIR", ""))
RUNINFO = DB_DIR / "runinfo" / "general.json"           # per run: beam energy, calibration file, dead modules
RUN_TYPE_OF = {
    24590: "0p73GeV_H2", 24603: "0p73GeV_H2",
    25228: "2p24GeV_H2", 25257: "2p24GeV_H2", 25260: "2p24GeV_H2", 25299: "2p24GeV_H2", 25320: "2p24GeV_H2",
    24698: "3p49GeV_H2", 24911: "3p49GeV_H2", 24937: "3p49GeV_H2",
    24386: "3p49GeV_H2_old", 24503: "3p49GeV_H2_old", 24512: "3p49GeV_H2_old",   # previous 3.49 GeV runs; bank of 3p49GeV_H2
    25800: "2p24GeV_Ta_X17", 26138: "2p24GeV_Ta_X17", 26230: "2p24GeV_Ta_X17",
}
M_PROTON, M_ELECTRON, U_AMU = 938.272, 0.51099895, 931.494   # MeV
TARGET_MASS = {"H2": M_PROTON, "Ta": 180.948 * U_AMU}        # elastic target (Ta for the X17 runs)
EBEAM_TOL = 1e-3                       # builder beam energy versus runinfo, relative

# --- scan: events, hits, noise reference
PHYSICS_BITS = 0x100                   # physics events: trigger bits == this (light-monitoring events fail it)
PHYSICS_MASK_ANY = 0xf00               # X17 runs: trigger bits & this != 0 (as template_extract.py) ...
PHYSICS_ANY = {"2p24GeV_Ta_X17"}       # ... because the exact bit starves these runs
LMS_BIT = 0x1000000                    # light-monitoring trigger bit, vetoed on both gates
NS_EXPECT, OVERFLOW, BLOCK = 100, 4095, 256   # samples per PbWO4 channel, ADC overflow, events per matched-filter block
PRE_LO, PRE_HI = 100.0, 200.0          # ns: a hit needs PRE_LO < t < PRE_HI (no time offsets)
N_PHASE = 256                          # sampling phases of the walk-out integral table (1.8 % peak to peak per module)
N_REF_EVENTS = 20000                   # first decoded events of an EVIO file that set its per-module noise reference
MIN_REF_WAVEFORMS = 20                 # modules with fewer reference waveforms take the median reference sigma
START_NSIGMA = 3.0                     # baseline farther than this many reference sigma from the reference starts the fit at it

# --- matched filter (copied from pileup_deconvolution_study.py, same values; comments there)
PED_NSAMPLES = 30       # baseline samples
UNDERSHOOT_ADC = 20.0   # a sample this far below the baseline, or any sample equal to 0, marks a bad waveform
KMAX = 6                # at most this many pulses per waveform
NSIGMA, MIN_AMP = 5.0, 5.0   # a pulse needs A >= max(NSIGMA sigma_A, MIN_AMP)
MIN_SEP_NS = 8.0        # smallest onset separation of two pulses
MERGE_NS = 6.5          # post-fit merge below this separation
R_MODEL = 0.02          # stage C: a shoulder needs A >= R_MODEL x every pulse within one kernel span (96 ns)
PEAK_SMOOTH_NS = 8.0    # stage A: Gaussian sigma (ns) of the amplitude smoothing before the peak search
PEAK_FRAC = 0.02        # stage A: a candidate's smoothed height >= PEAK_FRAC x the waveform's highest candidate
SHAPE_SIGMA_MIN = 0.005 # floor of shape_sigma (fraction of the peak) in the sigma(t) variance
A_MIN_SCALE = 20.0      # ADC; floor of the sqrt(a_band / A) scaling of shape_sigma
BASE_CLIP = 1.0         # free-baseline offset limited to +-1 sigma
RAW_NSIGMA = 5.0        # an accepted pulse needs raw max - baseline > RAW_NSIGMA sigma_ped within its own window
LATE_MIN_TAIL = 6       # a pulse near the window end needs this many samples after its peak inside the window
BIC_PENALTY = 1.0       # stage C: chi2 improvement must exceed BIC_PENALTY x 2 ln N_w
DT, NPH, L, LW, LEAD = 4.0, 4, 24, 32, 3              # ns per sample, phases, kernel / fit window samples, lead samples
LUT_DT, SIGMA_MIN, CAND_NSIGMA = 0.05, 0.3, 3.0        # template grid (ns), noise floor (ADC), candidate threshold (sigma)
F32 = np.float32

# --- clustering
CLUSTER_FIELDS = ("min_module_energy", "min_center_energy", "min_cluster_energy", "min_cluster_size",
                  "corner_conn", "split_iter", "least_split", "log_weight_thres")   # copied from the builder configuration
MIN_CENTER_E = 30.0                    # MeV, seed threshold (database default 10)
LEGACY_WINDOW = -1.0                   # seed_time_window of the legacy path (no time coincidence)

# --- pool selection (MeV at E_REF; times kE at other beam energies)
E_REF = 2239.51                        # reference beam energy (MeV)
SIGNAL_MIN = 500.0                     # events with exactly one cluster above this
SEED_MIN = 200.0                       # pool: seed module energy >= this
INNER_R = 200.0                        # pool: seed module centre radius < this (mm)
EXCLUDE_ROWCOL = (15, 20, 15, 20)      # pool: seed not in rows 15-20 and columns 15-20 (1-based, inclusive)

# --- module set (elastic-window method, no gain applied)
MIN_PER_MODULE = 30                    # events within mu_M +- WINDOW_NSIGMA sigma_M
WINDOW_NSIGMA = 2.0                    # elastic window half width, sigma_M
PREFIT_RANGE = (-0.15, 0.10)           # per-module Crystal Ball prefit range, fractions of mu_global (divided by kS)
MIN_FIT_EVENTS = 100                   # fewer events in the prefit range: median and the pool sigma instead of a fit
MOLLER_GAP_NSIGMA = 3.0                # module unused if its elastic-Moller gap < this x sigma_M

# --- fits
FIT_HALF = 800.0                       # spectrum range E' +- FIT_HALF kE
BIN = 20.0                             # bin width, x kS
FIX_HALF = 150.0                       # fixed-window Gaussian: +- FIX_HALF kS around the argmax
G2_NSIGMA, G2_MAX_ITER, G2_TOL = 2.0, 5, 0.01   # iterated Gaussian within 2 sigma
MOLLER_RES_A, MOLLER_RES_B = 0.025, 0.007        # Moller resolution sigma_E/E = sqrt((A / sqrt(E / GeV))^2 + B^2)
MOLLER_TIE_NSIGMA = 2.0                # two-peak fit: mu_M within E_ee(r98 .. r2) +- this x sigma_E(mu_M)
MOLLER_TIE_PCT = (2.0, 98.0)           # cluster-radius percentiles of the tie range
YIELD_NSIGMA = 2.5                     # peak yield: |E / E' - 1| <= this x s_ref
TAIL_NSIGMA = 3.0                      # tail fraction outside mu_CB +- this x sigma_CB
NPROC = 8                              # report: worker processes (one split each)

# --- gains (percrystal_calib.py solve, windows refitted in the first passes)
GAIN_SEED_FRAC = 0.64                  # mean seed energy fraction; sets the prior strength
GAIN_FREE_FRAC = 0.5                   # a module with a prior fraction below this is free (convergence test)
GAIN_DG_TOL, GAIN_TOL_SIGMA = 1e-3, 0.5   # converged when |dg| < max(GAIN_DG_TOL, GAIN_TOL_SIGMA sigma_g) for every free module
GAIN_MAX_ITER, GAIN_N_REFIT = 6, 3     # passes; the first GAIN_N_REFIT refit the module windows, later passes keep them

# --- position correction (position_corr.py map / fit --per-module; notes/poscorr_form_choice_20261106.md)
PC_SETTING = {"3p49GeV_H2": ("hybrid", 1000), "3p49GeV_H2_old": ("hybrid", 1000), "2p24GeV_H2": ("hybrid", 1000),
              "0p73GeV_H2": ("none", 0), "2p24GeV_Ta_X17": ("2p24GeV_H2", 0)}   # (form, own-polynomial threshold in clusters): hybrid = own polynomials + map, none = gains only, a run type = its hybrid correction borrowed (run that report first); best out-of-sample, poscorr_comparison.tsv
PC_FIT = {"3p49GeV_H2": "cb", "3p49GeV_H2_old": "cb", "2p24GeV_H2": "gauss", "0p73GeV_H2": "gauss"}   # map bin peak fitter
PC_R = (0.85, 1.12)                    # training window of r = E / E'
PC_BINS, PC_RANGE = 12, 0.6            # bilinear map: PC_BINS x PC_BINS bins over u, v in [-PC_RANGE, PC_RANGE]
PC_MIN_FIT, PC_MIN_N = 30, 300         # map bin: fitted from PC_MIN_FIT events, used from PC_MIN_N (others filled from neighbours)
PC_CLAMP = 0.6                         # hybrid: u, v clipped to +- this
PC_GRID, PC_GRID_MIN, PC_GRID_BINS = 6, 30, 10   # hybrid: grid bins per axis, events per grid bin, grid bins needed
PC_TRUNC_K = 2.0                       # hybrid: +- this sigma truncated mean per grid bin
PC_ACCEPT = (0.9, 1.1)                 # hybrid: own fit kept only if R stays in this range over the clamp square

# --- Moller rejection and non-linearity (old chain: --moller-reject kin, --nl-members on)
MOLLER_REJECT = False                  # kinematic Moller rejection in the arms that have it: event selection and gain fit (default off pending user confirmation, 2026-11-12)
MOLLER_KIN_N, MOLLER_KIN_M = 3.0, 2.0  # Moller if |E - E_ee| < N sigma and |E - E'| > M sigma (cluster polar angle)
NONLINEAR = False                      # calibration non-linearity factor in the arms that have it (clusterer's factor, per cluster; default off pending user confirmation, 2026-11-12)

# --- report arms: (gains + position correction, Moller rejection, non-linearity)
ARMS = {"none": (0, 0, 0), "gp": (1, 0, 0), "gpm": (1, 1, 0), "gpn": (1, 0, 1), "all": (1, 1, 1)}
ARM_LABEL = {"none": "none", "gp": "gains + position", "gpm": "+ Moller rejection", "gpn": "+ non-linearity", "all": "all on"}
DEFAULT_ARM = "all"                    # full gate scan (summary.tsv); the other arms: legacy and gate ARM_GATE
ARM_GATE = 12.0                        # ns
OUT_OF_SAMPLE = True                   # the default report also quotes sigma/E with gains and position correction fitted on the other half of the splits (even/odd cross fit)
OTHER_ARMS = [a for a in ARMS if a != DEFAULT_ARM]   # arms fitted for legacy and gate ARM_GATE only
CORR_Q = ("res_cb", "res_g2", "mu_ep")   # quantities kept per arm
OOS_Q = tuple(f"{q}_oos" for q in CORR_Q) if OUT_OF_SAMPLE else ()   # out-of-sample columns (default arm, legacy and every gate)
RATIO_EDGES = np.arange(0.6, 1.4 + 1e-9, 0.004)   # bins of E / E' stored in spectra.npz (summary pages)
NAN2 = (float("nan"),) * 2             # failed (mu, sigma) fit


class Bank:
    """Spline bank of one run type on the 0.05 ns grid, linear interpolation.  LUT row 0 is the bank's average spline (lead
    glass, PbWO4 without a valid row or dead); row[midx] maps modules to rows.  slut: the row's shape_sigma on the same grid
    (linear in the bank's t_grid, ends held), at least SHAPE_SIGMA_MIN; aband: the row's a_band (row 0: median over the
    bank's valid, not dead rows)."""


    def __init__(self, path, names, pwo4):
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


# scan =========================================================================
def bad_waveform(Y):
    """[n] bool: any sample == 0, or any sample < median(first PED_NSAMPLES) - UNDERSHOOT_ADC."""
    base = np.median(Y[:, :PED_NSAMPLES], axis=1)
    return (Y == 0).any(axis=1) | (Y < base[:, None] - UNDERSHOOT_ADC).any(axis=1)


def pedestal(Y):
    """Clipped mean of the first PED_NSAMPLES samples of each waveform: samples within 2 sigma (at least 1 ADC) of the median."""
    x = Y[:, :PED_NSAMPLES].astype(F32)
    med = np.median(x, axis=1)
    sig = np.maximum(1.4826 * np.median(np.abs(x - med[:, None]), axis=1), SIGMA_MIN)
    ok = np.abs(x - med[:, None]) <= np.maximum(2.0 * sig, 1.0)[:, None]
    return ((x * ok).sum(axis=1) / np.maximum(ok.sum(axis=1), 1)).astype(F32)


def noise_reference(X, mi, n_mod):
    """Per-module (b_ref, sigma_ref) from the pedestal samples X [n, PED_NSAMPLES] of waveforms of modules mi: median baseline and
    1.4826 median absolute deviation of the pooled residuals.  Modules with < MIN_REF_WAVEFORMS waveforms: b_ref NaN (no
    override) and the median sigma_ref."""
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
    return b_ref, sig_ref


def run_ref(mf, Y, midx, sig, b_ref):
    """MF.run with the noise sigma given per waveform and the baseline fit started at b_ref where the 30-sample baseline is farther
    than START_NSIGMA sigma from it.  Returns t0, a, valid, truncated [n, KMAX] and far [n] (start replaced by b_ref)."""
    Yf = Y.astype(F32)
    n, ns = Yf.shape
    ped = pedestal(Yf)
    far = np.abs(ped - b_ref) > START_NSIGMA * sig
    y = Yf - np.where(far, b_ref, ped).astype(F32)[:, None]
    clipped = (Yf >= OVERFLOW).any(axis=1)
    cand = ((y[:, :-2] + y[:, 1:-1] + y[:, 2:]) * (1.0 / 3.0)).max(axis=1) > np.maximum(CAND_NSIGMA * sig, 3.0)
    W = (Yf < OVERFLOW).astype(F32)
    T0, A, V = np.zeros((n, KMAX), F32), np.zeros((n, KMAX), F32), np.zeros((n, KMAX), bool)
    for masked in (False, True):
        sel = np.nonzero(cand & (clipped == masked))[0]
        if sel.size == 0:
            continue
        ys, ms, w = y[sel], midx[sel], W[sel] if masked else None
        if masked:
            ys = ys * w
        t0, a, v, _ = mf.find_pulses(ys, w, ms, sig[sel])
        if masked:
            v &= ~mf.bad_clipped(v, t0, w, ms, ns)[:, None]
        T0[sel], A[sel], V[sel] = t0, np.where(v, a, 0.0), v
    trunc = V & (np.ceil(T0 / DT).astype(np.int64) - LEAD + LW > ns)
    return T0.astype(np.float64), A.astype(np.float64), V, trunc, far


def run_config(run):
    """Beam energy (MeV), calibration file and dead modules of a run from runinfo/general.json (as template_extract's resolve_calib):
    "configurations" entries with from_run <= run override one another; dead modules from the latest "hycal_dead_modules" entry."""
    j = json.load(open(RUNINFO))
    beam, calib = None, ""
    for c in sorted((c for c in j["configurations"] if "from_run" in c), key=lambda c: c["from_run"]):
        if c["from_run"] <= run:
            beam, calib = c.get("beam_energy", beam), c.get("calibration", {}).get("file", calib)
    if beam is None or not calib:
        sys.exit(f"[ERROR] runinfo has no beam energy or calibration file for run {run}")
    d = max((e for e in j["hycal_dead_modules"] if e["from_run"] <= run), key=lambda e: e["from_run"])
    return dict(beam=float(beam), calib=calib, dead=list(d["dead_modules"]))


def build(run, calib):
    """prad2py pipeline of the run with the given calibration file."""
    b = det.PipelineBuilder()
    b.set_hycal_calib(calib)
    b.set_run_number(run)
    return b.build()


def template_columns(path, names):
    """t50 and iint_wa of the bank per module (median of the valid rows for modules without one)."""
    z = np.load(path, allow_pickle=False)
    ok = z["valid"] > 0
    row = {str(s): i for i, s in enumerate(z["names"]) if ok[i]}
    out = []
    for key in ("t50", "iint_wa"):
        v, med = z[key], float(np.median(z[key][ok]))
        out.append(np.array([v[row[n]] if n in row else med for n in names], dtype=np.float64))
    return out


def smooth_tri(x, order):
    """WaveAnalyzer triangular smoothing along the last axis (edge samples unsmoothed)."""
    n = x.shape[-1]
    val, wsum = x.copy(), np.ones(n)
    for j in range(1, order):
        w = 1.0 - j / (order + 1.0)
        val[..., j:n - j] += w * (x[..., :n - 2 * j] + x[..., 2 * j:])
        wsum[j:n - j] += 2.0 * w
    return val / wsum


def walk_integral(v, ratio, nbreak):
    """WaveAnalyzer peak integral of the smoothed rows v[R, N]: walk outward from the maximum; nbreak consecutive samples below
    ratio x peak (or negative) end the walk and are dropped."""
    R, N = v.shape
    ar = np.arange(R)
    pk = np.argmax(v, axis=1)
    h = v[ar, pk]
    cut = ratio * h
    integral = h.copy()
    for side in (-1, 1):
        run, pend, active = np.zeros(R, dtype=np.int64), np.zeros(R), np.ones(R, dtype=bool)
        for d in range(1, N):
            pos = pk + side * d
            inr = active & (pos >= 0) & (pos < N)
            if not inr.any():
                break
            vv = v[ar, np.clip(pos, 0, N - 1)]
            below = (vv < cut) | (vv < 0.0)
            mb = inr & below
            run = np.where(mb, run + 1, run)
            pend = np.where(mb, pend + vv, pend)
            active &= ~(mb & (run >= nbreak))
            ma = inr & ~below
            integral += np.where(ma, pend + vv, 0.0)
            pend = np.where(ma, 0.0, pend)
            run = np.where(ma, 0, run)
    return integral


def phase_table(bank):
    """Per bank row the walk-out integral (iint_wa of template_extract.py, not phase-averaged) at N_PHASE sampling phases / its mean."""
    tau = DT * (np.arange(64) - 12)[None, :] + (DT / N_PHASE) * np.arange(N_PHASE)[:, None]
    tab = np.array([walk_integral(smooth_tri(bank.shape(np.full(tau.shape, r), tau, np.float64), 2), 0.1, 2)
                    for r in range(len(bank.lut) // bank.n)])
    tab /= tab.mean(axis=1, keepdims=True)
    return np.hstack([tab, tab[:, :1]])


def scan(out, evios):
    t_start = time.monotonic()
    parsed = [re.search(r"prad_0*(\d+)\.evio\.(\d+)", os.path.basename(f)) for f in evios]
    if not all(parsed):
        sys.exit("[ERROR] file names must be prad_0<run>.evio.<split>")
    runs = {int(m.group(1)) for m in parsed}
    if len(runs) != 1:
        sys.exit(f"[ERROR] files of several runs: {sorted(runs)}")
    run, splits = runs.pop(), [int(m.group(2)) for m in parsed]
    rt = RUN_TYPE_OF.get(run)
    if rt is None:
        sys.exit(f"[ERROR] run {run} is not in RUN_TYPE_OF")
    gate = ((lambda b: b & PHYSICS_MASK_ANY != 0 and not b & LMS_BIT) if rt in PHYSICS_ANY
            else (lambda b: b == PHYSICS_BITS))
    cfg = run_config(run)
    calib = str(DB_DIR / cfg["calib"])
    core = build(run, calib)
    ebeam = float(core.run_cfg.Ebeam)
    if abs(ebeam / cfg["beam"] - 1.0) > EBEAM_TOL:
        sys.exit(f"[ERROR] builder beam energy {ebeam:g} MeV != runinfo {cfg['beam']:g} MeV")
    hycal = core.hycal
    mods = [hycal.module(i) for i in range(hycal.module_count())]
    names, n_mod = np.array([m.name for m in mods]), len(mods)
    pwo4 = np.array([bool(m.is_pwo4()) for m in mods])
    dead = set(cfg["dead"])
    if not dead <= set(names):
        sys.exit(f"[ERROR] dead modules not in the HyCal map: {sorted(dead - set(names))}")
    gain = np.array([m.energize(1.0) for m in mods])          # MeV per integral unit (energize is linear, checked)
    if any(abs(m.energize(1000.0) - 1000.0 * g) > 1e-6 * abs(1000.0 * g) + 1e-9 for m, g in zip(mods, gain)):
        sys.exit("[ERROR] Module.energize is not linear")
    bank_rt = rt.removesuffix("_old")                          # the old 3.49 GeV runs use the 3p49GeV_H2 bank (the old bank has no a_band)
    bank_path = BANK_ROOT / bank_rt / f"templates_{bank_rt}.npz"
    bank = Bank(bank_path, names, pwo4)
    mf = MF(bank)
    t50, iint = template_columns(bank_path, names)
    ptab = phase_table(bank)
    wave = dec.WaveAnalyzer(core.daq_cfg.wave_cfg)
    crate_map = {int(r.tag): int(r.crate) for r in core.daq_cfg.roc_tags}
    print(f"[setup] run {run} {rt}, {len(evios)} splits, beam {ebeam:g} MeV, calib {cfg['calib']}, "
          f"bank {bank_path.relative_to(PYDIR)}, dead {sorted(dead)}", flush=True)

    cache = {}

    def channels(fe):
        """(module index, channel) of every channel of an event that belongs to a live HyCal module."""
        for ri in range(fe.nrocs):
            roc = fe.roc(ri)
            crate = crate_map.get(roc.tag) if roc.present else None
            if crate is None:
                continue
            for sl in roc.present_slots():
                slot = roc.slot(sl)
                for c in slot.present_channels():
                    mi = cache.get((roc.tag, sl, c))
                    if mi is None:
                        m = hycal.module_by_daq(crate, sl, c)
                        mi = cache[(roc.tag, sl, c)] = m.index if (m is not None and m.is_hycal() and m.name not in dead) else -1
                    if mi >= 0:
                        yield mi, slot.channel(c)

    def events(path):
        """Decoded physics events of an EVIO file."""
        ch = dec.EvChannel()
        ch.set_config(core.daq_cfg)
        if ch.open_auto(str(path)) != dec.Status.success:
            sys.exit(f"[ERROR] cannot open {path}")
        while ch.read() == dec.Status.success:
            if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
                continue
            for i in range(ch.get_n_events()):
                d = ch.decode_event(i, with_ssp=False)
                if d["ok"]:
                    yield d["event"]

    def reference(path):
        """Noise reference of an EVIO file from the pedestal samples of its first N_REF_EVENTS events."""
        rows, mi_ = [], []
        for n, fe in enumerate(events(path), 1):
            for mi, cd in channels(fe):
                if pwo4[mi] and cd.nsamples == NS_EXPECT:
                    rows.append(np.asarray(cd.samples, dtype=np.uint16))
                    mi_.append(mi)
            if n >= N_REF_EVENTS:
                break
        if not rows:                                           # no waveform in the file: every module on the floor sigma
            return np.full(n_mod, np.nan, F32), np.full(n_mod, SIGMA_MIN, F32)
        Y, mi_ = np.stack(rows), np.array(mi_)
        good = ~bad_waveform(Y)
        return noise_reference(Y[good, :PED_NSAMPLES], mi_[good], n_mod)

    hits, rows, row_ev, row_mi, wa_buf = [], [], [], [], []   # per block: hit arrays; PbWO4 waveforms; WaveAnalyzer hits
    ev_number, ev_split = [], []
    n_bad, n_trunc, n_base_ref = (np.zeros(n_mod, np.int64) for _ in range(3))   # per module: bad waveforms, truncated main pulses, starts from b_ref
    n_fb = n_seen = 0

    def wa_hit(ev, mi, samples):
        for pk in wave.analyze(samples)[2]:
            if PRE_LO < pk.time < PRE_HI:
                wa_buf.append((ev, mi, pk.time, gain[mi] * pk.integral))

    def flush(b_ref, sig_ref):
        """Matched filter on the PbWO4 waveforms of the block, then its WaveAnalyzer hits."""
        nonlocal n_fb
        if rows:
            Y, ev, mi = np.stack(rows), np.array(row_ev), np.array(row_mi)
            rows.clear(), row_ev.clear(), row_mi.clear()
            bad = bad_waveform(Y)
            n_bad[:] += np.bincount(mi[bad], minlength=n_mod)
            Y, ev, mi = Y[~bad], ev[~bad], mi[~bad]
            t0, a, v, trunc, far = run_ref(mf, Y, bank.row[mi], sig_ref[mi], b_ref[mi])
            n_base_ref[:] += np.bincount(mi[far], minlength=n_mod)
            tm = v.any(axis=1) & trunc[np.arange(len(Y)), np.where(v, a, -1.0).argmax(axis=1)]   # main pulse's fit window cut by the waveform end
            n_trunc[:] += np.bincount(mi[tm], minlength=n_mod)
            t = t0 + t50[mi][:, None]
            r, j = np.nonzero(v & ~tm[:, None] & (t > PRE_LO) & (t < PRE_HI))
            x = (-(t0[r, j] % DT) / DT * N_PHASE) % N_PHASE          # table position of the onset phase
            i0 = x.astype(np.int64)
            ph = ptab[bank.row[mi[r]], i0] * (1 - x + i0) + ptab[bank.row[mi[r]], i0 + 1] * (x - i0)
            hits.append((ev[r], mi[r], t[r, j], gain[mi[r]] * iint[mi[r]] * ph * a[r, j], a[r, j], v.sum(1)[r]))
            fb = np.nonzero((Y >= OVERFLOW).any(axis=1) & ~v.any(axis=1))[0]   # clipped, no usable pulse
            n_fb += fb.size
            for q in fb:
                wa_hit(int(ev[q]), int(mi[q]), Y[q])
        if wa_buf:
            w = np.array(wa_buf, dtype=np.float64)
            wa_buf.clear()
            hits.append((w[:, 0].astype(np.int64), w[:, 1].astype(np.int64), w[:, 2], w[:, 3],
                         np.full(len(w), np.nan), np.zeros(len(w), dtype=np.int64)))

    for split, evio in zip(splits, evios):
        b_ref, sig_ref = reference(evio)                       # blocks never span two files: each has its own reference
        for fe in events(evio):
            n_seen += 1
            if n_seen % 20000 == 0:
                print(f"  [progress] split {split}: {n_seen} events, {time.monotonic() - t_start:.0f} s", flush=True)
            if not gate(int(fe.info.trigger_bits)):
                continue
            ev = len(ev_number)
            ev_number.append(int(fe.info.event_number))
            ev_split.append(split)
            for mi, cd in channels(fe):
                if cd.nsamples <= 0:
                    continue
                if pwo4[mi] and cd.nsamples == NS_EXPECT:
                    rows.append(np.asarray(cd.samples, dtype=np.uint16))
                    row_ev.append(ev)
                    row_mi.append(mi)
                elif bad_waveform(np.asarray(cd.samples, dtype=np.uint16)[None])[0]:
                    n_bad[mi] += 1
                else:
                    wa_hit(ev, mi, cd.samples)
            if len(ev_number) % BLOCK == 0:
                flush(b_ref, sig_ref)
        flush(b_ref, sig_ref)
    h_ev, h_mod, h_t, h_e, h_a, h_k = (np.concatenate([h[q] for h in hits]) if hits else np.zeros(0) for q in range(6))
    o = np.lexsort((-h_e, h_mod, h_ev))                        # per event and module, largest energy first
    h_ev, h_mod, h_t, h_e, h_a, h_k = (x[o] for x in (h_ev, h_mod, h_t, h_e, h_a, h_k))
    extra = np.r_[False, (h_ev[1:] == h_ev[:-1]) & (h_mod[1:] == h_mod[:-1])]   # not the module's largest in-window pulse
    wall = time.monotonic() - t_start
    with open(bank_path, "rb") as f:
        sha = hashlib.sha1(f.read()).hexdigest()
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    np.savez_compressed(out, run=run, run_type=rt, splits=np.array(splits), beam_mev=ebeam, calib=calib,
                        dead_modules=np.array(sorted(dead), dtype=str), bank=str(bank_path.relative_to(PYDIR)),
                        bank_sha1=sha, n_seen=n_seen, n_pass=len(ev_number), n_clip_fallback=n_fb,
                        mod_n_base_ref=n_base_ref, mod_n_bad=n_bad, mod_n_trunc=n_trunc, wall_s=wall,
                        mod_name=names, mod_pwo4=pwo4, ev_number=np.array(ev_number, dtype=np.int64),
                        ev_split=np.array(ev_split, dtype=np.int32), h_ev=h_ev.astype(np.int32),
                        h_mod=h_mod.astype(np.int16), h_t=h_t.astype(np.float32), h_e=h_e.astype(np.float32),
                        h_a=h_a.astype(np.float32), h_k=h_k.astype(np.int8), h_extra=extra)
    top = lambda c: ", ".join(f"{names[i]} {c[i]}" for i in np.argsort(-c)[:3])   # noqa: E731
    print(f"[done] {out}: {n_seen} events, {len(ev_number)} physics, {h_ev.size} in-window hits ({int(extra.sum())} extra), "
          f"{n_fb} clipped channels on WaveAnalyzer, {int(n_base_ref.sum())} baseline starts from b_ref; {wall:.0f} s\n"
          f"[bad waveforms] {n_bad.sum()} (most: {top(n_bad)}); [truncated main pulses] {n_trunc.sum()} (most: {top(n_trunc)})", flush=True)


# clustering (report) ==========================================================
_CL = {}                                                       # per worker process: (run, gates) -> clusterers, geometry
REC = ("ev", "seed", "seed_e", "e", "er", "x", "y", "pu", "clip", "mi", "me")   # per selected cluster


def clusterers(run, calib, gates):
    """Legacy and gated HyCalCluster instances of the run (builder fields CLUSTER_FIELDS, MIN_CENTER_E, its Geant4 profile and
    non-linearity) and the module geometry."""
    key = (run, tuple(gates))
    if key in _CL:
        return _CL[key]
    core = build(run, calib)
    hycal, base = core.hycal, core.hycal_cluster_cfg

    def make(window):
        cfg = det.HyCalClusterConfig()
        for f in CLUSTER_FIELDS:
            setattr(cfg, f, getattr(base, f))
        cfg.min_center_energy = MIN_CENTER_E
        cfg.seed_time_window = window
        cl = det.HyCalCluster(hycal)
        cl.set_config(base)                # installs the builder's Geant4 profile (a C++-only field) ...
        cl.set_config(cfg)                 # ... which this one keeps while setting every other field
        return cl

    cls = {"legacy": make(LEGACY_WINDOW)} | {f"g{g:g}": make(g) for g in gates}
    # profile check: three PbWO4 modules in a row with 1000/50/200 MeV share their energy differently with the Geant4 profile
    simple = det.HyCalCluster(hycal)
    simple.set_config(cls["legacy"].get_config())
    w = {(m.row, m.column): m.index for m in (hycal.module(i) for i in range(hycal.module_count())) if m.is_pwo4()}
    r0, c0 = next(k for k in sorted(w) if (k[0], k[1] + 1) in w and (k[0], k[1] + 2) in w)
    split = []
    for cl in (simple, cls["legacy"]):
        cl.clear()
        for dc, e in ((0, 1000.0), (1, 50.0), (2, 200.0)):
            cl.add_hit(w[(r0, c0 + dc)], e, 0.0)
        cl.form_clusters()
        split.append(sorted(round(x.hit.energy, 3) for x in cl.reconstruct_matched()))
        cl.clear()
    if split[0] == split[1]:
        sys.exit("[ERROR] the clusterer did not keep the builder's Geant4 profile")
    mods = [hycal.module(i) for i in range(hycal.module_count())]
    geo = dict(x=np.array([m.x for m in mods]), y=np.array([m.y for m in mods]), row=np.array([m.row for m in mods]),
               col=np.array([m.column for m in mods]), name=np.array([m.name for m in mods]),
               pwo4=np.array([bool(m.is_pwo4()) for m in mods]), sx=np.array([m.size_x for m in mods]),
               sy=np.array([m.size_y for m in mods]), hz=float(core.run_cfg.hycal_z), ebeam=float(core.run_cfg.Ebeam))
    _CL[key] = (cls, geo)
    return _CL[key]


def objects(v):
    """1-D object array of arrays (np.array would stack equal-length members into a matrix)."""
    a = np.empty(len(v), object)
    a[:] = v
    return a


def take(d, mask):
    return {q: v[mask] for q, v in d.items()}


def cat(ds, keys=None):
    return {q: np.concatenate([d[q] for d in ds]) for q in (keys or ds[0])}


def pool_modules(geo):
    """[n modules] bool: seed module centre radius < INNER_R and outside EXCLUDE_ROWCOL (1-based rows and columns, inclusive)."""
    rl, rh, cl, ch = EXCLUDE_ROWCOL
    row, col = geo["row"] + 1, geo["col"] + 1
    return (geo["x"] ** 2 + geo["y"] ** 2 < INNER_R ** 2) & ~((rl <= row) & (row <= rh) & (cl <= col) & (col <= ch))


def cluster_events(h, cls, geo, signal_min, seed_min):
    """Per path, the pool clusters (seed in pool_modules, seed energy >= seed_min) of the events of h with exactly one cluster
    above signal_min (fields REC): cluster energy e with the calibration's non-linearity factor, er without it (sum of the raw
    member energies me of modules mi), pile-up flag (a member module with >= 2 in-window pulses), clip flag (PbWO4 seed on
    WaveAnalyzer).  Also {path: number of such events}."""
    ev, pwo4, ok_mod = h["h_ev"], geo["pwo4"], pool_modules(geo)
    out = {p: {q: [] for q in REC} for p in cls}
    n_single = dict.fromkeys(cls, 0)
    if ev.size:
        lo = np.flatnonzero(np.r_[True, ev[1:] != ev[:-1]])
        hi = np.r_[lo[1:], ev.size]
        keep = np.add.reduceat(h["h_e"].astype(np.float64), lo) > 0.5 * signal_min   # a cluster above signal_min needs > half of it in hits
        mod, t, e, k, extra = (h[q].tolist() for q in ("h_mod", "h_t", "h_e", "h_k", "h_extra"))
        for s, f in zip(lo[keep].tolist(), hi[keep].tolist()):
            multi = {m for m, c in Counter(mod[s:f]).items() if c >= 2}
            wa = {m for m, kk in zip(mod[s:f], k[s:f]) if kk == 0}
            for p, cl in cls.items():
                cl.clear()
                for i in range(s, f):
                    if p != "legacy" or not extra[i]:
                        cl.add_hit(mod[i], e[i], t[i])
                cl.form_clusters()
                big = [r for r in cl.reconstruct_matched() if r.hit.energy > signal_min]
                if len(big) != 1:
                    continue
                r = big[0]
                n_single[p] += 1
                seed = r.cluster.center.index
                if not (ok_mod[seed] and r.cluster.center.energy >= seed_min):
                    continue
                mi = np.array([x.index for x in r.cluster.hits], np.int64)
                me = np.array([x.energy for x in r.cluster.hits])
                for q, v in zip(REC, (ev[s], seed, r.cluster.center.energy, r.hit.energy, me.sum(), r.hit.x, r.hit.y,
                                      bool(multi.intersection(mi.tolist())), seed in wa and pwo4[seed], mi, me)):
                    out[p][q].append(v)
    return {p: {q: objects(v) if q in ("mi", "me") else np.array(v) for q, v in d.items()} for p, d in out.items()}, n_single


# kinematics and fits (benchmark_hycal_timing_test_20261021.py) ================
def energy_loss(theta):
    """PhysicsTools::EnergyLoss: target and window material loss (MeV) at polar angle theta (degrees)."""
    c = np.cos(np.radians(theta))
    return (0.5 * 1.6 + 0.12 * 1.6 + 0.1 * 2.0 + 0.48 * 1.8) * np.where(c > 0.01, 1.0 / np.maximum(c, 0.01), 100.0)


def e_elastic(theta, ebeam, m_target):
    """Elastic E' = E M / (M + E (1 - cos theta)) minus the energy loss."""
    return ebeam * m_target / (m_target + ebeam * (1.0 - np.cos(np.radians(theta)))) - energy_loss(theta)


def e_moller(theta, ebeam):
    """Moller E_ee(theta) minus the energy loss; 0 where the denominator is not positive."""
    c2 = np.cos(np.radians(theta)) ** 2
    g = ebeam / M_ELECTRON
    den = (g + 1.0) - (g - 1.0) * c2
    return np.where(den > 0, M_ELECTRON * ((g + 1.0) + (g - 1.0) * c2) / np.where(den > 0, den, 1.0) - energy_loss(theta), 0.0)


def theta_deg(x, y, hz):
    """Polar angle (degrees) of a point on HyCal (no shower depth, no beam offset)."""
    return np.degrees(np.arctan(np.hypot(x, y) / hz))


def gauss(x, A, mu, sigma):
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def crystal_ball(x, A, mu, sigma, alpha, n):
    """Crystal Ball with a left power-law tail, alpha > 0, n > 1."""
    z = (np.asarray(x, dtype=float) - mu) / sigma
    a = abs(alpha)
    arg = n / a - a - z
    with np.errstate(over="ignore", invalid="ignore"):
        tail = (n / a) ** n * np.exp(-0.5 * a * a) * np.where(arg > 0, arg, 1e-30) ** (-n)
    return A * np.where(z > -alpha, np.exp(-0.5 * z * z), tail)


def hist(values, lo, hi, bw):
    edges = np.arange(lo, hi + bw, bw)
    counts, _ = np.histogram(values, bins=edges)
    return 0.5 * (edges[:-1] + edges[1:]), counts


def fit_curve(f, x, y, p0, lo, hi, **kw):
    """curve_fit within bounds: (parameters, covariance), (None, None) when the fit fails."""
    try:
        return curve_fit(f, x, y.astype(float), p0=p0, bounds=(lo, hi), **kw)
    except Exception:
        return None, None


def gauss_fixed(c, y, lo, hi, bw, kS):
    """Gaussian within +- max(FIX_HALF kS, 4 bins) of the argmax: (mu, sigma)."""
    if y.sum() == 0:
        return NAN2
    mu0 = float(c[int(np.argmax(y))])
    win = max(FIX_HALF * kS, 4.0 * bw)
    m = (c > mu0 - win) & (c < mu0 + win)
    if y[m].sum() < 20:
        return NAN2
    p, _ = fit_curve(gauss, c[m], y[m], [float(y[m].max()), mu0, max(bw, 60.0 * kS)], [0, lo, 1.0], [np.inf, hi, win])
    return NAN2 if p is None else (float(p[1]), abs(float(p[2])))


def gauss_2sigma(c, y, mu, sig):
    """Iterated Gaussian on the bins within mu +- G2_NSIGMA sigma, from (mu, sig), until mu and sigma change by < G2_TOL
    (at most G2_MAX_ITER refits; a failed refit keeps the last good fit)."""
    if not (math.isfinite(mu) and math.isfinite(sig) and sig > 0):
        return mu, sig
    for _ in range(G2_MAX_ITER):
        lo, hi = mu - G2_NSIGMA * sig, mu + G2_NSIGMA * sig
        m = (c >= lo) & (c <= hi)
        if m.sum() < 4 or y[m].sum() < 20:
            break
        p, _ = fit_curve(gauss, c[m], y[m], [float(y[m].max()), mu, sig], [0, lo, 1.0], [np.inf, hi, hi - lo])
        if p is None:
            break
        done = abs(p[1] - mu) < G2_TOL * abs(mu) and abs(abs(p[2]) - sig) < G2_TOL * sig
        mu, sig = float(p[1]), abs(float(p[2]))
        if done:
            break
    return mu, sig


def cb_single(c, y, lo, hi, bw, kS, mu0):
    """Crystal Ball on [mu0 - win, mu0 + win / 2], win = max(400 kS, 12 bins): (mu, sigma)."""
    win = max(400.0 * kS, 12.0 * bw)
    m = (c > mu0 - win) & (c < mu0 + 0.5 * win)
    if y[m].sum() < 50:
        return NAN2
    p, _ = fit_curve(crystal_ball, c[m], y[m], [float(y[m].max()), mu0, max(60.0 * kS, 3 * bw), 1.5, 2.0],
                     [0, lo - 200.0 * kS, 1.0, 0.1, 1.01], [np.inf, hi + 200.0 * kS, win, 20.0, 50.0], maxfev=10000)
    return NAN2 if p is None else (float(p[1]), abs(float(p[2])))


def cb_two_peak(c, y, lo, hi, bw, kS, mu0, sg, em0, em_range):
    """Elastic + Moller Crystal Balls (shared alpha, n), Moller mean seeded at em0; the single Crystal Ball if the Moller peak is not
    supported (em0 outside [lo, hi] or not below mu0 - sg, amplitude < 2 x its error, or mu_M off em_range).  (mu, sigma, mode)."""
    single = (*cb_single(c, y, lo, hi, bw, kS, mu0), "single")
    if not (math.isfinite(sg) and sg > 0 and math.isfinite(em0)) or not lo <= em0 <= hi or em0 > mu0 - sg:
        return single
    win = max(400.0 * kS, 12.0 * bw)
    m = (c > min(mu0 - win, em0 - 4.0 * sg)) & (c < mu0 + 0.5 * win)
    x, yy = c[m], y[m].astype(float)
    if yy.sum() < 50:
        return single
    near_el, near_mo = np.abs(x - mu0) < max(sg, bw), np.abs(x - em0) < max(sg, bw)
    a_el = float(yy[near_el].max()) if near_el.any() else float(yy.max())
    a_mo = float(yy[near_mo].mean()) if near_mo.any() else 1.0

    def two_cb(xx, Ael, muel, sgel, al, nn, Amo, dlt, sgmo):
        return crystal_ball(xx, Ael, muel, sgel, al, nn) + crystal_ball(xx, Amo, muel - dlt, sgmo, al, nn)

    p, cov = fit_curve(two_cb, x, yy, [a_el, mu0, sg, 1.5, 2.0, max(a_mo, 1.0), mu0 - em0, sg],
                       [0, mu0 - 3.0 * sg, 1.0, 0.1, 1.01, 0, 1.0 * sg, 0.5 * sg],
                       [np.inf, mu0 + 3.0 * sg, win, 20.0, 50.0, np.inf, (mu0 - em0) + 4.0 * sg, 2.0 * sg], maxfev=20000)
    if p is None or not np.all(np.isfinite(p)):
        return single
    vA = float(cov[5, 5]) if np.ndim(cov) == 2 else float("nan")
    eA = math.sqrt(vA) if (math.isfinite(vA) and vA >= 0) else float("nan")
    if not p[5] > 2.0 * eA:
        return single
    mu_m = p[1] - p[6]
    res = math.hypot(MOLLER_RES_A / math.sqrt(max(mu_m, 1.0) / 1000.0), MOLLER_RES_B)
    if not em_range[0] - MOLLER_TIE_NSIGMA * res * mu_m <= mu_m <= em_range[1] + MOLLER_TIE_NSIGMA * res * mu_m:
        return single
    return float(p[1]), abs(float(p[2])), "two_peak"


def module_cb(vals, lo, hi, seed, bw, sigma_glob, kS):
    """Per-module Crystal Ball on [lo, hi]: (mu, sigma), None unless lo < mu < hi and 0.25 <= sigma / sigma_glob <= 2."""
    c, y = hist(vals, lo, hi, bw)
    if y.sum() == 0:
        return None
    p, _ = fit_curve(crystal_ball, c, y, [float(y.max()), seed, 50.0 * kS, 1.5, 2.0],
                     [0, lo, 5.0 * kS, 0.1, 1.01], [np.inf, hi, hi - lo, 20.0, 50.0], maxfev=5000)
    if p is None:
        return None
    mu, sig = float(p[1]), abs(float(p[2]))
    ok = np.all(np.isfinite(p)) and lo < mu < hi and 0.25 * sigma_glob <= sig <= 2.0 * sigma_glob
    return (mu, sig) if ok else None


def pool_mass(rt):
    return TARGET_MASS["Ta" if "Ta" in rt else "H2"]


def module_set(seed, E, geo, m_target, kE, kS):
    """Elastic-window module set of a pool (seed modules, cluster energies): {module: E'} of the used modules, the mean E', and
    dict(win={module: (lo, hi) = mu_M +- WINDOW_NSIGMA sigma_M}, sigma_glob, mu_global)."""
    nan = float("nan")
    ebeam, hz, bw = geo["ebeam"], geo["hz"], BIN * kS
    sids = np.unique(seed)
    if sids.size == 0:
        return {}, nan, dict(win={}, sigma_glob=nan, mu_global=nan)
    th = theta_deg(geo["x"][sids], geo["y"][sids], hz)
    ep, ee = dict(zip(sids, e_elastic(th, ebeam, m_target))), dict(zip(sids, e_moller(th, ebeam)))
    ep_mean = float(np.mean(list(ep.values())))
    lo, hi = ep_mean - FIT_HALF * kE, ep_mean + FIT_HALF * kE
    c, y = hist(E, lo, hi, bw)
    mu_g, sig_g = gauss_fixed(c, y, lo, hi, bw, kS)
    mu_cb, sig_cb = cb_single(c, y, lo, hi, bw, kS, mu_g) if math.isfinite(mu_g) else NAN2
    sigma_glob = sig_cb if math.isfinite(sig_cb) else sig_g
    mu_global = mu_cb if math.isfinite(mu_cb) else mu_g
    if not (math.isfinite(sigma_glob) and math.isfinite(mu_global)):
        return {}, ep_mean, dict(win={}, sigma_glob=nan, mu_global=nan)
    gscale = mu_global / ep_mean
    used, win = {}, {}
    for s in sids:
        vals = E[seed == s]
        plo = mu_global * max(1.0 + PREFIT_RANGE[0] / kS, 0.5 * (1.0 + ee[s] / ep[s]))
        phi = mu_global * (1.0 + PREFIT_RANGE[1] / kS)
        sel = vals[(vals >= plo) & (vals <= phi)]
        if sel.size == 0:
            continue
        f = module_cb(vals, plo, phi, mu_global, bw, sigma_glob, kS) if sel.size >= MIN_FIT_EVENTS else None
        mu_m, sig_m = f if f else (float(np.median(sel)), sigma_glob)
        n_in = int(((vals >= mu_m - WINDOW_NSIGMA * sig_m) & (vals <= mu_m + WINDOW_NSIGMA * sig_m)).sum())
        if gscale * (ep[s] - ee[s]) < MOLLER_GAP_NSIGMA * sig_m or n_in < MIN_PER_MODULE:
            continue
        used[int(s)] = ep[s]
        win[int(s)] = (mu_m - WINDOW_NSIGMA * sig_m, mu_m + WINDOW_NSIGMA * sig_m)
    return used, ep_mean, dict(win=win, sigma_glob=sigma_glob, mu_global=mu_global)


def fit_path(pool, used, ep_mean, ref, geo, kE, kS):
    """Fits of the used-module spectrum on ref +- FIT_HALF kE: (results, spectrum); two-peak Crystal Ball (res_cb), Gaussian within
    2 sigma (res_g2), interquartile range / 1.349 (res_iqr), peak / E' - 1 (mu_ep), tail, pile-up and clip fractions, all in %."""
    bw, hz, ebeam = BIN * kS, geo["hz"], geo["ebeam"]
    u = np.isin(pool["seed"], list(used))
    E, xc, yc = pool["e"][u], pool["x"][u], pool["y"][u]
    if E.size == 0 or not math.isfinite(ref):
        return dict.fromkeys(("mu_fix", "sig_fix", "mu_g2", "sig_g2", "mu_cb", "sig_cb", "res_cb", "res_g2", "res_iqr",
                              "tail", "mu_ep", "pu", "clip"), float("nan")) | dict(n_used=int(E.size), cb_mode="none"), E
    lo, hi = ref - FIT_HALF * kE, ref + FIT_HALF * kE
    c, y = hist(E, lo, hi, bw)
    mu_f, sig_f = gauss_fixed(c, y, lo, hi, bw, kS)
    mu_2, sig_2 = gauss_2sigma(c, y, mu_f, sig_f)
    mu_cb, sig_cb, mode = (*NAN2, "none")
    if math.isfinite(mu_f):
        sc = mu_f / ep_mean
        em0 = float(e_moller(np.mean(theta_deg(xc, yc, hz)), ebeam)) * sc
        r_lo, r_hi = np.percentile(np.hypot(xc, yc), MOLLER_TIE_PCT)
        em_rng = tuple(float(e_moller(theta_deg(r, 0.0, hz), ebeam)) * sc for r in (r_hi, r_lo))
        mu_cb, sig_cb, mode = cb_two_peak(c, y, lo, hi, bw, kS, mu_f, sig_2, em0, em_rng)
    inr = E[(E >= lo) & (E <= hi)]
    q25, q50, q75 = np.percentile(inr, (25, 50, 75)) if inr.size else (np.nan,) * 3
    tail = float(np.mean(np.abs(inr - mu_cb) > TAIL_NSIGMA * sig_cb)) if inr.size and sig_cb > 0 else np.nan
    pct = lambda a, b: 100.0 * a / b if (b and math.isfinite(b) and b > 0) else float("nan")   # noqa: E731
    return dict(n_used=int(E.size), mu_fix=mu_f, sig_fix=sig_f, mu_g2=mu_2, sig_g2=sig_2, mu_cb=mu_cb, sig_cb=sig_cb,
                cb_mode=mode, res_cb=pct(sig_cb, mu_cb), res_g2=pct(sig_2, mu_2), res_iqr=pct((q75 - q25) / 1.349, q50),
                tail=100.0 * tail, mu_ep=pct(mu_cb, ep_mean) - 100.0, pu=100.0 * float(np.mean(pool["pu"][u])) if u.any() else np.nan,
                clip=100.0 * float(np.mean(pool["clip"][u])) if u.any() else np.nan), E


def peak_yield(pool, used, s_ref, common):
    """Pool events of the common modules within +- YIELD_NSIGMA s_ref of their module's E'."""
    if not (math.isfinite(s_ref) and s_ref > 0):
        return 0
    tgt = np.array([used.get(int(s), np.nan) if int(s) in common else np.nan for s in pool["seed"]])
    return int(np.sum(np.abs(pool["e"] / tgt - 1.0) <= YIELD_NSIGMA * s_ref))


def evaluate(pools, geo, mass, kE, kS):
    """Module set and fits of every path's pool on a fit range centred on the legacy pool's mean E'; peak yield = events within
    YIELD_NSIGMA s_ref of E' (s_ref: legacy Gaussian width) on the modules both paths use.  Returns {path: fields}, {path: spectrum}."""
    sets = {p: module_set(pl["seed"].astype(np.int64), pl["e"], geo, mass, kE, kS) for p, pl in pools.items()}
    ref = sets["legacy"][1]
    fits = {p: fit_path(pl, *sets[p][:2], ref, geo, kE, kS) for p, pl in pools.items()}
    s_ref = fits["legacy"][0]["sig_fix"] / fits["legacy"][0]["mu_fix"]
    out = {}
    for p, pl in pools.items():
        common = set(sets["legacy"][0]) & set(sets[p][0])
        out[p] = dict(fits[p][0], n_modules=len(sets[p][0]), eprime=sets[p][1], n_common=len(common),
                      yield_n=peak_yield(pl, sets[p][0], s_ref, common),
                      yield_legacy_n=peak_yield(pools["legacy"], sets["legacy"][0], s_ref, common),
                      modules=",".join(sorted(str(geo["name"][s]) for s in sets[p][0])))
    return out, {p: f[1] for p, f in fits.items()}


# per-crystal gains (percrystal_calib.py solve, windows refitted) ==============
def normal_equations(ptr, midx, mE, tgt, n_mod):
    """A^T A and A^T t for events given as CSR members (ptr, compact module index, member energy, target)."""
    counts = np.diff(ptr)
    N = np.zeros(n_mod * n_mod)
    for c in np.unique(counts):
        evs = np.flatnonzero(counts == c)
        step = max(1, 2_000_000 // int(c * c))
        for a in range(0, evs.size, step):
            pos = ptr[evs[a:a + step]][:, None] + np.arange(c)[None, :]
            m2, e2 = midx[pos], mE[pos]
            N += np.bincount((m2[:, :, None] * n_mod + m2[:, None, :]).ravel(),
                             weights=(e2[:, :, None] * e2[:, None, :]).ravel(), minlength=n_mod * n_mod)
    b = np.bincount(midx, weights=mE * np.repeat(tgt, counts), minlength=n_mod)
    return N.reshape(n_mod, n_mod), b


def solve_gains(pool, geo, m_target, kE, kS):
    """One gain per member module so that sum_i g_i E_i of the events inside their seed module's elastic window matches E' (prior
    lambda (g_i - g0)^2, g0 = mean E' / mu_global, lambda = MIN_PER_MODULE (GAIN_SEED_FRAC mean E')^2 / sigma_glob^2 of the g = 1
    pass); windows from module_set, refitted in the first GAIN_N_REFIT passes.  Modules never selected keep g0.  Returns
    dict(g, ecal = sum g_i E_i per cluster, rows, g0, iterations, converged, n_free) or None."""
    n_all = geo["x"].size
    seed, E = pool["seed"].astype(np.int64), pool["e"].astype(float)
    counts = np.array([m.size for m in pool["mi"]], np.int64)
    midx = np.concatenate(pool["mi"]).astype(np.int64)
    mE = np.concatenate(pool["me"]).astype(float)
    evid = np.repeat(np.arange(E.size), counts)
    ecal_of = lambda g: np.bincount(evid, weights=mE * g[midx], minlength=E.size)   # noqa: E731
    windows = lambda e: module_set(seed, e, geo, m_target, kE, kS)                    # noqa: E731
    used0, ep_mean, w0 = windows(E)
    sg0, mu0 = w0["sigma_glob"], w0["mu_global"]
    if not (math.isfinite(sg0) and math.isfinite(mu0) and sg0 > 0):
        return None
    g0, w = ep_mean / mu0, 1.0 / sg0 ** 2
    lam = MIN_PER_MODULE * (GAIN_SEED_FRAC * ep_mean) ** 2 * w
    g, converged = np.ones(n_all), False
    for it in range(1, GAIN_MAX_ITER + 1):
        ecal = ecal_of(g)
        if it <= GAIN_N_REFIT:
            used, _, wf = (used0, ep_mean, w0) if it == 1 else windows(ecal)
            sg = wf["sigma_glob"]
            if not (math.isfinite(sg) and math.isfinite(wf["mu_global"])):
                return None
            lo, hi, tg = (np.full(n_all, np.nan) for _ in range(3))
            for m_, (a, b) in wf["win"].items():
                lo[m_], hi[m_], tg[m_] = a, b, used[m_]
        sel = (ecal >= lo[seed]) & (ecal <= hi[seed])
        if not sel.any():
            return None
        msel = sel[evid]
        uniq, inv = np.unique(midx[msel], return_inverse=True)
        ptr_s = np.r_[0, np.cumsum(counts[sel])]
        Nm, b = normal_equations(ptr_s, inv, mE[msel], tg[seed[sel]], uniq.size)
        Nm, b = Nm * w, b * w
        D = np.diag(Nm).copy()
        A = Nm + lam * np.eye(uniq.size)
        gsol = np.linalg.solve(A, b + lam * g0)
        sigma_g = np.sqrt(np.diag(np.linalg.inv(A))) * (sg / sg0)
        pf = lam / (D + lam)
        gnew = np.full(n_all, g0)
        gnew[uniq] = gsol
        free = pf < GAIN_FREE_FRAC
        dg = np.abs(gnew[uniq] - g[uniq])[free]
        converged = bool(free.any() and (dg < np.maximum(GAIN_DG_TOL, GAIN_TOL_SIGMA * sigma_g[free])).all())
        g = gnew
        if converged:
            break
    n_seed = np.bincount(np.searchsorted(uniq, seed[sel]), minlength=uniq.size)
    n_member = np.bincount(inv, minlength=uniq.size)
    rows = [dict(module=str(geo["name"][m_]), g=gsol[j], sigma_g=sigma_g[j], n_seed=int(n_seed[j]),
                 n_member=int(n_member[j]), prior_frac=pf[j]) for j, m_ in enumerate(uniq.tolist())]
    return dict(g=g, ecal=ecal_of(g), rows=rows, g0=g0, iterations=it, converged=converged, n_free=int(free.sum()))


# position correction (position_corr.py map / fit --per-module) ================
def uv(pool, geo):
    """Impact position in the seed module: (centroid - seed centre) / module size."""
    s = pool["seed"].astype(np.int64)
    return (pool["x"] - geo["x"][s]) / geo["sx"][s], (pool["y"] - geo["y"][s]) / geo["sy"][s]


def bin_index(u, v, bins, rng):
    """Bin of (u, v) over [-rng, rng]^2; an exact 0 (single-module centroid) alternates over the two central bins."""
    w = 2 * rng / bins
    iu, iv = np.floor((u + rng) / w).astype(int), np.floor((v + rng) / w).astype(int)
    if bins % 2 == 0:
        idx = np.arange(u.size)
        iu = np.where(u == 0, bins // 2 - 1 + idx % 2, iu)
        iv = np.where(v == 0, bins // 2 - 1 + (idx // 2) % 2, iv)
    return iu, iv, (iu >= 0) & (iu < bins) & (iv >= 0) & (iv < bins)


def groups(flat):
    """{value: indices} of an integer array."""
    order = np.argsort(flat, kind="stable")
    ids, start = np.unique(flat[order], return_index=True)
    return dict(zip(ids.tolist(), np.split(order, start[1:])))


def peak_r(r, kind, tref, kS):
    """Peak of r with the report's fitters on r x tref (gauss: iterated Gaussian within 2 sigma; cb: Crystal Ball)."""
    x, lo, hi, bw = r * tref, tref * PC_R[0], tref * PC_R[1], BIN * kS
    c, y = hist(x, lo, hi, bw)
    mu, sig = gauss_fixed(c, y, lo, hi, bw, kS)
    if math.isfinite(mu):
        mu, sig = gauss_2sigma(c, y, mu, sig) if kind == "gauss" else cb_single(c, y, lo, hi, bw, kS, mu)
    return mu / tref, sig / tref


def fit_map(u, v, r, kind, tref, kS):
    """Bilinear map: per bin the peak of r (bins with >= PC_MIN_N events), the others the mean of their filled neighbours,
    normalised to an event-weighted mean of 1.  Returns (R, n, filled)."""
    B = PC_BINS
    iu, iv, inb = bin_index(u, v, B, PC_RANGE)
    R, n = np.full((B, B), np.nan), np.zeros((B, B), int)
    for b, idx in groups(iu[inb] * B + iv[inb]).items():
        i, j = divmod(b, B)
        n[i, j] = idx.size
        if idx.size >= max(PC_MIN_FIT, PC_MIN_N):
            mu, sig = peak_r(r[inb][idx], kind, tref, kS)
            if math.isfinite(mu) and math.isfinite(sig):
                R[i, j] = mu
    if not np.isfinite(R).any():
        sys.exit("[ERROR] position correction: no populated map bin")
    filled = ~np.isfinite(R)
    while np.isnan(R).any():
        new = R.copy()
        for i, j in zip(*np.nonzero(np.isnan(R))):
            blk = R[max(0, i - 1):i + 2, max(0, j - 1):j + 2]
            if np.isfinite(blk).any():
                new[i, j] = np.nanmean(blk)
        R = new
    return R / ((n * R).sum() / n.sum()), n, filled


def poly_r(c, u, v):
    """c0 (1 + c1 u^2 + c2 v^2 + c3 u^2 v^2 + c4 u^4 + c5 v^4 + c6 u + c7 v), u and v clipped to +- PC_CLAMP."""
    u, v = np.clip(u, -PC_CLAMP, PC_CLAMP), np.clip(v, -PC_CLAMP, PC_CLAMP)
    return c[0] * (1 + c[1] * u ** 2 + c[2] * v ** 2 + c[3] * u ** 2 * v ** 2 + c[4] * u ** 4 + c[5] * v ** 4
                   + c[6] * u + c[7] * v)


def trunc_mean(r, max_iter=20):
    """Iterated +- PC_TRUNC_K sigma truncated mean: (mu, sigma, n_in) or None."""
    k = PC_TRUNC_K
    ctrunc = 1 - 2 * k * math.exp(-k * k / 2) / math.sqrt(2 * math.pi) / math.erf(k / math.sqrt(2))
    mu = float(np.median(r))
    sig = 1.4826 * float(np.median(np.abs(r - mu)))
    if r.size < 5 or not sig > 0:
        return None
    for _ in range(max_iter):
        sel = np.abs(r - mu) <= k * sig
        if sel.sum() < 5:
            return None
        mu2, sig2 = float(r[sel].mean()), float(r[sel].std(ddof=1)) / math.sqrt(ctrunc)
        done = abs(mu2 - mu) < 1e-7 and abs(sig2 - sig) < 1e-4 * sig
        mu, sig = mu2, sig2
        if done:
            break
    return mu, sig, int((np.abs(r - mu) <= k * sig).sum())


def fit_module(u, v, r):
    """Own polynomial of one seed module: weighted least squares on the PC_GRID x PC_GRID truncated means, R averaging 1 over the
    module's clusters; None if too few bins or R leaves PC_ACCEPT."""
    cen = -PC_CLAMP + (np.arange(PC_GRID) + 0.5) * 2 * PC_CLAMP / PC_GRID
    iu, iv, inb = bin_index(u, v, PC_GRID, PC_CLAMP)
    pts = []
    for b, idx in groups(np.where(inb, iu * PC_GRID + iv, -1)).items():
        t = trunc_mean(r[idx]) if b >= 0 and idx.size >= PC_GRID_MIN else None
        if t:
            pts.append((cen[b // PC_GRID], cen[b % PC_GRID], t[0], t[1] / math.sqrt(t[2])))
    if len(pts) < PC_GRID_BINS:
        return None
    uc, vc, Rv, se = np.array(pts).T
    X = np.stack([np.ones_like(uc), uc ** 2, vc ** 2, uc ** 2 * vc ** 2, uc ** 4, vc ** 4, uc, vc], axis=-1)
    a = np.linalg.lstsq(X / se[:, None], Rv / se, rcond=None)[0]
    c = np.r_[a[0], a[1:] / a[0]]
    c[0] /= poly_r(c, u, v).mean()
    g = np.linspace(-PC_CLAMP, PC_CLAMP, 13)
    Rg = poly_r(c, *np.meshgrid(g, g))
    ok = np.isfinite(c).all() and PC_ACCEPT[0] <= Rg.min() and Rg.max() <= PC_ACCEPT[1]
    return c if ok else None


def fit_pc(pool, e, geo, rt, kE, kS):
    """Position correction of run type rt from a pool (seed, x, y) with energies e: bilinear map plus, for the hybrid form, own
    polynomials of the seed modules with >= PC_SETTING[rt][1] training clusters."""
    u, v = uv(pool, geo)
    s = pool["seed"].astype(np.int64)
    ep = e_elastic(theta_deg(geo["x"], geo["y"], geo["hz"]), geo["ebeam"], pool_mass(rt))
    r = e / ep[s]
    m = np.isfinite(u) & np.isfinite(v) & np.isfinite(r) & (r >= PC_R[0]) & (r <= PC_R[1])
    R, n, filled = fit_map(u[m], v[m], r[m], PC_FIT[rt], float(np.mean(ep[np.unique(s)])), kS)
    own = {}
    if PC_SETTING[rt][0] == "hybrid":
        for mod, idx in groups(s[m]).items():
            if idx.size >= PC_SETTING[rt][1]:
                c = fit_module(u[m][idx], v[m][idx], r[m][idx])
                if c is not None:
                    own[mod] = (c, idx.size)
    return dict(R=R, n=n, filled=filled, own=own, n_train=int(m.sum()))


def pc_eval(pc, pool, geo):
    """R per cluster: the seed module's own polynomial, else the map interpolated between bin centres."""
    u, v = uv(pool, geo)
    cu = -PC_RANGE + (np.arange(PC_BINS) + 0.5) * 2 * PC_RANGE / PC_BINS
    uc, vc = np.clip(u, cu[0], cu[-1]), np.clip(v, cu[0], cu[-1])
    iu = np.clip(np.searchsorted(cu, uc, side="right") - 1, 0, PC_BINS - 2)
    iv = np.clip(np.searchsorted(cu, vc, side="right") - 1, 0, PC_BINS - 2)
    tu, tv, G = (uc - cu[iu]) / (cu[iu + 1] - cu[iu]), (vc - cu[iv]) / (cu[iv + 1] - cu[iv]), pc["R"]
    R = ((1 - tu) * (1 - tv) * G[iu, iv] + tu * (1 - tv) * G[iu + 1, iv] + (1 - tu) * tv * G[iu, iv + 1]
         + tu * tv * G[iu + 1, iv + 1])
    for mod, (c, _) in pc["own"].items():
        m = pool["seed"] == mod
        R[m] = poly_r(c, u[m], v[m])
    return R


def write_pc(maps, rt, out, names, borrowed=None):
    """poscorr_map.tsv (bilinear maps) and poscorr_modules.tsv (own polynomials) of the switch combinations (moller, nonlinear)."""
    cu = -PC_RANGE + (np.arange(PC_BINS) + 0.5) * 2 * PC_RANGE / PC_BINS
    form = f"hybrid (own polynomials and bilinear map) of {borrowed}" if borrowed else PC_SETTING[rt][0]
    write_tsv(out / "poscorr_map.tsv", ("moller", "nonlinear", "iu", "iv", "u", "v", "n", "R", "filled"),
              [[m, n, i, j, f"{cu[i]:.3f}", f"{cu[j]:.3f}", pc["n"][i, j], f"{pc['R'][i, j]:.6f}", int(pc["filled"][i, j])]
               for (m, n), pc in maps.items() for i in range(PC_BINS) for j in range(PC_BINS)],
              f"# {rt}: form {form}, map bin fitter {PC_FIT[borrowed or rt]}; corrected E = E / R")
    write_tsv(out / "poscorr_modules.tsv", ("moller", "nonlinear", "module", "n") + tuple(f"c{k}" for k in range(8)),
              [[m, n, names[mod], cnt, *c] for (m, n), pc in maps.items() for mod, (c, cnt) in sorted(pc["own"].items())])


def read_pc(src, names):
    """Position corrections of a report directory src (poscorr_map.tsv and poscorr_modules.tsv): {(moller, nonlinear): pc}, own
    polynomials keyed by the module index in names."""
    a = np.loadtxt(src / "poscorr_map.tsv", skiprows=2, delimiter="\t", ndmin=2)     # moller, nonlinear, iu, iv, u, v, n, R, filled
    idx, maps = {str(n): i for i, n in enumerate(names)}, {}
    for m, n in {tuple(r) for r in a[:, :2].astype(int).tolist()}:
        b = a[(a[:, 0] == m) & (a[:, 1] == n)]
        R, cnt, filled = (np.zeros((PC_BINS, PC_BINS)) for _ in range(3))
        i, j = b[:, 2].astype(int), b[:, 3].astype(int)
        cnt[i, j], R[i, j], filled[i, j] = b[:, 6], b[:, 7], b[:, 8]
        maps[m, n] = dict(R=R, n=cnt.astype(int), filled=filled.astype(bool), own={}, n_train=int(cnt.sum()))
    for line in (src / "poscorr_modules.tsv").read_text().splitlines()[1:]:
        f = line.split("\t")                               # moller, nonlinear, module, n, c0..c7
        maps[int(f[0]), int(f[1])]["own"][idx[f[2]]] = (np.array(f[4:12], dtype=float), int(f[3]))
    return maps


# arms: Moller rejection, non-linearity, corrections ===========================
def switches(arm):
    """Effective (gains + position, Moller rejection, non-linearity) of an arm after the MOLLER_REJECT and NONLINEAR switches."""
    g, m, n = arm
    return g, int(m and MOLLER_REJECT), int(n and NONLINEAR)


def arm_name(a):
    """Arm label, marked when a switch turns one of its effects off."""
    return ARM_LABEL[a] + (" (switch off)" if switches(ARMS[a]) != ARMS[a] else "")


def moller_tag(pool, geo, mass):
    """Old chain's kinematic Moller tag (benchmark_hycal_timing_test_20261021.py moller_kin, lines 3297-3306): at the cluster's polar
    angle |E - E_ee| < MOLLER_KIN_N sigma and |E - E'| > MOLLER_KIN_M sigma, sigma = E' sqrt((A / sqrt(E' / GeV))^2 + B^2), E the
    corrected cluster energy."""
    th = theta_deg(pool["x"], pool["y"], geo["hz"])
    ep, ee = e_elastic(th, geo["ebeam"], mass), e_moller(th, geo["ebeam"])
    sig = ep * np.hypot(MOLLER_RES_A / np.sqrt(ep / 1000.0), MOLLER_RES_B)
    return (np.abs(pool["e"] - ee) < MOLLER_KIN_N * sig) & (np.abs(pool["e"] - ep) > MOLLER_KIN_M * sig)


def nl_view(pool, nonlinear):
    """Pool with (nonlinear) the clusterer's energy and members scaled by its factor f = e / er, the seed module's calibration
    nl1, nl2 as the old chain's --nl-members on (prad2det/src/HyCalCluster.cpp:457-473), or without it (e = er, raw members)."""
    if not nonlinear:
        return pool | {"e": pool["er"]}
    return pool | {"me": objects([w * f for w, f in zip(pool["me"], pool["e"] / pool["er"])])}


def apply_arm(pool, sw, corr, geo, mass):
    """Pool with the arm's corrected energy: gains, then position correction (corr[sw] = (gains, correction)), then Moller rejection."""
    g, m, n = sw
    pool = nl_view(pool, n)
    if g:
        gains, pc = corr[sw]
        e = np.array([gains[i] @ w for i, w in zip(pool["mi"], pool["me"])])
        pool = pool | {"e": e / pc_eval(pc, pool, geo)}
    return take(pool, ~moller_tag(pool, geo, mass)) if m else pool


def gains_job(a):
    pool, geo, mass, kE, kS = a
    g = solve_gains(pool, geo, mass, kE, kS)
    if g is None:                  # no module with an elastic window (too few clusters): gains 1
        print(f"[WARN] gain solve failed ({pool['e'].size} pool clusters), gains 1", flush=True)
        g = dict(g=np.ones(geo["x"].size), ecal=pool["e"], rows=[], g0=1.0, iterations=0, converged=False, n_free=0)
    return g


def pc_job(a):
    return fit_pc(*a)


def corrections(workers, rt, res, out_root, out, write=True):
    """Gains per run and position correction for every (Moller, non-linearity) setting of the arms, from the legacy pool of all
    splits as the old chain: gains, Moller tag with them (no position correction), gains refitted on the rejected pool,
    position correction fitted on it with these gains.  Writes gains.tsv and the poscorr tables (write; the halves of the cross
    fit are not written); returns {(1, moller, nonlinear): ({run: gains}, position correction)}."""
    geo, kE, kS, mass = res[0]["geo"], res[0]["kE"], res[0]["kS"], pool_mass(rt)
    runs = sorted({x["run"] for x in res})
    pools = {r: cat([x["pool"] for x in res if x["run"] == r]) for r in runs}
    combos = sorted({switches(a)[1:] for a in ARMS.values() if a[0]})                # (moller, nonlinear)
    solve = lambda ps: dict(zip(ps, workers.map(gains_job, [(p, geo, mass, kE, kS) for p in ps.values()])))   # noqa: E731
    train = {(0, n, r): nl_view(pools[r], n) for n in sorted({n for _, n in combos}) for r in runs}
    sol = solve(train)                                                                # gains of the unrejected pools
    if any(m for m, _ in combos):
        rej = {(1, n, r): take(p, ~moller_tag(p | {"e": sol[0, n, r]["ecal"], "x": geo["x"][p["seed"].astype(int)], "y": geo["y"][p["seed"].astype(int)]}, geo, mass)) for (_, n, r), p in train.items()}
        train |= rej
        sol |= solve(rej)
    cols = ("module", "g", "sigma_g", "n_seed", "n_member", "prior_frac")
    if write:
        write_tsv(out / "gains.tsv", ("moller", "nonlinear", "run") + cols,
                  [[m, n, r] + [x[q] for q in cols] for (m, n, r), g in sol.items() if (m, n) in combos for x in g["rows"]])
    for (m, n, r), g in sol.items():
        if write and (m, n) in combos:
            print(f"[gains] {rt} run {r} moller {m} nonlinear {n}: {len(g['rows'])} modules ({g['n_free']} free), g0 {g['g0']:.4f}, "
                  f"{g['iterations']} passes, converged {g['converged']}, {train[m, n, r]['e'].size} clusters", flush=True)
    if PC_SETTING[rt][0] == "none":      # gains only: R = 1
        flat = dict(R=np.ones((PC_BINS, PC_BINS)), n=np.zeros((PC_BINS, PC_BINS), int), filled=np.ones((PC_BINS, PC_BINS), bool), own={}, n_train=0)
        maps, borrowed = dict.fromkeys(combos, flat), None
    elif PC_SETTING[rt][0] == "hybrid":
        jobs = [(cat([train[m, n, r] for r in runs], ("seed", "x", "y")), np.concatenate([sol[m, n, r]["ecal"] for r in runs]),
                 geo, rt, kE, kS) for m, n in combos]
        maps, borrowed = dict(zip(combos, workers.map(pc_job, jobs))), None
    else:
        borrowed = PC_SETTING[rt][0]
        src = out_root / borrowed
        if not (src / "poscorr_modules.tsv").exists():
            sys.exit(f"[ERROR] {rt} borrows the position correction of {borrowed}: run its report first ({src})")
        maps = read_pc(src, geo["name"])
        if not set(combos) <= set(maps):
            sys.exit(f"[ERROR] {src} lacks the switch combinations {sorted(set(combos) - set(maps))} (switches changed since?)")
    if write:
        write_pc({c: maps[c] for c in combos}, rt, out, geo["name"], borrowed)
    for (m, n), pc in maps.items():
        if write:
            print(f"[poscorr] {rt} moller {m} nonlinear {n}{' (borrowed)' if borrowed else ''}: {pc['n_train']} training clusters, "
              f"{int((~pc['filled']).sum())}/{PC_BINS ** 2} map bins fitted, {len(pc['own'])} own modules", flush=True)
    return {(1, m, n): ({r: sol[m, n, r]["g"] for r in runs}, maps[m, n]) for m, n in combos}


# report =======================================================================
def split_job(job):
    """One split of one scan file.  corr None (first pass, legacy path only): the legacy pool with its members, input of the
    corrections.  Else corr = {switches: (gains of the run, position correction)}: the default arm for every path (variants all and
    nopu = pile-up-flagged clusters removed) and the other arms for legacy and gate ARM_GATE (variant all).  An optional fifth
    element oos (same form as corr, fitted on the other half of the splits) adds the default arm with it, columns *_oos."""
    fn, split, gates, corr = job[:4]
    oos = job[4] if len(job) > 4 else None
    z = np.load(fn, allow_pickle=False)
    run, rt = int(z["run"]), str(z["run_type"])
    cls, geo = clusterers(run, str(z["calib"]), gates)
    kE, kS, mass = round(geo["ebeam"] / E_REF, 6), round(math.sqrt(geo["ebeam"] / E_REF), 6), pool_mass(rt)
    evs, h_ev = np.flatnonzero(z["ev_split"] == split), z["h_ev"]
    lo, hi = (np.searchsorted(h_ev, [evs[0], evs[-1] + 1]) if evs.size else (0, 0))
    h = {q: z[q][lo:hi] for q in ("h_mod", "h_t", "h_e", "h_k", "h_extra")} | {"h_ev": h_ev[lo:hi]}
    t0 = time.monotonic()
    base, n_single = cluster_events(h, cls, geo, SIGNAL_MIN * kE, SEED_MIN * kE)
    t_cl = time.monotonic() - t0
    if corr is None:
        return dict(run=run, pool=base["legacy"], geo=geo, kE=kE, kS=kS)
    arm_pools = lambda arm, paths, c=corr: {p: apply_arm(base[p], switches(ARMS[arm]), c, geo, mass) for p in paths}   # noqa: E731
    arm_paths = [p for p in base if p in ("legacy", f"g{ARM_GATE:g}")]
    arms = {a: evaluate(arm_pools(a, arm_paths), geo, mass, kE, kS)[0] for a in OTHER_ARMS}
    if oos:
        arms["oos"] = evaluate(arm_pools(DEFAULT_ARM, list(base), oos), geo, mass, kE, kS)[0]   # every path: the gate choice is out-of-sample
    main = arm_pools(DEFAULT_ARM, base)
    rows = []
    for variant in ("all", "nopu"):
        pools = main if variant == "all" else {p: take(pl, ~pl["pu"].astype(bool)) for p, pl in main.items()}
        res, spec = evaluate(pools, geo, mass, kE, kS)
        for p, r in res.items():
            row = dict(run_type=rt, run=run, split=split, path=p, variant=variant, n_single=n_single[p],
                       n_pool=int(pools[p]["e"].size), **r)
            if variant == "all":
                row |= {f"{q}_{a}": res_a[p][q] for a, res_a in arms.items() if p in res_a for q in CORR_Q}
            rows.append(row)
        if variant == "all":
            spectra, ref = spec, res["legacy"]["eprime"]
            ratio = {p: np.histogram(spec[p] / ref, RATIO_EDGES)[0] for p in spec}      # E / E' of the used-module clusters
    print(f"  [split] run {run} split {split}: {evs.size} events, clustering {t_cl:.0f} s", flush=True)
    return dict(rows=rows, spectra=spectra, ratio=ratio, ref=ref, kE=kE, kS=kS)


COLS = ("run_type", "run", "split", "path", "variant", "n_single", "n_pool", "n_modules", "n_used", "eprime", "mu_fix",
        "sig_fix", "mu_g2", "sig_g2", "mu_cb", "sig_cb", "cb_mode", "res_cb", "res_g2", "res_iqr", "tail", "mu_ep", "pu",
        "clip", "n_common", "yield_n", "yield_legacy_n") + tuple(f"{q}_{a}" for a in OTHER_ARMS for q in CORR_Q) + OOS_Q
STATS = (("res_cb", "Crystal Ball sigma/E (%)"), ("res_g2", "Gaussian within 2 sigma, sigma/E (%)"),
         ("res_iqr", "interquartile range / 1.349 / median (%)"), ("mu_ep", "Crystal Ball peak / E' - 1 (%)"),
         ("tail", "outside 3 sigma of the Crystal Ball (%)"), ("pu", "pile-up flag (%)"),
         ("clip", "clipped seed on WaveAnalyzer (%)"))


def fmt(v):
    return v if isinstance(v, str) else str(v) if isinstance(v, (int, np.integer)) else f"{v:.6g}"


def write_tsv(path, header, rows, comment=None):
    """Tab-separated file: optional comment line, header, rows (floats as %.6g)."""
    with open(path, "w") as f:
        f.write((comment + "\n" if comment else "") + "\t".join(header) + "\n")
        f.writelines("\t".join(fmt(x) for x in r) + "\n" for r in rows)


def mean_se(v):
    """(mean, standard error) of the finite values."""
    v = np.array([x for x in v if math.isfinite(x)], dtype=float)
    if v.size == 0:
        return float("nan"), float("nan")
    return float(v.mean()), float(v.std(ddof=1) / math.sqrt(v.size)) if v.size > 1 else float("nan")


def label(p):
    return "legacy" if p == "legacy" else f"gate {p[1:]} ns"


def summarize(rows, paths):
    """Per path (variant all): mean +- standard error over splits of the default arm and, where present, the other arms; paired
    differences gated - legacy (d_cb, d_g2, and d_cb_<arm>), pile-up difference d_pu and the peak-yield ratio."""
    by = {}
    for r in rows:
        by.setdefault((r["path"], r["variant"]), {})[(r["run"], r["split"])] = r
    nan = float("nan")
    out = []
    for p in paths:
        a, leg, nopu = (by.get(q, {}) for q in ((p, "all"), ("legacy", "all"), (p, "nopu")))
        k = sorted(a)
        diff = lambda q: mean_se([a[x].get(q, nan) - leg[x].get(q, nan) for x in k if x in leg])   # noqa: E731
        s = dict(path=p, n_splits=len(k), n_used=sum(a[x].get("n_used", 0) for x in k),
                 n_two_peak=sum(a[x].get("cb_mode") == "two_peak" for x in k), d_cb=diff("res_cb"), d_g2=diff("res_g2"),
                 d_pu=mean_se([nopu[x].get("res_cb", nan) - a[x].get("res_cb", nan) for x in k if x in nopu]),
                 d_split={x: a[x].get("res_cb", nan) - leg[x].get("res_cb", nan) for x in k if x in leg})
        for q, _ in STATS:
            s[q] = mean_se([a[x].get(q, nan) for x in k])
        for arm in OTHER_ARMS:
            for q in CORR_Q:
                s[f"{q}_{arm}"] = mean_se([a[x].get(f"{q}_{arm}", nan) for x in k])
            s[f"d_cb_{arm}"] = diff(f"res_cb_{arm}")
        for q in OOS_Q:
            s[q] = mean_se([a[x].get(q, nan) for x in k])
        n_y, n_l = sum(a[x]["yield_n"] for x in k), sum(a[x]["yield_legacy_n"] for x in k)
        s["yield"] = 100.0 * n_y / n_l if n_l else nan
        out.append(s)
    return out


def md_table(head, rows):
    return "\n".join(["| " + " | ".join(head) + " |", "|" + "---|" * len(head)] + ["| " + " | ".join(r) + " |" for r in rows]) + "\n"


def write_tables(rows, summ, rt, runs, gates, out):
    """splits.tsv, modules.tsv, summary.tsv (default arm, every path; other arms for legacy and gate ARM_GATE) and summary.md."""
    write_tsv(out / "splits.tsv", COLS, [[r.get(c, float("nan")) for c in COLS] for r in rows])
    write_tsv(out / "modules.tsv", ("run", "split", "path", "variant", "n_modules", "modules"),
              [[r[c] for c in ("run", "split", "path", "variant", "n_modules", "modules")] for r in rows])
    keys = [q for q, _ in STATS] + ["d_cb", "d_g2", "d_pu"] + [f"{q}_{a}" for a in OTHER_ARMS for q in CORR_Q + ("d_cb",)] + list(OOS_Q)
    write_tsv(out / "summary.tsv", ["path", "n_splits", "n_used", "n_two_peak", "yield_pct"] + [f"{q}{e}" for q in keys for e in ("", "_se")],
              [[s["path"], s["n_splits"], s["n_used"], s["n_two_peak"], s["yield"]] + [x for q in keys for x in s[q]] for s in summ])
    pm = lambda m: f"{m[0]:.3f} ± {m[1]:.3f}" if math.isfinite(m[1]) else f"{m[0]:.3f}"   # noqa: E731
    shown = [s for s in summ if s["path"] in ("legacy", f"g{ARM_GATE:g}")]
    leg, gate = shown[0], shown[-1]
    val = lambda s, a, q: s[q] if a == DEFAULT_ARM else s[f"{q}_{a}"]   # noqa: E731
    on = ", ".join(n for n, f in zip(("gains and position correction", "Moller rejection", "non-linearity"), switches(ARMS[DEFAULT_ARM])) if f)
    with open(out / "summary.md", "w") as f:
        f.write(f"# {rt}: PbWO4 energy resolution, legacy versus time-gated clustering\n\nRuns {', '.join(map(str, runs))}; mean ± standard "
                f"error over {leg['n_splits']} splits; default arm = {on}.  The scan over the gates "
                f"{', '.join(f'{g:g}' for g in gates)} ns is in summary.tsv; columns are defined in README.md.\n\n"
                f"## Default arm, legacy and gate {ARM_GATE:g} ns (differences gated - legacy in percentage points)\n\n")
        oos = ["out-of-sample (even/odd cross fit), Crystal Ball sigma/E (%)"] if OOS_Q else []
        f.write(md_table(["path", "Crystal Ball sigma/E (%)"] + oos + ["Gaussian within 2 sigma (%)", "interquartile range / 1.349 (%)",
                          "Crystal Ball peak / E' - 1 (%)", "difference, Crystal Ball", "difference, Gaussian", "yield (%)",
                          "pile-up difference", "pile-up flag (%)", "clipped seed (%)", "two-peak fits"],
                         [[label(s["path"]), pm(s["res_cb"])] + [pm(s["res_cb_oos"]) for _ in oos] + [pm(s["res_g2"]), pm(s["res_iqr"]), pm(s["mu_ep"])]
                          + ["-" if s["path"] == "legacy" else pm(s[d]) for d in ("d_cb", "d_g2")]
                          + [f"{s['yield']:.2f}", pm(s["d_pu"]), pm(s["pu"]), pm(s["clip"]), f"{s['n_two_peak']}/{s['n_splits']}"]
                          for s in shown]))
        f.write("\n## Correction arms\n\nEvery arm after \"none\" has the gains and the position correction; Moller rejection and "
                "non-linearity are added as named.  Crystal Ball sigma/E (%) unless stated.\n\n")
        rows_ = [["legacy"] + [pm(val(leg, a, "res_cb")) for a in ARMS]]
        if gate is not leg:
            rows_ += [[label(gate["path"])] + [pm(val(gate, a, "res_cb")) for a in ARMS],
                      ["difference, gated - legacy (percentage points)"]
                      + [pm(gate["d_cb"] if a == DEFAULT_ARM else gate[f"d_cb_{a}"]) for a in ARMS]]
        rows_.append(["legacy, Crystal Ball peak / E' - 1 (%)"] + [pm(val(leg, a, "mu_ep")) for a in ARMS])
        f.write(md_table(["path"] + [arm_name(a) for a in ARMS], rows_))


def write_pdf(summ, spectra, rows, rt, gates, out):
    """report.pdf: page 1 legacy and every gate (resolution and peak position of the default arm and the other arms, peak yield,
    paired difference per split, pile-up flag systematic, table); page 2 the spectra and the resolution per split."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    plt.rcParams.update({
        "axes.labelsize": 16, "axes.labelweight": "bold", "axes.titlesize": 16, "axes.titleweight": "bold",
        "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 14, "axes.linewidth": 1.8,
        "xtick.major.width": 1.5, "ytick.major.width": 1.5, "xtick.major.size": 6, "ytick.major.size": 6,
        "lines.linewidth": 2.0, "figure.titlesize": 18, "figure.titleweight": "bold"})
    paths = [s["path"] for s in summ]
    x, keys, split_ticks = np.arange(len(paths)), sorted(summ[0]["d_split"]), lambda ax: ax.set_xticks(
        range(len(keys)), [f"{r}.{sp}" for r, sp in keys], rotation=60, fontsize=11)   # noqa: E731
    series = lambda ax, q, f, lab=None: ax.errorbar(x, [s[q][0] for s in summ], [s[q][1] for s in summ], fmt=f, capsize=3, label=lab)   # noqa: E731
    with PdfPages(out / "report.pdf") as pdf:
        fig, axs = plt.subplots(2, 3, figsize=(24, 14))
        for q, mk, lab in (("res_cb", "o", "two-peak Crystal Ball"), ("res_g2", "s", "Gaussian within 2 sigma"),
                           ("res_iqr", "^", "interquartile range / 1.349")):
            series(axs[0, 0], q, mk + "-", f"{lab}, {arm_name(DEFAULT_ARM)}")
        series(axs[0, 1], "mu_ep", "o-", arm_name(DEFAULT_ARM))
        if OOS_Q:
            series(axs[0, 0], "res_cb_oos", "x", f"two-peak Crystal Ball, {arm_name(DEFAULT_ARM)}, out-of-sample")
        for a in OTHER_ARMS:                                   # the other arms exist for legacy and gate ARM_GATE only
            series(axs[0, 0], f"res_cb_{a}", "D", f"two-peak Crystal Ball, {arm_name(a)}")
            series(axs[0, 1], f"mu_ep_{a}", "D", arm_name(a))
        axs[0, 2].plot(x, [s["yield"] for s in summ], "o-")
        series(axs[1, 1], "d_pu", "o-")
        for i, s in enumerate(summ):
            axs[1, 1].annotate(f"{s['pu'][0]:.2f} %", (i, s["d_pu"][0]), textcoords="offset points", xytext=(4, 6), fontsize=14)
        for s in summ[1:]:
            axs[1, 0].plot(range(len(keys)), [s["d_split"].get(k, np.nan) for k in keys], "o-", label=label(s["path"]))
        split_ticks(axs[1, 0])
        for ax, base, yl, title in ((axs[0, 0], None, "resolution sigma/E (%)", "resolution, mean ± standard error over splits"),
                                    (axs[0, 1], 0, "Crystal Ball peak / E' - 1 (%)", "calibration check"),
                                    (axs[0, 2], 100, "peak yield / legacy (%)", f"peak yield (within {YIELD_NSIGMA:g} legacy Gaussian widths of E',\ncommon modules)"),
                                    (axs[1, 0], 0, "Crystal Ball difference, gated - legacy\n(percentage points)", "paired difference per split (run.split)"),
                                    (axs[1, 1], 0, "Crystal Ball sigma/E without flagged\n- all (percentage points)",
                                     "pile-up flag systematic\n(label: flagged fraction)")):
            if base is not None:
                ax.axhline(base, color="k", lw=0.5)
            ax.set_ylabel(yl)
            ax.set_title(title)
            if ax is not axs[1, 0]:
                ax.set_xticks(x, [label(p) for p in paths])
            if ax in (axs[0, 0], axs[0, 1], axs[1, 0]):
                ax.legend(fontsize=11)
        axs[1, 2].axis("off")
        txt = ["path        sigma/E (%)      difference (percentage points)  yield (%)"]
        for s in summ:
            d = "-" if s["path"] == "legacy" else f"{s['d_cb'][0]:+.3f} ± {s['d_cb'][1]:.3f}"
            txt.append(f"{label(s['path']):<11} {s['res_cb'][0]:.3f} ± {s['res_cb'][1]:.3f}   {d:<31} {s['yield']:.2f}")
        axs[1, 2].text(0, 1, "\n".join(txt), va="top", family="monospace", fontsize=12)
        fig.suptitle(f"{rt}: legacy versus time-gated clustering, gates {', '.join(f'{g:g}' for g in gates)} ns "
                     f"({summ[0]['n_splits']} splits), {arm_name(DEFAULT_ARM)}")
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        pdf.savefig(fig)
        plt.close(fig)
        fig, axs = plt.subplots(1, 2, figsize=(24, 8.5))
        ref, kE, kS = spectra["ref"], spectra["kE"], spectra["kS"]
        if math.isfinite(ref):
            edges = np.arange(ref - FIT_HALF * kE, ref + FIT_HALF * kE + BIN * kS, BIN * kS)
            for p in paths:
                axs[0].hist(spectra[p], bins=edges, histtype="step", density=True, label=label(p))
            axs[0].axvline(ref, color="k", ls=":", lw=1)
        for p in paths:
            r = sorted((q for q in rows if q["path"] == p and q["variant"] == "all"), key=lambda q: (q["run"], q["split"]))
            axs[1].plot(range(len(r)), [q.get("res_cb", np.nan) for q in r], "o-", label=label(p))
        split_ticks(axs[1])
        axs[0].set_xlabel("cluster energy (MeV), used-module pool, all splits")
        axs[0].set_ylabel("density")
        axs[1].set_ylabel("Crystal Ball sigma/E (%)")
        for ax in axs:
            ax.legend()
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)


def report(scan_dir, gates, out_root=OUT_ROOT):
    files = sorted(Path(scan_dir).glob("*.npz"))
    if not files:
        sys.exit(f"[ERROR] no *.npz in {scan_dir}")
    jobs, run_of, types = [], {}, set()
    for fn in files:
        with np.load(fn, allow_pickle=False) as z:
            types.add(str(z["run_type"]))
            run_of[str(fn)] = int(z["run"])
            jobs += [(str(fn), int(s)) for s in z["splits"]]
    if len(types) != 1:
        sys.exit(f"[ERROR] {scan_dir} holds several run types: {sorted(types)}")
    rt, t0 = types.pop(), time.monotonic()
    out = Path(out_root) / rt
    out.mkdir(parents=True, exist_ok=True)
    paths = ["legacy"] + [f"g{g:g}" for g in gates]
    with mp.get_context("fork").Pool(NPROC) as workers:
        first = workers.map(split_job, [(fn, s, (), None) for fn, s in jobs], chunksize=1)
        corr = corrections(workers, rt, first, Path(out_root), out)
        oos = [None] * len(jobs)
        if OUT_OF_SAMPLE:       # even/odd cross fit: half = ordinal of the split within its run (split numbers are almost all even)
            by_run = {}
            for fn, s in jobs:
                by_run.setdefault(run_of[fn], []).append(s)
            half = [sorted(by_run[run_of[fn]]).index(s) % 2 for fn, s in jobs]
            fit = [corrections(workers, rt, [x for x, h in zip(first, half) if h == k], Path(out_root), out, write=False)
                   if k in half else None for k in (0, 1)]
            oos = [{a: (g[run_of[fn]], pc) for a, (g, pc) in f.items()} if f and run_of[fn] in next(iter(f.values()))[0] else None
                   for (fn, s), h in zip(jobs, half) for f in (fit[1 - h],)]
        res = workers.map(split_job, [(fn, s, tuple(gates), {a: (g[run_of[fn]], pc) for a, (g, pc) in corr.items()}, o)
                                      for (fn, s), o in zip(jobs, oos)], chunksize=1)
    rows = [r for x in res for r in x["rows"]]
    summ = summarize(rows, paths)
    spectra = {p: np.concatenate([x["spectra"][p] for x in res]) for p in paths}
    ref = [x["ref"] for x in res if math.isfinite(x["ref"])]
    spectra.update(ref=ref[0] if ref else np.nan, kE=res[0]["kE"], kS=res[0]["kS"])
    write_tables(rows, summ, rt, sorted(set(run_of.values())), gates, out)
    write_pdf(summ, spectra, rows, rt, gates, out)
    np.savez_compressed(out / "spectra.npz", edges=RATIO_EDGES, **{p: sum(x["ratio"][p] for x in res) for p in paths})
    print(f"[done] {rt}: {len(jobs)} splits, gates {', '.join(f'{g:g}' for g in gates)} ns, {time.monotonic() - t0:.0f} s -> {out}", flush=True)


# summary.pdf ==================================================================
RUN_TYPES = ("0p73GeV_H2", "2p24GeV_H2", "3p49GeV_H2", "3p49GeV_H2_old", "2p24GeV_Ta_X17")
TITLE = {"0p73GeV_H2": "0.73 GeV, hydrogen", "2p24GeV_H2": "2.24 GeV, hydrogen", "3p49GeV_H2": "3.49 GeV, hydrogen",
         "3p49GeV_H2_old": "3.49 GeV, hydrogen, previous runs", "2p24GeV_Ta_X17": "2.24 GeV, tantalum (X17 runs)"}
VARIANT = {"none": "none", "map8": "map 8 x 8", "map12": "map 12 x 12", "map16": "map 16 x 16", "map24": "map 24 x 24",
           "hybrid1000": "hybrid, own from 1000", "hybrid3000": "hybrid, own from 3000", "hybrid10000": "hybrid, own from 10000",
           "even4_global": "even4, global", "even6_global": "even6, global", "even4_pm": "even4, per module", "even6_pm": "even6, per module"}
YIELD_MIN = 99.0                       # percent of legacy: a gate must keep this peak yield to be chosen


def tsv_rows(path):
    """Rows of a tab-separated file (comment lines starting with # skipped) as dicts of strings."""
    return list(csv.DictReader((l for l in open(path) if not l.startswith("#")), delimiter="\t"))


def pc_variant(rt):
    """poscorr_comparison.tsv variant of the run type's position-correction setting."""
    form, own = PC_SETTING[rt]
    return pc_variant(form) if form in PC_SETTING else form if form == "none" else f"{form}{own}"


def run_table(out, gate_path):
    """Per run and path (legacy, gate_path): mean and standard error over the run's splits (variant all), plus the number of used
    clusters; flag = the run's out-of-sample sigma/E is more than 3 of its standard errors from the run-type mean over all splits.
    Writes <out>/runs.tsv and returns the rows."""
    keys = ("res_cb", "res_cb_oos", "mu_ep", "mu_ep_oos")
    sp = [r for r in tsv_rows(out / "splits.tsv") if r["variant"] == "all" and r["path"] in ("legacy", gate_path)]
    rows = []
    for path in ("legacy", gate_path):
        mean = mean_se([float(r["res_cb_oos"]) for r in sp if r["path"] == path])[0]
        for run in sorted({int(r["run"]) for r in sp}):
            d = [r for r in sp if r["path"] == path and int(r["run"]) == run]
            row = dict(run=run, path=path, n_splits=len(d), n_used=sum(int(float(r["n_used"])) for r in d),
                       **{k: mean_se([float(r[k]) for r in d]) for k in keys})
            row["flag"] = abs(row["res_cb_oos"][0] - mean) > 3 * row["res_cb_oos"][1]
            rows.append(row)
    write_tsv(out / "runs.tsv", ("run", "path", "n_splits") + tuple(f"{k}{e}" for k in keys for e in ("", "_se")) + ("n_used", "flag"),
              [[r["run"], r["path"], r["n_splits"]] + [x for k in keys for x in r[k]] + [r["n_used"], int(r["flag"])] for r in rows],
              "# per run: mean and standard error over the run's splits; flag = out-of-sample sigma/E more than 3 standard errors from the run-type mean")
    return rows


def write_summary(root):
    """summary.pdf: overview table with the best gate per run type, resolution progression, gate scans, position-correction
    variants, and one page per run type with the E / E' distribution of legacy and the best gate."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    plt.rcParams.update({"axes.labelsize": 16, "axes.labelweight": "bold", "axes.titlesize": 16, "axes.titleweight": "bold",
                         "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 14, "axes.linewidth": 1.8,
                         "lines.linewidth": 2.0, "figure.titlesize": 18, "figure.titleweight": "bold"})
    f = lambda d, k: float(d[k])                                                      # noqa: E731
    S = {rt: {d["path"]: d for d in tsv_rows(root / rt / "summary.tsv")} for rt in RUN_TYPES}
    gates = {rt: sorted((float(p[1:]), p) for p in S[rt] if p != "legacy") for rt in RUN_TYPES}
    best = {}
    for rt in RUN_TYPES:                                                              # smallest gate within one standard error of the minimum
        ok = [(g, p) for g, p in gates[rt] if f(S[rt][p], "yield_pct") >= YIELD_MIN and math.isfinite(f(S[rt][p], "res_cb_oos")) and math.isfinite(f(S[rt][p], "res_cb_oos_se"))]
        if not ok:
            sys.exit(f"[ERROR] {rt}: no gate keeps {YIELD_MIN:g} % of the legacy yield")
        lo = min(ok, key=lambda x: f(S[rt][x[1]], "res_cb_oos"))[1]
        best[rt] = next(x for x in ok if f(S[rt][x[1]], "res_cb_oos") <= f(S[rt][lo], "res_cb_oos") + f(S[rt][lo], "res_cb_oos_se"))
    pcr = [r for r in tsv_rows(root / "poscorr_comparison.tsv") if r["sample"] == "out"]
    pcv = lambda rt, var, path, q: next((float(r[q]) for r in pcr if (r["run_type"], r["variant"], r["path"]) == (rt, var, path)), np.nan)   # noqa: E731

    def form(rt):
        k, own = PC_SETTING[rt]
        return ("none (gains only)" if k == "none" else f"hybrid: own polynomial (>= {own} clusters) + {PC_BINS} x {PC_BINS} map"
                if k == "hybrid" else f"hybrid of {TITLE[k]} (polynomials + map)")
    sw = switches(ARMS[DEFAULT_ARM])
    cuts = (f"Cuts: one cluster above {SIGNAL_MIN:g} MeV per event (at 2.24 GeV, scaled with the beam energy); seed radius < {INNER_R:g} mm; seed energy >= "
            f"{SEED_MIN:g} MeV (scaled); seed outside rows {EXCLUDE_ROWCOL[0]}-{EXCLUDE_ROWCOL[1]} and columns {EXCLUDE_ROWCOL[2]}-{EXCLUDE_ROWCOL[3]};\n"
            f"modules with >= {MIN_PER_MODULE} events within {WINDOW_NSIGMA:g} sigma of their Crystal Ball peak and an elastic-Moller gap >= {MOLLER_GAP_NSIGMA:g} sigma; "
            f"elastic window: peak +- {WINDOW_NSIGMA:g} sigma.  Gains fitted per crystal; Moller rejection {'on' if sw[1] else 'off'}, non-linearity {'on' if sw[2] else 'off'}.")
    rule = (f"Best gate: the smallest gate whose out-of-sample Crystal Ball sigma/E is within one standard error of the minimum over the gates "
            f"with peak yield >= {YIELD_MIN:g} % of legacy.")
    q = lambda rt, p, k: f(S[rt][p], k)                                               # noqa: E731
    with PdfPages(root / "summary.pdf") as pdf:
        fig, ax = plt.subplots(figsize=(26, 10))
        ax.axis("off")
        head = ["run type", "corrections", "position-correction form", "best gate\n(ns)", "sigma/E (%),\nin-sample\nlegacy | best gate",
                "sigma/E (%),\nout-of-sample\nlegacy | best gate", "peak / E' - 1 (%),\nout-of-sample\nlegacy | best gate", "peak yield\n(%)",
                "events\nlegacy | best gate"]
        rows = [[TITLE[rt], "gains" + ("" if PC_SETTING[rt][0] == "none" else " + position"), form(rt), f"{best[rt][0]:g}",
                 f"{q(rt, 'legacy', 'res_cb'):.3f} | {q(rt, best[rt][1], 'res_cb'):.3f}",
                 f"{q(rt, 'legacy', 'res_cb_oos'):.3f} | {q(rt, best[rt][1], 'res_cb_oos'):.3f}",
                 f"{q(rt, 'legacy', 'mu_ep_oos'):+.3f} | {q(rt, best[rt][1], 'mu_ep_oos'):+.3f}", f"{q(rt, best[rt][1], 'yield_pct'):.2f}",
                 f"{q(rt, 'legacy', 'n_used'):.0f} | {q(rt, best[rt][1], 'n_used'):.0f}"] for rt in RUN_TYPES]
        t = ax.table(cellText=rows, colLabels=head, loc="upper center", cellLoc="center", colWidths=(0.14, 0.1, 0.29, 0.05, 0.08, 0.09, 0.11, 0.05, 0.09),
                     bbox=(0, 0.35, 1, 0.6))
        t.auto_set_font_size(False)
        t.set_fontsize(12)
        fig.suptitle("PbWO4 energy resolution: best configuration per run type")
        fig.text(0.02, 0.3, rule + "\n" + cuts, fontsize=12, va="top")
        pdf.savefig(fig)
        plt.close(fig)
        fig, ax = plt.subplots(figsize=(20, 9))
        steps = ("none", "gains", "gains + position", "best gate (out-of-sample)")
        for k, st in enumerate(steps):
            v = [q(rt, "legacy", "res_cb_none") if k == 0 else pcv(rt, "none", "legacy", "res_cb") if k == 1
                 else q(rt, "legacy", "res_cb_oos") if k == 2 else q(rt, best[rt][1], "res_cb_oos") for rt in RUN_TYPES]
            b = ax.bar(np.arange(len(RUN_TYPES)) + (k - 1.5) * 0.2, v, 0.19, label=st)
            ax.bar_label(b, fmt="%.2f", fontsize=12)
        ax.set_xticks(range(len(RUN_TYPES)), [TITLE[rt].replace(", ", "\n", 1) for rt in RUN_TYPES])
        ax.set_ylabel("Crystal Ball sigma/E (%)")
        ax.set_title("Resolution progression (legacy clustering up to the best gate; gains, position correction and gate out-of-sample)")
        ax.legend()
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
        fig, axs = plt.subplots(2, 3, figsize=(26, 14))
        for ax, rt in zip(axs.flat, RUN_TYPES):
            g = np.array([x[0] for x in gates[rt]])
            for k, ls in (("res_cb_oos", "o-"), ("res_cb", "s--")):
                y, e = (np.array([q(rt, p, c) for _, p in gates[rt]]) for c in (k, k + "_se"))
                ax.errorbar(g, y, e, fmt=ls, ms=4, capsize=2, label="out-of-sample" if k.endswith("oos") else "in-sample")
            ax.axhline(q(rt, "legacy", "res_cb_oos"), color="k", ls=":", label="legacy, out-of-sample")
            ax.axvline(best[rt][0], color="r", lw=1.5, label=f"chosen {best[rt][0]:g} ns")
            ax.set_title(TITLE[rt])
            ax.set_xlabel("time gate (ns)")
            ax.set_ylabel("Crystal Ball sigma/E (%)")
            ax.legend(fontsize=12)
        axs.flat[-1].axis("off")
        axs.flat[-1].text(0, 0.5, rule.replace(" whose", "\nwhose").replace(" of the minimum", "\nof the minimum"), fontsize=13)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
        per_run = {rt: run_table(root / rt, best[rt][1]) for rt in RUN_TYPES}
        fig, axs = plt.subplots(5, 2, figsize=(26, 28), gridspec_kw=dict(width_ratios=(3, 2)))
        for (ax, tab), rt in zip(axs, RUN_TYPES):
            rr = per_run[rt]
            for k, path in enumerate(("legacy", best[rt][1])):
                d = [r for r in rr if r["path"] == path]
                ax.errorbar(np.arange(len(d)) + (k - 0.5) * 0.15, [r["res_cb_oos"][0] for r in d], [r["res_cb_oos"][1] for r in d], fmt="os"[k],
                            capsize=3, color=f"C{k}", label=label(path))
                ax.plot([i + (k - 0.5) * 0.15 for i, r in enumerate(d) if r["flag"]], [r["res_cb_oos"][0] for r in d if r["flag"]], "x", color="r", ms=16, mew=3,
                        label="more than 3 standard errors from the run-type mean" if k == 0 else None)
            ax.set_xticks(range(len(d)), [str(r["run"]) for r in d])
            ax.set_xlabel("run number")
            ax.set_ylabel("out-of-sample sigma/E (%)")
            ax.set_title(TITLE[rt])
            ax.legend(fontsize=12)
            tab.axis("off")
            t = tab.table(cellText=[[str(r["run"]), r["path"] == "legacy" and "legacy" or f"gate {r['path'][1:]} ns", f"{r['res_cb_oos'][0]:.3f} ± {r['res_cb_oos'][1]:.3f}",
                                     f"{r['mu_ep_oos'][0]:+.3f}", f"{r['n_used']:,}", "yes" if r["flag"] else "no fit" if not math.isfinite(r["res_cb_oos"][0]) else ""] for r in rr],
                          colLabels=["run", "path", "sigma/E (%)", "peak / E' - 1 (%)", "events", "flagged"], loc="center")
            t.auto_set_font_size(False)
            t.set_fontsize(12)
            t.scale(1, 1.3)
        fig.suptitle("Resolution per run, out-of-sample")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
        fig, axs = plt.subplots(2, 3, figsize=(26, 14))
        for ax, rt in zip(axs.flat, RUN_TYPES):
            x = np.arange(len(VARIANT))
            for path, mk in (("legacy", "o"), (f"g{ARM_GATE:g}", "s")):
                ax.errorbar(x, [pcv(rt, v, path, "res_cb") for v in VARIANT], [pcv(rt, v, path, "res_cb_se") for v in VARIANT], fmt=mk,
                            capsize=3, label=label(path))
            ax.axvspan(list(VARIANT).index(pc_variant(rt)) - 0.4, list(VARIANT).index(pc_variant(rt)) + 0.4, color="gold", alpha=0.35, label="chosen")
            ax.set_xticks(x, list(VARIANT.values()), rotation=70, ha="right", fontsize=12)
            ax.set_ylabel("out-of-sample sigma/E (%)")
            ax.set_title(TITLE[rt])
            ax.legend(fontsize=12)
        axs.flat[-1].axis("off")
        fig.suptitle("Position-correction variants, gains on, out-of-sample")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
        for rt in RUN_TYPES:
            z = np.load(root / rt / "spectra.npz")
            c = 0.5 * (z["edges"][:-1] + z["edges"][1:])
            fig, ax = plt.subplots(figsize=(20, 11))
            inset = ax.inset_axes((0.56, 0.36, 0.42, 0.6))
            for p, col in (("legacy", "C0"), (best[rt][1], "C1")):
                y, name = z[p], label(p)
                for a in (ax, inset):
                    a.step(c, y, where="mid", color=col, label=f"{name}: {y.sum():,} clusters" if a is ax else None)
                m0 = c[np.argmax(y)]
                w = (c > m0 - 0.18) & (c < m0 + 0.09)
                fit, _ = fit_curve(crystal_ball, c[w], y[w], [y.max(), m0, 0.03, 1.5, 2.0], [0, 0.6, 1e-3, 0.1, 1.01], [np.inf, 1.4, 0.18, 20.0, 50.0], maxfev=10000)
                if fit is not None:
                    xs = np.linspace(0.6, 1.4, 800)
                    for a in (ax, inset):
                        a.plot(xs, crystal_ball(xs, *fit), "--", color=col, lw=1.5, label=f"Crystal Ball fit, sigma/E {100 * abs(fit[2]) / fit[1]:.2f} %" if a is ax else None)
                    if p == "legacy":
                        for a in (ax, inset):
                            a.axvspan(fit[1] - WINDOW_NSIGMA * abs(fit[2]), fit[1] + WINDOW_NSIGMA * abs(fit[2]), color="gray", alpha=0.2, label="elastic window" if a is ax else None)
            ax.set_yscale("log")
            ax.set_ylim(bottom=1)
            inset.set_xlim(0.85, 1.12)
            ax.set_xlabel("E / E' (selected clusters, best configuration)")
            ax.set_ylabel("clusters per bin (log scale)")
            ax.set_title(f"{TITLE[rt]}: energy distribution, legacy and gate {best[rt][0]:g} ns")
            ax.legend(loc="lower left")
            fig.text(0.02, 0.01, "Curves: Crystal Ball fit of the pooled histogram (the table quotes the mean over splits of the two-peak fit); inset: linear scale.", fontsize=12)
            fig.tight_layout(rect=(0, 0.03, 1, 1))
            pdf.savefig(fig)
            plt.close(fig)
    print(f"[done] {root / 'summary.pdf'}", flush=True)


def main():
    a = sys.argv[1:]
    if len(a) == 2 and a[0] == "summary":
        write_summary(Path(a[1]))
    elif len(a) >= 3 and a[0] == "scan":
        scan(a[1], a[2:])
    elif len(a) >= 3 and a[0] == "report":
        out_root = OUT_ROOT
        if a[1] == "--out":
            out_root, a = Path(a[2]), a[2:]
        gates = sorted({float(g) for g in a[2:]})
        if not gates or gates[0] <= 0:
            sys.exit("[ERROR] gates must be > 0 ns")
        report(a[1], gates, out_root)
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main()
