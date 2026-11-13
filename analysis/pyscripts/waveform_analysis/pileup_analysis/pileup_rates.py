#!/usr/bin/env python3
"""pileup_rates.py -- pile-up rates in the data, channel level and cluster level (one leading cluster per event).

    python3 waveform_analysis/pileup_analysis/pileup_rates.py scan <out.npz> <evio file> [<evio file> ...]
    python3 waveform_analysis/pileup_analysis/pileup_rates.py report <scan dir> <out dir>

scan    all EVIO files (splits) of one 0p73GeV_H2, 2p24GeV_H2, 3p49GeV_H2 or 2p24GeV_Ta_X17 run (other runs
        are skipped).  Events: trigger gate of the run type and LMS-bit veto only (no multiplicity
        or pulse-quality cut).  Every HyCal channel (PbWO4 and lead glass) with 100 samples is first checked by bad_waveform()
        (any sample equal to 0 or more than UNDERSHOOT_ADC below the baseline): flagged channels are removed and only counted
        per module (bad_count).  Every other channel whose
        raw max - median of the first 30 samples is >= PRECUT_ADC is fitted with the matched filter
        of pileup_deconvolution_study.py (class MF copied below; that docstring's Algorithm describes its
        three stages), using the run type's own spline bank
        (lead glass and modules without a valid row: the bank's average spline).  A channel is kept if its largest fitted amplitude is
        >= KEEP_ADC (10 ADC); k = number of fitted pulses, pile-up = k >= 2.  A kept channel whose main pulse (largest amplitude) has its
        fit window cut by the end of the waveform (trunc_main) stays in the raw k counts and k maps but is excluded from every
        E_main-based quantity (cut rates, slices, maps, clusters, cut extras, examples); it is counted per module (trunc_count).
        Saves the kept channels (with the fitted baseline offset dped and trunc_main),
        per-event counts, and two random reservoirs (k >= 2 and k == 1, A_main >= A_EX) of raw
        waveforms with their fits, and per module (up to EX_PER_MODULE per run) random pile-up (cut) channels and
        (up to EX_CLEAN_PER_MODULE) one-pulse channels above the E_main threshold, after the vetoes, for the hot-module and ring examples.
report  pool the <scan dir>/run_*.npz scans per run type -> pileup_rates.tsv, pileup_grid_<RT>.tsv,
        pileup_modules_<RT>.tsv, pileup_rate_map_<RT>.png, pileup_rate_slices_<RT>.png, pileup_slices.tsv,
        summary.txt and pileup_rates.pdf in <out dir>, plus examples.pdf (pile-up (cut) examples per run type; the X17 beam-hole
        ring first, then 16 spread examples per run type on 2 pages, the excluded-channels map page (bad waveforms, saturated and bad fits per module) and the N(k >= n) / N(kept channels after the vetoes) map pages for n = 2 ... 6, k >= 2 first, 5 pages)
        and hot_modules.pdf (the HOT_N modules with the highest cut rate per run type: summary table, then one page each).

Channel rate = N(k >= 2) / N(kept, A_main >= KEEP_ADC) +- binomial error.  Grid (k == 2): main = earlier pulse,
dt = t_later - t_earlier, ratio = A_later / A_earlier; cell = N(pairs in cell) / N(kept) in %.
Per-module k columns of pileup_modules_<RT>.tsv: rate_k2_raw = N(k >= 2) / N(all kept channels, before the vetoes); rate_k2 ... rate_k6 = N(k >= n) / N(kept channels after the vetoes and exclusions), the denominator of the k maps.
Cluster level (one leading cluster per event): in-window pulses have onset in the pre-window (130, 200) ns
(benchmark_hycal_timing_test_event_hits; time offsets off); seed = kept channel with the largest in-window pulse;
members = seed + modules within 1.5 x the mean module size of it with an in-window pulse; per member the used pulse
is the largest in-window one, the others are extras (in-window or not; r = A_extra / A_used).
Cut rates: a channel is pile-up (cut) when E_main = A_main x gain >= CUT_EFRAC (0.1) x E_beam of its run and another
pulse has A_other / A_main > CUT_R (0.1) and |t_other - t_main| >= CUT_SEP_NS (8 ns).  gain (MeV/ADC) = calibration factor of
the run (runinfo calibration file) x I_int of the module's spline template; modules without a factor are left out (counted in
summary.txt).  Channel cut rate = N(cut) / N(kept, E_main >= thr), per module likewise (modules with >= MIN_N_MOD such channels);
cluster cut rate = N(clusters with a member at E_used >= thr and an extra with r > 0.1) / N(clusters with seed E_used >= thr).
Slices: the cut rate in 6 slices of E_main / E_beam, [0.1, 0.3) ... [1.1, 1.3): N(cut channels in the slice) / N(all kept channels of the module after the vetoes and exclusions) per module (drawn for modules with >= 10 channels in the slice), pooled over the run type with N(all kept channels of the run type) (pileup_rate_slices_<RT>.png, pileup_slices.tsv).
Veto: a channel that is the only kept channel in its event is dropped from the cut rates, maps, slices, medians, dt histogram and
clustering (self-triggered bursts in a few modules); channels of an event with >= VETO_NKEPT (300) kept channels (array-wide pulses) are vetoed the same way;
raw k >= 2, grid and clipped rates keep both; the vetoed counts are reported.
Median dt (t_other - t_main, signed) and r of the cut channels' extra with the largest A_other, and of the largest-r extra that
triggered the cut in each cut cluster (all / in-window extras), are reported per run type.
Clipped rate (clipped(): any raw sample >= 4095) = N(clipped kept channels) / N(kept), cluster level
N(clusters with a clipped member) / N(clusters), and the k >= 2 fraction among clipped channels.
Per run (database/runinfo/general.json, run_config): beam energy (E_beam of the 0.1 E_beam threshold; stored in the scan as beam_energy_mev) and the
hycal_dead_modules list, whose modules are excluded from the scan (no channels, not in n_kept or clusters; stored as dead_modules); the report marks
them (dead column, hatched in the maps, count in summary.txt).
Standalone: numpy, scipy, matplotlib and the compiled prad2py decoder only.  Database files (DAQ config, hycal_map.json,
runinfo/general.json, calibration/*.json) are read from $PRAD2_DATABASE_DIR; source setup_env.sh first (it exports that
and puts prad2py on PYTHONPATH).
chi2/dof as in pileup_deconvolution_study.py (sigma(t) variance), with the module's per-file noise reference sigma_ref as noise sigma.
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")                     # one BLAS thread per process

import json
import re
import sys
import textwrap
import time
import warnings
from pathlib import Path

import matplotlib
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks

# plot style: large bold axis labels, large legends
matplotlib.rcParams.update({
    "axes.labelsize": 16, "axes.labelweight": "bold", "axes.titlesize": 16, "axes.titleweight": "bold",
    "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 14, "axes.linewidth": 1.8,
    "xtick.major.width": 1.5, "ytick.major.width": 1.5, "xtick.major.size": 6, "ytick.major.size": 6,
    "lines.linewidth": 2.0, "figure.titlesize": 18, "figure.titleweight": "bold"})

HERE = Path(__file__).resolve().parent
PYDIR = HERE.parent.parent                              # analysis/pyscripts
BANK_ROOT = PYDIR / "waveform_analysis" / "pulse_templates"
DB_DIR = Path(os.environ.get("PRAD2_DATABASE_DIR", "")) # database: always $PRAD2_DATABASE_DIR (setup_env.sh exports it)
HYCAL_MAP = DB_DIR / "hycal_map.json"
RUN_TYPE_OF = {
    24590: "0p73GeV_H2", 24603: "0p73GeV_H2",
    25228: "2p24GeV_H2", 25257: "2p24GeV_H2", 25260: "2p24GeV_H2", 25299: "2p24GeV_H2", 25320: "2p24GeV_H2",
    24698: "3p49GeV_H2", 24911: "3p49GeV_H2", 24937: "3p49GeV_H2",
    # previous 3p49GeV_H2 runs (replaced 2026-11-10)
    # 24386: "3p49GeV_H2", 24503: "3p49GeV_H2", 24512: "3p49GeV_H2",
    25800: "2p24GeV_Ta_X17", 26138: "2p24GeV_Ta_X17", 26230: "2p24GeV_Ta_X17",
}
RUN_TYPES = ["0p73GeV_H2", "2p24GeV_H2", "3p49GeV_H2", "2p24GeV_Ta_X17"]
PHYSICS_ANY = {"2p24GeV_Ta_X17"}       # trigger: (bits & 0xf00) != 0; others bits == 0x100
PHYSICS_BITS_EXACT, PHYSICS_MASK_ANY = 0x100, 0xf00
LMS_BIT = 0x1000000                    # events with this trigger bit are vetoed
VETO_NKEPT = 300                       # events with >= this many kept channels are array-wide pulses, not physics; vetoed like lone channels
NS_EXPECT, OVERFLOW, BLOCK = 100, 4095, 256
MAX_EVENTS = int(os.environ.get("PILEUP_MAX_EVENTS", "0"))   # test aid only: stop the scan after this many decoded events of the run (0 = all)
KEEP_ADC = 10.0   # kept channel: fitted main amplitude (baseline-subtracted) >= this; ≈ 9.5 MeV at the median gain; Møller partners and shower neighbours
PRECUT_ADC = 8.0   # matched filter only on channels with raw max - baseline >= this
UNDERSHOOT_ADC = 20.0   # a sample this far below the baseline, or any sample equal to 0, marks a bad waveform
SAT_NCLIP = 3     # a clipped main pulse is extrapolated by the fit (e.g. 12 000 ADC) and its wavy tail becomes a spurious multi-pulse fit
CHI2_MAX = 1000.0 # chi2/dof above this: distorted waveforms; cross-talk spikes (5–20) and unresolved pile-up (10–100s) stay
CHI_EDGES = np.logspace(-1, 5, 61)   # chi2/dof histogram bins (10 per decade); values outside are put in the end bins
REC = dict(midx=np.int16, ev=np.int32, k=np.int8, a=np.float32, t=np.float32, chi2=np.float32,
           ped=np.float32, sig=np.float32, nclip=np.int8, dped=np.float32, trunc_main=bool)          # per kept channel (ped includes dped)
A_EX, N_RES = 200.0, 60                # example reservoirs: A_main >= A_EX, N_RES per class and run
DT_EDGES = [8, 12, 16, 20, 24, 32, 48, 64, 96, 128, np.inf]
R_EDGES = [0, 0.05, 0.1, 0.2, 0.3, 0.5, 1, 2, np.inf]
MAP_LIM, MIN_N_MOD = 350.0, 100        # heat map range (mm); per-module cut and slice rates only for modules with >= MIN_N_MOD channels (raw k >= 2: > MIN_N_MOD kept)
PRE_LO, PRE_HI, CL_RADIUS = 130.0, 200.0, 1.5   # pre-window (ns, benchmark default; time offsets off); cluster radius / mean size
SLICE_MIN_N = 10   # per-module slice rate drawn only for modules with >= this many channels in the slice
SLICE_EDGES = np.round(np.arange(0.1, 1.31, 0.2), 1)   # E_main / E_beam slices of the cut rate: 6 slices 0.1-0.3 ... 1.1-1.3
CUT_EFRAC, CUT_R, CUT_SEP_NS = 0.1, 0.1, 8.0     # report: pile-up (cut) = E_main >= 0.1 E_beam, an other pulse with A_other / A_main > 0.1 and |dt| >= 8 ns
RUNINFO = DB_DIR / "runinfo" / "general.json"            # per-run beam energy, calibration file, dead modules
N_EX_PU, N_EX_1 = 20, 10
EX_PER_MODULE, EX_CLEAN_PER_MODULE = 3, 1   # per module and run: pile-up (cut) and one-pulse channels kept for examples.pdf / hot_modules.pdf
RING_RUN_TYPE, N_RING, N_RING_EX = "2p24GeV_Ta_X17", 12, 2   # beam-hole ring: the N_RING PbWO4 modules closest to the hole, N_RING_EX examples each
N_SPREAD, HOT_N, HOT_PU, HOT_CLEAN = 16, 8, 6, 2           # spread examples per run type; hot modules per run type, pile-up and one-pulse panels each

# --- matched filter: copied from pileup_deconvolution_study.py (same class MF and constants)
PED_NSAMPLES = 30
KMAX = 6              # at most this many pulses per waveform
NSIGMA, MIN_AMP = 5.0, 5.0   # a pulse needs A >= max(NSIGMA sigma_A, MIN_AMP) (stage A: prominence >= that)
MIN_SEP_NS = 8.0      # smallest onset separation of two pulses
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
DT, NPH, L, LW, LEAD = 4.0, 4, 24, 32, 3
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


def fit(mf, Y, m, sig, b_ref):
    """pileup_deconvolution_study.fit on channels Y (LUT rows m), noise sigma = sig [n] (the module's
    sigma_ref), baseline fit start b_ref [n] (run_ref); adds far (start replaced by b_ref).
    Pulses are returned in time order: t (NaN pad), a (0 pad), comp = a_k * shape_k [n, KMAX, ns]."""
    ped, sig, t0, a, v, W, dped, trunc, far = run_ref(mf, Y, m, OVERFLOW, sig, b_ref)
    ns = Y.shape[1]
    S = np.stack([mf.basis(m, t0[:, k].astype(F32), ns) for k in range(KMAX)], axis=1)
    model = ped[:, None] + np.einsum("nk,nks->ns", np.where(v, a, 0.0), S)
    mask = MF.windows(t0, v, ns) & (W > 0)
    noise = np.sqrt(sig.astype(np.float64) ** 2 + 1.0 / 12.0)          # + rounding variance
    var = mf.variance(m, t0.astype(F32), a, v, noise, ns)
    dof = mask.sum(1) - 2 * v.sum(1) - 1                  # 2 per pulse + the baseline offset
    with np.errstate(invalid="ignore", divide="ignore"):
        chi2 = np.where(v.any(1) & (dof > 0), (mask * (Y - model) ** 2 / var).sum(1) / dof, np.nan)
    o = np.argsort(np.where(v, t0, np.inf), axis=1)
    tk = lambda x: np.take_along_axis(x, o, axis=1)                   # noqa: E731
    vs = tk(v)
    comp = np.take_along_axis(np.where(v, a, 0.0)[:, :, None] * S, o[:, :, None], axis=1)
    return dict(ped=ped, sig=sig, k=v.sum(1), t=np.where(vs, tk(t0), np.nan), a=np.where(vs, tk(a), 0.0),
                comp=comp, chi2=chi2, nclip=clipped(Y), dped=dped, far=far,
                trunc_main=trunc[np.arange(Y.shape[0]), np.where(v, a, -1.0).argmax(1)])      # main pulse's fit window cut by the end of the waveform


def main_other(k, a, t):
    """Per channel: main pulse (largest amplitude), its amplitude, the other pulses with A_other / A_main > CUT_R and |dt| >= CUT_SEP_NS
    (mask [n, KMAX]), and the index of the largest of them."""
    n = k.size
    valid = np.arange(KMAX) < k[:, None]
    main = np.where(valid, a, -1.0).argmax(1)
    amain = a[np.arange(n), main].astype(np.float64)
    oth = valid & (np.arange(KMAX) != main[:, None]) & (a > CUT_R * amain[:, None]) \
        & (np.abs(t - t[np.arange(n), main][:, None]) >= CUT_SEP_NS)
    return main, amain, oth, np.where(oth, a, -1.0).argmax(1)


def reservoir(old, new, per):
    """Uniform random sample of up to `per` entries per module: the entries with the smallest random keys of old + new."""
    if old is not None:
        new = {k: np.concatenate([old[k], new[k]]) for k in new}
    if new["key"].size == 0:
        return old
    o = np.lexsort((new["key"], new["midx"]))
    m = new["midx"][o]
    i = np.arange(m.size)
    rank = i - np.maximum.accumulate(np.where(np.r_[True, m[1:] != m[:-1]], i, 0))
    return {k: v[o[rank < per]] for k, v in new.items()}


# scan =========================================================================
def scan(out, evios):
    from prad2py import dec, det                         # compiled decoder binding (PYTHONPATH from setup_env.sh)
    t_start = time.monotonic()
    parsed = [re.search(r"prad_0*(\d+)\.evio\.(\d+)", os.path.basename(f)) for f in evios]
    if not all(parsed):
        sys.exit("[ERROR] file names must be prad_0<run>.evio.<split>")
    runs = {int(m.group(1)) for m in parsed}
    if len(runs) != 1:
        sys.exit(f"[ERROR] files of several runs: {sorted(runs)}")
    run = runs.pop()
    splits = np.array([int(m.group(2)) for m in parsed], dtype=np.int64)
    rt = RUN_TYPE_OF.get(run)
    if rt is None:
        print(f"[skip] run {run} is not one of {RUN_TYPES}")
        return
    cfg = run_config(run)
    dead = set(cfg["dead"])                              # dead modules: no channels, no clusters
    if rt in PHYSICS_ANY:
        gate = lambda b: (b & PHYSICS_MASK_ANY) != 0          # noqa: E731
    else:
        gate = lambda b: b == PHYSICS_BITS_EXACT              # noqa: E731
    b = det.PipelineBuilder()                            # DAQ config + HyCal map of the run from $PRAD2_DATABASE_DIR
    b.set_run_number(run)
    core = b.build()
    hycal = core.hycal
    crate_map = {int(r.tag): int(r.crate) for r in core.daq_cfg.roc_tags}
    mods = [hycal.module(i) for i in range(hycal.module_count())]
    names = np.array([m.name for m in mods])
    pwo4 = np.array([bool(m.is_pwo4()) for m in mods])
    if not dead <= set(names):
        sys.exit(f"[ERROR] dead modules not in the HyCal map: {sorted(dead - set(names))}")
    bank = Bank(BANK_ROOT / rt / f"templates_{rt}.npz", names, pwo4)
    mf = MF(bank)
    rng = np.random.default_rng(run)
    rng_mod = np.random.default_rng(run + 1)                   # keys of the per-module reservoirs (separate: the other reservoirs stay unchanged)
    gain = module_gain(cfg, rt, names)                         # MeV per ADC; nan: no calibration
    athr = np.where(np.isfinite(gain), CUT_EFRAC * cfg["beam"] / gain, np.inf)   # ADC threshold of E_main >= CUT_EFRAC E_beam
    print(f"[setup] run {run} {len(evios)} splits {rt}; {len(mods)} modules; beam {cfg['beam']:g} MeV; dead {sorted(dead)}", flush=True)

    rec = {k: [] for k in REC}
    bad_count, trunc_count, digit_count = (np.zeros(len(mods), np.int64) for _ in range(3))   # per module: bad waveforms, kept channels with a truncated main pulse, digitised channels
    ex = {0: None, 1: None}                                    # reservoirs: 1 = k >= 2, 0 = k == 1
    pm = {0: None, 1: None}                                    # per-module reservoirs after the vetoes: 1 = pile-up (cut), 0 = k == 1 above threshold
    evnum, rows, evs, mis = [], [], [], []

    sig_ref = b_ref = None                                   # noise reference of the file being read
    ref_n, ref_b, ref_sig = [], [], []                       # per file: reference waveforms, b_ref, sigma_ref per module
    n_base_ref = np.zeros(len(mods), np.int64)                 # per module: waveforms whose baseline fit started at b_ref

    def flush():
        if not rows:
            return
        Y, ev, mi = np.stack(rows), np.array(evs), np.array(mis)
        rows.clear(), evs.clear(), mis.clear()
        bad = bad_waveform(Y)                                              # flagged before any cut: removed, only counted
        digit_count[:] += np.bincount(mi, minlength=len(mods))
        bad_count[:] += np.bincount(mi[bad], minlength=len(mods))
        Y, ev, mi = Y[~bad], ev[~bad], mi[~bad]
        x = np.median(Y[:, :PED_NSAMPLES], axis=1)
        s = np.nonzero(Y.max(axis=1) - x >= PRECUT_ADC)[0]
        f = fit(mf, Y[s], bank.row[mi[s]], sig_ref[mi[s]], b_ref[mi[s]])
        n_base_ref[:] += np.bincount(mi[s][f.pop("far")], minlength=len(mods))
        keep = f["a"].max(axis=1) >= KEEP_ADC
        s = s[keep]
        f = {k: v[keep] for k, v in f.items()}
        trunc_count[:] += np.bincount(mi[s][f["trunc_main"]], minlength=len(mods))
        for k, v in dict(midx=mi[s], ev=ev[s], **{k: f[k] for k in rec if k in f}).items():
            rec[k].append(v.astype(REC[k]))
        big = (f["a"].max(axis=1) >= A_EX) & ~f["trunc_main"]
        for c, sel in ((1, big & (f["k"] >= 2)), (0, big & (f["k"] == 1))):
            new = dict(samples=Y[s[sel]], midx=mi[s[sel]], evnum=np.array(evnum)[ev[s[sel]]],
                       key=rng.random(int(sel.sum())), **{k: v[sel] for k, v in f.items()})
            if ex[c] is not None:
                new = {k: np.concatenate([ex[c][k], new[k]]) for k in new}
            o = np.argsort(new["key"])[:N_RES]
            ex[c] = {k: v[o] for k, v in new.items()}
        _, inv, cnt = np.unique(ev[s], return_inverse=True, return_counts=True)
        nk = cnt[inv]                                                       # kept channels in the event (the block holds whole events)
        main, amain, oth, _ = main_other(f["k"], f["a"].astype(np.float32), f["t"].astype(np.float32))
        big = (amain >= athr[mi[s]]) & (nk > 1) & (nk < VETO_NKEPT) & ~f["trunc_main"]
        for c, sel in ((1, big & oth.any(1)), (0, big & (f["k"] == 1))):
            new = dict(samples=Y[s[sel]], midx=mi[s[sel]], evnum=np.array(evnum)[ev[s[sel]]], key=rng_mod.random(int(sel.sum())),
                       emain=(amain[sel] * gain[mi[s[sel]]]).astype(np.float32), k=f["k"][sel], **{q: f[q][sel].astype(np.float32) for q in ("ped", "sig", "t", "a", "chi2")})
            pm[c] = reservoir(pm[c], new, EX_PER_MODULE if c else EX_CLEAN_PER_MODULE)

    cache = {}
    def module_of(tag, crate, sl, c):
        """HyCal module index of a DAQ channel, -1 for channels that are not scanned (cached)."""
        mi = cache.get((tag, sl, c))
        if mi is None:
            m = hycal.module_by_daq(crate, sl, c)
            mi = cache[(tag, sl, c)] = m.index if (m is not None and m.is_hycal() and m.name not in dead) else -1
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
            return np.zeros(len(mods), np.int64), np.full(len(mods), np.nan, F32), np.full(len(mods), SIGMA_MIN, F32)
        return noise_reference(np.concatenate(X), np.concatenate(M), len(mods))

    n_seen = n_trig = n_lms = 0
    done = False
    for evio in evios:
        if done:
            break
        flush()                                                # blocks never span two files: each has its own noise reference
        n_ref, b_ref, sig_ref = read_reference(evio, np.ones(len(mods), bool))
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
    r = {k: (np.concatenate(v) if v else np.zeros((0, KMAX) if k in ("a", "t") else 0, REC[k])) for k, v in rec.items()}
    n_pass = len(evnum)
    pu = r["k"] >= 2
    ev_out = dict(ev_number=np.array(evnum, dtype=np.int64),
                  ev_n_kept=np.bincount(r["ev"], minlength=n_pass))
    wall = time.monotonic() - t_start
    exo = {f"ex{c}_{k}": v for c in (0, 1) if ex[c] is not None for k, v in ex[c].items()}
    exo.update({f"pm{c}_{k}": v for c in (0, 1) if pm[c] is not None for k, v in pm[c].items()})
    sr = np.array(ref_sig)[np.array(ref_n) >= MIN_REF_WAVEFORMS]
    sr = sr if sr.size else np.array([SIGMA_MIN])
    print(f"[noise reference] sigma_ref median {np.median(sr):.3f} ADC (16-84 %: {np.percentile(sr, 16):.3f}-{np.percentile(sr, 84):.3f}); "
          f"{int((np.array(ref_n) < MIN_REF_WAVEFORMS).sum())} module-files on the median sigma_ref; {int(n_base_ref.sum())} baseline starts from b_ref", flush=True)
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    np.savez_compressed(out, run=run, splits=splits, run_type=rt, n_seen=n_seen, n_trig=n_trig, n_lms=n_lms,
                        n_pass=n_pass, wall_s=wall, beam_energy_mev=cfg["beam"], dead_file=cfg["dead_src"],
                        mod_n_ref=np.array(ref_n), mod_b_ref=np.array(ref_b), mod_sig_ref=np.array(ref_sig), mod_n_base_ref=n_base_ref, dead_modules=np.array(sorted(dead), dtype=str), bad_count=bad_count, trunc_count=trunc_count, digit_count=digit_count, mod_name=names, mod_pwo4=pwo4,
                        **{f"ch_{k}": v for k, v in r.items()}, **ev_out, **exo)
    print(f"[done] {out}: {n_seen} events, {n_pass} after trigger + LMS, {r['k'].size} channels kept, "
          f"{int(pu.sum())} with k >= 2; {wall:.0f} s", flush=True)


# report =======================================================================
def binom(k, n):
    r = k / max(n, 1)
    return r, np.sqrt(r * (1 - r) / max(n, 1))


def neighbours(names):
    """Boolean matrix [seed, module]: centre within CL_RADIUS x the mean of the two module sizes."""
    geo = {m["n"]: m["geo"] for m in json.load(open(HYCAL_MAP))}
    x, y, sz = (np.array([geo[n][q] for n in names]) for q in ("x", "y", "sx"))
    dist = np.hypot(x[:, None] - x, y[:, None] - y)
    return dist <= CL_RADIUS * 0.5 * (sz[:, None] + sz)


def clusters(z, nbr, athr):
    """Leading cluster of each event of one scan file (see the module docstring); athr = per-module ADC threshold of the E_main cut; returns the counts."""
    k, a, t, ev, mi = z["ch_k"], z["ch_a"].astype(np.float64), z["ch_t"], z["ch_ev"], z["ch_midx"]
    valid = np.arange(KMAX) < k[:, None]
    win = valid & (t > PRE_LO) & (t < PRE_HI)
    used = np.where(win, a, -1.0).argmax(1)                    # largest in-window pulse
    a_used = np.where(win.any(1), a[np.arange(k.size), used], 0.0)
    el = np.nonzero(win.any(1))[0]                             # channels with an in-window pulse
    o = el[np.lexsort((-a_used[el], ev[el]))]                  # per event, largest first
    seed = o[np.r_[True, ev[o][1:] != ev[o][:-1]]]
    seed_of = np.full(z["n_pass"], -1)
    seed_of[ev[seed]] = mi[seed]
    mem = el[nbr[seed_of[ev[el]], mi[el]]]                     # members of the leading cluster
    extra = valid[mem] & (np.arange(KMAX) != used[mem][:, None])
    r = a[mem] / a_used[mem][:, None]
    hi = (a_used >= athr[mi]) & ~z["ch_trunc_main"]            # E_used >= threshold, main pulse not truncated
    cut = (extra & (r > CUT_R) & (np.abs(t[mem] - t[mem, used[mem]][:, None]) >= CUT_SEP_NS)).any(1) & hi[mem]
    fl = {"clip": z["ch_nclip"][mem] > 0, "cut": cut, "cutin": cut & (extra & win[mem] & (r > CUT_R)).any(1)}
    out = dict(n_ev_kept=np.unique(ev).size, n_cl=seed.size, n_cl_emain=int(hi[seed].sum()))
    trig = extra & (r > CUT_R) & (np.abs(t[mem] - t[mem, used[mem]][:, None]) >= CUT_SEP_NS) & hi[mem][:, None]
    for name, m in (("all", trig), ("in", trig & win[mem])):   # per cluster: the extra with the largest r that triggered the cut
        rr = np.where(m, r, -1.0)
        j = rr.argmax(1)
        h = np.nonzero(rr[np.arange(mem.size), j] > 0)[0]
        h = h[np.lexsort((-rr[h, j[h]], ev[mem][h]))]
        h = h[np.r_[True, ev[mem][h][1:] != ev[mem][h][:-1]]] if h.size else h
        out[f"dt_{name}"] = t[mem[h], j[h]] - t[mem[h], used[mem[h]]]
        out[f"r_{name}"] = rr[h, j[h]]
    for f, v in fl.items():
        out[f] = (np.bincount(ev[mem], weights=v, minlength=z["n_pass"])[ev[seed]] > 0).sum()
    return out


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


def load(scan_dir):
    """Run scans pooled per run type, reduced per run to what the report uses: k, midx, the chi2
    of k == 1 channels with A <= 1000 ADC, (dt, ratio) of k == 2 channels, event counts, reservoirs."""
    files = sorted(Path(scan_dir).glob("run_*.npz"))
    if not files:
        sys.exit(f"[ERROR] no run_*.npz in {scan_dir}")
    D, nbr, G = {}, None, {}
    for fn in files:
        z = np.load(fn, allow_pickle=False)
        rt, run = str(z["run_type"]), int(z["run"])
        d = D.setdefault(rt, dict(splits=0, runs=set(), wall=[], n_pass=0, n_seen=0, n_ev_arr=0, n_lms=0, cl=dict(
                                      n_ev_kept=0, n_cl=0, dt_all=[], r_all=[], dt_in=[], r_in=[],
                                      clip=0, cut=0, cutin=0, n_cl_emain=0), bad=0, trunc=0, digit=0, sat=0, bfit=0, n_pre=0, chi_h=np.zeros((2, len(CHI_EDGES) - 1), np.int64), ch={}, ex={0: {}, 1: {}}, pm={0: {}, 1: {}}, beam={}, dead=set()))
        d["splits"] += len(z["splits"])
        d["runs"].add(run)
        d["beam"][run] = float(z["beam_energy_mev"])
        d["dead"] |= set(map(str, z["dead_modules"]))
        d["wall"].append(float(z["wall_s"]))
        for q in ("bad", "trunc", "digit"):
            d[q] = d[q] + z[f"{q}_count"]
        for k in ("n_pass", "n_seen", "n_lms"):
            d[k] += int(z[k])
        names, pwo4 = z["mod_name"], z["mod_pwo4"]
        nbr = neighbours(names) if nbr is None else nbr
        chi = z["ch_chi2"]
        sat, bfit = z["ch_nclip"] >= SAT_NCLIP, chi > CHI2_MAX             # report-level exclusions (nan chi2 passes): only counted per module
        for key, f in (("sat", sat), ("bfit", bfit)):
            d[key] = d[key] + np.bincount(z["ch_midx"][f], minlength=names.size)
        d["n_pre"] += chi.size
        for i, sel in enumerate((z["ch_k"] == 1, z["ch_k"] >= 2)):           # chi2/dof of all kept channels, before the exclusions
            d["chi_h"][i] += np.histogram(np.clip(chi[sel & np.isfinite(chi)], CHI_EDGES[0], CHI_EDGES[-1]), CHI_EDGES)[0]
        zc = {q: z[q][~(sat | bfit)] for q in z.files if q.startswith("ch_")}
        if run not in G:
            G[run] = module_gain(run_config(run), rt, names)
        gain = G[run]
        athr = CUT_EFRAC * d["beam"][run] / gain                 # ADC threshold of E_main >= CUT_EFRAC E_beam; nan: no calibration
        athr = np.where(np.isfinite(athr), athr, np.inf)
        nev = np.bincount(zc["ch_ev"], minlength=z["ev_n_kept"].size)     # kept channels per event after the exclusions
        nk_ev = nev[zc["ch_ev"]]
        lone, arr = nk_ev == 1, nk_ev >= VETO_NKEPT              # only kept channel in its event; array-wide event
        veto = lone | arr
        d["n_ev_arr"] += int((nev >= VETO_NKEPT).sum())
        c = clusters(dict({q: zc[q][~veto] for q in ("ch_k", "ch_a", "ch_t", "ch_ev", "ch_midx", "ch_nclip", "ch_trunc_main")}, n_pass=int(z["n_pass"])), nbr, athr)
        cl = d["cl"]
        for key in ("n_ev_kept", "clip", "cut", "cutin", "n_cl_emain"):
            cl[key] += int(c[key])
        cl["n_cl"] += int(c["n_cl"])
        for key in ("dt_all", "r_all", "dt_in", "r_in"):
            cl[key].append(c[key])
        k, a, t = zc["ch_k"], zc["ch_a"], zc["ch_t"]
        two, one = k == 2, (k == 1) & (a[:, 0] <= 1000)
        valid = np.arange(KMAX) < k[:, None]
        main = np.where(valid, a, -1.0).argmax(1)
        amain = a[np.arange(k.size), main].astype(np.float64)
        oth = valid & (np.arange(KMAX) != main[:, None]) & (a > CUT_R * amain[:, None]) \
            & (np.abs(t - t[np.arange(k.size), main][:, None]) >= CUT_SEP_NS)
        cal = (amain >= athr[zc["ch_midx"]]) & ~zc["ch_trunc_main"]  # E_main >= threshold, main pulse not truncated (False: no calibration)
        cut0 = cal & oth.any(1)
        cutm = np.nonzero(cut0 & ~veto)[0]
        j = np.where(oth, a, -1.0).argmax(1)[cutm]               # the extra with the largest A_other
        for key, v in dict(k=k, amain=amain, big=cal & ~veto, eb=np.where(veto | zc["ch_trunc_main"], np.nan, amain * gain[zc["ch_midx"]] / d["beam"][run]), nocal=~np.isfinite(gain[zc["ch_midx"]]) & ~veto, veto=lone, cut=cut0 & ~veto, cut_veto=cut0 & lone, arr=arr, cut_arr=cut0 & arr,
                           cut_dt=t[cutm, j] - t[cutm, main[cutm]], cut_r=a[cutm, j] / amain[cutm], midx=zc["ch_midx"], keep=~veto, dsig=zc["ch_dped"] / zc["ch_sig"], nclip=zc["ch_nclip"], chi2_k1=zc["ch_chi2"][one], dt=t[two, 1] - t[two, 0],
                           ratio=a[two, 1] / a[two, 0]).items():
            d["ch"].setdefault(key, []).append(v)
        for c in (0, 1):
            e = {key[4:]: z[key] for key in z.files if key.startswith(f"ex{c}_")}
            if e:
                good = (clipped(e["samples"]) < SAT_NCLIP) & ~(e["chi2"] > CHI2_MAX)
                e = {key: v[good] for key, v in e.items()}
                e["run"] = np.full(e["key"].size, run)
                for key, v in e.items():
                    d["ex"][c].setdefault(key, []).append(v)
            e = {key[4:]: z[key] for key in z.files if key.startswith(f"pm{c}_")}
            if e:
                good = (clipped(e["samples"]) < SAT_NCLIP) & ~(e["chi2"] > CHI2_MAX)
                e = {key: v[good] for key, v in e.items()}
                e["run"] = np.full(e["key"].size, run)
                for key, v in e.items():
                    d["pm"][c].setdefault(key, []).append(v)
    for d in D.values():
        d["dead"] = np.isin(names, sorted(d["dead"]))             # dead in any run of the type
        d["gain"] = np.concatenate([G[r] for r in sorted(d["runs"])])
        for key in ("dt_all", "r_all", "dt_in", "r_in"):
            d["cl"][key] = np.concatenate(d["cl"][key])
        d["ch"] = {k: np.concatenate(v) for k, v in d["ch"].items()}
        d["ex"] = {c: {k: np.concatenate(v) for k, v in e.items()} for c, e in d["ex"].items()}
        d["pm"] = {c: {k: np.concatenate(v) for k, v in e.items()} for c, e in d["pm"].items()}
    return files, D, names, pwo4


def k_rates(c, nm):
    """[KMAX - 1, nm] per-module N(k >= n) / N(kept channels after the vetoes and flags) for n = 2 ... KMAX (row n - 2; truncated channels included);
    nan for modules with fewer than MIN_N_MOD such channels."""
    n = np.bincount(c["midx"][c["keep"]], minlength=nm)
    cnt = np.stack([np.bincount(c["midx"][c["keep"] & (c["k"] >= j)], minlength=nm) for j in range(2, KMAX + 1)])
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(n >= MIN_N_MOD, cnt / n, np.nan)


def grid(d):
    """k == 2 pairs per (ratio, dt) cell: counts, % of kept channels and binomial error."""
    c = d["ch"]
    cnt = np.histogram2d(c["ratio"], c["dt"], bins=[R_EDGES, DT_EDGES])[0]
    n = c["k"].size
    p = cnt / max(n, 1)
    return cnt, 100 * p, 100 * np.sqrt(p * (1 - p) / max(n, 1))


def lbl(edges):
    return [f"{a:g}-{b:g}" if np.isfinite(b) else f">={a:g}" for a, b in zip(edges[:-1], edges[1:])]


def ex_title(e, i, rt, names):
    k = int(e["k"][i])
    t, a = e["t"][i][:k], e["a"][i][:k]
    s = (f"{names[e['midx'][i]]}  {rt}  run {e['run'][i]} ev {e['evnum'][i]}\n"
         f"k = {k}: t = {', '.join(f'{x:.1f}' for x in t)} ns; A = {', '.join(f'{x:.0f}' for x in a)} ADC\n")
    if k >= 2:
        cl = cell(e, i)
        s += (f"dt = {t[1] - t[0]:.1f} ns, ratio A2/A1 = {a[1] / a[0]:.3f}; cell dt {lbl(DT_EDGES)[cl[0]]} ns, "
              f"r {lbl(R_EDGES)[cl[1]]}\n" if cl else f"dt = {t[1] - t[0]:.1f} ns, ratio A2/A1 = {a[1] / a[0]:.3f} (k >= 3 cell); ")
    return s + f"chi2/dof = {e['chi2'][i]:.2f}, ped = {e['ped'][i]:.1f}"


def wrap(text, width):
    """Plot text with every line wrapped at width characters (the plot-style sizes make long titles and footnotes wide)."""
    return "\n".join(textwrap.fill(p, width, break_on_hyphens=False) if p else "" for p in text.split("\n"))


def page_fit(fig, title, note):
    """tight_layout that leaves room for the page title (figure.titlesize, bold) and the footnote (14 pt), by line count."""
    h = fig.get_size_inches()[1] * 72
    fig.tight_layout(rect=(0, 0.01 + (note.count("\n") + 1) * 17 / h, 1, 1 - ((title.count("\n") + 1) * 26 + 10) / h))


def draw_panel(ax, e, i, title):
    """One example: raw samples, pedestal, total fit and the fitted pulses."""
    tt = DT * np.arange(NS_EXPECT)
    ped = e["ped"][i]
    ax.plot(tt, e["samples"][i], "k.", ms=3, label="data")
    ax.axhline(ped, color="0.6", lw=0.8)
    ax.plot(tt, ped + e["comp"][i].sum(0), "r-", lw=1, label="fit")
    for c, j in zip(("tab:blue", "tab:orange", "tab:green", "tab:purple", "tab:brown", "tab:pink"), range(int(e["k"][i]))):
        ax.plot(tt, ped + e["comp"][i][j], "--", color=c, lw=0.8)
    ax.set_xlim(0, DT * NS_EXPECT)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("t (ns)")
    ax.tick_params(labelsize=12)


def draw_examples(pdf, plt, picks, names, title):
    for p0 in range(0, len(picks), 10):
        fig, axs = plt.subplots(2, 5, figsize=(30, 13))
        for ax in axs.flat:
            ax.axis("off")
        for ax, (rt, e, i) in zip(axs.flat, picks[p0:p0 + 10]):
            ax.axis("on")
            draw_panel(ax, e, i, ex_title(e, i, rt, names))
        ttl = wrap(f"{title} ({p0 + 1}-{min(p0 + 10, len(picks))} of {len(picks)})", 130)
        fig.suptitle(ttl)
        note = wrap("black: raw samples; grey: pedestal; red: total fit; dashed: fitted pulses in time "
                    "order (blue, orange, green, purple, brown, pink).  t = onset (ns after readout start), A in ADC.", 160)
        fig.text(0.01, 0.005, note, fontsize=14)
        page_fit(fig, ttl, note)
        pdf.savefig(fig)
        plt.close(fig)


def cell(e, i):
    """(dt bin, ratio bin) index into lbl(DT_EDGES), lbl(R_EDGES) of a k == 2 example, else None."""
    if e["k"][i] != 2:
        return None
    return (int(np.digitize(e["t"][i, 1] - e["t"][i, 0], DT_EDGES)) - 1,
            int(np.digitize(e["a"][i, 1] / e["a"][i, 0], R_EDGES)) - 1)


def far_order(cells):
    """Cells (dt bin, ratio bin) in farthest-point order, starting at the first."""
    cells, order = list(cells), []
    while cells:
        d = [min((abs(x[0] - y[0]) + abs(x[1] - y[1]) for y in order), default=0) for x in cells]
        order.append(cells.pop(int(np.argmax(d))))
    return order


def pick_examples(D, c, n):
    """k == 1 (c = 0): round-robin over run types.  Pile-up (c = 1): one example per populated
    (dt, ratio) cell of k == 2 (plus one k >= 3 cell) before any cell repeats; cells are taken in
    farthest-point order in (dt bin, ratio bin) starting at the shortest dt / smallest ratio, so
    short and long dt and small and large ratios appear early; within a cell run types alternate,
    each in random-key order."""
    pool = {}                                              # cell -> run type -> examples
    for rt in RUN_TYPES:
        e = D.get(rt, {}).get("ex", {}).get(c, {})
        if not e or e["key"].size == 0:
            continue
        for i in np.argsort(e["key"]):
            pool.setdefault(cell(e, i) if c else None, {}).setdefault(rt, []).append((rt, e, int(i)))
    if c == 0:
        per = list(pool.get(None, {}).values())
    else:
        order = far_order(sorted(x for x in pool if x is not None))
        order.insert(min(3, len(order)), None)             # the k >= 3 cell
        per = [_alternate(pool[cl]) for cl in order if cl in pool]
        n_cells = len(per)
        out = [lst[0] for lst in per][:n]                  # first pass: one per cell
        rest = [x for lst in per for x in lst[1:]]
        return out + rest[:max(n - n_cells, 0)]
    out = []
    while len(out) < n and any(per):
        for lst in per:
            if lst and len(out) < n:
                out.append(lst.pop(0))
    return out


def _alternate(by_rt):
    """Examples of one cell, alternating run types."""
    lsts, out = [list(v) for v in by_rt.values()], []
    while any(lsts):
        out += [lst.pop(0) for lst in lsts if lst]
    return out


def template_class(rt, names, pwo4):
    """Per module, the template the run type's bank gives it: own elastic (pool of elastic-energy pulses), own general (pool of
    general pulses; only if the bank has a pool_general key), single (one pulse), or average (lead glass, no valid row, fallback)."""
    z = np.load(BANK_ROOT / rt / f"templates_{rt}.npz", allow_pickle=False)
    row = {str(n): i for i, n in enumerate(z["names"])}
    general = z["pool_general"] > 0 if "pool_general" in z.files else np.zeros(len(row), dtype=bool)
    out = []
    for nm, p in zip(names, pwo4):
        r = row.get(str(nm))
        if not p or r is None or z["valid"][r] <= 0 or z["dead"][r] > 0 or z["fallback"][r] > 0:
            out.append("average")
        else:
            out.append("single" if z["single"][r] > 0 else "own general" if general[r] else "own elastic")
    return np.array(out)


def add_components(e, bank):
    """Fitted pulses e["comp"] [n, KMAX, samples] = amplitude x the module's unit shape at the fitted onset (as fit())."""
    tt = DT * np.arange(NS_EXPECT, dtype=F32)
    m = bank.row[e["midx"]]
    e["comp"] = np.stack([e["a"][:, j, None] * bank.shape(m[:, None], tt[None, :] - np.nan_to_num(e["t"][:, j])[:, None])
                          for j in range(KMAX)], axis=1)


def cut_fields(e):
    """Main amplitude, largest extra amplitude, dt = t_extra - t_main and r = A_extra / A_main of pile-up (cut) examples."""
    main, amain, _, j = main_other(e["k"], e["a"], e["t"])
    n = np.arange(e["k"].size)
    return dict(amain=amain, aex=e["a"][n, j], dt=e["t"][n, j] - e["t"][n, main], r=e["a"][n, j] / amain)


def cut_title(e, i, names):
    return (f"{names[e['midx'][i]]}  run {e['run'][i]}  event {e['evnum'][i]}\n"
            f"k = {int(e['k'][i])}; A_main {e['amain'][i]:.0f}, A_extra {e['aex'][i]:.0f} ADC; Δt = {e['dt'][i]:+.1f} ns; r = {e['r'][i]:.2f}\n"
            f"E_main {e['emain'][i]:.0f} MeV; chi2/dof {e['chi2'][i]:.2f}")


def clean_title(e, i, names):
    return (f"{names[e['midx'][i]]}  run {e['run'][i]}  event {e['evnum'][i]}  (one pulse, no pile-up)\n"
            f"A_main {e['a'][i, 0]:.0f} ADC; E_main {e['emain'][i]:.0f} MeV; chi2/dof {e['chi2'][i]:.2f}")


def draw_grid(pdf, plt, items, title, note):
    """One page of 2 x 4 panels; items = [(examples, index, panel title)]."""
    fig, axs = plt.subplots(2, 4, figsize=(26, 12))
    for ax in axs.flat:
        ax.axis("off")
    for ax, (e, i, text) in zip(axs.flat, items):
        ax.axis("on")
        draw_panel(ax, e, i, text)
    ttl = wrap(title, 130)
    fig.suptitle(ttl)
    note = wrap(note + "  Black: raw samples; grey: pedestal; red: total fit; dashed: fitted pulses in time order "
                "(blue, orange, green, purple, brown, pink).  t = onset (ns after readout start), A in ADC.", 160)
    fig.text(0.01, 0.005, note, fontsize=14)
    page_fit(fig, ttl, note)
    pdf.savefig(fig)
    plt.close(fig)


def ring_modules(names, pwo4, geo):
    """The N_RING PbWO4 modules closest to the beam-hole centre."""
    idx = np.nonzero(pwo4)[0]
    return idx[np.argsort([np.hypot(geo[names[i]]["x"], geo[names[i]]["y"]) for i in idx], kind="stable")][:N_RING]


def pick_spread(e, n):
    """n pile-up (cut) examples spread over the (|Δt|, r) cells (cells in far_order, one example per cell per round) and over modules
    (within a cell, a module not yet shown comes first)."""
    cl = list(zip(np.digitize(np.abs(e["dt"]), DT_EDGES) - 1, np.digitize(e["r"], R_EDGES) - 1))
    pool = {}
    for i in np.argsort(e["key"]):
        pool.setdefault(cl[i], []).append(int(i))
    order, shown, out = far_order(sorted(pool)), set(), []
    while len(out) < n and any(pool.values()):
        for c in order:
            if pool[c] and len(out) < n:
                i = next((q for q in pool[c] if e["midx"][q] not in shown), pool[c][0])
                pool[c].remove(i)
                shown.add(e["midx"][i])
                out.append(i)
    return out


def write_examples(out, D, rts, names, pwo4, geo, mc, mod_map):
    """examples.pdf and hot_modules.pdf from the per-module reservoirs (see the module docstring)."""
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    nbr = neighbours(names)
    rts = sorted((rt for rt in rts if D[rt]["pm"][1]), key=lambda rt: rt != RING_RUN_TYPE)   # the ring run type first; scans without per-module examples are skipped
    hot, cls, thr = {}, {}, {}
    for rt in rts:
        pm = D[rt]["pm"]
        bank = Bank(BANK_ROOT / rt / f"templates_{rt}.npz", names, pwo4)
        for c in (0, 1):
            if pm[c]:
                add_components(pm[c], bank)
        pm[1].update(cut_fields(pm[1]))
        cls[rt] = template_class(rt, names, pwo4)
        thr[rt] = float(np.median([CUT_EFRAC * b for b in D[rt]["beam"].values()]))
        nb, nc = mc[rt]
        with np.errstate(invalid="ignore", divide="ignore"):
            rate = np.where(pwo4 & (nb >= MIN_N_MOD), nc / nb, np.nan)
        order = np.argsort(-np.nan_to_num(rate, nan=-1.0), kind="stable")
        hot[rt] = [int(i) for i in order if np.isfinite(rate[i])][:HOT_N]
    sel = (f"E_main >= {CUT_EFRAC:g} E_beam and another pulse with A_other / A_main > {CUT_R:g}, |Δt| >= {CUT_SEP_NS:g} ns; "
           f"single-channel and array-wide vetoes applied")
    with PdfPages(out / "examples.pdf") as pdf:
        for rt in rts:
            e = D[rt]["pm"][1]
            if rt == RING_RUN_TYPE:
                ring = ring_modules(names, pwo4, geo)
                items = []
                for m in ring:
                    ii = np.nonzero(e["midx"] == m)[0]
                    items += [(e, i, cut_title(e, i, names)) for i in ii[np.argsort(e["key"][ii])][:N_RING_EX]]
                miss = [names[m] for m in ring if not (e["midx"] == m).any()]
                npg = -(-len(items) // 8)
                for p in range(npg):
                    draw_grid(pdf, plt, items[8 * p:8 * p + 8],
                              f"{rt}: pile-up (cut) in the beam-hole ring, the {N_RING} PbWO4 modules closest to the hole ({', '.join(names[ring])}); "
                              f"up to {N_RING_EX} per module (page {p + 1} of {npg})\n{sel}",
                              "No example in the sample: " + (", ".join(miss) if miss else "none") + ".")
            for p, i0 in enumerate(range(0, N_SPREAD, 8)):
                picks = pick_spread(e, N_SPREAD)[i0:i0 + 8]
                if picks:
                    draw_grid(pdf, plt, [(e, i, cut_title(e, i, names)) for i in picks],
                              f"{rt}: pile-up (cut) examples spread over modules and over the (|Δt|, r) grid (page {p + 1} of {N_SPREAD // 8})\n{sel}",
                              "Cells of |Δt| and r in turn, a module not yet shown first.")
        cols = [(f"bad waveform (any sample 0 or more than {UNDERSHOOT_ADC:g} ADC below the baseline)", "bad"), (f"saturated (clipped samples >= {SAT_NCLIP})", "sat"),
                (f"bad fit (chi-squared per degree of freedom > {CHI2_MAX:g})", "bfit")]
        fig, axs = plt.subplots(len(RUN_TYPES), 3, figsize=(21, 28), layout="constrained")
        for j, (lab, q) in enumerate(cols):
            allv = np.concatenate([np.where(D[rt]["dead"], np.nan, D[rt][q].astype(float))[pwo4] for rt in rts])
            vmax = float(np.nanmax(allv)) or 1.0                          # one colour scale per column
            for ax, rt in zip(axs[:, j], RUN_TYPES):
                if rt not in D:
                    ax.axis("off")
                    continue
                pc = mod_map(ax, np.where(D[rt]["dead"], np.nan, D[rt][q].astype(float)), vmax, D[rt]["dead"])
                ax.set_title(f"{rt}: {lab}", fontsize=11)
            fig.colorbar(pc, ax=axs[:, j], shrink=0.5, location="bottom").ax.tick_params(labelsize=12)
        fig.suptitle("excluded channels per module: bad waveform, saturated, bad fit (count per PbWO4 module; one colour scale per column; hatched: dead)")
        pdf.savefig(fig)
        plt.close(fig)
        pages = []
        pages += [(f"N(k >= {n}) / N(kept channels after the vetoes and flags) per PbWO4 module (%, colour scale up to the 99th percentile)",
                   {rt: 100 * k_rates(D[rt]["ch"], names.size)[n - 2] for rt in D}, 99) for n in range(2, KMAX + 1)]
        for title, vals, pct in pages:
            vals = {rt: np.where(D[rt]["dead"], np.nan, v) for rt, v in vals.items()}
            allv = np.concatenate([v[pwo4] for v in vals.values()])
            allv = allv[np.isfinite(allv)]
            vmax = (float(np.percentile(allv, pct) if pct else allv.max()) if allv.size else 0.0) or 1.0   # shared colour scale; 1 if all zero
            fig, axs = plt.subplots(2, 2, figsize=(16, 14), layout="constrained")
            for ax, rt in zip(axs.flat, RUN_TYPES):
                if rt not in vals:
                    ax.axis("off")
                    continue
                pc = mod_map(ax, vals[rt], vmax, D[rt]["dead"])
                ax.set_title(rt, fontsize=14)
            fig.suptitle(wrap(f"{title}; max {vmax:.3g}; " + (f"grey: fewer than {MIN_N_MOD} kept channels, " if pct else "") + "hatched: dead", 105))
            cb = fig.colorbar(pc, ax=axs, shrink=0.7)
            cb.ax.tick_params(labelsize=14)
            pdf.savefig(fig)
            plt.close(fig)
    with PdfPages(out / "hot_modules.pdf") as pdf:
        rows = []
        for rt in rts:
            nb, nc = mc[rt]
            for rank, i in enumerate(hot[rt], 1):
                g = geo[names[i]]
                dn = names[nbr[i] & D[rt]["dead"]]
                rows.append([rt, rank, names[i], f"{g['x']:.0f}", f"{g['y']:.0f}", int(max(abs(g["x"]), abs(g["y"])) / g["sx"]),
                             nb[i], nc[i], "%.2f ± %.2f" % tuple(100 * v for v in binom(nc[i], nb[i])), cls[rt][i],
                             ", ".join(dn) if dn.size else "none"])
        fig, ax = plt.subplots(figsize=(24, 16))
        ax.axis("off")
        tb = ax.table(cellText=rows, colLabels=["run type", "rank", "module", "x (mm)", "y (mm)", "ring", "N(E_main >= thr)", "N(cut)",
                                                "cut rate (%)", "template", "dead neighbours"], loc="upper center", cellLoc="center")
        tb.auto_set_font_size(False)
        tb.set_fontsize(12)
        tb.auto_set_column_width(range(11))
        tb.scale(1, 1.5)
        ax.set_title(wrap(f"The {HOT_N} PbWO4 modules with the highest pile-up (cut) rate per run type (modules with >= {MIN_N_MOD} channels at E_main >= {CUT_EFRAC:g} E_beam)", 110))
        foot = wrap(f"Pile-up (cut): {sel}.\nRing = floor(max(|x|, |y|) / module size); ring 1 = the beam-hole ring.\n"
                    "Template: own elastic / own general (module's own spline from elastic / general pulses), single (one pulse), average (the bank's average spline).\n"
                    "Dead neighbours: dead modules within 1.5 x the mean module size.", 170)
        fig.text(0.02, 0.02, foot, fontsize=14)
        pdf.savefig(fig)
        plt.close(fig)
        for rt in rts:
            nb, nc = mc[rt]
            for i in hot[rt]:
                e1, e0 = D[rt]["pm"][1], D[rt]["pm"][0]
                u = np.nonzero(e1["midx"] == i)[0]
                c = np.nonzero(e0["midx"] == i)[0] if e0 else np.zeros(0, dtype=int)
                items = [(e1, k, cut_title(e1, k, names)) for k in u[np.argsort(e1["key"][u])][:HOT_PU]]
                items += [(e0, k, clean_title(e0, k, names)) for k in c[np.argsort(e0["key"][c])][:HOT_CLEAN]]
                g = geo[names[i]]
                dn = names[nbr[i] & D[rt]["dead"]]
                r, er = binom(nc[i], nb[i])
                draw_grid(pdf, plt, items,
                          f"{rt}: {names[i]} at x = {g['x']:.0f}, y = {g['y']:.0f} mm, ring {int(max(abs(g['x']), abs(g['y'])) / g['sx'])}; "
                          f"N(cut) / N(E_main >= {thr[rt]:.0f} MeV) = {nc[i]} / {nb[i]} = {100 * r:.2f} ± {100 * er:.2f} %; template: {cls[rt][i]}; "
                          f"dead neighbours: {', '.join(dn) if dn.size else 'none'}\n{HOT_PU} pile-up (cut) examples, then {HOT_CLEAN} one-pulse examples of the module; {sel}",
                          "")


def report(scan_dir, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.collections import PatchCollection
    from matplotlib.colors import LogNorm
    from matplotlib.patches import Rectangle
    files, D, names, pwo4 = load(scan_dir)
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    geo = {m["n"]: m["geo"] for m in json.load(open(HYCAL_MAP))}
    rts = [rt for rt in RUN_TYPES if rt in D]
    nm = names.size
    summ, lines = [], [
        f"pileup_rates.py report: {len(files)} run scans ({scan_dir})",
        f"MF: nsigma {NSIGMA:g}, min_amp {MIN_AMP:g}, kmax {KMAX}, r_model {R_MODEL:g}, MIN_SEP_NS = {MIN_SEP_NS:g} ns, "
        f"clipped samples weight 0; kept channel: A_main >= {KEEP_ADC:g} ADC (MF only if raw max - baseline >= {PRECUT_ADC:g} ADC); "
        f"note: the raw k >= 2 counts and the k maps now use the {KEEP_ADC:g} ADC denominator",
        f"bad waveforms (any sample 0 or more than {UNDERSHOOT_ADC:g} ADC below the baseline) are removed before the matched filter and excluded from everything; "
        "channels whose main pulse is truncated by the end of the waveform stay in the raw k counts and k maps but are excluded from every E_main-based quantity",
        "channels that are the only kept channel in their event are vetoed (self-triggered bursts in a few modules), as are channels of events with "
        f">= {VETO_NKEPT} kept channels (array-wide pulses); the cut rates, maps and clusters exclude them, raw k >= 2 numbers do not",
        f"saturated channels (nclip >= {SAT_NCLIP}) and bad fits (chi2/dof > {CHI2_MAX:g}) are excluded like bad waveforms (cut rates, clusters, slices, medians, maps, examples; only counted per module)",
        "events: trigger gate + LMS veto only; channels: PbWO4 and lead glass (lead glass fitted with the bank's average PbWO4 spline)",
        f"cluster: seed = kept channel with the largest pulse with onset in ({PRE_LO:g}, {PRE_HI:g}) ns (as the benchmark's pre-window, time offsets off); "
        f"members = seed + modules within {CL_RADIUS:g} x the mean module size that have an in-window pulse; used pulse = largest in-window one, "
        f"other pulses = extras (in-window or not; r = A_extra / A_used)",
        f"per-module rates (pileup_modules_*.tsv, maps): raw k >= 2 rate only for modules with n_kept > {MIN_N_MOD}; cut rate only for modules with >= {MIN_N_MOD} channels at E_main >= threshold",
        ""]
    thr_rt = {rt: sorted(CUT_EFRAC * b for b in D[rt]["beam"].values()) for rt in rts}      # per-run E_main thresholds (MeV)
    mev = lambda v: f"{v[0]:.1f}" if v[0] == v[-1] else f"{v[0]:.1f}-{v[-1]:.1f}"           # noqa: E731
    thr_txt = ", ".join(f"{rt} {mev(thr_rt[rt])} MeV" for rt in rts)
    lines.insert(1, f"pile-up (cut): E_main >= {CUT_EFRAC:g} E_beam of the run (runinfo beam_energy; {thr_txt}; E_main = A_main x gain, gain = calibration factor of the run "
                    f"(database/runinfo/general.json -> database/calibration/*.json) x I_int of the module's spline template; modules without a factor are left out) "
                    f"and another fitted pulse with A_other / A_main > {CUT_R:g} and |dt| >= {CUT_SEP_NS:g} ns "
                    f"(suppresses noise-driven secondaries); cut rates = N(pile-up (cut)) / N(kept channels with E_main >= thr); cluster cut rates = clusters with a member at E_used >= thr and an extra pulse with r > {CUT_R:g}, over leading clusters with seed E_used >= thr; raw (uncut) rates are comparison columns")
    rect = {i: Rectangle((geo[names[i]]["x"] - geo[names[i]]["sx"] / 2, geo[names[i]]["y"] - geo[names[i]]["sy"] / 2),
                         geo[names[i]]["sx"], geo[names[i]]["sy"]) for i in np.nonzero(pwo4)[0]}
    ns_ = len(SLICE_EDGES) - 1
    sl_lbl = [f"{a:.1f}–{b:.1f}" for a, b in zip(SLICE_EDGES[:-1], SLICE_EDGES[1:])]

    def mod_rate(nc, n):                                # per-module cut rate (fraction) and error, nan unless >= MIN_N_MOD channels
        with np.errstate(invalid="ignore", divide="ignore"):
            r = np.where(n >= MIN_N_MOD, nc / n, np.nan)
            return r, np.sqrt(r * (1 - r) / n)

    def mod_map(ax, rate, vmax, dead, vmin=None):                  # PbWO4 module map of rate (%), grey where nan, hatched where dead; returns the collection
        has = [i for i in rect if np.isfinite(rate[i])]
        ax.add_collection(PatchCollection([rect[i] for i in rect if i not in has and not dead[i]], facecolor="0.8", edgecolor="0.6", lw=0.2))
        ax.add_collection(PatchCollection([rect[i] for i in rect if dead[i]], facecolor="white", edgecolor="0.3", lw=0.2, hatch="////"))
        pc = PatchCollection([rect[i] for i in has], cmap="viridis", edgecolor="0.3", lw=0.2)
        pc.set_array(rate[has])
        if vmin is None:
            pc.set_clim(0, vmax)
        else:
            pc.set_norm(LogNorm(vmin, vmax))                # logarithmic scale; the caller puts zero rates at vmin
        ax.add_collection(pc)
        ax.set_xlim(-MAP_LIM, MAP_LIM)
        ax.set_ylim(-MAP_LIM, MAP_LIM)
        ax.set_aspect("equal")
        return pc

    mc, sl, nkeep = {}, {}, {}                                     # per run type: module counts (E_main >= thr: n, n_cut); [module, slice] counts
    for rt in rts:
        c = D[rt]["ch"]
        mc[rt] = tuple(np.bincount(c["midx"][c[q]], minlength=nm) for q in ("big", "cut"))
        si = np.digitize(c["eb"], SLICE_EDGES) - 1      # nan -> len(edges): outside
        ok = (si >= 0) & (si < ns_)
        sl[rt] = [np.bincount(c["midx"][m] * ns_ + si[m], minlength=nm * ns_).reshape(nm, ns_) for m in (ok, ok & c["cut"])]
        nkeep[rt] = np.bincount(c["midx"][c["keep"]], minlength=nm)      # kept channels per module after the vetoes and exclusions, every energy
    vmax = lambda rates: max((float(np.nanmax(100 * v[pwo4])) for v in rates if np.isfinite(v[pwo4]).any()), default=0.0)   # noqa: E731
    map_vmax = vmax(mod_rate(*mc[rt][::-1])[0] for rt in rts)        # shared colour scales: the four maps; all slice maps

    def sl_rate(rt, i):                                 # per-module slice rate (%): N(cut in slice) / N(kept, every energy); nan below SLICE_MIN_N channels in the slice
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where((sl[rt][0][:, i] >= SLICE_MIN_N) & (nkeep[rt] > 0), 100 * sl[rt][1][:, i] / nkeep[rt], np.nan)
    sl_all = np.concatenate([sl_rate(rt, i)[pwo4 & ~D[rt]["dead"]] for rt in rts for i in range(ns_)])
    sl_nz = sl_all[np.isfinite(sl_all) & (sl_all > 0)]
    sl_vmin, sl_vmax = max(float(np.percentile(sl_nz, 1)), 1e-4), float(np.percentile(sl_nz, 99))   # shared log colour scale (%)
    print(f"[report] slice colour scale: vmin {sl_vmin:.4g} %, vmax {sl_vmax:.4g} %")
    with PdfPages(out / "pileup_rates.pdf") as pdf:
        fig_pages = []
        for rt in rts:
            d = D[rt]
            c = d["ch"]
            n, k = c["k"].size, c["k"]
            npu = int((k >= 2).sum())
            r_ch, e_ch = binom(npu, n)
            cl = d["cl"]
            ncl = cl["n_cl"]
            cr = {key: binom(cl[key], ncl) for key in ("clip",)}
            ncl5 = cl["n_cl_emain"]
            cr.update({key: binom(cl[key], ncl5) for key in ("cut", "cutin")})
            med = lambda v: float(np.median(v)) if len(v) else np.nan        # noqa: E731
            cdt, cr_ = c["cut_dt"], c["cut_r"]
            ext = dict(ch_cut_dt_med=med(cdt), ch_cut_absdt_med=med(np.abs(cdt)), ch_cut_r_med=med(cr_),
                       ch_cut_frac_before=float((cdt < 0).mean()) if cdt.size else np.nan,
                       cl_cut_dt_med=med(cl["dt_all"]), cl_cut_r_med=med(cl["r_all"]),
                       cl_cutin_dt_med=med(cl["dt_in"]), cl_cutin_r_med=med(cl["r_in"]))
            kd = [int((k == i).sum()) for i in range(1, KMAX + 1)]
            pw = pwo4[c["midx"]]
            npu_pw = int(((k >= 2) & pw).sum())
            r_pw, e_pw = binom(npu_pw, int(pw.sum()))
            clip = c["nclip"] > 0
            n_cl, n_cl_pu = int(clip.sum()), int((clip & (k >= 2)).sum())
            r_cl, e_cl = binom(n_cl, n)
            r_clpu, e_clpu = binom(n_cl_pu, n_cl)
            big = c["big"]
            dead = d["dead"]
            thr, gmed = float(np.median(thr_rt[rt])), float(np.nanmedian(d["gain"]))
            n_nocal_mod, n_nocal = int((np.isnan(d["gain"]) & ~np.tile(dead, len(d["runs"]))).sum() // len(d["runs"])), int(c["nocal"].sum())
            nveto, ncutveto = int(c["veto"].sum()), int(c["cut_veto"].sum())
            narr, ncutarr = int(c["arr"].sum()), int(c["cut_arr"].sum())
            nemain, ncut = int(big.sum()), int(c["cut"].sum())
            r_cut, e_cut = binom(ncut, nemain)
            dq = np.percentile(c["dsig"], [16, 50, 84]) if n else np.full(3, np.nan)
            chi_med = float(np.nanmedian(c["chi2_k1"])) if c["chi2_k1"].size else np.nan
            cnt, pct, err = grid(d)
            summ.append(dict(run_type=rt, runs=",".join(map(str, sorted(d["runs"]))), n_splits=d["splits"],
                             n_events=d["n_pass"], n_kept=n, n_k1=kd[0], n_k2=kd[1], n_k3=kd[2], n_k4=kd[3],
                             rate_ch=r_ch, err_ch=e_ch, n_kept_pbwo4=int(pw.sum()), rate_ch_pbwo4=r_pw,
                             err_ch_pbwo4=e_pw, n_clusters=ncl, n_ev_no_inwin=cl["n_ev_kept"] - ncl,
                             n_clipped=n_cl, rate_clipped=r_cl, err_clipped=e_cl,
                             n_cl_clipped=cl["clip"], rate_cl_clipped=cr["clip"][0], err_cl_clipped=cr["clip"][1],
                             n_clipped_pileup=n_cl_pu, rate_pileup_clipped=r_clpu, err_pileup_clipped=e_clpu,
                             chi2dof_med_k1_a1000=chi_med, n_single_veto=nveto, n_cut_single_veto=ncutveto, n_ev_array_veto=d["n_ev_arr"], n_array_veto=narr, n_cut_array_veto=ncutarr, thr_mev=thr, median_gain=gmed, n_nocal_modules=n_nocal_mod, n_nocal_channels=n_nocal, n_emain=nemain, n_cut=ncut, rate_ch_cut=r_cut, err_ch_cut=e_cut,
                             n_dead_modules=int(dead.sum()), rate_cl_cut=cr["cut"][0], err_cl_cut=cr["cut"][1], rate_cl_cutin=cr["cutin"][0], err_cl_cutin=cr["cutin"][1],
                             n_cl_seed_emain=ncl5, **ext))
            lines += [f"[{rt}] runs {sorted(d['runs'])}, {d['splits']} splits, {d['n_seen']} events seen, "
                      f"{d['n_lms']} LMS-vetoed, {d['n_pass']} after trigger + LMS",
                      f"  channel level, PRIMARY (cut): {ncut}/{nemain} = {100 * r_cut:.3f} +- {100 * e_cut:.3f} % (N(kept, E_main >= {thr:.1f} MeV) = {nemain})",
                      f"  dead modules (runinfo hycal_dead_modules, excluded from the scan; dead in any run of the type): {int(dead.sum())}: {', '.join(names[dead])}",
                      f"  E_main >= {CUT_EFRAC:g} E_beam = {mev(thr_rt[rt])} MeV (beam energies {', '.join(f'{b:g}' for b in sorted(set(d['beam'].values())))}); median gain {gmed:.4f} MeV/ADC (equivalent ADC threshold {thr / gmed:.0f} at the median); "
                      f"{n_nocal_mod} modules without a calibration factor, {n_nocal} kept non-veto channels in them left out of the cut rates",
                      f"  channel level, cut extras (largest A_other per channel; dt = t_other - t_main): median dt {ext['ch_cut_dt_med']:.1f} ns, "
                      f"median |dt| {ext['ch_cut_absdt_med']:.1f} ns, median r {ext['ch_cut_r_med']:.3f}, before the main pulse {100 * ext['ch_cut_frac_before']:.1f} %",
                      f"  cluster level, PRIMARY (cut): {ncl5} leading clusters with seed E_used >= {thr:.1f} MeV (of {ncl}); "
                      f"clusters with a member at E_used >= {thr:.1f} MeV and an extra with r > {CUT_R:g}, |dt| >= {CUT_SEP_NS:g} ns: "
                      f"all extras {100 * cr['cut'][0]:.3f} +- {100 * cr['cut'][1]:.3f} %, in-window extras {100 * cr['cutin'][0]:.3f} +- {100 * cr['cutin'][1]:.3f} %",
                      f"  cluster level, cut extras (largest-r triggering extra per cluster; dt = t_extra - t_used): all extras median dt {ext['cl_cut_dt_med']:.1f} ns, "
                      f"median r {ext['cl_cut_r_med']:.3f}; in-window extras median dt {ext['cl_cutin_dt_med']:.1f} ns, median r {ext['cl_cutin_r_med']:.3f}",
                      f"  single-channel veto: {nveto} of {n} kept channels vetoed ({100 * nveto / max(n, 1):.2f} %), {ncutveto} of them would be pile-up (cut)",
                      f"  array-wide veto: {d['n_ev_arr']} events, {narr} kept channels vetoed, {ncutarr} of them would be pile-up (cut)",
                      f"  channel level, raw (comparison): {npu}/{n} = {100 * r_ch:.3f} +- {100 * e_ch:.3f} % with k >= 2 "
                      f"(PbWO4 only: {npu_pw}/{int(pw.sum())} = {100 * r_pw:.3f} +- {100 * e_pw:.3f} %)",
                      f"  k distribution (k = 1 ... {KMAX}): {kd} -> k = 1 fraction {kd[0] / max(n, 1):.4f}; k = {KMAX} channels: {kd[-1]}",
                      f"  bad waveforms: {int(d['bad'].sum())} of {int(d['digit'].sum())} digitised channels ({100 * d['bad'].sum() / max(d['digit'].sum(), 1):.3f} %); "
                      f"modules with the most: {', '.join(f'{names[i]} {int(d['bad'][i])}' for i in np.argsort(-d['bad'], kind='stable')[:3])}",
                      *[f"  {lab}: {int(d[q].sum())} of {int(d['n_pre'])} kept channels before the exclusion ({100 * d[q].sum() / max(d['n_pre'], 1):.3f} %); "
                        f"modules with the most: {', '.join(f'{names[i]} {int(d[q][i])}' for i in np.argsort(-d[q], kind='stable')[:3])}"
                        for lab, q in ((f"saturated (nclip >= {SAT_NCLIP})", "sat"), (f"bad fit (chi2/dof > {CHI2_MAX:g})", "bfit"))],
                      f"  channels with a truncated main pulse: {int(d['trunc'].sum())} of {n} kept ({100 * d['trunc'].sum() / max(n, 1):.3f} %), excluded from the E_main-based quantities",
                      f"  baseline offset dped / sigma of the kept channels: median {dq[1]:.3f}, 16th-84th percentile {dq[0]:.3f} ... {dq[2]:.3f}; "
                      f"fraction clipped at ±{BASE_CLIP:g} sigma: {100 * float((np.abs(c['dsig']) >= BASE_CLIP - 1e-3).mean() if n else np.nan):.3f} %",
                      f"  leading clusters: {ncl} of {cl['n_ev_kept']} events with a kept channel "
                      f"({d['n_pass']} after trigger + LMS; {cl['n_ev_kept'] - ncl} without an in-window pulse, "
                      f"{100 * (cl['n_ev_kept'] - ncl) / max(cl['n_ev_kept'], 1):.2f} %)",
                      f"  clipped (any sample >= {OVERFLOW}): channel level {n_cl}/{n} = {100 * r_cl:.4f} +- {100 * e_cl:.4f} %, "
                      f"cluster level {cl['clip']}/{ncl} = {100 * cr['clip'][0]:.4f} +- {100 * cr['clip'][1]:.4f} %; "
                      f"pile-up (k >= 2) among clipped channels: {n_cl_pu}/{n_cl} = {100 * r_clpu:.3f} +- {100 * e_clpu:.3f} %",
                      f"  median chi2/dof of k = 1 channels with A <= 1000 ADC: {chi_med:.3f}",
                      f"  scan wall time per run: median {np.median(d['wall']):.0f} s, max {max(d['wall']):.0f} s", ""]
            with open(out / f"pileup_grid_{rt}.tsv", "w") as fh:
                fh.write("ratio_lo\tratio_hi\tdt_lo_ns\tdt_hi_ns\tn_pairs\tpct_of_kept\terr_pct\n")
                for (i, j), v in np.ndenumerate(cnt):
                    fh.write(f"{R_EDGES[i]:g}\t{R_EDGES[i + 1]:g}\t{DT_EDGES[j]:g}\t{DT_EDGES[j + 1]:g}\t"
                             f"{int(v)}\t{pct[i, j]:.5f}\t{err[i, j]:.5f}\n")
            fig, ax = plt.subplots(figsize=(16, 8))
            ax.imshow(pct, cmap="Greens", vmin=0, vmax=max(pct.max(), 1e-9), aspect="auto")
            for (i, j), v in np.ndenumerate(pct):
                ax.text(j, i, f"{v:.3f} ± {err[i, j]:.3f}", ha="center", va="center", fontsize=11,
                        color="white" if v > 0.7 * pct.max() else "black")
            ax.set_xticks(range(len(DT_EDGES) - 1), lbl(DT_EDGES))
            ax.set_yticks(range(len(R_EDGES) - 1), lbl(R_EDGES))
            ax.set_xlabel("time gap dt = t_later - t_earlier (ns)")
            ax.set_ylabel("amplitude ratio A_later / A_earlier")
            ttl = wrap(f"{rt}: k = 2 channels per (ratio, dt) cell, % of all {n} kept channels "
                       f"(sum {pct.sum():.3f} % = k = 2 rate; all k >= 2: {100 * r_ch:.3f} %)", 100)
            ax.set_title(ttl)
            note = wrap(f"cell = N(k = 2 pairs in cell) / N(kept channels, A_main >= {KEEP_ADC:g} ADC) ± binomial "
                        f"error (%); MIN_SEP_NS = {MIN_SEP_NS:g} ns, so no pair has dt < {DT_EDGES[0]} ns.", 130)
            fig.text(0.01, 0.005, note, fontsize=14)
            page_fit(fig, ttl, note)
            fig_pages.append(fig)

            nk = np.bincount(c["midx"], minlength=nm)
            np_ = np.bincount(c["midx"][k >= 2], minlength=nm)
            with np.errstate(invalid="ignore", divide="ignore"):
                rate = np.where(nk > 0, np_ / nk, np.nan)
                rerr = np.sqrt(rate * (1 - rate) / nk)
            nclip = np.bincount(c["midx"][clip], minlength=nm)
            ok = nk > MIN_N_MOD
            rate, rerr = np.where(ok, rate, np.nan), np.where(ok, rerr, np.nan)
            nb, nc = mc[rt]
            rc, rce = mod_rate(nc, nb)
            okc = np.isfinite(rc)
            med, p90 = (np.percentile(100 * rc[okc & pwo4], q) if (okc & pwo4).any() else np.nan for q in (50, 90))
            lines.insert(len(lines) - 1, f"  PbWO4 modules with n_kept > {MIN_N_MOD}: {int((ok & pwo4).sum())} of {int((pwo4 & ~dead).sum())} (dead excluded); "
                         f"with >= {MIN_N_MOD} channels at E_main >= thr (cut rate): {int((okc & pwo4).sum())}, "
                         f"per-module cut rate median {med:.2f} %, p90 {p90:.2f} %")
            sc = sl[rt][1].sum(0)
            lines.insert(len(lines) - 1, f"  pooled cut rate (%) in E_main/E_beam slices {sl_lbl[0]} ... {sl_lbl[-1]} = N(cut channels with E_main in the slice) / N(all {int(nkeep[rt].sum())} kept channels of the run type, every energy; "
                         f"per-module values: / kept channels of the module; drawn for >= {SLICE_MIN_N} channels in the slice): "
                         + ", ".join(f"{100 * binom(a, nkeep[rt].sum())[0]:.3f}" for a in sc))
            # rate_k2_raw, err_k2_raw: N(k >= 2) / N(all kept channels, before the vetoes); rate_k2 ... rate_k6: N(k >= n) / N(kept channels after the vetoes and exclusions), as in the k maps
            with open(out / f"pileup_modules_{rt}.tsv", "w") as fh:
                fh.write("module\tis_pbwo4\tdead\tn_kept\tn_pileup\trate_k2_raw\terr_k2_raw\tn_clipped\trate_clipped\tn_emain\tn_cut\trate_cut\terr_cut\trate_k2\trate_k3\trate_k4\trate_k5\trate_k6\tn_saturated\tn_bad_fit\n")
                rk = k_rates(c, nm)
                for i in range(nm):
                    v = (f"{rate[i]:.5f}", f"{rerr[i]:.5f}", f"{nclip[i] / nk[i]:.5f}") if ok[i] else ("", "", "")
                    u = (f"{rc[i]:.5f}", f"{rce[i]:.5f}") if okc[i] else ("", "")
                    fh.write(f"{names[i]}\t{int(pwo4[i])}\t{int(dead[i])}\t{nk[i]}\t{np_[i]}\t{v[0]}\t{v[1]}\t{nclip[i]}\t{v[2]}\t{nb[i]}\t{nc[i]}\t{u[0]}\t{u[1]}\t" + "\t".join(f"{x:.5f}" for x in rk[:, i]) + f"\t{d['sat'][i]}\t{d['bfit'][i]}\n")
            fig, ax = plt.subplots(figsize=(11, 10), layout="constrained")
            pc = mod_map(ax, 100 * rc, map_vmax, dead)
            cb = fig.colorbar(pc, ax=ax)
            cb.set_label(wrap(f"pile-up (cut) rate (% of channels with E_main >= {CUT_EFRAC:g} E_beam ({thr:.0f} MeV))\nshared scale, max {map_vmax:.1f} %", 36),
                         fontsize=16, fontweight="bold")
            cb.ax.tick_params(labelsize=14)
            ax.set_xlabel("x (mm)")
            ax.set_ylabel("y (mm)")
            ax.set_title(wrap(f"{rt}: pile-up rate per PbWO4 module, cuts E_main >= {CUT_EFRAC:g} E_beam ({thr:.0f} MeV), r > {CUT_R:g}, |dt| >= {CUT_SEP_NS:g} ns\n"
                              f"(modules with >= {MIN_N_MOD} channels at E_main >= thr; grey: fewer, hatched: dead;\nchannels alone in their event vetoed)", 62))
            fig.savefig(out / f"pileup_rate_map_{rt}.png", dpi=150)
            fig_pages.append(fig)

        with open(out / "pileup_slices.tsv", "w") as fh:
            fh.write("run_type\tslice_lo\tslice_hi\tn\tn_cut\trate\terr\n")
            for rt in rts:
                for i, c_ in enumerate(sl[rt][1].sum(0)):          # n = all kept channels of the run type (every energy)
                    fh.write(f"{rt}\t{SLICE_EDGES[i]:.1f}\t{SLICE_EDGES[i + 1]:.1f}\t{nkeep[rt].sum()}\t{c_}\t%.6f\t%.6f\n" % binom(c_, nkeep[rt].sum()))
        for rt in rts:
            fig, axs = plt.subplots(2, 3, figsize=(20, 13), layout="constrained")
            for i, ax in enumerate(axs.flat):
                c_ = sl[rt][1][:, i]
                r_ = sl_rate(rt, i)
                pc = mod_map(ax, np.where(r_ == 0, sl_vmin, r_), sl_vmax, dead, sl_vmin)     # zero rates at the floor
                ax.set_title(f"E_main/E_beam {sl_lbl[i]}: pooled {100 * binom(c_.sum(), nkeep[rt].sum())[0]:.3f} % ({int(c_.sum())} of {int(nkeep[rt].sum())} kept)", fontsize=12)
            fig.suptitle(wrap(f"{rt}: pile-up (cut) rate per PbWO4 module in E_main/E_beam slices (r > {CUT_R:g}, |dt| >= {CUT_SEP_NS:g} ns; grey: < {SLICE_MIN_N} channels in the slice, hatched: dead); "
                              f"rate = N(cut channels with E_main in the slice) / N(all kept channels of the module, every energy), pooled: of the run type; "
                              f"one logarithmic colour scale for all slices and run types, {sl_vmin:.2g}-{sl_vmax:.2g} %", 105))
            cb = fig.colorbar(pc, ax=axs, shrink=0.7, extend="both")
            cb.set_label("pile-up (cut) rate in the slice / all kept channels of the module (%); 0 drawn at the floor", fontsize=14, fontweight="bold")
            cb.ax.tick_params(labelsize=14)
            fig.savefig(out / f"pileup_rate_slices_{rt}.png", dpi=110)
            fig_pages.append(fig)

        cols = list(summ[0]) if summ else []
        with open(out / "pileup_rates.tsv", "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for s in summ:
                fh.write("\t".join(v if isinstance(v, str) else f"{v:.6g}" for v in s.values()) + "\n")
        (out / "summary.txt").write_text("\n".join(lines) + "\n")

        def pct(s, x, f=3):
            return f"{100 * s['rate_' + x]:.{f}f} ± {100 * s['err_' + x]:.{f}f}"
        pages = [
            ("Pile-up (cut) rates: channel level (matched filter)",
             ["run type", "splits", "events", "N(E_main >= thr)", "thr (MeV)", "channel cut rate (%)", "median dt (ns)", "median r", "kept ch.", "raw k >= 2 rate (%)", "k = 1/2/3/4"],
             [[s["run_type"], s["n_splits"], s["n_events"], s["n_emain"], f"{s['thr_mev']:.1f}", pct(s, "ch_cut"), f"{s['ch_cut_dt_med']:.1f}", f"{s['ch_cut_r_med']:.2f}", s["n_kept"], pct(s, "ch"),
               f"{s['n_k1']}/{s['n_k2']}/{s['n_k3']}/{s['n_k4']}"] for s in summ],
             "\n".join(lines[1:4]) + "\nchannel cut rate = N(pile-up (cut)) / N(kept with E_main >= 0.1 E_beam); raw = N(k >= 2) / N(kept); ± binomial errors.\n"
             "median dt, r: of the pile-up (cut) channels' extra pulse with the largest A_other; dt = t_other - t_main (negative = before the main pulse), r = A_other / A_main; more in summary.txt."),
            ("Pile-up (cut) rates: cluster level and clipped",
             ["run type", "clusters, seed E >= thr", "cluster cut (%)", "median dt (ns)", "median r", "cluster cut in-window (%)", "median dt in-window (ns)", "median r in-window",
              "clipped ch. (%)", "clipped clusters (%)", "pile-up in clipped (%)", "chi2/dof med."],
             [[s["run_type"], s["n_cl_seed_emain"], pct(s, "cl_cut"), f"{s['cl_cut_dt_med']:.1f}", f"{s['cl_cut_r_med']:.2f}",
               pct(s, "cl_cutin"), f"{s['cl_cutin_dt_med']:.1f}", f"{s['cl_cutin_r_med']:.2f}",
               pct(s, "clipped", 4), pct(s, "cl_clipped", 4), f"{100 * s['rate_pileup_clipped']:.1f} ± {100 * s['err_pileup_clipped']:.1f}",
               f"{s['chi2dof_med_k1_a1000']:.2f}"] for s in summ],
             f"cluster cut = N(clusters with a member at E_used >= 0.1 E_beam and an extra pulse with r > {CUT_R:g}, |dt| >= {CUT_SEP_NS:g} ns) / "
             f"N(leading clusters with seed E_used >= 0.1 E_beam); in-window: the extra pulse is inside the pre-window.\n"
             "median dt, r: per cut cluster, the triggering extra with the largest r = A_extra / A_used; dt = t_extra - t_used. \n"
             "clipped = any sample >= 4095 (channel / cluster level; pile-up = k >= 2 among clipped channels).  chi2/dof med.: k = 1 channels with A <= 1000 ADC."),
        ]
        for title, hdr, cells, foot in pages:
            fig, ax = plt.subplots(figsize=(24, 10))
            ax.axis("off")
            tb = ax.table(cellText=cells, colLabels=hdr, loc="upper center", cellLoc="center")
            tb.auto_set_font_size(False)
            tb.set_fontsize(12)
            tb.auto_set_column_width(range(len(hdr)))
            tb.scale(1, 2)
            ax.set_title(title)
            fig.text(0.02, 0.03, wrap(foot, 170), fontsize=14)
            pdf.savefig(fig)
            plt.close(fig)
        fig, axs = plt.subplots(2, 2, figsize=(16, 10), sharex=True, layout="constrained")
        for ax, rt in zip(axs.flat, rts):
            v = D[rt]["ch"]["cut_dt"]
            ax.hist(v, bins=np.arange(-100, 101, 4), color="C0")
            ax.axvline(np.median(v), color="C3", label=f"median {np.median(v):.1f} ns")
            ax.set_yscale("log")
            ax.set_title(f"{rt}: {v.size} cut channels")
            ax.legend()
        for ax in axs[1]:
            ax.set_xlabel("dt = t_other - t_main (ns)")
        for ax in axs[:, 0]:
            ax.set_ylabel("channels / 4 ns")
        fig.suptitle(wrap(f"pile-up (cut) channels: dt of the extra pulse with the largest A_other (E_main >= {CUT_EFRAC:g} E_beam, r > {CUT_R:g}, |dt| >= {CUT_SEP_NS:g} ns)", 100))
        pdf.savefig(fig)
        plt.close(fig)
        fig, axs = plt.subplots(2, 2, figsize=(16, 10), sharex=True, layout="constrained")
        for ax, rt in zip(axs.flat, rts):
            for h, lab in zip(D[rt]["chi_h"], ("one pulse (k = 1)", "two or more pulses (k >= 2)")):
                ax.stairs(h, CHI_EDGES, label=lab)
            ax.axvline(CHI2_MAX, color="C3", ls="--", label=f"exclusion threshold {CHI2_MAX:g}")
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_title(rt)
            ax.legend()
        for ax in axs[1]:
            ax.set_xlabel("chi-squared per degree of freedom")
        for ax in axs[:, 0]:
            ax.set_ylabel("kept channels per bin")
        fig.suptitle(wrap("chi-squared per degree of freedom of all kept channels before the saturation and bad-fit exclusions (10 bins per decade; values outside 0.1-1e5 in the end bins)", 100))
        pdf.savefig(fig)
        plt.close(fig)
        for fig in fig_pages:
            pdf.savefig(fig)
            plt.close(fig)
        draw_examples(pdf, plt, pick_examples(D, 1, N_EX_PU), names,
                      f"pile-up examples (k >= 2, A_main >= {A_EX:g} ADC)")
        draw_examples(pdf, plt, pick_examples(D, 0, N_EX_1), names,
                      f"no pile-up examples (k = 1, A_main >= {A_EX:g} ADC)")
    write_examples(out, D, rts, names, pwo4, geo, mc, mod_map)
    print("\n".join(lines))
    print(f"[report] wrote {out}/: pileup_rates.pdf examples.pdf hot_modules.pdf pileup_rates.tsv summary.txt pileup_grid_*.tsv "
          "pileup_modules_*.tsv pileup_rate_map_*.png pileup_rate_slices_*.png pileup_slices.tsv")


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
