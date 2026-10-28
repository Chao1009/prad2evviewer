#!/usr/bin/env python3
"""spline_template_extract.py -- sampled (spline) per-module PbWO4 pulse templates.

Plan: plans/plan_20261024.md section 6.  Four subcommands:

  select  one EVIO split -> compressed npz of clean single-pulse candidates
          (this is the only expensive pass; one process per split)
  build   pool candidates of a run type, phase-align, fit the spline
  ampcorr compare the run types' f(A) (amplitude-dependent integral factor)
  check   held-out gates G1-4..G1-9 (+ the 24386 PASS 1 macro cross-check,
          G1-1..G1-3) -> gates.json / gates.tsv / detail TSVs / PNGs

`select` uses the PASS 1 selection of Waveform_Wassim/best_ref_shape_fit.C
(`pass1_evaluate`, ported line by line) plus the vetoes PASS 1 lacks.  The
low-threshold matched filter is used ONLY to count pulses, to read the onset
and to seed (t0, A); it is never a chi2 selection.

    python3 spline_template_extract.py select \\
        /path/prad_025320.evio.00050 \\
        --out output/templates_spline/2p24GeV_H2/pulses/025320_00050.npz
    # opt-in (defaults reproduce the above exactly): --per-module-best 50 keeps the 50
    #   highest-E candidates per module (E = mf_amp * I_int * cal, MeV; --calib-file;
    #   npz mf_energy) instead of the first-come cap; with --emax-frac 1.2 (E < 1.2 E_beam)
    #   --mf-strict (mf_npk == 1) --max-modules-per-event 20 --exclude-trigger-bits 0x1000000
    #   --mf-bank output/templates_spline/<slug>/templates_<slug>.npz

    python3 spline_template_extract.py build \\
        --pulses "output/templates_spline/2p24GeV_H2/pulses/*.npz" \\
        --run-type 2p24GeV_H2 --h50 output/pulse_templates_h50.dat \\
        --out output/templates_spline/2p24GeV_H2/templates_2p24GeV_H2
    # -> <stem>.npz (spline_templates_v1, loads with mf_pulse.build_bank),
    #    <stem>.root, <stem>.pdf and summary.tsv next to them.  20 % of the
    #    candidates are held out by `holdout_mask` (for `check`).
    #    --bands 100,250,500,1000,2000,4000: f(A) from these bands (n knots; npz
    #    band_edges / shape_band / r_band / n_band) instead of lo / mid / top.
    #    --pool-order energy (best-K candidate files): pool = top n-max by mf_energy,
    #    hold-out = the next ~n-max*frac/(1-frac) (npz holdout_keys, used by `check`).
    # opt-in: --lambda-select heldout (per-module lambda from the held-out chi2, grid
    #   0..1e4; summary.tsv / meta get lam, lam_flag, n_ho, ho_chi2_*); --converge-df X
    #   (convergence threshold on max |f_k - f_(k-1)|, default 1e-3).
    # opt-in, with --fit-points direct / --lambda-select heldout: --shape-sigma S (adds
    #   (S*ref)^2 to the point variance), --lambda-chi2-sigma full (held-out chi2 with the
    #   timing + shape terms), --lambda-min X, --lambda-fallback Y (grid floor / fallback).
    # opt-in: --osc-escalate (module template with spurious wiggles: lambda x10 and rebuild
    #   until none, run-type average at 1e4; adds osc_* / lam_final to summary.tsv, lam_module).

    # opt-in Wassim single-pulse reference (best_ref_shape_fit.C, PDF section 4):
    #   select --wassim-ref  keeps per module ONLY the best valid PASS 1 candidate of the
    #     split (betterCandidate order; clip veto + PASS 1 only, no mf / other cuts);
    #   build --single-pulse --pulses "<those npz>"  one bank row per module = the best
    #     candidate over all files, R(k) = (W - bg_mean) / max(W - bg_mean) as a not-a-knot cubic
    #     interpolating spline (exp. tail only to reach 400 ns; no f(A); fallback 2 = average)

    python3 spline_template_extract.py ampcorr \\
        --templates output/templates_spline/*/templates_*.npz \\
        --out output/templates_spline/ampcorr_summary.json
    # f(A): one global correction in every .npz if the run types agree (within 2
    #   combined SE), else each keeps its own (keys ampcorr_logA/ampcorr_f/ampcorr_meta)
    # --per-module-f OUTDIR (banks built with --bands): writes COPIES
    #   OUTDIR/<slug>/templates_<slug>.npz whose ampcorr_f is a per-module table
    #   [module, ln A grid] (a module keeps its own band ratio r_b / r_mod where its
    #   stat SE < 0.3 %, else the run-type knot); the inputs are not modified.

    python3 spline_template_extract.py check \\
        --pulses "output/templates_spline/2p24GeV_H2/pulses/*.npz" \\
        --templates output/templates_spline/2p24GeV_H2/templates_2p24GeV_H2.npz \\
        --h50 output/pulse_templates_h50.dat \\
        --out output/templates_spline/2p24GeV_H2/check/
    # --amp-bins E0,...,En: the G1-7 bins (default 100,300,1000,2500).
    # cross-check (24386 split 00010) adds: --pass1-macro-out p1.root --raw-root
    #   prad_024386.00010_raw.root --macro-tsv p1_cands.tsv --evio-cands <nocap.npz>
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import heapq
import json
import math
import os
import re
import subprocess
import sys
import time
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

# =========================================================== constants / cuts ==
SCHEMA = "spline_cand_v2"   # v2: no PASS 1 shape cut, p1_shape_ok stored
NS_EXPECT = 100                 # samples per PbWO4 readout
DT_NS = 4.0                     # sample period (macro SAMPLE_DT_NS)
BLOCK = 256                     # events per batched pass (plot_pileup_examples)

# --- cut values (module level so build/check and the PDF can report them) ----
CLIP_ADC = 4095                 # cut 1: no raw sample >= this anywhere
MF_NSIGMA = 3.0                 # low-threshold mf (count / onset / seed only)
MF_MIN_AMP = 3.0
MF_KMAX = 4
MF_R_MODEL = 0.02
# Cut 2 as literally specified ("exactly one mf pulse") rejects ~99 % of the
# otherwise clean pulses: at nsigma 3 / min_amp 3 the two-tau bank (1-2 % short
# in the tail) fakes 2-3 extra pulses of 2.5-5 % of the main amplitude, 55-150 ns
# after it (25320, 3000 events: 2410 secondaries, ratio median 3.0 %, 95th pct
# 5.0 %).  A mf pulse therefore counts as a SECONDARY only if its amplitude is
# >= MF_SEC_RATIO of the largest mf pulse (the plan's own pile-up definition,
# D6: 10 %).  mf_npk (all pulses) and mf_sec_ratio are stored so `build` can
# re-apply the literal rule (mf_npk == 1) or any tighter ratio.
MF_SEC_RATIO = 0.10
MF_ONSET_MIN_NS = 60.0          # cut 3: mf onset (t0) >= 60 ns ...
FIRST_THR_MIN = 20              # ... and first_thr >= sample 20
PEAK_MAX_SAMPLE = 59            # cut 4: contained (>= 40 samples of tail)
A_STORE_MIN = 50.0              # cut 5: stored if A = raw peak - ped >= 50
BAND_B0 = (300.0, 2500.0)       # amplitude bands (tags only; build selects)
BAND_LO = (100.0, 400.0)
BAND_HI = (400.0, 2500.0)
BAND_MID = (400.0, 1000.0)      # f(A) bands (build only): lo / mid / top, non-overlapping
BAND_TOP = (1000.0, 2500.0)
P1_THR = 2.0                    # cut 6/7: PASS 1 (macro defaults)
P1_MIN_RISING = 4
P1_MIN_FALLING = 10
P1_BG_SAMPLES = 20
P1_BG_MAX_ABS = 2.0
P1_BG_MAX_RMS = 1.0
P1_BG_MAX_ROUGH = 1.0
P1_SEC_FRAC = 0.15
P1_SEC_EXCL = 6
P1_T_LO = 150.0                 # macro time window; NOT a select cut (see
P1_T_HI = 170.0                 # peak_hist), kept for pass1_evaluate.valid_p1
PED_RMS_MAX = 1.5               # cut 8: WaveAnalyzer rms
PED_SLOPE_MAX = 0.05            # cut 8: |slope| of 20 pre-pulse samples
TAIL_OFF0, TAIL_OFF1 = 40, 50   # cut 9: samples [peak+40, min(peak+50, 99)]
# The plan's tail limit (|mean| <= max(2 ADC, 2*rms) AND <= 0.5 % of A) rejects
# EVERY pulse: the real PbWO4 pulse is still at 1.4-2.7 % of A (median 2.0 %,
# independent of A) 160-200 ns after the peak (25320, 1500 events).  The tail
# cut is therefore |mean| <= max(TAIL_ABS_ADC, TAIL_RMS_MULT*rms) +
# TAIL_FRAC_A*A with TAIL_FRAC_A = 4 % (about twice the median real tail); it
# only removes late pile-up / baseline excursions.  tail_dev is stored so
# `build` can re-cut.
TAIL_ABS_ADC = 2.0
TAIL_RMS_MULT = 2.0
TAIL_FRAC_A = 0.04
PHYSICS_BITS_EXACT = 0x100      # C.PHYSICS_TRIGGER_BITS (user decision)
PHYSICS_MASK_ANY = 0xf00        # C.PHYSICS_TRIGGER_MASK (--trigger physics-any)

# Cut-flow columns, in the order they are APPLIED.  The first N_SEQ are
# sequential (each counts the channels that passed it and every earlier one);
# the rest are informational and not part of the chain.
CUTFLOW_SEQ = [
    "nseen",            # PbWO4, nsamples == 100, event passed the trigger gate
    "c1_clip_ok",       # cut 1: no raw sample >= 4095
    "c5_amp_ge50",      # cut 5: A = raw peak - WaveAnalyzer ped >= 50
    "c4_contained",     # cut 4: peak sample <= 59
    "c3a_first_thr",    # cut 3 (part): first_thr >= 20 (pre-pulse window exists)
    "c8_pedestal",      # cut 8: WA rms <= 1.5 and |slope| <= 0.05
    "c9_tail",          # cut 9: tail back at baseline
    "c2_mf_no_secondary",  # cut 2: no mf pulse >= MF_SEC_RATIO of the largest
    "c3b_mf_onset",     # cut 3 (part): every mf pulse onset >= 60 ns
                        #   (== passes all cuts)
    "stored",           # after the per-module caps
]
CUTFLOW_INFO = [
    "info_p1_shape_ok",        # of the stored candidates: the old cut 6 (PASS 1 shape)
                               #   would have passed (rising>=4, falling>=10, no secondary)
    "info_tier2",              # of the all-cut survivors: tier 2
    "info_b0",                 # ... in band b0 [300, 2500]
    "info_tier2_b0",           # ... both
    "info_fail_falling_plateau_only",  # see _plateau_only (not in the chain)
    "info_wa_npeaks_ne1",      # stored with WaveAnalyzer npeaks != 1
    "info_mf_npk_eq1",         # of the all-cut survivors: literal cut 2 (mf_npk == 1)
]
CUTFLOW_NAMES = CUTFLOW_SEQ + CUTFLOW_INFO
N_SEQ = len(CUTFLOW_SEQ)
CUT_ORDER_TEXT = (
    "nseen > c1 clip veto > c5 A>=50 (WaveAnalyzer ped) > c4 contained "
    "(peak<=59) > c3a first_thr>=20 > c8 WA rms<=1.5 & |slope|<=0.05 > "
    "c9 tail (|mean| <= max(2,2*rms)+4% A) > c2 no mf secondary >= 10% of the "
    "largest mf pulse (literal 'exactly one pulse' is stored as mf_npk and "
    "info_mf_npk_eq1) > c3b every mf pulse onset >= 60 ns > per-module caps "
    "(separately for band-b0 candidates and for all others, each first come, "
    "so a module stores <= 2*cap).  The PASS 1 shape cut (rising>=4, falling>=10 "
    "strict, no secondary >15% outside +-6) is NOT applied; the stored "
    "p1_shape_ok flag records it so `build` can re-apply it.  The (expensive) mf runs only on channels "
    "that passed everything before it.  Tier>=1 is implied by c3a (tier 0 = no pre-pulse "
    "window); tier 2 vs 1 is stored as a flag, never cut.")


# ======================================================== pure functions ====
# Port of Waveform_Wassim/best_ref_shape_fit.C (WW/best_ref_shape_fit.C),
# namespace HyCalBestRef.  Line numbers below refer to that file.

def pass1_evaluate(wf, thr=P1_THR, min_rising=P1_MIN_RISING,
                   min_falling=P1_MIN_FALLING, bg_samples=P1_BG_SAMPLES,
                   bg_max_abs=P1_BG_MAX_ABS, bg_max_rms=P1_BG_MAX_RMS,
                   bg_max_rough=P1_BG_MAX_ROUGH, sec_frac=P1_SEC_FRAC,
                   sec_excl=P1_SEC_EXCL, t_lo=P1_T_LO, t_hi=P1_T_HI):
    """PASS 1 `evaluateCandidate` on ONE pedestal-subtracted waveform.

    Returned keys: amp, peak, first_thr, rising, falling, falling_nonstrict,
    has_secondary, sec_max, bg_ok, bg_mean, bg_rms, bg_maxabs, roughness, tier
    (2 strict, 1 loose, 0 no pre-pulse window), strict, shape_ok (rising,
    falling, secondary, bg window; NO time window), time_ok (peak*4 ns in
    [t_lo, t_hi]), valid_p1 (the macro's `Candidate.valid`: amp > thr, time
    window, first_thr found, bg window, rising, falling, no secondary).

    rising/falling/has_secondary/bg_* are evaluated whenever amp > thr, also
    where the macro would already have returned (time window, bg window); only
    `valid_p1` reproduces the macro's early-return chain.  Sums are accumulated
    sequentially in the macro's order (bit-identical to the vectorised
    `pass1_evaluate_batch`).  `falling_nonstrict` (run of wf[i+1] <= wf[i]) is
    an addition used only for the plateau diagnostic.
    """
    wf = np.asarray(wf, dtype=np.float64)
    ns = wf.size
    w = [float(v) for v in wf]
    # :400-403  main positive pulse = std::max_element (FIRST maximum)
    peak = 0
    amp = w[0]
    for i in range(1, ns):
        if w[i] > amp:
            amp = w[i]
            peak = i
    out = dict(amp=amp, peak=peak, first_thr=-1, rising=0, falling=0,
               falling_nonstrict=0, has_secondary=False, sec_max=0.0,
               bg_ok=False, bg_mean=float("nan"), bg_rms=float("nan"),
               bg_maxabs=float("nan"), roughness=float("nan"), tier=0,
               strict=False, shape_ok=False,
               time_ok=bool(peak * DT_NS >= t_lo and peak * DT_NS <= t_hi),
               valid_p1=False)
    if amp <= thr:                                   # :406-407
        return out

    # :188-199 firstThresholdCrossing: first i <= peak with wf[i] > thr
    ft = -1
    for i in range(peak + 1):
        if w[i] > thr:
            ft = i
            break
    out["first_thr"] = ft

    # :213-238 longestAscendingRun: samples above thr, strictly increasing
    # (the first_thr sample counts), run restarts at 1 after a non-increase;
    # the best run anywhere in [first_thr, peak]; 0 if peak <= first_thr.
    rising = 0
    if ft >= 0 and peak > ft:
        run = 0
        for i in range(ft, peak + 1):
            if w[i] <= thr:
                run = 0
                continue
            if i == ft or w[i] > w[i - 1]:
                run += 1
            else:
                run = 1
            if run > rising:
                rising = run
    out["rising"] = rising

    # :240-282 measurePrePulseBackground: exactly npre samples
    # [ft-npre, ft-1]; invalid if ft < npre (:248)
    if ft >= bg_samples:
        b0 = ft - bg_samples
        b1 = ft - 1
        mean = 0.0
        for i in range(b0, b1 + 1):
            mean += w[i]
        mean /= bg_samples
        s2 = 0.0
        d2 = 0.0
        largest = 0.0
        for i in range(b0, b1 + 1):
            c = w[i] - mean
            s2 += c * c
            largest = max(largest, abs(w[i]))
            if i > b0:
                d = w[i] - w[i - 1]
                d2 += d * d
        rms = math.sqrt(s2 / bg_samples)
        rough = math.sqrt(d2 / (bg_samples - 1)) if bg_samples > 1 else 0.0
        out.update(bg_ok=True, bg_mean=mean, bg_rms=rms, bg_maxabs=largest,
                   roughness=rough)
        # :284-295 isStrictBackground (all <=); :454-460 tier 2 / 1
        strict = (largest <= bg_max_abs and rms <= bg_max_rms
                  and rough <= bg_max_rough)
        out["strict"] = bool(strict)
        out["tier"] = 2 if strict else 1

    # :297-318 consecutiveDescendingAfterPeak (strict decrease, from the peak)
    fall = 0
    if 0 <= peak < ns - 1:
        for i in range(peak, ns - 1):
            if w[i + 1] < w[i]:
                fall += 1
            else:
                break
    out["falling"] = fall
    fn = 0
    if 0 <= peak < ns - 1:
        for i in range(peak, ns - 1):
            if w[i + 1] <= w[i]:
                fn += 1
            else:
                break
    out["falling_nonstrict"] = fn

    # :320-358 findSecondPeak: 3-point local max (y > prev && y >= next) with
    # |i - peak| > excl; second_peak_amp = running max (starting at 0) of
    # ALL local maxima visited up to and including the first qualifying one
    # (y > thr && y > frac*amp), where the macro returns.
    sec_max = 0.0
    has_sec = False
    for i in range(1, ns - 1):
        if abs(i - peak) <= sec_excl:
            continue
        y = w[i]
        if not (y > w[i - 1] and y >= w[i + 1]):
            continue
        sec_max = max(sec_max, y)
        if y > thr and y > sec_frac * amp:
            has_sec = True
            break
    out["has_secondary"] = has_sec
    out["sec_max"] = sec_max

    shape_ok = (out["bg_ok"] and rising >= min_rising and fall >= min_falling
                and not has_sec)
    out["shape_ok"] = bool(shape_ok)
    # :372-490: validity = amp > thr (:406), time window (:416-420),
    # first_thr >= 0 (:425-426), bg window (:437-445), rising (:468),
    # descending (:478), no second peak (:483) [max_amplitude off]
    out["valid_p1"] = bool(shape_ok and out["time_ok"])
    return out


def _first_false(ok):
    """Per row: index of the first False in ok[n, m] (m if none)."""
    m = ok.shape[1]
    f = np.argmin(ok, axis=1)
    return np.where(ok.all(axis=1), m, f)


def pass1_evaluate_batch(W, thr=P1_THR, min_rising=P1_MIN_RISING,
                         min_falling=P1_MIN_FALLING, bg_samples=P1_BG_SAMPLES,
                         bg_max_abs=P1_BG_MAX_ABS, bg_max_rms=P1_BG_MAX_RMS,
                         bg_max_rough=P1_BG_MAX_ROUGH, sec_frac=P1_SEC_FRAC,
                         sec_excl=P1_SEC_EXCL, t_lo=P1_T_LO, t_hi=P1_T_HI):
    """Vectorised `pass1_evaluate` over W[n, ns] (pedestal-subtracted); returns
    a dict of arrays with the same keys, bit-identical to the scalar version
    (the sums are accumulated column by column in the scalar's order)."""
    W = np.ascontiguousarray(W, dtype=np.float64)
    n, ns = W.shape
    rows = np.arange(n)
    idx = np.arange(ns)[None, :]
    peak = np.argmax(W, axis=1)                       # first maximum
    amp = W[rows, peak]
    live = amp > thr

    # first_thr: first i <= peak with wf > thr
    above = (W > thr) & (idx <= peak[:, None])
    ft = np.where(above.any(axis=1), np.argmax(above, axis=1), -1)
    ft = np.where(live, ft, -1)
    has_ft = ft >= 0

    # rising: per-column state machine over i in [ft, peak]
    run = np.zeros(n, dtype=np.int64)
    best = np.zeros(n, dtype=np.int64)
    for i in range(ns):
        act = has_ft & (i >= ft) & (i <= peak)
        over = W[:, i] <= thr
        inc = (ft == i) | (W[:, i] > W[:, i - 1])
        newrun = np.where(over, 0, np.where(inc, run + 1, 1))
        run = np.where(act, newrun, run)
        best = np.maximum(best, np.where(act, run, 0))
    rising = np.where(has_ft & (peak > ft), best, 0)

    # background window [ft - npre, ft - 1]
    bg_ok = has_ft & (ft >= bg_samples)
    ftc = np.where(bg_ok, ft, bg_samples)             # safe index base
    cols = (ftc - bg_samples)[:, None] + np.arange(bg_samples)[None, :]
    B = W[rows[:, None], cols]
    mean = np.zeros(n)
    for j in range(bg_samples):
        mean = mean + B[:, j]
    mean = mean / bg_samples
    s2 = np.zeros(n)
    d2 = np.zeros(n)
    largest = np.zeros(n)
    for j in range(bg_samples):
        c = B[:, j] - mean
        s2 = s2 + c * c
        largest = np.maximum(largest, np.abs(B[:, j]))
        if j > 0:
            d = B[:, j] - B[:, j - 1]
            d2 = d2 + d * d
    bg_rms = np.sqrt(s2 / bg_samples)
    rough = np.sqrt(d2 / (bg_samples - 1)) if bg_samples > 1 else np.zeros(n)
    strict = bg_ok & (largest <= bg_max_abs) & (bg_rms <= bg_max_rms) \
        & (rough <= bg_max_rough)
    tier = np.where(bg_ok, np.where(strict, 2, 1), 0)
    nan = np.full(n, np.nan)
    bg_mean = np.where(bg_ok, mean, nan)
    bg_rms = np.where(bg_ok, bg_rms, nan)
    bg_maxabs = np.where(bg_ok, largest, nan)
    rough = np.where(bg_ok, rough, nan)

    # falling: strict and non-strict descent from the peak
    def _fall(strictly):
        d = (W[:, 1:] < W[:, :-1]) if strictly else (W[:, 1:] <= W[:, :-1])
        ok = d | (np.arange(ns - 1)[None, :] < peak[:, None])
        f = _first_false(ok)                           # first non-descent i
        return np.where(live, f - peak, 0)
    falling = _fall(True)
    falling_ns = _fall(False)

    # secondary
    y = W[:, 1:-1]
    lm = (y > W[:, :-2]) & (y >= W[:, 2:])
    out_excl = np.abs(np.arange(1, ns - 1)[None, :] - peak[:, None]) > sec_excl
    vals = np.where(lm & out_excl, y, 0.0)
    qual = lm & out_excl & (y > thr) & (y > sec_frac * amp[:, None])
    has_sec = qual.any(axis=1)
    fq = np.where(has_sec, np.argmax(qual, axis=1), ns - 3)
    cm = np.maximum.accumulate(np.maximum(vals, 0.0), axis=1)
    sec_max = cm[rows, fq]

    shape_ok = live & bg_ok & (rising >= min_rising) \
        & (falling >= min_falling) & ~has_sec
    time_ok = (peak * DT_NS >= t_lo) & (peak * DT_NS <= t_hi)
    return dict(
        amp=amp, peak=peak, first_thr=ft,
        rising=np.where(live, rising, 0), falling=falling,
        falling_nonstrict=falling_ns,
        has_secondary=live & has_sec, sec_max=np.where(live, sec_max, 0.0),
        bg_ok=live & bg_ok, bg_mean=bg_mean, bg_rms=bg_rms,
        bg_maxabs=bg_maxabs, roughness=rough,
        tier=tier, strict=strict, shape_ok=shape_ok, time_ok=time_ok,
        valid_p1=shape_ok & time_ok)


def prepulse_slope(W, first_thr, npre=P1_BG_SAMPLES):
    """Least-squares slope (per sample) of the npre samples before first_thr;
    rows need first_thr >= npre."""
    rows = np.arange(W.shape[0])
    cols = (first_thr - npre)[:, None] + np.arange(npre)[None, :]
    B = W[rows[:, None], cols]
    x = np.arange(npre, dtype=np.float64) - (npre - 1) / 2.0
    return (B * x[None, :]).sum(axis=1) / (x * x).sum()


def tail_mean(W, peak, off0=TAIL_OFF0, off1=TAIL_OFF1):
    """Mean of W over samples [peak+off0, min(peak+off1, ns-1)] per row
    (pedestal-subtracted input); rows need peak + off0 <= ns - 1."""
    n, ns = W.shape
    idx = np.arange(ns)[None, :]
    m = (idx >= (peak + off0)[:, None]) & \
        (idx <= np.minimum(peak + off1, ns - 1)[:, None])
    return (W * m).sum(axis=1) / m.sum(axis=1)


def _plateau_only(P, ok_mask):
    """Rows (within ok_mask) that fail falling >= 10 ONLY because of a plateau:
    pass rising/secondary/bg window, strict falling < 10 but the non-strict
    descent (allowing an equal consecutive pair) is >= 10."""
    return (ok_mask & P["bg_ok"] & (P["rising"] >= P1_MIN_RISING)
            & ~P["has_secondary"] & (P["falling"] < P1_MIN_FALLING)
            & (P["falling_nonstrict"] >= P1_MIN_FALLING))


# ================================================================== select ==
_SPLIT_PAT = re.compile(r"prad_0*(\d+)\.evio\.(\d+)", re.IGNORECASE)


def _git_short():
    try:
        return subprocess.check_output(
            ["git", "-C", str(HERE), "rev-parse", "--short", "HEAD"],
            stderr=subprocess.DEVNULL, text=True).strip()
    except Exception:                                   # noqa: BLE001
        return "unknown"


def cut_values():
    """All cut constants as a JSON-able dict (for metadata and the PDF)."""
    g = globals()
    keys = ["CLIP_ADC", "MF_NSIGMA", "MF_MIN_AMP", "MF_KMAX", "MF_R_MODEL",
            "MF_SEC_RATIO", "MF_ONSET_MIN_NS", "FIRST_THR_MIN", "PEAK_MAX_SAMPLE",
            "A_STORE_MIN", "BAND_B0", "BAND_LO", "BAND_HI", "P1_THR",
            "P1_MIN_RISING", "P1_MIN_FALLING", "P1_BG_SAMPLES",
            "P1_BG_MAX_ABS", "P1_BG_MAX_RMS", "P1_BG_MAX_ROUGH", "P1_SEC_FRAC",
            "P1_SEC_EXCL", "PED_RMS_MAX", "PED_SLOPE_MAX", "TAIL_OFF0",
            "TAIL_OFF1", "TAIL_ABS_ADC", "TAIL_RMS_MULT", "TAIL_FRAC_A",
            "NS_EXPECT"]
    return {k: (list(g[k]) if isinstance(g[k], tuple) else g[k]) for k in keys}


class Selector:
    """Block-batched candidate selection (the plot_pileup_examples.Scanner
    pattern): rows of PbWO4 raw waveforms are stacked per BLOCK events and
    cut with numpy; the mf runs only on the survivors of every other cut."""

    def __init__(self, p, rec, per_module_cap, run, split, best_k=0,
                 mf_strict=False, energy=None, emax_mev=None, max_modules=0,
                 module_h=50.0, t0_off=None):
        self.p, self.rec, self.cap = p, rec, per_module_cap
        # opt-in extras (all off = the original first-come selection)
        self.best_k, self.mf_strict = int(best_k), bool(mf_strict)
        self.energy, self.emax_mev = energy, emax_mev   # energy(midx, t0, amp) -> MeV
        self.max_modules, self.module_h = int(max_modules), float(module_h)
        self.t0_off = t0_off                # --mf-bank: per-module bank-axis offset (ns) or None
        self.heaps = [[] for _ in range(p.hycal.module_count())]
        self.seq = 0
        self.n_no_energy = self.n_above_emax = self.n_rej_mf_strict = 0
        self.n_rows_mult_rej = 0
        self.event_counts = []              # per flush: channels with A >= module_h per event
        self.run, self.split = run, split
        self.nm = p.hycal.module_count()
        self.mods = [p.hycal.module(i) for i in range(self.nm)]
        self.names = np.array([m.name for m in self.mods])
        self.mids = np.array([m.id for m in self.mods], dtype=np.int64)
        self.overflow = int(p.wave_ana.cfg.overflow)
        self.cutflow = np.zeros((self.nm, len(CUTFLOW_NAMES)), dtype=np.int64)
        self.peak_hist = np.zeros(NS_EXPECT, dtype=np.int64)
        # cap counters per module and class (0: band b0, 1: the rest)
        self.n_stored = np.zeros((self.nm, 2), dtype=np.int64)
        self.n_wa_ne1 = 0
        self.rows, self.meta = [], []
        self.out = {k: [] for k in (
            "event", "event_number", "trigger_bits", "midx", "samples", "ped",
            "rms", "amp", "peak_sample", "first_thr", "tier", "rising",
            "falling", "sec_max", "p1_shape_ok", "slope", "tail_dev", "mf_t0",
            "mf_amp", "mf_npk", "mf_sec_ratio", "mf_t0_min", "wa_integral",
            "wa_npeaks")}
        if energy is not None:
            self.out["mf_energy"] = []
        self._best_done = False
        self.ci = {n: i for i, n in enumerate(CUTFLOW_NAMES)}

    def add(self, row, ev, evnum, trig, midx):
        self.rows.append(row)
        self.meta.append((ev, evnum, trig, midx))

    def _count(self, name, midx_rows):
        np.add.at(self.cutflow[:, self.ci[name]], midx_rows, 1)

    def flush(self):
        if not self.rows:
            return
        ci = self.ci
        Y = np.stack(self.rows)                       # uint16 [n, 100]
        meta = np.array(self.meta, dtype=np.int64)    # ev, evnum, trig, midx
        self.rows, self.meta = [], []
        pre = None
        if self.max_modules > 0:
            # --max-modules-per-event: WaveAnalyzer on ALL rows (clipped too),
            # reused below; drop every row of an event with > N channels at
            # A = raw max - ped >= module_h
            ana = self.p.wave_ana.analyze
            nr0 = Y.shape[0]
            ped_a, rms_a = np.empty(nr0), np.empty(nr0)
            npk_a = np.empty(nr0, dtype=np.int64)
            int_a = np.full(nr0, np.nan)
            for k in range(nr0):
                pd, rm, pks = ana(Y[k])
                ped_a[k], rms_a[k], npk_a[k] = pd, rm, len(pks)
                if len(pks) == 1:
                    int_a[k] = pks[0].integral
            hit = (Y.max(axis=1) - ped_a) >= self.module_h
            _, inv = np.unique(meta[:, 0], return_inverse=True)
            cnt = np.bincount(inv, weights=hit).astype(np.int64)
            self.event_counts.append(cnt)
            rej = cnt > self.max_modules
            kr = ~rej[inv]
            self.n_rows_mult_rej += int((~kr).sum())
            Y, meta = Y[kr], meta[kr]
            pre = (ped_a[kr], rms_a[kr], npk_a[kr], int_a[kr])
            if Y.shape[0] == 0:
                return
        mid = meta[:, 3]
        self._count("nseen", mid)

        # c1 clipping veto (exact)
        ok = Y.max(axis=1) < CLIP_ADC
        self._count("c1_clip_ok", mid[ok])
        sel = np.nonzero(ok)[0]
        if sel.size == 0:
            return

        # WaveAnalyzer pedestal / rms / integral on every unclipped channel
        # (cheap: ~10 us/channel) so that A is exact
        n1 = sel.size
        if pre is None:
            ana = self.p.wave_ana.analyze
            ped = np.empty(n1)
            rms = np.empty(n1)
            npk = np.empty(n1, dtype=np.int64)
            integ = np.full(n1, np.nan)
            for k, r in enumerate(sel):
                pd, rm, pks = ana(Y[r])
                ped[k], rms[k], npk[k] = pd, rm, len(pks)
                if len(pks) == 1:
                    integ[k] = pks[0].integral
        else:
            ped, rms, npk, integ = (x[sel] for x in pre)
        Y1 = Y[sel]
        mid1 = mid[sel]
        peak1 = np.argmax(Y1, axis=1)
        A1 = Y1[np.arange(n1), peak1].astype(np.float64) - ped

        # c5 amplitude >= 50; peak-sample histogram (clip + amplitude cuts)
        k5 = A1 >= A_STORE_MIN
        self._count("c5_amp_ge50", mid1[k5])
        self.peak_hist += np.bincount(peak1[k5], minlength=NS_EXPECT)
        # c4 contained
        k4 = k5 & (peak1 <= PEAK_MAX_SAMPLE)
        self._count("c4_contained", mid1[k4])
        s4 = np.nonzero(k4)[0]
        if s4.size == 0:
            return

        # PASS 1 on pedestal-subtracted floats (samples - WaveAnalyzer ped)
        Wf = Y1[s4].astype(np.float64) - ped[s4][:, None]
        P = pass1_evaluate_batch(Wf)
        m4 = mid1[s4]
        # plateau-only falling failures (informational, not in the chain)
        pl = _plateau_only(P, np.ones(s4.size, dtype=bool))
        self._count("info_fail_falling_plateau_only", m4[pl])
        k6 = (P["rising"] >= P1_MIN_RISING) & (P["falling"] >= P1_MIN_FALLING) \
            & ~P["has_secondary"]
        # PASS 1 shape (old cut 6) is only recorded, not applied (p1_shape_ok)
        k3a = P["bg_ok"] & (P["first_thr"] >= FIRST_THR_MIN)
        self._count("c3a_first_thr", m4[k3a])

        # c8 pedestal quality, c9 tail (only meaningful where bg window exists)
        slope = np.full(s4.size, np.nan)
        tdev = np.full(s4.size, np.nan)
        g = np.nonzero(k3a)[0]
        if g.size:
            slope[g] = prepulse_slope(Wf[g], P["first_thr"][g])
            tdev[g] = tail_mean(Wf[g], P["peak"][g])
        r4 = rms[s4]
        k8 = k3a & (r4 <= PED_RMS_MAX) & (np.abs(slope) <= PED_SLOPE_MAX)
        self._count("c8_pedestal", m4[k8])
        lim = np.maximum(TAIL_ABS_ADC, TAIL_RMS_MULT * r4)
        k9 = k8 & (np.abs(tdev) <= lim + TAIL_FRAC_A * P["amp"])
        self._count("c9_tail", m4[k9])
        s9 = np.nonzero(k9)[0]
        if s9.size == 0:
            return

        # mf on the survivors only: cut 2 (exactly one pulse), cut 3 (onset)
        rows9 = s4[s9]                                  # indices into Y1/sel
        Ym = Y1[rows9]
        res = self.rec.run(Ym, mid1[rows9], self.overflow)
        npk_mf = res.npk
        nr = rows9.size
        t0 = np.full(nr, np.nan)            # onset of the largest mf pulse
        amf = np.full(nr, np.nan)           # its amplitude
        sec = np.zeros(nr)                  # largest other / largest
        t0min = np.full(nr, np.nan)         # earliest mf onset
        st = np.searchsorted(res.p_ch, np.arange(nr + 1))
        for k in range(nr):
            lo, hi = st[k], st[k + 1]
            if hi == lo:
                continue
            a_k = res.p_amp[lo:hi]
            j = int(np.argmax(a_k))
            t0[k], amf[k] = res.p_t0[lo + j], a_k[j]
            t0min[k] = res.p_t0[lo:hi].min()
            if hi - lo > 1:
                sec[k] = np.delete(a_k, j).max() / a_k[j]
        m9 = mid1[rows9]
        t0_bank = t0                        # the mf bank's own onset axis (for the energy)
        if self.t0_off is not None:         # --mf-bank: onsets on the h50 axis (t0 - off)
            t0 = t0 - self.t0_off[m9]
            t0min = t0min - self.t0_off[m9]
        one = (npk_mf >= 1) & (sec < MF_SEC_RATIO)
        self._count("c2_mf_no_secondary", m9[one])
        fin = one & (t0min >= MF_ONSET_MIN_NS)
        self._count("c3b_mf_onset", m9[fin])
        f = np.nonzero(fin)[0]
        if f.size == 0:
            return

        # opt-in extras on the all-cut survivors (default: f unchanged)
        if self.mf_strict:                  # exactly one mf pulse, after c3b
            g1 = npk_mf[f] == 1
            self.n_rej_mf_strict += int((~g1).sum())
            f = f[g1]
        Ef = None
        if self.energy is not None:         # E (MeV) of the survivors only
            Ef = np.full(nr, np.nan)
            if f.size:
                Ef[f] = self.energy(m9[f], t0_bank[f], amf[f])
                fe = np.isfinite(Ef[f])
                self.n_no_energy += int((~fe).sum())
                f = f[fe]
                if self.emax_mev is not None:
                    lo_e = Ef[f] < self.emax_mev
                    self.n_above_emax += int((~lo_e).sum())
                    f = f[lo_e]
        if f.size == 0:
            return

        tier9 = P["tier"][s9]
        A9 = P["amp"][s9]
        b0 = (A9 >= BAND_B0[0]) & (A9 <= BAND_B0[1])
        t2 = tier9 == 2

        def fields(ix):                     # every output field of candidates ix
            r9 = rows9[ix]
            sk = sel[r9]                    # indices into Y/meta
            d = dict(
                event=meta[sk, 0], event_number=meta[sk, 1],
                trigger_bits=meta[sk, 2], midx=meta[sk, 3], samples=Y[sk],
                ped=ped[r9], rms=rms[r9], amp=A9[ix],
                peak_sample=P["peak"][s9][ix], first_thr=P["first_thr"][s9][ix],
                tier=tier9[ix], rising=P["rising"][s9][ix],
                falling=P["falling"][s9][ix], sec_max=P["sec_max"][s9][ix],
                p1_shape_ok=k6[s9][ix], slope=slope[s9][ix],
                tail_dev=tdev[s9][ix], mf_t0=t0[ix], mf_amp=amf[ix],
                mf_npk=npk_mf[ix], mf_sec_ratio=sec[ix], mf_t0_min=t0min[ix],
                wa_integral=integ[r9], wa_npeaks=npk[r9])
            if Ef is not None:
                d["mf_energy"] = Ef[ix]
            return d

        if self.best_k > 0:
            # per-module best-K by energy: no cap; counts of the stored set
            # (stored / info_p1_shape_ok / info_wa_npeaks_ne1) come at the end
            self._count("info_mf_npk_eq1", m9[f][npk_mf[f] == 1])
            self._count("info_tier2", m9[f][t2[f]])
            self._count("info_b0", m9[f][b0[f]])
            self._count("info_tier2_b0", m9[f][(t2 & b0)[f]])
            self._push_best(f, m9[f], Ef[f], fields)
            return

        # per-module caps (first come, separate for band b0 and the rest)
        keep = []
        for j in f:
            mi, cl = int(m9[j]), 0 if b0[j] else 1
            if self.n_stored[mi, cl] < self.cap:
                self.n_stored[mi, cl] += 1
                keep.append(j)
        keep = np.asarray(keep, dtype=np.int64)
        self._count("info_mf_npk_eq1", m9[f][npk_mf[f] == 1])
        self._count("info_p1_shape_ok", m9[keep][k6[s9][keep]])
        self._count("info_tier2", m9[f][t2[f]])
        self._count("info_b0", m9[f][b0[f]])
        self._count("info_tier2_b0", m9[f][(t2 & b0)[f]])
        if keep.size == 0:
            return
        self._count("stored", m9[keep])
        o = self.out
        d = fields(keep)
        for k, v in d.items():
            o[k].append(v)
        ne1 = d["wa_npeaks"] != 1
        self.n_wa_ne1 += int(ne1.sum())
        self._count("info_wa_npeaks_ne1", m9[keep][ne1])

    def _push_best(self, f, mf_, Ee, fields):
        """Per-module bounded min-heap of size best_k keyed by (E, -seq): on
        equal E the earlier candidate wins.  Payload = one tuple of every
        output field (the samples row copied)."""
        K, heaps = self.best_k, self.heaps
        # a candidate not above the root of a full heap can never enter (the
        # root only grows), so only the others are materialised
        adm = np.array([len(heaps[m]) < K or e > heaps[m][0][0]
                        for m, e in zip(mf_.tolist(), Ee.tolist())], dtype=bool)
        seqs = self.seq + np.nonzero(adm)[0]
        self.seq += f.size
        if not adm.any():
            return
        d = fields(f[adm])
        self._pkeys = list(d)
        cols = [d[k] for k in self._pkeys]
        si = self._pkeys.index("samples")
        for i, (m, e) in enumerate(zip(mf_[adm].tolist(), Ee[adm].tolist())):
            row = [c[i] for c in cols]
            row[si] = row[si].copy()
            h = heaps[m]
            (heapq.heappush if len(h) < K else heapq.heappushpop)(
                h, (e, -int(seqs[i]), tuple(row)))

    def kept(self):
        if self.best_k > 0:
            return sum(len(h) for h in self.heaps)
        return int(self.n_stored.sum())

    def _finalize_best(self):
        """Heap contents in event order -> self.out; stored-set cutflow."""
        self._best_done = True
        items = sorted((it for h in self.heaps for it in h), key=lambda t: -t[1])
        if not items:
            return
        o = self.out
        for j, k in enumerate(self._pkeys):
            col = [it[2][j] for it in items]
            o[k].append(np.stack(col) if k == "samples" else np.array(col))
        mi, ci, cf = o["midx"][0], self.ci, self.cutflow
        cf[:, ci["stored"]] = np.bincount(mi, minlength=self.nm)
        cf[:, ci["info_p1_shape_ok"]] = np.bincount(
            mi[o["p1_shape_ok"][0]], minlength=self.nm)
        ne1 = o["wa_npeaks"][0] != 1
        cf[:, ci["info_wa_npeaks_ne1"]] = np.bincount(mi[ne1], minlength=self.nm)
        self.n_wa_ne1 = int(ne1.sum())

    def arrays(self):
        if self.best_k > 0 and not self._best_done:
            self._finalize_best()
        o = self.out
        dt = dict(event=np.int64, event_number=np.int64, trigger_bits=np.int64,
                  midx=np.int64, samples=np.uint16, ped=np.float64,
                  rms=np.float64, amp=np.float64, peak_sample=np.int16,
                  first_thr=np.int16, tier=np.int8, rising=np.int16,
                  falling=np.int16, sec_max=np.float64, p1_shape_ok=np.bool_,
                  slope=np.float64,
                  tail_dev=np.float64, mf_t0=np.float64, mf_amp=np.float64,
                  mf_npk=np.int16, mf_sec_ratio=np.float64,
                  mf_t0_min=np.float64, wa_integral=np.float64,
                  wa_npeaks=np.int16)
        if self.energy is not None:
            dt["mf_energy"] = np.float64
        res = {}
        for k, d in dt.items():
            if o[k]:
                res[k] = np.concatenate(o[k]).astype(d)
            else:
                res[k] = np.zeros((0, NS_EXPECT) if k == "samples" else (0,), d)
        mi = res["midx"]
        n = mi.size
        res["name"] = self.names[mi] if n else np.array([], dtype="<U8")
        res["module_id"] = self.mids[mi] if n else np.zeros(0, np.int64)
        A = res["amp"]
        res["band_b0"] = (A >= BAND_B0[0]) & (A <= BAND_B0[1])
        res["band_lo"] = (A >= BAND_LO[0]) & (A <= BAND_LO[1])
        res["band_hi"] = (A >= BAND_HI[0]) & (A <= BAND_HI[1])
        return res


# ---------------------------------------------- Wassim single-pulse reference --
# Opt-in `select --wassim-ref` (and `build --single-pulse`): the single best PASS 1
# candidate per module of a split (Waveform_Wassim/best_ref_shape_fit.C, the macro's
# betterCandidate order), nothing else.  Selector above is not involved.
WASSIM_CUTFLOW = [
    "nseen",              # PbWO4, nsamples == 100, event passed the trigger gate
    "c1_clip_ok",         # no raw sample >= 4095
    "c_time_ok",          # raw argmax (== argmax of W) in the macro window 150-170 ns
    "c_amp_gt_thr",       # A = max(samples - ped) > 2 ADC (counted among time_ok rows)
    "c_valid_p1",         # macro Candidate.valid (pass1_evaluate_batch valid_p1)
    "info_tier2_valid",   # of the valid rows: tier 2 (strict background)
    "stored",             # 1 if the module has a best candidate in this split
]
WASSIM_N_SEQ = 5
WASSIM_CUT_ORDER_TEXT = (
    "nseen > clip veto (no raw sample >= 4095) > raw argmax in 150-170 ns (same row as "
    "argmax of samples - ped, evaluated before the pedestal so that WaveAnalyzer runs "
    "only on these rows) > A = max(samples - ped) > 2 ADC > valid_p1 = pass1_evaluate_batch "
    "(macro defaults, incl. rising >= 4, falling >= 10, no secondary > 15 %, 20-sample "
    "pre-pulse window) > per module only the best valid candidate of the split by the "
    "macro's betterCandidate order (tier desc, amp desc, roughness asc, bg_rms asc; "
    "strict comparisons: first come on exact ties).  No mf, no c3a/c8/c9/c2/c3b, no "
    "A >= 50, no containment, no multiplicity / energy extras.")
WASSIM_NEW_FIELDS = ("bg_mean", "bg_rms", "bg_maxabs", "roughness")


def _wassim_better(new, old):
    """Macro betterCandidate: tier desc, amp desc, roughness asc, bg_rms asc (strict)."""
    for x, y, hi in zip(new, old, (True, True, False, False)):
        if x != y:
            return (x > y) if hi else (x < y)
    return False


class WassimSelector:
    """Block-batched `select --wassim-ref`: rows of PbWO4 raw waveforms are stacked per
    BLOCK events; per module only the best valid_p1 candidate of the split is kept."""

    def __init__(self, p, run, split):
        self.p, self.run, self.split = p, run, split
        self.nm = p.hycal.module_count()
        self.mods = [p.hycal.module(i) for i in range(self.nm)]
        self.names = np.array([m.name for m in self.mods])
        self.mids = np.array([m.id for m in self.mods], dtype=np.int64)
        self.cutflow = np.zeros((self.nm, len(WASSIM_CUTFLOW)), dtype=np.int64)
        self.peak_hist = np.zeros(NS_EXPECT, dtype=np.int64)
        self.ci = {n: i for i, n in enumerate(WASSIM_CUTFLOW)}
        self.rows, self.meta = [], []
        self.best = {}                      # midx -> (key tuple, payload dict)

    def add(self, row, ev, evnum, trig, midx):
        self.rows.append(row)
        self.meta.append((ev, evnum, trig, midx))

    def _count(self, name, midx_rows):
        np.add.at(self.cutflow[:, self.ci[name]], midx_rows, 1)

    def kept(self):
        return len(self.best)

    def flush(self):
        if not self.rows:
            return
        Y = np.stack(self.rows)                       # uint16 [n, 100]
        meta = np.array(self.meta, dtype=np.int64)    # ev, evnum, trig, midx
        self.rows, self.meta = [], []
        mid = meta[:, 3]
        self._count("nseen", mid)
        ok = Y.max(axis=1) < CLIP_ADC
        self._count("c1_clip_ok", mid[ok])
        sel = np.nonzero(ok)[0]
        if sel.size == 0:
            return
        Y1 = Y[sel]
        peak1 = np.argmax(Y1, axis=1)
        # peak_hist: unclipped rows with (max - mean of the first 20 samples) >= 50
        # (the WaveAnalyzer pedestal is only evaluated on the time_ok rows)
        px = Y1.max(axis=1).astype(np.float64) - Y1[:, :P1_BG_SAMPLES].mean(axis=1)
        self.peak_hist += np.bincount(peak1[px >= A_STORE_MIN], minlength=NS_EXPECT)
        tok = (peak1 * DT_NS >= P1_T_LO) & (peak1 * DT_NS <= P1_T_HI)
        self._count("c_time_ok", mid[sel][tok])
        s = np.nonzero(tok)[0]
        if s.size == 0:
            return
        ana = self.p.wave_ana.analyze
        ns_ = s.size
        ped, rms = np.empty(ns_), np.empty(ns_)
        npk = np.empty(ns_, dtype=np.int64)
        integ = np.full(ns_, np.nan)
        Ys = Y1[s]
        for k in range(ns_):
            pd, rm, pks = ana(Ys[k])
            ped[k], rms[k], npk[k] = pd, rm, len(pks)
            if len(pks) == 1:
                integ[k] = pks[0].integral
        W = Ys.astype(np.float64) - ped[:, None]
        P = pass1_evaluate_batch(W)
        msel = mid[sel][s]
        self._count("c_amp_gt_thr", msel[P["amp"] > P1_THR])
        v = np.nonzero(P["valid_p1"])[0]
        self._count("c_valid_p1", msel[v])
        if v.size == 0:
            return
        self._count("info_tier2_valid", msel[v][P["tier"][v] == 2])
        mv = msel[v]
        o = np.lexsort((v, P["bg_rms"][v], P["roughness"][v], -P["amp"][v],
                        -P["tier"][v], mv))
        vs = v[o]
        mo = mv[o]
        first = np.concatenate([[True], mo[1:] != mo[:-1]])
        for j in vs[first]:
            mi = int(msel[j])
            key = (int(P["tier"][j]), float(P["amp"][j]), float(P["roughness"][j]),
                   float(P["bg_rms"][j]))
            cur = self.best.get(mi)
            if cur is not None and not _wassim_better(key, cur[0]):
                continue
            sk = sel[s[j]]
            self.best[mi] = (key, dict(
                event=meta[sk, 0], event_number=meta[sk, 1], trigger_bits=meta[sk, 2],
                midx=meta[sk, 3], samples=Y[sk].copy(), ped=ped[j], rms=rms[j],
                amp=P["amp"][j], peak_sample=P["peak"][j], first_thr=P["first_thr"][j],
                tier=P["tier"][j], rising=P["rising"][j], falling=P["falling"][j],
                sec_max=P["sec_max"][j], p1_shape_ok=P["shape_ok"][j],
                bg_mean=P["bg_mean"][j], bg_rms=P["bg_rms"][j],
                bg_maxabs=P["bg_maxabs"][j], roughness=P["roughness"][j],
                wa_integral=integ[j], wa_npeaks=npk[j]))

    def arrays(self):
        dt = dict(event=np.int64, event_number=np.int64, trigger_bits=np.int64,
                  midx=np.int64, samples=np.uint16, ped=np.float64, rms=np.float64,
                  amp=np.float64, peak_sample=np.int16, first_thr=np.int16,
                  tier=np.int8, rising=np.int16, falling=np.int16, sec_max=np.float64,
                  p1_shape_ok=np.bool_, bg_mean=np.float64, bg_rms=np.float64,
                  bg_maxabs=np.float64, roughness=np.float64,
                  wa_integral=np.float64, wa_npeaks=np.int16)
        items = [self.best[m][1] for m in sorted(self.best)]
        n = len(items)
        res = {}
        for k, d in dt.items():
            if n:
                res[k] = np.array([it[k] for it in items]).astype(d)
            else:
                res[k] = np.zeros((0, NS_EXPECT) if k == "samples" else (0,), d)
        # not computed in this mode
        for k in ("slope", "tail_dev", "mf_t0", "mf_amp", "mf_sec_ratio", "mf_t0_min"):
            res[k] = np.full(n, np.nan)
        res["mf_npk"] = np.full(n, -1, dtype=np.int16)
        mi = res["midx"]
        res["name"] = self.names[mi] if n else np.array([], dtype="<U8")
        res["module_id"] = self.mids[mi] if n else np.zeros(0, np.int64)
        A = res["amp"]
        res["band_b0"] = (A >= BAND_B0[0]) & (A <= BAND_B0[1])
        res["band_lo"] = (A >= BAND_LO[0]) & (A <= BAND_LO[1])
        res["band_hi"] = (A >= BAND_HI[0]) & (A <= BAND_HI[1])
        if n:
            self.cutflow[mi, self.ci["stored"]] = 1
        return res


def _wassim_select_finish(a, sel, arr, mods, wcfg, wave_cfg, run, split, evio, out_path,
                          trig_rule, n_seen, n_trig, wall):
    """Meta, npz and terminal summary of `select --wassim-ref`."""
    pm = np.array([m.is_pwo4() for m in mods])
    cf = sel.cutflow[pm].sum(axis=0)
    cut_tot = {n: int(v) for n, v in zip(WASSIM_CUTFLOW, cf)}
    st = os.stat(evio)
    meta = dict(
        schema=SCHEMA, mode="wassim_ref", evio=os.path.abspath(evio),
        evio_size=st.st_size, evio_mtime=datetime.fromtimestamp(st.st_mtime).isoformat(),
        run=run, split=split, trigger=a.trigger, trigger_rule=trig_rule,
        max_events=a.max_events,
        wassim=dict(
            ref="Waveform_Wassim/best_ref_shape_fit.C (PASS 1 evaluateCandidate / "
                "betterCandidate), Hycal_prad_wf_analysis.pdf section 4",
            params=dict(thr=P1_THR, min_rising=P1_MIN_RISING, min_falling=P1_MIN_FALLING,
                        bg_samples=P1_BG_SAMPLES, bg_max_abs=P1_BG_MAX_ABS,
                        bg_max_rms=P1_BG_MAX_RMS, bg_max_rough=P1_BG_MAX_ROUGH,
                        sec_frac=P1_SEC_FRAC, sec_excl=P1_SEC_EXCL,
                        t_window_ns=[P1_T_LO, P1_T_HI], clip_adc=CLIP_ADC),
            better_candidate="tier desc, amp desc, roughness asc, bg_rms asc; strict "
                             "comparisons (first come on exact ties); one best valid "
                             "candidate per module and split",
            cutflow_total=cut_tot,
            cutflow_note="sequential columns nseen > c1_clip_ok > c_time_ok (raw argmax "
                         "in the window) > c_amp_gt_thr (counted among the time_ok rows; "
                         "the pedestal is only evaluated there) > c_valid_p1; "
                         "info_tier2_valid = tier 2 among the valid rows; stored = 1 per "
                         "module with a candidate; per module in the npz cutflow"),
        n_events_seen=n_seen, n_events_after_trigger=n_trig,
        n_candidates=int(arr["midx"].size),
        cut_order=WASSIM_CUT_ORDER_TEXT, cut_values=cut_values(),
        cutflow_names=WASSIM_CUTFLOW, n_cutflow_sequential=WASSIM_N_SEQ,
        peak_hist_note="raw argmax of unclipped channels with (max - mean of the first 20 "
                       "samples) >= 50 (proxy pedestal; WaveAnalyzer runs only on the "
                       "150-170 ns rows), no amplitude / containment cut",
        wave_config=dict(
            wave_cfg, ped_nsamples=int(wcfg.ped_nsamples),
            overflow=int(wcfg.overflow), clk_mhz=float(wcfg.clk_mhz),
            peak_nsigma=float(wcfg.peak_nsigma),
            min_peak_height=float(wcfg.min_peak_height),
            ped_flatness=float(wcfg.ped_flatness),
            ped_max_iter=int(wcfg.ped_max_iter)),
        amplitude_def="A = max(samples - WaveAnalyzer pedestal)",
        pedestal_def="WaveAnalyzer pedestal / rms (wave_ana.analyze); W = samples - ped; "
                     "bg_mean = mean of the 20 samples before first_thr (stored)",
        peak_def="first argmax of the samples (macro max_element)",
        not_computed="slope, tail_dev, mf_* (NaN, mf_npk -1); wa_* are the WaveAnalyzer "
                     "single-peak integral / peak count of the stored candidate",
        wall_s=round(wall, 1), events_per_s=round(n_seen / wall, 1),
        git=_git_short(), created=datetime.now().isoformat(timespec="seconds"),
        python=sys.version.split()[0], numpy=np.__version__, argv=sys.argv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_name(out_path.name + ".tmp")
    with open(tmp, "wb") as fh:
        np.savez_compressed(
            fh, run=np.int64(run), split=np.int64(split), **arr,
            cutflow_names=np.array(WASSIM_CUTFLOW), cutflow=sel.cutflow,
            peak_hist=sel.peak_hist, meta=np.array(json.dumps(meta)))
    os.replace(tmp, out_path)
    print(f"[done] {out_path}  {n_seen} events ({n_trig} after trigger) in "
          f"{wall:.0f} s = {n_seen / wall:.0f} ev/s  ({wall / max(n_seen, 1):.4f} s/event)")
    print("[cutflow] (Wassim single-pulse reference, PbWO4 totals)")
    for nme, v in cut_tot.items():
        print(f"  {nme:34s} {v:>10d}")
    t2 = int((arr["tier"] == 2).sum())
    print(f"[stored] {arr['midx'].size} reference candidates (one per module), "
          f"{t2} tier-2, {int(pm.sum())} PbWO4 modules", flush=True)
    return 0


def _trigger_gate(spec, bits_file):
    """Returns (callable(fadc_evt) -> bool, rule text)."""
    if spec == "physics":
        return (lambda fe: int(fe.info.trigger_bits) == PHYSICS_BITS_EXACT,
                f"trigger_bits == 0x{PHYSICS_BITS_EXACT:x} (exact; "
                "C.PHYSICS_TRIGGER_BITS, SSP0 only)")
    if spec == "physics-any":
        return (lambda fe: (int(fe.info.trigger_bits) & PHYSICS_MASK_ANY) != 0,
                f"(trigger_bits & 0x{PHYSICS_MASK_ANY:x}) != 0 "
                "(any physics trigger; C.PHYSICS_TRIGGER_MASK)")
    if spec == "rawsum":
        from fit_pulse_template import _resolve_trigger_event_types
        types = _resolve_trigger_event_types("SSP_RawSum", bits_file)
        return (lambda fe: int(fe.info.trigger_type) in types,
                "trigger_type in {" + ",".join(f"0x{t:02x}" for t in
                                               sorted(types)) + "} (SSP_RawSum)")
    return (lambda fe: True, "all events (no trigger gate)")


def _mf_bank_offsets(path, nmod):
    """Per-module offset (ns) from the --mf-bank onset axis to the h50 onset axis
    (t0_h50 = t0_bank - off), the convention of check's `t0_seed_sampled = mf_t0 +
    onset_shift_ns`: onset_shift_ns of the module's row (fallback 0/1); a fallback-2
    row and a module absent from the bank take d_typ = median over the fallback-0/1
    rows (+ the row's own shift for fallback 2)."""
    z = np.load(path, allow_pickle=False)
    sh = z["onset_shift_ns"].astype(np.float64)
    own = z["fallback"] <= 1
    d_typ = float(np.median(sh[own])) if own.any() else 0.0
    off = np.full(nmod, d_typ)
    for r, m in enumerate(z["module_index"]):
        if 0 <= int(m) < nmod:
            off[int(m)] = sh[r] if own[r] else d_typ + sh[r]
    return off


def cmd_select(a):
    import _common as C
    from _common import dec
    import mf_pulse

    t_start = time.monotonic()
    evio = str(a.evio)
    mt = _SPLIT_PAT.search(os.path.basename(evio))
    if not mt:
        sys.exit(f"[ERROR] cannot read run/split from {evio} "
                 "(expected prad_0RRRRR.evio.SSSSS)")
    run, split = int(mt.group(1)), int(mt.group(2))
    out_path = Path(a.out)
    if out_path.suffix != ".npz":
        sys.exit("[ERROR] --out must end in .npz")
    wass = bool(a.wassim_ref)
    if wass:
        bad = [n for n, on in (
            ("--per-module-best", a.per_module_best != 0),
            ("--mf-strict", a.mf_strict), ("--max-modules-per-event", a.max_modules_per_event != 0),
            ("--emax-frac", a.emax_frac is not None),
            ("--exclude-trigger-bits", a.exclude_trigger_bits != 0),
            ("--mf-bank", bool(a.mf_bank))) if on]
        if bad:
            sys.exit(f"[ERROR] --wassim-ref cannot be combined with {', '.join(bad)}")

    p = C.setup_pipeline(evio_path=evio, run_num=run)
    mods = [p.hycal.module(i) for i in range(p.hycal.module_count())]
    pwo4 = [bool(m.is_pwo4()) for m in mods]
    wcfg = p.wave_ana.cfg
    wave_cfg = dict(smooth_order=int(wcfg.smooth_order),
                    int_tail_ratio=float(wcfg.int_tail_ratio),
                    tail_break_n=int(wcfg.tail_break_n))
    h50 = Path(a.h50)
    if not h50.is_absolute():
        h50 = HERE / h50
    column = mf_pulse.resolve_column(str(h50), "auto", prefer="25320")
    mf_bank_path = mf_bank_sha1 = mf_bank_kind = ""
    t0_off = None
    if wass:                                    # no matched filter in this mode
        bank = info = None
    elif a.mf_bank:
        mf_bank_path = os.path.abspath(a.mf_bank)
        bank, info = mf_pulse.build_bank([m.name for m in mods], pwo4, mf_bank_path,
                                         wave_cfg)
        mf_bank_sha1, mf_bank_kind = info["sha1"], info["kind"]
        t0_off = _mf_bank_offsets(mf_bank_path, len(mods))
    else:
        bank, info = mf_pulse.build_bank([m.name for m in mods], pwo4, str(h50),
                                         wave_cfg, column=column)
    rec = None if wass else mf_pulse.MFReconstructor(
        bank, nsigma=MF_NSIGMA, min_amp=MF_MIN_AMP, kmax=MF_KMAX,
        r_model=MF_R_MODEL, ped_nsamples=int(wcfg.ped_nsamples))

    bits_file = Path(os.environ.get("PRAD2_DATABASE_DIR", "")) / "trigger_bits.json"
    gate, trig_rule = _trigger_gate(a.trigger, bits_file)
    if wass:
        print(f"[setup] run {run} split {split}  trigger: {trig_rule}\n"
              "[setup] --wassim-ref: single best valid_p1 candidate per module "
              "(best_ref_shape_fit.C betterCandidate order); no mf, no cuts but the clip "
              "veto and PASS 1", flush=True)
    else:
        print(f"[setup] run {run} split {split}  trigger: {trig_rule}\n"
              + (f"[setup] mf: h50 {h50.name} column {column}" if not a.mf_bank else
                 f"[setup] mf: --mf-bank {mf_bank_path} ({mf_bank_kind}, sha1 "
                 f"{mf_bank_sha1[:12]}, ampcorr {info['ampcorr']})")
              + f"  nsigma {MF_NSIGMA:g} "
              f"min_amp {MF_MIN_AMP:g} kmax {MF_KMAX} r_model {MF_R_MODEL:g}",
              flush=True)

    extras = (a.per_module_best > 0 or a.emax_frac is not None or a.mf_strict
              or a.max_modules_per_event > 0 or a.exclude_trigger_bits != 0
              or bool(a.mf_bank))
    energy = emax_mev = ebeam = None
    calib_path = calib_sha1 = ""
    if a.per_module_best < 0 or a.max_modules_per_event < 0:
        sys.exit("[ERROR] --per-module-best and --max-modules-per-event must be >= 0")
    if a.emax_frac is not None and not a.emax_frac > 0:
        sys.exit("[ERROR] --emax-frac must be > 0")
    if a.per_module_best > 0 or a.emax_frac is not None:
        calib_path = a.calib_file
        if not calib_path:
            import module_timing_dump as mtd
            if run not in mtd.RUN_CFG:
                sys.exit(f"[ERROR] run {run} is not in module_timing_dump.RUN_CFG: "
                         "pass --calib-file")
            calib_path = os.path.join(mtd.CALIB_DIR, mtd.RUN_CFG[run][1])
        cal = np.full(len(mods), np.nan)
        by_name = {m["name"]: float(m.get("factor", 0.0))
                   for m in json.load(open(calib_path))}
        for i, m in enumerate(mods):
            fc = by_name.get(m.name, 0.0)
            if fc > 0:
                cal[i] = fc
        calib_sha1 = _sha1(calib_path)
        ebeam = float(p.geo.Ebeam)
        if a.emax_frac is not None:
            emax_mev = a.emax_frac * ebeam
        energy = lambda mi_, t0_, amp_: (                    # noqa: E731
            mf_pulse.pulse_energy_units(bank, mi_, t0_, amp_, kappa=1.0) * cal[mi_])
        print(f"[setup] energy E = mf_amp * I_int * cal (MeV): calib {calib_path}, "
              f"E_beam {ebeam:.1f} MeV" + ("" if emax_mev is None else
              f", E < {a.emax_frac:g} * E_beam = {emax_mev:.1f} MeV"), flush=True)

    sel = WassimSelector(p, run, split) if wass else Selector(
        p, rec, a.per_module_cap, run, split, best_k=a.per_module_best,
                   mf_strict=a.mf_strict, energy=energy, emax_mev=emax_mev,
                   max_modules=a.max_modules_per_event,
                   module_h=a.module_count_height_min, t0_off=t0_off)
    ch = dec.EvChannel()
    ch.set_config(p.cfg)
    if ch.open_auto(evio) != dec.Status.success:
        sys.exit(f"[ERROR] cannot open {evio}")

    cache = {}                                  # (roc tag, slot, ch) -> midx|-1
    n_seen = n_trig = n_blk = n_excl_bits = 0
    done = False
    next_prog = 20000
    while not done and ch.read() == dec.Status.success:
        if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
            continue
        for i in range(ch.get_n_events()):
            d = ch.decode_event(i, with_ssp=False)
            if not d["ok"]:
                continue
            ev_idx = n_seen
            n_seen += 1
            fe = d["event"]
            if a.max_events and n_seen >= a.max_events:
                done = True
            if not gate(fe):
                if done:
                    break
                continue
            n_trig += 1
            evnum = int(fe.info.event_number)
            trig = int(fe.info.trigger_bits)
            if a.exclude_trigger_bits and (trig & a.exclude_trigger_bits):
                n_excl_bits += 1
                if done:
                    break
                continue
            for ri in range(fe.nrocs):
                roc = fe.roc(ri)
                if not roc.present:
                    continue
                crate = p.crate_map.get(roc.tag)
                if crate is None:
                    continue
                for s in roc.present_slots():
                    slot = roc.slot(s)
                    for c in slot.present_channels():
                        key = (roc.tag, s, c)
                        mi = cache.get(key)
                        if mi is None:
                            m = p.hycal.module_by_daq(crate, s, c)
                            mi = (m.index if (m is not None and m.is_hycal()
                                              and m.is_pwo4()) else -1)
                            cache[key] = mi
                        if mi < 0:
                            continue
                        cd = slot.channel(c)
                        if cd.nsamples != NS_EXPECT:
                            continue
                        sel.add(np.asarray(cd.samples, dtype=np.uint16),
                                ev_idx, evnum, trig, mi)
            n_blk += 1
            if n_blk >= BLOCK:
                sel.flush()
                n_blk = 0
            if n_seen >= next_prog:
                el = time.monotonic() - t_start
                print(f"  [progress] {n_seen} events ({n_trig} gated)  stored "
                      f"{sel.kept()}  {n_seen / el:.0f} ev/s  "
                      f"{el:.0f} s", flush=True)
                next_prog += 20000
            if done:
                break
    sel.flush()
    wall = time.monotonic() - t_start

    arr = sel.arrays()
    if not a.max_events and n_seen == 0:
        sys.exit("[ERROR] no events decoded")
    if wass:
        return _wassim_select_finish(a, sel, arr, mods, wcfg, wave_cfg, run, split, evio,
                                     out_path, trig_rule, n_seen, n_trig, wall)
    n_analyzed = n_trig - n_excl_bits
    ex = None
    if extras:
        ex = dict(
            mode=(f"per_module_best K={a.per_module_best}, ranked by mf_energy"
                  if a.per_module_best > 0 else "per_module_cap (first come)"),
            per_module_best=a.per_module_best, emax_frac=a.emax_frac,
            mf_strict=bool(a.mf_strict),
            max_modules_per_event=a.max_modules_per_event,
            module_count_height_min=a.module_count_height_min,
            exclude_trigger_bits=a.exclude_trigger_bits,
            n_excl_bits=n_excl_bits, n_events_analyzed=n_analyzed,
            n_rej_mf_strict=sel.n_rej_mf_strict,
            **({} if not a.mf_bank else dict(
                mf_bank_ampcorr=info["ampcorr"],
                mf_t0_axis="h50 onset axis: mf_t0 / mf_t0_min = bank onset - "
                           "onset_shift_ns of the module (fallback-2 rows: + the "
                           "median of the fallback-0/1 rows; modules absent from the "
                           "bank: that median), as check does; the energy uses the "
                           "bank's own axis")))
        if energy is not None:
            ex.update(
                energy_def="mf_energy = mf_amp * I_int(" + ("h50 bank" if not a.mf_bank
                                                           else "--mf-bank (x f(A) if it carries ampcorr)")
                           + (", mf_t0" if not a.mf_bank else ", onset on the bank axis = mf_t0 + module offset")
                           + ") * cal_factor "
                           "(MeV; mf_pulse.pulse_energy_units kappa=1 times the "
                           "per-module `factor` of calib_file; factor <= 0 or "
                           "missing -> NaN, dropped)",
                calib_file=os.path.abspath(calib_path), calib_sha1=calib_sha1,
                e_beam_mev=ebeam, e_max_mev=emax_mev,
                n_no_energy=sel.n_no_energy, n_above_emax=sel.n_above_emax)
        if a.max_modules_per_event > 0:
            nev = np.concatenate(sel.event_counts + [np.zeros(
                max(n_analyzed - sum(c.size for c in sel.event_counts), 0),
                dtype=np.int64)])
            n_rej_ev = int((nev > a.max_modules_per_event).sum())
            ex.update(
                mult_mean=float(nev.mean()) if nev.size else float("nan"),
                mult_p99=float(np.percentile(nev, 99)) if nev.size else float("nan"),
                n_events_mult_rejected=n_rej_ev,
                frac_events_mult_rejected=(n_rej_ev / nev.size if nev.size
                                           else float("nan")),
                n_rows_mult_rejected=sel.n_rows_mult_rej)
    st = os.stat(evio)
    meta = dict(
        schema=SCHEMA, evio=os.path.abspath(evio), evio_size=st.st_size,
        evio_mtime=datetime.fromtimestamp(st.st_mtime).isoformat(),
        run=run, split=split, trigger=a.trigger, trigger_rule=trig_rule,
        max_events=a.max_events, per_module_cap=a.per_module_cap,
        per_module_cap_rule=("not applied: per-module best K by energy"
                             if a.per_module_best > 0 else
                             "cap applied separately to band-b0 candidates and "
                             "to all other candidates, first come in event order"),
        p1_shape_def="p1_shape_ok = rising>=4 & falling>=10 (strict) & no "
                     "secondary >15% outside +-6 (old cut 6; stored, not applied)",
        n_events_seen=n_seen, n_events_after_trigger=n_trig,
        n_candidates=int(arr["midx"].size),
        n_stored_wa_npeaks_ne1=sel.n_wa_ne1,
        cut_order=CUT_ORDER_TEXT, cut_values=cut_values(),
        cutflow_names=CUTFLOW_NAMES, n_cutflow_sequential=len(CUTFLOW_SEQ),
        cutflow_note=(
            "sequential columns count channels passing that cut and all "
            "earlier ones; info_* columns are counted on their stated subset "
            "and are not part of the chain.  info_fail_falling_plateau_only: "
            "channels passing clip, A>=50 and contained, with a bg window, "
            "rising>=4, no secondary, strict falling<10 but non-strict "
            "(equal pair allowed) falling>=10"),
        peak_hist_note=("raw argmax sample of channels passing clip + "
                        "A>=50 (before the containment cut)"),
        wave_config=dict(
            wave_cfg, ped_nsamples=int(wcfg.ped_nsamples),
            overflow=int(wcfg.overflow), clk_mhz=float(wcfg.clk_mhz),
            peak_nsigma=float(wcfg.peak_nsigma),
            min_peak_height=float(wcfg.min_peak_height),
            ped_flatness=float(wcfg.ped_flatness),
            ped_max_iter=int(wcfg.ped_max_iter)),
        amplitude_def="A = raw peak sample - WaveAnalyzer pedestal",
        pedestal_def="WaveAnalyzer pedestal / rms (wave_ana.analyze)",
        peak_def="first argmax of the raw samples (PASS 1 max_element)",
        tail_def=("|mean(samples - ped)| over [peak+40, min(peak+50, 99)] <= "
                  "max(2, 2*rms) + 0.04*A (the plan's 0.5 % of A rejects "
                  "every real pulse: the tail is ~2 % of A there); tail_dev "
                  "stores the signed mean"),
        slope_def="least-squares slope of the 20 samples [first_thr-20, "
                  "first_thr-1]",
        mf=dict(**(dict(h50=str(h50), column=column) if not a.mf_bank else
                   dict(h50=mf_bank_path, sha1=mf_bank_sha1, kind=mf_bank_kind)),
                nsigma=MF_NSIGMA, min_amp=MF_MIN_AMP, kmax=MF_KMAX, r_model=MF_R_MODEL,
                sec_ratio_used=MF_SEC_RATIO,
                onset_def="mf_t0 / mf_amp = onset (ns) / amplitude of the "
                          "largest mf pulse; mf_t0_min = earliest mf onset; "
                          "mf_npk = number of mf pulses; mf_sec_ratio = "
                          "largest other / largest amplitude",
                use="pulse count, onset and (t0, A) seed only; never a chi2 "
                    "selection"),
        **({} if ex is None else dict(selection_extras=ex)),
        wall_s=round(wall, 1), events_per_s=round(n_seen / wall, 1),
        git=_git_short(), created=datetime.now().isoformat(timespec="seconds"),
        python=sys.version.split()[0], numpy=np.__version__,
        argv=sys.argv)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_name(out_path.name + ".tmp")
    with open(tmp, "wb") as fh:
        np.savez_compressed(
            fh, run=np.int64(run), split=np.int64(split), **arr,
            cutflow_names=np.array(CUTFLOW_NAMES), cutflow=sel.cutflow,
            peak_hist=sel.peak_hist, meta=np.array(json.dumps(meta)))
    os.replace(tmp, out_path)

    # ------------------------------------------------------------ summary --
    pm = np.array([m.is_pwo4() for m in mods])
    cf = sel.cutflow[pm].sum(axis=0)
    print(f"[done] {out_path}  {n_seen} events ({n_trig} after trigger) in "
          f"{wall:.0f} s = {n_seen / wall:.0f} ev/s")
    print("[cutflow] (PbWO4 totals, in applied order)")
    for nme, v in zip(CUTFLOW_NAMES, cf):
        print(f"  {nme:34s} {int(v):>10d}")
    ct2 = (arr["tier"] == 2) & arr["band_b0"]
    cb0 = np.bincount(arr["midx"][arr["band_b0"]], minlength=len(mods))
    cnt = np.bincount(arr["midx"][ct2], minlength=len(mods))
    print(f"[stored] {arr['midx'].size} candidates, {int(arr['band_b0'].sum())} b0 "
          f"({int(ct2.sum())} tier-2 b0, {int(arr['p1_shape_ok'].sum())} p1_shape_ok); "
          f"modules with >= 20 b0: {int((cb0 >= 20).sum())}, >= 50 b0: "
          f"{int((cb0 >= 50).sum())}, >= 50 tier-2 b0: {int((cnt >= 50).sum())}; "
          f"peak-sample histogram mode = {int(sel.peak_hist.argmax())} "
          f"({int(sel.peak_hist.max())} ch)", flush=True)
    if ex is not None:
        print(f"[extras] mode: {ex['mode']}")
        if a.exclude_trigger_bits:
            print(f"[extras] exclude-trigger-bits 0x{a.exclude_trigger_bits:x}: "
                  f"{n_excl_bits} of {n_trig} gated events rejected "
                  f"({n_analyzed} analysed)")
        if a.max_modules_per_event > 0:
            print(f"[extras] max-modules-per-event {a.max_modules_per_event} "
                  f"(A >= {a.module_count_height_min:g} ADC): per-event count mean "
                  f"{ex['mult_mean']:.2f}, 99th pct {ex['mult_p99']:.1f}; "
                  f"{ex['n_events_mult_rejected']} of {n_analyzed} events rejected "
                  f"({100 * ex['frac_events_mult_rejected']:.2f} %), "
                  f"{sel.n_rows_mult_rej} channel rows dropped")
        if a.mf_strict:
            print(f"[extras] mf-strict: {sel.n_rej_mf_strict} c3b survivors "
                  "rejected (mf_npk != 1)")
        if energy is not None:
            print(f"[extras] energy: {sel.n_no_energy} dropped (non-finite E: no "
                  f"calibration factor), "
                  + ("" if emax_mev is None else
                     f"{sel.n_above_emax} rejected above E_max = {emax_mev:.1f} MeV "
                     f"({a.emax_frac:g} * E_beam {ebeam:.1f}), ")
                  + f"calib {os.path.basename(calib_path)}")
        if a.per_module_best > 0:
            nk = np.bincount(arr["midx"], minlength=len(mods))[pm]
            print(f"[extras] best-K {a.per_module_best}: {arr['midx'].size} stored, "
                  f"{int((nk >= a.per_module_best).sum())} of {int(pm.sum())} modules "
                  f"full, {int((nk == 0).sum())} empty")


# =================================================================== build ==
# Plan section 6.2 with the decisions of the build brief (tail extension, 100-bin
# ROOT references, T_end = 400 ns).  Everything below is pure numpy/scipy except
# `cmd_build` (geometry via _common.setup_pipeline) and the writers.

TS = DT_NS * np.arange(NS_EXPECT)       # sample times in the readout (ns)
SIG_MIN = 0.3                           # ADC floor on the pedestal rms (mf SIGMA_MIN)
WIN = (-20.0, 160.0)                    # alignment / chi2 window, ns after onset
POOL_TAU_MIN = -20.0                    # pooled points: every sample with tau >= this
BIN_W = 0.5                             # ns
T_END = 400.0                           # stored template extent (ns after onset)
T_AXIS_END = 416.0                      # spline domain end on the PRE-shift axis
                                        # (so the shifted spline still covers 400)
T_LEFT = -20.0
TAIL_BIN_MIN_PTS = 5                    # T_data: last bin with >= this many points
TAIL_FIT_FROM = 140.0                   # exponential fitted on [140, T_data]
TAIL_PSEUDO_DT = 4.0
TAIL_TAU_DEFAULT = 170.0
SIGMA_T_NS = 0.15                       # timing-jitter term in the per-bin sigma floor
SMOOTH_LAMBDA = 10000.0                   # penalty lam * int f''(t)^2 dt  (f in peak units); value chosen by the lambda scan rule recorded in the plan progress log
ONSET_FRAC = 0.005                      # onset = last time before peak <= 0.5 % of max
DROP_FRAC = 0.20
N_IT = 3                                # alignment iterations (0, 1, 2)
MAX_IT = 10                             # --max-iterations default (only with --converge-tol)
FREEZE_IT = 3                           # from this iteration on the pool is frozen to iteration FREEZE_IT-1's keep mask
CONV_DF = 1e-3                          # convergence: max |f_k - f_(k-1)| on T_GRID (unit-peak shapes)
CONV_CHI2 = 0.01                        # convergence: chi2/dof not decreasing by more than this fraction
FIT_POINTS = "bins"                     # build --fit-points (set in cmd_build); bins | direct
ITER_COLS = ("dt0_med_ns", "chi2_dof", "max_dshape", "rise_10_90_ns", "width_99_ns", "n_pool")
GRID_DT = 0.1
T_GRID = np.round(-8.0 + GRID_DT * np.arange(int(round((T_END + 8.0) / GRID_DT)) + 1), 6)
LUT_DT = 0.05
LAM_GRID = (0.0, 1e1, 1e2, 1e3, 1e4)    # build --lambda-select heldout
LAM_FALLBACK = 1e2                      # < LAM_MIN_HO held-out pulses
LAM_MIN_HO = 5
LAM_MIN = 0.0                           # build --lambda-min: held-out grid = LAM_GRID values >= this
LAM_CHI2_SIGMA = "noise"                # build --lambda-chi2-sigma (set in cmd_build): noise | full
SHAPE_SIGMA = 0.0                       # build --shape-sigma (set in cmd_build): relative shape term of the direct point sigma
OSC_PROM = 0.002                        # build --osc-escalate: find_peaks prominence / main-peak value
OSC_WIN = (2.0, 150.0)                  # ns: wiggle search window
OSC_EXCL = 3.0                          # ns: |t - t_peak| excluded
OSC_PEAK_WIN = (3.0, 30.0)              # ns: the main peak is the largest value in here
HOLDOUT_FRAC = 0.2
BAND_MIN = 20                           # pulses needed for a lo / hi band
RELAXED_BAND = (100.0, 2500.0)
RELAX_TIER_RULE = {0: "tier>=1 & b0 [300,2500]", 1: "tier>=1 & [100,2500]",
                   2: "run-type average"}
ROOT_NBINS = 100
ROOT_ONSET_BIN = 20
ROOT_DT_BIN = 4.0
ROOT_MODULE_TYPE = 2                    # EventData.h MOD_PbWO4 (hycal.module_type in _raw.root)


def _mix64(x):
    """splitmix64 finaliser on a uint64 array."""
    x = x ^ (x >> np.uint64(30))
    x = x * np.uint64(0xBF58476D1CE4E5B9)
    x = x ^ (x >> np.uint64(27))
    x = x * np.uint64(0x94D049BB133111EB)
    return x ^ (x >> np.uint64(31))


def holdout_mask(run, split, event, midx, seed=1, frac=HOLDOUT_FRAC):
    """True for the candidates held out of `build` (and used by `check`).

    A deterministic hash of (seed, run, split, event, midx), independent of file
    order and of how many candidates exist; `event` is the candidate's
    `event_number` (global physics event number).  Arguments broadcast.
    """
    a = np.broadcast_arrays(*[np.asarray(v) for v in (run, split, event, midx)])
    h = np.full(a[0].shape, int(seed), dtype=np.uint64)
    with np.errstate(over="ignore"):
        h = _mix64(h + np.uint64(0x9E3779B97F4A7C15))
        for v in a:
            h = _mix64(h ^ v.astype(np.int64).astype(np.uint64))
    u = (h >> np.uint64(11)).astype(np.float64) * (1.0 / 9007199254740992.0)
    return u < frac


def holdout_energy(D, n_max, frac):
    """--pool-order energy hold-out (replaces holdout_mask).  Per module, the
    eligible candidates (tier >= 1, band b0, finite mf_energy) are ranked by
    mf_energy descending (stable on the candidate index): ranks [0, n_max) are
    the fit pool (never held out), the next n_hold = int(n_max * frac / (1 -
    frac) + 0.5) ranks are held out for `check`, the rest is not held out.
    Returns (bool mask, n_hold)."""
    n_hold = int(n_max * frac / (1.0 - frac) + 0.5)
    hold = np.zeros(D["midx"].size, dtype=bool)
    E = D["mf_energy"]
    idx = np.nonzero((D["tier"] >= 1) & D["band_b0"] & np.isfinite(E))[0]
    if idx.size == 0:
        return hold, n_hold
    order = idx[np.lexsort((idx, -E[idx], D["midx"][idx]))]
    ms = D["midx"][order]
    start = np.searchsorted(ms, ms, side="left")
    rank = np.arange(order.size) - start
    hold[order[(rank >= n_max) & (rank < n_max + n_hold)]] = True
    return hold, n_hold


# -------------------------------------------------------- shape evaluators --
class TwoTauShape:
    """Unit-peak two-tau (h50) shape; eval(tau) is 0 for tau <= 0."""

    def __init__(self, tr, tf):
        import mf_pulse
        self.tr, self.tf = float(tr), float(tf)
        self.tm = mf_pulse.tmax_of(self.tr, self.tf)
        self._tpl = mf_pulse.template

    def eval(self, tau):
        return self._tpl(tau, self.tr, self.tf, self.tm)


class LutShape:
    """Sampled shape: values v on t = 0, 0.05, ..., T_END; linear interpolation
    (what SampledTemplateBank does), 0 for tau <= 0 or tau >= T_END."""

    def __init__(self, v):
        self.v = np.asarray(v, dtype=np.float64)
        self.n = self.v.size

    @classmethod
    def from_spline(cls, t, c, k=3):
        from scipy.interpolate import BSpline
        tg = LUT_DT * np.arange(int(round(T_END / LUT_DT)) + 1)
        return cls(BSpline(t, c, k)(tg))

    def eval(self, tau):
        tau = np.asarray(tau, dtype=np.float64)
        x = tau * (1.0 / LUT_DT)
        i = np.clip(np.floor(x), 0, self.n - 2).astype(np.int64)
        f = x - i
        v = self.v[i] * (1.0 - f) + self.v[i + 1] * f
        return np.where((tau > 0.0) & (tau < T_END), v, 0.0)


# ------------------------------------------------------------- alignment --
_NCOL = 52                              # sample columns kept around the seed


def _align_columns(t_seed):
    base = np.clip(np.floor((t_seed - 26.0) / DT_NS).astype(np.int64), 0,
                   NS_EXPECT - _NCOL)
    return base[:, None] + np.arange(_NCOL)[None, :]


def _num_den(Ys, cols, t0, shape, win=WIN):
    tau = TS[cols] - t0[:, None]
    msk = (tau >= win[0]) & (tau <= win[1])
    m = shape.eval(tau) * msk
    return (Ys * m).sum(1), (m * m).sum(1)


def align_pulses(Y, sig, t_seed, shape, halfwidth, win=WIN):
    """Fit (A, t0) per pulse against `shape` (anything with .eval(tau)).

    Y[n,100] pedestal-subtracted, sig[n] per-pulse sigma (ADC), t_seed[n] the
    start (ns).  Closed-form A = sum(y m)/sum(m^2) over the samples with
    tau = t - t0 in `win`; t0 searched within +-halfwidth of the seed in three
    nested passes (0.25 ns, 0.05 ns, 0.01 ns steps, each around the previous
    best), i.e. the result sits on a 0.01 ns grid.  Returns dict(t0, A,
    chi2dof (sigma = sig, window `win`, dof = n_window - 2), nwin, ok).
    """
    n = Y.shape[0]
    t_seed = np.asarray(t_seed, dtype=np.float64)
    cols = _align_columns(t_seed)
    Ys = np.take_along_axis(Y, cols, axis=1)
    best_t = t_seed.copy()
    best_s = np.full(n, -np.inf)

    def sweep(center, offs):
        nonlocal best_t, best_s
        bt = center.copy()
        bs = np.full(n, -np.inf)
        for o in offs:
            t0 = center + o
            num, den = _num_den(Ys, cols, t0, shape, win)
            s = np.where(den > 1e-12, num * np.abs(num) / np.maximum(den, 1e-12),
                         -np.inf)
            up = s > bs
            bs = np.where(up, s, bs)
            bt = np.where(up, t0, bt)
        best_t, best_s = bt, bs

    sweep(t_seed, np.arange(-halfwidth, halfwidth + 1e-9, 0.25))
    sweep(best_t, np.arange(-0.25, 0.25 + 1e-9, 0.05))
    sweep(best_t, np.arange(-0.05, 0.05 + 1e-9, 0.01))
    t0 = np.round(best_t, 6)
    return _chi2_at(Ys, cols, sig, t0, shape, win)


def chi2_at(Y, sig, t0, shape, win=WIN):
    """A and chi2/dof at fixed t0 (the final-template quality measure)."""
    cols = _align_columns(np.asarray(t0, dtype=np.float64))
    return _chi2_at(np.take_along_axis(Y, cols, axis=1), cols, sig,
                    np.asarray(t0, dtype=np.float64), shape, win)


def _chi2_at(Ys, cols, sig, t0, shape, win):
    tau = TS[cols] - t0[:, None]
    msk = (tau >= win[0]) & (tau <= win[1])
    m = shape.eval(tau) * msk
    den = (m * m).sum(1)
    A = np.where(den > 1e-12, (Ys * m).sum(1) / np.maximum(den, 1e-12), 0.0)
    r = (Ys - A[:, None] * m) * msk
    nwin = msk.sum(1)
    chi2 = (r * r).sum(1) / (sig * sig)
    dof = np.maximum(nwin - 2, 1)
    ok = (A > 0) & np.isfinite(chi2) & (nwin > 10)
    return dict(t0=t0, A=A, chi2dof=np.where(ok, chi2 / dof, np.inf), nwin=nwin,
                ok=ok)


# ----------------------------------------------------------- pool and bins --
def pool_points(Y, sig, t0, A, sel, ref):
    """Flat (tau, y/A, sigma/A, y/A detrended to the bin centre) of the pulses
    in `sel` for every sample with tau >= POOL_TAU_MIN.

    The within-bin detrend subtracts ref(tau) - ref(bin centre): on the rise a
    0.5 ns bin spans +-6 % of the peak, which would otherwise dominate the
    per-bin MAD.  `ref` is the template the pulses were aligned with."""
    idx = np.flatnonzero(sel)
    tau = TS[None, :] - t0[idx, None]
    keep = tau >= POOL_TAU_MIN
    a = A[idx][:, None]
    y = (Y[idx] / a)[keep]
    sg = np.broadcast_to((sig[idx][:, None] / a), tau.shape)[keep]
    tau = tau[keep]
    b = np.floor((tau - POOL_TAU_MIN) / BIN_W).astype(np.int64)
    ctr = POOL_TAU_MIN + (b + 0.5) * BIN_W
    yd = y - (ref.eval(tau) - ref.eval(ctr))
    return dict(tau=tau, y=y, sg=sg, yd=yd, bin=b)


def _grouped_median(vals, b, nb):
    """Per-bin median of vals (bins b); returns (median[nb], count[nb])."""
    cnt = np.bincount(b, minlength=nb)
    order = np.lexsort((vals, b))
    vs = vals[order]
    start = np.cumsum(cnt) - cnt
    ok = cnt > 0
    lo = np.where(ok, start + (cnt - 1) // 2, 0)
    hi = np.where(ok, start + cnt // 2, 0)
    med = 0.5 * (vs[lo] + vs[hi])
    return np.where(ok, med, np.nan), cnt


def bin_stats(pp, ref):
    """0.5 ns bins: centre, n, median, sigma = max(1.4826 MAD/sqrt(n), floor);
    floor = sqrt(noise^2 + (SIGMA_T_NS * ref'(x))^2)/sqrt(n) with noise the mean
    sigma_pulse/A of the bin (so n = 1 gives sigma_pulse/A plus the timing term).
"""
    b = pp["bin"]
    nb = int(b.max()) + 1
    med, cnt = _grouped_median(pp["yd"], b, nb)
    dev = np.abs(pp["yd"] - np.nan_to_num(med)[b])
    mad, _ = _grouped_median(dev, b, nb)
    ok = cnt > 0
    ctr = POOL_TAU_MIN + (np.arange(nb) + 0.5) * BIN_W
    msg = np.bincount(b, weights=pp["sg"], minlength=nb) / np.maximum(cnt, 1)
    slope = (ref.eval(ctr + 0.05) - ref.eval(ctr - 0.05)) / 0.1
    rn = np.sqrt(np.maximum(cnt, 1))
    floor = np.sqrt(msg ** 2 + (SIGMA_T_NS * slope) ** 2) / rn
    sigma = np.maximum(1.4826 * np.nan_to_num(mad) / rn, floor)
    return dict(x=ctr[ok], n=cnt[ok], y=med[ok], s=sigma[ok])


# ------------------------------------------------------------------ spline --
def knot_plan(n_pool, scale=1):
    """Interior knots (ns, pre-shift axis).  Plan 6.2 step 5 up to 80 ns, then
    8 ns on (80, 160] and 16 ns on (160, 400]; `scale` doubles every spacing
    (Schoenberg-Whitney retries).  One knot at -10 ns holds the pre-onset flat."""
    s0 = 0.5 if n_pool >= 100 else 1.0
    segs = [(-1.0, 10.0, s0), (10.0, 20.0, 1.0), (20.0, 40.0, 2.0),
            (40.0, 80.0, 4.0), (80.0, 160.0, 8.0), (160.0, 400.0, 16.0)]
    k = [-10.0, -1.0]
    for lo, hi, st in segs:
        st *= scale
        k.extend(np.arange(lo + st, hi + 1e-9, st))
    k = np.unique(np.round(np.asarray(k), 6))
    return k[(k > T_LEFT) & (k < T_AXIS_END)]


def tail_extension(bs):
    """T_data, weighted exponential a*exp(-t/tau) on [140, T_data], pseudo-points."""
    from scipy.optimize import curve_fit
    big = bs["n"] >= TAIL_BIN_MIN_PTS
    t_data = float(bs["x"][big].max()) if big.any() else float(bs["x"].max())
    sel = (bs["x"] >= TAIL_FIT_FROM) & (bs["x"] <= t_data)
    a_s = tau_s = np.nan
    if sel.sum() >= 5:
        x, y, s = bs["x"][sel], bs["y"][sel], bs["s"][sel]
        try:
            pos = y > 0
            p0 = (0.02, TAIL_TAU_DEFAULT)
            if pos.sum() >= 3:
                c1 = np.polyfit(x[pos] - TAIL_FIT_FROM, np.log(y[pos]), 1)
                if c1[0] < 0:
                    p0 = (math.exp(c1[1]), -1.0 / c1[0])

            def f(t, a, tau):
                return a * np.exp(-(t - TAIL_FIT_FROM) / tau)

            popt, _ = curve_fit(f, x, y, p0=p0, sigma=s, absolute_sigma=True,
                                bounds=([0.0, 20.0], [1.0, 2000.0]), maxfev=2000)
            tau_s = float(popt[1])
            a_s = float(popt[0]) * math.exp(TAIL_FIT_FROM / tau_s)   # a at t = 0
        except Exception:                                           # noqa: BLE001
            a_s = tau_s = np.nan
    if not np.isfinite(tau_s):
        tail = bs["x"] >= t_data - 20.0
        yv = float(np.median(bs["y"][tail]))
        tau_use = TAIL_TAU_DEFAULT
        a_use = max(yv, 0.0) * math.exp(t_data / tau_use)
    else:
        tau_use, a_use = tau_s, a_s
    tx = np.arange(t_data + TAIL_PSEUDO_DT, T_AXIS_END + 1e-9, TAIL_PSEUDO_DT)
    wmed = float(np.median(1.0 / bs["s"]))
    return dict(T_data=t_data, tau_slow=tau_s, a_slow=a_s, tau_used=tau_use,
                a_used=a_use, tx=tx, ty=a_use * np.exp(-tx / tau_use),
                w_pseudo=wmed)


_PEN_CACHE = {}


def _penalty(t):
    """Gram matrix G = int B_i''(t) B_j''(t) dt of the cubic basis on knots t
    (0.05 ns Riemann sum), cached per knot vector."""
    from scipy.interpolate import BSpline
    key = t.tobytes()
    if key not in _PEN_CACHE:
        n = t.size - 4
        g = np.arange(T_LEFT, T_AXIS_END, 0.05)
        d2 = BSpline(t, np.eye(n), 3).derivative(2)(g)
        _PEN_CACHE[key] = (d2.T @ d2) * 0.05
        if len(_PEN_CACHE) > 8:
            _PEN_CACHE.pop(next(iter(_PEN_CACHE)))
    return _PEN_CACHE[key]


def fit_spline(bs, n_pool, pp=None, ref=None):
    """Weighted cubic B-spline through the bin medians + the tail pseudo-points
    (weights 1/sigma_bin), penalised by SMOOTH_LAMBDA * int f''^2.

    The trigger-synchronous readout puts most pulses of a module at ONE
    sub-sample phase (t0 mod 4 ns; Rayleigh R ~ 0.8-0.9), so the 0.5 ns bins
    are filled only near that phase and an unpenalised least-squares spline
    (make_lsq_spline) is free to oscillate by tens of percent (or is singular)
    between the populated phases; the tiny penalty only acts where the data
    leave the spline undetermined.  If the solve fails every knot spacing is
    doubled (the plan's Schoenberg-Whitney retry).  Returns dict(t, c, scale,
    tail, n_knots) on the pre-shift axis, or None.

    FIT_POINTS == "direct" (build --fit-points direct; needs pp, ref): for tau <=
    T_data the pooled points themselves are fitted (x = tau, y = y/A, not the
    bin-detrended value) with sigma = sqrt((sigma_pulse/A)^2 + (SIGMA_T_NS *
    ref'(tau))^2) instead of the 0.5 ns bin medians (SHAPE_SIGMA > 0 adds
    (SHAPE_SIGMA * ref(tau))^2 under the root); tail extension, pseudo-points
    and their weight are the bins' (tail_extension(bs))."""
    from scipy.interpolate import BSpline
    tail = tail_extension(bs)
    if FIT_POINTS == "direct" and pp is not None:
        dm = pp["tau"] <= tail["T_data"]
        xd, yd_ = pp["tau"][dm], pp["y"][dm]
        slope = (ref.eval(xd + 0.05) - ref.eval(xd - 0.05)) / 0.1
        var = pp["sg"][dm] ** 2 + (SIGMA_T_NS * slope) ** 2
        if SHAPE_SIGMA > 0:
            var = var + (SHAPE_SIGMA * ref.eval(xd)) ** 2
        sgd = np.sqrt(var)
        x = np.concatenate([xd, tail["tx"]])
        y = np.concatenate([yd_, tail["ty"]])
        w = np.concatenate([1.0 / sgd, np.full(tail["tx"].size, tail["w_pseudo"])])
    else:
        d = bs["x"] <= tail["T_data"]
        x = np.concatenate([bs["x"][d], tail["tx"]])
        y = np.concatenate([bs["y"][d], tail["ty"]])
        w = np.concatenate([1.0 / bs["s"][d], np.full(tail["tx"].size, tail["w_pseudo"])])
    o = np.argsort(x, kind="stable")
    x, y, w = x[o], y[o], w[o]
    for scale in (1, 2, 4):
        ik = knot_plan(n_pool, scale)
        t = np.concatenate([[T_LEFT] * 4, ik, [T_AXIS_END] * 4])
        n = t.size - 4
        try:
            Bx = BSpline.design_matrix(x, t, 3).toarray()
            A = Bx.T @ (Bx * (w * w)[:, None]) + SMOOTH_LAMBDA * _penalty(t)
            c = np.linalg.solve(A, Bx.T @ (w * w * y))
        except (ValueError, np.linalg.LinAlgError):
            continue
        if np.all(np.isfinite(c)):
            return dict(t=t, c=c, scale=scale, tail=tail, n_knots=int(ik.size))
    return None


def finalize_spline(t, c, shift=True):
    """Unit peak and onset at 0.  Returns dict(t, c) on the shifted axis plus
    t_on (the shift, 0 if shift=False) and the pre-normalisation peak."""
    from scipy.interpolate import BSpline
    spl = BSpline(t, c, 3)
    g = np.arange(-15.0, 80.0, 0.01)
    v = spl(g)
    ip = int(np.argmax(v))
    vmax = float(v[ip])
    t_on = 0.0
    if shift:
        below = np.flatnonzero(v[:ip] <= ONSET_FRAC * vmax)
        if below.size:
            j = int(below[-1])
            v0, v1 = v[j], v[j + 1]
            fr = (ONSET_FRAC * vmax - v0) / (v1 - v0) if v1 > v0 else 0.0
            t_on = float(g[j] + fr * 0.01)
        else:
            t_on = float(g[0])
    return dict(t=t - t_on, c=c / vmax, t_on=t_on, vmax=vmax)


def spline_values(t, c, tt):
    from scipy.interpolate import BSpline
    return BSpline(t, c, 3)(tt)


def derived_from_spline(t, c, ctx, with_iint=True):
    """t_pk, t50, iint_wa (phase-averaged walk integral), area -- on the shifted
    axis; the same recipe as SampledTemplateBank (shape 0 outside (0, 400)).
    `build` passes with_iint=False (iint_wa = NaN) and fills the integral from
    bank_integrals() instead; G1-9 keeps the direct-spline value."""
    import mf_pulse
    from scipy.interpolate import BSpline
    spl = BSpline(t, c, 3)

    def shp(tau):
        return np.where((tau > 0) & (tau < T_END), spl(np.clip(tau, t[3], t[-4])), 0.0)

    g = np.arange(0.0, 60.0, 0.01)
    v = shp(g)
    ip = int(np.argmax(v))
    t_pk = float(g[ip])
    i1 = int(np.argmax(v[:ip + 1] >= 0.5))
    i0 = max(i1 - 1, 0)
    fr = (0.5 - v[i0]) / (v[i1] - v[i0]) if v[i1] > v[i0] else 0.0
    t50 = float(g[i0] + fr * 0.01)
    area = float(shp(DT_NS * np.arange(int(T_END / DT_NS) + 1)).sum())
    if not with_iint:
        return dict(t_pk=t_pk, t50=t50, iint_wa=float("nan"), area=area)
    u = DT_NS * np.arange(mf_pulse.NIPH) / mf_pulse.NIPH
    iint = mf_pulse.unit_integral_fn(shp, u, ctx["smooth_order"],
                                     ctx["int_tail_ratio"], ctx["tail_break_n"])
    return dict(t_pk=t_pk, t50=t50, iint_wa=float(iint.mean()), area=area)


def bank_integrals(splines, ctx, return_bank=False):
    """Energy integral exactly as the bank defines it: the shapes (list of
    (t, c) cubic B-splines) go into one in-memory SampledTemplateBank (same LUT,
    same walk-out code as mf_pulse.build_bank on the written .npz) and the result
    is I_cont (= I_mean, the mean of the exact discrete walk-out over the 256
    sub-sample phases) and the peak-to-peak of the 256-phase table over I_cont."""
    import mf_pulse
    from scipy.interpolate import BSpline
    rows = [mf_pulse._ClippedSpline(BSpline(np.asarray(t, dtype=np.float64),
                                            np.asarray(c, dtype=np.float64), 3))
            for t, c in splines]
    n = len(rows)
    bank = mf_pulse.SampledTemplateBank(
        rows, [f"s{i}" for i in range(n)], [True] * n,
        smooth_order=ctx["smooth_order"], int_tail_ratio=ctx["int_tail_ratio"],
        tail_break_n=ctx["tail_break_n"], T_end=T_END)
    tab = bank.I_tab[:, :mf_pulse.NIPH].astype(np.float64)
    out = (np.asarray(bank.I_cont, dtype=np.float64).copy(),
           (tab.max(axis=1) - tab.min(axis=1)) / bank.I_cont)
    return out + (bank,) if return_bank else out


# ------------------------------------------------------------- one shape --
class PulseSet:
    """Pedestal-subtracted pulses of one pool (all arrays length n)."""

    def __init__(self, samples, ped, rms, t_seed):
        self.Y = samples.astype(np.float64) - ped[:, None]
        self.sig = np.maximum(rms, SIG_MIN)
        self.seed = np.asarray(t_seed, dtype=np.float64)
        self.n = self.Y.shape[0]


def _rise_width(f):
    """(10-90 % rise time, width of the region >= 0.99 of the peak), ns, of a
    finalised shape f on T_GRID; linear interpolation walking outward from the
    peak (NaN if a level is not crossed inside the grid)."""
    ip = int(np.argmax(f))
    pk = float(f[ip])

    def cross(level, step):
        j = ip
        while 0 <= j + step < f.size and f[j] >= level:
            j += step
        if f[j] >= level:
            return np.nan
        k = j - step                      # the point nearer the peak (>= level)
        return float(T_GRID[j] + (level - f[j]) / (f[k] - f[j]) * (T_GRID[k] - T_GRID[j]))

    rise = cross(0.9 * pk, -1) - cross(0.1 * pk, -1)
    width = cross(0.99 * pk, 1) - cross(0.99 * pk, -1)
    return rise, width


def build_shape(ps, ref0, *, seed_hw=3.0, later_hw=2.0, nit=N_IT,
                drop_first=False, shift=True, rng=None, conv_tol=None,
                max_it=MAX_IT, diag=False):
    """Alignment / pooling / spline iterations (plan 6.2 steps 2-6) for one set.

    Iteration 0 aligns against `ref0` (anything with .eval; the h50 two-tau for a
    module); every later iteration against the previous spline.  From iteration 1
    on (or 0 with drop_first) the worst DROP_FRAC of ALL pulses by chi2/dof
    (sigma = rms, no model floor) are left out of the pool; from iteration
    FREEZE_IT (3) on the pool is frozen to iteration 2's keep mask (and the
    current alignment's ok), later iterations only realign it.  Returns None if
    the spline fit fails, else a dict with the final template and diagnostics.

    nit iterations are run; with conv_tol (ns) the loop then continues while not
    converged and fewer than max_it iterations are done.  Converged (evaluated
    for every iteration k >= FREEZE_IT against k - 1): median |t0 - t_seed| over
    the pool < conv_tol, max |f_k - f_(k-1)| on T_GRID < CONV_DF and the median
    chi2/dof not decreasing by more than CONV_CHI2 (fraction).  diag: also return
    hist (one row per iteration: ITER_COLS) and converged.
    """
    n = ps.n
    ref = ref0
    t_seed = ps.seed
    hw = seed_hw
    N, shifts, chi2_it, hist = [], [], [], []
    prev = None
    keep_frozen = None
    f_prev = None
    converged = False
    it = 0
    while True:
        t_in = t_seed
        al = align_pulses(ps.Y, ps.sig, t_seed, ref, hw)
        ok = al["ok"]
        chi2_it.append(float(np.median(al["chi2dof"][ok])) if ok.any() else np.nan)
        if keep_frozen is not None:
            keep = keep_frozen & ok
        elif it >= 1 or drop_first:
            nkeep = max(int(round((1.0 - DROP_FRAC) * n)), 1)
            order = np.argsort(np.where(ok, al["chi2dof"], np.inf), kind="stable")
            keep = np.zeros(n, dtype=bool)
            keep[order[:nkeep]] = True
            keep &= ok
        else:
            keep = ok.copy()
        if it == FREEZE_IT - 1:
            keep_frozen = keep.copy()
        if keep.sum() < 5:
            return None
        pp = pool_points(ps.Y, ps.sig, al["t0"], al["A"], keep, ref)
        bs = bin_stats(pp, ref)
        sp = fit_spline(bs, int(keep.sum()), pp, ref)
        if sp is None:
            return None
        fin = finalize_spline(sp["t"], sp["c"], shift=shift)
        N.append(int(keep.sum()))
        shifts.append(fin["t_on"])
        if diag:
            f = spline_values(fin["t"], fin["c"], T_GRID)
            dt0 = (float(np.median(np.abs(al["t0"] - t_in)[keep])) if it >= 1
                   else np.nan)
            chg = float(np.max(np.abs(f - f_prev))) if it >= 1 else np.nan
            hist.append((dt0, chi2_it[-1], chg, *_rise_width(f), int(keep.sum())))
            f_prev = f
            if it >= FREEZE_IT and conv_tol is not None:
                c1, c2 = chi2_it[-2], chi2_it[-1]
                converged = bool(dt0 < conv_tol and chg < CONV_DF
                                 and (c1 - c2) / c1 <= CONV_CHI2)
        prev = dict(al=al, keep=keep, pp=pp, bs=bs, sp=sp, fin=fin)
        ref = LutShape.from_spline(fin["t"], fin["c"])
        t_seed = al["t0"] + fin["t_on"]       # t0 on the new (shifted) axis
        hw = later_hw
        it += 1
        if it >= nit and (conv_tol is None or converged or it >= max_it):
            break
    # final (t0, A, chi2/dof) of ALL pulses against the final template
    q = align_pulses(ps.Y, ps.sig, t_seed, ref, 0.5)
    prev.update(N=N, shifts=shifts, chi2_it=chi2_it, final=q, ref=ref)
    if diag:
        prev.update(hist=np.array(hist, dtype=np.float64), converged=converged)
    return prev


def band_shape(ps, ref):
    """Amplitude-band template (plan 6.2 step 7, steps 2-6 repeated): iteration 0
    aligns the band pool against the module's final template `ref` (no drop);
    iteration 1 aligns against the band's own iteration-0 spline and drops the
    worst 20 % against it.  No onset shift: the band stays on the module's axis."""
    return build_shape(ps, ref, seed_hw=3.0, nit=2, drop_first=False, shift=False)


# ------------------------------------------------------------ candidates --
BUILD_FIELDS = ("event_number", "midx", "samples", "ped", "rms", "amp", "tier",
                "band_b0", "mf_t0", "p1_shape_ok")


def _sha1(path):
    h = hashlib.sha1()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def load_candidates(files, fields=BUILD_FIELDS):
    """Concatenate the candidate files (order of `files`).  Returns a dict with
    the requested fields plus `run`, `split` (per candidate), `cutflow` (summed),
    `cutflow_names`, `peak_hist` (summed), `metas`, `inputs`."""
    out = {k: [] for k in fields}
    runs, splits, metas, inputs = [], [], [], []
    cutflow = peak_hist = names0 = None
    for f in files:
        z = np.load(f, allow_pickle=False)
        meta = json.loads(str(z["meta"]))
        if meta.get("schema") != SCHEMA:
            sys.exit(f"[ERROR] {f}: schema {meta.get('schema')} != {SCHEMA}")
        n = int(z["midx"].size)
        for k in fields:
            if k not in z.files:
                sys.exit(f"[ERROR] {f}: no `{k}` array (--pool-order energy needs "
                         "candidate files from `select` with --per-module-best or "
                         "--emax-frac, which store mf_energy)")
            out[k].append(z[k])
        runs.append(np.full(n, int(z["run"]), dtype=np.int64))
        splits.append(np.full(n, int(z["split"]), dtype=np.int64))
        names = [str(x) for x in z["cutflow_names"]]
        if names0 is None:
            names0, cutflow = names, z["cutflow"].astype(np.int64)
            peak_hist = z["peak_hist"].astype(np.int64)
        elif names != names0:
            sys.exit(f"[ERROR] {f}: cutflow columns differ from the first file")
        else:
            cutflow += z["cutflow"]
            peak_hist += z["peak_hist"]
        st = os.stat(f)
        inputs.append(dict(path=os.path.abspath(f), size=st.st_size,
                           mtime=datetime.fromtimestamp(st.st_mtime).isoformat(),
                           sha1=_sha1(f), run=int(z["run"]), split=int(z["split"]),
                           n_candidates=n))
        metas.append(meta)
    D = {k: np.concatenate(v) for k, v in out.items()}
    D.update(run=np.concatenate(runs), split=np.concatenate(splits),
             cutflow=cutflow, cutflow_names=names0, peak_hist=peak_hist,
             metas=metas, inputs=inputs)
    return D


def _wave_ctx(metas):
    """The wave config (smooth_order, int_tail_ratio, tail_break_n) of the
    candidate files; they must agree."""
    w0 = {k: metas[0]["wave_config"][k]
          for k in ("smooth_order", "int_tail_ratio", "tail_break_n")}
    for m in metas[1:]:
        w = {k: m["wave_config"][k] for k in w0}
        if w != w0:
            sys.exit(f"[ERROR] candidate files disagree on the wave config: "
                     f"{w0} vs {w}")
    return w0


# ---------------------------------------------------------- per module ----
def _row_arrays(fin, bs, sp):
    """shape / shape_err / resid on T_GRID (float32, NaN where undefined)."""
    t, c, t_on, vmax = fin["t"], fin["c"], fin["t_on"], fin["vmax"]
    shape = spline_values(t, c, T_GRID)
    x = bs["x"] - t_on
    d = bs["x"] <= sp["tail"]["T_data"]
    xd = x[d]
    err = np.interp(T_GRID, xd, bs["s"][d] / vmax, left=np.nan, right=np.nan)
    rs = (bs["y"][d] / vmax - spline_values(t, c, xd))
    res = np.interp(T_GRID, xd, rs, left=np.nan, right=np.nan)
    return (shape.astype(np.float32), err.astype(np.float32),
            res.astype(np.float32), xd, rs, bs["s"][d] / vmax, bs["n"][d],
            bs["y"][d] / vmax)


def _pad(rows, n):
    out = np.full((len(rows), n), np.nan)
    for i, r in enumerate(rows):
        out[i, :len(r)] = r
    return out


def build_lambda_select(ctx, ps, ref0, idx):
    """--lambda-select heldout: build_shape for each LAM_GRID value >= LAM_MIN on the
    same pool; the held-out pulses (ctx['hold']) of the module are aligned against each
    final template (seed = mf_t0 + cumulative onset shift, +-3 ns) and the lambda with
    the lowest mean chi2/dof over the pulses ok for every lambda is kept (chi2 with
    sigma = rms; LAM_CHI2_SIGMA == "full": _chi2_full at the alignment's A, t0).
    Fewer than LAM_MIN_HO held-out pulses: LAM_FALLBACK.  Leaves SMOOTH_LAMBDA = the
    chosen value.  Returns (r, lam, flag, n_ho, ho_chi2[len(LAM_GRID)], NaN for the
    values not tried); r is None if the build failed."""
    global SMOOTH_LAMBDA
    D = ctx["D"]
    grid = tuple(x for x in LAM_GRID if x >= LAM_MIN)
    hc = np.full(len(LAM_GRID), np.nan)
    h = idx[ctx["hold"][idx]]
    if h.size < LAM_MIN_HO:
        SMOOTH_LAMBDA = LAM_FALLBACK
        return build_shape(ps, ref0, **ctx["iter_kw"]), LAM_FALLBACK, 1, 0, hc
    hs = PulseSet(D["samples"][h], D["ped"][h], D["rms"][h], D["mf_t0"][h])
    rs, als = [], []
    for lam in grid:
        SMOOTH_LAMBDA = lam
        r = build_shape(ps, ref0, **ctx["iter_kw"])
        rs.append(r)
        als.append(None if r is None else align_pulses(
            hs.Y, hs.sig, hs.seed + float(np.sum(r["shifts"])), r["ref"], 3.0))
    good = [i for i, al in enumerate(als) if al is not None]
    if not good:
        return None, np.nan, 1, 0, hc
    ok = np.all([als[i]["ok"] for i in good], axis=0)
    n_ho = int(ok.sum())
    if n_ho == 0:                       # nothing comparable: fixed fallback
        k = grid.index(LAM_FALLBACK) if LAM_FALLBACK in grid \
            and rs[grid.index(LAM_FALLBACK)] is not None else good[0]
        flag = 1
    else:
        hg = np.full(len(grid), np.nan)
        for i in good:
            c2 = als[i]["chi2dof"] if LAM_CHI2_SIGMA == "noise" else _chi2_full(
                hs.Y, hs.sig, als[i]["t0"], als[i]["A"], rs[i]["ref"])
            hg[i] = float(np.mean(c2[ok]))
        k = int(np.argmin(np.where(np.isnan(hg), np.inf, hg)))
        for i, x in enumerate(grid):
            hc[LAM_GRID.index(x)] = hg[i]
        flag = 0
    SMOOTH_LAMBDA = grid[k]
    return rs[k], grid[k], flag, n_ho, hc


def _phase_strat(E, t0, n_max, nb=8):
    """--pool-order energy-phase: indices of n_max candidates taken in descending E
    with at most ceil(n_max/nb) per bin of the sub-sample phase t0 mod DT_NS (nb equal
    bins); if the bins cannot supply n_max, the rest is filled by descending E.
    Returned sorted."""
    o = np.argsort(-E, kind="stable")
    b = np.minimum((np.mod(t0[o], DT_NS) / (DT_NS / nb)).astype(np.int64), nb - 1)
    cap = -(-n_max // nb)
    rank = np.zeros(o.size, dtype=np.int64)          # running count within the bin
    for k in range(nb):
        m = b == k
        rank[m] = np.arange(int(m.sum()))
    take = o[rank < cap][:n_max]
    if take.size < n_max:
        take = np.concatenate([take, o[rank >= cap][:n_max - take.size]])
    return np.sort(take)


def _osc_count(f):
    """--osc-escalate: number of spurious extrema (maxima + minima, prominence >=
    OSC_PROM * peak, in OSC_WIN, farther than OSC_EXCL from the main peak) of the finalised
    shape f on T_GRID; the main peak is the largest f in OSC_PEAK_WIN."""
    from scipy.signal import find_peaks
    f = np.asarray(f, dtype=np.float64)
    mp = (T_GRID >= OSC_PEAK_WIN[0]) & (T_GRID <= OSC_PEAK_WIN[1])
    ip = int(np.flatnonzero(mp)[np.argmax(f[mp])])
    pk, t_pk = f[ip], T_GRID[ip]
    w = np.flatnonzero((T_GRID >= OSC_WIN[0]) & (T_GRID <= OSC_WIN[1]))
    n = 0
    for sg in (f[w], -f[w]):
        j = w[find_peaks(sg, prominence=OSC_PROM * pk)[0]]
        n += int((np.abs(T_GRID[j] - t_pk) > OSC_EXCL).sum())
    return n


def process_module(ctx, midx, idx):
    """One PbWO4 module: the D4 chain, three iterations, amplitude bands.
    idx = candidate indices (into ctx['D']) of this module."""
    global SMOOTH_LAMBDA
    D, a = ctx["D"], ctx["args"]
    heldout = a.lambda_select == "heldout"
    amp, tier = D["amp"][idx], D["tier"][idx]
    held = ctx["hold"][idx]
    p0 = (~held) & (tier >= 1) & D["band_b0"][idx]
    p1 = (~held) & (tier >= 1) & (amp >= RELAXED_BAND[0]) & (amp <= RELAXED_BAND[1])
    chain = []
    if p0.sum() >= a.n_min:
        chain.append((0, p0, BAND_B0))
    if p1.sum() >= a.n_min_relaxed:
        chain.append((1, p1, RELAXED_BAND))
    rng = np.random.default_rng([int(a.seed), int(midx)])
    name = ctx["names"][midx]
    ref0 = TwoTauShape(*ctx["tau"].get(name, ctx["tau_agg"]))
    res = dict(midx=midx, fallback=2, valid=1, n_cand=int(p1.sum()),
               band=(np.nan, np.nan), n_pool=0, N=(0, 0, 0), kept_gidx=None,
               amp_b0=amp[p0].astype(np.float64),     # H of the f(A) r_mod waveforms
               **(dict(lam=np.nan, lam_flag=2, n_ho=0, ho_chi2=np.full(len(LAM_GRID), np.nan))
                  if heldout else {}),
               **(dict(osc_n_before=-1, lam_final=np.nan, osc_escalated=0, osc_fallback=0)
                  if a.osc_escalate else {}))
    for fb, mask, band in chain:
        pool = idx[mask]
        n_cand = int(pool.size)
        if pool.size > a.n_max:
            rnd = np.sort(rng.choice(pool, a.n_max, replace=False))
            if a.pool_order == "energy":        # rng stream kept; top n_max by E
                o = np.argsort(-D["mf_energy"][pool], kind="stable")
                pool = np.sort(pool[o[:a.n_max]])
            elif a.pool_order == "energy-phase":   # top by E, <= ceil(n_max/8) per phase bin
                pool = pool[_phase_strat(D["mf_energy"][pool], D["mf_t0"][pool], a.n_max)]
            else:
                pool = rnd
        ps = PulseSet(D["samples"][pool], D["ped"][pool], D["rms"][pool],
                      D["mf_t0"][pool])
        if heldout:
            r, lam, lam_flag, n_ho, ho_chi2 = build_lambda_select(ctx, ps, ref0, idx)
        else:
            r = build_shape(ps, ref0, **ctx["iter_kw"])
        if r is None:
            continue
        if a.osc_escalate:                  # wiggles: lambda x10 and rebuild the same pool
            n0 = n = _osc_count(spline_values(r["fin"]["t"], r["fin"]["c"], T_GRID))
            steps = 0
            while n > 0 and SMOOTH_LAMBDA < 1e4:
                SMOOTH_LAMBDA = min(1e1 if SMOOTH_LAMBDA <= 0 else SMOOTH_LAMBDA * 10, 1e4)
                r2 = build_shape(ps, ref0, **ctx["iter_kw"])
                steps += 1
                if r2 is None:
                    break
                r = r2
                n = _osc_count(spline_values(r["fin"]["t"], r["fin"]["c"], T_GRID))
            res.update(osc_n_before=n0, osc_escalated=steps)
            if n > 0:                       # still wiggly: run-type average (fallback 2)
                res["osc_fallback"] = 1
                break
            res["lam_final"] = SMOOTH_LAMBDA
        fin, bs, sp = r["fin"], r["bs"], r["sp"]
        if len(r["N"]) > 3:                 # fixed-shape arrays: iterations 0, 1, last
            r["N"] = [r["N"][i] for i in (0, 1, -1)]
            r["chi2_it"] = [r["chi2_it"][i] for i in (0, 1, -1)]
        shape, err, resid, bx, brs, bsig, bn, by = _row_arrays(fin, bs, sp)
        keep = r["keep"] & r["final"]["ok"]
        tail = sp["tail"]
        npl = min(3000, r["pp"]["tau"].size)
        sel = np.sort(rng.choice(r["pp"]["tau"].size, npl, replace=False))
        t0f = r["final"]["t0"]
        ph = np.exp(2j * np.pi * t0f[keep] / DT_NS)
        d = derived_from_spline(fin["t"], fin["c"], ctx["wave"], with_iint=False)
        res.update(
            n_p1ok=int(D["p1_shape_ok"][pool].sum()),
            fallback=fb, valid=1, band=band, n_cand=n_cand, n_pool=int(pool.size),
            **(dict(lam=lam, lam_flag=lam_flag, n_ho=n_ho, ho_chi2=ho_chi2)
               if heldout else {}),
            N=tuple(r["N"]), chi2=float(np.median(r["final"]["chi2dof"][keep])),
            shape=shape, err=err, resid=resid, t=fin["t"], c=fin["c"],
            shift=float(r["shifts"][-1]), shift0=float(r["shifts"][0]),
            vmax=fin["vmax"], T_data=tail["T_data"], tau_slow=tail["tau_slow"],
            a_slow=tail["a_slow"], tau_used=tail["tau_used"], a_used=tail["a_used"],
            knot_scale=sp["scale"], n_knots=sp["n_knots"], **d,
            rise_empty=int(30 - np.unique(np.floor(bs["x"][(bs["x"] > 0) & (bs["x"] < 15)]
                                                   / BIN_W)).size),
            phase_R=float(abs(ph.mean())) if ph.size else np.nan,
            pts=(np.float32(r["pp"]["tau"][sel] - fin["t_on"]),
                 np.float32(r["pp"]["y"][sel] / fin["vmax"])),
            bins=(bx, by, bsig, bn, brs),
            kept_gidx=pool[keep], kept_t0=t0f[keep], ref=r["ref"],
            cum_shift=float(np.sum(r["shifts"])), chi2_it=tuple(r["chi2_it"]),
            **({} if "hist" not in r else dict(
                iter_hist=r["hist"], n_iterations=r["hist"].shape[0],
                converged=(-1 if ctx["iter_kw"]["conv_tol"] is None
                           else int(r["converged"])))))
        # amplitude bands against the module's own final template
        cand = ~held & (tier >= 1)
        # lo / hi (G1-9) first, then the f(A) bands (default mid / top, or --bands):
        # the lo / hi draws are unchanged
        for key, lo, hi, incl in (("lo", *BAND_LO, False), ("hi", *BAND_HI, True),
                                  *[(b["key"], b["lo"], b["hi"], b["incl"])
                                    for b in ctx["fa_bands"] if b["key"] != "lo"]):
            m = cand & (amp >= lo) & ((amp <= hi) if incl else (amp < hi))
            res["n_" + key] = int(m.sum())
            if m.sum() < BAND_MIN:
                continue
            bp = idx[m]
            if bp.size > a.n_max:
                bp = np.sort(rng.choice(bp, a.n_max, replace=False))
            bps = PulseSet(D["samples"][bp], D["ped"][bp], D["rms"][bp],
                           D["mf_t0"][bp] + res["cum_shift"])
            br = band_shape(bps, r["ref"])
            if br is None:
                continue
            res["shape_" + key] = spline_values(br["fin"]["t"], br["fin"]["c"],
                                                T_GRID).astype(np.float32)
            res["spl_" + key] = (br["fin"]["t"], br["fin"]["c"])
            res["amp_" + key] = D["amp"][bp]
        break
    if heldout or a.osc_escalate:
        SMOOTH_LAMBDA = float(a.smooth_lambda)      # the run-type average uses --smooth-lambda
    return res


def run_average(ctx, results):
    """Run-type average: <= 20 aligned pulses from each fallback-0 module (from
    every fallback <= 1 module if fewer than 50 of fallback 0), the module
    iterations with the median module shape as the starting reference."""
    D = ctx["D"]
    rng = np.random.default_rng([int(ctx["args"].seed), 987654321])
    use = [r for r in results if r["fallback"] == 0]
    rule = "fallback 0"
    if len(use) < 50:
        use = [r for r in results if r["fallback"] <= 1 and r["kept_gidx"] is not None]
        rule = "fallback <= 1 (< 50 of fallback 0)"
    if not use:
        return None, rule, 0
    g, t0 = [], []
    for r in use:
        n = r["kept_gidx"].size
        sel = np.sort(rng.choice(n, min(20, n), replace=False))
        g.append(r["kept_gidx"][sel])
        t0.append(r["kept_t0"][sel])
    g, t0 = np.concatenate(g), np.concatenate(t0)
    shp = np.nanmedian(np.stack([r["shape"] for r in use]).astype(np.float64), axis=0)
    tl = LUT_DT * np.arange(int(round(T_END / LUT_DT)) + 1)
    ref0 = LutShape(np.interp(tl, T_GRID, np.nan_to_num(shp)))
    ps = PulseSet(D["samples"][g], D["ped"][g], D["rms"][g], t0)
    r = build_shape(ps, ref0, seed_hw=1.0)
    return r, rule, int(g.size)


# --------------------------------------------------------------- writers ----
def _root_hist(shape_fn):
    tau = (np.arange(ROOT_NBINS) - ROOT_ONSET_BIN) * ROOT_DT_BIN
    return np.where(tau > 0, shape_fn(tau), 0.0)


def write_root(path, ctx, rows, avg_fn, meta_json):
    import uproot
    t0 = time.monotonic()
    nr = len(rows)
    cols = {k: np.zeros(nr, dtype=np.int32) for k in (
        "module_id", "module_type", "nseen", "nused", "polarity", "valid",
        "nsamples", "peak_sample", "onset_bin", "fallback", "n_iter")}
    cols["best_event"] = np.full(nr, -1, dtype=np.int64)
    for k in ("band_lo", "band_hi", "chi2_ndf", "t_pk_ns", "t50_ns", "iint_wa",
              "dt_bin_ns", "ref_rms_pre"):
        cols[k] = np.zeros(nr)
    edges = np.arange(ROOT_NBINS + 1, dtype=np.float64)
    with uproot.recreate(path) as f:
        for i, r in enumerate(rows):
            fn = r["shape_fn"]
            h = _root_hist(fn)
            f[f"refs/ref_{r['module_id']}"] = (h, edges)
            cols["module_id"][i] = r["module_id"]
            cols["module_type"][i] = ROOT_MODULE_TYPE
            cols["nseen"][i] = r["n_cand"]
            cols["nused"][i] = r["N"][2]
            cols["polarity"][i] = 1
            cols["valid"][i] = r["valid"]
            cols["nsamples"][i] = ROOT_NBINS
            cols["peak_sample"][i] = ROOT_ONSET_BIN + int(round(r["t_pk"] / ROOT_DT_BIN))
            cols["onset_bin"][i] = ROOT_ONSET_BIN
            cols["fallback"][i] = r["fallback"]
            cols["n_iter"][i] = r.get("n_iterations", 3) if r["fallback"] < 2 else 0
            cols["band_lo"][i], cols["band_hi"][i] = r["band"]
            cols["chi2_ndf"][i] = r.get("chi2", np.nan)     # median of the kept 80 % of pulses
            cols["t_pk_ns"][i] = r["t_pk"]
            cols["t50_ns"][i] = r["t50"]
            cols["iint_wa"][i] = r["iint_wa"]
            cols["dt_bin_ns"][i] = ROOT_DT_BIN
            cols["ref_rms_pre"][i] = r.get("rms_pre", np.nan)
        f["refs/ref_avg"] = (_root_hist(avg_fn), edges)
        # mktree + extend: a dict assignment would write an RNTuple in uproot 5.7
        f.mktree("ref_summary", {k: v.dtype for k, v in cols.items()})
        f["ref_summary"].extend(cols)
        f["meta"] = meta_json
    print(f"[root] {path}  {nr} references  ({time.monotonic() - t0:.1f} s)",
          flush=True)


def write_tsv(path, rows, cutflow_names, fa_keys=None, lam=False, osc=False):
    """fa_keys (--bands): the n_mid / n_top / r_lo / r_mid / r_top columns are
    replaced by n_<band> / r_<band> columns of the generalised f(A) bands."""
    leg = fa_keys is None
    cols = ["name", "module_id", "module_index", "fallback", "valid", "tier_rule",
            "n_cand", "n_pool", "N0", "N1", "N2", "band_lo", "band_hi", "chi2_ndf",
            "t_pk", "t50", "iint_wa", "area", "onset_shift_ns", "onset_shift_it0",
            "T_data", "tau_slow", "a_slow", "knot_scale", "n_knots",
            "rise_empty_bins", "phase_R", "n_lo", "n_hi"] \
        + (["n_mid", "n_top"] if leg else []) \
        + ["lohi_maxdiff", "lohi_iint_ratio", "iint_phase_pp"] \
        + (["r_lo", "r_mid", "r_top"] if leg else []) + ["r_mod"] \
        + ([] if leg else ["n_" + k for k in fa_keys] + ["r_" + k for k in fa_keys]) \
        + ["n_p1ok"] + ["cf_" + c for c in cutflow_names[:N_SEQ]] \
        + (["lam", "lam_flag", "n_ho"] + ["ho_chi2_l" + ("0" if x == 0 else f"1e{round(math.log10(x))}") for x in LAM_GRID]
           if lam else []) \
        + (["osc_n_before", "lam_final", "osc_escalated", "osc_fallback"] if osc else [])

    def fmt(v):
        if isinstance(v, (float, np.floating)):
            return "nan" if not np.isfinite(v) else f"{v:.6g}"
        return str(v)

    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            v = [r["name"], r["module_id"], r["midx"], r["fallback"], r["valid"],
                 RELAX_TIER_RULE[r["fallback"]], r["n_cand"], r["n_pool"],
                 *r["N"], *r["band"], r.get("chi2", np.nan), r["t_pk"], r["t50"],
                 r["iint_wa"], r["area"], r["cum_shift"], r["shift0"], r["T_data"],
                 r["tau_slow"], r["a_slow"], r.get("knot_scale", 0),
                 r.get("n_knots", 0), r.get("rise_empty", -1),
                 r.get("phase_R", np.nan), r.get("n_lo", 0), r.get("n_hi", 0),
                 *([r.get("n_mid", 0), r.get("n_top", 0)] if leg else []),
                 r.get("lohi_maxdiff", np.nan), r.get("lohi_iint_ratio", np.nan),
                 r["iint_pp"],
                 *([r.get("r_lo", np.nan), r.get("r_mid", np.nan), r.get("r_top", np.nan)]
                   if leg else []),
                 r.get("r_mod", np.nan),
                 *([] if leg else [r.get("n_" + k, 0) for k in fa_keys]
                   + [r.get("r_" + k, np.nan) for k in fa_keys]),
                 r.get("n_p1ok", -1),
                 *[int(x) for x in r["cutflow"][:N_SEQ]],
                 *([r["lam"], r["lam_flag"], r["n_ho"], *r["ho_chi2"]] if lam else []),
                 *([r["osc_n_before"], r["lam_final"], r["osc_escalated"], r["osc_fallback"]]
                   if osc else [])]
            fh.write("\t".join(fmt(x) for x in v) + "\n")


# ------------------------------------------------------------------- PDF ----
class _Pages:
    """One reused figure: per module a main/tail/residual/band page."""

    def __init__(self):
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec
        self.plt = plt
        self.fig = plt.figure(figsize=(11, 8.5), dpi=90)
        gs = GridSpec(3, 2, figure=self.fig, height_ratios=[3.2, 2.2, 1.5],
                      left=0.07, right=0.98, top=0.93, bottom=0.03,
                      hspace=0.34, wspace=0.18)
        self.a_main = self.fig.add_subplot(gs[0, 0])
        self.a_tail = self.fig.add_subplot(gs[0, 1])
        self.a_res = self.fig.add_subplot(gs[1, 0])
        sub = gs[1, 1].subgridspec(2, 1, height_ratios=[2, 1], hspace=0.08)
        self.a_band = self.fig.add_subplot(sub[0])
        self.a_diff = self.fig.add_subplot(sub[1], sharex=self.a_band)
        self.a_txt = self.fig.add_subplot(gs[2, :])
        self.axes = [self.a_main, self.a_tail, self.a_res, self.a_band,
                     self.a_diff, self.a_txt]

    def module(self, r, ctx):
        from scipy.interpolate import BSpline
        fig = self.fig
        for ax in self.axes:
            ax.cla()
        am, at, ar, ab, ad, ax_t = self.axes
        name = r["name"]
        tt = np.arange(-10.0, 80.0, 0.1)
        spl = BSpline(r["t"], r["c"], 3)
        h50 = TwoTauShape(*ctx["tau"].get(name, ctx["tau_agg"]))
        fb = r["fallback"]
        fig.suptitle(f"{name}  (module id {r['module_id']})   fallback {fb}: "
                     f"{(SINGLE_TIER_RULE if ctx.get('single') else RELAX_TIER_RULE)[fb]}"
                     f"   run type {ctx['run_type']}",
                     fontsize=11)
        if r.get("pts") is not None:
            px, py = r["pts"]
            am.scatter(px, py, s=2, c="0.55", lw=0, rasterized=True,
                       label=f"pooled points ({px.size})")
            bx, by, bsg, bn, brs = r["bins"]
            e = np.interp(tt, bx, bsg, left=np.nan, right=np.nan)
            am.fill_between(tt, spl(tt) - e, spl(tt) + e, color="C3", alpha=0.3,
                            lw=0, rasterized=True, label="+- bin sigma")
            am.plot(bx, by, ".", ms=3, c="C0", label="bin medians", rasterized=True)
        am.plot(tt, spl(tt), "C3", lw=1.2, label="spline")
        am.plot(tt, h50.eval(tt), "k--", lw=0.9, label="h50 two-tau")
        am.set_xlim(-10, 60)
        am.set_ylim(-0.05, 1.1)
        am.set_xlabel("t - onset (ns)")
        am.set_ylabel("y / A")
        am.legend(fontsize=7, loc="upper right")
        am.grid(alpha=0.3)
        # tail, log y
        t2 = np.arange(0.0, 300.0, 0.5)
        if r.get("bins") is not None:
            bx, by, bsg, bn, brs = r["bins"]
            p = by > 0
            at.errorbar(bx[p], by[p], bsg[p], fmt=".", ms=2, c="C0", lw=0.5,
                        rasterized=True)
            te = t2[t2 >= r["T_data"] - r["shift"]]
            at.plot(te, r["a_used"] * np.exp(-(te + r["shift"]) / r["tau_used"])
                    / r["vmax"], "C2--", lw=1.2,
                    label=f"exp ext. tau={r['tau_used']:.0f} ns")
            at.axvline(r["T_data"] - r["shift"], c="0.5", lw=0.7, ls=":")
        at.plot(t2, spl(t2), "C3", lw=1.0)
        at.plot(t2, h50.eval(t2), "k--", lw=0.8)
        at.set_yscale("log")
        at.set_ylim(2e-3, 1.5)
        at.set_xlim(0, 300)
        at.set_xlabel("t - onset (ns)")
        at.set_title("tail (log), dotted = T_data", fontsize=8)
        if at.get_legend_handles_labels()[0]:
            at.legend(fontsize=7, loc="upper right")
        at.grid(alpha=0.3, which="both")
        # residuals
        if r.get("bins") is not None:
            bx, by, bsg, bn, brs = r["bins"]
            ar.errorbar(bx, 100 * brs, 100 * bsg, fmt=".", ms=2, c="C0", lw=0.4,
                        rasterized=True)
            lim = float(np.clip(np.nanpercentile(np.abs(100 * brs), 99) * 1.3, 0.3, 6))
            ar.set_ylim(-lim, lim)
        else:
            ar.text(0.5, 0.5, "no own pooled data (run-type average)", ha="center",
                    va="center", transform=ar.transAxes, fontsize=9, color="0.4")
        ar.axhline(0, c="k", lw=0.6)
        ar.set_xlim(-20, 260)
        ar.set_xlabel("t - onset (ns)")
        ar.set_ylabel("bin median - spline (% of peak)")
        ar.grid(alpha=0.3)
        # amplitude bands
        t3 = np.arange(0.0, 40.0, 0.1)
        ab.plot(t3, spl(t3), "k", lw=1.0, label="module")
        lo, hi = r.get("shape_lo"), r.get("shape_hi")
        if lo is not None:
            ab.plot(t3, np.interp(t3, T_GRID, lo), "C0", lw=0.9,
                    label=f"b_lo 100-400 (n={r['n_lo']})")
        if hi is not None:
            ab.plot(t3, np.interp(t3, T_GRID, hi), "C3", lw=0.9,
                    label=f"b_hi 400-2500 (n={r['n_hi']})")
        if lo is not None and hi is not None:
            ad.plot(t3, 100 * (np.interp(t3, T_GRID, hi) - np.interp(t3, T_GRID, lo)),
                    "C1", lw=0.9)
        ab.legend(fontsize=7, loc="lower right")
        ab.set_ylabel("shape")
        ab.tick_params(labelbottom=False)
        ad.axhline(0, c="k", lw=0.5)
        ad.set_ylabel("hi - lo (%)", fontsize=8)
        ad.set_xlabel("t - onset (ns)")
        ab.grid(alpha=0.3)
        ad.grid(alpha=0.3)
        # text
        ax_t.axis("off")
        cf = r["cutflow"]
        cfn = ctx["cutflow_names"]
        items = [f"{n} {int(v)}" for n, v in zip(cfn[:N_SEQ], cf[:N_SEQ])]
        cft = "  ".join(items[:6]) + "\n          " + "  ".join(items[6:])
        lines = [
            f"N per iteration (0,1,2): {r['N'][0]} / {r['N'][1]} / {r['N'][2]}   "
            f"pool {r['n_pool']} of {r['n_cand']} candidates   band "
            f"[{r['band'][0]:.0f}, {r['band'][1]:.0f}]   median chi2/dof (kept 80 %) "
            f"{r.get('chi2', float('nan')):.1f}",
            f"t_pk {r['t_pk']:.2f} ns   t50 {r['t50']:.2f} ns   iint_wa "
            f"{r['iint_wa']:.4f} (phase p-p {100 * r['iint_pp']:.2f} %)   area {r['area']:.3f}   onset shift total "
            f"{r['cum_shift']:.2f} (last it {r['shift']:.2f}, it0 {r['shift0']:.2f}) ns   knots x{r.get('knot_scale', 0)}"
            f" ({r.get('n_knots', 0)})",
            f"T_data {r['T_data']:.0f} ns   tail tau_slow {r['tau_slow']:.0f} ns   "
            f"a_slow {r['a_slow']:.3f}   empty rise bins {r.get('rise_empty', -1)}   "
            f"phase R {r.get('phase_R', float('nan')):.2f}   lo/hi max|d| "
            f"{r.get('lohi_maxdiff', float('nan')) * 100:.2f} %   I_int hi/lo "
            f"{r.get('lohi_iint_ratio', float('nan')):.4f}",
            "cut flow: " + cft,
        ]
        if fb == 2:
            lines.insert(0, "FALLBACK 2: this module has no own template; the "
                            "run-type average is shown")
        ax_t.text(0.0, 1.0, "\n".join(lines), va="top", ha="left", fontsize=7.5,
                  family="monospace", transform=ax_t.transAxes)
        return fig


def cover_pages(pdf, ctx, rows, avg, avg_rule, n_avg, args):
    import matplotlib.pyplot as plt
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Rectangle
    fb = np.array([r["fallback"] for r in rows])
    ncand = np.array([r["n_cand"] for r in rows], dtype=float)
    fig = plt.figure(figsize=(11, 8.5))
    ax = fig.add_axes([0.06, 0.08, 0.60, 0.84])
    g = ctx["geom"]
    ax.add_collection(PatchCollection(
        [Rectangle((x - sx / 2, y - sy / 2), sx, sy)
         for x, y, sx, sy in g["glass"]], facecolor="none", edgecolor="0.8", lw=0.3))
    rect = [Rectangle((x - sx / 2, y - sy / 2), sx, sy) for x, y, sx, sy in g["rows"]]
    val = np.log10(np.maximum(ncand, 1))
    own = fb < 2
    pc = PatchCollection([p for p, o in zip(rect, own) if o], cmap="viridis",
                         edgecolor="none")
    pc.set_array(val[own])
    pc.set_clim(0, max(np.log10(max(ncand.max(), 10)), 1))
    ax.add_collection(pc)
    ax.add_collection(PatchCollection(
        [Rectangle((x - sx / 2, y - sy / 2), sx, sy) for (x, y, sx, sy), o
         in zip(g["rows"], own) if not o],
        facecolor="0.8", edgecolor="0.45", hatch="xxxx", lw=0))
    sel = [p for p, f in zip(g["rows"], fb) if f == 1]
    if sel:
        ax.add_collection(PatchCollection(
            [Rectangle((x - sx / 2, y - sy / 2), sx, sy) for x, y, sx, sy in sel],
            facecolor="none", edgecolor="w", hatch="////", lw=0))
    ax.set_xlim(g["xlim"])
    ax.set_ylim(g["ylim"])
    ax.set_aspect("equal")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    cb = fig.colorbar(pc, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("log10(candidates in the pool)")
    nfb = [int((fb == k).sum()) for k in range(3)]
    ax.set_title(f"{ctx['run_type']}: coverage; white hatch = fallback 1, grey "
                 f"cross-hatch = fallback 2 (average)", fontsize=9)
    txt = [f"Run type {ctx['run_type']}",
           f"runs {sorted(set(int(x) for x in ctx['D']['run']))}",
           f"{len(ctx['files'])} candidate files, {ctx['D']['midx'].size} candidates",
           ("single-pulse reference per module (Wassim eq. 18): no hold-out"
            if ctx.get("single") else
            f"held out for check: {int(ctx['hold'].sum())} "
            f"({100 * ctx['hold'].mean():.1f} %, hash seed {args.seed})"
            if args.pool_order == "random" else
            f"held out for check: {int(ctx['hold'].sum())} "
            f"({100 * ctx['hold'].mean():.1f} %; energy rule, next "
            f"{ctx['n_hold']} per module after the top "
            f"{args.n_max})"),
           "", f"modules (PbWO4 rows): {len(rows)}",
           (f"  fallback 0 (single reference, tier 2): {nfb[0]}" if ctx.get("single") else
            f"  fallback 0 (tier>=1, b0, >= {args.n_min}): {nfb[0]}"),
           (f"  fallback 1 (single reference, tier 1): {nfb[1]}" if ctx.get("single") else
            f"  fallback 1 (tier>=1, [100,2500], >= {args.n_min_relaxed}): {nfb[1]}"),
           f"  fallback 2 (run-type average): {nfb[2]}", "",
           f"pool N (candidates), fb 0: median "
           f"{np.median(ncand[fb == 0]) if nfb[0] else float('nan'):.0f}",
           f"average template: {n_avg} pulses from",
           f"  {avg_rule}",
           f"  T_data {avg['sp']['tail']['T_data']:.0f} ns, tau_slow "
           f"{avg['sp']['tail']['tau_slow']:.0f} ns"]
    pk = ctx.get("p1ok")
    if pk and pk["n_pool"]:
        txt += [f"PASS 1 shape ok (pooled): {pk['n_ok']}/{pk['n_pool']}"
                f" = {100 * pk['n_ok'] / pk['n_pool']:.1f} %"]
    ac = ctx.get("ampcorr")
    if ac and ac.get("f") is not None:
        txt += ["", f"f(A) ({ac['n_modules']} modules, PCHIP through the band knots):"]
        txt += [f"  {b:3s} A {k['A']:.0f}  f {k['f']:.4f} +- {k['se']:.4f} "
                f"(n {k['n_modules']})" for b, k in ac["knots"].items()]
    chi = np.array([r.get("chi2", np.nan) for r in rows])
    ok = fb < 2
    txt += ["", f"median chi2/dof (final, kept 80 %): {np.nanmedian(chi[ok]):.1f}"]
    fig.text(0.70, 0.92, "\n".join(txt), va="top", ha="left", fontsize=8,
             family="monospace")
    pdf.savefig(fig)
    plt.close(fig)

    # page 2: cut flow, peak histogram, distributions
    D = ctx["D"]
    fig = plt.figure(figsize=(11, 8.5))
    gs = fig.add_gridspec(2, 2, left=0.07, right=0.97, top=0.93, bottom=0.07,
                          hspace=0.3, wspace=0.25)
    at = fig.add_subplot(gs[0, 0])
    at.axis("off")
    cf = ctx["cutflow_tot"]
    nms = ctx["cutflow_names"]
    lines = ["Cut flow, PbWO4 channels summed over all files", ""]
    prev = None
    for n, v in zip(nms, cf):
        lines.append(f"{n:34s} {int(v):>12d}"
                     + (f"  {100.0 * v / max(prev, 1):6.1f} %" if (prev and n in nms[1:N_SEQ]) else ""))
        if n in nms[:N_SEQ]:
            prev = v
    at.text(0, 1, "\n".join(lines), va="top", family="monospace", fontsize=7)
    ap = fig.add_subplot(gs[0, 1])
    ap.bar(np.arange(NS_EXPECT) * DT_NS, D["peak_hist"], width=DT_NS * 0.9)
    ap.axvspan(P1_T_LO, P1_T_HI, color="C1", alpha=0.2, label="PASS 1 window 150-170")
    ap.axvline(PEAK_MAX_SAMPLE * DT_NS, c="r", ls="--", lw=0.8, label="containment cut")
    ap.set_xlabel("peak sample time (ns)")
    ap.set_ylabel("channels")
    ap.set_title("raw argmax (clip + A>=50 channels)", fontsize=8)
    ap.legend(fontsize=7)
    a1 = fig.add_subplot(gs[1, 0])
    a1.hist(np.log10(np.maximum(ncand, 1)), bins=30, color="C0")
    a1.set_xlabel("log10(candidates in the pool)")
    a1.set_ylabel("modules")
    a2 = fig.add_subplot(gs[1, 1])
    tau = np.array([r["tau_slow"] for r in rows if r["fallback"] < 2])
    a2.hist(tau[np.isfinite(tau)], bins=40, color="C2")
    a2.set_xlabel("tail tau_slow (ns), own templates")
    a2.set_ylabel("modules")
    pdf.savefig(fig)
    plt.close(fig)


def write_pdf(path, ctx, rows, avg, avg_rule, n_avg, args):
    from matplotlib.backends.backend_pdf import PdfPages
    t0 = time.monotonic()
    pg = _Pages()
    with PdfPages(path) as pdf:
        cover_pages(pdf, ctx, rows, avg, avg_rule, n_avg, args)
        for i, r in enumerate(rows):
            fig = pg.module(r, ctx)
            pdf.savefig(fig, dpi=90)
            if i % 200 == 0:
                print(f"  [pdf] {i}/{len(rows)}  {time.monotonic() - t0:.0f} s",
                      flush=True)
    print(f"[pdf] {path}  ({time.monotonic() - t0:.0f} s)", flush=True)


class _PngShim:
    """Stands in for PdfPages: saves each figure as cover_<n>.png."""

    def __init__(self, d):
        self.d, self.n = Path(d), 0

    def savefig(self, fig, **kw):
        self.n += 1
        fig.savefig(self.d / f"cover_{self.n}.png", dpi=90)


def write_pngs(png_dir, ctx, rows, want, cover_args=None):
    pg = _Pages()
    Path(png_dir).mkdir(parents=True, exist_ok=True)
    if cover_args is not None:
        cover_pages(_PngShim(png_dir), ctx, rows, *cover_args)
    for r in rows:
        if r["name"] in want:
            pg.module(r, ctx).savefig(Path(png_dir) / f"{r['name']}.png", dpi=110)


# ---------------------------------------------------------------- command ----
def _iter_kw(a):
    """build_shape kwargs of the module templates: {} for the default
    (3 iterations, no convergence test), else the --iterations options."""
    if a.iterations == N_IT and a.converge_tol is None:
        return {}
    return dict(nit=a.iterations, conv_tol=a.converge_tol, max_it=a.max_iterations,
                diag=True)


def _print_iter_log(ctx, results):
    """Per-iteration diagnostics of the module templates (non-default
    --iterations / --converge-tol only)."""
    a, kw = ctx["args"], ctx["iter_kw"]
    res = [r for r in results if "iter_hist" in r]
    if not res:
        return
    if a.modules:
        for r in res:
            for k, h in enumerate(r["iter_hist"]):
                print(f"  [iter] {ctx['names'][r['midx']]} it {k}: dt0_med {h[0]:.4f} ns  "
                      f"chi2/dof {h[1]:.3f}  max|df| {h[2]:.2e}  rise10-90 {h[3]:.3f} ns  "
                      f"width>=0.99 {h[4]:.3f} ns  N {int(h[5])}"
                      + (f"  converged {r['converged']}" if k == r["iter_hist"].shape[0] - 1
                         else ""), flush=True)
    m = max(r["iter_hist"].shape[0] for r in res)
    print(f"[iter] module templates: --iterations {a.iterations}"
          + ("" if a.converge_tol is None else
             f" --converge-tol {a.converge_tol:g} ns --max-iterations {a.max_iterations}")
          + f"; {len(res)} modules; medians over the modules that reached each iteration",
          flush=True)
    print("[iter]  it  n_mod  " + "  ".join(f"{c:>13s}" for c in ITER_COLS), flush=True)
    for k in range(m):
        H = np.array([r["iter_hist"][k] for r in res if r["iter_hist"].shape[0] > k])
        with np.errstate(all="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)       # all-NaN column
            med = np.nanmedian(H, axis=0) if np.isfinite(H).any() else np.full(6, np.nan)
        print(f"[iter] {k:3d}  {H.shape[0]:5d}  " + "  ".join(f"{v:13.5g}" for v in med),
              flush=True)
    if a.converge_tol is not None:
        cv = np.array([r["converged"] for r in res])
        ni = np.array([r["iter_hist"].shape[0] for r in res])
        print(f"[iter] converged {int((cv == 1).sum())}/{len(res)}; stopped at "
              f"--max-iterations {a.max_iterations} without converging "
              f"{int(((cv == 0) & (ni >= a.max_iterations)).sum())}; iterations run "
              f"min/median/max {ni.min()}/{np.median(ni):.0f}/{ni.max()}", flush=True)


def cmd_build(a):
    global SMOOTH_LAMBDA, FIT_POINTS, CONV_DF, SHAPE_SIGMA, LAM_CHI2_SIGMA, LAM_MIN, LAM_FALLBACK
    if a.single_pulse:
        return cmd_build_single(a)
    import _common as C
    import mf_pulse
    SMOOTH_LAMBDA = float(a.smooth_lambda)       # read by fit_spline and recorded in the meta
    FIT_POINTS = a.fit_points                    # read by fit_spline
    CONV_DF = float(a.converge_df)               # read by build_shape, recorded in the meta
    SHAPE_SIGMA = float(a.shape_sigma)           # read by fit_spline (direct points)
    LAM_CHI2_SIGMA = a.lambda_chi2_sigma         # read by build_lambda_select
    LAM_MIN = float(a.lambda_min)                # read by build_lambda_select
    LAM_FALLBACK = float(a.lambda_fallback)      # read by build_lambda_select, recorded in the meta
    if a.iterations < FREEZE_IT or a.max_iterations < FREEZE_IT:
        sys.exit(f"[ERROR] --iterations and --max-iterations must be >= {FREEZE_IT}")
    if a.converge_tol is not None and not a.converge_tol > 0:
        sys.exit("[ERROR] --converge-tol must be > 0")
    if not a.converge_df > 0:
        sys.exit("[ERROR] --converge-df must be > 0")
    if SHAPE_SIGMA < 0:
        sys.exit("[ERROR] --shape-sigma must be >= 0")
    if not any(x >= LAM_MIN for x in LAM_GRID):
        sys.exit(f"[ERROR] --lambda-min {LAM_MIN:g} leaves no value of {list(LAM_GRID)}")
    t_start = time.monotonic()
    files = sorted({f for pat in a.pulses.split(",") for f in glob.glob(pat.strip())})
    if not files:
        sys.exit(f"[ERROR] no files match {a.pulses}")
    stem = Path(a.out)
    stem.parent.mkdir(parents=True, exist_ok=True)
    energy_pool = a.pool_order in ("energy", "energy-phase")
    D = load_candidates(files, BUILD_FIELDS + ("mf_energy",) if energy_pool
                        else BUILD_FIELDS)
    wave = _wave_ctx(D["metas"])
    runs = sorted(set(int(x) for x in D["run"]))
    print(f"[load] {len(files)} files, runs {runs}, {D['midx'].size} candidates "
          f"({time.monotonic() - t_start:.1f} s)", flush=True)
    if energy_pool:
        hold, n_hold = holdout_energy(D, a.n_max, a.holdout_frac)
        print(f"[pool] --pool-order {a.pool_order}: fit pool = top {a.n_max} by mf_energy "
              + ("(<= ceil(n_max/8) per 0.5 ns phase bin of mf_t0) " if a.pool_order == "energy-phase" else "")
              + 
              f"per module; held out for check = the next {n_hold} per module "
              f"(tier>=1, band b0, finite mf_energy): {int(hold.sum())} candidates",
              flush=True)
    else:
        n_hold = None
        hold = holdout_mask(D["run"], D["split"], D["event_number"], D["midx"],
                            seed=a.seed, frac=a.holdout_frac)

    p = C.setup_pipeline(evio_path="", run_num=runs[0])
    mods = [p.hycal.module(i) for i in range(p.hycal.module_count())]
    names = [m.name for m in mods]
    pwo4 = [bool(m.is_pwo4()) for m in mods]
    rows_idx = [i for i, f in enumerate(pwo4) if f]
    if a.modules:
        want = set(s.strip() for s in a.modules.split(","))
        rows_idx = [i for i in rows_idx if names[i] in want]
    h50 = Path(a.h50)
    if not h50.is_absolute():
        h50 = HERE / h50
    column = mf_pulse.resolve_column(str(h50), "auto", prefer="25320")
    tau_rows, h50_sha1, agg_json = mf_pulse.load_h50(str(h50), column)
    tau_agg = mf_pulse.load_aggregate(agg_json)
    print(f"[setup] h50 {h50.name} column {column}: {len(tau_rows)} modules with "
          f"own (tau_r, tau_f), aggregate {tau_agg}", flush=True)

    glass = [(m.x, m.y, m.size_x, m.size_y) for m in mods
             if m.is_hycal() and not m.is_pwo4()]
    geo_rows = [(mods[i].x, mods[i].y, mods[i].size_x, mods[i].size_y)
                for i in rows_idx]
    allg = glass + geo_rows
    xs = [x - sx / 2 for x, y, sx, sy in allg] + [x + sx / 2 for x, y, sx, sy in allg]
    ys = [y - sy / 2 for x, y, sx, sy in allg] + [y + sy / 2 for x, y, sx, sy in allg]
    pm = np.array(pwo4)
    cutflow_tot = D["cutflow"][pm].sum(axis=0)
    ctx = dict(D=D, args=a, iter_kw=_iter_kw(a), hold=hold, n_hold=n_hold, names=names, tau=tau_rows, tau_agg=tau_agg,
               band_edges=a.bands,
               fa_bands=FA_BANDS_DEFAULT if a.bands is None else fa_bands_from_edges(a.bands),
               wave=wave, run_type=a.run_type, files=files,
               cutflow_names=D["cutflow_names"], cutflow_tot=cutflow_tot,
               geom=dict(glass=glass, rows=geo_rows, xlim=(min(xs), max(xs)),
                         ylim=(min(ys), max(ys))))

    order = np.argsort(D["midx"], kind="stable")
    ms = D["midx"][order]
    lo = np.searchsorted(ms, np.arange(len(mods)), side="left")
    hi = np.searchsorted(ms, np.arange(len(mods)), side="right")
    results = []
    t_mod = time.monotonic()
    for k, mi in enumerate(rows_idx):
        r = process_module(ctx, mi, order[lo[mi]:hi[mi]])
        results.append(r)
        if (k + 1) % 100 == 0 or k + 1 == len(rows_idx):
            nf = np.bincount([x["fallback"] for x in results], minlength=3)
            print(f"  [modules] {k + 1}/{len(rows_idx)}  fallback {nf.tolist()}  "
                  f"{time.monotonic() - t_mod:.0f} s", flush=True)

    if ctx["iter_kw"]:
        _print_iter_log(ctx, results)
    if a.lambda_select == "heldout":
        lam = np.array([r["lam"] for r in results])
        fl = np.bincount([r["lam_flag"] for r in results], minlength=3)
        print("[lambda] heldout: modules per lambda " + "  ".join(
            f"{x:g}:{int((lam == x).sum())}" for x in sorted(set(LAM_GRID) | {LAM_FALLBACK}))
            + f";  flag 0/1/2 = {fl[0]}/{fl[1]}/{fl[2]}", flush=True)
    if a.osc_escalate:
        ev = [r for r in results if r["osc_n_before"] > 0]
        st = np.bincount([r["osc_escalated"] for r in ev], minlength=5)
        print(f"[osc] escalate: modules flagged before {len(ev)}; escalated by 1/2/3 steps "
              f"{st[1]}/{st[2]}/{st[3]}" + (f" (4 steps {st[4]})" if st[4:].sum() else "")
              + f"; osc_fallback {sum(r['osc_fallback'] for r in results)}", flush=True)
    t_avg = time.monotonic()
    av, avg_rule, n_avg = run_average(ctx, results)
    if av is None:
        sys.exit("[ERROR] no module with an own template: cannot build the "
                 "run-type average")
    avg_shape, avg_err, avg_res = _row_arrays(av["fin"], av["bs"], av["sp"])[:3]
    avg_fin = av["fin"]
    avg_der = derived_from_spline(avg_fin["t"], avg_fin["c"], wave, with_iint=False)
    avg_tail = av["sp"]["tail"]
    print(f"[average] {n_avg} pulses from {avg_rule}: N {av['N']} "
          f"chi2/dof med {np.median(av['final']['chi2dof']):.1f}  T_data "
          f"{avg_tail['T_data']:.0f}  tau_slow {avg_tail['tau_slow']:.0f}  "
          f"({time.monotonic() - t_avg:.0f} s)", flush=True)
    avg_t, avg_c = avg_fin["t"], avg_fin["c"]
    ctx["avg_shape_fn"] = (lambda tt: np.where(
        (tt > 0) & (tt < T_END),
        spline_values(avg_t, avg_c, np.clip(tt, avg_t[3], avg_t[-4])), 0.0))

    # fallback-2 rows and the per-row derived fields
    rows = []
    for r in results:
        mi = r["midx"]
        r["name"], r["module_id"] = names[mi], int(mods[mi].id)
        r["cutflow"] = D["cutflow"][mi]
        if r["fallback"] == 2:
            r.update(shape=avg_shape, err=avg_err, resid=avg_res,
                     t=avg_fin["t"], c=avg_fin["c"], shift=float(av["shifts"][-1]),
                     shift0=float(av["shifts"][0]),
                     cum_shift=float(np.sum(av["shifts"])), vmax=avg_fin["vmax"],
                     T_data=avg_tail["T_data"], tau_slow=avg_tail["tau_slow"],
                     a_slow=avg_tail["a_slow"], tau_used=avg_tail["tau_used"],
                     a_used=avg_tail["a_used"], pts=None, bins=None, **avg_der)
            r["N"] = (0, 0, 0)
        d = r.get("shape_lo"), r.get("shape_hi")
        if d[0] is not None and d[1] is not None:
            m = (T_GRID >= 0) & (T_GRID <= 20)
            r["lohi_maxdiff"] = float(np.nanmax(np.abs(d[1][m].astype(float)
                                                       - d[0][m])))
        # pre-onset noise of the pooled points (ROOT ref_rms_pre)
        if r.get("pts") is not None:
            px, py = r["pts"]
            q = px < -2.0
            r["rms_pre"] = float(np.sqrt(np.mean(py[q] ** 2))) if q.any() else np.nan
        sp_t, sp_c = r["t"], r["c"]
        r["shape_fn"] = (lambda tt, t=sp_t, c=sp_c:
                         np.where((tt > 0) & (tt < T_END),
                                  spline_values(t, c, np.clip(tt, t[3], t[-4])), 0.0))
        rows.append(r)
    ctx.update(p=p, rows=rows, av=av, avg_rule=avg_rule, n_avg=n_avg, avg_fin=avg_fin,
               avg_shape=avg_shape, avg_der=avg_der, avg_tail=avg_tail, h50=h50,
               h50_sha1=h50_sha1, column=column, agg_json=agg_json, mods=mods,
               pwo4=pwo4, runs=runs, t_start=t_start, stem=stem)
    return finish_build(ctx)


# ------------------------------------------- build --single-pulse (Wassim ref) --
# One best valid PASS 1 candidate per module (over all input files, the macro's
# betterCandidate order) -> its reference R(k) = (W(k) - bg_mean) / max_k(W - bg_mean)
# (Wassim eq. 18, 100 samples at 4 ns) as a not-a-knot cubic INTERPOLATING spline (the
# TSpline3 analogue; no smoothing, no pooling, no alignment iterations), in the bank
# format.  The exponential tail extension exists only because the bank needs T_END.
SINGLE_TIER_RULE = {0: "single best valid_p1 candidate, tier 2 (strict bg)",
                    1: "single best valid_p1 candidate, tier 1 (loose bg)",
                    2: "run-type average"}
SINGLE_FIELDS = ("event_number", "midx", "samples", "ped", "rms", "amp", "tier",
                 "peak_sample", "first_thr", "bg_mean", "bg_rms", "bg_maxabs",
                 "roughness")
SINGLE_AVG_STEP = 5                     # run-type average: T_GRID every 0.5 ns


def _best_wassim(D):
    """{module index: row of D} of the best candidate over ALL rows (all files in
    order), the macro's betterCandidate order: tier desc, amp desc, roughness asc,
    bg_rms asc; exact ties go to the first row (strict comparisons)."""
    n = int(D["midx"].size)
    if n == 0:
        return {}
    o = np.lexsort((np.arange(n), D["bg_rms"], D["roughness"], -D["amp"],
                    -D["tier"].astype(np.int64), D["midx"]))
    ms = D["midx"][o]
    first = np.concatenate([[True], ms[1:] != ms[:-1]])
    return {int(D["midx"][i]): int(i) for i in o[first]}


def _interp_spline(x, y):
    """Natural cubic interpolating spline (knots / coefficients of the BSpline)."""
    from scipy.interpolate import make_interp_spline
    s = make_interp_spline(np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64),
                           k=3, bc_type="not-a-knot")
    return s.t, s.c


def wassim_template(samples, ped, bg_mean, bg_rms, ref0):
    """Bank-format template of one reference pulse.  Returns a dict with the final
    spline (fin: t, c; unit spline peak, onset at 0), the onset / normalisation
    bookkeeping, the tail extension and `shift` = onset on the template axis minus the
    onset of a two-tau (`ref0`) fit to the pulse (what `check` adds to the h50-axis
    mf_t0)."""
    W = samples.astype(np.float64) - ped
    R = W - bg_mean
    R = R / float(R.max())                              # eq. 18
    ip = int(np.argmax(R))
    below = np.flatnonzero(R[:ip] <= ONSET_FRAC)
    t_est = float(TS[0])
    if below.size:
        j = int(below[-1])
        fr = (ONSET_FRAC - R[j]) / (R[j + 1] - R[j]) if R[j + 1] > R[j] else 0.0
        t_est = float(TS[j] + DT_NS * fr)
    x0 = TS - t_est                                     # onset near 0 (finalize_spline window)
    f0 = finalize_spline(*_interp_spline(x0, R), shift=True)    # onset of the data spline
    x1 = x0 - f0["t_on"]
    y1 = R / f0["vmax"]                                 # unit spline peak
    tail = tail_extension(dict(x=x1, y=y1, s=np.ones(x1.size),
                               n=np.full(x1.size, TAIL_BIN_MIN_PTS)))
    t, c = _interp_spline(np.concatenate([x1, tail["tx"]]),
                          np.concatenate([y1, tail["ty"]]))
    fin = finalize_spline(t, c, shift=True)
    t_on_ro = t_est + f0["t_on"] + fin["t_on"]          # onset (readout time) on the final axis
    al = align_pulses((W - bg_mean)[None, :], np.array([max(bg_rms, SIG_MIN)]),
                      np.array([t_on_ro]), ref0, 16.0)
    shift = float(t_on_ro - al["t0"][0]) if al["ok"][0] else float("nan")
    return dict(fin=fin, tail=tail, x=x1 - fin["t_on"], y=y1 / fin["vmax"],
                t_on0=float(f0["t_on"]), t_on1=float(fin["t_on"]),
                peak_over_sample=float(f0["vmax"] * fin["vmax"]), shift=shift,
                t_on_ro=float(t_on_ro))


def finalize_integrals_single(ctx):
    """iint_wa / iint_phase_pp of every row from one in-memory SampledTemplateBank
    (bank_integrals); no f(A), no bands."""
    rows = ctx["rows"]
    t0 = time.monotonic()
    shapes, own = [], {}
    for r in rows:
        if r["fallback"] < 2:
            own[r["midx"]] = len(shapes)
            shapes.append((r["t"], r["c"]))
    avg_i = len(shapes)
    shapes.append((ctx["avg_fin"]["t"], ctx["avg_fin"]["c"]))
    I, pp = bank_integrals(shapes, ctx["wave"])
    for r in rows:
        k = own.get(r["midx"], avg_i)
        r["iint_wa"], r["iint_pp"] = float(I[k]), float(pp[k])
    ctx["avg_der"]["iint_wa"] = float(I[avg_i])
    ctx["ampcorr"] = dict(rule="none: single-pulse reference bank (f = 1)", scope="none",
                          run_type=ctx["run_type"], n_modules=0, knots={}, f=None,
                          logA=None, per_run_type=None)
    n_own = len(own)
    ctx["p1ok"] = dict(n_pool=n_own, n_ok=n_own)
    print(f"[integrals] {len(shapes)} shapes through one in-memory bank "
          f"({time.monotonic() - t0:.0f} s); phase p-p median "
          f"{100 * np.median([r['iint_pp'] for r in rows]):.2f} %, max "
          f"{100 * max(r['iint_pp'] for r in rows):.2f} %", flush=True)


def write_tsv_single(path, rows, cutflow_names):
    cols = ["name", "module_id", "module_index", "fallback", "valid", "tier_rule",
            "n_cand", "run", "split", "event_number", "amp", "tier", "peak_sample",
            "first_thr", "roughness", "bg_rms", "bg_mean", "bg_maxabs",
            "spline_peak_over_sample", "t_pk", "t50", "iint_wa", "iint_phase_pp", "area",
            "onset_shift_ns", "T_data", "tau_slow", "tau_used"] \
        + ["cf_" + c for c in cutflow_names]
    nan = np.nan

    def fmt(v):
        if isinstance(v, (float, np.floating)):
            return "nan" if not np.isfinite(v) else f"{v:.6g}"
        return str(v)

    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            c = r.get("cand")
            cv = ((c["run"], c["split"], c["event_number"], c["amp"], c["tier"],
                   c["peak_sample"], c["first_thr"], c["roughness"], c["bg_rms"],
                   c["bg_mean"], c["bg_maxabs"]) if c else
                  (-1, -1, -1, nan, -1, -1, -1, nan, nan, nan, nan))
            v = [r["name"], r["module_id"], r["midx"], r["fallback"], r["valid"],
                 SINGLE_TIER_RULE[r["fallback"]], r["n_cand"], *cv,
                 r.get("peak_over_sample", nan), r["t_pk"], r["t50"], r["iint_wa"],
                 r["iint_pp"], r["area"], r["cum_shift"], r["T_data"], r["tau_slow"],
                 r["tau_used"], *[int(x) for x in r["cutflow"]]]
            fh.write("\t".join(fmt(x) for x in v) + "\n")


def _meta_json_single(ctx):
    """_meta_json with the rules that do not apply replaced."""
    m = json.loads(_meta_json(ctx))
    m["mode"] = "single_pulse_wassim"
    m["cuts"]["build"]["single_pulse"] = True
    m["amplitude_def"] = ("none fitted: R(k) = (W(k) - bg_mean) / max_k(W(k) - bg_mean) of "
                          "the module's reference candidate, W = samples - WaveAnalyzer ped")
    m["pedestal_def"] = ("WaveAnalyzer pedestal (select); bg_mean = mean of the 20 samples "
                         "before first_thr (select, stored)")
    m["single_pulse"] = dict(
        reference="Waveform_Wassim/best_ref_shape_fit.C, Hycal_prad_wf_analysis.pdf section "
                  "4, eq. 18",
        selection="per module the best valid_p1 candidate over all input files in the "
                  "macro's betterCandidate order (tier desc, amp desc, roughness asc, "
                  "bg_rms asc; ties: first file / first row)",
        template="not-a-knot cubic interpolating spline (scipy make_interp_spline, "
                 "bc_type not-a-knot, as TSpline3 with no end conditions) through R(k) at t = 4k ns, k = 0..99; "
                 "no smoothing, no averaging, no alignment iterations, no pool",
        bank_axis="the data spline is shifted so its 0.5 % onset is at 0 and divided by "
                  "its fine-grid (0.01 ns) maximum (finalize_spline): the bank's unit peak "
                  "is the SPLINE maximum, not the maximum sample; peak_over_sample = spline "
                  "maximum / maximum sample (>= 1, overshoot between samples) is in "
                  "summary.tsv (spline_peak_over_sample)",
        tail_extension="the readout ends ~250 ns after the onset; a*exp(-t/tau) fitted "
                       "(tail_extension, unweighted) on the samples in [140 ns after onset, "
                       "T_data] and added as pseudo-points every 4 ns up to 416 ns, only "
                       "because the bank format needs T_END = 400 ns; the spline stays "
                       "interpolating through them",
        onset_shift_ns="onset of the template axis minus the onset of a two-tau (h50) fit "
                       "to the reference pulse (align_pulses, +-8 ns): check seeds "
                       "t0 = mf_t0 (h50 axis) + onset_shift_ns",
        average="run-type average (fallback 2): " + ctx["avg_rule"],
        no_ampcorr="no f(A): the npz has no ampcorr_* keys (mf_pulse uses f = 1)",
        n_modules_with_reference=int(sum(r["fallback"] < 2 for r in ctx["rows"])),
        n_tier2=int(sum(r["fallback"] == 0 for r in ctx["rows"])))
    m["fallback_rule"] = ("fallback 0: own single reference, tier 2; fallback 1: own single "
                          "reference, tier 1; fallback 2: module without any valid candidate "
                          "-> run-type average (valid stays 1, as in the normal build)")
    m["holdout_rule"] = ("none: the bank is one candidate per module; check's holdout_mask "
                         "(seed, frac) applies to the candidate files given to check")
    for k in ("smoothing_rule", "iteration_rule", "band_rule", "tail_rule",
              "fit_points_rule"):
        m.pop(k, None)
    m["average_rule"] = ctx["avg_rule"]
    return json.dumps(m)


def cmd_build_single(a):
    import _common as C
    import mf_pulse
    bad = [n for n, on in (
        ("--pool-order", a.pool_order != "random"),
        ("--lambda-select", a.lambda_select != "fixed"), ("--osc-escalate", a.osc_escalate),
        ("--bands", a.bands is not None), ("--converge-tol", a.converge_tol is not None),
        ("--iterations", a.iterations != N_IT), ("--fit-points", a.fit_points != "bins"),
        ("--max-iterations", a.max_iterations != MAX_IT)) if on]
    if bad:
        sys.exit(f"[ERROR] build --single-pulse cannot be combined with {', '.join(bad)}")
    t_start = time.monotonic()
    files = sorted({f for pat in a.pulses.split(",") for f in glob.glob(pat.strip())})
    if not files:
        sys.exit(f"[ERROR] no files match {a.pulses}")
    for f in files:
        m0 = json.loads(str(np.load(f, allow_pickle=False)["meta"]))
        if m0.get("mode") != "wassim_ref":
            sys.exit(f"[ERROR] {f}: mode {m0.get('mode')!r} != 'wassim_ref' (build "
                     "--single-pulse needs `select --wassim-ref` files)")
    stem = Path(a.out)
    stem.parent.mkdir(parents=True, exist_ok=True)
    D = load_candidates(files, SINGLE_FIELDS)
    wave = _wave_ctx(D["metas"])
    runs = sorted(set(int(x) for x in D["run"]))
    print(f"[load] {len(files)} files, runs {runs}, {D['midx'].size} reference "
          f"candidates ({time.monotonic() - t_start:.1f} s)", flush=True)

    p = C.setup_pipeline(evio_path="", run_num=runs[0])
    mods = [p.hycal.module(i) for i in range(p.hycal.module_count())]
    names = [m.name for m in mods]
    pwo4 = [bool(m.is_pwo4()) for m in mods]
    rows_idx = [i for i, f in enumerate(pwo4) if f]
    if a.modules:
        want = set(s.strip() for s in a.modules.split(","))
        rows_idx = [i for i in rows_idx if names[i] in want]
    h50 = Path(a.h50)
    if not h50.is_absolute():
        h50 = HERE / h50
    column = mf_pulse.resolve_column(str(h50), "auto", prefer="25320")
    tau_rows, h50_sha1, agg_json = mf_pulse.load_h50(str(h50), column)
    tau_agg = mf_pulse.load_aggregate(agg_json)
    print(f"[setup] h50 {h50.name} column {column} (onset_shift_ns reference only)",
          flush=True)

    glass = [(m.x, m.y, m.size_x, m.size_y) for m in mods
             if m.is_hycal() and not m.is_pwo4()]
    geo_rows = [(mods[i].x, mods[i].y, mods[i].size_x, mods[i].size_y) for i in rows_idx]
    allg = glass + geo_rows
    xs = [x - sx / 2 for x, y, sx, sy in allg] + [x + sx / 2 for x, y, sx, sy in allg]
    ys = [y - sy / 2 for x, y, sx, sy in allg] + [y + sy / 2 for x, y, sx, sy in allg]
    pm = np.array(pwo4)
    ctx = dict(D=D, args=a, iter_kw={}, hold=np.zeros(D["midx"].size, dtype=bool),
               n_hold=None, names=names, tau=tau_rows, tau_agg=tau_agg, band_edges=None,
               fa_bands=FA_BANDS_DEFAULT, wave=wave, run_type=a.run_type, files=files,
               cutflow_names=D["cutflow_names"], cutflow_tot=D["cutflow"][pm].sum(axis=0),
               single=True,
               geom=dict(glass=glass, rows=geo_rows, xlim=(min(xs), max(xs)),
                         ylim=(min(ys), max(ys))))

    best = _best_wassim(D)
    n_per_mod = np.bincount(D["midx"], minlength=len(mods))
    nan_row = np.full(T_GRID.size, np.nan, dtype=np.float32)
    results = []
    for mi in rows_idx:
        res = dict(midx=mi, fallback=2, valid=1, n_cand=0, band=(np.nan, np.nan),
                   n_pool=0, N=(0, 0, 0), cand=None)
        g = best.get(mi)
        if g is not None:
            ref0 = TwoTauShape(*tau_rows.get(names[mi], tau_agg))
            w = wassim_template(D["samples"][g], float(D["ped"][g]), float(D["bg_mean"][g]),
                                float(D["bg_rms"][g]), ref0)
            fin, tail = w["fin"], w["tail"]
            d = derived_from_spline(fin["t"], fin["c"], wave, with_iint=False)
            amp = float(D["amp"][g])
            bx, by = w["x"], w["y"]
            bsg = np.full(bx.size, max(float(D["bg_rms"][g]), SIG_MIN) / amp
                          / w["peak_over_sample"])
            brs = by - spline_values(fin["t"], fin["c"], bx)
            res.update(
                fallback=0 if int(D["tier"][g]) == 2 else 1, n_cand=int(n_per_mod[mi]),
                band=(amp, amp), n_pool=1, N=(1, 1, 1),
                cand=dict(run=int(D["run"][g]), split=int(D["split"][g]),
                          event_number=int(D["event_number"][g]), amp=amp,
                          tier=int(D["tier"][g]), peak_sample=int(D["peak_sample"][g]),
                          first_thr=int(D["first_thr"][g]),
                          roughness=float(D["roughness"][g]), bg_rms=float(D["bg_rms"][g]),
                          bg_mean=float(D["bg_mean"][g]),
                          bg_maxabs=float(D["bg_maxabs"][g])),
                shape=spline_values(fin["t"], fin["c"], T_GRID).astype(np.float32),
                err=nan_row, resid=nan_row, t=fin["t"], c=fin["c"], shift=w["t_on1"],
                shift0=w["t_on0"], cum_shift=w["shift"], vmax=fin["vmax"],
                peak_over_sample=w["peak_over_sample"], T_data=tail["T_data"],
                tau_slow=tail["tau_slow"], a_slow=tail["a_slow"],
                tau_used=tail["tau_used"], a_used=tail["a_used"], knot_scale=0,
                n_knots=int(fin["t"].size), n_iterations=0, **d,
                pts=(np.float32(bx), np.float32(by)),
                bins=(bx, by, bsg, np.ones(bx.size), brs))
            res["rms_pre"] = float(np.sqrt(np.mean(by[bx < -2.0] ** 2))) \
                if (bx < -2.0).any() else np.nan
        results.append(res)
    own_rows = [r for r in results if r["fallback"] < 2]
    if not own_rows:
        sys.exit("[ERROR] no module with a reference candidate: cannot build the "
                 "run-type average")
    sh = np.array([r["cum_shift"] for r in own_rows])
    sh_med = float(np.median(sh[np.isfinite(sh)])) if np.isfinite(sh).any() else 0.0
    for r in own_rows:
        if not np.isfinite(r["cum_shift"]):
            r["cum_shift"] = sh_med                 # two-tau fit failed: median shift
    n_t2 = sum(r["fallback"] == 0 for r in results)
    print(f"[single] {len(own_rows)}/{len(results)} modules with a reference "
          f"({n_t2} tier 2, {len(own_rows) - n_t2} tier 1), "
          f"{len(results) - len(own_rows)} without (run-type average); two-tau fit "
          f"failed for {int((~np.isfinite(sh)).sum())}; onset_shift_ns median "
          f"{sh_med:.2f} ns; spline peak / max sample median "
          f"{np.median([r['peak_over_sample'] for r in own_rows]):.4f}, max "
          f"{max(r['peak_over_sample'] for r in own_rows):.4f}", flush=True)

    # run-type average: per-sample median of the own shapes (exponential beyond 400 ns)
    tau_av = float(np.median([r["tau_used"] for r in own_rows]))
    shp = np.nanmedian(np.stack([r["shape"] for r in own_rows]).astype(np.float64), axis=0)
    xa = T_GRID[::SINGLE_AVG_STEP]
    ya = np.nan_to_num(shp[::SINGLE_AVG_STEP])
    xe = np.arange(xa[-1] + 0.5, T_AXIS_END + 1e-9, 0.5)
    avg_fin = finalize_spline(*_interp_spline(
        np.concatenate([xa, xe]),
        np.concatenate([ya, ya[-1] * np.exp(-(xe - xa[-1]) / tau_av)])), shift=True)
    avg_der = derived_from_spline(avg_fin["t"], avg_fin["c"], wave, with_iint=False)
    avg_shape = spline_values(avg_fin["t"], avg_fin["c"], T_GRID).astype(np.float32)
    avg_tail = dict(T_data=float(np.median([r["T_data"] for r in own_rows])),
                    tau_slow=float(np.nanmedian([r["tau_slow"] for r in own_rows])),
                    a_slow=float(np.nanmedian([r["a_slow"] for r in own_rows])),
                    tau_used=tau_av,
                    a_used=float(np.median([r["a_used"] for r in own_rows])))
    avg_rule = ("per-sample median of the own modules' single-pulse shapes (T_GRID every "
                "0.5 ns, exponential tail beyond 400 ns), not-a-knot cubic interpolating "
                "spline, re-finalised to unit peak / onset 0")
    av = dict(sp=dict(tail=avg_tail), N=(0,), shifts=[avg_fin["t_on"]],
              final=dict(chi2dof=np.array([np.nan])), fin=avg_fin)
    avg_t, avg_c = avg_fin["t"], avg_fin["c"]
    ctx["avg_shape_fn"] = (lambda tt: np.where(
        (tt > 0) & (tt < T_END),
        spline_values(avg_t, avg_c, np.clip(tt, avg_t[3], avg_t[-4])), 0.0))
    print(f"[average] {len(own_rows)} module shapes, onset shift {avg_fin['t_on']:.4f} ns",
          flush=True)

    rows = []
    for r in results:
        mi = r["midx"]
        r["name"], r["module_id"] = names[mi], int(mods[mi].id)
        r["cutflow"] = D["cutflow"][mi]
        if r["fallback"] == 2:
            r.update(shape=avg_shape, err=nan_row, resid=nan_row, t=avg_fin["t"],
                     c=avg_fin["c"], shift=float(avg_fin["t_on"]),
                     shift0=float(avg_fin["t_on"]), cum_shift=float(avg_fin["t_on"]),
                     vmax=avg_fin["vmax"], T_data=avg_tail["T_data"],
                     tau_slow=avg_tail["tau_slow"], a_slow=avg_tail["a_slow"],
                     tau_used=avg_tail["tau_used"], a_used=avg_tail["a_used"],
                     pts=None, bins=None, **avg_der)
        sp_t, sp_c = r["t"], r["c"]
        r["shape_fn"] = (lambda tt, t=sp_t, c=sp_c:
                         np.where((tt > 0) & (tt < T_END),
                                  spline_values(t, c, np.clip(tt, t[3], t[-4])), 0.0))
        rows.append(r)
    ctx.update(p=p, rows=rows, av=av, avg_rule=avg_rule, n_avg=len(own_rows),
               avg_fin=avg_fin, avg_shape=avg_shape, avg_der=avg_der, avg_tail=avg_tail,
               h50=h50, h50_sha1=h50_sha1, column=column, agg_json=agg_json, mods=mods,
               pwo4=pwo4, runs=runs, t_start=t_start, stem=stem)
    return finish_build(ctx)


def _meta_json(ctx):
    a, D = ctx["args"], ctx["D"]
    cuts = dict(select=cut_values(), build=dict(
        n_min=a.n_min, n_min_relaxed=a.n_min_relaxed, n_max=a.n_max, seed=a.seed,
        holdout_frac=a.holdout_frac, band_b0=list(BAND_B0), band_relaxed=list(RELAXED_BAND),
        band_lo=list(BAND_LO), band_hi=list(BAND_HI), band_mid=list(BAND_MID),
        band_top=list(BAND_TOP), band_min_pulses=BAND_MIN,
        **({} if ctx["band_edges"] is None else dict(fa_band_edges=list(ctx["band_edges"]))),
        align_window_ns=list(WIN), pool_tau_min_ns=POOL_TAU_MIN, bin_ns=BIN_W,
        drop_frac=DROP_FRAC, n_iterations=N_IT, sigma_min_adc=SIG_MIN,
        sigma_t_ns=SIGMA_T_NS, smooth_lambda=SMOOTH_LAMBDA, onset_frac=ONSET_FRAC,
        **({} if not ctx["iter_kw"] else dict(iteration_config=dict(
            iterations=a.iterations, converge_tol_ns=a.converge_tol,
            max_iterations=a.max_iterations, freeze_from_iteration=FREEZE_IT,
            conv_max_dshape=CONV_DF, conv_chi2_frac=CONV_CHI2,
            applies_to="module templates (process_module); band templates (2 "
                       "iterations) and the run-type average unchanged",
            fixed_arrays_note="n_iter (rows, 3) and chi2_iter (rows, 4) hold "
                              "iterations 0, 1 and the LAST iteration (chi2_iter: "
                              "then final); iter_hist [rows, M, 6] (columns "
                              f"{list(ITER_COLS)}, NaN-padded), n_iterations, "
                              "converged (-1 without --converge-tol) hold all"))),
        **({} if a.fit_points == "bins" else dict(fit_points=a.fit_points)),
        **({} if not SHAPE_SIGMA > 0 else dict(shape_sigma=SHAPE_SIGMA)),
        **({} if not a.osc_escalate else dict(osc_escalate=dict(
            prominence=OSC_PROM, window=list(OSC_WIN), excl=OSC_EXCL,
            peak_window=list(OSC_PEAK_WIN),
            rule="per module template: count the maxima + minima of the shape on t_grid in "
                 "window with prominence >= prominence * (largest value in peak_window), "
                 "excluding |t - t_peak| <= excl; while count > 0 and lambda < 1e4: "
                 "lambda = 1e1 if lambda <= 0 else 10 * lambda (cap 1e4) and rebuild the "
                 "same pool (full re-iteration); still > 0 at 1e4 (or the rebuild failed): "
                 "the run-type average (fallback 2); lam_final = lambda used"))),
        **({} if a.lambda_select != "heldout" else dict(lambda_select=dict(
            mode="heldout", grid=[x for x in LAM_GRID if x >= LAM_MIN],
            fallback_lambda=LAM_FALLBACK, min_heldout=LAM_MIN_HO,
            chi2_sigma=LAM_CHI2_SIGMA,
            rule="per module, the lambda of the grid with the lowest mean chi2/dof "
                 "(sigma per chi2_sigma) of the module's held-out pulses aligned against the "
                 "final template (pulses ok for every lambda; ties: first); fewer than "
                 "min_heldout held-out pulses: fallback_lambda; null lam = fallback 2 "
                 "(run-type average, smooth_lambda); lists aligned to the npz rows",
            lam=[None if not np.isfinite(r["lam_final" if a.osc_escalate else "lam"])
                 else float(r["lam_final" if a.osc_escalate else "lam"])
                 for r in ctx["rows"]],
            lam_flag=[int(r["lam_flag"]) for r in ctx["rows"]],
            n_ho=[int(r["n_ho"]) for r in ctx["rows"]])))))
    return json.dumps(dict(
        schema="spline_templates_v1", run_type=ctx["run_type"], runs=ctx["runs"],
        splits=[[i["run"], i["split"]] for i in D["inputs"]], inputs=D["inputs"],
        modules_subset=a.modules or None, cuts=cuts,
        trigger_rule=D["metas"][0].get("trigger_rule"),
        trigger_note="see trigger_rule (select files)",
        pedestal_def="WaveAnalyzer pedestal/rms (wave_ana.analyze); y = raw - ped",
        amplitude_def="A fitted per pulse (closed form, least squares against the "
                      "current template over tau in [-20, 160] ns); y/A pooled",
        wave_config=ctx["wave"], mf_seed=D["metas"][0].get("mf"),
        h50=dict(path=str(ctx["h50"]), sha1=ctx["h50_sha1"], column=ctx["column"],
                 agg_json=ctx["agg_json"], aggregate=list(ctx["tau_agg"])),
        onset_rule="spline shifted so the last time before the peak at which it is "
                   "<= 0.5 % of its (fine-grid, 0.01 ns) maximum is t = 0; divided "
                   "by that maximum; BSpline(t) is the unit-peak template for "
                   "0 < t < 400 ns (0 outside; t < 0 holds the pre-onset spline)",
        grid_convention="t_grid = -8..400 ns step 0.1 on the shifted axis; shape, "
                        "shape_err (per-bin sigma, interpolated, NaN beyond T_data), "
                        "resid (bin median - spline, fraction of peak, NaN where no "
                        "bin) are float32",
        tail_rule="no pin at the end; T_data = last 0.5 ns bin with >= 5 points; "
                  "a*exp(-t/tau) weighted-fitted to the bin medians on [140, T_data] "
                  "ns; pseudo-points from it every 4 ns beyond T_data to 416 ns "
                  "(pre-shift axis) with the median bin weight; knots 0.5 or 1 ns "
                  "on [-1,10], 1 ns (10,20], 2 ns (20,40], 4 ns (40,80], 8 ns "
                  "(80,160], 16 ns (160,400] (+ one knot at -10 ns), doubling on "
                  "Schoenberg-Whitney failure",
        smoothing_rule="penalised least squares: sum((y-f)^2/sigma^2) + lam * "
                       f"int f''(t)^2 dt with lam = {SMOOTH_LAMBDA:g} (f in peak units, t in ns); "
                       "needed because the sub-sample phase of a module's pulses "
                       "is concentrated (phase_R); per-bin sigma floor "
                       "sqrt(noise^2 + (0.15 ns * slope)^2)/sqrt(n); within-bin "
                       "detrend (y - (ref(tau) - ref(bin centre)))",
        iteration_rule="iteration 0 aligns against the h50 two-tau (+-3 ns of the "
                       "mf seed), pools all pulses; iterations 1, 2 align against "
                       "the previous sampled template (+-2 ns of the previous t0 "
                       "on the new axis), drop the worst 20 % of ALL pulses by "
                       "chi2/dof (sigma = max(rms, 0.3)), re-pool",
        band_rule="b_lo/b_hi: held-in tier>=1 pulses of the band, two iterations "
                  "(plan 6.2 step 7): iteration 0 aligns them against the module's "
                  "final template (+-3 ns) and pools all; iteration 1 aligns against "
                  "the band's own iteration-0 spline (+-2 ns) and drops the worst "
                  "20 % by chi2/dof against it; same axis as the module (no onset "
                  "shift), same binning and spline",
        average_rule=ctx["avg_rule"], average_pulses=ctx["n_avg"],
        fallback_rule="fallback 0: >= n_min held-in candidates with tier >= 1 in "
                      "b0 [300, 2500]; fallback 1: >= n_min_relaxed with tier >= 1 in "
                      "[100, 2500]; fallback 2: run-type average (the tier stays "
                      "stored per candidate; select has no PASS 1 shape cut, "
                      "p1_shape_ok is stored)",
        integral_rule="iint_wa = mean over 256 sub-sample phases of the exact "
                      "discrete walk-out integral on the template (mf_pulse "
                      "SampledTemplateBank.I_cont, same code); iint_phase_pp = "
                      "(max - min) of the 256-phase table / iint_wa",
        ampcorr=ctx["ampcorr"],
        ampcorr_note="the .npz keys ampcorr_logA / ampcorr_f / ampcorr_meta are "
                     "authoritative (the `ampcorr` subcommand may replace the "
                     "per-run-type values by one global correction in the .npz; "
                     "this copy is the per-run-type value at build time)",
        p1_shape_ok_pooled=ctx["p1ok"],
        holdout_rule=("holdout_mask(run, split, event_number, midx, seed, frac): "
                      "splitmix64 hash of the five integers; True = held out"
                      if a.pool_order == "random" else
                      "energy: per module the eligible candidates (tier>=1, band "
                      "b0, finite mf_energy) sorted by mf_energy descending "
                      "(stable on the candidate index); ranks [0, n_max) = fit "
                      "pool (never held out), the next holdout_n_per_module = "
                      "int(n_max*frac/(1-frac)+0.5) ranks are held out, the "
                      "rest is not held out; the list is the npz array "
                      "holdout_keys (run, split, event_number, module_index)"),
        **({} if a.pool_order == "random" else dict(
            pool_order=a.pool_order, holdout_n_per_module=ctx["n_hold"],
            pool_rule="fit pool = top n_max by mf_energy of the module's "
                      "eligible candidates (the rng draw is still made)"
                      + ("; energy-phase: taken in descending mf_energy with at most "
                         "ceil(n_max/8) per bin of mf_t0 mod 4 ns (8 equal bins), the "
                         "rest filled by mf_energy" if a.pool_order == "energy-phase" else ""))),
        holdout_seed=a.seed, holdout_frac=a.holdout_frac,
        root_layout=dict(nbins=ROOT_NBINS, onset_bin=ROOT_ONSET_BIN,
                         dt_bin_ns=ROOT_DT_BIN,
                         module_type=f"{ROOT_MODULE_TYPE} = EventData.h MOD_PbWO4 "
                                     "(what hycal.module_type holds in _raw.root)"),
        **({} if a.fit_points == "bins" else dict(
            fit_points_rule="direct: the spline is fitted to the pooled points "
                            "(tau <= T_data; y/A, sigma = sqrt((sigma_pulse/A)^2 + "
                            "(0.15 ns * ref'(tau))^2)) instead of the 0.5 ns bin "
                            "medians; tail extension from the bin stats as in bins")),
        n_rows=len(ctx["rows"]), git=_git_short(),
        created=datetime.now().isoformat(timespec="seconds"),
        python=sys.version.split()[0], numpy=np.__version__,
        scipy=__import__("scipy").__version__, argv=sys.argv))


AMPCORR_RULE = ("E = A * I_int * f(A); f = PCHIP (monotone cubic) in ln A through the "
                "three band knots (lo 100-400, mid 400-1000, top 1000-2500; knot at the "
                "band's median amplitude), stored as a dense table (8 points per "
                "interval, knots included) so the linear np.interp of mf_pulse is "
                "smooth; constant outside the outer knots; f_b = <C++ walk-out integral "
                "/ mf energy> of noise-free band-template waveforms fitted with the "
                "module template, relative to the same for the module template; "
                "16 phases; knots chained within modules: f_lo = median(r_lo / r_mod), "
                "f_mid = f_lo * median(r_mid / r_lo), f_top = f_mid * median(r_top / "
                "r_mid) (chosen on the production banks, see plan progress log)")
FA_BANDS = (("lo", BAND_LO), ("mid", BAND_MID), ("top", BAND_TOP))
FA_BANDS_DEFAULT = [dict(key=k, lo=lo, hi=hi, incl=(k == "top")) for k, (lo, hi) in FA_BANDS]


def fa_bands_from_edges(edges):
    """--bands E0,...,En -> the n f(A) bands [E_k, E_k+1) (the last one closed),
    keyed "E_k-E_k+1" (that label is the knot's identity in ampcorr_meta)."""
    return [dict(key=_bin_label(lo, hi), lo=float(lo), hi=float(hi),
                 incl=(k == len(edges) - 2))
            for k, (lo, hi) in enumerate(zip(edges[:-1], edges[1:]))]


def ampcorr_rule(edges):
    """AMPCORR_RULE for the default bands, its generalisation for --bands."""
    if edges is None:
        return AMPCORR_RULE
    return (f"E = A * I_int * f(A); f = PCHIP (monotone cubic) in ln A through the band "
            f"knots (bands {', '.join(_bin_label(lo, hi) for lo, hi in zip(edges[:-1], edges[1:]))} "
            f"ADC; a band with no usable module has no knot; knot at the band's median "
            f"amplitude), stored as a dense table (8 points per interval, knots included) "
            f"so the linear np.interp of mf_pulse is smooth; constant outside the outer "
            f"knots; f_b = <C++ walk-out integral / mf energy> of noise-free band-template "
            f"waveforms fitted with the module template, relative to the same for the "
            f"module template; 16 phases; knots chained within modules: f_b0 = "
            f"median(r_b0 / r_mod), f_bk = f_b(k-1) * median(r_bk / r_b(k-1)) over the "
            f"modules having both bands")
FA_PER_INTERVAL = 8                     # dense-table points per knot interval (17 for 3 knots)


def fa_table(knots):
    """knots {band: dict(A, f)} (the bands that exist) -> (logA, f) lists of the
    stored f(A) table, or None if fewer than two usable knots.  3 knots: PCHIP in
    ln A through them, sampled at FA_PER_INTERVAL equal ln A steps per interval (the
    knots are table points; 17 points for 3 knots); 2 knots: the two points."""
    pts = sorted((float(k["A"]), float(k["f"])) for k in knots.values()
                 if np.isfinite([k["A"], k["f"]]).all() and k["A"] > 0)
    if len(pts) < 2 or np.any(np.diff([p[0] for p in pts]) <= 0):
        return None
    la = np.log([p[0] for p in pts])
    fv = np.array([p[1] for p in pts])
    if len(pts) == 2:
        return la.tolist(), fv.tolist()
    from scipy.interpolate import PchipInterpolator
    pc = PchipInterpolator(la, fv)
    grid = np.concatenate([np.linspace(la[i], la[i + 1], FA_PER_INTERVAL + 1)[:-1]
                           for i in range(len(la) - 1)] + [la[-1:]])
    return grid.tolist(), pc(grid).tolist()
FA_NPH = 16                             # phases of the f(A) waveforms (as G1-6)
FA_MF = dict(nsigma=5.0, min_amp=5.0, kmax=4, r_model=0.3)   # production working point


def _median_se(x):
    x = np.asarray(x, dtype=np.float64)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return float("nan"), float("nan"), 0
    se = 1.2533 * x.std(ddof=1) / np.sqrt(x.size) if x.size > 1 else float("nan")
    return float(np.median(x)), float(se), int(x.size)


def _fa_ratios(ctx, bank, sets, mod_rows, own, band_i, A_b, H_b0):
    """r_b / r_mod: noise-free integer waveforms ped 200 + H shape_b(4 i - t0),
    t0 = 140 + phi, FA_NPH phases, of the band template b (H = A_b[b]; rows of
    sets[b]) or of the module template (H = H_b0; rows of mod_rows), set in the
    rows as r_<b> / r_mod;  r = phase mean of Q_raw / E_mf, Q_raw = C++ wave_ana
    single-peak integral, E_mf = mf_pulse energy (module template, no ampcorr,
    FA_MF) of the single fitted pulse.  Waveforms without exactly one peak (C++ or
    mf) are skipped and counted; a module with none left has NaN.  Returns the
    skip counts."""
    import mf_pulse
    ana = ctx["p"].wave_ana.analyze
    ped_n = int(ctx["p"].wave_ana.cfg.ped_nsamples)
    mf = mf_pulse.MFReconstructor(bank, ped_nsamples=ped_n, **FA_MF)
    phi = DT_NS * np.arange(FA_NPH) / FA_NPH
    skipped = {}
    todo = [(b["key"], sets[b["key"]], A_b[b["key"]]) for b in ctx["fa_bands"]
            if sets.get(b["key"])]
    todo.append(("mod", mod_rows, H_b0))
    for key, rr, H in todo:
        nm = len(rr)
        own_k = np.array([own[r["midx"]] for r in rr], dtype=np.int64)
        gen = own_k if key == "mod" else np.array(
            [band_i[(r["midx"], key)] for r in rr], dtype=np.int64)
        gk = np.repeat(gen, FA_NPH)
        t0 = 140.0 + np.tile(phi, nm)
        S = bank.shape(gk[:, None], DT_NS * np.arange(NS_EXPECT)[None, :] - t0[:, None],
                       dtype=np.float64)
        Y = np.rint(200.0 + H * S).astype(np.uint16)
        Q = np.full(Y.shape[0], np.nan)
        for j in range(Y.shape[0]):
            pks = ana(Y[j])[2]
            if len(pks) == 1:
                Q[j] = pks[0].integral
        res = mf.run(Y, np.repeat(own_k, FA_NPH))
        E = np.full(Y.shape[0], np.nan)
        one = np.nonzero(res.npk == 1)[0]
        if one.size:
            pos = np.searchsorted(res.p_ch, one)       # p_ch is sorted by channel
            E[one] = mf_pulse.pulse_energy_units(
                bank, own_k[one // FA_NPH], res.p_t0[pos], res.p_amp[pos], kappa=1.0)
        skipped[key] = dict(cpp_not_one_peak=int(np.isnan(Q).sum()),
                            mf_not_one_peak=int(np.isnan(E).sum()), n=int(Y.shape[0]))
        ratio = (Q / E).reshape(nm, FA_NPH)
        okp = np.isfinite(ratio)
        rb = np.where(okp.any(axis=1), np.where(okp, ratio, 0.0).sum(axis=1)
                      / np.maximum(okp.sum(axis=1), 1), np.nan)
        for r, v in zip(rr, rb):
            r["r_" + key] = float(v)
    return skipped


def finalize_integrals(ctx):
    """iint_wa / iint_phase_pp of every row, the lo / hi band integrals (G1-9),
    r_lo / r_mid / r_top / r_mod (Q_raw / E_mf of noise-free waveforms, _fa_ratios)
    and the run type's f(A) -- all from one in-memory SampledTemplateBank
    (bank_integrals)."""
    rows, av = ctx["rows"], ctx["av"]
    t0 = time.monotonic()
    shapes, own, avg_i = [], {}, None
    for r in rows:
        if r["fallback"] < 2:
            own[r["midx"]] = len(shapes)
            shapes.append((r["t"], r["c"]))
    avg_i = len(shapes)
    shapes.append((ctx["avg_fin"]["t"], ctx["avg_fin"]["c"]))
    band_i = {}
    for r in rows:
        for key in ("lo", "hi", *[b["key"] for b in ctx["fa_bands"] if b["key"] != "lo"]):
            if "spl_" + key in r:
                band_i[(r["midx"], key)] = len(shapes)
                shapes.append(r["spl_" + key])
    I, pp, bank = bank_integrals(shapes, ctx["wave"], return_bank=True)
    for r in rows:
        k = own.get(r["midx"], avg_i)
        r["iint_wa"], r["iint_pp"] = float(I[k]), float(pp[k])
        for key in ("lo", "hi"):
            if (r["midx"], key) in band_i:
                r["iint_" + key] = float(I[band_i[(r["midx"], key)]])
        if "iint_lo" in r and "iint_hi" in r:
            r["lohi_iint_ratio"] = r["iint_hi"] / r["iint_lo"]
    ctx["avg_der"]["iint_wa"] = float(I[avg_i])
    own_rows = [r for r in rows if r["fallback"] < 2]
    fa_keys = [b["key"] for b in ctx["fa_bands"]]
    sets = {b: [r for r in own_rows if "spl_" + b in r] for b in fa_keys}
    mod_rows = [r for r in own_rows if any("spl_" + b in r for b in fa_keys)]
    ac = dict(rule=ampcorr_rule(ctx["band_edges"]), scope="per_run_type",
              run_type=ctx["run_type"],
              bands=(dict(lo=list(BAND_LO), mid=list(BAND_MID), top=list(BAND_TOP),
                          hi_g19_only=list(BAND_HI)) if ctx["band_edges"] is None else
                     dict({b["key"]: [b["lo"], b["hi"]] for b in ctx["fa_bands"]},
                          hi_g19_only=list(BAND_HI))), band_min_pulses=BAND_MIN,
              n_modules=len(mod_rows), knots={}, f=None, logA=None, per_run_type=None,
              se_rule="1.2533 * std(ddof=1) / sqrt(n) of the per-module ratios",
              integral_def="I = mean over 256 sub-sample phases of the exact discrete "
                           "walk-out integral on the template (mf_pulse "
                           "SampledTemplateBank.I_cont)")
    if mod_rows:
        A_b = {b: float(np.median(np.concatenate([r["amp_" + b] for r in sets[b]])))
               for b in fa_keys if sets[b]}
        H_b0 = float(np.median(np.concatenate([r["amp_b0"] for r in mod_rows])))
        skipped = _fa_ratios(ctx, bank, sets, mod_rows, own, band_i, A_b, H_b0)
        ac.update(H_b0=H_b0, n_skipped_waveforms=skipped,
                  fa_mf=dict(FA_MF), fa_phases=FA_NPH,
                  fa_note=(
                      "per-module r_b = mean over phases of Q_raw / E_mf (b = lo, "
                      "mid, top, mod); f_lo = median over the modules having lo "
                      "(and a finite r_mod) of r_lo / r_mod; f_mid = f_lo * median "
                      "over the modules having lo and mid of r_mid / r_lo; f_top = "
                      "f_mid * median over the modules having mid and top of r_top / "
                      "r_mid (chained within-module steps; chosen on the production "
                      "banks, see plan progress log); SEs of the steps are combined "
                      "in quadrature of relative errors, n_modules = size of the "
                      "step's module set; npz r_lo / r_mid / r_top / r_mod hold the r_b"
                      if ctx["band_edges"] is None else
                      "per-module r_b = mean over phases of Q_raw / E_mf (b = the "
                      "bands, mod); f_b0 = median over the modules having b0 (and a "
                      "finite r_mod) of r_b0 / r_mod; f_bk = f_b(k-1) * median over the "
                      "modules having b(k-1) and bk of r_bk / r_b(k-1) (chained "
                      "within-module steps; a band with no such module ends the "
                      "chain); SEs of the steps are combined in quadrature of relative "
                      "errors, n_modules = size of the step's module set; npz r_band "
                      "[nb, nm] and r_mod hold the r_b"))
        fin = lambda r, *ks: all(np.isfinite(r.get("r_" + k, np.nan)) for k in ks)
        chain = tuple((b, "mod" if k == 0 else fa_keys[k - 1], None if k == 0 else fa_keys[k - 1])
                      for k, b in enumerate(fa_keys))
        for b, den, prev in chain:
            if b not in A_b:
                print(f"[WARN] f(A): no module with a usable {b} band", flush=True)
                break
            rb = [r for r in sets[b] if fin(r, b, den)]
            if not rb:
                print(f"[WARN] f(A): no module with a usable {b} band", flush=True)
                break
            step, step_se, n_b = _median_se([r["r_" + b] / r["r_" + den] for r in rb])
            if prev is None:
                f_b, se_b = step, step_se
            else:
                fp = ac["knots"][prev]
                f_b = fp["f"] * step
                se_b = f_b * float(np.hypot(fp["se"] / fp["f"], step_se / step))
            ac["knots"][b] = dict(A=A_b[b], f=f_b, se=se_b, n_modules=n_b)
    tab = fa_table(ac["knots"]) if ac["knots"] else None
    if tab is None:
        print("[WARN] fewer than two usable f(A) band knots: no f(A) stored", flush=True)
    else:
        if len(ac["knots"]) < len(fa_keys):
            print(f"[WARN] f(A) from {len(ac['knots'])} band knots only "
                  f"({', '.join(ac['knots'])})", flush=True)
        ac.update(logA=tab[0], f=tab[1], n_table=len(tab[0]),
                  per_run_type=dict(knots=ac["knots"], n_modules=ac["n_modules"]))
    ctx["ampcorr"] = ac
    ctx["p1ok"] = dict(n_pool=int(sum(r["n_pool"] for r in own_rows)),
                       n_ok=int(sum(r["n_p1ok"] for r in own_rows)))
    print(f"[integrals] {len(shapes)} shapes through one in-memory bank "
          f"({time.monotonic() - t0:.0f} s); phase p-p median "
          f"{100 * np.median([r['iint_pp'] for r in rows]):.2f} %, max "
          f"{100 * max(r['iint_pp'] for r in rows):.2f} %", flush=True)
    if ac["knots"]:
        print("[ampcorr] knots: " + ";  ".join(
            f"{b} A {k['A']:.0f} f {k['f']:.4f} +- {k['se']:.4f} (n {k['n_modules']})"
            for b, k in ac["knots"].items()) + f";  {ac['n_modules']} modules, "
            f"table {ac.get('n_table', 0)} points", flush=True)
        print(f"[ampcorr] f_b chained ({', '.join(f'r_{b} / r_{d}' for b, d, _ in chain)}); "
              f"H_b0 {ac['H_b0']:.0f}; skipped waveforms "
              f"(C++ / mf not one peak): "
              + ", ".join(f"{k} {v['cpp_not_one_peak']} / {v['mf_not_one_peak']} of {v['n']}"
                          for k, v in ac["n_skipped_waveforms"].items()), flush=True)
    pk = ctx["p1ok"]
    print(f"[p1_shape_ok] {pk['n_ok']}/{pk['n_pool']} pooled pulses "
          f"({100 * pk['n_ok'] / max(pk['n_pool'], 1):.1f} %) pass the PASS 1 shape "
          f"cut (informational)", flush=True)


def finish_build(ctx):
    import mf_pulse
    a, D, rows, stem = ctx["args"], ctx["D"], ctx["rows"], ctx["stem"]
    single = bool(ctx.get("single"))
    (finalize_integrals_single if single else finalize_integrals)(ctx)
    nr, nt = len(rows), T_GRID.size
    g = lambda k, dt=float: np.array([r.get(k, np.nan) for r in rows], dtype=dt)   # noqa: E731
    nan_row = np.full(nt, np.nan, dtype=np.float32)
    nk = max(r["t"].size for r in rows)
    nc = max(r["c"].size for r in rows)
    avg_fin = ctx["avg_fin"]
    cf_sel = np.array([r["cutflow"] for r in rows])
    arr = dict(
        names=np.array([r["name"] for r in rows]),
        module_id=np.array([r["module_id"] for r in rows], dtype=np.int64),
        module_index=np.array([r["midx"] for r in rows], dtype=np.int64),
        valid=np.array([r["valid"] for r in rows], dtype=np.int8),
        fallback=np.array([r["fallback"] for r in rows], dtype=np.int8),
        tier_rule=np.array([(SINGLE_TIER_RULE if single else RELAX_TIER_RULE)[r["fallback"]]
                            for r in rows]),
        band_lo=np.array([r["band"][0] for r in rows]),
        band_hi=np.array([r["band"][1] for r in rows]),
        n_cand=np.array([r["n_cand"] for r in rows], dtype=np.int32),
        n_pool=np.array([r["n_pool"] for r in rows], dtype=np.int32),
        n_iter=np.array([r["N"] for r in rows], dtype=np.int32).reshape(nr, 3),
        chi2_ndf=g("chi2"),
        chi2_iter=np.array([list(r["chi2_it"]) + [r["chi2"]] if "chi2_it" in r
                            else [np.nan] * 4 for r in rows]),
        t_grid=T_GRID,
        shape=np.stack([r["shape"] for r in rows]),
        shape_err=np.stack([r["err"] for r in rows]),
        resid=np.stack([r["resid"] for r in rows]),
        spline_k=np.int64(3),
        spline_t=_pad([r["t"] for r in rows], nk),
        spline_c=_pad([r["c"] for r in rows], nc),
        spline_nt=np.array([r["t"].size for r in rows], dtype=np.int32),
        spline_nc=np.array([r["c"].size for r in rows], dtype=np.int32),
        avg_spline_t=avg_fin["t"], avg_spline_c=avg_fin["c"],
        onset_shift_ns=g("cum_shift"),     # total: final axis vs the mf seed axis (fb 2: the average's, vs the modules' final axes)
        t_pk=g("t_pk"), t50=g("t50"), iint_wa=g("iint_wa"),
        iint_phase_pp=g("iint_pp"), iint_lo=g("iint_lo"), iint_hi=g("iint_hi"),
        r_mod=g("r_mod"),
        area=g("area"),   # r_*: Q_raw / E_mf of the f(A) waveforms (_fa_ratios), not iint ratios
        T_data=g("T_data"), tau_slow=g("tau_slow"), a_slow=g("a_slow"),
        knot_scale=np.array([r.get("knot_scale", 0) for r in rows], dtype=np.int32),
        phase_R=g("phase_R"), rise_empty_bins=np.array(
            [r.get("rise_empty", -1) for r in rows], dtype=np.int32),
        shape_lo=np.stack([r.get("shape_lo", nan_row) for r in rows]),
        shape_hi=np.stack([r.get("shape_hi", nan_row) for r in rows]),
        n_lo=np.array([r.get("n_lo", 0) for r in rows], dtype=np.int32),
        n_hi=np.array([r.get("n_hi", 0) for r in rows], dtype=np.int32),
        lohi_maxdiff=g("lohi_maxdiff"), lohi_iint_ratio=g("lohi_iint_ratio"),
        shape_avg=ctx["avg_shape"], T_end=np.float64(T_END),
        cutflow_names=np.array(ctx["cutflow_names"]), cutflow=cf_sel,
        meta=np.array(_meta_json_single(ctx) if single else _meta_json(ctx)))
    if ctx["band_edges"] is None:
        # default bands: the lo / mid / top arrays (the lo / hi ones above are G1-9's)
        arr.update(
            r_lo=g("r_lo"), r_mid=g("r_mid"), r_top=g("r_top"),
            shape_mid=np.stack([r.get("shape_mid", nan_row) for r in rows]),
            shape_top=np.stack([r.get("shape_top", nan_row) for r in rows]),
            n_mid=np.array([r.get("n_mid", 0) for r in rows], dtype=np.int32),
            n_top=np.array([r.get("n_top", 0) for r in rows], dtype=np.int32))
    else:
        # --bands: generalised layout, band k = [band_edges[k], band_edges[k + 1])
        keys = [b["key"] for b in ctx["fa_bands"]]
        arr.update(
            band_edges=np.array(ctx["band_edges"], dtype=np.float64),
            shape_band=np.stack([np.stack([r.get("shape_" + k, nan_row) for r in rows])
                                 for k in keys]),
            r_band=np.stack([g("r_" + k) for k in keys]),
            n_band=np.array([[r.get("n_" + k, 0) for r in rows] for k in keys],
                            dtype=np.int32),
            A_med_band=np.array([[float(np.median(r["amp_" + k])) if "amp_" + k in r
                                  else np.nan for r in rows] for k in keys]))
    ac = ctx["ampcorr"]
    if ac["f"] is not None:
        arr.update(ampcorr_logA=np.array(ac["logA"]), ampcorr_f=np.array(ac["f"]),
                   ampcorr_meta=np.array(json.dumps(ac)))
    if ctx["iter_kw"]:
        mw = max(a.iterations, a.max_iterations) if a.converge_tol is not None \
            else a.iterations
        ih = np.full((nr, mw, len(ITER_COLS)), np.nan, dtype=np.float32)
        for i, r in enumerate(rows):
            if "iter_hist" in r:
                ih[i, :r["iter_hist"].shape[0]] = r["iter_hist"]
        arr.update(iter_hist=ih,
                   n_iterations=np.array([r.get("n_iterations", 0) for r in rows],
                                         dtype=np.int16),
                   converged=np.array([r.get("converged", -1) for r in rows],
                                      dtype=np.int8))
    if a.lambda_select == "heldout":
        arr["lam_module"] = g("lam_final" if a.osc_escalate else "lam")   # lambda used per row (NaN: fallback 2)
    elif a.osc_escalate:
        arr["lam_module"] = g("lam_final")
    if a.pool_order != "random":
        hm = ctx["hold"]
        arr["holdout_keys"] = np.stack(
            [D["run"][hm], D["split"][hm], D["event_number"][hm], D["midx"][hm]],
            axis=1).astype(np.int64)
    out_npz = stem.with_suffix(".npz")
    tmp = out_npz.with_name(out_npz.name + ".tmp")
    with open(tmp, "wb") as fh:
        np.savez_compressed(fh, **arr)
    os.replace(tmp, out_npz)
    print(f"[npz] {out_npz}  {out_npz.stat().st_size / 1e6:.1f} MB", flush=True)

    # --- verification: the file must load as a SampledTemplateBank ------------
    mods, names, pwo4 = ctx["mods"], [m.name for m in ctx["mods"]], ctx["pwo4"]
    bank, info = mf_pulse.build_bank(names, pwo4, str(out_npz), ctx["wave"])
    mi = arr["module_index"]
    tab = bank.I_tab[mi, :mf_pulse.NIPH].astype(np.float64)
    dpp = np.abs((tab.max(axis=1) - tab.min(axis=1)) / bank.I_cont[mi]
                 - arr["iint_phase_pp"])
    print(f"[bank] f(A): {bank.ampcorr_str()};  max |phase p-p diff| {dpp.max():.2e}",
          flush=True)
    di = np.abs(bank.I_cont[mi] - arr["iint_wa"]) / arr["iint_wa"]
    dp = np.abs(bank.t_pk[mi] - arr["t_pk"])
    print(f"[bank] build_bank OK: n_module {info['n_module']}  n_fallback "
          f"{info['n_fallback']}  T_end {bank.T_end:g}  max |I_cont/iint_wa - 1| "
          f"{di.max():.2e}  max |t_pk diff| {dp.max():.3f} ns", flush=True)

    # --- ROOT, tsv, pdf --------------------------------------------------------
    meta_json = str(arr["meta"])
    write_root(str(stem.with_suffix(".root")), ctx, rows, ctx["avg_shape_fn"],
               meta_json)
    if single:
        write_tsv_single(str(stem.parent / "summary.tsv"), rows, ctx["cutflow_names"])
    else:
        write_tsv(str(stem.parent / "summary.tsv"), rows, ctx["cutflow_names"],
                  None if ctx["band_edges"] is None else [b["key"] for b in ctx["fa_bands"]],
                  lam=a.lambda_select == "heldout", osc=a.osc_escalate)
    print(f"[tsv] {stem.parent / 'summary.tsv'}", flush=True)
    if a.png_dir:
        want = set(s.strip() for s in a.png_modules.split(",")) if a.png_modules \
            else {r["name"] for r in rows[:4]}
        write_pngs(a.png_dir, ctx, rows, want,
                   (ctx["av"], ctx["avg_rule"], ctx["n_avg"], a))
    if not a.no_pdf:
        write_pdf(str(stem.with_suffix(".pdf")), ctx, rows, ctx["av"],
                  ctx["avg_rule"], ctx["n_avg"], a)

    # --- terminal summary ------------------------------------------------------
    fb = arr["fallback"]
    print(f"[summary] rows {nr}: fallback 0/1/2 = "
          f"{int((fb == 0).sum())}/{int((fb == 1).sum())}/{int((fb == 2).sum())}")
    if single:
        print(f"[done] total {time.monotonic() - ctx['t_start']:.0f} s", flush=True)
        return 0
    nu = arr["n_pool"][fb < 2]
    if nu.size:
        print(f"[summary] pool size (after n_max): min {nu.min()} median "
              f"{np.median(nu):.0f} max {nu.max()};  n_cand median "
              f"{np.median(arr['n_cand'][fb < 2]):.0f}")
        ci = arr["chi2_iter"][fb < 2]
        print("[summary] median chi2/dof: it0 (two-tau ref) "
              f"{np.nanmedian(ci[:, 0]):.1f}  it1 {np.nanmedian(ci[:, 1]):.1f}  "
              f"{'it2' if not ctx['iter_kw'] else 'it_last'} {np.nanmedian(ci[:, 2]):.1f}  final (kept) "
              f"{np.nanmedian(ci[:, 3]):.1f}")
        ts = arr["tau_slow"][fb < 2]
        ts = ts[np.isfinite(ts)]
        print(f"[summary] tau_slow: n {ts.size}, percentiles 10/50/90 = "
              f"{np.percentile(ts, 10):.0f}/{np.median(ts):.0f}/"
              f"{np.percentile(ts, 90):.0f} ns;  T_data median "
              f"{np.median(arr['T_data'][fb < 2]):.0f} ns")
    print(f"[done] total {time.monotonic() - ctx['t_start']:.0f} s", flush=True)
    return 0


# ================================================================ ampcorr ==
AMPCORR_NSIG = 2.0


def _rewrite_npz(path, arrays):
    """Write `arrays` to `path` through a temp file + rename."""
    out = Path(path)
    tmp = out.with_name(out.name + ".tmp")
    with open(tmp, "wb") as fh:
        np.savez_compressed(fh, **arrays)
    os.replace(tmp, out)


PM_SE_MAX = 0.003                       # --per-module-f: own band ratio kept below this SE
PM_GRID_STEP = 0.02                     # --per-module-f: ln A step of the common table grid


def _pm_nmax(z):
    """The build's per-template pulse cap n_max (band and module templates are fit
    on at most this many pulses; n_band holds the count before the cap)."""
    return int(json.loads(str(z["meta"])).get("cuts", {}).get("build", {})
               .get("n_max", 200))


def _pm_cells(z):
    """Per-module f(A) inputs of a `build --bands` bank: (labels, f_own, A_own, x)
    with f_own = r_b / r_mod [nb, nm] (NaN where either is missing or the module
    has no own template), A_own = the module's median amplitude in band b,
    x = 1 / min(n_b, n_max) + 1 / min(n_pool, n_max) (the variance weight of the
    stat model)."""
    if "band_edges" not in z:
        sys.exit("[ERROR] --per-module-f needs banks built with --bands (band_edges, "
                 "r_band, n_band, A_med_band)")
    be = np.asarray(z["band_edges"], dtype=np.float64)
    labels = [_bin_label(lo, hi) for lo, hi in zip(be[:-1], be[1:])]
    with np.errstate(invalid="ignore", divide="ignore"):
        f_own = np.asarray(z["r_band"]) / np.asarray(z["r_mod"])[None, :]
        nb = np.minimum(np.asarray(z["n_band"], dtype=np.float64), _pm_nmax(z))
        npl = np.minimum(np.asarray(z["n_pool"], dtype=np.float64), _pm_nmax(z))
        x = 1.0 / np.where(nb > 0, nb, np.nan) + 1.0 / npl[None, :]
    A_own = np.asarray(z["A_med_band"], dtype=np.float64)
    ok = np.isfinite(f_own) & np.isfinite(A_own) & (A_own > 0) & np.isfinite(x)
    return labels, np.where(ok, f_own, np.nan), np.where(ok, A_own, np.nan), \
        np.where(ok, x, np.nan)


def _pm_se_model(Z):
    """Stat model of a module's band ratio, calibrated on the banks Z: the step
    ratio s = r_b / r_(b-1) (r_(-1) = r_mod) over the modules of all banks has a
    robust scatter sigma^2 = s0^2 + c^2 (1 / n_b + 1 / n_ref), n_ref = n_(b-1)
    (n_pool for r_mod), every n capped at the build's n_max; c and s0 from a least-squares line through the binned
    (6 quantile bins of x per band and bank) robust variances (1.4826 MAD).  s0 is
    the n-independent scatter (real module differences and what the band template
    does not average out); the stat SE of a module's r_b / r_mod is then
    c sqrt(1 / n_b + 1 / n_pool) (capped).  Returns dict(c, s0, n_bins)."""
    X, Y = [], []
    for z in Z:
        nmax = _pm_nmax(z)
        rb = np.asarray(z["r_band"])
        nb = np.minimum(np.asarray(z["n_band"], dtype=np.float64), nmax)
        rm = np.asarray(z["r_mod"])
        npool = np.minimum(np.asarray(z["n_pool"], dtype=np.float64), nmax)
        for b in range(rb.shape[0]):
            den, nd = (rm, npool) if b == 0 else (rb[b - 1], nb[b - 1])
            m = np.isfinite(rb[b]) & np.isfinite(den) & (nb[b] > 0) & (nd > 0)
            if m.sum() < 30:
                continue
            st = rb[b][m] / den[m]
            x = 1.0 / nb[b][m] + 1.0 / nd[m]
            d = st - np.median(st)
            q = np.quantile(x, np.linspace(0.0, 1.0, 7))
            for lo, hi in zip(q[:-1], q[1:]):
                k = (x >= lo) & (x <= hi)
                if k.sum() < 5:
                    continue
                X.append(x[k].mean())
                Y.append((1.4826 * np.median(np.abs(d[k] - np.median(d[k])))) ** 2)
    if len(X) < 3:
        sys.exit("[ERROR] --per-module-f: too few modules to calibrate the stat model")
    coef = np.linalg.lstsq(np.vstack([np.ones(len(X)), X]).T, np.array(Y), rcond=None)[0]
    if not coef[1] > 0:
        sys.exit(f"[ERROR] --per-module-f: stat model c^2 = {coef[1]:.3g} <= 0")
    return dict(c=float(np.sqrt(coef[1])), s0=float(np.sqrt(max(coef[0], 0.0))),
                n_bins=len(X))


def _pm_table(z, am, own_knots, model):
    """Per-module f(A) of one bank: (logA[ng], f[nm, ng], f_rt[ng], own[nb, nm], se[nb,
    nm], stats).  Cell (module, band) keeps the module's own r_b / r_mod (knot at its
    median band amplitude; the same number as the chain f_b0 = r_b0 / r_mod, f_bk =
    f_b(k-1) r_bk / r_b(k-1) of the run-type rule for a module that has every band
    of the chain) where SE = c sqrt(1 / n_b + 1 / n_pool) < PM_SE_MAX; a band where
    it does not is set from the module's nearest own band b' (lower one on a tie) by
    the run-type step, f = f_own(b') f_rt(b) / f_rt(b'), at the run-type knot A of b,
    so the module keeps its level; a module with no own cell takes the run-type row.  A module's row is PCHIP in ln A through its
    knots (constant outside, one knot: constant), sampled on a common grid: the run
    type's table am['logA'] plus a uniform ln A step PM_GRID_STEP over the span of
    all knots (the own knots of the modules sit at different A, so the 8-point
    run-type grid would cut their kinks)."""
    from scipy.interpolate import PchipInterpolator
    labels, f_own, A_own, x = _pm_cells(z)
    nbd, nm = f_own.shape
    se = model["c"] * np.sqrt(x)
    own = np.isfinite(se) & (se < PM_SE_MAX)
    rt = {b: own_knots[b] for b in labels if b in own_knots}
    kn = []                                  # per module: (lnA[], f[])
    for m in range(nm):
        ob = np.array([b for b in np.nonzero(own[:, m])[0] if labels[b] in rt], int)
        pts = []
        for b, lab in enumerate(labels):
            if own[b, m]:
                pts.append((A_own[b, m], f_own[b, m]))
            elif lab in rt:
                if ob.size:
                    b1 = ob[np.argmin(np.abs(ob - b))]
                    pts.append((rt[lab]["A"], f_own[b1, m] * rt[lab]["f"]
                                / rt[labels[b1]]["f"]))
                else:
                    pts.append((rt[lab]["A"], rt[lab]["f"]))
        kn.append((np.log([p[0] for p in pts]), np.array([p[1] for p in pts])))
    la_rt, f_rt = np.array(am["logA"]), np.array(am["f"])
    lo = min([la_rt[0]] + [k[0][0] for k in kn if k[0].size])
    hi = max([la_rt[-1]] + [k[0][-1] for k in kn if k[0].size])
    ug = np.linspace(lo, hi, int(np.ceil((hi - lo) / PM_GRID_STEP)) + 1)
    ug = ug[np.abs(ug[:, None] - la_rt[None, :]).min(axis=1) > 1e-3]
    grid = np.sort(np.concatenate([la_rt, ug]))
    F = np.empty((nm, grid.size))
    for m, (la, fv) in enumerate(kn):
        if la.size == 0:
            F[m] = np.interp(grid, la_rt, f_rt)
        elif la.size == 1:
            F[m] = fv[0]
        else:
            F[m] = PchipInterpolator(la, fv)(np.clip(grid, la[0], la[-1]))
    # interpolation error of the stored table at the modules' knots
    errs = np.array([np.abs(np.interp(la, grid, F[m]) - fv).max()
                     for m, (la, fv) in enumerate(kn) if la.size])
    err = errs.max() if errs.size else 0.0
    has = np.isfinite(f_own)
    stats = dict(
        labels=labels, se_max=PM_SE_MAX, se_model=model, n_modules=nm,
        n_cells_with_ratio=int(has.sum()), n_cells_own=int(own.sum()),
        frac_own_of_cells_with_ratio=float(own.sum() / max(has.sum(), 1)),
        per_band={lab: dict(n_with_ratio=int(has[b].sum()), n_own=int(own[b].sum()),
                            frac_own=float(own[b].sum() / max(has[b].sum(), 1)),
                            se_median=float(np.nanmedian(se[b])) if has[b].any() else None,
                            own_f_median=float(np.median(f_own[b][own[b]]))
                            if own[b].any() else None,
                            own_f_std=float(np.std(f_own[b][own[b]], ddof=1))
                            if own[b].sum() > 1 else None,
                            run_type_f=rt[lab]["f"] if lab in rt else None)
                  for b, lab in enumerate(labels)},
        n_grid=int(grid.size), table_interp_err_at_knots=float(err),
        table_interp_err_median=float(np.median(errs)) if errs.size else 0.0,
        n_modules_all_run_type=int((~own.any(axis=0)).sum()))
    return grid, F, np.interp(grid, la_rt, f_rt), own.astype(np.int8), se, stats


def cmd_ampcorr(a):
    """Compare the run types' band knots (f_b, b = lo / mid / top, or the E_k-E_k+1
    labels of a `build --bands` bank: any number of knots; the files' band labels
    are pooled); one global f(A) if consistent.

    Rule: for each band b present (finite positive SE) in at least two files (all
    of them if a band is in every file), over the files having it, the SE-weighted
    mean fbar_b = sum(f_i / se_i^2) / sum(1 / se_i^2), se_bar = (sum 1 /
    se_i^2)^-1/2; run type i is consistent if |f_i - fbar_b| <= 2 sqrt(se_i^2 +
    se_bar^2).  If every run type passes for every band it has (and at least two
    bands qualify), ONE global correction is written into every file: knots (A_b =
    median over the files having b, f_b = fbar_b; a band a run type lacks is still
    in its global knots), scope "global"; otherwise each file keeps (or is reset to)
    its own per-run-type knots, scope "per_run_type".  A band in fewer than two
    files is left out (warning).  One file: nothing to compare, scope
    "per_run_type".  The stored table is fa_table(knots).

    --per-module-f OUTDIR: the inputs stay as they are; every bank is copied to
    OUTDIR/<slug>/templates_<slug>.npz with ampcorr_f [module, ln A grid] (the module
    rows of the bank file), ampcorr_f_rt (the run-type f on that grid, for modules
    without a row), ampcorr_pm_own / ampcorr_pm_se [band, module] and
    ampcorr_meta.per_module = true (_pm_table).  The consistency test is run and
    reported as above; the knots a module falls back to are always its run type's."""
    files = [str(Path(f)) for f in a.templates]
    Z, own, rules = [], [], []
    for f in files:
        z = dict(np.load(f, allow_pickle=False))
        if "ampcorr_meta" not in z:
            sys.exit(f"[ERROR] {f}: no ampcorr_meta (built before f(A)?)")
        am = json.loads(str(z["ampcorr_meta"]))
        o = am.get("per_run_type")
        if not o or "knots" not in o:
            sys.exit(f"[ERROR] {f}: ampcorr_meta has no per-run-type band knots "
                     f"(built before the three-band f(A)?)")
        meta = json.loads(str(z["meta"]))
        Z.append(z)
        own.append(dict(o, run_type=meta.get("run_type", Path(f).stem)))
        rules.append(am.get("rule", AMPCORR_RULE))
    n = len(Z)
    # band labels of all files, in order of the mean knot amplitude
    labels = {}
    for o in own:
        for b, k in o["knots"].items():
            labels.setdefault(b, []).append(k["A"])
    bands = sorted(labels, key=lambda b: float(np.mean(labels[b])))
    stat, usable = {}, []
    ok_all = n > 1
    for q in bands:
        idx = [i for i, o in enumerate(own) if q in o["knots"]]
        if len(idx) < n:
            print(f"[WARN] ampcorr: band {q} missing in "
                  f"{[o['run_type'] for o in own if q not in o['knots']]}; "
                  f"weighted mean and test over {len(idx)} file(s)", flush=True)
        if n > 1 and len(idx) < 2:
            print(f"[WARN] ampcorr: band {q} in fewer than two files; band left out",
                  flush=True)
            continue
        if not idx:
            continue
        f = np.array([own[i]["knots"][q]["f"] for i in idx])
        se = np.array([own[i]["knots"][q]["se"] for i in idx])
        if n > 1 and not (np.isfinite(se).all() and (se > 0).all()):
            print(f"[WARN] ampcorr: non-finite or zero SE in f_{q}: {se}; band left out",
                  flush=True)
            continue
        if n > 1:
            w = 1.0 / se ** 2
            fbar, sebar = float((w * f).sum() / w.sum()), float(w.sum() ** -0.5)
            dev = (f - fbar) / np.sqrt(se ** 2 + sebar ** 2)
            ok = np.abs(dev) <= AMPCORR_NSIG
            ok_all = ok_all and bool(ok.all())
        else:
            fbar, sebar, dev, ok = float(f[0]), float(se[0]), np.zeros(1), np.ones(1, bool)
        usable.append(q)
        stat[q] = dict(files=idx, f=f.tolist(), se=se.tolist(), wmean=fbar,
                       wmean_se=sebar, dev_in_combined_se=dev.tolist(),
                       consistent=ok.tolist())
    if n > 1 and len(usable) < 2:
        ok_all = False
    for i, o in enumerate(own):
        print(f"[ampcorr] {o['run_type']:16s} n {o['n_modules']:4d}  " + "  ".join(
            f"f_{b} {k['f']:.4f} +- {k['se']:.4f} @ {k['A']:.0f}"
            for b, k in o["knots"].items())
              + ("  dev " + " ".join(
                  f"{b} {stat[b]['dev_in_combined_se'][stat[b]['files'].index(i)]:+.2f}"
                  for b in usable if i in stat[b]["files"]) if n > 1 else ""))
    glob_ = n > 1 and ok_all
    if a.per_module_f and glob_:
        print("[WARN] ampcorr: run types consistent, but --per-module-f falls back to "
              "the per-run-type knots (scope per_run_type)", flush=True)
        glob_ = False
    pm_model = _pm_se_model(Z) if a.per_module_f else None
    if pm_model:
        print(f"[ampcorr] per-module stat model: SE(r_b / r_mod) = {pm_model['c']:.4f} "
              f"sqrt(1/n_b + 1/n_pool) (n-independent scatter s0 {pm_model['s0']:.4f}, "
              f"{pm_model['n_bins']} bins); own value kept below SE {PM_SE_MAX:g}",
              flush=True)
    if n > 1:
        print("[ampcorr] weighted mean " + ", ".join(
            f"f_{b} {stat[b]['wmean']:.4f} +- {stat[b]['wmean_se']:.4f}" for b in usable)
              + f"; all within {AMPCORR_NSIG:g} combined SE: {ok_all} -> scope "
              f"{'global' if glob_ else 'per_run_type'}", flush=True)
    else:
        print("[ampcorr] one file: nothing to compare; scope per_run_type", flush=True)
    if glob_:
        if len(usable) < len(bands):
            print(f"[WARN] ampcorr: global f(A) from {len(usable)} bands only "
                  f"({', '.join(usable)})", flush=True)
        gl_knots = {b: dict(A=float(np.median([own[i]["knots"][b]["A"]
                                               for i in stat[b]["files"]])),
                            f=stat[b]["wmean"], se=stat[b]["wmean_se"],
                            n_modules=int(sum(own[i]["knots"][b]["n_modules"]
                                              for i in stat[b]["files"])))
                    for b in usable}
        gl = dict(knots=gl_knots, n_modules=int(sum(o["n_modules"] for o in own)))
    summary = dict(
        rule=rules[0],
        consistency_rule=f"per band q in ({', '.join(bands)}) present in at least two run "
                         f"types: weighted mean over the run types having it of f_q "
                         f"with weights 1/se^2, se_bar = (sum 1/se^2)^-1/2; every run "
                         f"type must satisfy |f_i - fbar| <= {AMPCORR_NSIG:g} "
                         f"sqrt(se_i^2 + se_bar^2) for every band it has; then ONE "
                         f"global f (fbar_q; A_q = medians over the run types having q) "
                         f"in every file, else each keeps its own",
        scope="global" if glob_ else "per_run_type",
        note=("single file: consistency not testable" if n == 1 else
              ("all run types consistent" if glob_ else "at least one run type deviates")),
        run_types=[o["run_type"] for o in own], per_run_type=own, stat=stat,
        global_values=gl if glob_ else None, files=[], created=datetime.now().isoformat(
            timespec="seconds"), argv=sys.argv)
    for f, z, o, rule in zip(files, Z, own, rules):
        v = gl if glob_ else o
        tab = fa_table(v["knots"])
        if tab is None:
            print(f"[WARN] ampcorr: {f}: fewer than two usable knots; file left as is",
                  flush=True)
            continue
        am = dict(rule=rule, scope="global" if glob_ else "per_run_type",
                  run_type=o["run_type"], n_modules=v["n_modules"], knots=v["knots"],
                  f=tab[1], logA=tab[0], n_table=len(tab[0]),
                  per_run_type={k: o[k] for k in ("knots", "n_modules")},
                  se_rule="1.2533 * std(ddof=1) / sqrt(n) of the per-module ratios "
                          "(global: SE of the weighted mean)",
                  note=summary["note"],
                  consistency=dict(rule=summary["consistency_rule"],
                                   run_types=summary["run_types"], stat=stat))
        z["ampcorr_logA"] = np.array(am["logA"])
        z["ampcorr_f"] = np.array(am["f"])
        dest = f
        if a.per_module_f:
            slug = Path(f).stem.removeprefix("templates_")
            dest = str(Path(a.per_module_f) / slug / f"templates_{slug}.npz")
            if Path(dest).resolve() == Path(f).resolve():
                sys.exit(f"[ERROR] --per-module-f: {f} would overwrite its own input")
            grid, F, f_rt, own_c, se_c, pst = _pm_table(z, am, v["knots"], pm_model)
            Path(dest).parent.mkdir(parents=True, exist_ok=True)
            am["per_module"] = True
            am["per_module_rule"] = (
                "ampcorr_f[nm, ngrid] rows = the module rows of this file; module "
                "cell (band b): own r_b / r_mod at A = the module's median amplitude "
                "in b if SE = c sqrt(1/n_b + 1/n_pool) < se_max (c from the pooled "
                "regression of the robust step scatter on 1/n_b + 1/n_ref, s0 the "
                "n-independent part; n capped at n_max), else chained from the module's "
                "nearest own band b' by the run-type step: f_own(b') f_rt(b) / "
                "f_rt(b') at the run-type knot A (no own cell: the run-type row); PCHIP in ln A per "
                "module, constant outside its knots; ampcorr_f_rt = the run-type "
                "f(A) (the 1-D table of ampcorr_meta f / logA) on the same grid; "
                "ampcorr_pm_own / ampcorr_pm_se [band, module]")
            am["per_module_stats"] = pst
            z["ampcorr_logA"] = grid
            z["ampcorr_f"] = F
            z["ampcorr_f_rt"] = f_rt
            z["ampcorr_pm_own"] = own_c
            z["ampcorr_pm_se"] = se_c
            print(f"[ampcorr] per-module {slug}: {pst['n_cells_own']} of "
                  f"{pst['n_cells_with_ratio']} cells with a ratio keep their own value "
                  f"({100 * pst['frac_own_of_cells_with_ratio']:.1f} %); " + ", ".join(
                      f"{b} {100 * d['frac_own']:.0f} %" for b, d in pst["per_band"].items())
                  + f"; table {F.shape[1]} points, knot interpolation error "
                  f"{pst['table_interp_err_at_knots']:.1e}", flush=True)
        z["ampcorr_meta"] = np.array(json.dumps(am))
        meta = json.loads(str(z["meta"]))
        meta["ampcorr"] = am
        z["meta"] = np.array(json.dumps(meta))
        _rewrite_npz(dest, z)
        summary["files"].append(dict(path=os.path.abspath(dest), sha1=_sha1(dest),
                                     scope=am["scope"], knots=am["knots"],
                                     logA=am["logA"], f=am["f"]))
        if a.per_module_f:
            summary["files"][-1].update(per_module=True, source=os.path.abspath(f),
                                        per_module_stats=pst)
        print(f"[ampcorr] rewrote {dest}: scope {am['scope']}, knots " + ", ".join(
            f"{b} f {k['f']:.4f} @ {k['A']:.0f}" for b, k in v["knots"].items())
              + f"; table {len(tab[0])} points", flush=True)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w") as fh:
        json.dump(_jsonable(summary), fh, indent=1)
    print(f"[ampcorr] summary {a.out}", flush=True)
    return 0


# =================================================================== check ==
# Plan sections 6.4-6.6.  Gates G1-4 .. G1-9 run on the held-out 20 % of the
# candidate files; G1-1 .. G1-3 (cross-check against the PASS 1 macro on a
# `_raw.root` split) run when --pass1-macro-out / --raw-root / --macro-tsv are
# given.  Everything is written to --out: gates.json, gates.tsv, per-gate detail
# TSVs and one PNG per gate.

CHK_AMP_BINS = ((100.0, 300.0), (300.0, 1000.0), (1000.0, 2500.0))
CHK_MF_FRAC = 0.03                      # the mf model-error term (r_model)
CHK_RISE_TAU = (0.0, 6.0)               # G1-5 rise window (ns after onset)
CHK_MIN_N = 30                          # G1-5 / G1-4: smaller groups are "insufficient"
CHK_MIN_MOD_PULSES = 10                 # G1-7 per-module RMS
CHK_MAXDIFF_WIN = (0.0, 20.0)           # G1-9 window (ns)
CHK_RNG_SEED = 20261024                 # held-out subsample / noise-floor halves
CHK_FIELDS = ("event_number", "midx", "samples", "ped", "rms", "amp", "tier",
              "mf_t0", "wa_integral")
CHK_COVERAGE_RULE = {"2p24GeV_H2": 1100, "3p49GeV_H2": 1100, "2p24GeV_C": 1000}
RAW_STRIDE = 200                        # hycal.samples is uint16[nch][200]
RAW_PED_N = 20                          # PASS 1 pedestal: mean of raw samples 0-19


def _jsonable(o):
    """numpy -> plain python; NaN / inf -> None."""
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return _jsonable(o.tolist())
    if isinstance(o, (np.integer, np.bool_)):
        return o.item()
    if isinstance(o, (float, np.floating)):
        v = float(o)
        return v if math.isfinite(v) else None
    return o


def _fmt(v):
    if isinstance(v, (float, np.floating)):
        return "nan" if not math.isfinite(v) else f"{v:.6g}"
    return str(v)


def _write_tab(path, header, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(_fmt(v) for v in r) + "\n")


def _bin_label(lo, hi):
    return f"{lo:g}-{hi:g}"


def _med(x):
    x = np.asarray(x, dtype=np.float64)
    x = x[np.isfinite(x)]
    return float(np.median(x)) if x.size else float("nan")


class Gates:
    """Collects the gate records (id, rule, result, numbers) and writes them."""

    def __init__(self):
        self.rec = {}

    def add(self, gid, title, rule, result, numbers, summary):
        assert result in ("PASS", "FAIL", "REPORT", "N/A")
        self.rec[gid] = dict(id=gid, title=title, rule=rule, result=result,
                             summary=summary, numbers=numbers)
        print(f"[{gid}] {result:6s} {title}: {summary}", flush=True)

    def write(self, out, header):
        order = sorted(self.rec, key=lambda g: [int(x) if x.isdigit() else 0
                                                for x in re.split(r"[-.]", g[1:])])
        doc = dict(header, gates=[self.rec[g] for g in order])
        with open(out / "gates.json", "w") as fh:
            json.dump(_jsonable(doc), fh, indent=1)
        _write_tab(out / "gates.tsv", ["gate", "result", "title", "rule", "summary"],
                   [[g, self.rec[g]["result"], self.rec[g]["title"],
                     self.rec[g]["rule"], self.rec[g]["summary"]] for g in order])


# ------------------------------------------------------------ loading ------
def _load_check_inputs(files, seed, frac, noise_mods, noise_band, cap, rng,
                       hold_keys=None):
    """One pass over the candidate files.  Returns (held, ntot, nheld, pool,
    metas): `held` = dict of the held-out candidates (CHK_FIELDS + run, split),
    `pool` = non-held-out tier>=1 candidates of the modules in `noise_mods` in
    amplitude band `noise_band` (for the G1-9 noise floor).  hold_keys = set of
    (run, split, event_number, midx) tuples (templates npz `holdout_keys` of a
    --pool-order energy build) replaces holdout_mask."""
    held, pool, metas = [], [], []
    ntot = 0
    nm_set = np.array(sorted(noise_mods), dtype=np.int64)
    for f in files:
        D = load_candidates([f], CHK_FIELDS)
        metas.append(D["metas"][0])
        if hold_keys is None:
            h = holdout_mask(D["run"], D["split"], D["event_number"], D["midx"],
                             seed=seed, frac=frac)
        else:
            r0, s0 = int(D["run"][0]) if D["run"].size else 0, \
                int(D["split"][0]) if D["split"].size else 0
            h = np.array([(r0, s0, e, m) in hold_keys for e, m in
                          zip(D["event_number"].tolist(), D["midx"].tolist())],
                         dtype=bool)
        ntot += h.size
        held.append({k: D[k][h] for k in CHK_FIELDS + ("run", "split")})
        if nm_set.size:
            q = (~h) & (D["tier"] >= 1) & np.isin(D["midx"], nm_set) \
                & (D["amp"] >= noise_band[0]) & (D["amp"] <= noise_band[1])
            pool.append({k: D[k][q] for k in ("midx", "samples", "ped", "rms", "mf_t0")})
        del D
    H = {k: np.concatenate([x[k] for x in held]) for k in held[0]}
    nheld = int(H["midx"].size)
    P = ({k: np.concatenate([x[k] for x in pool]) for k in pool[0]}
         if pool else None)
    # per-module subsample of the held-out set (deterministic)
    keep = np.ones(nheld, dtype=bool)
    if cap and cap > 0:
        order = np.argsort(H["midx"], kind="stable")
        ms = H["midx"][order]
        u, st, cnt = np.unique(ms, return_index=True, return_counts=True)
        for s0, c in zip(st, cnt):
            if c > cap:
                drop = rng.choice(c, c - cap, replace=False)
                keep[order[s0 + drop]] = False
    H = {k: v[keep] for k, v in H.items()}
    return H, ntot, nheld, P, metas


# ----------------------------------------------------------- fit pass -------
def _chi2_mf(Y, sig, t0, A, shape, win=WIN):
    """chi2/dof with the mf variance sigma^2 + (0.03 * model)^2."""
    cols = _align_columns(t0)
    Ys = np.take_along_axis(Y, cols, axis=1)
    tau = TS[cols] - t0[:, None]
    msk = (tau >= win[0]) & (tau <= win[1])
    mdl = A[:, None] * shape.eval(tau)
    var = sig[:, None] ** 2 + (CHK_MF_FRAC * mdl) ** 2
    chi2 = (((Ys - mdl) ** 2 / var) * msk).sum(1)
    return chi2 / np.maximum(msk.sum(1) - 2, 1)


def _chi2_full(Y, sig, t0, A, shape, win=WIN):
    """chi2/dof with the variance sigma^2 + (A*SIGMA_T_NS*f'(tau))^2 +
    (A*SHAPE_SIGMA*f(tau))^2 over the `win` samples (dof = n_window - 2); the
    held-out measure of --lambda-chi2-sigma full."""
    cols = _align_columns(t0)
    Ys = np.take_along_axis(Y, cols, axis=1)
    tau = TS[cols] - t0[:, None]
    msk = (tau >= win[0]) & (tau <= win[1])
    f = shape.eval(tau)
    slope = (shape.eval(tau + 0.05) - shape.eval(tau - 0.05)) / 0.1
    var = sig[:, None] ** 2 + (A[:, None] * SIGMA_T_NS * slope) ** 2 \
        + (A[:, None] * SHAPE_SIGMA * f) ** 2
    chi2 = (((Ys - A[:, None] * f) ** 2 / var) * msk).sum(1)
    return chi2 / np.maximum(msk.sum(1) - 2, 1)


def _rise_resid(Y, t0, A, shape):
    """Mean over the samples with 0 < tau < 6 ns of (y / A - shape(tau)), in
    units of the peak (NaN if the pulse has no such sample)."""
    cols = _align_columns(t0)
    Ys = np.take_along_axis(Y, cols, axis=1)
    tau = TS[cols] - t0[:, None]
    msk = (tau > CHK_RISE_TAU[0]) & (tau < CHK_RISE_TAU[1])
    r = (Ys / np.where(A > 0, A, np.nan)[:, None] - shape.eval(tau)) * msk
    n = msk.sum(1)
    return np.where(n > 0, r.sum(1) / np.maximum(n, 1), np.nan), n


def fit_heldout(H, z, bank_s, tau_rows, tau_agg, names):
    """Fit (A, t0) of every held-out pulse with its module's sampled template and
    its h50 two-tau.  Returns per-pulse arrays (aligned with H) in a dict with
    sub-dicts 's' and 't'."""
    n = H["midx"].size
    row_of = {int(m): r for r, m in enumerate(z["module_index"])}
    out = {k: {f: np.full(n, np.nan) for f in
               ("t0", "A", "chi2_rms", "chi2_mf", "rise", "rise_n", "iint")}
           for k in ("s", "t")}
    fb = np.full(n, -1, dtype=np.int8)
    shift = z["onset_shift_ns"].astype(np.float64)    # total: final axis vs the mf seed axis
    # a fallback-2 row holds the run-type average, whose total shift is relative to
    # the modules' own final axes: the module's own offset (not stored) is replaced
    # by the median of the fallback-0/1 modules and the search window is widened
    own = z["fallback"] <= 1
    d_typ = float(np.median(shift[own])) if own.any() else 0.0
    bound = {"s": np.zeros(n, dtype=bool), "t": np.zeros(n, dtype=bool)}
    order = np.argsort(H["midx"], kind="stable")
    ms = H["midx"][order]
    u, st, cnt = np.unique(ms, return_index=True, return_counts=True)
    for m, s0, c in zip(u, st, cnt):
        r = row_of.get(int(m))
        if r is None:
            continue
        idx = order[s0:s0 + c]
        fb[idx] = int(z["fallback"][r])
        Y = H["samples"][idx].astype(np.float64) - H["ped"][idx][:, None]
        sig = np.maximum(H["rms"][idx], SIG_MIN)
        seed = H["mf_t0"][idx]
        hw_s = 3.0 if own[r] else 6.0
        sh_s = shift[r] if own[r] else d_typ + shift[r]
        for key, shp, sd, hw in (
                ("s", LutShape(bank_s.lut[m].astype(np.float64)), seed + sh_s, hw_s),
                ("t", TwoTauShape(*tau_rows.get(names[m], tau_agg)), seed, 3.0)):
            al = align_pulses(Y, sig, sd, shp, hw)
            ok = al["ok"]
            bound[key][idx] = np.abs(al["t0"] - sd) > hw - 0.05
            o = out[key]
            g = lambda v: np.where(ok, v, np.nan)                # noqa: E731
            o["t0"][idx], o["A"][idx] = g(al["t0"]), g(al["A"])
            o["chi2_rms"][idx] = g(al["chi2dof"])
            o["chi2_mf"][idx] = g(_chi2_mf(Y, sig, al["t0"], al["A"], shp))
            rr, rn = _rise_resid(Y, al["t0"], al["A"], shp)
            o["rise"][idx], o["rise_n"][idx] = g(rr), rn
    return out, fb, bound


def _iint_at(bank, midx, t0, A, with_f=True):
    """Q_template = mf_pulse.pulse_energy_units(bank, m, t0, A, kappa=1) = A * I_int
    (* f(A) for a bank carrying ampcorr; with_f=False: A * bank.iint(m, t0) only);
    NaN where the fit failed."""
    import mf_pulse
    q = np.full(midx.size, np.nan)
    ok = np.isfinite(t0) & np.isfinite(A)
    if ok.any():
        if with_f:
            q[ok] = mf_pulse.pulse_energy_units(bank, midx[ok], t0[ok], A[ok], kappa=1.0)
        else:
            q[ok] = A[ok] * bank.iint(midx[ok], t0[ok])
    return q


def _bin_stats(r, gm, amp, bins=CHK_AMP_BINS):
    """Median and SE (1.2533 std / sqrt(n)) of the ratio r in the amplitude bins
    (default the three CHK_AMP_BINS; the last bin is closed)."""
    out = []
    for lo, hi in bins:
        m = gm & (amp >= lo) & (amp < hi if (lo, hi) != bins[-1] else amp <= hi) \
            & np.isfinite(r)
        n = int(m.sum())
        out.append(dict(bin=_bin_label(lo, hi), n=n,
                        median=_med(r[m]) if n else float("nan"),
                        se=float(1.2533 * np.std(r[m], ddof=1) / np.sqrt(n)) if n > 1
                        else float("nan")))
    return out


def _flatness(bins, nsig=2.0):
    """Every pairwise difference of the bin medians (bins with n >= CHK_MIN_N)
    within nsig x its combined SE."""
    have = [d for d in bins if d["n"] >= CHK_MIN_N and math.isfinite(d["se"])]
    pairs = []
    for i in range(len(have)):
        for j in range(i + 1, len(have)):
            a_, b_ = have[i], have[j]
            diff = b_["median"] - a_["median"]
            se = math.hypot(a_["se"], b_["se"])
            pairs.append(dict(a=a_["bin"], b=b_["bin"], diff_pct=100 * diff,
                              se_pct=100 * se, n_sigma=diff / se if se > 0 else float("inf"),
                              ok=bool(abs(diff) <= nsig * se)))
    return dict(rule=f"every pairwise difference of bin medians <= {nsig:g} x combined SE",
                n_bins=len(have), pairs=pairs,
                flat=bool(len(have) >= 2 and all(p_["ok"] for p_ in pairs)))


# --------------------------------------------------------------- gates ------
def _groups(fb):
    return {"fb01": (fb >= 0) & (fb <= 1), "fb2": fb == 2, "all": fb >= 0}


def gate_g14(G, H, F, fb, outd, png):
    amp = H["amp"]
    res, rows = {}, []
    for gname, gm in _groups(fb).items():
        gb = []
        for lo, hi in CHK_AMP_BINS:
            m = gm & (amp >= lo) & (amp < hi if hi < 2500 else amp <= hi)
            ms = m & np.isfinite(F["s"]["chi2_rms"]) & np.isfinite(F["t"]["chi2_rms"])
            d = dict(bin=_bin_label(lo, hi), n=int(ms.sum()))
            if ms.sum():
                for k in ("rms", "mf"):
                    a, b = _med(F["s"]["chi2_" + k][ms]), _med(F["t"]["chi2_" + k][ms])
                    d["sampled_" + k], d["two_tau_" + k] = a, b
                    d["ratio_" + k] = a / b if b > 0 else float("nan")
            gb.append(d)
            rows.append([gname, d["bin"], d["n"]] + [d.get(k, float("nan")) for k in
                         ("sampled_rms", "two_tau_rms", "ratio_rms",
                          "sampled_mf", "two_tau_mf", "ratio_mf")])
        res[gname] = gb
    ev = [d for d in res["fb01"] if d["n"] >= CHK_MIN_N]
    passed = all(d["ratio_rms"] <= 0.5 for d in ev)
    for d in res["fb01"]:
        d["status"] = ("insufficient" if d["n"] < CHK_MIN_N else
                       ("ok" if d["ratio_rms"] <= 0.5 else "FAIL"))
    result = "N/A" if not ev else ("PASS" if passed else "FAIL")
    summ = "; ".join(f"{d['bin']}: {d['sampled_rms']:.1f}/{d['two_tau_rms']:.1f}="
                     f"{d['ratio_rms']:.2f} (n={d['n']})" for d in res["fb01"]
                     if d["n"]) or "no fallback-0/1 held-out pulses"
    G.add("G1-4", "held-out chi2/dof",
          "median chi2/dof (sigma = rms) of the sampled template <= 0.5 x the two-tau "
          f"h50 value in every amplitude bin (fallback-0/1 modules; bins with n < "
          f"{CHK_MIN_N} are insufficient); fallback-2 and all reported",
          result, res, "fb01 sampled/two-tau (ratio): " + summ)
    _write_tab(outd / "g14_chi2_bins.tsv",
               ["group", "amp_bin", "n", "sampled_chi2_rms", "two_tau_chi2_rms",
                "ratio_rms", "sampled_chi2_mf", "two_tau_chi2_mf", "ratio_mf"], rows)
    if png:
        fig, ax = png.subplots(1, 2, figsize=(10, 3.8))
        for a_, k, ttl in ((ax[0], "rms", "sigma = rms"), (ax[1], "mf", "sigma^2 + (0.03 model)^2")):
            x = np.arange(len(CHK_AMP_BINS))
            for j, (gname, c) in enumerate((("fb01", "C0"), ("fb2", "C1"))):
                v = [d.get("sampled_" + k, np.nan) for d in res[gname]]
                w = [d.get("two_tau_" + k, np.nan) for d in res[gname]]
                a_.bar(x + 0.36 * j - 0.27, v, 0.17, color=c, label=f"sampled {gname}")
                a_.bar(x + 0.36 * j - 0.09, w, 0.17, color=c, alpha=0.4,
                       label=f"two-tau {gname}")
            a_.set_xticks(x)
            a_.set_xticklabels([_bin_label(*b) for b in CHK_AMP_BINS])
            a_.set_yscale("log")
            a_.set_xlabel("candidate amplitude (ADC)")
            a_.set_ylabel("median chi2/dof")
            a_.set_title(ttl)
        ax[0].legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(outd / "g14_chi2.png", dpi=110)
        png.close(fig)


def gate_g15(G, H, F, fb, outd, png):
    amp = H["amp"]
    good = amp >= CHK_AMP_BINS[0][0]
    res, rows = {}, []
    for gname, gm in _groups(fb).items():
        gq = {}
        for key, lab in (("s", "sampled"), ("t", "two_tau")):
            t0, rs = F[key]["t0"], F[key]["rise"]
            m = gm & good & np.isfinite(rs)
            q = np.floor(np.mod(t0, DT_NS)).astype(np.int64)
            lst = []
            for k in range(4):
                mk = m & (q == k)
                n = int(mk.sum())
                d = dict(quartile=f"[{k},{k + 1}) ns", n=n)
                if n:
                    v = rs[mk] * 100.0
                    d["mean_pct"] = float(v.mean())
                    d["se_pct"] = float(v.std(ddof=1) / math.sqrt(n)) if n > 1 else float("nan")
                    d["status"] = ("insufficient" if n < CHK_MIN_N else
                                   ("ok" if abs(d["mean_pct"]) < 0.5 else "FAIL"))
                else:
                    d["mean_pct"] = d["se_pct"] = float("nan")
                    d["status"] = "insufficient"
                lst.append(d)
                rows.append([gname, lab, k, n, d["mean_pct"], d["se_pct"], d["status"]])
            gq[lab] = lst
        res[gname] = gq
    ev = [d for d in res["fb01"]["sampled"] if d["status"] != "insufficient"]
    result = "N/A" if not ev else ("PASS" if all(d["status"] == "ok" for d in ev) else "FAIL")
    summ = "; ".join(f"{d['quartile']}: {d['mean_pct']:+.2f}+-{d['se_pct']:.2f}% n={d['n']}"
                     f"{' (insufficient)' if d['status'] == 'insufficient' else ''}"
                     for d in res["fb01"]["sampled"])
    G.add("G1-5", "sub-sample rise",
          "mean residual (y/A - template) on 0 < tau < 6 ns per pulse, grouped by "
          "the fitted phase (t0 mod 4 ns) in four 1 ns quartiles: |mean| < 0.5 % of "
          f"peak in every quartile with N >= {CHK_MIN_N} pulses (fallback-0/1 "
          "modules, amp >= 100); N counts pulses; fallback-2, all and the two-tau "
          "(informational) are in the json", result, res,
          "fb01 sampled: " + summ)
    _write_tab(outd / "g15_quartiles.tsv",
               ["group", "template", "phase_quartile", "n_pulses", "mean_resid_pct",
                "se_pct", "status"], rows)
    if png:
        fig, ax = png.subplots(1, 2, figsize=(10, 3.8))
        for j, (gname, c) in enumerate((("fb01", "C0"), ("fb2", "C1"))):
            d = res[gname]["sampled"]
            ax[0].errorbar(np.arange(4) + 0.1 * j, [x["mean_pct"] for x in d],
                           [x["se_pct"] for x in d], fmt="o", color=c,
                           label=f"sampled {gname}")
        d = res["fb01"]["two_tau"]
        ax[0].errorbar(np.arange(4) + 0.3, [x["mean_pct"] for x in d],
                       [x["se_pct"] for x in d], fmt="s", color="gray",
                       label="two-tau fb01 (own t0)")
        ax[0].axhspan(-0.5, 0.5, color="g", alpha=0.12)
        ax[0].set_xticks(range(4))
        ax[0].set_xticklabels(["[0,1)", "[1,2)", "[2,3)", "[3,4)"])
        ax[0].set_xlabel("fitted phase t0 mod 4 ns")
        ax[0].set_ylabel("mean rise residual (% of peak)")
        ax[0].legend(fontsize=7)
        t0 = F["s"]["t0"][good & np.isfinite(F["s"]["t0"])]
        ax[1].hist(np.mod(t0, DT_NS), bins=np.arange(0, 4.01, 0.25), color="C0")
        ax[1].set_xlabel("t0 mod 4 ns (sampled fit)")
        ax[1].set_ylabel("held-out pulses")
        fig.tight_layout()
        fig.savefig(outd / "g15_rise.png", dpi=110)
        png.close(fig)


G16_NPH = 16                           # phases phi = 0, 0.25, ..., 3.75 ns


def _g16_run(p, mi, names, bank, amp, nph=G16_NPH):
    """Noise-free integer waveforms ped 200 + amp * shape(4 i - t0), t0 = 140 + phi,
    phi = 4 k / nph ns, through wave_ana.analyze.  Per module: the phase mean of
    pk.integral / amp against the phase mean of bank.iint(m, t0) (a constant for
    the continuous sampled bank) and the per-phase deviations."""
    ana = p.wave_ana.analyze
    nm = mi.size
    Q = np.full((nm, nph), np.nan)
    R = np.zeros((nm, nph))
    nbad = 0
    for j in range(nph):
        phi = DT_NS * j / nph
        t0 = 140.0 + phi
        S = bank.shape(mi[:, None], DT_NS * np.arange(NS_EXPECT)[None, :] - t0,
                       dtype=np.float64)
        Yw = np.rint(200.0 + amp * S).astype(np.uint16)
        R[:, j] = bank.iint(mi, np.full(nm, t0))
        for k in range(nm):
            _, _, pks = ana(Yw[k])
            if len(pks) != 1:
                nbad += 1
            else:
                Q[k, j] = pks[0].integral / amp
    dph = Q / R - 1.0                                   # per phase (NaN: no single peak)
    qm, rm = Q.mean(axis=1), R.mean(axis=1)             # NaN if any phase failed
    dmean = qm / rm - 1.0

    def stats(d):
        a_ = np.where(np.isfinite(d), np.abs(d), 9.9)
        f = a_[a_ < 9.9]
        return dict(max_abs_dev=float(f.max()) if f.size else float("nan"),
                    p99_abs_dev=float(np.percentile(f, 99)) if f.size else float("nan"),
                    median_abs_dev=float(np.median(f)) if f.size else float("nan"))

    pm_over = np.where(np.isfinite(dmean), np.abs(dmean), 9.9) > 0.005
    pp_mod = np.where(np.isfinite(dph), np.abs(dph), 9.9).max(axis=1)
    ref_const = bool(np.all(R == R[:, :1]))
    worst_i = np.argsort(-np.where(np.isfinite(dmean), np.abs(dmean), 9.9))[:10]
    st = dict(amplitude=amp, n_modules=int(nm), n_phases=nph, n_waveforms=int(Q.size),
              n_no_single_peak=int(nbad), bank_iint_phase_independent=ref_const,
              phase_mean=dict(stats(dmean), n_modules_over_0p5pct=int(pm_over.sum())),
              per_phase=dict(stats(dph.reshape(-1)),
                             n_modules_over_0p5pct=int((pp_mod > 0.005).sum())),
              worst_phase_mean=[dict(module=names[i], dev=float(dmean[i]),
                                     q_mean=float(qm[i]), bank_iint_mean=float(rm[i]))
                                for i in worst_i])
    recs = [(names[k], DT_NS * j / nph, float(Q[k, j]), float(R[k, j]), float(dph[k, j]))
            for k in range(nm) for j in range(nph)]
    return st, recs, dmean, dph


def gate_g16(G, p, mods_pwo4, bank_s, bank_t, outd, png):
    """Phase-mean C++ integral vs the bank's I_int on noise-free integer waveforms
    (16 phases; amplitude 1000 gates the sampled bank, amplitude 3000 and the
    per-phase deviations are informational: the continuous I_int does not step
    with the sub-sample phase, the C++ integral does)."""
    mi = np.array([i for i, f in enumerate(mods_pwo4) if f], dtype=np.int64)
    names = [p.hycal.module(int(i)).name for i in mi]
    res, rows, dm, dp_ = {}, [], {}, {}
    for lab, bank in (("sampled", bank_s), ("two_tau", bank_t)):
        st, recs, dmean, dph = _g16_run(p, mi, names, bank, 1000.0)
        res[lab] = st
        rows += [[lab, r[0], r[1], r[2], r[3], r[4]] for r in recs]
        dm[lab], dp_[lab] = dmean, dph
    res["sampled"]["amplitude_3000"] = _g16_run(p, mi, names, bank_s, 3000.0)[0]
    s_ = res["sampled"]
    bad = s_["phase_mean"]["n_modules_over_0p5pct"] + s_["n_no_single_peak"]
    result = "PASS" if bad == 0 else "FAIL"
    summ = "; ".join(
        f"{k}: phase-mean max {res[k]['phase_mean']['max_abs_dev'] * 100:.3f}% p99 "
        f"{res[k]['phase_mean']['p99_abs_dev'] * 100:.3f}% over 0.5%: "
        f"{res[k]['phase_mean']['n_modules_over_0p5pct']}; per-phase max "
        f"{res[k]['per_phase']['max_abs_dev'] * 100:.3f}% p99 "
        f"{res[k]['per_phase']['p99_abs_dev'] * 100:.3f}% (modules with a phase over "
        f"0.5%: {res[k]['per_phase']['n_modules_over_0p5pct']}); no-peak "
        f"{res[k]['n_no_single_peak']}" for k in ("sampled", "two_tau"))
    a3 = s_["amplitude_3000"]
    summ += (f"; sampled at amp 3000: phase-mean max {a3['phase_mean']['max_abs_dev'] * 100:.3f}%"
             f", over 0.5%: {a3['phase_mean']['n_modules_over_0p5pct']}")
    res["note"] = ("waveforms are rounded to integers (ped 200 + 1000 shape(4 i - t0), "
                   "t0 = 140 + phi, phi = 0, 0.25, ..., 3.75 ns); the binding casts "
                   "float input to integers by truncation, so an unrounded variant "
                   "cannot be run.  Gate: sampled bank, amplitude 1000, mean over the "
                   "16 phases of pk.integral/1000 vs bank.iint (continuous I_int) per "
                   "module; the per-phase deviations (the step the continuous mode "
                   "removes) and the two-tau bank (iint = phase table, compared at "
                   "each t0 and as the phase mean) are informational")
    G.add("G1-6", "I_int vs the C++ integral",
          "sampled bank: phase mean (16 phases, amplitude 1000, ped 200, integer "
          "waveforms) of pk.integral/1000 vs bank.iint within 0.5 % for every PbWO4 "
          "module; per-phase deviations and the two-tau bank informational",
          result, res, summ)
    _write_tab(outd / "g16_modules.tsv",
               ["bank", "module", "phi_ns", "pk_integral_over_1000", "bank_iint", "dev"], rows)
    if png:
        fig, ax = png.subplots(figsize=(6, 3.8))
        for k, c in (("sampled", "C0"), ("two_tau", "C3")):
            for d, ls in ((dm[k], "-"), (dp_[k].reshape(-1), ":")):
                v = np.abs(d[np.isfinite(d)]) * 100
                ax.hist(np.clip(v, 1e-4, None), bins=np.logspace(-4, 1, 60),
                        histtype="step", color=c, ls=ls,
                        label=f"{k} {'phase mean' if ls == '-' else 'per phase'}")
        ax.axvline(0.5, color="k", ls="--")
        ax.set_xscale("log")
        ax.set_xlabel("|pk.integral / bank.iint - 1| (%), amplitude 1000")
        ax.set_ylabel("module (x phase)")
        ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(outd / "g16_iint.png", dpi=110)
        png.close(fig)


def gate_g17(G, H, F, fb, bank_s, bank_t, names, outd, png, amp_bins=CHK_AMP_BINS):
    """amp_bins: the bins of the flatness gate and the bin tables (default
    CHK_AMP_BINS); the per-module RMS (A > 100 / A > 300), the pre-f(A) absolute
    conditions, kappa and the two-tau check always use CHK_AMP_BINS."""
    amp, midx = H["amp"], H["midx"]
    qraw = H["wa_integral"]
    ratio = {}
    for key, bank in (("s", bank_s), ("t", bank_t)):
        q = _iint_at(bank, midx, F[key]["t0"], F[key]["A"])      # includes f(A) (sampled)
        ratio[key] = q / qraw
    ratio_nof = _iint_at(bank_s, midx, F["s"]["t0"], F["s"]["A"], with_f=False) / qraw
    res, rows = {}, []
    for gname, gm in _groups(fb).items():
        gb = {}
        for key, lab in (("s", "sampled"), ("t", "two_tau")):
            lst = _bin_stats(ratio[key], gm, amp, amp_bins)
            for d in lst:
                rows.append([gname, lab, d["bin"], d["n"], d["median"], d["se"]])
            # per-module statistics of the ratio.  within_rms: RMS of the ratio
            # about the module mean, median over modules (pulse-level scatter;
            # amp100 is only the literal within-module scatter number);
            # spread: RMS across modules of the per-module MEDIAN ratio (the
            # note section 6.1 "per-module RMS"), A > 300, >= 10 pulses.
            # stat_expected: sampling error of those medians,
            # sqrt(mean(within_i^2 / n_i)) * sqrt(pi / 2); spread_corrected removes
            # it in quadrature
            rmsv = {}
            for tag, ampmin, minn in (("amp100", CHK_AMP_BINS[0][0], CHK_MIN_MOD_PULSES),
                                      ("amp300", CHK_AMP_BINS[1][0], CHK_MIN_MOD_PULSES)):
                mm = gm & (amp >= ampmin) & np.isfinite(ratio[key])
                within, meds, nn = [], [], []
                if mm.any():
                    o = np.argsort(midx[mm], kind="stable")
                    mids, r_ = midx[mm][o], ratio[key][mm][o]
                    uu, ss, cc = np.unique(mids, return_index=True, return_counts=True)
                    for s0, c in zip(ss, cc):
                        if c >= minn:
                            within.append(float(np.std(r_[s0:s0 + c])))
                            meds.append(float(np.median(r_[s0:s0 + c])))
                            nn.append(int(c))
                meds = np.array(meds)
                spread = (float(np.sqrt(np.mean((meds - meds.mean()) ** 2)))
                          if meds.size else float("nan"))
                stat = (float(np.sqrt(np.mean(np.array(within) ** 2 / np.array(nn))) * 1.2533)
                        if meds.size else float("nan"))
                rmsv[tag] = dict(
                    n_modules=len(within),
                    within_module_rms_median=float(np.median(within)) if within else float("nan"),
                    module_median_spread_rms=spread,
                    module_median_stat_expected=stat,
                    module_median_spread_corrected=(float(np.sqrt(max(spread ** 2 - stat ** 2, 0.0)))
                                                    if meds.size else float("nan")),
                    median_of_module_medians=float(np.median(meds)) if meds.size else float("nan"))
                if tag == "amp100":
                    gb.setdefault("_modrms_" + lab, within)
                if tag == "amp300":
                    gb.setdefault("_modmed_" + lab, meds.tolist())
            gb[lab] = dict(bins=lst, per_module_rms=rmsv)
        gb["sampled_no_fA"] = dict(bins=_bin_stats(ratio_nof, gm, amp, amp_bins))
        for d in gb["sampled_no_fA"]["bins"]:
            rows.append([gname, "sampled_no_fA", d["bin"], d["n"], d["median"], d["se"]])
        gb["flatness_with_fA"] = _flatness(gb["sampled"]["bins"])
        gb["flatness_without_fA"] = _flatness(gb["sampled_no_fA"]["bins"])
        res[gname] = gb
    modrms = {g: {k: v for k, v in d.items() if k.startswith("_mod")} for g, d in res.items()}
    for g in res:
        for k in list(res[g]):
            if k.startswith("_mod"):
                del res[g][k]
    s_fb = res["fb01"]["sampled"]
    b = s_fb["bins"]

    def leg_bins(key, gname):                    # the CHK_AMP_BINS stats of a group
        if amp_bins == CHK_AMP_BINS:
            return res[gname]["sampled" if key == "s" else "two_tau"]["bins"]
        return _bin_stats(ratio[key], _groups(fb)[gname], amp)

    bl = leg_bins("s", "fb01")
    ok_med = [abs(d["median"] - 1.0) <= 0.005 for d in bl[1:] if d["n"] >= CHK_MIN_N]
    have = [d for d in bl[1:] if d["n"] >= CHK_MIN_N]
    b2b = abs(bl[1]["median"] - bl[2]["median"]) if len(have) == 2 else float("nan")
    allmed = [d["median"] for d in bl if d["n"] >= CHK_MIN_N]
    rmsm = s_fb["per_module_rms"]["amp300"]["module_median_spread_rms"]
    rms_stat = s_fb["per_module_rms"]["amp300"]["module_median_stat_expected"]
    rms_corr = s_fb["per_module_rms"]["amp300"]["module_median_spread_corrected"]
    cond = [all(ok_med) and len(ok_med) == 2, b2b <= 0.005 if math.isfinite(b2b) else False,
            rmsm <= 0.004 if math.isfinite(rmsm) else False]
    flat = res["fb01"]["flatness_with_fA"]
    flat0 = res["fb01"]["flatness_without_fA"]
    result = "N/A" if flat["n_bins"] < 2 else ("PASS" if flat["flat"] else "FAIL")
    kap = 1.0 / bl[1]["median"] if bl[1]["n"] and math.isfinite(bl[1]["median"]) else float("nan")
    ba = leg_bins("s", "all")
    kap_all = 1.0 / ba[1]["median"] if ba[1]["n"] else float("nan")
    res["kappa_candidate_sampled_300_1000_fb01"] = kap
    res["kappa_candidate_sampled_300_1000_all"] = kap_all
    res["conditions_note"] = ("informational (the pre-f(A) absolute criteria); the gate "
                              "is the flatness rule with f(A) on fallback-0/1 modules")
    res["ampcorr"] = bank_s.ampcorr_str()
    res["ampcorr_meta"] = dict(getattr(bank_s, "ampcorr_meta", {}) or {})
    kn = res["ampcorr_meta"].get("knots") or {}
    res["conditions"] = dict(median_300_1000_and_1000_2500_within_0p005=cond[0],
                             bin_to_bin_300_1000_vs_1000_2500=b2b,
                             bin_to_bin_ok=cond[1],
                             max_bin_spread_all_three_bins=(max(allmed) - min(allmed)
                                                            if len(allmed) > 1 else float("nan")),
                             per_module_spread_rms_amp300=rmsm, per_module_rms_ok=cond[2],
                             per_module_stat_expected_amp300=rms_stat,
                             per_module_spread_quadrature_corrected_amp300=rms_corr,
                             within_module_rms_median_amp100=s_fb["per_module_rms"]["amp100"][
                                 "within_module_rms_median"],
                             per_module_rms_note="gate = RMS across modules of the per-module "
                                                 "median ratio at A > 300 (note section 6.1 "
                                                 "definition, today 0.6 %); the literal "
                                                 "within-module scatter (~1.8 % for both "
                                                 "banks, pulse-level noise) is reported "
                                                 "beside it")
    tt = leg_bins("t", "fb01")
    if amp_bins != CHK_AMP_BINS:
        medc = [d["median"] for d in b if d["n"] >= CHK_MIN_N]
        res["amp_bins"] = [list(x) for x in amp_bins]
        res["conditions"]["max_bin_spread_amp_bins"] = (max(medc) - min(medc)
                                                        if len(medc) > 1 else float("nan"))
    res["two_tau_reproduces_deficit_100_1000"] = bool(
        all(0.975 <= d["median"] <= 0.996 for d in tt[:2] if d["n"] >= CHK_MIN_N)) \
        if all(d["n"] >= CHK_MIN_N for d in tt[:2]) else None
    G.add("G1-7", "the deficit",
          f"median Q_template/Q_raw flat across "
          f"{' / '.join(_bin_label(*x) for x in amp_bins)} ADC "
          "(fallback-0/1 modules, bins with >= 30 pulses): every pairwise difference of "
          "the bin medians <= 2 x its combined SE (SE = 1.2533 std / sqrt(n)); "
          "Q_template = pulse_energy_units(bank, m, t0_fit, A_fit, kappa=1) = A I_int "
          "f(A), Q_raw = candidate wa_integral; the same without f(A), the module-median "
          "spread (RMS across modules, >= 10 pulses, A > 300), the older absolute "
          "criteria (medians within 1 +- 0.005 at >= 300 ADC, bin-to-bin <= 0.5 %, "
          "spread <= 0.4 %) and the two-tau medians are reported beside it",
          result, res,
          "fb01 sampled median with f(A): " + ", ".join(
              f"{d['bin']}: {d['median']:.4f}+-{d['se']:.4f} (n={d['n']})" for d in b)
          + "; without f(A): " + ", ".join(
              f"{d['median']:.4f}+-{d['se']:.4f}"
              for d in res["fb01"]["sampled_no_fA"]["bins"])
          + "; pairwise diffs with f(A) (%): " + ", ".join(
              f"{q['a']} vs {q['b']} {q['diff_pct']:+.2f} ({q['n_sigma']:+.1f} se)"
              for q in flat["pairs"]) + f" flat={flat['flat']}; without: " + ", ".join(
              f"{q['diff_pct']:+.2f} ({q['n_sigma']:+.1f} se)" for q in flat0["pairs"])
          + f" flat={flat0['flat']}"
          + ("; f(A) band knots (" + str(res["ampcorr_meta"].get("scope")) + "): " + ", ".join(
              f"{b} A={k['A']:.0f} f={k['f']:.4f}+-{k['se']:.4f} (n_mod {k['n_modules']})"
              for b, k in kn.items()) if kn else "")
          + f"; module-median spread RMS {rmsm * 100:.2f} % (stat expected {rms_stat * 100:.2f} %, "
          f"quadrature-corrected {rms_corr * 100:.2f} %; within-module scatter "
          f"{s_fb['per_module_rms']['amp100']['within_module_rms_median'] * 100:.2f} %); two-tau fb01: "
          + ", ".join(f"{d['median']:.4f}" for d in tt) + f"; kappa candidate {kap:.4f}")
    _write_tab(outd / "g17_ratio_bins.tsv", ["group", "template", "amp_bin", "n", "median_ratio", "se"], rows)
    mrows = []
    for gname in ("fb01", "fb2"):
        for lab in ("sampled", "two_tau"):
            for v in modrms[gname].get("_modrms_" + lab, []):
                mrows.append([gname, lab, "within_module_rms_amp100", v])
            for v in modrms[gname].get("_modmed_" + lab, []):
                mrows.append([gname, lab, "module_median_ratio_amp300", v])
    _write_tab(outd / "g17_module_stats.tsv", ["group", "template", "quantity", "value"], mrows)
    if png:
        fig, ax = png.subplots(1, 2, figsize=(10, 3.8))
        e = np.geomspace(amp_bins[0][0], amp_bins[-1][1], 12)
        for key, lab, c in (("s", "sampled", "C0"), ("t", "two-tau", "C3")):
            m = (fb >= 0) & (fb <= 1) & np.isfinite(ratio[key])
            xc, md, lo_, hi_ = [], [], [], []
            for a0, a1 in zip(e[:-1], e[1:]):
                mm = m & (amp >= a0) & (amp < a1)
                if mm.sum() >= 10:
                    v = ratio[key][mm]
                    xc.append(np.sqrt(a0 * a1))
                    md.append(np.median(v))
                    lo_.append(np.percentile(v, 16))
                    hi_.append(np.percentile(v, 84))
            ax[0].plot(xc, md, "o-", color=c, label=lab + (" (f(A))" if key == "s" else ""))
            ax[0].fill_between(xc, lo_, hi_, color=c, alpha=0.15)
            if key == "s":
                xn, mn = [], []
                for a0, a1 in zip(e[:-1], e[1:]):
                    mm = m & (amp >= a0) & (amp < a1) & np.isfinite(ratio_nof)
                    if mm.sum() >= 10:
                        xn.append(np.sqrt(a0 * a1))
                        mn.append(np.median(ratio_nof[mm]))
                ax[0].plot(xn, mn, "s--", color="C2", label="sampled, no f(A)")
        ax[0].axhspan(0.995, 1.005, color="g", alpha=0.12)
        ax[0].set_xscale("log")
        ax[0].set_xlabel("candidate amplitude (ADC)")
        ax[0].set_ylabel("Q_template / Q_raw  (median, 16-84 %)")
        ax[0].legend()
        for lab, c in (("sampled", "C0"), ("two_tau", "C3")):
            v = modrms["fb01"].get("_modmed_" + lab, [])
            if v:
                ax[1].hist(v, bins=np.linspace(0.94, 1.04, 41), histtype="step", color=c,
                           label=lab)
        ax[1].axvline(1.0, color="k", ls="--")
        ax[1].set_xlabel("per-module median Q_template / Q_raw (A > 300, fb01)")
        ax[1].set_ylabel("modules")
        ax[1].legend()
        fig.tight_layout()
        fig.savefig(outd / "g17_deficit.png", dpi=110)
        png.close(fig)


def gate_g18(G, z, run_type, outd, png):
    fbs = z["fallback"].astype(int)
    n = np.bincount(fbs, minlength=3)
    n01 = int(n[0] + n[1])
    need = CHK_COVERAGE_RULE.get(run_type)
    result = "REPORT" if need is None else ("PASS" if n01 >= need else "FAIL")
    G.add("G1-8", "coverage",
          "fallback <= 1 for >= 1100 modules (2.24 H2, 3.49 H2), >= 1000 (2.24 C); "
          "other run types reported", result,
          dict(run_type=run_type, n_rows=int(fbs.size), n_fallback0=int(n[0]),
               n_fallback1=int(n[1]), n_fallback2=int(n[2]), n_fallback_le1=n01,
               threshold=need),
          f"{run_type}: fallback 0/1/2 = {n[0]}/{n[1]}/{n[2]}  (<=1: {n01}, "
          f"threshold {need})")
    _write_tab(outd / "g18_coverage.tsv",
               ["module", "module_id", "fallback", "band_lo", "band_hi", "n_cand", "n_pool"],
               [[z["names"][r], z["module_id"][r], fbs[r], z["band_lo"][r], z["band_hi"][r],
                 z["n_cand"][r], z["n_pool"][r]] for r in range(fbs.size)])
    if png:
        fig, ax = png.subplots(figsize=(4.5, 3.5))
        ax.bar(["0", "1", "2"], n, color=["C2", "C1", "C3"])
        ax.set_xlabel("fallback status")
        ax.set_ylabel("modules")
        fig.tight_layout()
        fig.savefig(outd / "g18_coverage.png", dpi=110)
        png.close(fig)


def _lohi_shift(hi, lo, win=CHK_MAXDIFF_WIN, span=1.0, step=0.02):
    """Best-fit shift d minimising sum_t (hi(t) - lo(t + d))^2 on `win`
    (d > 0: the high band is EARLIER than the low band); also the max |diff| on
    the window before and after the shift."""
    m = (T_GRID >= win[0]) & (T_GRID <= win[1])
    t = T_GRID[m]
    hv = hi[m].astype(np.float64)
    best = (np.inf, 0.0)
    for d in np.arange(-span, span + 1e-9, step):
        v = np.interp(t + d, T_GRID, lo.astype(np.float64))
        s = float(((hv - v) ** 2).sum())
        if s < best[0]:
            best = (s, float(d))
    d0 = best[1]
    after = float(np.abs(hv - np.interp(t + d0, T_GRID, lo.astype(np.float64))).max())
    before = float(np.abs(hv - lo[m]).max())
    return d0, before, after


def gate_g19(G, z, P, bank_s, wave, noise_rows, rng, outd, png, n_max):
    """Amplitude bands, the lo/hi time shift and the half-vs-half noise floor."""
    global SMOOTH_LAMBDA
    lam0 = SMOOTH_LAMBDA
    lam_row = json.loads(str(z["meta"])).get("cuts", {}).get("build", {}) \
        .get("lambda_select", {}).get("lam")      # per-module lambda (--lambda-select heldout)
    both = np.isfinite(z["lohi_maxdiff"]) & np.isfinite(z["lohi_iint_ratio"])
    nb = int(both.sum())
    md = z["lohi_maxdiff"][both]
    ir = z["lohi_iint_ratio"][both]
    ok = (md <= 0.01) & (np.abs(ir - 1.0) <= 0.005)
    frac = float(ok.mean()) if nb else float("nan")
    result = "N/A" if nb == 0 else ("PASS" if frac >= 0.9 else "FAIL")
    res = dict(n_modules_with_both_bands=nb, n_rows=int(z["names"].size),
               frac_pass=frac, n_fail_maxdiff=int((md > 0.01).sum()),
               n_fail_iint=int((np.abs(ir - 1.0) > 0.005).sum()),
               maxdiff_median=_med(md), maxdiff_p90=float(np.percentile(md, 90)) if nb else float("nan"),
               iint_ratio_dev_median=_med(np.abs(ir - 1.0)),
               iint_ratio_median=_med(ir),
               iint_ratio_median_se=(float(1.2533 * ir.std(ddof=1) / math.sqrt(nb))
                                     if nb > 1 else float("nan")),
               iint_ratio_note="signed median of iint_hi / iint_lo over the modules; "
                               "its standard error is the sampling error of that "
                               "median (a bias well beyond it is not noise)")
    # best-fit time shift between shape_hi and shape_lo
    rows_both = np.nonzero(both)[0]
    sh = np.array([_lohi_shift(z["shape_hi"][r], z["shape_lo"][r]) for r in rows_both]) \
        if nb else np.zeros((0, 3))
    if nb:
        res["time_shift_hi_vs_lo"] = dict(
            definition="hi(t) ~ lo(t + d) on [0, 20] ns, d > 0: hi band earlier",
            median_ns=float(np.median(sh[:, 0])), p16_ns=float(np.percentile(sh[:, 0], 16)),
            p84_ns=float(np.percentile(sh[:, 0], 84)),
            frac_positive=float((sh[:, 0] > 0).mean()),
            frac_abs_gt_0p1ns=float((np.abs(sh[:, 0]) > 0.1).mean()),
            maxdiff_before_median=float(np.median(sh[:, 1])),
            maxdiff_after_shift_median=float(np.median(sh[:, 2])),
            frac_maxdiff_le1pct_after_shift=float((sh[:, 2] <= 0.01).mean()))
    # noise floor: half-vs-half of the hi pool, same procedure as build
    nf_rows = []
    if P is not None and noise_rows.size:
        shift = z["onset_shift_ns"].astype(np.float64)
        o = np.argsort(P["midx"], kind="stable")
        pm = P["midx"][o]
        for r in noise_rows:
            m = int(z["module_index"][r])
            SMOOTH_LAMBDA = lam0 if lam_row is None or lam_row[r] is None \
                else float(lam_row[r])
            lo_i, hi_i = np.searchsorted(pm, m, "left"), np.searchsorted(pm, m, "right")
            idx = o[lo_i:hi_i]
            if idx.size > n_max:
                idx = rng.choice(idx, n_max, replace=False)
            if idx.size < 2 * BAND_MIN:
                continue
            idx = rng.permutation(idx)
            h = idx.size // 2
            ref = LutShape(bank_s.lut[m].astype(np.float64))
            sh_ = []
            for part in (idx[:h], idx[h:2 * h]):
                ps = PulseSet(P["samples"][part], P["ped"][part], P["rms"][part],
                              P["mf_t0"][part] + shift[r])
                br = band_shape(ps, ref)
                if br is None:
                    break
                sh_.append((spline_values(br["fin"]["t"], br["fin"]["c"], T_GRID),
                            derived_from_spline(br["fin"]["t"], br["fin"]["c"], wave)["iint_wa"]))
            if len(sh_) < 2:
                continue
            w = (T_GRID >= CHK_MAXDIFF_WIN[0]) & (T_GRID <= CHK_MAXDIFF_WIN[1])
            mdiff = float(np.abs(sh_[0][0][w] - sh_[1][0][w]).max())
            nhi, nlo = int(z["n_hi"][r]), int(z["n_lo"][r])
            # 1/n scaling to the real lo-vs-hi comparison of this module
            scale = math.sqrt((1.0 / max(min(nhi, n_max), 1) + 1.0 / max(min(nlo, n_max), 1))
                              * min(nhi, n_max) / 4.0) if nlo >= BAND_MIN else float("nan")
            nf_rows.append([z["names"][r], int(idx.size), h, mdiff, mdiff * scale,
                            float(z["lohi_maxdiff"][r]) if both[r] else float("nan"),
                            abs(sh_[0][1] / sh_[1][1] - 1.0)])
    SMOOTH_LAMBDA = lam0
    if nf_rows:
        a = np.array([x[1:] for x in nf_rows], dtype=np.float64)
        res["noise_floor_halves"] = dict(
            n_modules=len(nf_rows),
            definition="b_hi pool of the module (held-in, tier>=1, 400-2500 ADC, <= n_max "
                       "pulses) split into two random halves, each built with band_shape "
                       "(the build code path); max |half1 - half2| on [0, 20] ns",
            median_pool=float(np.median(a[:, 0])), median_half=float(np.median(a[:, 1])),
            median_maxdiff_half_vs_half=float(np.median(a[:, 2])),
            p84_maxdiff_half_vs_half=float(np.percentile(a[:, 2], 84)),
            median_predicted_lohi_noise=_med(a[:, 3]),
            median_actual_lohi_maxdiff_same_modules=_med(a[:, 4]),
            predicted_note="half-vs-half scaled by sqrt((1/n_hi + 1/n_lo) n_hi / 4) "
                           "(1/n noise scaling to the lo-vs-hi comparison)",
            median_iint_ratio_dev_half_vs_half=_med(a[:, 5]))
    G.add("G1-9", "amplitude dependence",
          "max |shape_hi - shape_lo| on [0, 20] ns <= 1 % of peak and |iint ratio - 1| "
          "<= 0.5 % for >= 90 % of the modules with both bands; the noise floor of "
          "the max-difference and the lo/hi time shift are reported", result, res,
          f"{nb} modules with both bands, {frac * 100:.1f} % pass; lohi_maxdiff median "
          f"{res['maxdiff_median'] * 100:.2f} %"
          + (f"; noise floor (half vs half) median {res['noise_floor_halves']['median_maxdiff_half_vs_half'] * 100:.2f} % "
             f"(n={len(nf_rows)}), predicted lo-vs-hi noise "
             f"{res['noise_floor_halves']['median_predicted_lohi_noise'] * 100:.2f} %" if nf_rows else "")
          + f"; iint_hi/iint_lo median {res['iint_ratio_median']:.4f} +- {res['iint_ratio_median_se']:.4f}"
          + (f"; lo/hi time shift median {res['time_shift_hi_vs_lo']['median_ns']:+.3f} ns, "
             f"max diff after shift {res['time_shift_hi_vs_lo']['maxdiff_after_shift_median'] * 100:.2f} %"
             if nb else ""))
    _write_tab(outd / "g19_modules.tsv",
               ["module", "fallback", "n_lo", "n_hi", "lohi_maxdiff", "lohi_iint_ratio",
                "shift_hi_vs_lo_ns", "maxdiff_before", "maxdiff_after_shift"],
               [[z["names"][r], z["fallback"][r], z["n_lo"][r], z["n_hi"][r],
                 z["lohi_maxdiff"][r], z["lohi_iint_ratio"][r], *sh[k]]
                for k, r in enumerate(rows_both)])
    _write_tab(outd / "g19_noise_floor.tsv",
               ["module", "n_pool", "n_half", "maxdiff_half_vs_half", "predicted_lohi_noise",
                "actual_lohi_maxdiff", "iint_ratio_dev_half_vs_half"], nf_rows)
    if png:
        fig, ax = png.subplots(1, 2, figsize=(10, 3.8))
        if nb:
            ax[0].hist(md * 100, bins=np.linspace(0, 6, 49), histtype="step", color="C0",
                       label="lohi_maxdiff (npz)")
        if nf_rows:
            ax[0].hist(np.array([x[3] for x in nf_rows]) * 100, bins=np.linspace(0, 6, 49),
                       histtype="step", color="C3", label="half vs half")
            ax[0].hist(np.array([x[4] for x in nf_rows]) * 100, bins=np.linspace(0, 6, 49),
                       histtype="step", color="C2", label="half vs half, scaled to lo/hi")
        ax[0].axvline(1.0, color="k", ls="--")
        ax[0].set_xlabel("max |diff| on [0, 20] ns (% of peak)")
        ax[0].set_ylabel("modules")
        ax[0].legend(fontsize=7)
        if nb:
            ax[1].hist(sh[:, 0], bins=np.linspace(-1, 1, 41), color="C0")
        ax[1].set_xlabel("best-fit shift of hi vs lo (ns; > 0 = hi earlier)")
        ax[1].set_ylabel("modules")
        fig.tight_layout()
        fig.savefig(outd / "g19_bands.png", dpi=110)
        png.close(fig)


# ------------------------------------------------------ cross-check mode ----
def _raw_events_key(f):
    cyc = [int(m.group(1)) for k in f.keys() if (m := re.fullmatch(r"events;(\d+)", k))]
    return "events;%d" % max(cyc)


def raw_chunks(path, step=5000, stop=None):
    """uproot reader of a `_raw.root` split (port of the scratch read_raw_24386.py):
    yields, per block of `step` entries, flat per-channel arrays entry, ich,
    module_id, module_type, ns and samples (uint16 [N, 200], fixed stride; the
    first ns are data).  Tree `events`, highest cycle."""
    import awkward as ak
    import uproot
    br = ["hycal.nch", "hycal.module_id", "hycal.module_type", "hycal.nsamples",
          "hycal.samples"]
    with uproot.open(path) as f:
        t = f[_raw_events_key(f)]
        n_ent = t.num_entries if stop is None else min(stop, t.num_entries)
        for s0 in range(0, n_ent, step):
            s1 = min(s0 + step, n_ent)
            a = t.arrays(br, entry_start=s0, entry_stop=s1)
            nch = ak.to_numpy(a["hycal.nch"]).astype(np.int64)
            tot = int(nch.sum())
            if tot == 0:
                continue
            smp = ak.to_numpy(ak.flatten(a["hycal.samples"], axis=1))
            assert smp.shape == (tot, RAW_STRIDE), smp.shape
            first = np.cumsum(nch) - nch
            yield dict(
                entry=np.repeat(np.arange(s0, s1, dtype=np.int64), nch),
                ich=np.arange(tot, dtype=np.int64) - np.repeat(first, nch),
                module_id=ak.to_numpy(ak.flatten(a["hycal.module_id"])).astype(np.int64),
                module_type=ak.to_numpy(ak.flatten(a["hycal.module_type"])).astype(np.int64),
                ns=ak.to_numpy(ak.flatten(a["hycal.nsamples"])).astype(np.int64),
                samples=smp.astype(np.uint16))


def raw_ped(smp, ns):
    """PASS 1 patched-macro pedestal: mean and population rms of raw samples 0-19
    (min(20, ns)); smp[N, >=20] uint16-valued."""
    np_ = np.minimum(RAW_PED_N, ns)
    assert (np_ == RAW_PED_N).all(), "channel with fewer than 20 samples"
    h = smp[:, :RAW_PED_N].astype(np.float64)
    ped = h.sum(axis=1) / RAW_PED_N
    sp2 = (h * h).sum(axis=1) / RAW_PED_N
    return ped, np.sqrt(np.maximum(0.0, sp2 - ped * ped))


def p1_eval_rows(smp, ns, ped):
    """pass1_evaluate_batch on pedestal-subtracted rows, grouped by nsamples."""
    keys = ("amp", "peak", "first_thr", "rising", "falling", "has_secondary", "sec_max",
            "bg_ok", "bg_mean", "bg_rms", "bg_maxabs", "roughness", "tier", "time_ok",
            "valid_p1")
    n = smp.shape[0]
    out = {k: np.zeros(n, dtype=np.float64) for k in keys}
    for nsv in np.unique(ns):
        m = np.nonzero(ns == nsv)[0]
        if nsv <= 0:
            continue
        W = smp[m, :nsv].astype(np.float64) - ped[m][:, None]
        P = pass1_evaluate_batch(W)
        for k in keys:
            out[k][m] = P[k]
    return out


class RawScan:
    """Full pass over a `_raw.root` split: PASS 1 evaluation of every channel;
    keeps the flat per-row table, the debug-TSV rows (entry < tsv_entries) and the
    raw samples of requested (entry, ich) rows."""

    FIELDS = ("entry", "ich", "module_id", "module_type", "ns", "ped", "ped_rms", "amp",
              "peak", "first_thr", "rising", "falling", "has_secondary", "sec_max",
              "bg_ok", "bg_mean", "bg_rms", "bg_maxabs", "roughness", "tier",
              "time_ok", "valid_p1")

    def __init__(self, path, want_keys, tsv_entries=20000):
        self.tab = {k: [] for k in self.FIELDS}
        self.want = {}
        wk = np.asarray(sorted(want_keys), dtype=np.int64)
        t0 = time.monotonic()
        nrow = 0
        for ch in raw_chunks(path):
            ped, prms = raw_ped(ch["samples"], ch["ns"])
            P = p1_eval_rows(ch["samples"], ch["ns"], ped)
            rec = dict(entry=ch["entry"], ich=ch["ich"], module_id=ch["module_id"],
                       module_type=ch["module_type"], ns=ch["ns"], ped=ped, ped_rms=prms, **P)
            for k in self.FIELDS:
                self.tab[k].append(rec[k])
            if wk.size:
                key = ch["entry"] * 100000 + ch["ich"]
                hit = np.nonzero(np.isin(key, wk))[0]
                for i in hit:
                    self.want[int(key[i])] = (ch["samples"][i].copy(), float(ped[i]),
                                              float(prms[i]), int(ch["ns"][i]))
            nrow += ch["entry"].size
        self.t = {k: np.concatenate(v) for k, v in self.tab.items()}
        del self.tab
        print(f"[raw] {nrow} channel rows, {int(self.t['entry'].max()) + 1} entries, "
              f"PASS 1 evaluated in {time.monotonic() - t0:.0f} s", flush=True)


def _best_per_module(T):
    """The macro's betterCandidate applied to all valid rows in order; returns
    {module_id: row index}."""
    v = np.nonzero(T["valid_p1"] > 0.5)[0]
    order = np.lexsort((v, T["bg_rms"][v], T["roughness"][v], -T["amp"][v],
                        -T["tier"][v], T["module_id"][v]))
    vs = v[order]
    first = np.concatenate([[True], T["module_id"][vs][1:] != T["module_id"][vs][:-1]])
    return {int(T["module_id"][i]): int(i) for i in vs[first]}


def gate_g11(G, T, tsv_path, summ, outd):
    """G1-1: per-(entry, channel) field comparison and best-per-module."""
    M = np.loadtxt(tsv_path, delimiter="\t", skiprows=1, dtype=np.float64)
    cols = ["entry", "ich", "module_id", "module_type", "valid", "tier", "rising",
            "falling", "has_secondary", "amplitude", "peak_sample", "first_thr", "bg_mean",
            "bg_rms", "bg_maxabs", "bg_rough", "peak_in_window", "ped", "ped_rms",
            "second_peak"]
    M = {c: M[:, i] for i, c in enumerate(cols)}
    kM = M["entry"].astype(np.int64) * 100000 + M["ich"].astype(np.int64)
    n_ent = int(M["entry"].max()) + 1
    sel = np.nonzero(T["entry"] < n_ent)[0]
    kP = T["entry"][sel] * 100000 + T["ich"][sel]
    oM, oP = np.argsort(kM, kind="stable"), np.argsort(kP, kind="stable")
    same_rows = kM.size == kP.size and np.array_equal(kM[oM], kP[oP])
    res = dict(n_entries_compared=n_ent, n_rows_macro=int(kM.size), n_rows_python=int(kP.size),
               rows_identical_keys=bool(same_rows))
    if not same_rows:
        G.add("G1-1", "PASS 1 port fidelity", "exact agreement", "FAIL", res,
              f"row sets differ: macro {kM.size}, python {kP.size}")
        return
    Mi = {c: M[c][oM] for c in cols}
    Pi = {k: T[k][sel][oP] for k in RawScan.FIELDS}
    thr, minr, minf = P1_THR, P1_MIN_RISING, P1_MIN_FALLING
    s1 = (Pi["amp"] > thr) & (Pi["time_ok"] > 0.5)
    s2 = s1 & (Pi["first_thr"] >= 0)
    s3 = s2 & (Pi["bg_ok"] > 0.5)
    s4 = s3 & (Pi["rising"] >= minr)
    s5 = s4 & (Pi["falling"] >= minf)
    allr = np.ones(kP.size, dtype=bool)

    def close(a, b):
        return np.abs(a - b) <= 1e-6 * np.maximum(1.0, np.abs(b))

    spec = [("valid", allr, Pi["valid_p1"], Mi["valid"], False),
            ("amplitude", allr, Pi["amp"], Mi["amplitude"], True),
            ("peak_sample", allr, Pi["peak"], Mi["peak_sample"], False),
            ("peak_in_window", allr, Pi["time_ok"], Mi["peak_in_window"], False),
            ("ped", allr, Pi["ped"], Mi["ped"], True),
            ("ped_rms", allr, Pi["ped_rms"], Mi["ped_rms"], True),
            ("first_thr", s1, Pi["first_thr"], Mi["first_thr"], False),
            ("tier", s2, Pi["tier"], Mi["tier"], False),
            ("bg_mean", s3, Pi["bg_mean"], Mi["bg_mean"], True),
            ("bg_rms", s3, Pi["bg_rms"], Mi["bg_rms"], True),
            ("bg_maxabs", s3, Pi["bg_maxabs"], Mi["bg_maxabs"], True),
            ("bg_rough", s3, Pi["roughness"], Mi["bg_rough"], True),
            ("rising", s3, Pi["rising"], Mi["rising"], False),
            ("falling", s4, Pi["falling"], Mi["falling"], False),
            ("has_secondary", s5, Pi["has_secondary"], Mi["has_secondary"], False),
            ("second_peak", s5, Pi["sec_max"], Mi["second_peak"], True)]
    fields, rows, nbad_tot = {}, [], 0
    for nm, st, a, b, flt in spec:
        bad = st & ~(close(a, b) if flt else (a == b))
        fields[nm] = dict(n_compared=int(st.sum()), n_mismatch=int(bad.sum()))
        nbad_tot += int(bad.sum())
        for i in np.nonzero(bad)[0][:200]:
            rows.append([nm, int(kP[oP][i] // 100000), int(kP[oP][i] % 100000),
                         float(a[i]), float(b[i])])
    res["fields"] = fields
    _write_tab(outd / "x_g11_mismatch.tsv", ["field", "entry", "ich", "python", "macro"], rows)
    # best per module over the whole split
    bp = _best_per_module(T)
    sm = summ
    mism, agree, tie, rows_b = [], 0, 0, []
    mod_ids = set(bp) | {int(m) for m, v in zip(sm["module_id"], sm["valid"]) if v == 1}
    kT = T["entry"] * 100000 + T["ich"]
    srt = np.argsort(kT, kind="stable")
    ksorted = kT[srt]
    for mid in sorted(mod_ids):
        i = int(np.nonzero(sm["module_id"] == mid)[0][0]) if (sm["module_id"] == mid).any() else -1
        mv = bool(i >= 0 and sm["valid"][i] == 1)
        pi = bp.get(mid)
        if not mv and pi is None:
            continue
        if mv and pi is not None:
            me, mc = int(sm["best_event"][i]), int(sm["best_channel"][i])
            pe, pc = int(T["entry"][pi]), int(T["ich"][pi])
            if (me, mc) == (pe, pc):
                agree += 1
                continue
            j = srt[np.searchsorted(ksorted, me * 100000 + mc)]
            ident = all(T[k][j] == T[k][pi] for k in ("tier", "amp", "roughness", "bg_rms")) \
                and T["valid_p1"][j] > 0.5
            tie += int(ident)
            mism.append(dict(module_id=mid, macro=(me, mc), python=(pe, pc),
                             tie=bool(ident),
                             macro_tier_amp=(float(T["tier"][j]), float(T["amp"][j])),
                             python_tier_amp=(float(T["tier"][pi]), float(T["amp"][pi]))))
        else:
            mism.append(dict(module_id=mid, macro_valid=mv, python_valid=pi is not None,
                             tie=False))
    res["best_per_module"] = dict(
        n_modules_macro_valid=int((sm["valid"] == 1).sum()), n_modules_python_valid=len(bp),
        n_agree=agree, n_disagree=len(mism), n_disagree_documented_ties=int(tie),
        n_disagree_not_tie=len(mism) - int(tie), disagreements=mism[:50])
    _write_tab(outd / "x_g11_best.tsv", ["module_id", "kind", "macro", "python", "tie"],
               [[m["module_id"], "valid_mismatch" if "macro_valid" in m else "different_pick",
                 m.get("macro", m.get("macro_valid")), m.get("python", m.get("python_valid")),
                 m["tie"]] for m in mism])
    result = "PASS" if (nbad_tot == 0 and res["best_per_module"]["n_disagree_not_tie"] == 0) \
        else "FAIL"
    G.add("G1-1", "PASS 1 port fidelity",
          "pass1_evaluate agrees with the macro debug TSV in every field (compared only "
          "where the macro reached that stage; valid on all rows) for entries < "
          f"{n_ent}, and the macro betterCandidate ranking over the split picks the "
          "macro's best (entry, channel) per module (identical except exact ties)",
          result, res,
          f"{int(kP.size)} rows, mismatches " + ", ".join(
              f"{k}:{v['n_mismatch']}" for k, v in fields.items() if v["n_mismatch"])
          + (" none" if nbad_tot == 0 else "")
          + f"; best per module: {agree} agree, {len(mism)} disagree "
            f"({tie} exact ties)")


def _selector_veto_check(p, rec, rows, midx_list):
    """Run `Selector.flush` (the select code path, no cap) on single waveforms;
    returns, per row, the first sequential cut that rejects it ('' if none)."""
    sel = Selector(p, rec, 10 ** 9, 0, 0)
    for k, (row, mi) in enumerate(zip(rows, midx_list)):
        sel.add(np.asarray(row, dtype=np.uint16), k, k, 0x100, int(mi))
    sel.flush()
    seq = CUTFLOW_SEQ[:CUTFLOW_SEQ.index("stored")]
    out = []
    for mi in midx_list:
        cf = sel.cutflow[int(mi)]
        why = ""
        for nm in seq:
            if cf[sel.ci[nm]] == 0:
                why = nm
                break
        out.append(why)
    return out


def gate_g12_g13(G, p, rec, z, bank_s, summ, raw, cands_path, f_root, outd, png):
    import uproot
    mods = [p.hycal.module(i) for i in range(p.hycal.module_count())]
    midx_of_id = {int(m.id): i for i, m in enumerate(mods)}
    row_of = {int(m): r for r, m in enumerate(z["module_index"])}
    use = [i for i in range(summ["module_id"].size)
           if summ["valid"][i] == 1 and summ["module_type"][i] == ROOT_MODULE_TYPE
           and int(summ["module_id"][i]) in midx_of_id and mods[midx_of_id[int(summ["module_id"][i])]].is_pwo4()]
    recs = []
    for i in use:
        mid = int(summ["module_id"][i])
        key = int(summ["best_event"][i]) * 100000 + int(summ["best_channel"][i])
        w = raw.want.get(key)
        if w is None:
            continue
        recs.append(dict(module_id=mid, midx=midx_of_id[mid], name=mods[midx_of_id[mid]].name,
                         entry=int(summ["best_event"][i]), ich=int(summ["best_channel"][i]),
                         samples=w[0][:NS_EXPECT], ped=w[1], ped_rms=w[2], ns=w[3],
                         amp=float(summ["best_amplitude_mV"][i])))
    nrec = len(recs)
    # ---- G1-2
    if cands_path:
        zc = np.load(cands_path, allow_pickle=False)
        cs = zc["samples"]
        look = {cs[i].tobytes(): i for i in range(cs.shape[0])}
        for r in recs:
            r["match"] = look.get(r["samples"].astype(np.uint16).tobytes(), -1)
        veto = _selector_veto_check(p, rec, [r["samples"] for r in recs],
                                    [r["midx"] for r in recs])
        for r, v in zip(recs, veto):
            r["veto"] = v
        matched = [r for r in recs if r["match"] >= 0]
        un = [r for r in recs if r["match"] < 0]
        by_cut = {}
        for r in un:
            by_cut[r["veto"] or "passes_all_vetoes"] = by_cut.get(r["veto"] or "passes_all_vetoes", 0) + 1
        n_pass = len(matched) + by_cut.get("passes_all_vetoes", 0)
        matched_fail_veto = sum(1 for r in matched if r["veto"])
        frac = len(matched) / n_pass if n_pass else float("nan")
        res = dict(n_modules=nrec, n_matched=len(matched), n_unmatched=len(un),
                   matched_but_fail_a_veto=matched_fail_veto,
                   unmatched_failing_cut=by_cut, n_pass_all_vetoes=n_pass,
                   match_fraction_of_passing=frac,
                   unmatched_label="passes_all_vetoes = not in EVIO candidates (trigger "
                                   "!= 0x100 or event mismatch)",
                   candidates_file=str(cands_path), n_candidates=int(cs.shape[0]))
        G.add("G1-2", "EVIO vs macro candidates",
              "the macro's best waveform per module (matched by the 100 raw samples) is "
              "among the EVIO candidates for >= 95 % of the modules that pass the added "
              "vetoes", "PASS" if frac >= 0.95 else "FAIL", res,
              f"{nrec} modules, {len(matched)} matched, {len(un)} unmatched "
              f"({', '.join(f'{k}: {v}' for k, v in sorted(by_cut.items()))}); "
              f"matched / passing vetoes = {frac * 100:.1f} %")
        _write_tab(outd / "x_g12_modules.tsv",
                   ["module_id", "module", "entry", "ich", "amp", "matched", "first_failed_cut"],
                   [[r["module_id"], r["name"], r["entry"], r["ich"], r["amp"],
                     int(r["match"] >= 0), r["veto"] or "none"] for r in recs])
    # ---- G1-3
    rows3, chi_all, chi_mf = [], [], []
    h = uproot.open(f_root)
    for r in recs:
        row = row_of.get(r["midx"])
        if row is None:
            continue
        shp = LutShape(bank_s.lut[r["midx"]].astype(np.float64))
        y = (r["samples"].astype(np.float64) - r["ped"])[None, :]
        sig = np.array([max(r["ped_rms"], SIG_MIN)])
        pk = int(np.argmax(r["samples"]))
        seed = np.array([DT_NS * pk - float(z["t_pk"][row])])
        al = align_pulses(y, sig, seed, shp, 8.0)
        c2 = float(al["chi2dof"][0]) if al["ok"][0] else float("nan")
        c2m = float(_chi2_mf(y, sig, al["t0"], al["A"], shp)[0]) if al["ok"][0] else float("nan")
        # PASS 1 reference histogram vs our spline
        try:
            href = h[f"refs/ref_{r['module_id']}"].values().astype(np.float64)[None, :NS_EXPECT]
        except Exception:                                           # noqa: BLE001
            href = None
        rms_a = rms_b = float("nan")
        if href is not None:
            pr = int(np.argmax(href[0]))
            a2 = align_pulses(href, np.array([1.0]), np.array([DT_NS * pr - float(z["t_pk"][row])]),
                              shp, 8.0)
            if a2["ok"][0]:
                t0r, Ar = a2["t0"][0], a2["A"][0]
                tau = TS - t0r
                res_ = (href[0] - Ar * shp.eval(tau)) * 100.0
                for lo_, tgt in ((0.0, "a"), (WIN[0], "b")):
                    mm = (tau >= lo_) & (tau <= WIN[1])
                    v = float(np.sqrt(np.mean(res_[mm] ** 2)))
                    if tgt == "a":
                        rms_a = v
                    else:
                        rms_b = v
        clipped = bool(r["samples"].max() >= CLIP_ADC)
        rows3.append([r["module_id"], r["name"], int(z["fallback"][row]), r["amp"], int(clipped),
                      c2, c2m, rms_a, rms_b, float(al["t0"][0]), int(r.get("match", -2) >= 0)])
    a3 = np.array([[x[2], x[3], x[4], x[5], x[6], x[7], x[8]] for x in rows3], dtype=np.float64)
    if a3.size:
        fb3, amp3, clip3, c2, c2m, ra, rb = a3.T
        unc = clip3 == 0
        res = dict(
            n_modules=int(a3.shape[0]),
            chi2_ndf_median=_med(c2), chi2_ndf_median_fb01=_med(c2[fb3 <= 1]),
            chi2_ndf_median_fb2=_med(c2[fb3 == 2]), chi2_ndf_median_unclipped=_med(c2[unc]),
            chi2_ndf_median_mfsigma=_med(c2m), chi2_ndf_median_mfsigma_fb01=_med(c2m[fb3 <= 1]),
            chi2_definition="sigma = max(pre-pulse rms of samples 0-19, 0.3), ndf = n_window - 2, "
                            "window [-20, 160] ns; mfsigma = sigma^2 + (0.03 model)^2",
            amp_median=_med(amp3), n_clipped=int((clip3 > 0).sum()),
            spline_vs_pass1_rms_pct_median=_med(ra),
            spline_vs_pass1_rms_pct_median_fb01=_med(ra[fb3 <= 1]),
            spline_vs_pass1_rms_pct_median_window_m20=_med(rb),
            spline_vs_pass1_rms_pct_p84=float(np.nanpercentile(ra, 84)),
            pass1_rms_definition="RMS of (ref sample - A spline(t - t0)) / max-sample on the 4 ns "
                                 "grid for 0 <= tau <= 160 ns (window_m20: -20 <= tau <= 160), "
                                 "(A, t0) fitted to the macro's refs/ref_<id> TH1D")
        ok = (res["chi2_ndf_median"] <= 2.0) and (res["spline_vs_pass1_rms_pct_median"] <= 1.0)
        G.add("G1-3", "shape agreement within noise",
              "median chi2/ndf of the macro's best waveforms against our spline <= 2 and "
              "median RMS difference between our spline and PASS 1's single-waveform "
              "template (4 ns grid, 0-160 ns) <= 1 % of peak",
              "PASS" if ok else "FAIL", res,
              f"n={a3.shape[0]} chi2/ndf median {res['chi2_ndf_median']:.1f} "
              f"(fb01 {res['chi2_ndf_median_fb01']:.1f}, unclipped "
              f"{res['chi2_ndf_median_unclipped']:.1f}, mf-sigma "
              f"{res['chi2_ndf_median_mfsigma']:.1f}); spline-vs-PASS1 RMS median "
              f"{res['spline_vs_pass1_rms_pct_median']:.2f} %")
    else:
        G.add("G1-3", "shape agreement within noise", "see plan", "N/A", {}, "no waveforms")
    _write_tab(outd / "x_g13_modules.tsv",
               ["module_id", "module", "fallback", "amp", "clipped", "chi2_ndf", "chi2_ndf_mfsigma",
                "rms_vs_pass1_pct_0_160", "rms_vs_pass1_pct_m20_160", "t0_fit", "evio_match"],
               rows3)
    if png and a3.size:
        fig, ax = png.subplots(1, 2, figsize=(10, 3.8))
        ax[0].hist(np.clip(a3[:, 3], 0.1, None), bins=np.logspace(-1, 4, 50), histtype="step",
                   label="sigma = pre-pulse rms")
        ax[0].hist(np.clip(a3[:, 4], 0.1, None), bins=np.logspace(-1, 4, 50), histtype="step",
                   label="mf sigma")
        ax[0].axvline(2, color="k", ls="--")
        ax[0].set_xscale("log")
        ax[0].set_xlabel("chi2/ndf of the macro best waveform vs our spline")
        ax[0].set_ylabel("modules")
        ax[0].legend(fontsize=7)
        ax[1].hist(a3[:, 5][np.isfinite(a3[:, 5])], bins=np.linspace(0, 5, 51))
        ax[1].axvline(1, color="k", ls="--")
        ax[1].set_xlabel("RMS(spline - PASS 1 template), % of peak")
        ax[1].set_ylabel("modules")
        fig.tight_layout()
        fig.savefig(outd / "x_g13_shapes.png", dpi=110)
        png.close(fig)


# ----------------------------------------------------------------- driver ---
def cmd_check(a):
    import _common as C
    import mf_pulse
    t_start = time.monotonic()
    files = sorted({f for pat in a.pulses.split(",") for f in glob.glob(pat.strip())})
    if not files:
        sys.exit(f"[ERROR] no files match {a.pulses}")
    xc = [a.pass1_macro_out, a.raw_root, a.macro_tsv]
    if any(xc) and not all(xc):
        sys.exit("[ERROR] --pass1-macro-out, --raw-root and --macro-tsv go together")
    outd = Path(a.out)
    outd.mkdir(parents=True, exist_ok=True)
    png = None
    if not a.no_png:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as png                     # noqa: F811
    tpl = Path(a.templates)
    z = np.load(tpl, allow_pickle=False)
    tmeta = json.loads(str(z["meta"]))
    run_type = tmeta.get("run_type")
    frac = float(tmeta.get("holdout_frac", HOLDOUT_FRAC))
    warn = []
    hold_keys = None
    if "holdout_keys" in z.files:               # --pool-order energy build
        hold_keys = set(map(tuple, z["holdout_keys"].tolist()))
        print(f"[check] templates carry holdout_keys: using the build's energy-mode "
              f"held-out list ({len(hold_keys)} candidates, "
              f"pool_order {tmeta.get('pool_order')}) instead of holdout_mask",
              flush=True)
    if hold_keys is None and \
            int(tmeta.get("holdout_seed", a.holdout_seed)) != int(a.holdout_seed):
        warn.append(f"--holdout-seed {a.holdout_seed} differs from the build's "
                    f"{tmeta.get('holdout_seed')}: the held-out set overlaps the "
                    "training set")
        print("[WARN] " + warn[-1], flush=True)
    n_max = int(tmeta.get("cuts", {}).get("build", {}).get("n_max", 200))
    global SMOOTH_LAMBDA, FIT_POINTS, SHAPE_SIGMA  # the noise-floor band rebuilds use the build's fit settings
    bld = tmeta.get("cuts", {}).get("build", {})
    SMOOTH_LAMBDA = float(bld.get("smooth_lambda", SMOOTH_LAMBDA))
    FIT_POINTS = bld.get("fit_points", "bins")
    SHAPE_SIGMA = float(bld.get("shape_sigma", 0.0))
    rng = np.random.default_rng(CHK_RNG_SEED)
    z0 = np.load(files[0], allow_pickle=False)
    run0 = int(z0["run"])
    p = C.setup_pipeline(evio_path="", run_num=run0)
    nmod = p.hycal.module_count()
    mods = [p.hycal.module(i) for i in range(nmod)]
    names = [m.name for m in mods]
    pwo4 = [bool(m.is_pwo4()) for m in mods]

    # noise-floor sample of modules (G1-9): both bands present and a big enough pool
    both = np.isfinite(z["lohi_maxdiff"]) & (z["n_hi"] >= 2 * BAND_MIN) & (z["n_lo"] >= BAND_MIN)
    cand_rows = np.nonzero(both)[0]
    noise_rows = np.sort(rng.choice(cand_rows, min(a.noise_modules, cand_rows.size), replace=False)) \
        if cand_rows.size else np.zeros(0, dtype=np.int64)
    noise_mids = z["module_index"][noise_rows].tolist()

    H, ntot, nheld, P, metas = _load_check_inputs(files, a.holdout_seed, frac, noise_mids,
                                                  BAND_HI, a.max_per_module, rng, hold_keys)
    wave = _wave_ctx(metas)
    print(f"[load] {len(files)} files: {ntot} candidates, {nheld} held out "
          + (f"(seed {a.holdout_seed}, frac {frac:g}), " if hold_keys is None else
             "(build's energy-mode list), ") + f"{H['midx'].size} used "
          f"(<= {a.max_per_module} per module)  ({time.monotonic() - t_start:.0f} s)", flush=True)

    h50 = Path(a.h50)
    if not h50.is_absolute():
        h50 = HERE / h50
    column = mf_pulse.resolve_column(str(h50), "auto", prefer="25320")
    tau_rows, h50_sha1, agg_json = mf_pulse.load_h50(str(h50), column)
    tau_agg = mf_pulse.load_aggregate(agg_json)
    bank_s, info_s = mf_pulse.build_bank(names, pwo4, str(tpl), wave)
    bank_t, info_t = mf_pulse.build_bank(names, pwo4, str(h50), wave, column=column)
    print(f"[banks] sampled {tpl.name} ({info_s['n_fallback']} fallback names), two-tau "
          f"{h50.name} column {column} ({len(info_t['fallback'])} aggregate)", flush=True)

    G = Gates()
    t1 = time.monotonic()
    F, fb, bound = fit_heldout(H, z, bank_s, tau_rows, tau_agg, names)
    nfit_fail = {k: int(np.isnan(F[k]["A"]).sum()) for k in F}
    nbound = {k: {g: int((bound[k] & m).sum()) for g, m in _groups(fb).items()} for k in bound}
    print(f"[fit] {H['midx'].size} held-out pulses fitted with both templates in "
          f"{time.monotonic() - t1:.0f} s (failed fits: {nfit_fail}; at the t0 search "
          f"boundary: {nbound})", flush=True)
    fbc = {g: int(m.sum()) for g, m in _groups(fb).items()}
    gate_g14(G, H, F, fb, outd, png)
    gate_g15(G, H, F, fb, outd, png)
    gate_g16(G, p, pwo4, bank_s, bank_t, outd, png)
    gate_g17(G, H, F, fb, bank_s, bank_t, names, outd, png,
             CHK_AMP_BINS if a.amp_bins is None else tuple(zip(a.amp_bins[:-1], a.amp_bins[1:])))
    gate_g18(G, z, run_type, outd, png)
    gate_g19(G, z, P, bank_s, wave, noise_rows, rng, outd, png, n_max)

    if a.raw_root:
        import uproot
        macro = uproot.open(a.pass1_macro_out)
        sm = macro["ref_summary"].arrays(library="np")
        want = {int(e) * 100000 + int(c) for e, c, v in
                zip(sm["best_event"], sm["best_channel"], sm["valid"]) if v == 1}
        raw = RawScan(a.raw_root, want)
        gate_g11(G, raw.t, a.macro_tsv, sm, outd)
        wcfg = p.wave_ana.cfg
        rec = mf_pulse.MFReconstructor(
            bank_t, nsigma=MF_NSIGMA, min_amp=MF_MIN_AMP, kmax=MF_KMAX,
            r_model=MF_R_MODEL, ped_nsamples=int(wcfg.ped_nsamples))
        gate_g12_g13(G, p, rec, z, bank_s, sm, raw, a.evio_cands, a.pass1_macro_out, outd, png)

    header = dict(
        tool="spline_template_extract.py check", schema="spline_check_v1",
        run_type=run_type, templates=str(tpl.resolve()), templates_sha1=info_s["sha1"],
        h50=str(h50), h50_sha1=h50_sha1, h50_column=column,
        pulses=[os.path.abspath(f) for f in files],
        holdout=dict(seed=a.holdout_seed, frac=frac, build_seed=tmeta.get("holdout_seed"),
                     **({} if hold_keys is None else dict(
                         source="templates npz holdout_keys (build --pool-order energy)",
                         n_keys=len(hold_keys))),
                     n_candidates=ntot, n_heldout=nheld, n_used=int(H["midx"].size),
                     max_per_module=a.max_per_module,
                     subsample_note=("held-out pulses are subsampled to at most "
                                     f"{a.max_per_module} per module for G1-4/5/7 (deterministic, "
                                     f"seed {CHK_RNG_SEED}); 0 = no subsampling"),
                     n_by_fallback_group=fbc, failed_fits=nfit_fail,
                     fits_at_t0_search_boundary=nbound),
        warnings=warn, fit=dict(window_ns=list(WIN), sigma="max(rms, 0.3)",
                                t0_seed_sampled="mf_t0 + onset_shift_ns (total shift of the final axis vs the mf seed axis; +-3 ns); "
                                                "fallback-2 rows: + median over fallback-0/1 rows "
                                                "(+-6 ns)",
                                t0_seed_two_tau="mf_t0 (+-3 ns)"),
        wall_s=round(time.monotonic() - t_start, 1), git=_git_short(),
        created=datetime.now().isoformat(timespec="seconds"), argv=sys.argv)
    G.write(outd, header)
    print(f"[done] {outd / 'gates.json'}  {time.monotonic() - t_start:.0f} s", flush=True)
    return 0


# ==================================================================== main ==
def _edges_arg(txt):
    """'E0,E1,...,En' -> tuple of strictly increasing positive floats (n >= 1 bands)."""
    try:
        e = tuple(float(x) for x in txt.split(","))
    except ValueError:
        raise argparse.ArgumentTypeError(f"not a comma-separated list of numbers: {txt!r}")
    if len(e) < 2 or e[0] <= 0 or any(b <= a for a, b in zip(e[:-1], e[1:])):
        raise argparse.ArgumentTypeError("need at least two strictly increasing "
                                         f"positive edges: {txt!r}")
    return e


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("select", help="EVIO split -> candidate npz")
    s.add_argument("evio", help="EVIO split, e.g. prad_025320.evio.00050")
    s.add_argument("--out", required=True, help="output .npz")
    s.add_argument("--trigger", choices=("physics", "physics-any", "rawsum", "all"),
                   default="physics",
                   help="physics = trigger_bits == 0x100 (exact); physics-any = "
                        "(trigger_bits & 0xf00) != 0; rawsum = SSP_RawSum "
                        "event type; all = no gate")
    s.add_argument("--max-events", type=int, default=0,
                   help="stop after this many decoded events (0 = all)")
    s.add_argument("--per-module-cap", type=int, default=150,
                   help="max stored per module and split, applied separately "
                        "to band-b0 candidates and to the rest (first come)")
    s.add_argument("--h50", default="output/pulse_templates_h50.dat")
    s.add_argument("--mf-bank", default="", metavar="PATH",
                   help="matched-filter bank: a spline_templates_v1 templates_<slug>.npz "
                        "(default '' = the h50 two-tau bank of --h50).  Used for the "
                        "MFReconstructor (same nsigma / min_amp / kmax / r_model) and for "
                        "the energy estimate E, where mf_pulse.pulse_energy_units also "
                        "applies the bank's f(A) when the npz carries ampcorr.  Path, sha1 "
                        "and kind are recorded in meta['mf'].  mf_t0 / mf_t0_min are "
                        "converted to the h50 onset axis (t0 - onset_shift_ns of the module, "
                        "as check does), so build/check seeds keep their convention")
    s.add_argument("--per-module-best", type=int, default=0, metavar="K",
                   help="keep, per module and split, the K all-cut candidates with "
                        "the highest energy estimate E (see --calib-file) instead of "
                        "the first-come --per-module-cap (0 = off, cap used); the stored "
                        "set is written in event order, ties go to the earlier candidate; stored / info_p1_shape_ok / "
                        "info_wa_npeaks_ne1 are counted on the final set; adds the "
                        "npz array mf_energy")
    s.add_argument("--emax-frac", type=float, default=None, metavar="F",
                   help="keep only candidates with E < F * E_beam (E_beam = runinfo "
                        "Ebeam of the run; also applied before the cap); E = mf_amp * "
                        "I_int(MF bank: h50 or --mf-bank, MF onset on that bank's axis) * calibration factor (MeV); adds the "
                        "npz array mf_energy")
    s.add_argument("--calib-file", default="", metavar="PATH",
                   help="calibration JSON (list of {name, factor}) for E; factor <= 0 "
                        "or missing -> E is NaN and the candidate is dropped.  Default: "
                        "the file module_timing_dump.RUN_CFG gives for the run")
    s.add_argument("--mf-strict", action="store_true",
                   help="also require exactly one matched-filter pulse (mf_npk == 1), "
                        "applied after c3b; counted separately (n_rej_mf_strict).  With "
                        "the h50 bank nearly every A >= 100 ADC pulse has mf_npk >= 2 "
                        "(template mismatch); use --mf-bank")
    s.add_argument("--max-modules-per-event", type=int, default=0, metavar="N",
                   help="reject every channel of an event with more than N PbWO4 "
                        "channels at A >= --module-count-height-min (0 = off; the LMS "
                        "pulses light up the whole calorimeter).  Mirrors "
                        "fit_pulse_template.py --max-modules-per-event / "
                        "--module-count-height-min (same option names, event rejected "
                        "when the count > N; the threshold there defaults to its --height-min "
                        "(20 ADC), here --module-count-height-min 50) but counts only "
                        "PbWO4 channels with 100 samples, uses A = raw max - "
                        "WaveAnalyzer pedestal instead of the first WaveAnalyzer peak "
                        "height, and applies NO quality check (pile-up, overflow and "
                        "clipped channels count, unlike fit_pulse_template's default "
                        "--module-count-require-clean)")
    s.add_argument("--module-count-height-min", type=float, default=50.0, metavar="H",
                   help="A threshold (ADC) of --max-modules-per-event (default 50)")
    s.add_argument("--exclude-trigger-bits", type=lambda x: int(x, 0), default=0,
                   metavar="MASK",
                   help="reject events with (trigger_bits & MASK) != 0, under every "
                        "--trigger choice (e.g. 0x1000000 = LMS bit); base 0, "
                        "0x... accepted; counted as n_excl_bits")
    s.add_argument("--wassim-ref", action="store_true",
                   help="Wassim single-pulse reference mode (best_ref_shape_fit.C): per PbWO4 "
                        "row only the clip veto (no sample >= 4095) and PASS 1 valid_p1 "
                        "(macro defaults incl. the 150-170 ns peak window) on samples - "
                        "WaveAnalyzer ped; keeps per module ONLY the best valid candidate of "
                        "the split (macro betterCandidate: tier desc, amp desc, roughness asc, "
                        "bg_rms asc; first come on ties).  No mf, no other cuts; adds bg_mean "
                        "/ bg_rms / bg_maxabs / roughness to the npz, meta mode 'wassim_ref'.  "
                        "Not combinable with --per-module-best, --mf-strict, "
                        "--max-modules-per-event, --emax-frac, --exclude-trigger-bits, "
                        "--mf-bank; --per-module-cap is ignored.  Input of build "
                        "--single-pulse")
    s.set_defaults(fn=cmd_select)

    b = sub.add_parser("build", help="candidate npz files -> templates npz/root/pdf")
    b.add_argument("--pulses", required=True,
                   help="glob of candidate files (quote it; several globs may be "
                        "comma-separated), e.g. "
                        "'output/templates_spline/<slug>/pulses/*.npz'")
    b.add_argument("--run-type", required=True, help="slug, e.g. 2p24GeV_H2")
    b.add_argument("--h50", default="output/pulse_templates_h50.dat")
    b.add_argument("--out", required=True,
                   help="output stem: <stem>.npz/.root/.pdf and summary.tsv in "
                        "its directory")
    b.add_argument("--n-min", type=int, default=20,
                   help="candidates (tier>=1, band b0, held-in) for fallback 0")
    b.add_argument("--n-min-relaxed", type=int, default=20,
                   help="candidates (tier>=1, band [100,2500]) for fallback 1")
    b.add_argument("--n-max", type=int, default=200)
    b.add_argument("--seed", type=int, default=1)
    b.add_argument("--holdout-frac", type=float, default=HOLDOUT_FRAC)
    b.add_argument("--pool-order", choices=("random", "energy", "energy-phase"), default="random",
                   help="random = n-max drawn at random, hash hold-out (default).  "
                        "energy = the fit pool is the top n-max by the stored mf_energy "
                        "(candidate files from select --per-module-best/--emax-frac) of "
                        "the module's tier>=1 band-b0 candidates; the hold-out for check "
                        "is the NEXT int(n-max*frac/(1-frac)+0.5) by energy (13 for "
                        "n-max 50, frac 0.2), replacing the hash hold-out, and is saved "
                        "in the npz as holdout_keys (run, split, event_number, module).  "
                        "energy-phase = as energy, but at most ceil(n-max/8) pool pulses per "
                        "0.5 ns bin of the sub-sample phase (mf_t0 mod 4 ns), taken by energy")
    b.add_argument("--smooth-lambda", type=float, default=SMOOTH_LAMBDA,
                   help="spline roughness penalty lam * int f''(t)^2 dt (peak units)")
    b.add_argument("--iterations", type=int, default=N_IT, metavar="N",
                   help="alignment / pooling / spline iterations of the MODULE templates "
                        f"(default {N_IT}, >= {N_IT}).  The worst 20 %% drop is recomputed "
                        "in iterations 1 and 2; from iteration 3 the pool is frozen and "
                        "later iterations only realign it.  With N != 3 or --converge-tol "
                        "adds the npz arrays iter_hist / n_iterations / converged and a "
                        "per-iteration log")
    b.add_argument("--converge-tol", type=float, default=None, metavar="T_NS",
                   help="after --iterations keep iterating until converged or "
                        "--max-iterations; converged (evaluated for every iteration "
                        "k >= 3 against k-1) = median |t0 - t_seed| over the pool < T_NS "
                        "AND max |f_k - f_(k-1)| on the 0.1 ns grid < 1e-3 (unit-peak "
                        "shapes) AND median chi2/dof not decreasing by more than 1 %%.  "
                        "Consistency only, no smoothness criterion (default: off, exactly "
                        "--iterations iterations)")
    b.add_argument("--converge-df", type=float, default=CONV_DF, metavar="X",
                   help="convergence threshold on max |f_k - f_(k-1)| (unit-peak shapes, "
                        f"default {CONV_DF:g}); used with --converge-tol")
    b.add_argument("--lambda-select", choices=("fixed", "heldout"), default="fixed",
                   help="fixed = one global --smooth-lambda (default).  heldout = per module "
                        "the lambda in {0, 1e1, 1e2, 1e3, 1e4} with the lowest mean chi2/dof "
                        "of the module's held-out pulses (1e2 if fewer than 5); adds lam / "
                        "lam_flag / n_ho / ho_chi2_* to summary.tsv, lam / lam_flag / n_ho to the npz meta")
    b.add_argument("--lambda-chi2-sigma", choices=("noise", "full"), default=LAM_CHI2_SIGMA,
                   help="--lambda-select heldout: held-out chi2 variance.  noise = sigma = rms "
                        "(default).  full = sigma^2 + (A*0.15 ns*f')^2 + (A*shape_sigma*f)^2 "
                        "per sample, f = that lambda's final template")
    b.add_argument("--lambda-min", type=float, default=LAM_MIN, metavar="X",
                   help="--lambda-select heldout: the grid is the values >= X (default 0 = all)")
    b.add_argument("--lambda-fallback", type=float, default=LAM_FALLBACK, metavar="Y",
                   help="--lambda-select heldout: lambda used with fewer than 5 (or no "
                        f"comparable) held-out pulses (default {LAM_FALLBACK:g})")
    b.add_argument("--shape-sigma", type=float, default=SHAPE_SIGMA, metavar="S",
                   help="--fit-points direct: adds (S*ref(tau))^2 (ref = alignment reference, "
                        "unit peak) to the point variance (default 0 = off).  Ignored with "
                        "--fit-points bins, whose bin sigma (MAD) already holds the shape "
                        "variation; recorded in the meta when > 0")
    b.add_argument("--osc-escalate", action="store_true",
                   help="module templates with spurious wiggles (extrema of prominence >= 0.2 %% of "
                        "the peak outside +-3 ns of it, 2-150 ns): lambda x10 and rebuild the pool "
                        "until none; still wiggly at 1e4: the run-type average (osc_* columns)")
    b.add_argument("--max-iterations", type=int, default=MAX_IT, metavar="M",
                   help=f"iteration cap with --converge-tol (default {MAX_IT})")
    b.add_argument("--fit-points", choices=("bins", "direct"), default="bins",
                   help="bins = spline fitted to the 0.5 ns bin medians (default).  "
                        "direct = fitted to the pooled points themselves (tau <= T_data; "
                        "y/A with sigma = sqrt((sigma_pulse/A)^2 + (0.15 ns * ref'(tau))^2), "
                        "weight 1/sigma); same lambda, knots, tail extension and "
                        "pseudo-points; the bin stats still give the tail, T_data and the "
                        "shape_err / resid outputs.  Applies to module, band and run-average templates")
    b.add_argument("--bands", type=_edges_arg, default=None, metavar="E0,E1,...,En",
                   help="f(A) amplitude bands [E_k, E_k+1) in ADC (default: the "
                        "lo / mid / top bands 100-400 / 400-1000 / 1000-2500); the npz "
                        "gets band_edges / shape_band / r_band / n_band / A_med_band "
                        "instead of the lo / mid / top arrays")
    b.add_argument("--single-pulse", action="store_true",
                   help="Wassim single-pulse reference bank from `select --wassim-ref` "
                        "files: per module the best candidate over all files (betterCandidate "
                        "order, ties: first file); template R(k) = (W - bg_mean) / max_k(W - "
                        "bg_mean), a not-a-knot cubic interpolating spline through the 100 "
                        "samples (no smoothing / pool / alignment iterations), exponential "
                        "tail only to reach T_END = 400 ns, bank unit peak = spline "
                        "maximum; modules without a candidate: run-type average "
                        "(fallback 2); no f(A) (no ampcorr keys), no holdout_keys, band "
                        "arrays NaN / 0.  Same outputs; not combinable with --pool-order, "
                        "--lambda-select, --osc-escalate, --bands, --iterations, "
                        "--converge-tol, --max-iterations, --fit-points; --n-min / --n-max / "
                        "--smooth-lambda / --holdout-frac have no effect")
    b.add_argument("--no-pdf", action="store_true")
    b.add_argument("--modules", default="",
                   help="debug: comma-separated module names (e.g. W100,W200); "
                        "only these rows are built")
    b.add_argument("--png-dir", default="",
                   help="debug: also save the pages of --png-modules (default the "
                        "first four rows) as PNG here")
    b.add_argument("--png-modules", default="")
    b.set_defaults(fn=cmd_build)

    m = sub.add_parser("ampcorr", help="compare the run types' f(A); write one "
                                       "global f(A) if consistent")
    m.add_argument("--templates", nargs="+", required=True,
                   help="templates_<slug>.npz of the run types (rewritten in place)")
    m.add_argument("--out", required=True, help="summary .json")
    m.add_argument("--per-module-f", default=None, metavar="OUTDIR",
                   help="banks built with --bands: leave --templates untouched, write "
                        "copies OUTDIR/<slug>/templates_<slug>.npz with a per-module "
                        "ampcorr_f [module, ln A grid]")
    m.set_defaults(fn=cmd_ampcorr)

    q = sub.add_parser("check", help="held-out gates G1-4..G1-9 and the PASS 1 "
                                     "cross-check G1-1..G1-3")
    q.add_argument("--pulses", required=True,
                   help="glob of candidate files (the files `build` used; all splits "
                        "of the run type)")
    q.add_argument("--templates", required=True, help="templates_<slug>.npz")
    q.add_argument("--h50", default="output/pulse_templates_h50.dat")
    q.add_argument("--out", required=True, help="output directory")
    q.add_argument("--holdout-seed", type=int, default=1,
                   help="must equal the `build --seed` (warns otherwise)")
    q.add_argument("--max-per-module", type=int, default=100,
                   help="subsample the held-out pulses to this many per module for "
                        "G1-4/5/7 (0 = all)")
    q.add_argument("--noise-modules", type=int, default=100,
                   help="modules for the G1-9 half-vs-half noise-floor estimate")
    q.add_argument("--amp-bins", type=_edges_arg, default=None, metavar="E0,E1,...,En",
                   help="amplitude bin edges (ADC) of G1-7 (default: CHK_AMP_BINS "
                        "100-300 / 300-1000 / 1000-2500); the last bin is closed")
    q.add_argument("--no-png", action="store_true")
    q.add_argument("--pass1-macro-out", default="",
                   help="cross-check: the patched PASS 1 macro ROOT output")
    q.add_argument("--raw-root", default="", help="cross-check: the _raw.root split it ran on")
    q.add_argument("--macro-tsv", default="", help="cross-check: the macro debug TSV")
    q.add_argument("--evio-cands", default="",
                   help="cross-check G1-2: uncapped select npz of the same split")
    q.set_defaults(fn=cmd_check)

    a = ap.parse_args(argv)
    rc = a.fn(a)
    sys.exit(rc or 0)


if __name__ == "__main__":
    main()
