#!/usr/bin/env python3
"""benchmark_hycal_timing_test_20261021.py — legacy vs. timing-coincident HyCal clustering.

Runs WaveAnalyzer once per HyCal channel and feeds the same peak set
into two HyCalCluster instances:

  legacy : seed_time_window ≤ 0, one ModuleHit per module — the
           largest-integral peak within [pre_lo, pre_hi] ns.
  gated  : seed_time_window = W, every in-window peak pushed as a
           separate ModuleHit; HyCalCluster applies the seed-anchored
           coincidence cut during BFS.

Headline selection: events with exactly one reconstructed cluster of
energy > --signal-min MeV (default 500·kE, kE = E_beam/2239.51).  Per path the script counts
matches, histograms the cluster energy, and fits a Gaussian.

With --calibrate the script also fits per-seed-module gain corrections
on an inner-ring high-statistics sample (radius < --inner-radius mm,
seed energy ≥ --min-seed-energy MeV, ≥ --min-per-module events) and
re-fits the recalibrated spectrum with both Gaussian and Crystal-Ball
shapes — the headline observable is then σ_E/E, comparable across the
two paths under matched per-module calibration.  See the technical
note in docs/technical_notes/hycal_clustering for the methodology.

This 20261021 test version adds two things on top of benchmark_hycal_timing.py
(see plan_benchmark_hycal_timing_test_20261021.md):

  LMS gain correction (--lms-corr, --lms-ref).  Per-batch, per-PbWO4-module
    LMS factors are read from the C++ replay_gainCorr output
    (prad_NNNNNN_gain_corr.root, via uproot, imported lazily).  The factor
    is the alpha-normalised (LMS2+LMS3)/2 reference, renormalised to the
    run's first batch (--lms-ref first-batch; ref-run keeps the stored
    scale).  It multiplies every HyCal pulse integral before energize()
    in both the legacy and gated paths; the batch is looked up by the
    TI event number.  --lms-corr off gives all-1.0 factors and does not
    need uproot.

  Elastic-window gain refinement (--cal-method elastic-window, default).
    Per seed module: Crystal-Ball pre-fit of the LMS-corrected cluster
    energy on mu_global·[1-0.15/kS, 1+0.10/kS] (mu_global = CB mu of the pool
    spectrum) with a Moller-aware lower edge, fixed window
    mu_M +- nsigma*sigma_M, in-window median -> gain so that the median
    lands on the elastic energy E'(theta_M) of the run's target (--target
    or the RUN_TARGET map: H2 = e-p, C / Ta = e-nucleus recoil kinematics;
    Python port of PhysicsTools::ExpectedEnergy incl. EnergyLoss).  Modules
    whose elastic-Moller gap is below 3 sigma_M are not used.
    Single pass: the window is fixed in pre-refinement energy, so a second
    pass would return gain 1 exactly.  --cal-method cb-peak / median are
    the unchanged old refinements (use --target-peak, 2 iterations).

  All beam energies.  Thresholds, bin width, histogram range and
  --target-peak default from Ebeam (kE = Ebeam/2239.51, kS = sqrt(kE)):
  500·kE signal cut, 500..3500·kE range, 200·kE min seed energy, 20·kS MeV
  bins, target = Ebeam; peak-shape constants scale with kS.  Explicit CLI
  values override.  The effective calibration file's energy tag (0p7, 2p1,
  2p2, 3p5, x17) is checked against Ebeam: a mismatching runinfo default is
  an error (exit 2), a mismatching explicit --hc-calib-file only a warning.

Outputs:
  <out>.tsv               — summary + histogram bins
  <out>.png               — energy spectrum + Gaussian fits
  <out>_calibrated.png    — recalibrated spectrum + CB fit (--calibrate only)
  <out>_clusters.tsv      — per-event table (--write-cluster-tsv)
  <out>_lms.tsv / _lms.png — per-batch LMS factor summary (png: not for off)
  <out>_window.tsv        — per-module elastic-window table
                            (--calibrate --cal-method elastic-window)
"""

from __future__ import annotations

import argparse
import bisect
import math
import os
import sys
import time
from dataclasses import dataclass, field

import numpy as np
import matplotlib.pyplot as plt

import _common as C
from prad2py import dec, det  # noqa: E402

try:
    from scipy.optimize import curve_fit
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False


def gauss(x, A, mu, sigma):
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def crystal_ball(x, A, mu, sigma, alpha, n):
    """Crystal-Ball PDF (left-side power-law tail).  α > 0, n > 1."""
    z = (np.asarray(x, dtype=float) - mu) / sigma
    abs_a = abs(alpha)
    A_coef = (n / abs_a) ** n * np.exp(-0.5 * abs_a ** 2)
    B_coef = n / abs_a - abs_a
    gauss_part = np.exp(-0.5 * z * z)
    arg = B_coef - z
    arg = np.where(arg > 0, arg, 1e-30)
    tail_part = A_coef * arg ** (-n)
    return A * np.where(z > -alpha, gauss_part, tail_part)


# ---- elastic kinematics: double-precision port of PhysicsTools.cpp:284-332 ----
M_PROTON, M_ELECTRON, DEG2RAD = 938.272, 0.51099895, 3.14159265 / 180.0
U_AMU = 931.494                      # atomic mass unit (MeV)
TARGET_MASS = {'H2': M_PROTON,       # e-p
               'C':  12.0 * U_AMU,   # e-12C elastic
               'Ta': 180.948 * U_AMU}  # e-Ta elastic
# Elastic target per run (documented map; --target overrides).  The 0.73 and
# 3.49 GeV runs are e-p runs per the user.
RUN_TARGET = {**{r: 'H2' for r in (24386, 24503, 24512, 24578, 24590,
                                   24603, 25228, 25257, 25260, 25299,
                                   25320)},
              25308: 'C',
              **{r: 'Ta' for r in (25800, 26138, 26230)}}

E_REF = 2239.51                      # reference beam energy (MeV) of the
                                     # default (2.24 GeV) thresholds/shapes
# Calibration-file energy tags (basename substring -> beam energy, MeV).
CALIB_TAGS = (("x17", 2239.5), ("0p7", 728.9), ("2p2", 2239.5),
              ("2p1", 2108.8), ("3p5", 3488.4))
CALIB_MATCH_TOL = 0.03               # 2p1 (2108.8) vs 2p2 (2239.5) differ by 5.8%

# Elastic-window refinement tunables (module constants, not flags).
PREFIT_RANGE   = (-0.15, 0.10)   # per-module CB pre-fit range, fractions of mu_global (at E_REF;
                                 # divided by kS at other energies)
MIN_FIT_EVENTS = 100             # min events in the pre-fit range for a CB fit
MOLLER_GAP_NSIGMA = 3.0          # exclude a module if its elastic-Moller gap < this × sigma_M


def energy_loss(theta_deg):
    """PhysicsTools::EnergyLoss — target/window material loss (MeV)."""
    c = math.cos(theta_deg * DEG2RAD)
    sec = 1.0 / c if c > 0.01 else 100.0
    return (0.5 * 1.6 + 0.12 * 1.6 + 0.1 * 2.0 + 0.48 * 1.8) * sec


def expected_energy_elastic(theta_deg, ebeam, m_target):
    """Elastic recoil kinematics E·M/(M+E(1-cosθ)) minus EnergyLoss.
    C++ PhysicsTools has no nucleus mode, so for C / Ta the same
    EnergyLoss as the "ep" mode is applied."""
    cos_t = math.cos(theta_deg * DEG2RAD)
    e = ebeam * m_target / (m_target + ebeam * (1.0 - cos_t))
    return e - energy_loss(theta_deg)


def expected_energy_ep(theta_deg, ebeam):
    """PhysicsTools::ExpectedEnergy(..., "ep"): elastic e-p minus EnergyLoss."""
    return expected_energy_elastic(theta_deg, ebeam, M_PROTON)


def expected_energy_ee(theta_deg, ebeam):
    """PhysicsTools::ExpectedEnergy(..., "ee"): Moller minus EnergyLoss;
    0 if the denominator is not positive."""
    cos_t = math.cos(theta_deg * DEG2RAD)
    gamma = ebeam / M_ELECTRON
    num = (gamma + 1.0) + (gamma - 1.0) * cos_t * cos_t
    den = (gamma + 1.0) - (gamma - 1.0) * cos_t * cos_t
    if den <= 0:
        return 0.0
    return M_ELECTRON * num / den - energy_loss(theta_deg)


def module_theta_deg(x, y, hycal_z):
    """Polar angle of a module centre, as physics_calib.cpp: no shower
    depth, no beam offset."""
    return math.degrees(math.atan(math.hypot(x, y) / hycal_z))


def fit_peak_cb(values, e_lo, e_hi, bin_width, mu0_hint=None, escale=1.0):
    """Crystal-Ball fit of values in [e_lo, e_hi].  Returns
    (mu, sigma, alpha, n, A, n_in_fit, edges, counts).  μ seeded from
    `mu0_hint` (else histogram argmax); a wide ±400·escale MeV window
    includes the low-side tail.  `escale` (kS) scales the peak-shape
    constants with the beam energy (1.0 at 2.24 GeV)."""
    edges = np.arange(e_lo, e_hi + bin_width, bin_width)
    counts, _ = np.histogram(values, bins=edges)
    centres = 0.5 * (edges[:-1] + edges[1:])
    nan = float('nan')
    if not HAVE_SCIPY or counts.sum() == 0:
        return nan, nan, nan, nan, nan, 0, edges, counts

    imax = int(np.argmax(counts))
    mu0  = float(mu0_hint) if mu0_hint is not None else float(centres[imax])
    win  = max(400.0 * escale, 12.0 * bin_width)
    mask = (centres > mu0 - win) & (centres < mu0 + 0.5 * win)
    x = centres[mask]
    y = counts[mask].astype(float)
    if y.sum() < 50:
        return nan, nan, nan, nan, nan, int(y.sum()), edges, counts
    try:
        p0     = [float(y.max()), mu0, max(60.0 * escale, bin_width * 3), 1.5, 2.0]
        bounds = ([0,        e_lo - 200.0 * escale, 1.0,  0.1,  1.01],
                  [np.inf,   e_hi + 200.0 * escale, win,  20.0, 50.0])
        popt, _ = curve_fit(crystal_ball, x, y, p0=p0, bounds=bounds,
                            maxfev=10000)
        A, mu, sig, alpha, n = popt
        return (float(mu), abs(float(sig)), float(alpha), float(n),
                float(A), int(y.sum()), edges, counts)
    except Exception as e:
        sys.stderr.write(f"[WARN] crystal-ball fit failed: {e}\n")
        return nan, nan, nan, nan, nan, int(y.sum()), edges, counts


def per_module_peak(values, bin_width, escale=1.0):
    """Per-module peak position for the gain refinement.  The narrow
    ±200 MeV window keeps the per-module radiative tail from pulling
    the fit away from the elastic position; CB lets us include a
    little of the tail without bias.  Tries in order: Crystal-Ball →
    narrow Gaussian → histogram argmax → median (each fall-through
    only on the previous step's failure)."""
    vals = np.asarray(values, dtype=float)
    if len(vals) == 0:
        return float('nan')
    med = float(np.median(vals))
    if len(vals) < 10 or not HAVE_SCIPY:
        return med

    e_lo, e_hi = med - 600.0 * escale, med + 600.0 * escale
    edges  = np.arange(e_lo, e_hi + bin_width, bin_width)
    counts, _ = np.histogram(vals, bins=edges)
    if counts.sum() < 10:
        return med
    centres = 0.5 * (edges[:-1] + edges[1:])
    mu0 = float(centres[int(np.argmax(counts))])

    win  = 200.0 * escale
    mask = (centres > mu0 - win) & (centres < mu0 + 0.5 * win)
    x, y = centres[mask], counts[mask].astype(float)
    if y.sum() < 10:
        return mu0

    try:
        popt, _ = curve_fit(
            crystal_ball, x, y,
            p0     = [float(y.max()), mu0, 50.0 * escale, 1.5, 2.0],
            bounds = ([0, e_lo, 5.0 * escale, 0.1, 1.01],
                      [np.inf, e_hi, win, 20.0, 50.0]),
            maxfev = 5000)
        return float(popt[1])
    except Exception:
        pass
    try:
        popt, _ = curve_fit(
            gauss, x, y,
            p0     = [float(y.max()), mu0, 50.0 * escale],
            bounds = ([0, e_lo, 5.0 * escale], [np.inf, e_hi, win]),
            maxfev = 2000)
        return float(popt[1])
    except Exception:
        pass
    return mu0


def _select_pool(events, *, inner_r, min_seed_E, exclude_rowcol):
    """Calibration pool shared by both refinements: seed radius < inner_r,
    seed energy ≥ min_seed_E, and outside the exclude_rowcol rectangle."""
    inner_r2 = inner_r * inner_r
    excl = exclude_rowcol

    def keep(e):
        _sid, x, y, row, col, se, _ce, *_ = e
        if (x * x + y * y) >= inner_r2 or se < min_seed_E:
            return False
        if excl is not None and excl[0] <= row <= excl[1] and excl[2] <= col <= excl[3]:
            return False
        return True

    return [e for e in events if keep(e)]


def fit_per_module_gains(events, *, target_peak, inner_r, min_seed_E,
                         min_per_mod, n_iter=2, exclude_rowcol=None,
                         method='median', bin_width=20.0, escale=1.0):
    """Per-seed-module gain refinement on the inner-ring high-statistics
    sample.  See the technical note for the methodology.

    `events`         list of (seed_id, x, y, row, col, seed_E, cluster_E,
                     f_eff); f_eff is not used here.
    `exclude_rowcol` (R_LO, R_HI, C_LO, C_HI): drop seeds whose row/col
                     falls inside the rectangle (used to peel the
                     leakage-prone ring around the beam hole).
    `method`         'median' (cheap, biased low by the radiative tail)
                     or 'cb-peak' (per_module_peak()).

    Returns (gains, raw, corrected, n_modules) where `raw` and
    `corrected` are aligned numpy arrays of cluster energies for the
    events whose seed module had ≥ min_per_mod entries."""
    pool = _select_pool(events, inner_r=inner_r, min_seed_E=min_seed_E,
                        exclude_rowcol=exclude_rowcol)
    if not pool:
        return {}, np.array([]), np.array([]), 0

    seed_ids = sorted({e[0] for e in pool})
    gains = {sid: 1.0 for sid in seed_ids}

    for _ in range(n_iter):
        per_mod = {sid: [] for sid in seed_ids}
        for sid, _x, _y, _r, _c, _se, ce, *_ in pool:
            per_mod[sid].append(ce * gains[sid])
        for sid, vals in per_mod.items():
            if len(vals) < min_per_mod:
                continue
            ref = (per_module_peak(vals, bin_width, escale) if method == 'cb-peak'
                   else float(np.median(vals)))
            if ref and ref > 0 and not math.isnan(ref):
                gains[sid] *= target_peak / ref

    counts = {sid: 0 for sid in seed_ids}
    for sid, *_ in pool:
        counts[sid] += 1
    converged = {sid: g for sid, g in gains.items()
                 if counts[sid] >= min_per_mod}
    raw  = np.array([e[6] for e in pool if e[0] in converged])
    corr = np.array([e[6] * converged[e[0]] for e in pool if e[0] in converged])
    return converged, raw, corr, len(converged)


def fit_module_cb(vals, lo, hi, seed, bin_width, sigma_glob, escale=1.0):
    """Crystal-Ball fit of one module's energies on [lo, hi].  Returns
    (mu, sigma, alpha, n_cb, ok).  `ok` needs a clean fit with finite
    parameters, lo < mu < hi and 0.25·sigma_glob ≤ sigma ≤ 2·sigma_glob."""
    nan = float('nan')
    if not HAVE_SCIPY:
        return nan, nan, nan, nan, False
    edges = np.arange(lo, hi + bin_width, bin_width)
    counts, _ = np.histogram(vals, bins=edges)
    if counts.sum() == 0:
        return nan, nan, nan, nan, False
    centres = 0.5 * (edges[:-1] + edges[1:])
    y = counts.astype(float)
    try:
        popt, _ = curve_fit(
            crystal_ball, centres, y,
            p0     = [float(y.max()), seed, 50.0 * escale, 1.5, 2.0],
            bounds = ([0, lo, 5.0 * escale, 0.1, 1.01],
                      [np.inf, hi, hi - lo, 20.0, 50.0]),
            maxfev = 5000)
    except Exception:
        return nan, nan, nan, nan, False
    _A, mu, sig, alpha, n_cb = popt
    sig = abs(float(sig))
    ok = (bool(np.all(np.isfinite(popt))) and lo < mu < hi
          and 0.25 * sigma_glob <= sig <= 2.0 * sigma_glob)
    return float(mu), sig, float(alpha), float(n_cb), ok


def fit_elastic_window_gains(events, *, ebeam, hycal_z, m_target, inner_r,
                             min_seed_E, min_per_mod, exclude_rowcol, nsigma,
                             bin_width, escale=1.0, kE=1.0):
    """Elastic-window per-seed-module gain refinement (single pass).

    Pool as fit_per_module_gains.  E'(theta_M) is the elastic energy for
    target mass `m_target` (expected_energy_elastic).  The pool spectrum is
    fitted globally (Gaussian, then CB, on Ē' ± 800·kE MeV): mu_global =
    CB mu (Gaussian mu if the CB fit fails), sigma_glob = CB sigma
    (Gaussian if it fails).  Per module: CB pre-fit of the (LMS-corrected,
    pre-refinement) cluster energy on
    [mu_global·max(1+LO, (1+E_ee/E')/2), mu_global·(1+HI)] with
    (LO, HI) = PREFIT_RANGE/escale (constant multiples of sigma_exp; the
    Moller-aware lower edge; both measured from mu_global), seeded at
    mu_global; a failed or under-populated fit falls back to centre =
    median in that range and sigma = sigma_glob.  The window
    mu_M ± nsigma·sigma_M is fixed (not clipped); the gain is
    E'(theta_M) / median(in-window events).
    A module is not used (exclude_reason) if
      "moller_gap"        gap_meas < MOLLER_GAP_NSIGMA·sigma_M, with
                          gap_meas = (mu_global/Ē')·(E'_el,M - E_ee,M) (the
                          scale factor puts the gap in measured-energy
                          units, like sigma_M);
      "few_in_window"     n_in_window < min_per_mod;
      "no_prefit_events"  no events in the pre-fit range.
    There is no iteration: the selected event set is fixed in
    pre-refinement energy, so a second pass would return gain 1 exactly.
    `min_per_mod` counts in-window events.  flag_mu_off =
    |mu_M - mu_global| > 2·sigma_M.

    `events` are (seed_id, x, y, row, col, seed_E, cluster_E, f_eff).
    Returns (gains, raw, corr, n_modules, rows, sigma_glob, eprime_mean,
    mu_global); the global offset is mu_global/eprime_mean - 1.
    raw/corr hold all pool events of used modules (uncorrected / times
    gain); `rows` is the per-module table (list of dicts, no `path`)."""
    nan = float('nan')
    pool = _select_pool(events, inner_r=inner_r, min_seed_E=min_seed_E,
                        exclude_rowcol=exclude_rowcol)
    if not pool:
        return {}, np.array([]), np.array([]), 0, [], nan, nan, nan

    per_mod = {}
    for sid, x, y, row, col, _se, ce, feff in pool:
        d = per_mod.setdefault(sid, dict(x=x, y=y, row=row, col=col,
                                         ce=[], ce_nolms=[]))
        d['ce'].append(ce)
        d['ce_nolms'].append(ce / feff)
    seed_ids = sorted(per_mod)

    theta = {s: module_theta_deg(per_mod[s]['x'], per_mod[s]['y'], hycal_z)
             for s in seed_ids}
    eprime = {s: expected_energy_elastic(theta[s], ebeam, m_target)
              for s in seed_ids}
    emoll  = {s: expected_energy_ee(theta[s], ebeam) for s in seed_ids}
    eprime_mean = float(np.mean([eprime[s] for s in seed_ids]))

    # Global mu / sigma from the whole pool (same two-step as the post-cal fit).
    all_ce = np.array([e[6] for e in pool])
    g_lo, g_hi = eprime_mean - 800.0 * kE, eprime_mean + 800.0 * kE
    mu_g, sig_g, *_ = fit_peak(all_ce, g_lo, g_hi, bin_width, escale)
    mu_cb_g, sig_cb_g = fit_peak_cb(all_ce, g_lo, g_hi, bin_width,
                                    mu0_hint=mu_g, escale=escale)[:2]
    sigma_glob = sig_cb_g if not math.isnan(sig_cb_g) else sig_g
    mu_global = mu_cb_g if not math.isnan(mu_cb_g) else mu_g
    if math.isnan(sigma_glob) or math.isnan(mu_global):
        sys.stderr.write("[ERROR] elastic-window: global mu/sigma fit failed; "
                         "skipping the refinement for this path.\n")
        return {}, np.array([]), np.array([]), 0, [], nan, eprime_mean, nan

    lo_f, hi_f = PREFIT_RANGE[0] / escale, PREFIT_RANGE[1] / escale
    gscale = mu_global / eprime_mean
    gains, rows = {}, []
    for sid in seed_ids:
        d = per_mod[sid]
        vals = np.asarray(d['ce'], dtype=float)
        vals_nl = np.asarray(d['ce_nolms'], dtype=float)
        ep, ee = eprime[sid], emoll[sid]
        lo = mu_global * max(1.0 + lo_f, 0.5 * (1.0 + ee / ep))
        hi = mu_global * (1.0 + hi_f)
        sel = vals[(vals >= lo) & (vals <= hi)]
        n_prefit = len(sel)

        row = dict(module=sid, row=d['row'], col=d['col'],
                   r_mm=math.hypot(d['x'], d['y']), theta_deg=theta[sid],
                   eprime=ep, emoller=ee, prefit_lo=lo, prefit_hi=hi,
                   n_total=len(vals), n_prefit=n_prefit,
                   mu=nan, sigma=nan, sigma_source='none',
                   win_lo=nan, win_hi=nan, n_in_window=0, n_in_nolms=0,
                   frac_lms=nan, frac_nolms=nan, median_in=nan,
                   target=ep, gain=nan, flag_mu_off=0, used=0,
                   moller_gap_nsigma=nan, exclude_reason='no_prefit_events')
        rows.append(row)
        if n_prefit == 0:
            continue

        mu_m, sig_m, source = nan, nan, 'global'
        if n_prefit >= MIN_FIT_EVENTS:
            fmu, fsig, _a, _n, ok = fit_module_cb(vals, lo, hi, mu_global,
                                                  bin_width, sigma_glob, escale)
            if ok:
                mu_m, sig_m, source = fmu, fsig, 'module'
        if source == 'global':
            mu_m, sig_m = float(np.median(sel)), sigma_glob

        w_lo, w_hi = mu_m - nsigma * sig_m, mu_m + nsigma * sig_m
        in_win = vals[(vals >= w_lo) & (vals <= w_hi)]
        n_in = len(in_win)
        n_in_nl = int(np.count_nonzero((vals_nl >= w_lo) & (vals_nl <= w_hi)))
        gap_meas = gscale * (ep - ee)
        row.update(mu=mu_m, sigma=sig_m, sigma_source=source,
                   win_lo=w_lo, win_hi=w_hi, n_in_window=n_in,
                   n_in_nolms=n_in_nl,
                   frac_lms=n_in / len(vals), frac_nolms=n_in_nl / len(vals),
                   flag_mu_off=int(abs(mu_m - mu_global) > 2.0 * sig_m),
                   moller_gap_nsigma=gap_meas / sig_m)
        if gap_meas < MOLLER_GAP_NSIGMA * sig_m:
            row['exclude_reason'] = 'moller_gap'
            continue
        if n_in < min_per_mod or n_in == 0:
            row['exclude_reason'] = 'few_in_window'
            continue
        med = float(np.median(in_win))
        if not (med > 0 and math.isfinite(med)):
            row['exclude_reason'] = 'bad_median'
            continue
        row.update(median_in=med, gain=ep / med, used=1, exclude_reason='')
        gains[sid] = ep / med

    raw  = np.array([e[6] for e in pool if e[0] in gains])
    corr = np.array([e[6] * gains[e[0]] for e in pool if e[0] in gains])
    return gains, raw, corr, len(gains), rows, sigma_glob, eprime_mean, mu_global


def fit_peak(values, e_lo, e_hi, bin_width, escale=1.0):
    """Gaussian fit around the histogram argmax in a ±150·escale MeV window.
    Returns (mu, sigma, A, n_in_fit, edges, counts); mu/sigma NaN if
    the fit failed."""
    edges  = np.arange(e_lo, e_hi + bin_width, bin_width)
    counts, _ = np.histogram(values, bins=edges)
    centres = 0.5 * (edges[:-1] + edges[1:])
    nan = float('nan')
    if not HAVE_SCIPY or counts.sum() == 0:
        return nan, nan, nan, 0, edges, counts

    mu0  = float(centres[int(np.argmax(counts))])
    win  = max(150.0 * escale, 4.0 * bin_width)
    mask = (centres > mu0 - win) & (centres < mu0 + win)
    x, y = centres[mask], counts[mask].astype(float)
    if y.sum() < 20:
        return nan, nan, nan, int(y.sum()), edges, counts
    try:
        popt, _ = curve_fit(
            gauss, x, y,
            p0     = [float(y.max()), mu0, max(bin_width, 60.0 * escale)],
            bounds = ([0, e_lo, 1.0], [np.inf, e_hi, win]))
        A, mu, sig = popt
        return float(mu), abs(float(sig)), float(A), int(y.sum()), edges, counts
    except Exception as e:
        sys.stderr.write(f"[WARN] gauss fit failed: {e}\n")
        return nan, nan, nan, int(y.sum()), edges, counts


# ---- LMS gain-drift correction ----------------------------------------------

def _repo_database_dir():
    return os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                         "..", "..", "database"))


def resolve_lms_path(arg, run_num):
    """--lms-corr value -> gain_corr.root path, or None for 'off'.
    'auto' = <PRAD2_DATABASE_DIR>/gain_factor/gain_correction/
    prad_<run:06d>_gain_corr.root (repo database/ without the env var)."""
    if arg == 'off':
        return None
    if arg != 'auto':
        return arg
    db = os.environ.get("PRAD2_DATABASE_DIR") or _repo_database_dir()
    return os.path.join(db, "gain_factor", "gain_correction",
                        f"prad_{run_num:06d}_gain_corr.root")


def _require_uproot():
    try:
        import uproot  # noqa: F401  (lazy: --lms-corr off must not need it)
    except ImportError:
        sys.stderr.write(
            "[ERROR] uproot is required to read the LMS gain_corr file.\n"
            "        Install with:  python3 -m pip install uproot\n"
            "        or run with --lms-corr off.\n")
        sys.exit(2)
    return uproot


@dataclass
class LmsTable:
    """Per-batch, per-module LMS factors, indexed by HyCal module.index."""
    path: str
    ref_mode: str
    ref_run: int
    event_start: list
    event_end: list
    unix_time: list
    n_lms_events: list
    factors: np.ndarray                  # (n_batches, n_modules)
    corrected: np.ndarray                # bool, same shape: factor not defaulted
    pwo4_index: list                     # module.index of PbWO4 rows
    names: dict = field(default_factory=dict)   # module.index -> name
    factors_list: list = field(default_factory=list)

    def __post_init__(self):
        self.factors_list = [row.tolist() for row in self.factors]

    @property
    def n_batches(self):
        return len(self.event_start)

    def batch_of(self, ev):
        """Last batch with event_start <= ev, clipped to [0, n-1] (GetCorr)."""
        b = bisect.bisect_right(self.event_start, ev) - 1
        return 0 if b < 0 else (self.n_batches - 1 if b >= self.n_batches else b)

    def n_uncorrected(self, b):
        idx = self.pwo4_index
        return int(len(idx) - np.count_nonzero(self.corrected[b, idx]))

    def uncorrected_names(self, batches):
        idx = np.asarray(self.pwo4_index, dtype=int)
        bad = ~np.all(self.corrected[np.ix_(list(batches), idx)], axis=0) \
            if len(batches) else np.zeros(len(idx), dtype=bool)
        return [self.names[int(i)] for i in idx[bad]]


def load_lms_gain_corr(path, hycal, ref_mode):
    """Read replay_gainCorr output into an LmsTable.  path=None -> 'off'
    (one batch of 1.0, no uproot).  Row i of gain_corr_W is module id
    1001+i; columns are LMS1, LMS2, LMS3 and (LMS2+LMS3)/2 is the
    reference.  A column is valid if > 0 and != 1.0 (1.0 = failed fit).
    first-batch: A_b/A_0 (both columns, else the single common column,
    else 1.0 "uncorrected").  ref-run: stored scale, 0 -> 1.0 per column,
    as Replay.cpp."""
    n_mod = hycal.module_count()
    if path is None:
        return LmsTable(path="off", ref_mode=ref_mode, ref_run=0,
                        event_start=[0], event_end=[0], unix_time=[0],
                        n_lms_events=[0], factors=np.ones((1, n_mod)),
                        corrected=np.ones((1, n_mod), dtype=bool),
                        pwo4_index=[], names={})

    uproot = _require_uproot()
    with uproot.open(path) as fr:
        t = fr["gain_corr"]
        if t.num_entries == 0:
            print(f"[ERROR] LMS gain_corr file has no batches: {path}\n"
                  f"        (use --lms-corr off to run uncorrected)",
                  file=sys.stderr)
            sys.exit(2)
        g = {k: t[k].array(library="np") for k in
             ("event_num_start", "event_num_end", "n_lms_events",
              "ref_run", "gain_corr_W")}
        # older files have no unix_time branch (C++ loader: has_unix_time)
        g["unix_time"] = (t["unix_time"].array(library="np")
                          if "unix_time" in t.keys()
                          else np.zeros(t.num_entries, dtype=np.int64))
    order = np.argsort(g["event_num_start"], kind="stable")
    start = np.asarray(g["event_num_start"])[order]
    assert np.all(np.diff(start) > 0), "event_num_start not strictly increasing"
    gc = np.asarray(g["gain_corr_W"], dtype=np.float64)[order]
    gc = gc.reshape(len(start), -1, 3)
    nb = len(start)

    factors   = np.ones((nb, n_mod))
    corrected = np.zeros((nb, n_mod), dtype=bool)
    pwo4_index, names = [], {}
    for i in range(gc.shape[1]):
        mod = hycal.module_by_id(1001 + i)
        if mod is None:
            continue
        mi = mod.index
        pwo4_index.append(mi)
        names[mi] = mod.name
        c1, c2 = gc[:, i, 1], gc[:, i, 2]
        v1 = (c1 > 0) & (c1 != 1.0)
        v2 = (c2 > 0) & (c2 != 1.0)
        if ref_mode == 'ref-run':
            a1 = np.where(c1 > 0, c1, 1.0)
            a2 = np.where(c2 > 0, c2, 1.0)
            factors[:, mi] = 0.5 * (a1 + a2)
            corrected[:, mi] = v1 | v2
            continue
        ok1, ok2 = v1 & v1[0], v2 & v2[0]      # valid in batch b AND batch 0
        both, one1, one2 = ok1 & ok2, ok1 & ~ok2, ok2 & ~ok1
        f = np.ones(nb)
        with np.errstate(divide='ignore', invalid='ignore'):  # masked below
            f[both] = ((c1 + c2) / (c1[0] + c2[0]))[both]
            f[one1] = (c1 / c1[0])[one1]
            f[one2] = (c2 / c2[0])[one2]
        factors[:, mi] = f
        corrected[:, mi] = both | one1 | one2
    return LmsTable(path=path, ref_mode=ref_mode,
                    ref_run=int(np.asarray(g["ref_run"])[order][0]),
                    event_start=start.astype(np.int64).tolist(),
                    event_end=np.asarray(g["event_num_end"])[order].astype(np.int64).tolist(),
                    unix_time=np.asarray(g["unix_time"])[order].astype(np.int64).tolist(),
                    n_lms_events=np.asarray(g["n_lms_events"])[order].astype(np.int64).tolist(),
                    factors=factors, corrected=corrected,
                    pwo4_index=pwo4_index, names=names)


def cluster_lms_feff(cluster, fvec):
    """Energy-weighted effective LMS factor of a cluster:
    Σ E_h / Σ (E_h / f_h).  E_c / f_eff ≈ the no-LMS cluster energy."""
    num = den = 0.0
    for h in cluster.hits:
        num += h.energy
        den += h.energy / fvec[h.index]
    return num / den if den > 0 else 1.0


def lms_batch_stats(lms, b):
    """Factor statistics (median, p05, p95, min, max) over the PbWO4
    modules that were actually corrected in batch b."""
    idx = np.asarray(lms.pwo4_index, dtype=int)
    if len(idx) == 0:
        return (1.0,) * 5
    f = lms.factors[b, idx][lms.corrected[b, idx]]
    if len(f) == 0:
        return (1.0,) * 5
    return (float(np.median(f)), float(np.percentile(f, 5)),
            float(np.percentile(f, 95)), float(f.min()), float(f.max()))


def write_lms_tsv(path, lms, n_ev_per_batch, touched):
    with open(path, "w") as f:
        f.write("# benchmark_hycal_timing LMS factor summary\n")
        f.write(f"# lms_corr={lms.path} ref_mode={lms.ref_mode} "
                f"ref_run={lms.ref_run} n_batches={lms.n_batches} "
                f"n_lms_events_batch0={lms.n_lms_events[0]}\n")
        f.write(f"# event_range={lms.event_start[0]}..{lms.event_end[-1]} "
                f"batches_touched={','.join(map(str, touched))} "
                f"n_uncorrected_module_batches="
                f"{sum(lms.n_uncorrected(b) for b in touched)}\n")
        f.write("batch_id\tevent_num_start\tevent_num_end\tunix_time"
                "\tminutes_since_batch0\tn_lms_events\tn_events_this_job"
                "\tf_median\tf_p05\tf_p95\tf_min\tf_max\n")
        for b in range(lms.n_batches):
            med, p05, p95, fmin, fmax = lms_batch_stats(lms, b)
            mins = (lms.unix_time[b] - lms.unix_time[0]) / 60.0
            f.write(f"{b}\t{lms.event_start[b]}\t{lms.event_end[b]}"
                    f"\t{lms.unix_time[b]}\t{mins:.2f}\t{lms.n_lms_events[b]}"
                    f"\t{n_ev_per_batch[b]}\t{med:.5f}\t{p05:.5f}\t{p95:.5f}"
                    f"\t{fmin:.5f}\t{fmax:.5f}\n")


def plot_lms_factors(path, lms, n_ev_per_batch):
    st = np.array([lms_batch_stats(lms, b) for b in range(lms.n_batches)])
    mins = (np.asarray(lms.unix_time, dtype=float) - lms.unix_time[0]) / 60.0
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.fill_between(mins, st[:, 1], st[:, 2], color='#1f77b4', alpha=0.25,
                    label='5–95 % modules')
    ax.plot(mins, st[:, 0], 'o-', color='#1f77b4', label='median')
    ax.plot(mins, st[:, 3], '--', color='#888', lw=1, label='min / max')
    ax.plot(mins, st[:, 4], '--', color='#888', lw=1)
    step = (mins[1] - mins[0]) if len(mins) > 1 else 1.0
    for b in range(lms.n_batches):
        if n_ev_per_batch[b] > 0:
            ax.axvspan(mins[b] - 0.5 * step, mins[b] + 0.5 * step,
                       color='#ff7f0e', alpha=0.15, lw=0,
                       label='batches in this job' if b == min(
                           i for i, n in enumerate(n_ev_per_batch) if n > 0)
                       else None)
    ax.set_xlabel('minutes since batch 0')
    ax.set_ylabel('LMS factor')
    ax.set_title(f"LMS factor per batch ({lms.ref_mode})\n"
                 f"{os.path.basename(lms.path)}", fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(loc='best', fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def write_window_tsv(path, rows_by_path, args, eprime_mean, sigma_glob,
                     mu_global, kin_info):
    cols = [("path", "{}"), ("module", "{}"), ("module_id", "{}"),
            ("row", "{}"), ("col", "{}"), ("r_mm", "{:.2f}"),
            ("theta_deg", "{:.4f}"), ("eprime_MeV", "{:.3f}"),
            ("emoller_MeV", "{:.3f}"), ("prefit_lo_MeV", "{:.2f}"),
            ("prefit_hi_MeV", "{:.2f}"), ("n_total", "{}"),
            ("n_prefit", "{}"), ("mu_M_MeV", "{:.2f}"),
            ("sigma_M_MeV", "{:.2f}"), ("sigma_source", "{}"),
            ("win_lo_MeV", "{:.2f}"), ("win_hi_MeV", "{:.2f}"),
            ("n_in_window", "{}"), ("frac_in_window_lms", "{:.4f}"),
            ("frac_in_window_nolms", "{:.4f}"),
            ("median_in_window_MeV", "{:.3f}"), ("target_MeV", "{:.3f}"),
            ("gain", "{:.5f}"), ("flag_mu_off", "{}"), ("used", "{}"),
            ("moller_gap_nsigma", "{:.2f}"), ("exclude_reason", "{}")]
    keys = ["path", "module", "module_id", "row", "col", "r_mm", "theta_deg",
            "eprime", "emoller", "prefit_lo", "prefit_hi", "n_total",
            "n_prefit", "mu", "sigma", "sigma_source", "win_lo", "win_hi",
            "n_in_window", "frac_lms", "frac_nolms", "median_in", "target",
            "gain", "flag_mu_off", "used", "moller_gap_nsigma",
            "exclude_reason"]
    with open(path, "w") as f:
        f.write("# benchmark_hycal_timing elastic-window per-module table\n")
        f.write(f"# window_nsigma={args.window_nsigma} "
                f"prefit_range={kin_info['prefit']} "
                f"min_fit_events={MIN_FIT_EVENTS} "
                f"moller_gap_nsigma={MOLLER_GAP_NSIGMA} "
                f"min_per_module(in-window)={args.min_per_module}\n")
        f.write(f"# target={kin_info['target']} kinematics={kin_info['label']} "
                f"mass_MeV={kin_info['mass']:.3f} Ebeam_MeV={kin_info['ebeam']:.2f}\n")
        f.write("# eprime_mean_MeV " + " ".join(
            f"{lab}={em:.2f}(sigma_glob={sg:.1f} mu_global={mg:.2f} "
            f"offset={(mg / em - 1) * 100:+.2f}%)"
            for lab, em, sg, mg in [(k, eprime_mean[k], sigma_glob[k],
                                     mu_global[k]) for k in rows_by_path])
            + "\n")
        f.write("\t".join(c for c, _ in cols) + "\n")
        for label, rows in rows_by_path.items():
            for r in rows:
                r = dict(r, path=label, module_id=r['module'])
                r['module'] = _mod_name(r['module_id'])
                f.write("\t".join(fmt.format(r[k])
                                  for (c, fmt), k in zip(cols, keys)) + "\n")


def _mod_name(sid):
    # PrimEx IDs: 1-900 → "G<n>", 1001-2152 → "W<n>"
    if 1 <= sid <= 900:
        return f"G{sid}"
    if 1001 <= sid <= 2152:
        return f"W{sid - 1000}"
    return str(sid)


def calib_tag_energy(path):
    """Beam energy (MeV) implied by a calibration-file basename, or None."""
    base = os.path.basename(path or "").lower()
    for tag, e in CALIB_TAGS:
        if tag in base:
            return e
    return None


def calib_candidates(ebeam):
    """Files in the database calibration dir whose tag matches `ebeam`."""
    db = os.environ.get("PRAD2_DATABASE_DIR") or _repo_database_dir()
    cdir = os.path.join(db, "calibration")
    try:
        names = sorted(os.listdir(cdir))
    except OSError:
        return cdir, []
    out = []
    for n in names:
        e = calib_tag_energy(n)
        if n.endswith(".json") and e is not None and \
                abs(e / ebeam - 1.0) < CALIB_MATCH_TOL:
            out.append(os.path.join(cdir, n))
    return cdir, out


def check_calib_file(calib, explicit, ebeam):
    """Compare the effective calibration file's energy tag with Ebeam.
    Returns 'match' / 'mismatch' / 'unknown-tag'.  A runinfo-default file
    that does not match is fatal (exit 2); an explicit --hc-calib-file
    only gets a warning banner."""
    e_tag = calib_tag_energy(calib)
    if e_tag is None:
        status = 'unknown-tag'
    else:
        status = ('match' if abs(e_tag / ebeam - 1.0) < CALIB_MATCH_TOL
                  else 'mismatch')
    if status == 'match':
        return status
    why = ("has no recognised energy tag (0p7/2p1/2p2/3p5/x17)"
           if status == 'unknown-tag'
           else f"is tagged ~{e_tag:.0f} MeV")
    cdir, cands = calib_candidates(ebeam)
    cand_txt = ("\n".join(f"            {c}" for c in cands)
                if cands else f"            (none found in {cdir})")
    if explicit:
        bar = "!" * 72
        sys.stderr.write(
            f"{bar}\n"
            f"!! WARNING: calibration file does not match the beam energy\n"
            f"!!   Ebeam         : {ebeam:.1f} MeV\n"
            f"!!   --hc-calib-file: {calib}\n"
            f"!!   this file {why}\n"
            f"!!   Continuing because the file was given explicitly.\n"
            f"!!   Matching candidates in {cdir}:\n"
            + "\n".join(f"!!     {c}" for c in cands or ["(none)"]) + "\n"
            f"{bar}\n")
        sys.stderr.flush()
        return status
    sys.stderr.write(
        f"[ERROR] runinfo calibration file does not match the beam energy.\n"
        f"        Ebeam : {ebeam:.1f} MeV\n"
        f"        file  : {calib}\n"
        f"        this file {why}.\n"
        f"        Pass --hc-calib-file with one of the matching candidates:\n"
        f"{cand_txt}\n")
    sys.exit(2)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    C.add_common_args(ap)
    ap.add_argument("--window", type=float, default=8.0,
                    help="seed_time_window (ns) for the new method (default 8).")
    ap.add_argument("--pre-window", type=float, nargs=2, default=(130.0, 200.0),
                    metavar=("LO", "HI"),
                    help="Pre-window applied to peaks before clustering for "
                         "both paths (default 130 200 ns — runinfo for 24386).")
    ap.add_argument("--signal-min", type=float, default=None,
                    help="Single-cluster energy cut in MeV (default "
                         "500·kE, kE = Ebeam/2239.51).")
    ap.add_argument("--bin-width", type=float, default=None,
                    help="Energy histogram bin width in MeV (default "
                         "20·kS, kS = sqrt(Ebeam/2239.51)).")
    ap.add_argument("--e-range", type=float, nargs=2, default=None,
                    metavar=("LO", "HI"),
                    help="Energy histogram range in MeV (default "
                         "500·kE 3500·kE).")
    ap.add_argument("--calibrate", action="store_true",
                    help="After the main loop, free per-seed-module gain "
                         "constants on the inner-ring high-statistics sample "
                         "(target peak = Ebeam by default), then refit the "
                         "Gaussian to compare σ_E between methods under "
                         "matched per-module calibration.")
    ap.add_argument("--target-peak", type=float, default=None,
                    help="Target μ for the median / cb-peak gain refinement "
                         "(default Ebeam in MeV; ignored by elastic-window).")
    ap.add_argument("--inner-radius", type=float, default=200.0,
                    help="Restrict per-module gain fits to seeds within this "
                         "radius (mm) from beam (default 200).")
    ap.add_argument("--min-seed-energy", type=float, default=None,
                    help="Skip events whose seed energy is below this when "
                         "computing per-module gains (default 200·kE MeV).")
    ap.add_argument("--min-per-module", type=int, default=30,
                    help="Minimum events per module before a gain is fit "
                         "(default 30).")
    ap.add_argument("--exclude-rowcol", type=int, nargs=4, default=None,
                    metavar=("R_LO", "R_HI", "C_LO", "C_HI"),
                    help="Drop seed modules with row in [R_LO, R_HI] AND "
                         "column in [C_LO, C_HI] from the per-module "
                         "calibration sample (e.g. --exclude-rowcol 16 19 "
                         "16 19 removes the 4x4 PbWO4 block surrounding the "
                         "beam hole, where leakage biases the per-module mean).")
    ap.add_argument("--cal-method",
                    choices=['elastic-window', 'cb-peak', 'median'],
                    default='elastic-window',
                    help="Per-module gain refinement: 'elastic-window' "
                         "(default; CB pre-fit, fixed ±--window-nsigma window, "
                         "in-window median -> elastic E'(theta_M) of --target; "
                         "ignores --target-peak), 'median' (legacy, biased by the "
                         "radiative tail) or 'cb-peak' (Crystal-Ball fitted "
                         "peak per module, robust against the tail; falls "
                         "back to a narrow-window Gaussian and then the "
                         "histogram argmax for low-statistics modules).  "
                         "median / cb-peak use --target-peak, 2 iterations.")
    ap.add_argument("--target", choices=['H2', 'C', 'Ta'], default=None,
                    help="Elastic target for elastic-window (H2 = e-p, C / Ta "
                         "= e-nucleus recoil kinematics).  Default: from the "
                         "per-run map RUN_TARGET; an unmapped run needs this "
                         "flag with elastic-window (exit 2 otherwise).")
    ap.add_argument("--window-nsigma", type=float, default=2.0,
                    help="Half-width of the elastic-window refinement window "
                         "in units of the per-module CB sigma (default 2).")
    ap.add_argument("--lms-corr", default='auto', metavar="{auto,off,PATH}",
                    help="LMS gain correction from replay_gainCorr output: "
                         "'auto' = $PRAD2_DATABASE_DIR/gain_factor/"
                         "gain_correction/prad_<run>_gain_corr.root, 'off' = "
                         "no correction (uproot not needed), or an explicit "
                         "file.  A missing file is an error (exit 2).")
    ap.add_argument("--lms-ref", choices=['first-batch', 'ref-run'],
                    default='first-batch',
                    help="LMS normalisation: 'first-batch' = A_b/A_0 "
                         "(in-run drift only, default) or 'ref-run' = the "
                         "stored scale, as the C++ replay.")
    ap.add_argument("--write-cluster-tsv", action="store_true",
                    help="Write <out_path>_clusters.tsv with one row per "
                         "accepted single-cluster event (columns: event_num, "
                         "config_name, seed_module_id, seed_x_mm, seed_y_mm, "
                         "seed_row, seed_col, seed_energy_MeV, "
                         "cluster_E_raw_MeV, cluster_E_cal_MeV).  Requires "
                         "--calibrate (the event lists are only populated in "
                         "that mode).")
    args = ap.parse_args(argv)

    if args.write_cluster_tsv and not args.calibrate:
        sys.stderr.write(
            "[ERROR] --write-cluster-tsv requires --calibrate to be set.\n"
            "        Re-run with both flags: --calibrate --write-cluster-tsv\n"
        )
        sys.exit(2)

    run_for_lms = (args.run_num if args.run_num > 0
                   else C.extract_run_number(args.evio_path))
    target = args.target or RUN_TARGET.get(run_for_lms)
    if args.calibrate and args.cal_method == 'elastic-window' and target is None:
        sys.stderr.write(
            f"[ERROR] elastic-window needs the elastic target: run "
            f"{run_for_lms} is not in RUN_TARGET and --target was not given.\n"
            f"        Pass --target {{H2,C,Ta}} (or use --cal-method "
            f"cb-peak / median).\n")
        sys.exit(2)
    m_target = TARGET_MASS[target] if target else None
    kin_label = ('ep' if target == 'H2' else 'eA') if target else 'n/a'
    lms_path = resolve_lms_path(args.lms_corr, run_for_lms)
    if lms_path is not None:
        if args.lms_corr == 'auto' and run_for_lms <= 0:
            sys.stderr.write("[ERROR] --lms-corr auto: cannot determine the "
                             "run number; pass --run-num, a PATH, or "
                             "--lms-corr off.\n")
            sys.exit(2)
        if not os.path.isfile(lms_path):
            sys.stderr.write(
                f"[ERROR] LMS gain_corr file not found: {lms_path}\n"
                "        Produce it with replay_gainCorr on ALL splits of the run "
                "(plan §8), e.g.:\n"
                "          DB=/scratch/$USER/db_b1600; cp -R $REPO/database $DB; "
                "mkdir -p /scratch/$USER/lms_<run>\n"
                "          PRAD2_DATABASE_DIR=$DB $REPO/build/bin/"
                "prad2ana_replay_gainCorr -o /scratch/$USER/lms_<run> -b 1600 "
                "-r 24827 -j 4 $RAW/prad_<run>.evio.*\n"
                "          cp $DB/gain_factor/gain_correction/"
                "prad_<run>_gain_corr.root <the path above>\n"
                "        or re-run with --lms-corr off (or --lms-corr PATH).\n")
            sys.exit(2)
        _require_uproot()

    p = C.setup_pipeline(
        evio_path     = args.evio_path,
        max_events    = args.max_events,
        run_num       = args.run_num,
        gem_ped_file  = args.gem_ped_file,
        gem_cm_file   = args.gem_cm_file,
        hc_calib_file = args.hc_calib_file,
        daq_config    = args.daq_config,
        gem_map_file  = args.gem_map_file,
        hc_map_file   = args.hc_map_file,
    )
    # ---- beam-energy scaling of the default thresholds / shapes ------------
    ebeam = float(p.geo.Ebeam)
    kE = round(ebeam / E_REF, 6)
    kS = round(math.sqrt(ebeam / E_REF), 6)
    auto = dict(signal_min=args.signal_min is None,
                e_range=args.e_range is None,
                min_seed_energy=args.min_seed_energy is None,
                bin_width=args.bin_width is None,
                target_peak=args.target_peak is None)
    if auto['signal_min']:      args.signal_min = 500.0 * kE
    if auto['e_range']:         args.e_range = (500.0 * kE, 3500.0 * kE)
    if auto['min_seed_energy']: args.min_seed_energy = 200.0 * kE
    if auto['bin_width']:       args.bin_width = 20.0 * kS
    if auto['target_peak']:     args.target_peak = ebeam
    src = lambda k: 'auto' if auto[k] else 'cli'
    prefit = (PREFIT_RANGE[0] / kS, PREFIT_RANGE[1] / kS)
    print(f"[setup] Ebeam {ebeam:.2f} MeV  kE={kE:g}  kS={kS:g}  "
          f"hycal_z {float(p.geo.hycal_z):.1f} mm", flush=True)
    print(f"[setup] signal-min {args.signal_min:.1f} ({src('signal_min')})  "
          f"e-range {args.e_range[0]:.1f} {args.e_range[1]:.1f} "
          f"({src('e_range')})  min-seed-energy {args.min_seed_energy:.1f} "
          f"({src('min_seed_energy')})", flush=True)
    print(f"[setup] bin-width {args.bin_width:.3f} ({src('bin_width')})  "
          f"target-peak {args.target_peak:.1f} ({src('target_peak')})  "
          f"pre-fit fractions ({prefit[0]:.4f}, {prefit[1]:.4f})", flush=True)
    if target:
        print(f"[setup] target {target}  kinematics {kin_label}  "
              f"mass {m_target:.3f} MeV", flush=True)
    else:
        print("[setup] target unresolved (not needed for this --cal-method)",
              flush=True)

    # ---- calibration file vs beam energy ------------------------------------
    calib_file = args.hc_calib_file or str(p.geo.energy_calib_file)
    calib_match = check_calib_file(calib_file, bool(args.hc_calib_file), ebeam)
    print(f"[setup] calib file {calib_file}  ({calib_match})", flush=True)

    pre_lo, pre_hi = args.pre_window
    W = args.window
    print(f"[setup] pre-window  : [{pre_lo}, {pre_hi}] ns", flush=True)
    print(f"[setup] seed window : ±{W} ns (new method only)", flush=True)
    print(f"[setup] signal cut  : single cluster E > {args.signal_min:.0f} MeV",
          flush=True)
    lms = load_lms_gain_corr(lms_path, p.hycal, args.lms_ref)
    print(f"[setup] LMS corr    : {lms.path}  ({lms.n_batches} batch(es), "
          f"ref {args.lms_ref})", flush=True)
    lms_tag = "LMS: off" if lms_path is None else "LMS: on"
    n_ev_per_batch = [0] * lms.n_batches

    # Two clusterers sharing the geometry; differ only in seed_time_window.
    base_cfg = p.hc_clusterer.get_config()
    def _make_cluster(seed_window):
        cfg = det.HyCalClusterConfig()
        for f in ("min_module_energy", "min_center_energy", "min_cluster_energy",
                  "min_cluster_size", "corner_conn", "split_iter",
                  "least_split", "log_weight_thres"):
            setattr(cfg, f, getattr(base_cfg, f))
        cfg.seed_time_window = seed_window
        cl = det.HyCalCluster(p.hycal)
        cl.set_config(cfg)
        return cl
    cl_legacy = _make_cluster(-1.0)
    cl_new    = _make_cluster(W)

    # Bulk: single-cluster energies for the bulk Gaussian fit.
    # Per-event: (seed_id, x, y, row, col, seed_E, cluster_E, f_eff) —
    # populated only when --calibrate is set, so the memory cost is opt-in.
    energies_legacy, energies_new = [], []
    events_legacy,   events_new   = [], []
    n_legacy = n_new = n_phys = n_kept = n_read = n_files_open = 0

    ch = dec.EvChannel()
    ch.set_config(p.cfg)

    t0 = time.monotonic()
    try:
        for fpath in p.evio_files:
            if ch.open_auto(fpath) != dec.Status.success:
                print(f"[WARN] skip (cannot open): {fpath}", flush=True)
                continue
            n_files_open += 1
            print(f"[file {n_files_open}/{len(p.evio_files)}] {fpath}",
                  flush=True)

            done = False
            while ch.read() == dec.Status.success:
                n_read += 1
                if not ch.scan():
                    continue
                if ch.get_event_type() != dec.EventType.Physics:
                    continue

                for i in range(ch.get_n_events()):
                    decoded = ch.decode_event(i, with_ssp=False)
                    if not decoded["ok"]:
                        continue
                    n_phys += 1
                    fadc_evt = decoded["event"]
                    if fadc_evt.info.trigger_bits != C.PHYSICS_TRIGGER_BITS:
                        if args.max_events > 0 and n_phys >= args.max_events:
                            done = True; break
                        continue
                    n_kept += 1
                    b = lms.batch_of(int(fadc_evt.info.event_number))
                    fvec = lms.factors_list[b]
                    n_ev_per_batch[b] += 1

                    cl_legacy.clear()
                    cl_new.clear()

                    # One waveform pass; both clusterers consume the same peaks.
                    for ri in range(fadc_evt.nrocs):
                        roc = fadc_evt.roc(ri)
                        if not roc.present:
                            continue
                        crate = p.crate_map.get(roc.tag)
                        if crate is None:
                            continue
                        for s in roc.present_slots():
                            slot = roc.slot(s)
                            for c in slot.present_channels():
                                mod = p.hycal.module_by_daq(crate, s, c)
                                if mod is None or not mod.is_hycal():
                                    continue
                                cd = slot.channel(c)
                                if cd.nsamples <= 0:
                                    continue
                                _, _, peaks = p.wave_ana.analyze(cd.samples)
                                f = fvec[mod.index]

                                best_int  = -1.0
                                best_pk   = None
                                for pk in peaks:
                                    if pk.time <= pre_lo or pk.time >= pre_hi:
                                        continue
                                    cl_new.add_hit(mod.index,
                                                   mod.energize(pk.integral * f),
                                                   float(pk.time))
                                    if pk.integral > best_int:
                                        best_int = pk.integral
                                        best_pk  = pk
                                if best_pk is not None:
                                    cl_legacy.add_hit(mod.index,
                                                      mod.energize(best_pk.integral * f),
                                                      float(best_pk.time))

                    cl_legacy.form_clusters()
                    cl_new.form_clusters()
                    # reconstruct_matched gives us seed energy + center
                    # module without a second lookup; needed only by the
                    # --calibrate branch but cheap enough to use always.
                    rec_l = cl_legacy.reconstruct_matched()
                    rec_n = cl_new.reconstruct_matched()
                    big_l = [r for r in rec_l if r.hit.energy > args.signal_min]
                    big_n = [r for r in rec_n if r.hit.energy > args.signal_min]
                    if len(big_l) == 1:
                        r = big_l[0]
                        n_legacy += 1
                        energies_legacy.append(r.hit.energy)
                        if args.calibrate:
                            sm = p.hycal.module_by_id(r.hit.center_id)
                            events_legacy.append((r.hit.center_id, sm.x, sm.y,
                                                  sm.row, sm.column,
                                                  r.cluster.center.energy,
                                                  r.hit.energy,
                                                  cluster_lms_feff(r.cluster, fvec)))
                    if len(big_n) == 1:
                        r = big_n[0]
                        n_new += 1
                        energies_new.append(r.hit.energy)
                        if args.calibrate:
                            sm = p.hycal.module_by_id(r.hit.center_id)
                            events_new.append((r.hit.center_id, sm.x, sm.y,
                                               sm.row, sm.column,
                                               r.cluster.center.energy,
                                               r.hit.energy,
                                               cluster_lms_feff(r.cluster, fvec)))

                    if args.max_events > 0 and n_phys >= args.max_events:
                        done = True; break

                if done:
                    break
                if n_phys > 0 and n_phys % 5000 == 0:
                    elapsed = time.monotonic() - t0
                    rate = n_phys / max(elapsed, 1e-3)
                    print(f"[progress] {n_phys} physics  legacy={n_legacy} "
                          f"new={n_new}  ({rate:.1f} ev/s)", flush=True)

            ch.close()
            if done:
                break
    finally:
        elapsed = time.monotonic() - t0

    # ---- summary + Gaussian fits ------------------------------------------
    e_lo, e_hi = args.e_range
    mu_l, sig_l, A_l, nf_l, edges, counts_l = fit_peak(
        energies_legacy, e_lo, e_hi, args.bin_width, kS)
    mu_n, sig_n, A_n, nf_n, _,     counts_n = fit_peak(
        energies_new,    e_lo, e_hi, args.bin_width, kS)
    centres = 0.5 * (edges[:-1] + edges[1:])

    rate = n_phys / max(elapsed, 1e-3)
    print()
    print("=" * 72)
    print("benchmark_hycal_timing — legacy vs. seed-time-gated clustering")
    print("=" * 72)
    print(f"  files {n_files_open}/{len(p.evio_files)},  "
          f"records {n_read},  physics {n_phys},  triggered {n_kept}")
    print(f"  elapsed {elapsed:.1f} s  ({rate:.1f} ev/s)")
    print(f"  selection: single cluster E > {args.signal_min:.0f} MeV")
    print()
    print(f"               | legacy        | gated (W={W:g} ns)")
    print(f"  events kept  | {n_legacy:>13d} | {n_new:>13d}    "
          f"({(n_new-n_legacy)/max(n_legacy,1)*100:+.1f} % vs. legacy)")
    if HAVE_SCIPY:
        print(f"  Gauss μ (MeV)| {mu_l:>13.1f} | {mu_n:>13.1f}")
        print(f"  Gauss σ (MeV)| {sig_l:>13.1f} | {sig_n:>13.1f}")
        print(f"  N in fit win | {nf_l:>13d} | {nf_n:>13d}")
    else:
        print("  (scipy not available — fit skipped)")

    # ---- LMS summary (stdout) ----------------------------------------------
    touched = [bi for bi in range(lms.n_batches) if n_ev_per_batch[bi] > 0]
    print()
    print(f"  LMS correction  : {lms.path}  (ref {args.lms_ref})")
    if lms_path is not None:
        fr = [lms_batch_stats(lms, bi) for bi in touched]
        print(f"    batches touched : "
              + ", ".join(f"{bi}({n_ev_per_batch[bi]} ev)" for bi in touched))
        if fr:
            print(f"    factor range    : {min(x[3] for x in fr):.4f} .. "
                  f"{max(x[4] for x in fr):.4f}  (touched batches, corrected modules)")
        unc = lms.uncorrected_names(touched)
        print(f"    uncorrected     : {len(unc)} module(s)"
              + (f"  {unc[:10]}" if unc else ""))

    # ---- Optional: per-seed-module calibration on the inner ring ----------
    cal = None
    rows_by_path = {}
    eprime_by_path, sigma_glob_by_path, mu_global_by_path = {}, {}, {}
    target_ref = args.target_peak
    if args.calibrate:
        ew = args.cal_method == 'elastic-window'
        print()
        print(f"  --calibrate mode  : target μ = "
              + ("E'(θ_M) per module (elastic-window)" if ew
                 else f"{args.target_peak:.1f} MeV")
              + f", inner r < {args.inner_radius:.0f} mm,")
        print(f"                      seed E ≥ {args.min_seed_energy:.0f} MeV, "
              f"≥ {args.min_per_module} events/module")

        excl = tuple(args.exclude_rowcol) if args.exclude_rowcol else None
        if excl:
            print(f"                      exclude rows {excl[0]}..{excl[1]} "
                  f"AND cols {excl[2]}..{excl[3]}")
        print(f"                      per-module reference: {args.cal_method}")
        if ew:
            ew_kw = dict(ebeam=float(p.geo.Ebeam), hycal_z=float(p.geo.hycal_z),
                         m_target=m_target, escale=kS, kE=kE,
                         inner_r=args.inner_radius,
                         min_seed_E=args.min_seed_energy,
                         min_per_mod=args.min_per_module, exclude_rowcol=excl,
                         nsigma=args.window_nsigma, bin_width=args.bin_width)
            print(f"                      Ebeam {ew_kw['ebeam']:.2f} MeV, "
                  f"hycal_z {ew_kw['hycal_z']:.1f} mm, window ±"
                  f"{args.window_nsigma:g}σ, target {target} ({kin_label}, "
                  f"M {m_target:.1f} MeV)")
            gains_l, raw_l, corr_l, nmod_l, rows_l, sg_l, ep_l, mg_l = \
                fit_elastic_window_gains(events_legacy, **ew_kw)
            gains_n, raw_n, corr_n, nmod_n, rows_n, sg_n, ep_n, mg_n = \
                fit_elastic_window_gains(events_new, **ew_kw)
            rows_by_path = {'legacy': rows_l, 'new': rows_n}
            eprime_by_path = {'legacy': ep_l, 'new': ep_n}
            sigma_glob_by_path = {'legacy': sg_l, 'new': sg_n}
            mu_global_by_path = {'legacy': mg_l, 'new': mg_n}
            fin = [v for v in (ep_l, ep_n) if not math.isnan(v)]
            if fin:
                target_ref = float(np.mean(fin))
        else:
            gains_l, raw_l, corr_l, nmod_l = fit_per_module_gains(
                events_legacy,
                target_peak=args.target_peak, inner_r=args.inner_radius,
                min_seed_E=args.min_seed_energy, min_per_mod=args.min_per_module,
                exclude_rowcol=excl, method=args.cal_method,
                bin_width=args.bin_width, escale=kS)
            gains_n, raw_n, corr_n, nmod_n = fit_per_module_gains(
                events_new,
                target_peak=args.target_peak, inner_r=args.inner_radius,
                min_seed_E=args.min_seed_energy, min_per_mod=args.min_per_module,
                exclude_rowcol=excl, method=args.cal_method,
                bin_width=args.bin_width, escale=kS)

        # Fit window wider than the post-cal peak so a poorly-calibrated
        # pre-cal histogram still falls inside.
        cal_lo, cal_hi = target_ref - 800 * kE, target_ref + 800 * kE
        bw = args.bin_width
        mu_l_pre,  sig_l_pre,  _, _, _,         counts_l_pre  = fit_peak(raw_l,  cal_lo, cal_hi, bw, kS)
        mu_n_pre,  sig_n_pre,  _, _, _,         counts_n_pre  = fit_peak(raw_n,  cal_lo, cal_hi, bw, kS)
        mu_l_post, sig_l_post, A_l_post, _, edges_cal, counts_l_post = fit_peak(corr_l, cal_lo, cal_hi, bw, kS)
        mu_n_post, sig_n_post, A_n_post, _, _,         counts_n_post = fit_peak(corr_n, cal_lo, cal_hi, bw, kS)
        cal_centres = 0.5 * (edges_cal[:-1] + edges_cal[1:])

        cb_l = fit_peak_cb(corr_l, cal_lo, cal_hi, bw, mu0_hint=mu_l_post, escale=kS)
        cb_n = fit_peak_cb(corr_n, cal_lo, cal_hi, bw, mu0_hint=mu_n_post, escale=kS)
        mu_cb_l_post, sig_cb_l_post, alpha_l_post, n_l_post, A_cb_l_post, *_ = cb_l
        mu_cb_n_post, sig_cb_n_post, alpha_n_post, n_n_post, A_cb_n_post, *_ = cb_n

        # ALSO fit CB on pre-cal spectra (symmetric with post for comparison)
        cb_l_pre = fit_peak_cb(raw_l, cal_lo, cal_hi, bw, mu0_hint=mu_l_pre, escale=kS)
        cb_n_pre = fit_peak_cb(raw_n, cal_lo, cal_hi, bw, mu0_hint=mu_n_pre, escale=kS)
        mu_cb_l_pre, sig_cb_l_pre, alpha_l_pre, n_l_pre, A_cb_l_pre, *_ = cb_l_pre
        mu_cb_n_pre, sig_cb_n_pre, alpha_n_pre, n_n_pre, A_cb_n_pre, *_ = cb_n_pre

        def _ratio(s, m):
            return s / m if (m and m > 0 and not math.isnan(m)) else float('nan')

        cal = dict(
            n_modules_l=nmod_l, n_modules_n=nmod_n,
            n_events_l=len(raw_l), n_events_n=len(raw_n),
            mu_l_pre=mu_l_pre, sig_l_pre=sig_l_pre,
            mu_n_pre=mu_n_pre, sig_n_pre=sig_n_pre,
            mu_l_post=mu_l_post, sig_l_post=sig_l_post,
            mu_n_post=mu_n_post, sig_n_post=sig_n_post,
            A_l_post=A_l_post, A_n_post=A_n_post,   # Gaussian amplitudes (plot only)
            # Post-cal CB (renamed from mu_cb_l / sig_cb_l etc.)
            mu_cb_l_post=mu_cb_l_post, sig_cb_l_post=sig_cb_l_post,
            alpha_l_post=alpha_l_post, n_l_post=n_l_post, A_cb_l_post=A_cb_l_post,
            mu_cb_n_post=mu_cb_n_post, sig_cb_n_post=sig_cb_n_post,
            alpha_n_post=alpha_n_post, n_n_post=n_n_post, A_cb_n_post=A_cb_n_post,
            # NEW: Pre-cal CB
            mu_cb_l_pre=mu_cb_l_pre, sig_cb_l_pre=sig_cb_l_pre,
            alpha_l_pre=alpha_l_pre, n_l_pre=n_l_pre, A_cb_l_pre=A_cb_l_pre,
            mu_cb_n_pre=mu_cb_n_pre, sig_cb_n_pre=sig_cb_n_pre,
            alpha_n_pre=alpha_n_pre, n_n_pre=n_n_pre, A_cb_n_pre=A_cb_n_pre,
            edges=edges_cal, centres=cal_centres,
            counts_l_pre=counts_l_pre, counts_n_pre=counts_n_pre,
            counts_l_post=counts_l_post, counts_n_post=counts_n_post,
        )

        sigE_l_g,      sigE_n_g      = _ratio(sig_l_post, mu_l_post),         _ratio(sig_n_post, mu_n_post)
        sigE_l_cb_post, sigE_n_cb_post = _ratio(sig_cb_l_post, mu_cb_l_post), _ratio(sig_cb_n_post, mu_cb_n_post)
        sigE_l_cb_pre,  sigE_n_cb_pre  = _ratio(sig_cb_l_pre,  mu_cb_l_pre),  _ratio(sig_cb_n_pre,  mu_cb_n_pre)
        # Keep old names as aliases for backward-compat downstream (if used)
        sigE_l_cb, sigE_n_cb = sigE_l_cb_post, sigE_n_cb_post
        print()
        print(f"                  | legacy        | gated (W={W:g} ns)")
        print(f"  inner events    | {len(raw_l):>13d} | {len(raw_n):>13d}")
        print(f"  modules used    | {nmod_l:>13d} | {nmod_n:>13d}")
        print(f"  μ before  (MeV) | {mu_l_pre:>13.1f} | {mu_n_pre:>13.1f}")
        print(f"  σ before  (MeV) | {sig_l_pre:>13.1f} | {sig_n_pre:>13.1f}")
        print(f"  μ  after  (MeV) | {mu_l_post:>13.1f} | {mu_n_post:>13.1f}  "
              f"(target {target_ref:.1f})")
        print(f"  σ  Gaussian     | {sig_l_post:>13.1f} | {sig_n_post:>13.1f}")
        print(f"  σ_E / E (G)     | {sigE_l_g*100:>12.2f}% | {sigE_n_g*100:>12.2f}%")
        print(f"  μ  CB post(MeV) | {mu_cb_l_post:>13.1f} | {mu_cb_n_post:>13.1f}")
        print(f"  σ  CB post(MeV) | {sig_cb_l_post:>13.1f} | {sig_cb_n_post:>13.1f}  "
              f"(α_l={alpha_l_post:.2f} n_l={n_l_post:.1f}  α_n={alpha_n_post:.2f} n_n={n_n_post:.1f})")
        print(f"  σ_E / E (CB po) | {sigE_l_cb_post*100:>12.2f}% | {sigE_n_cb_post*100:>12.2f}%  "
              f"(Δ = {(sigE_n_cb_post - sigE_l_cb_post)*100:+.2f} pp)")
        print(f"  μ  CB pre (MeV) | {mu_cb_l_pre:>13.1f} | {mu_cb_n_pre:>13.1f}")
        print(f"  σ  CB pre (MeV) | {sig_cb_l_pre:>13.1f} | {sig_cb_n_pre:>13.1f}  "
              f"(α_l={alpha_l_pre:.2f} n_l={n_l_pre:.1f}  α_n={alpha_n_pre:.2f} n_n={n_n_pre:.1f})")
        print(f"  σ_E / E (CB pr) | {sigE_l_cb_pre*100:>12.2f}% | {sigE_n_cb_pre*100:>12.2f}%  "
              f"(Δ = {(sigE_n_cb_pre - sigE_l_cb_pre)*100:+.2f} pp)")

        if ew:
            print()
            print(f"  elastic-window    : window ±{args.window_nsigma:g}σ, "
                  f"in-window min/module {args.min_per_module}, "
                  f"target mean {target_ref:.1f} MeV")
            for label, rws in rows_by_path.items():
                us = [r for r in rws if r['used']]
                n_in = sum(r['n_in_window'] for r in us)
                n_nl = sum(r['n_in_nolms'] for r in us)
                n_tot = sum(r['n_total'] for r in us)
                print(f"    {label:<6}: used {len(us)}/{len(rws)} modules, "
                      f"fits module/global "
                      f"{sum(r['sigma_source'] == 'module' for r in us)}/"
                      f"{sum(r['sigma_source'] == 'global' for r in us)}, "
                      f"flag_mu_off {sum(r['flag_mu_off'] for r in us)}, "
                      f"moller-excluded "
                      f"{sum(r['exclude_reason'] == 'moller_gap' for r in rws)}, "
                      f"σ_glob {sigma_glob_by_path[label]:.1f} MeV, "
                      f"μ_global {mu_global_by_path[label]:.1f} MeV "
                      f"(offset {(mu_global_by_path[label] / eprime_by_path[label] - 1) * 100:+.2f}% vs Ē')")
                if n_tot:
                    print(f"            pooled in-window fraction "
                          f"{n_in / n_tot:.3f} with LMS, "
                          f"{n_nl / n_tot:.3f} without")

    # ---- Optional: per-event cluster TSV ------------------------------------
    if args.write_cluster_tsv and cal is not None:
        # gains_l / gains_n come from the chosen --cal-method above.
        # events_X lists: (seed_id, x, y, row, col, seed_E, cluster_E, f_eff)
        clust_path = (args.out_path[:-4] if args.out_path.endswith('.tsv')
                      else args.out_path) + "_clusters.tsv"
        n_rows_l = n_rows_n = 0
        with open(clust_path, "w") as cf:
            cf.write("event_num\tconfig_name\tseed_module_id\t"
                     "seed_x_mm\tseed_y_mm\tseed_row\tseed_col\t"
                     "seed_energy_MeV\tcluster_E_raw_MeV\tcluster_E_cal_MeV\n")
            for ev_idx, (sid, x, y, row, col, se, ce, *_) in enumerate(events_legacy):
                g = gains_l.get(sid, float('nan'))
                ce_cal = ce * g if not math.isnan(g) else float('nan')
                cf.write(f"{ev_idx}\tlegacy\t{_mod_name(sid)}\t"
                         f"{x:.3f}\t{y:.3f}\t{row}\t{col}\t"
                         f"{se:.3f}\t{ce:.3f}\t{ce_cal:.3f}\n")
                n_rows_l += 1
            for ev_idx, (sid, x, y, row, col, se, ce, *_) in enumerate(events_new):
                g = gains_n.get(sid, float('nan'))
                ce_cal = ce * g if not math.isnan(g) else float('nan')
                cf.write(f"{ev_idx}\tgated\t{_mod_name(sid)}\t"
                         f"{x:.3f}\t{y:.3f}\t{row}\t{col}\t"
                         f"{se:.3f}\t{ce:.3f}\t{ce_cal:.3f}\n")
                n_rows_n += 1
        print(f"\n  wrote: {clust_path}  "
              f"(N rows: {n_rows_l} legacy + {n_rows_n} gated)")

    # ---- TSV output --------------------------------------------------------
    tsv_path = args.out_path if args.out_path.endswith('.tsv') \
               else args.out_path + ".tsv"
    with open(tsv_path, "w") as f:
        f.write("# benchmark_hycal_timing\n")
        f.write(f"# evio_path={args.evio_path} max_events={args.max_events}\n")
        f.write(f"# pre_window=[{pre_lo},{pre_hi}] ns  seed_time_window={W} ns\n")
        f.write(f"# selection: single cluster E > {args.signal_min} MeV\n")
        f.write(f"# n_phys={n_phys} n_kept={n_kept} elapsed_s={elapsed:.1f}\n")
        f.write(f"# lms_corr={lms.path} lms_ref={args.lms_ref} "
                f"batches_touched={','.join(map(str, touched))} n_uncorrected_module_batches="
                f"{sum(lms.n_uncorrected(bi) for bi in touched)}\n")
        f.write(f"# resolved Ebeam={ebeam:.2f} kE={kE:g} kS={kS:g} "
                f"signal_min={args.signal_min:.2f}({src('signal_min')}) "
                f"e_range={args.e_range[0]:.2f},{args.e_range[1]:.2f}"
                f"({src('e_range')}) "
                f"min_seed_energy={args.min_seed_energy:.2f}"
                f"({src('min_seed_energy')}) "
                f"bin_width={args.bin_width:.4f}({src('bin_width')}) "
                f"target_peak={args.target_peak:.2f}({src('target_peak')}) "
                f"prefit_fractions=({prefit[0]:.4f},{prefit[1]:.4f})\n")
        f.write(f"# kinematics target={target} label={kin_label} "
                f"mass_MeV={(m_target if m_target else float('nan')):.3f} "
                f"calib_file={calib_file} calib_match={calib_match}\n")
        f.write(f"# cal_method={args.cal_method} "
                f"target_mean_MeV={target_ref:.2f} "
                f"window_nsigma={args.window_nsigma}\n")
        if mu_global_by_path:
            f.write("# window_mu_global " + " ".join(
                f"{lab}={mg:.2f}MeV(offset={(mg / eprime_by_path[lab] - 1) * 100:+.2f}%)"
                for lab, mg in mu_global_by_path.items()) + "\n")
        f.write("\n")
        f.write("path\tn_events\tgauss_mu_MeV\tgauss_sigma_MeV\tn_in_fit\n")
        f.write(f"legacy\t{n_legacy}\t{mu_l:.2f}\t{sig_l:.2f}\t{nf_l}\n")
        f.write(f"new\t{n_new}\t{mu_n:.2f}\t{sig_n:.2f}\t{nf_n}\n")
        f.write("\n")
        f.write("# energy histogram bin centres (MeV) and counts\n")
        f.write("e_centre\tcount_legacy\tcount_new\n")
        for c, l, n in zip(centres, counts_l, counts_n):
            f.write(f"{c:.1f}\t{l}\t{n}\n")
        if cal is not None:
            f.write("\n")
            f.write(f"# --calibrate inner ring (r<{args.inner_radius:.0f} mm,"
                    f" seed E≥{args.min_seed_energy:.0f} MeV,"
                    f" min/mod={args.min_per_module},"
                    f" target μ={target_ref:.1f} MeV)\n")
            f.write("path\tn_modules\tn_inner_events\tmu_pre_MeV\tsigma_pre_MeV"
                    "\tmu_post_MeV\tsigma_post_MeV\tsigma_over_E_post"
                    "\tmu_cb_MeV\tsigma_cb_MeV\talpha_cb\tn_cb"
                    "\tsigma_over_E_cb"
                    "\tmu_cb_pre_MeV\tsigma_cb_pre_MeV\talpha_cb_pre\tn_cb_pre"
                    "\tsigma_over_E_cb_pre\n")
            for path, nmod, ne, mu_pre, sig_pre, mu_post, sig_post, \
                mu_cb, sig_cb, a_cb, n_cb, \
                mu_cb_pre, sig_cb_pre, a_cb_pre, n_cb_pre in [
                ('legacy', cal['n_modules_l'], cal['n_events_l'],
                 cal['mu_l_pre'], cal['sig_l_pre'],
                 cal['mu_l_post'], cal['sig_l_post'],
                 cal['mu_cb_l_post'], cal['sig_cb_l_post'], cal['alpha_l_post'], cal['n_l_post'],
                 cal['mu_cb_l_pre'], cal['sig_cb_l_pre'], cal['alpha_l_pre'], cal['n_l_pre']),
                ('new', cal['n_modules_n'], cal['n_events_n'],
                 cal['mu_n_pre'], cal['sig_n_pre'],
                 cal['mu_n_post'], cal['sig_n_post'],
                 cal['mu_cb_n_post'], cal['sig_cb_n_post'], cal['alpha_n_post'], cal['n_n_post'],
                 cal['mu_cb_n_pre'], cal['sig_cb_n_pre'], cal['alpha_n_pre'], cal['n_n_pre']),
            ]:
                ratio = (sig_post / mu_post if (mu_post and mu_post > 0 and not math.isnan(mu_post))
                         else float('nan'))
                ratio_cb = (sig_cb / mu_cb if (mu_cb and mu_cb > 0 and not math.isnan(mu_cb))
                            else float('nan'))
                ratio_cb_pre = (sig_cb_pre / mu_cb_pre if (mu_cb_pre and mu_cb_pre > 0 and not math.isnan(mu_cb_pre))
                                else float('nan'))
                f.write(f"{path}\t{nmod}\t{ne}\t{mu_pre:.2f}\t{sig_pre:.2f}"
                        f"\t{mu_post:.2f}\t{sig_post:.2f}\t{ratio:.4f}"
                        f"\t{mu_cb:.2f}\t{sig_cb:.2f}\t{a_cb:.3f}\t{n_cb:.3f}"
                        f"\t{ratio_cb:.4f}"
                        f"\t{mu_cb_pre:.2f}\t{sig_cb_pre:.2f}\t{a_cb_pre:.3f}\t{n_cb_pre:.3f}"
                        f"\t{ratio_cb_pre:.4f}\n")
            f.write("\n# inner-ring histogram bin centres + counts (pre/post recal)\n")
            f.write("e_centre\tlegacy_pre\tlegacy_post\tnew_pre\tnew_post\n")
            for c, lp, lq, np_, nq in zip(cal['centres'],
                                          cal['counts_l_pre'], cal['counts_l_post'],
                                          cal['counts_n_pre'], cal['counts_n_post']):
                f.write(f"{c:.1f}\t{lp}\t{lq}\t{np_}\t{nq}\n")
    print(f"\n  wrote: {tsv_path}")

    # ---- LMS per-batch summary + elastic-window module table ---------------
    base_path = tsv_path[:-4]
    lms_tsv = base_path + "_lms.tsv"
    write_lms_tsv(lms_tsv, lms, n_ev_per_batch, touched)
    print(f"  wrote: {lms_tsv}")
    if lms_path is not None:
        try:
            lms_png = base_path + "_lms.png"
            plot_lms_factors(lms_png, lms, n_ev_per_batch)
            print(f"  wrote: {lms_png}")
        except Exception as lms_plot_exc:
            sys.stderr.write(f"[WARN] LMS PNG plot failed: {lms_plot_exc}\n")
            try:
                plt.close('all')
            except Exception:
                pass
    if rows_by_path:
        win_tsv = base_path + "_window.tsv"
        write_window_tsv(win_tsv, rows_by_path, args, eprime_by_path,
                         sigma_glob_by_path, mu_global_by_path,
                         dict(prefit=f"({prefit[0]:.4f},{prefit[1]:.4f})",
                              target=target, label=kin_label, mass=m_target,
                              ebeam=ebeam))
        print(f"  wrote: {win_tsv}")

    # ---- Headline plot (safe-guard against matplotlib crashes) -------------
    png_path = args.out_path if args.out_path.endswith('.png') \
               else args.out_path + ".png"
    try:
        fig, ax = plt.subplots(figsize=(9, 6))
        ax.step(centres, counts_l, where='mid', color='#888',
                lw=1.5, label=f'legacy   N={n_legacy}')
        ax.step(centres, counts_n, where='mid', color='#1f77b4',
                lw=1.5, label=f'new W={W:g} ns   N={n_new}')

        if HAVE_SCIPY and not math.isnan(mu_l):
            x_fit = np.linspace(mu_l - 4 * sig_l, mu_l + 4 * sig_l, 200)
            ax.plot(x_fit, gauss(x_fit, A_l, mu_l, sig_l),
                    color='#888', ls='--', lw=1.0,
                    label=f'   μ={mu_l:.0f} σ={sig_l:.0f} MeV')
        if HAVE_SCIPY and not math.isnan(mu_n):
            x_fit = np.linspace(mu_n - 4 * sig_n, mu_n + 4 * sig_n, 200)
            ax.plot(x_fit, gauss(x_fit, A_n, mu_n, sig_n),
                    color='#1f77b4', ls='--', lw=1.0,
                    label=f'   μ={mu_n:.0f} σ={sig_n:.0f} MeV')

        ax.axvline(args.signal_min, color='#d62728', lw=0.8, ls=':',
                   label=f'cut E > {args.signal_min:.0f} MeV')
        ax.set_xlabel('cluster energy (MeV)')
        ax.set_ylabel('events / bin')
        ax.set_title(f"Single-cluster energy: legacy vs seed-time-gated\n"
                     f"run window [{pre_lo:.0f}, {pre_hi:.0f}] ns, "
                     f"{n_phys} physics events, {n_kept} triggered  [{lms_tag}]")
        ax.grid(alpha=0.3)
        ax.legend(loc='upper right', fontsize=9)
        ax.set_xlim(e_lo, e_hi)
        fig.tight_layout()
        fig.savefig(png_path, dpi=130)
        plt.close(fig)
        print(f"  wrote: {png_path}")
    except Exception as plot_exc:
        sys.stderr.write(
            f"[WARN] headline PNG plot failed: {plot_exc}\n"
            f"       TSV written to {tsv_path} is still valid.\n"
        )
        try:
            plt.close('all')  # cleanup any dangling figures
        except Exception:
            pass

    # ---- Optional: calibration before/after plot --------------------------
    if cal is not None:
        try:
            cal_png = png_path[:-4] + "_calibrated.png"
            fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), sharey=True)

            for ax, label, c_pre, c_post, mu_pre, sig_pre, mu_post, sig_post, \
                mu_cb, sig_cb, alpha_cb, n_cb, A_cb, A_g, base_col in [
                (axes[0], 'legacy',
                 cal['counts_l_pre'], cal['counts_l_post'],
                 cal['mu_l_pre'], cal['sig_l_pre'],
                 cal['mu_l_post'], cal['sig_l_post'],
                 cal['mu_cb_l_post'], cal['sig_cb_l_post'], cal['alpha_l_post'],
                 cal['n_l_post'], cal['A_cb_l_post'], cal['A_l_post'], '#888'),
                (axes[1], f'new (W={W:g} ns)',
                 cal['counts_n_pre'], cal['counts_n_post'],
                 cal['mu_n_pre'], cal['sig_n_pre'],
                 cal['mu_n_post'], cal['sig_n_post'],
                 cal['mu_cb_n_post'], cal['sig_cb_n_post'], cal['alpha_n_post'],
                 cal['n_n_post'], cal['A_cb_n_post'], cal['A_n_post'], '#1f77b4'),
            ]:
                ax.step(cal['centres'], c_pre, where='mid', color=base_col,
                        alpha=0.45, lw=1.3,
                        label=f'pre   μ={mu_pre:.0f}  σ_G={sig_pre:.0f} MeV')
                ax.step(cal['centres'], c_post, where='mid', color=base_col,
                        lw=1.8,
                        label=f'post  μ={mu_post:.0f}  σ_G={sig_post:.0f} MeV')
                # Crystal-ball curve overlay (guard against numeric overflow)
                if not math.isnan(mu_cb) and A_cb > 0:
                    try:
                        xfine = np.linspace(mu_cb - 6 * sig_cb,
                                            mu_cb + 4 * sig_cb, 400)
                        with np.errstate(over='ignore', invalid='ignore'):
                            yfine = crystal_ball(xfine, A_cb, mu_cb, sig_cb,
                                                 alpha_cb, n_cb)
                        # Filter out NaN/inf from the overflow
                        yfine = np.where(np.isfinite(yfine), yfine, 0.0)
                        ax.plot(xfine, yfine, color='#d62728', lw=1.4,
                                label=f'CB    μ={mu_cb:.0f}  σ={sig_cb:.0f} MeV  '
                                      f'(α={alpha_cb:.2f} n={n_cb:.1f})')
                    except Exception as cb_exc:
                        sys.stderr.write(
                            f"[WARN] CB curve overlay failed for {label}: {cb_exc}\n"
                        )
                # Gaussian fit of the post spectrum: dashed, same colour as CB
                if not (math.isnan(mu_post) or math.isnan(sig_post)
                        or math.isnan(A_g)) and A_g > 0:
                    xg = np.linspace(mu_post - 3 * sig_post,
                                     mu_post + 3 * sig_post, 300)
                    ax.plot(xg, gauss(xg, A_g, mu_post, sig_post),
                            color='#d62728', lw=1.4, ls='--',
                            label=f'Gauss μ={mu_post:.0f}  σ={sig_post:.0f} MeV '
                                  f'(dashed)')
                ax.axvline(target_ref, color='#444', ls=':', lw=0.8)
                ax.set_xlabel('cluster energy (MeV)')
                ax.set_title(label)
                ax.grid(alpha=0.3)
                ax.legend(loc='upper right', fontsize=8.5)
            axes[0].set_ylabel('events / bin')
            excl_str = (f",  exclude rows [{args.exclude_rowcol[0]},"
                        f"{args.exclude_rowcol[1]}]×cols "
                        f"[{args.exclude_rowcol[2]},{args.exclude_rowcol[3]}]"
                        if args.exclude_rowcol else "")
            sigEcb_l = (cal['sig_cb_l_post']/cal['mu_cb_l_post']*100
                        if cal['mu_cb_l_post'] and not math.isnan(cal['mu_cb_l_post'])
                        else float('nan'))
            sigEcb_n = (cal['sig_cb_n_post']/cal['mu_cb_n_post']*100
                        if cal['mu_cb_n_post'] and not math.isnan(cal['mu_cb_n_post'])
                        else float('nan'))
            fig.suptitle(
                f"Inner-ring per-module gain refinement "
                f"(r < {args.inner_radius:.0f} mm,  "
                f"seed E ≥ {args.min_seed_energy:.0f} MeV,  "
                f"≥ {args.min_per_module} ev/mod{excl_str})  "
                f"[{lms_tag}, {args.cal_method}]\n"
                f"σ_E/E (Gauss):  legacy = {cal['sig_l_post']/cal['mu_l_post']*100:.2f} %   "
                f"new = {cal['sig_n_post']/cal['mu_n_post']*100:.2f} %     "
                f"(CB):  legacy = {sigEcb_l:.2f} %   new = {sigEcb_n:.2f} %",
                fontsize=10.5)
            fig.tight_layout()
            fig.savefig(cal_png, dpi=130)
            plt.close(fig)
            print(f"  wrote: {cal_png}")
        except Exception as cal_plot_exc:
            sys.stderr.write(
                f"[WARN] calibrated PNG plot failed: {cal_plot_exc}\n"
                f"       TSV written to {tsv_path} is still valid.\n"
            )
            try:
                plt.close('all')
            except Exception:
                pass

    return 0


if __name__ == "__main__":
    sys.exit(main())
