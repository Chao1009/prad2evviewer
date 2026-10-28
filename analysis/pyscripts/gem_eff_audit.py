#!/usr/bin/env python3
"""
gem_eff_audit.py — offline LOO audit of the GEM tracking-efficiency monitor.

Purpose
-------
For each test detector D in {0,1,2,3}, the OTHER 3 GEMs + HyCal define a
straight-line anchor track; we then project the anchor to D and count
whether D recorded a hit at the predicted position.  D itself contributes
nothing to the anchor (matching candidate set excludes D, fit excludes D),
so the test is genuinely unbiased.  Three LOO variants are run side-by-side
on every event so they can be compared:

  * loo              GEM-seeded anchor (HyCal + each OTHER GEM hit drawn as a
                     seed line; lowest-χ² survivor wins).  Fit through
                     HyCal + 3 OTHER GEMs.  Vertex-agnostic — accepts any
                     straight track that lights up the 3 OTHER GEMs.
  * loo-target-in    Same GEM-seeded matching, but the fit additionally
                     includes (target_x, target_y, target_z) as a weighted
                     measurement with σ = `--sigma-target-x/-y/-z`
                     (≈ beam spot).
                     Pulls the line toward the target → kills upstream halo
                     and beam-gas tracks; the χ² gate then rejects anything
                     that doesn't actually point back to (0,0,0).
  * loo-target-seed  Single-seed variant: the seed line is
                     (target → HyCal cluster); no GEM-pair seeding.
                     Fit through HyCal + 3 OTHER GEMs (target only seeds,
                     never enters the fit).  Cheapest of the three.

Match definition (applied to the anchor in every variant):
  - all 3 OTHER detectors found a hit within match_nsigma · σ_total of the
    seed-line projection (σ_total = sqrt(σ_HC@gem² + σ_GEM²)),
  - weighted line fit's χ²/dof ≤ max_chi2,
  - per-detector residual against the FIT line within match_nsigma · σ_GEM[d].

Test (applied at D after the anchor passes):
  - D has a hit within match_nsigma · σ_GEM[D] of the projection.
  - Increment denominator[D] (anchor was good); increment numerator[D]
    if the test hit was found.

Event filter: ≥1 HyCal cluster AND ≥3 GEM detectors with at least one
cluster — minimum to make any anchor possible.

Usage
-----
    python analysis/pyscripts/gem_eff_audit.py <evio_path> <out_dir> \\
        [--match-nsigma 3.0] [--max-chi2 3.5] [--sigma-gem MM] \\
        [--sigma-target-x MM] [--sigma-target-y MM] [--sigma-target-z MM] \\
        [--min-cluster-energy 500] [--max-hits-per-det 50] \\
        [--n-dets 4] [--n-event-plots 2] [--max-events N]
        (plus the common run/database overrides, see --help)

`<evio_path>` accepts a glob (quote it), directory, or single split.
`<out_dir>` is created if missing; the script writes `efficiency.png`,
`anchor_chi2.png` (two-row anchor-quality plot — χ²/dof distribution +
cumulative acceptance, NOT detector efficiency) and, per LOO mode,
`residuals_<variant>.png`, `eff_ineff_<variant>.png` and up to
--n-event-plots eff / ineff `event_<variant>_*.png` displays.  Detector
efficiency numbers go to stdout in the text summary.
"""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

import _common as C
from _common import det


def find_closest(hits: Sequence[Tuple[float, float, float]],
                 pred_x: float, pred_y: float, cut: float
                 ) -> Tuple[int, float]:
    """Return (index, dr) of the lab-frame hit closest to (pred_x, pred_y)
    within `cut`, or (-1, +inf) if none qualify."""
    best_idx = -1
    best_d2 = cut * cut
    for i, h in enumerate(hits):
        dx = h[0] - pred_x
        dy = h[1] - pred_y
        d2 = dx * dx + dy * dy
        if d2 < best_d2:
            best_d2 = d2
            best_idx = i
    if best_idx < 0:
        return -1, float("inf")
    return best_idx, math.sqrt(best_d2)


# The seed-line projection window is `σ_total = sqrt(σ_HC@gem² + σ_GEM²)`
# because that's the uncertainty before the fit refines the line.  The
# post-fit gate below uses just σ_GEM[d] because by then HyCal is already
# in the fit and the projection at any plane is dominated by the GEM
# resolution.
def fit_residuals_within_window(
        fit_params: Tuple[float, float, float, float],
        matched: Sequence[bool],
        cand: Sequence[Optional[Tuple[float, float, float]]],
        params: "TrackingParams",
        ) -> bool:
    ax, bx, ay, by = fit_params
    for d, m in enumerate(matched):
        if not m or cand[d] is None:
            continue
        h = cand[d]
        pred_x = ax + bx * h[2]
        pred_y = ay + by * h[2]
        cut = params.match_nsigma * params.gem_pos_res[d]
        dx = h[0] - pred_x
        dy = h[1] - pred_y
        if dx * dx + dy * dy > cut * cut:
            return False
    return True


# ---- Tracking-efficiency runs ----

@dataclass
class TrackingParams:
    match_nsigma:      float
    max_chi2:          float
    max_hits_per_det:  int
    gem_pos_res:       List[float]   # per-detector (mm), len 4

@dataclass
class TrackResult:
    chi2_per_dof: float
    matched:      List[bool]                      # len n_dets
    cand:         List[Optional[Tuple[float, float, float]]]
    fit:          Tuple[float, float, float, float]   # ax, bx, ay, by


def _try_seed(seed_d: int, seed_idx: int,
              hcx: float, hcy: float, hcz: float,
              sigma_hc: float,
              hits_by_det: List[List[Tuple[float, float, float]]],
              gem_z: Sequence[float],
              candidate_dets: Sequence[int],
              params: TrackingParams,
              min_match: int,
              target: Tuple[float, float, float] = (0.0, 0.0, 0.0),
              target_in_fit: bool = False,
              sigma_target: Tuple[float, float, float] = (1.0, 1.0, 20.0),
              cap_candidates: bool = False,
              diag: Optional[Dict[str, List[int]]] = None,
              test_d: int = -1,
              ) -> Optional[TrackResult]:
    """Build one track, like the C++ trySeed(S, ...).  With seed_d >= 0 the
    seed line runs from HyCal through hit `seed_idx` of GEM `seed_d`, which
    counts as matched (seed_d must be in `candidate_dets`).  With seed_d < 0
    it runs from `target` to HyCal and nothing is pre-matched; such tracks
    point back to the target, which rejects upstream halo / vacuum
    interactions.  The other `candidate_dets` are matched within
    match_nsigma · σ_total of the seed line; at least `min_match` matched
    detectors are needed for a fit.

    If `target_in_fit`, append `target` to the weighted fit as a soft
    "track originated at target" constraint.  σ_target_z (the target's
    longitudinal extent) couples to the transverse measurement at
    z = target_z via the track slope:
        σ_x_eff² = σ_target_x² + (bx_est · σ_target_z)²
        σ_y_eff² = σ_target_y² + (by_est · σ_target_z)²
    where the slope estimate comes from the (target → HyCal cluster)
    line; for central tracks this leaves the σ_x,y contributions
    unchanged, while peripheral tracks (large |HyCal x,y|) get the
    expected σ_z lever-arm widening.

    `diag` (with `test_d`) collects the per-stage counters n_3matched,
    n_pass_chi2 and n_pass_resid."""
    tx, ty, tz = target
    n_dets = len(hits_by_det)
    matched = [False] * n_dets
    cand: List[Optional[Tuple[float, float, float]]] = [None] * n_dets
    if seed_d < 0:
        ax_s, bx_s, ay_s, by_s = det.seed_line(tx, ty, tz, hcx, hcy, hcz)
    else:
        if seed_d not in candidate_dets:
            return None
        hits_s = hits_by_det[seed_d]
        if seed_idx >= len(hits_s):
            return None
        g0 = hits_s[seed_idx]
        ax_s, bx_s, ay_s, by_s = det.seed_line(hcx, hcy, hcz,
                                               g0[0], g0[1], g0[2])
        matched[seed_d] = True
        cand[seed_d] = g0

    # σ_HC is 0 at the target and sigma_hc at HyCal; at z_gem the projected
    # positional uncertainty scales with the lever arm from the target.
    lever_hc = (hcz - tz) if hcz != tz else 1.0
    cap = params.max_hits_per_det if cap_candidates else 0
    for d in candidate_dets:
        if matched[d]:
            continue
        zd = gem_z[d]
        pred_x = ax_s + bx_s * zd
        pred_y = ay_s + by_s * zd
        s_hc_at_gem = sigma_hc * abs((zd - tz) / lever_hc)
        s_gem = params.gem_pos_res[d]
        s_total = math.sqrt(s_hc_at_gem * s_hc_at_gem + s_gem * s_gem)
        cut = params.match_nsigma * s_total
        # The C++ findClosest caps candidates at gem_eff_max_hits_per_det in
        # every mode (also the test-detector probe in _record_loo); here
        # only the target-seeded mode does.
        cand_hits = hits_by_det[d][:cap] if cap > 0 else hits_by_det[d]
        idx, _ = find_closest(cand_hits, pred_x, pred_y, cut)
        if idx >= 0:
            matched[d] = True
            cand[d] = cand_hits[idx]

    if sum(matched) < min_match:
        return None
    if diag is not None:
        diag["n_3matched"][test_d] += 1

    # Weighted fit of HyCal + every matched GEM (+ target if requested).
    # x and y use the same per-point σ except at the target: σ_target_z
    # couples differently into σ_x_eff and σ_y_eff via the slope.
    z_arr: List[float]  = [hcz]
    x_arr: List[float]  = [hcx]
    y_arr: List[float]  = [hcy]
    w_x_arr: List[float] = [1.0 / (sigma_hc * sigma_hc)]
    w_y_arr: List[float] = [1.0 / (sigma_hc * sigma_hc)]
    if target_in_fit:
        # Slope estimate from (target → HyCal cluster) — used to widen the
        # target's transverse σ by the σ_z lever arm.  Falls back to zero
        # slope when HyCal sits on the target plane.
        if hcz != tz:
            bx_est = (hcx - tx) / (hcz - tz)
            by_est = (hcy - ty) / (hcz - tz)
        else:
            bx_est = by_est = 0.0
        stx, sty, stz = sigma_target
        sx_eff_sq = stx * stx + (bx_est * stz) ** 2
        sy_eff_sq = sty * sty + (by_est * stz) ** 2
        z_arr.append(tz)
        x_arr.append(tx)
        y_arr.append(ty)
        w_x_arr.append(1.0 / sx_eff_sq)
        w_y_arr.append(1.0 / sy_eff_sq)
    for d in range(n_dets):
        if not matched[d] or cand[d] is None:
            continue
        h = cand[d]
        z_arr.append(h[2])
        x_arr.append(h[0])
        y_arr.append(h[1])
        s_gem = params.gem_pos_res[d]
        w_x_arr.append(1.0 / (s_gem * s_gem))
        w_y_arr.append(1.0 / (s_gem * s_gem))

    fit = det.fit_weighted_line(z_arr, x_arr, y_arr, w_x_arr, w_y_arr)
    if fit is None:
        return None
    ax, bx, ay, by, chi2_per_dof = fit
    if chi2_per_dof > params.max_chi2:
        return None
    if diag is not None:
        diag["n_pass_chi2"][test_d] += 1
    if not fit_residuals_within_window((ax, bx, ay, by), matched, cand, params):
        return None
    if diag is not None:
        diag["n_pass_resid"][test_d] += 1
    return TrackResult(chi2_per_dof, matched, cand, (ax, bx, ay, by))


def best_track(hcx: float, hcy: float, hcz: float, sigma_hc: float,
               hits_by_det: List[List[Tuple[float, float, float]]],
               gem_z: Sequence[float], seed_dets: Sequence[int],
               candidate_dets: Sequence[int], params: TrackingParams,
               min_match: int,
               target: Tuple[float, float, float] = (0.0, 0.0, 0.0),
               target_in_fit: bool = False,
               sigma_target: Tuple[float, float, float] = (1.0, 1.0, 20.0),
               ) -> Optional[TrackResult]:
    """Try every (seed_d, seed_idx) pair.  Return the lowest-χ²/dof
    TrackResult that passes the chi²_max gate, or None.  The target args
    are forwarded to _try_seed unchanged."""
    best: Optional[TrackResult] = None
    for seed_d in seed_dets:
        if seed_d >= len(hits_by_det):
            continue
        n_seeds = min(len(hits_by_det[seed_d]), params.max_hits_per_det)
        for seed_idx in range(n_seeds):
            r = _try_seed(seed_d, seed_idx, hcx, hcy, hcz, sigma_hc,
                          hits_by_det, gem_z, candidate_dets, params,
                          min_match, target, target_in_fit, sigma_target)
            if r is None:
                continue
            if best is None or r.chi2_per_dof < best.chi2_per_dof:
                best = r
    return best


# Per-event accumulator — one set per LOO variant.  Each variant runs N times
# per HyCal cluster (one per test detector), so the denominator is per-detector
# (the "anchor exists" counter for that detector excluded), and the numerator
# is per-detector too.  chi2_list aggregates anchor χ²/dof across all test
# detectors.

@dataclass
class CapturedEvent:
    """One full LOO test snapshot, kept for the eff/ineff single-event plots.
    Lab-frame coordinates; positions are mm.  `matched_hits[d]` is the hit
    selected from GEM d when it was an anchor (None for the test detector
    and any non-anchor GEM); `test_hit` is the (matched) hit at the test
    detector if the LOO call counted as efficient, else None."""
    test_d:         int
    is_matched:     bool
    chi2_per_dof:   float
    hcx:            float
    hcy:            float
    hcz:            float
    target_x:       float
    target_y:       float
    target_z:       float
    all_hits:       List[List[Tuple[float, float, float]]]
    matched_hits:   List[Optional[Tuple[float, float, float]]]
    test_hit:       Optional[Tuple[float, float, float]]
    pred_x:         float
    pred_y:         float
    pred_z:         float
    fit:            Tuple[float, float, float, float]
    match_nsigma:   float
    sigma_gem_test: float


@dataclass
class LooStats:
    name:           str
    n_attempted:    List[int]         = field(default_factory=lambda: [0]*4)
    n_matched:      List[int]         = field(default_factory=lambda: [0]*4)
    chi2_list:      List[float]       = field(default_factory=list)
    residuals_x:    List[List[float]] = field(default_factory=lambda: [[],[],[],[]])
    residuals_y:    List[List[float]] = field(default_factory=lambda: [[],[],[],[]])
    # Predicted hit position at the test detector, in detector-local (mm),
    # split into "matched" (efficiency map) and "unmatched" (inefficiency
    # map) per test detector.
    eff_local_x:    List[List[float]] = field(default_factory=lambda: [[],[],[],[]])
    eff_local_y:    List[List[float]] = field(default_factory=lambda: [[],[],[],[]])
    ineff_local_x:  List[List[float]] = field(default_factory=lambda: [[],[],[],[]])
    ineff_local_y:  List[List[float]] = field(default_factory=lambda: [[],[],[],[]])
    # Up to N captured single-event snapshots per bucket (eff/ineff), for
    # event-display plots.  Capacity gated by `--n-event-plots`.
    captured_eff:   List["CapturedEvent"] = field(default_factory=list)
    captured_ineff: List["CapturedEvent"] = field(default_factory=list)


def _record_loo(stats: "LooStats", test_d: int,
                track: Optional[TrackResult],
                hits_by_det: List[List[Tuple[float, float, float]]],
                gem_z: Sequence[float],
                params: TrackingParams,
                test_xform,
                hc_lab: Optional[Tuple[float, float, float]] = None,
                target_lab: Optional[Tuple[float, float, float]] = None,
                n_event_plots: int = 0) -> None:
    """Bookkeeping for a single LOO test attempt.  If `track` is None the
    anchor failed (3-GEM fit didn't pass χ²/residual gates); we don't
    increment any counter — denominator only counts events where the
    OTHER 3 GEMs delivered a clean anchor.  `test_xform` is the test
    detector's DetectorTransform, used to convert the predicted hit
    position (lab) to detector-local for the efficiency / inefficiency
    histograms.

    When `hc_lab` and `target_lab` are supplied and `n_event_plots > 0`,
    the first `n_event_plots` eff and `n_event_plots` ineff outcomes are
    snapshotted into `stats.captured_eff` / `stats.captured_ineff` for
    the single-event display plots."""
    if track is None:
        return
    stats.chi2_list.append(track.chi2_per_dof)
    ax, bx, ay, by = track.fit
    zd = gem_z[test_d]
    pred_x = ax + bx * zd
    pred_y = ay + by * zd
    s_gem = params.gem_pos_res[test_d]
    cut = params.match_nsigma * s_gem
    idx, _ = find_closest(hits_by_det[test_d], pred_x, pred_y, cut)
    stats.n_attempted[test_d] += 1
    # Predicted hit position in test-detector local coords (for the heatmap).
    pred_lx, pred_ly, _ = test_xform.lab_to_local(pred_x, pred_y, zd)
    if idx >= 0:
        stats.n_matched[test_d] += 1
        h = hits_by_det[test_d][idx]
        stats.residuals_x[test_d].append(h[0] - pred_x)
        stats.residuals_y[test_d].append(h[1] - pred_y)
        stats.eff_local_x[test_d].append(pred_lx)
        stats.eff_local_y[test_d].append(pred_ly)
    else:
        stats.ineff_local_x[test_d].append(pred_lx)
        stats.ineff_local_y[test_d].append(pred_ly)

    if hc_lab is not None and target_lab is not None and n_event_plots > 0:
        is_matched = (idx >= 0)
        bucket = stats.captured_eff if is_matched else stats.captured_ineff
        if len(bucket) < n_event_plots:
            test_hit = hits_by_det[test_d][idx] if is_matched else None
            bucket.append(CapturedEvent(
                test_d         = test_d,
                is_matched     = is_matched,
                chi2_per_dof   = track.chi2_per_dof,
                hcx            = hc_lab[0],
                hcy            = hc_lab[1],
                hcz            = hc_lab[2],
                target_x       = target_lab[0],
                target_y       = target_lab[1],
                target_z       = target_lab[2],
                all_hits       = [list(hd) for hd in hits_by_det],
                matched_hits   = list(track.cand),
                test_hit       = test_hit,
                pred_x         = pred_x,
                pred_y         = pred_y,
                pred_z         = zd,
                fit            = track.fit,
                match_nsigma   = params.match_nsigma,
                sigma_gem_test = s_gem,
            ))


# ---- Main loop ----

def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    C.add_common_args(ap)
    ap.add_argument("--match-nsigma",     type=float, default=3.0,
                    help="Matching window in σ_total (mirrors monitor "
                         "default).")
    ap.add_argument("--max-chi2",         type=float, default=3.5,
                    help="χ²/dof gate for accepting an anchor fit.  Default "
                         "3.5 — picks up real tracks while cutting the "
                         "anchor χ² tail of pile-up combinations.  The "
                         "anchor_chi2.png plot sweeps a range around this.")
    ap.add_argument("--max-hits-per-det", type=int,   default=50,
                    help="Cap seed/match candidates per detector (mirrors "
                         "AppState::gem_eff_max_hits_per_det = 50).")
    ap.add_argument("--min-cluster-energy", type=float, default=500.0,
                    help="Minimum HyCal cluster energy (MeV) to enter the "
                         "denominator.  Mirrors AppState::"
                         "gem_eff_min_cluster_energy.  Default 500 MeV.")
    ap.add_argument("--sigma-gem",        type=float, default=None,
                    help="Override σ_GEM (mm) for ALL detectors — useful for "
                         "absorbing residual alignment into a 'working' σ "
                         "before per-detector calibration is applied.  "
                         "Default: use reconstruction_config.json:matching:"
                         "gem_pos_res.")
    ap.add_argument("--sigma-target-x",   type=float, default=None,
                    help="σ_target_x (mm) — transverse beam-spot size in x. "
                         "Used by loo-target-in only.  Default: from "
                         "reconstruction_config.json:matching:target_pos_res.")
    ap.add_argument("--sigma-target-y",   type=float, default=None,
                    help="σ_target_y (mm) — transverse beam-spot size in y. "
                         "Default: from reconstruction_config.json.")
    ap.add_argument("--sigma-target-z",   type=float, default=None,
                    help="σ_target_z (mm) — target longitudinal extent.  "
                         "Couples to σ_x_eff and σ_y_eff via the track "
                         "slope at the target plane.  Default: from "
                         "reconstruction_config.json.")
    ap.add_argument("--n-dets",           type=int,   default=4,
                    help="Number of GEM detectors to consider (default 4).")
    ap.add_argument("--n-event-plots",    type=int,   default=2,
                    help="Save up to N eff and N ineff single-event display "
                         "plots per LOO mode.  Each plot shows GEM hits + "
                         "HyCal cluster + fit line on the X-Y (top-down) "
                         "and Z-Y (side) planes in lab coords.  0 disables. "
                         "Default 2 → up to 12 single-event PNGs total.")
    args = ap.parse_args(argv)

    out_dir = Path(args.out_path)
    out_dir.mkdir(parents=True, exist_ok=True)

    p = C.setup_pipeline_from_args(args)

    (pr_A, pr_B, pr_C), gem_pos_res, target_pos_res = C.load_matching_config(p)
    if args.sigma_gem is not None:
        gem_pos_res = [args.sigma_gem] * 4
    cfg_tgt_x, cfg_tgt_y, cfg_tgt_z = target_pos_res
    sigma_target_x = (args.sigma_target_x
                      if args.sigma_target_x is not None else cfg_tgt_x)
    sigma_target_y = (args.sigma_target_y
                      if args.sigma_target_y is not None else cfg_tgt_y)
    sigma_target_z = (args.sigma_target_z
                      if args.sigma_target_z is not None else cfg_tgt_z)
    sigma_target = (sigma_target_x, sigma_target_y, sigma_target_z)
    print(f"[setup] HC σ(E) = sqrt(({pr_A:.3f}/√E_GeV)²"
          f"+({pr_B:.3f}/E_GeV)²+{pr_C:.3f}²) mm")
    print(f"[setup] σ_GEM   = {gem_pos_res} mm"
          + ("  (overridden)" if args.sigma_gem is not None else ""))
    print(f"[setup] match_nsigma={args.match_nsigma}  "
          f"max_chi2={args.max_chi2}  max_hits_per_det={args.max_hits_per_det}")

    n_dets = min(args.n_dets, 4, p.gem_sys.get_n_detectors())
    gem_z = list(p.geo.gem_z[:n_dets])
    # Per-detector half-extents (mm) for the local-coord heatmaps.
    det_half: List[Tuple[float, float]] = []
    _det_cfgs = p.gem_sys.get_detectors()
    for d in range(n_dets):
        dc = _det_cfgs[d]
        det_half.append((dc.plane_x.size * 0.5, dc.plane_y.size * 0.5))
    while len(det_half) < 4:
        det_half.append((300.0, 300.0))
    print(f"[setup] GEM lab z = {gem_z}")
    print(f"[setup] target    = ({p.geo.target_x:g}, {p.geo.target_y:g}, "
          f"{p.geo.target_z:g}) mm  (from runinfo)")

    params = TrackingParams(
        match_nsigma     = args.match_nsigma,
        max_chi2         = args.max_chi2,
        max_hits_per_det = args.max_hits_per_det,
        gem_pos_res      = gem_pos_res,
    )
    print(f"[setup] min HyCal cluster E = {args.min_cluster_energy:.0f} MeV")

    overridden_target = any(a is not None for a in (
        args.sigma_target_x, args.sigma_target_y, args.sigma_target_z))
    print(f"[setup] σ_target = (x={sigma_target_x:.2f}, "
          f"y={sigma_target_y:.2f}, z={sigma_target_z:.2f}) mm  "
          f"(only used by loo-target-in"
          f"{', overridden' if overridden_target else ''})")

    loo            = LooStats("loo (GEM-seeded, HyCal + 3 OTHER GEMs in fit)")
    loo_target_in  = LooStats(
        f"loo-target-in (GEM-seeded, target+HyCal+3 OTHER GEMs in fit, "
        f"σ_target=(x={sigma_target_x:.2f},y={sigma_target_y:.2f},"
        f"z={sigma_target_z:.2f}) mm)")
    loo_target_seed = LooStats(
        f"loo-target-seed (target=({p.geo.target_x:g},{p.geo.target_y:g},"
        f"{p.geo.target_z:g})→HyCal seed, HyCal + 3 OTHER GEMs in fit)")
    loo_modes = (loo, loo_target_in, loo_target_seed)
    # Per-test_d stage counters for the loo-target-seed mode (matches the
    # server's default loo_mode).  Lets us pinpoint where Python and the
    # C++ server's anchor counts diverge.
    diag_target_seed: Dict[str, List[int]] = {
        "n_call":        [0]*4,   # target-seeded anchor attempts for test_d
        "n_3matched":    [0]*4,   # all 3 candidate dets matched in seed window
        "n_pass_chi2":   [0]*4,   # passed χ²/dof gate (after 3-match)
        "n_pass_resid":  [0]*4,   # passed per-det residual gate (= denominator)
    }
    n_events_filter_pass = n_used = 0
    target_lab_pt = (p.geo.target_x, p.geo.target_y, p.geo.target_z)

    stats = C.LoopStats()
    for fadc_evt, ssp_evt in C.iter_physics_events(
            p, stats, with_ssp=True, max_events=args.max_events,
            accept=C.passes_physics_trigger):
        hc_raw = C.reconstruct_hycal(p, fadc_evt)
        if not hc_raw:
            continue
        hc_lab = [(*C.hycal_to_lab(p, h), float(h.energy)) for h in hc_raw]

        C.reconstruct_gem(p, ssp_evt)
        hits_by_det: List[List[Tuple[float, float, float]]] = [
            [] for _ in range(n_dets)
        ]
        for d in range(n_dets):
            xform = p.gem_xforms[d]
            for g in p.gem_sys.get_hits(d):
                x, y, z = xform.to_lab(g.x, g.y)
                hits_by_det[d].append((x, y, z))

        # Event filter — at least 3 GEM detectors must each have
        # ≥1 cluster: a track needs 4 points (HyCal + ≥3 GEMs) to be a
        # non-trivial 4-parameter fit (dof = 2N − 4 ≥ 2).
        n_dets_with_hits = sum(1 for h in hits_by_det if h)
        if n_dets_with_hits < 3:
            continue
        n_events_filter_pass += 1

        for hcx, hcy, hcz, energy in hc_lab:
            if hcz <= 0:
                continue
            if energy < args.min_cluster_energy:
                continue
            sigma_hc = C.hycal_pos_resolution(pr_A, pr_B, pr_C, energy)
            n_used += 1

            hc_lab_pt = (hcx, hcy, hcz)

            for test_d in range(n_dets):
                others = [d for d in range(n_dets) if d != test_d]
                seeds_for_loo = [d for d in others if hits_by_det[d]]
                test_xform = p.gem_xforms[test_d]

                # --- loo / loo-target-in: GEM-seeded ---
                if seeds_for_loo:
                    for mode, target_in_fit in ((loo, False),
                                                (loo_target_in, True)):
                        track = best_track(hcx, hcy, hcz, sigma_hc,
                            hits_by_det, gem_z,
                            seed_dets=seeds_for_loo,
                            candidate_dets=others,
                            params=params, min_match=3,
                            target=target_lab_pt,
                            target_in_fit=target_in_fit,
                            sigma_target=sigma_target)
                        _record_loo(mode, test_d, track,
                                    hits_by_det, gem_z, params, test_xform,
                                    hc_lab=hc_lab_pt,
                                    target_lab=target_lab_pt,
                                    n_event_plots=args.n_event_plots)

                # --- loo-target-seed: target→HyCal seed, no target in fit ---
                diag_target_seed["n_call"][test_d] += 1
                track = _try_seed(-1, -1, hcx, hcy, hcz, sigma_hc,
                                  hits_by_det, gem_z, others, params,
                                  min_match=3, target=target_lab_pt,
                                  cap_candidates=True,
                                  diag=diag_target_seed, test_d=test_d)
                _record_loo(loo_target_seed, test_d, track,
                            hits_by_det, gem_z, params, test_xform,
                            hc_lab=hc_lab_pt,
                            target_lab=target_lab_pt,
                            n_event_plots=args.n_event_plots)

    print(f"[done] {stats.n_phys} physics events, "
          f"{n_events_filter_pass} pass filter (≥1 HyCal + ≥3 GEMs), "
          f"{n_used} HyCal clusters used, {stats.elapsed:.1f}s")

    print()
    print("loo-target-seed per-stage breakdown (per test_d):")
    print(f"  {'test_d':>8} {'n_call':>10} {'n_3matched':>12} "
          f"{'n_pass_chi2':>12} {'n_pass_resid':>13}")
    for d in range(4):
        print(f"  {d:>8} "
              f"{diag_target_seed['n_call'][d]:>10} "
              f"{diag_target_seed['n_3matched'][d]:>12} "
              f"{diag_target_seed['n_pass_chi2'][d]:>12} "
              f"{diag_target_seed['n_pass_resid'][d]:>13}")
    print()

    write_outputs(out_dir, loo_modes, params,
                  stats.n_phys, n_events_filter_pass, n_used,
                  args.min_cluster_energy, det_half, gem_z)
    return 0


# ---- Reporting ----

def _eff(num: int, den: int) -> str:
    if den == 0:
        return "—"
    return f"{100.0 * num / den:5.1f}%  ({num}/{den})"


def write_outputs(out_dir: Path,
                  loo_modes: Sequence["LooStats"],
                  params: TrackingParams,
                  n_phys: int, n_events_filter_pass: int, n_used: int,
                  min_cluster_energy: float,
                  det_half: Sequence[Tuple[float, float]],
                  gem_z: Sequence[float]
                  ) -> None:

    # ---- text summary (stdout) --------------------------------------------
    print()
    print("GEM tracking-efficiency audit (leave-one-out)")
    print("=============================================")
    print()
    print(f"physics events processed : {n_phys}")
    print(f"events passing filter    : {n_events_filter_pass}  "
          f"(≥1 HyCal cluster AND ≥3 GEM detectors with hits)")
    print(f"HyCal clusters considered: {n_used}  "
          f"(after E ≥ {min_cluster_energy:.0f} MeV cut)")
    print(f"match_nsigma             : {params.match_nsigma}")
    print(f"max_chi2_per_dof         : {params.max_chi2}")
    print(f"max_hits_per_det         : {params.max_hits_per_det}")
    print(f"σ_GEM (mm)               : {params.gem_pos_res}")
    print()
    print("Each LOO variant runs once per (HyCal cluster, test detector D):")
    print("  - the OTHER 3 GEMs anchor a line fit (all 3 matched in seed")
    print("    window, χ²/dof ≤ max gate, per-det residual within "
          f"{params.match_nsigma}·σ_GEM),")
    print("  - the fit is projected to D and a hit is searched within "
          f"{params.match_nsigma}·σ_GEM[D].")
    print("Denominator = anchors that succeeded; numerator = anchors with a "
          "hit at D.")
    print()

    for mode in loo_modes:
        print(f"--- {mode.name} ---")
        for d in range(4):
            print(f"  GEM{d}: {_eff(mode.n_matched[d], mode.n_attempted[d])}")
        if mode.chi2_list:
            arr = sorted(mode.chi2_list)
            med  = arr[len(arr)//2]
            p90  = arr[int(0.9 * len(arr))]
            p99  = arr[int(0.99 * len(arr))]
            print(f"  anchor χ²/dof: median={med:.3f}  "
                  f"90%={p90:.3f}  99%={p99:.3f}  (n={len(arr)})")
        print()

    # ---- plots (matplotlib optional) --------------------------------------
    plt = C.import_pyplot()
    if plt is None:
        print("[plot] matplotlib not available; skipping PNGs")
        return

    _plot_efficiency_bars(plt, loo_modes, out_dir / "efficiency.png")
    _plot_anchor_chi2(plt, loo_modes, params,
                      out_dir / "anchor_chi2.png")
    for mode in loo_modes:
        slug = mode.name.split()[0].replace("/", "_")
        _plot_residuals(plt, mode, out_dir / f"residuals_{slug}.png")
        _plot_eff_ineff_local(plt, mode, det_half,
                              out_dir / f"eff_ineff_{slug}.png")
        # Single-event displays: a few eff and ineff snapshots per mode.
        for kind, bucket in (("eff",   mode.captured_eff),
                              ("ineff", mode.captured_ineff)):
            for k, ev in enumerate(bucket):
                out_path = (out_dir
                            / f"event_{slug}_{kind}_{k:02d}_test{ev.test_d}.png")
                _plot_single_event(plt, ev, gem_z, slug, kind, k, out_path)


def _plot_efficiency_bars(plt, modes: Sequence[LooStats], out: Path) -> None:
    """One panel per GEM (2x2 grid).  Each panel stacks one horizontal
    progress bar per LOO variant: full width = 100%, filled portion =
    detector efficiency, annotation = "XX.X%  (num/den)".  The bar IS the
    efficiency — denominator goes into the count text on the right."""
    n_dets = 4
    n_modes = len(modes)
    fig, axes = plt.subplots(2, 2, figsize=(11, 6.5))
    axes_flat = axes.flatten()
    bar_h = 0.65
    for d in range(n_dets):
        ax = axes_flat[d]
        # Plot bottom-up so loo (mode 0) ends up on top
        ys = list(range(n_modes - 1, -1, -1))
        for i, m in enumerate(modes):
            n_match = m.n_matched[d]
            n_att   = m.n_attempted[d]
            eff     = (100.0 * n_match / n_att) if n_att > 0 else 0.0
            y       = ys[i]
            ax.barh(y, 100.0, height=bar_h, color="lightgray", alpha=0.4,
                    edgecolor="gray", linewidth=0.5, zorder=1)
            ax.barh(y, eff, height=bar_h, color=f"C{i}", alpha=0.85,
                    zorder=2)
            ax.text(-2, y, m.name.split()[0], va="center", ha="right",
                    fontsize=9, color=f"C{i}", fontweight="bold")
            ann = (f"{eff:5.1f}%  ({n_match}/{n_att})"
                   if n_att > 0 else "no anchors")
            ax.text(102, y, ann, va="center", ha="left", fontsize=9)
        ax.set_xlim(0, 140)
        ax.set_ylim(-0.6, n_modes - 0.4)
        ax.set_yticks([])
        ax.set_xticks([0, 25, 50, 75, 100])
        ax.set_xlabel("efficiency (%)" if d >= 2 else "")
        ax.set_title(f"GEM{d}", fontsize=11, fontweight="bold")
        ax.grid(True, axis="x", alpha=0.3, zorder=0)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle("LOO per-detector efficiency", fontsize=12)
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"[plot] {out}")


def _plot_anchor_chi2(plt, modes: Sequence[LooStats],
                      params: TrackingParams, out: Path) -> None:
    """Two-row anchor-quality plot for the LOO line fits (HyCal + 3 OTHER
    GEMs).  Top row: per-variant χ²/dof histograms.  Bottom row: cumulative
    fraction of anchors retained as the χ²/dof cut sweeps.  Both rows show
    *anchor* statistics — they describe the quality of the 3-GEM-anchor
    sample, NOT detector efficiency (which is the GEM0..GEM3 numbers in the
    text summary)."""
    if not any(m.chi2_list for m in modes):
        return
    fig, (ax_hist, ax_acc) = plt.subplots(2, 1, figsize=(9, 8),
                                          sharex=True)

    # ---- top: χ²/dof histograms ------------------------------------------
    bins = np.linspace(0, max(params.max_chi2, 1.0) * 1.1, 60)
    for i, m in enumerate(modes):
        if not m.chi2_list:
            continue
        label = f"{m.name.split()[0]}  (n={len(m.chi2_list)})"
        ax_hist.hist(m.chi2_list, bins=bins, color=f"C{i}", alpha=0.55,
                     label=label)
    ax_hist.axvline(params.max_chi2, color="k", ls="--", lw=1,
                    label=f"current cut = {params.max_chi2}")
    ax_hist.set_ylabel("LOO anchors accepted")
    ax_hist.set_title("LOO anchor χ²/dof distribution "
                      "(HyCal + 3 OTHER GEMs)")
    ax_hist.grid(True, alpha=0.3)
    ax_hist.legend(loc="best", fontsize=9)

    # ---- bottom: cumulative anchor acceptance ----------------------------
    cuts = np.linspace(0.5, max(params.max_chi2, 1.0) * 1.5, 60)
    for i, m in enumerate(modes):
        if not m.chi2_list:
            continue
        arr = np.asarray(m.chi2_list)
        accepted_at_cut = np.array(
            [(arr <= c).sum() for c in cuts], dtype=float)
        total = float(len(arr))
        ax_acc.plot(cuts, 100.0 * accepted_at_cut / total, color=f"C{i}",
                    lw=1.6, label=m.name.split()[0])
    ax_acc.axvline(params.max_chi2, color="k", ls="--", lw=1,
                   label=f"current cut = {params.max_chi2}")
    ax_acc.set_xlabel("χ²/dof cut")
    ax_acc.set_ylabel("anchors retained (%)  [NOT detector efficiency]")
    ax_acc.set_title("Cumulative anchor acceptance vs χ²/dof gate")
    ax_acc.grid(True, alpha=0.3)
    ax_acc.legend(loc="best", fontsize=9)

    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"[plot] {out}")


def _plot_eff_ineff_local(plt, loo: LooStats,
                          det_half: Sequence[Tuple[float, float]],
                          out: Path) -> None:
    """2 rows × 4 columns of local-coord 2D histograms for one LOO variant.
    Top row: efficiency map — predicted hit position (in detector-local mm)
    for every LOO test where the test detector also recorded a hit.
    Bottom row: inefficiency map — predicted positions where it didn't.
    Columns are GEM 0..3.  Axes are shared within columns (x) and within
    rows (y) so the eight panels pack tight; colorbars sit at the right
    edge of each row.

    Color scale is log per row, with a single LogNorm shared across all
    four detectors so colors compare directly across GEMs."""
    has_eff   = any(loo.eff_local_x[d]   for d in range(4))
    has_ineff = any(loo.ineff_local_x[d] for d in range(4))
    if not (has_eff or has_ineff):
        return
    from matplotlib.colors import LogNorm
    fig, axes = plt.subplots(
        2, 4, figsize=(13, 11), sharex="col", sharey="row",
        gridspec_kw={"wspace": 0.02, "hspace": 0.10,
                     "left": 0.07, "right": 0.94,
                     "top": 0.95, "bottom": 0.05})
    row_meta = (
        ("efficiency",   "expected and detected",     "viridis",
         lambda d: loo.eff_local_x[d],   lambda d: loo.eff_local_y[d]),
        ("inefficiency", "expected but not detected", "magma",
         lambda d: loo.ineff_local_x[d], lambda d: loo.ineff_local_y[d]),
    )
    for row, (kind, blurb, cmap, get_x, get_y) in enumerate(row_meta):
        # Pre-pass: compute per-detector bins and the global max count for
        # this row so all four panels share one LogNorm and colors compare
        # directly across GEMs.
        panel_specs = []
        row_max = 0
        for d in range(4):
            xmax, ymax = (det_half[d] if d < len(det_half) else (300.0, 300.0))
            bins_x = np.linspace(-xmax, xmax, 60)
            bins_y = np.linspace(-ymax, ymax, 60)
            xs, ys = get_x(d), get_y(d)
            if xs:
                H, _, _ = np.histogram2d(xs, ys, bins=[bins_x, bins_y])
                if H.size:
                    row_max = max(row_max, int(H.max()))
            panel_specs.append((bins_x, bins_y, xmax, ymax, xs, ys))
        # vmin=1 masks empty bins (count 0) — they show through as the
        # axes background, not as the lowest cmap color.  vmax floored at 2
        # so vmin<vmax even when a row has at most one count anywhere.
        norm = LogNorm(vmin=1, vmax=max(row_max, 2))
        last_im = None
        for d in range(4):
            bins_x, bins_y, xmax, ymax, xs, ys = panel_specs[d]
            ax = axes[row, d]
            n = len(xs)
            if n > 0:
                _, _, _, im = ax.hist2d(xs, ys, bins=[bins_x, bins_y],
                                        cmap=cmap, norm=norm)
                last_im = im
            else:
                ax.text(0.5, 0.5, "no data", ha="center", va="center",
                        transform=ax.transAxes)
            ax.set_aspect("equal")
            ax.set_xlim(-xmax, xmax)
            ax.set_ylim(-ymax, ymax)
            ax.set_title(f"GEM{d}  (n={n})", fontsize=10)
            if row == 1:
                ax.set_xlabel("local x (mm)")
            if d == 0:
                ax.set_ylabel(f"{kind}\n({blurb})\nlocal y (mm)",
                              fontsize=9)
            ax.tick_params(labelsize=8)
            ax.grid(True, alpha=0.2)
        # one shared colorbar per row, riding on the rightmost panel
        if last_im is not None:
            cax = fig.add_axes([0.945, axes[row, -1].get_position().y0,
                                0.012, axes[row, -1].get_position().height])
            fig.colorbar(last_im, cax=cax, label="counts (log)")
    fig.suptitle(f"LOO predicted hit positions at the test detector  "
                 f"[{loo.name.split()[0]}]", fontsize=11, y=0.98)
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"[plot] {out}")


def _draw_event_primitives(ax, ev: "CapturedEvent", ix: int, iy: int,
                           line_h, line_v, labels: bool) -> None:
    """Draw one snapshot's GEM hits, anchor matches, HyCal cluster, target,
    fit line (line_h, line_v), prediction and test hit on `ax`, using lab
    coordinates ix / iy (0 = x, 1 = y, 2 = z) as the horizontal / vertical
    axes.  Legend labels are attached only when `labels`."""
    lab = (lambda s: s) if labels else (lambda s: None)
    hc = (ev.hcx, ev.hcy, ev.hcz)
    tg = (ev.target_x, ev.target_y, ev.target_z)
    pr = (ev.pred_x, ev.pred_y, ev.pred_z)
    for d, hits in enumerate(ev.all_hits):
        for k, h in enumerate(hits):
            ax.plot(h[ix], h[iy], 'o', color=f"C{d}", markersize=4,
                    alpha=0.4, label=lab(f"GEM{d} hits") if k == 0 else None)
    first_anchor = True
    for d in range(len(ev.all_hits)):
        m = ev.matched_hits[d]
        if d == ev.test_d or m is None:
            continue
        ax.plot(m[ix], m[iy], 's', color=f"C{d}", markersize=12, mfc='none',
                mew=2, label=lab("anchor match") if first_anchor else None)
        first_anchor = False
    ax.plot(hc[ix], hc[iy], '*', color='black', markersize=16,
            label=lab('HyCal cluster'))
    ax.plot(tg[ix], tg[iy], 'x', color='black', markersize=12, mew=2,
            label=lab('target'))
    ax.plot(line_h, line_v, '--', color='red', lw=1.2, alpha=0.7,
            label=lab('fit'))
    ax.plot(pr[ix], pr[iy], 'P', color='red', markersize=14, mfc='none',
            mew=2, label=lab(f'pred @ GEM{ev.test_d}'))
    if ev.test_hit is not None:
        ax.plot(ev.test_hit[ix], ev.test_hit[iy], 'D', color='red',
                markersize=10, label=lab(f'GEM{ev.test_d} hit'))


def _plot_single_event(plt, ev: "CapturedEvent",
                       gem_z: Sequence[float],
                       mode_short: str, kind: str, idx: int,
                       out: Path) -> None:
    """One LOO snapshot rendered in lab coords on two panels:

      * X-Y (top-down): all GEM hits per detector, anchor matches squared,
        HyCal star, target ×, fit-line projection, predicted position at
        the test detector (red plus), and the test-detector hit (red
        diamond) when the LOO call was efficient.  A dotted circle of
        radius `match_nsigma · σ_GEM[test_d]` around the prediction
        visualizes the matching window.
      * Z-Y (side):  same primitives, with z on the horizontal axis so the
        full lever arm from target through GEMs to HyCal is visible.

    GEMs are color-coded C0..C3; the test detector's data is highlighted
    in red regardless of detector index."""
    fig, (ax_xy, ax_zy) = plt.subplots(1, 2, figsize=(14, 6))

    a_x, b_x, a_y, b_y = ev.fit
    zs_ref = [ev.target_z, ev.hcz] + list(gem_z)
    z_min  = min(zs_ref) - 50.0
    z_max  = max(zs_ref) + 50.0
    z_arr  = np.linspace(z_min, z_max, 200)
    line_x = a_x + b_x * z_arr
    line_y = a_y + b_y * z_arr

    # ----- X-Y (top-down) -------------------------------------------------
    _draw_event_primitives(ax_xy, ev, 0, 1, line_x, line_y, True)
    cut_radius = ev.match_nsigma * ev.sigma_gem_test
    circle = plt.Circle((ev.pred_x, ev.pred_y), cut_radius,
                        facecolor='none', edgecolor='red',
                        linestyle=':', lw=1.0, alpha=0.6,
                        label=f'{ev.match_nsigma:g}σ match window')
    ax_xy.add_patch(circle)
    ax_xy.set_xlabel('lab x (mm)')
    ax_xy.set_ylabel('lab y (mm)')
    ax_xy.set_title('X-Y (top-down)')
    ax_xy.set_aspect('equal', adjustable='datalim')
    ax_xy.grid(True, alpha=0.3)
    ax_xy.legend(loc='best', fontsize=7, framealpha=0.85)

    # ----- Z-Y (side view) ------------------------------------------------
    _draw_event_primitives(ax_zy, ev, 2, 1, z_arr, line_y, False)
    # Dotted verticals at the GEM planes — orient the eye to where each
    # detector lives along z.
    for d, zd in enumerate(gem_z):
        ax_zy.axvline(zd, color=f"C{d}", ls=':', lw=0.6, alpha=0.4)
    ax_zy.axvline(ev.hcz, color='black', ls=':', lw=0.6, alpha=0.4)
    ax_zy.set_xlabel('lab z (mm)')
    ax_zy.set_ylabel('lab y (mm)')
    ax_zy.set_title('Z-Y (side view)')
    ax_zy.grid(True, alpha=0.3)

    flag = "matched" if ev.is_matched else "unmatched"
    title = (f"[{mode_short}] {kind} event #{idx} — "
             f"test_d=GEM{ev.test_d} ({flag})  "
             f"χ²/dof={ev.chi2_per_dof:.2f}")
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"[plot] {out}")


def _plot_residuals(plt, loo: LooStats, out: Path) -> None:
    if not any(loo.residuals_x):
        return
    fig, axes = plt.subplots(2, 4, figsize=(16, 6), sharex='row')
    for d in range(4):
        for row, (label, arr) in enumerate(
            (("Δx (mm)", loo.residuals_x[d]),
             ("Δy (mm)", loo.residuals_y[d]))):
            ax = axes[row, d]
            if not arr:
                ax.text(0.5, 0.5, "no data", ha="center", va="center",
                        transform=ax.transAxes)
                ax.set_xlabel(label)
                ax.set_title(f"GEM{d}")
                continue
            data = np.asarray(arr)
            lo, hi = np.percentile(data, [1, 99])
            pad = max(1.0, (hi - lo) * 0.1)
            bins = np.linspace(lo - pad, hi + pad, 60)
            ax.hist(data, bins=bins, color=f"C{d}", alpha=0.8)
            med = float(np.median(data))
            std = float(np.std(data))
            ax.axvline(med, color="k", ls="--", lw=0.8,
                       label=f"med={med:.2f}\nσ={std:.2f}")
            ax.set_xlabel(label)
            ax.set_title(f"GEM{d}  (n={len(data)})")
            ax.legend(loc="best", fontsize=8)
            ax.grid(True, alpha=0.3)
    fig.suptitle(f"LOO residuals: hit − projected at test detector  "
                 f"[{loo.name.split()[0]}]")
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"[plot] {out}")


if __name__ == "__main__":
    sys.exit(main())
