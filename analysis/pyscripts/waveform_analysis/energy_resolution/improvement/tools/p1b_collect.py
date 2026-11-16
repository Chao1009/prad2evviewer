#!/usr/bin/env python3
"""p1b_collect.py: AMP P1b budget table AMP_P1b_excess.tsv (variance in percent^2 of the cluster energy, robust-sigma basis) from the P1b tables.
Primary numbers: report arms on the common elastic clusters (AMP_P1b_arms*.tsv, includes the gain solve), cross-check: direct cluster-level D
(AMP_P1b_tables/direct/direct.tsv, gain-solve free), prediction: hit-level fit (AMP_P1b_tables/budget/predict.tsv)."""
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
import energy_resolution_benchmark as bm        # noqa: E402

IMP = HERE.parent
T = IMP / "AMP_P1b_tables"
con = {r["contrast"]: (float(r["d_var_pct2"]), float(r["d_var_se_pct2"])) for r in bm.tsv_rows(IMP / "AMP_P1b_arms_contrasts.tsv")}
dr = {r["variant"]: (float(r["var_pct2"]), float(r["se_pct2"])) for r in bm.tsv_rows(T / "direct" / "direct.tsv")}
pr = {r["component"]: float(r["cluster_sigma_pct"]) for r in bm.tsv_rows(T / "budget" / "predict.tsv")}
fit = {r["column"]: r for r in bm.tsv_rows(T / "s50" / "fit.tsv")}
MODEL = 0.15                                        # relative modelling uncertainty of a predicted component (fit weighting, binning, split to split)
tot, tot_se = con["integral -> raw (total)"]
rows = []


def add(cause, v, se, how, recover, src):
    rows.append([cause, v, se, how, recover, src])


ph, ph_se = con["phase table (vs integral)"]
am, am_se = con["amplitude table after phase table"]
ph, am = -ph, -am                                   # recovered variance (positive)
D = dr["rawpk: all members"][0]
add("total excess raw over integral (robust sigma, 19 k elastic clusters)", tot, tot_se, "arms", "-", "AMP_P1b_arms.tsv")
add("total, direct cluster-level D (variance, mean removed)", D, dr["rawpk: all members"][1], "direct", "-", "direct.tsv")
add("1 sampling-phase curve of the raw peak against the integral (per-module 16-bin table)", ph, ph_se, "arms: raw minus raw with table", "calibration (P5 M1-type table)", "AMP_P1b_arms_contrasts.tsv")
add("1 (direct)", D - dr["rawpk: phase table"][0], np.nan, "direct", "calibration", "direct.tsv")
add("2 amplitude-dependent bias g(A) of raw / integral (swing 0.9 %), after the phase table", am, am_se, "arms", "calibration (small)", "AMP_P1b_arms_contrasts.tsv")
add("2 (direct)", dr["rawpk: phase table"][0] - dr["rawpk: phase + amplitude tables"][0], np.nan, "direct", "calibration", "direct.tsv")
kn, kn_se = con["pile-up hits only (raw in kn)"]
add("3 pile-up members (h_k >= 2 or extra pulse), raw only in those hits, no tables", kn, kn_se, "arms", "estimator / fitter (not calibration)", "AMP_P1b_arms_contrasts.tsv")
add("3 (direct)", D - dr["rawpk: single-pulse members only (pile-up members u=0)"][0], np.nan, "direct", "estimator", "direct.tsv")
a, c, n = float(fit["rawpk"]["a_adc"]), float(fit["rawpk"]["c"]), float(fit["rawpk"]["n_adc"])
for lab, k in (("4 single-sample fluctuation ~ 1 / sqrt(A) (a / A term of the hit scatter)", "stochastic a/A"), ("5 amplitude-proportional shape term (c, hit scatter 0.6 %)", "constant c"),
               ("6 white pedestal noise of one sample (n / A term, fitted n)", "white noise n/A")):
    v = pr[k] ** 2
    add(lab, v, MODEL * v, "predicted: hit-level fit propagated to clusters", "intrinsic to the single sample", "budget/predict.tsv")
v = (pr["noise in excess of the matched filter (psd sqrt(1 - 1/sumsq)), all members"]) ** 2
add("6 (analytic) white noise from the pedestal rms, 0.84 ADC excess over the matched filter", v, 0.2 * v, "predicted: pedestal statistics", "intrinsic", "budget/predict.tsv")
add("7 cross terms of the member shifts after both tables (intra-cluster correlation; direct)", dr["rawpk: phase + amplitude tables: cross terms = variance of D minus the independent part"][0],
    dr["rawpk: phase + amplitude tables: cross terms = variance of D minus the independent part"][1], "direct", "none: consistent with zero (phase harmonic fit removes 0.0005)", "direct.tsv")
add("arms total minus direct total: gain solve and position-map cost of the noisier estimator (not isolated)", tot - D, tot_se, "difference", "statistics of the calibration (more runs)", "")
add("remainder after phase and amplitude tables (raw, direct; all independent hit scatter)", dr["rawpk: phase + amplitude tables"][0], dr["rawpk: phase + amplitude tables"][1], "direct", "-", "direct.tsv")
known = ph + max(am, 0) + kn + pr["stochastic a/A"] ** 2 + pr["constant c"] ** 2 + pr["white noise n/A"] ** 2
add("sum of causes 1 to 6 (arms for 1-3, 2 taken as >= 0)", known, np.sqrt(ph_se ** 2 + kn_se ** 2 + (MODEL * (pr["stochastic a/A"] ** 2 + pr["constant c"] ** 2)) ** 2), "sum", "-", "")
add("unexplained: total minus sum of causes 1 to 6", tot - known, np.sqrt(tot_se ** 2 + ph_se ** 2 + kn_se ** 2 + (MODEL * (pr["stochastic a/A"] ** 2 + pr["constant c"] ** 2)) ** 2), "difference", "overlaps the arms-minus-direct row (calibration solve), not intra-cluster correlation", "")
bm.write_tsv(IMP / "AMP_P1b_excess.tsv", ("cause", "variance_pct2", "se_pct2", "method", "recoverable_by", "source"), rows)
print(open(IMP / "AMP_P1b_excess.tsv").read())
