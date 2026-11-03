# PbWO4 energy resolution: R5 without gain refinement (R5nr) as the working reference

**Date:** 2026-10-29
**Branch:** `waveform_jz` (R5nr code state: commit `71e049f`, benchmark sha1 08b19869; S0 state in `output/r5_20261028/configs/git_state.txt`)
**Plan:** `analysis/pyscripts/plans/plan_r5_20261028.md` (S0–S6; this note is S7)
**Previous notes:** `analysis/pyscripts/notes/PbWO4_energy_resolution_20261028.md` (§8 seed scan, §9 merge, §10 Møller, §11 Geant4 profile), `analysis/pyscripts/notes/PbWO4_energy_resolution_summary.md` (§1.4, §4)

Paths are relative to `analysis/pyscripts/` in `/Users/L00338853/work/PRad/prad2evviewer` unless stated; `$O` = `output/r5_20261028`. R4 = `output/resolution_fit_20261027_best50/r4_splits` (e-p/C) and `output/resolution_fit_20261028_ta/n200/splits` (Ta). R5nr = `$O/scan/r5nr_splits` (118 splits = 93 e-p/C + 25 Ta). Refined R5 = `$O/scan/r5_splits` (same recipe with the gain refinement applied; reference arm only). Fit convention from 10-29 (user decision 1): the headline is the CB core σ for every run type, and the iterative ±2σ Gaussian is a secondary column. Where the two disagree (|mean(d_CB − d_G)| > 2 SE or > 0.05 pp, paired), the row is flagged CB≠G and no fit is chosen. Δ in pp, paired per split, mean ± SE. Ta runs are reported per run (25800, 26138, 26230) and never pooled.

## 1. Status

- **User decision (2026-10-29), verbatim:** "Move on with R5 no refinement". R5nr is the working reference from today. R4 becomes the previous reference.
- **Promotion check (plan S6):** criteria 1, 2 and 6 PASS, 3 DISCUSS, 4 and 5 FAIL (`$O/tables/promotion.md`). Under the plan's rule R5nr is "reported, not promoted"; the user proceeds with it anyway. The change of reference is therefore a user decision and is not a pass of the criteria (§7).
- **Headline (legacy, CB).** On CB, e-p and C are within 0.011 pp of R4: Δ −0.000 ± 0.005 (3.49 GeV H2), −0.006 ± 0.003 (2.24 GeV H2), +0.011 ± 0.012 (2.24 GeV C). At 0.73 GeV, CB +0.442 ± 0.019 and Gauss −0.264 ± 0.008; the CB-vs-Gaussian question there is open (§8). Ta is worse on all three runs, by +0.046 / +0.145 / +0.152 pp. The 2-term a/√E ⊕ b fit is unphysical (b² < 0, χ² 203 for 1 dof).

## 2. What R5 changes relative to R4

- **Clustering seed** `min_center_energy` 10 → 30 MeV, applied in-process from `$O/configs/reconstruction_config_r5.json` with no `database/` write (20261028 §8).
- **log_weight_thres** 3.6 → 5.2, main's database value (standing rule of 10-29: always use main's current database values). Runinfo is main's merged `database/runinfo/general.json`.
- **Møller rejection** `--moller-reject kin` (n = 3, m = 2) is applied both in the sample and in the gains. The per-crystal solve drops Møller-tagged clusters (`percrystal_calib.py --moller-reject kin`), and `--moller-reject-gains on` is set (20261028 §10).
- **Geant4 cluster profile** `--cluster-profile geant4` (`output/cluster_profile_20261028/report.md`).
- **Fit.** Two-peak CB with χ²/dof selection against the Gaussian (`--cb-two-peak --fit-select chi2`, the benchmark defaults since `be24414`). The headline is the CB core σ.
- **Re-derived inputs.** Per-crystal gains at the R5 settings, Møller-rejected (`$O/gains/<run>/`). Position maps re-derived per run type; R5nr uses the refine-off maps `$O/poscorr/maps_norefine/<slug>/poscorr_r5nr_<slug>_all.json`.
- **Gain refinement off** (`--gain-refine off`, contract change of 10-29, confirmed by the user at 13:40). The online seed-indexed refinement is still measured (`result_window.tsv:gain`) but applied as 1, so post = pre. R4 values are post-refinement. No second gain pass was run.
- **Unchanged from R4:** template banks (best-50 for H2 and C, n200 for Ta), the 16 ns gate and 100–200 ns pre-window, the exclusion rectangle, the ±2σ Gaussian and the 118 splits.

## 3. σ_E/E, R5nr vs R4

Legacy path, % (mean ± SE over splits); Δ = R5nr − R4 in pp (`$O/tables/r3_r4_r5nr.md`; per split in `$O/ablation_nr/compare/r5_vs_r4{,_ta}_splits.tsv`):

| Point | n | CB R4 | **CB R5nr** | Δ CB | Gauss R4 | Gauss R5nr | Δ Gauss | flag |
|---|---:|---|---|---|---|---|---|---|
| 3.49 GeV H2 | 27 | 1.494 ± 0.008 | **1.494 ± 0.008** | −0.000 ± 0.005 | 1.573 ± 0.009 | 1.541 ± 0.005 | −0.032 ± 0.007 | CB≠G |
| 2.24 GeV H2 | 42 | 1.807 ± 0.005 | **1.800 ± 0.004** | −0.006 ± 0.003 | 1.855 ± 0.004 | 1.845 ± 0.003 | −0.010 ± 0.002 | CB≠G |
| 2.24 GeV C | 8 | 1.789 ± 0.016 | **1.799 ± 0.009** | +0.011 ± 0.012 | 1.861 ± 0.007 | 1.843 ± 0.007 | −0.018 ± 0.005 | CB≠G |
| 0.73 GeV H2 | 16 | 2.974 ± 0.020 | **3.416 ± 0.009** | +0.442 ± 0.019 | 3.714 ± 0.010 | 3.450 ± 0.011 | −0.264 ± 0.008 | CB≠G |
| Ta 25800 | 8 | 1.983 ± 0.036 | **2.028 ± 0.031** | +0.046 ± 0.028 | 2.018 ± 0.030 | 2.047 ± 0.031 | +0.028 ± 0.021 | |
| Ta 26138 | 4 | 1.949 ± 0.108 | **2.094 ± 0.046** | +0.145 ± 0.111 | 2.035 ± 0.053 | 2.114 ± 0.039 | +0.080 ± 0.051 | CB≠G |
| Ta 26230 | 13 | 1.981 ± 0.070 | **2.132 ± 0.055** | +0.152 ± 0.073 | 2.057 ± 0.054 | 2.126 ± 0.059 | +0.069 ± 0.047 | CB≠G |

- **Gated, CB R5nr:** 1.493 / 1.798 / 1.797 / 3.415 % (3.49 H2 / 2.24 H2 / C / 0.73), within 0.002 of legacy. Ta 1.951 / 2.040 / 1.964 %, Δ vs R4 +0.043 ± 0.038 / +0.149 ± 0.101 / +0.075 ± 0.053.
- **Preferred fit (χ² selection, information only), legacy:** 1.494 / 1.800 / 1.799 / 3.440 % (0.73: CB in 5 splits, Gauss in 11); Ta 2.028 / 2.094 / 2.148 % (26230: CB 8, Gauss 5).
- **Projection vs measurement.** The plan expected about −0.04 / −0.03 / −0.02 pp at 3.49 H2 / 2.24 H2 / C from the single-ingredient studies, and the f08 bars of the 10-28 deck projected −0.059 / −0.046 / −0.033 pp from seed 30 alone (20261028 §8). On CB the measured change is −0.000 / −0.006 / +0.011 pp; on the Gaussian it is −0.032 / −0.010 / −0.018 pp. §5 attributes the difference.
- Figure: `notes/figures/r5_20261029/resolution_fit_r5nr.png` (copy of `$O/figures/resolution_fit_r5nr.png`; legacy CB, R4 hollow, R5nr filled; the label "R5" there means R5nr).

### 3.1 Refinement-on reference arm and the cost of turning it off

Legacy CB, refined R5 (benchmark d91af9b7, refinement applied) and the cost Δ(R5nr − refined) (`$O/tables/r5nr_table.md`). The last column is refined − R4 (`promotion.md`):

| Point | n | CB refined R5 | Δ(R5nr − refined) | Δ(refined − R4) |
|---|---:|---|---|---|
| 3.49 GeV H2 | 27 | 1.472 ± 0.007 | +0.021 ± 0.003 | −0.021 ± 0.004 |
| 2.24 GeV H2 | 42 | 1.773 ± 0.005 | +0.027 ± 0.002 | −0.033 ± 0.003 |
| 2.24 GeV C | 8 | 1.774 ± 0.012 | +0.025 ± 0.004 (CB≠G) | −0.014 ± 0.010 |
| 0.73 GeV H2 | 16 | 3.374 ± 0.010 | +0.042 ± 0.003 | +0.399 ± 0.018 |
| Ta 25800 | 8 | 1.959 ± 0.053 | +0.070 ± 0.024 | −0.024 ± 0.048 |
| Ta 26138 | 4 | 2.004 ± 0.022 | +0.090 ± 0.033 | +0.055 ± 0.100 |
| Ta 26230 | 13 | 2.064 ± 0.050 | +0.069 ± 0.031 (CB≠G) | +0.083 ± 0.072 |

- Turning the refinement off costs +0.021 to +0.027 pp on e-p and C, +0.042 at 0.73 GeV and +0.07 to +0.09 pp on Ta. With the refinement on, R5 is below R4 by more than 1 SE on CB at 3.49 and 2.24 GeV H2 and C; 0.73 GeV is still worse (+0.399 ± 0.018), and Ta is within 1.2 SE of R4 on all three runs.
- **S1 residual, refine-off** (`$O/tables/window_gain_check_nr.md`; diagnostic only, no second pass under the contract): the legacy pass_abs criterion (|median − 1| ≤ 0.002 and rms ≤ 0.005) holds for 4 of 7 groups. It fails for 0.73 GeV H2 (median 1.0010, rms 0.0057), Ta 26138 (1.0025, 0.0051) and Ta 26230 (1.0015, 0.0051). pass_rel holds for all 7: no group is worse than R4's own residual by more than 0.001.

## 4. a/√E ⊕ b on the three H2 points

2-term σ²/E² = a²/E + b², weighted with the per-point SE (1 dof), E = 0.727 / 2.236 / 3.484 GeV, legacy (`$O/fit/out_r5.txt`, ±2σ mode). The gated fits agree with legacy to ≤ 0.007 in a and b. C and Ta are not fitted.

| Arm | Series | a (%·√GeV) | b (%) | χ² (1 dof) | 3-term a, b, c |
|---|---|---|---|---|---|
| R4 | CB | 2.484 ± 0.020 | 0.700 ± 0.031 | 3.58 | c² < 0 |
| R5nr | CB (headline) | 2.934 ± 0.011 | −0.709 ± 0.021 (b² < 0) | 202.98 | 2.11, 0.86, 1.60 |
| R5nr | Gauss | 2.858 ± 0.011 | −0.347 ± 0.039 (b² < 0) | 331.65 | 2.14, 0.92, 1.58 |
| R5nr | preferred (info) | 2.924 ± 0.013 | −0.692 ± 0.026 (b² < 0) | 210.59 | 2.08, 0.88, 1.65 |
| refined R5 (ref.) | CB | 2.888 ± 0.011 | −0.665 ± 0.022 (b² < 0) | 206.88 | 2.06, 0.86, 1.60 |

- **Every R5 2-term fit is unphysical.** R4's CB fit (b = 0.70 %) was physical. In R5nr the CB points at 2.24 and 3.49 GeV move by ≤ 0.007 pp from R4, while the 0.73 GeV CB point moves from 2.974 to 3.416 %; the 2-term form does not fit the three together (χ² ≈ 200). The Gaussian series is unphysical in R4 as well (b² < 0, χ² 587).
- **The 3-term fits have c² > 0**, but they have 3 parameters on 3 points (0 dof). They are an exact interpolation, not a test.
- The fit depends on the open 0.73 GeV question (§8). No a, b value for R5nr is quoted as a result.

## 5. Ablation on the 41 reference splits

Each variant is the R5nr command with one change (A3: config and maps). Δ = variant − R5nr, CB, legacy, pp; Δ > 0 means the ingredient helps R5nr (`$O/tables/ablation_nr.md`; 205/205 jobs rc = 0). n = 4 per e-p/C point and 8 / 4 / 13 for Ta. † = CB≠G.

| Variant (change vs R5nr) | 3.49 H2 | 2.24 H2 | C | 0.73 H2 | Ta 25800 | Ta 26138 | Ta 26230 |
|---|---|---|---|---|---|---|---|
| A1 Møller rejection off in the sample | +0.000 ± 0.000 | −0.008 ± 0.003† | −0.005 ± 0.002† | −0.156 ± 0.011† | +0.014 ± 0.009 | +0.001 ± 0.001 | −0.001 ± 0.005 |
| A2 seed 30 → 10 | +0.009 ± 0.003 | +0.012 ± 0.003 | +0.010 ± 0.002 | +0.016 ± 0.005 | −0.007 ± 0.006 | +0.010 ± 0.015 | −0.076 ± 0.038† |
| A3 log_weight 3.6 with R3 maps | −0.041 ± 0.005† | −0.015 ± 0.007 | −0.013 ± 0.003 | +0.008 ± 0.003 | +0.015 ± 0.008† | −0.025 ± 0.030 | +0.001 ± 0.013 |
| A4 legacy fit, single CB | +0.000 ± 0.000 | +0.000 ± 0.000 | +0.000 ± 0.000 | +0.000 ± 0.000 | +0.000 ± 0.000 | +0.000 ± 0.000 | +0.000 ± 0.000 |
| A5 SimpleProfile | +0.000 ± 0.002† | −0.004 ± 0.001† | −0.001 ± 0.001† | +0.003 ± 0.003 | −0.007 ± 0.016† | −0.034 ± 0.035 | −0.058 ± 0.034† |

- **Seed 30 helps,** but less than at R4 settings: +0.009 to +0.016 pp on e-p and C, against −0.033 to −0.059 pp validated on the odd splits at R4 settings (20261028 §8). Ta is mixed.
- **log_weight 5.2 costs** 0.041 pp at 3.49 GeV (Gauss 0.018; CB≠G), 0.015 at 2.24 GeV H2 and 0.013 on C. It is the largest ablated cost on e-p and C. A3 also swaps in the R3 maps (derived at 3.6 with seed 10 and R3-era gains), a second-order confound.
- **Møller rejection in the sample** matters only at 0.73 GeV, and there the two fits disagree. Removing it lowers CB by 0.156 ± 0.011 and raises the Gaussian by 0.234 ± 0.017 (Δ(CB − G) −0.390 ± 0.023).
- **Two-peak CB / χ² selection:** Δ CB = +0.000 at every point. With Møller rejection on, the single and two-peak CB give the same σ.
- **Geant4 profile at fixed R5nr gains and maps:** |Δ| ≤ 0.004 pp on e-p and C, and SimpleProfile is slightly better on Ta (−0.007 / −0.034 / −0.058). The −0.029 to −0.047 pp gain measured at R4 settings (`output/cluster_profile_20261028/report.md`) does not appear at R5 settings. Not investigated.
- **Consistency** (Σ of −Δ over A1–A5 vs the measured R5nr − R4 on the same splits, CB). The two agree at 2.24 GeV H2 (+0.014 ± 0.006 vs +0.014 ± 0.016), C, Ta 26138 and Ta 26230. They differ (> 2 SE) at 3.49 GeV (+0.032 ± 0.004 vs +0.003 ± 0.009), 0.73 GeV (+0.129 ± 0.022 vs +0.449 ± 0.029) and Ta 25800 (−0.015 ± 0.027 vs +0.046 ± 0.028). The steps not ablated (gain and map re-derivation, refinement off, runinfo) and interactions between ingredients carry the remainder; they are not separated. On the Gaussian at 0.73 GeV the sum is close but differs at 3 SE (−0.262 ± 0.021 vs −0.241 ± 0.017).

## 6. Gain choice (user decision 2)

Arm `r4gains` = R5nr on the 41 splits with the R4 gains (`--gain-file` = R4 tables). Δ CB (r4gains − R5nr gains), legacy: +0.003 ± 0.003 (3.49 H2), −0.008 ± 0.005 (2.24 H2), −0.009 ± 0.004 (C), +0.025 ± 0.005 (0.73), +0.024 ± 0.011 / +0.032 ± 0.008 / −0.014 ± 0.028 (Ta). Mean of the 7 per-point Δ: **+0.008 ± 0.005 pp** (SE = √(Σ SE_i²)/7). The rule fixed before the results adopts R4 gains only if the mean is < 0 and more than 1 SE below zero, so the **Møller-rejected R5 gains are kept** (`$O/tables/ablation_nr.md`, "Gain selection"). The refined `r4gains` arm (`$O/gainsel/r4gains`) was abandoned at 10/41 splits when the contract changed.

## 7. Promotion check and decision

| Criterion (plan S6) | Result | Basis |
|---|---|---|
| 1 S2–S4 verifications | PASS | 0 failures (dump 119, maps 118, S4nr 118); ε_peak within ±0.5 % on e-p/C; S1 residual reported, not converged (§3.1) |
| 2 chain reproduction | PASS | S0.10 identical to `b_maincfg` on 4 splits; refine-off rerun of 024590.00080 identical |
| 3 headline e-p/C ≤ R4 by ≥ 1 SE | DISCUSS | 2.24 GeV H2 passes on CB; 3.49 H2, C and 0.73 do not on CB but do on the Gaussian, and all are CB≠G, so they go to discussion (decision 1) |
| 4 Ta per run within 1 SE of R4 | FAIL | all six rows (3 runs × legacy/gated) above +1 SE (§3) |
| 5 a/√E ⊕ b physical, χ² ≤ R4 + 2 | FAIL | b² < 0, χ² 202.98 vs R4 3.58 (legacy); 205.30 vs 2.72 (gated) |
| 6 ablation table complete | PASS | 205/205, 41 splits per variant, consistency section present |

Under the plan, R5nr is reported, not promoted. The user's decision of 2026-10-29 makes R5nr the working reference regardless. R4 stays documented as the previous reference (summary §4). The presentation (`output/presentation_20261028/`) is not updated: S7 updates it only on promotion.

## 8. Open items

- **0.73 GeV: CB or Gaussian (open, not decided).** R5nr gives CB 3.416 ± 0.009 and Gauss 3.450 ± 0.011. Against R4, CB is +0.442 ± 0.019 and Gauss −0.264 ± 0.008 (Δ(CB − G) +0.705 ± 0.022). R4's CB (2.974) was a single CB fitted over the Møller peak, with χ²/ndf ~50 (20261028 §10). The ablated ingredients account for +0.129 ± 0.022 of the +0.449 ± 0.029 CB rise on the 4 reference splits; on the Gaussian the sum is close (−0.262 vs −0.241) but differs at 3 SE. The headline, criterion 3 and the a/√E ⊕ b fit all depend on this choice.
- **Ta regression.** CB worse than R4 on all three runs (legacy +0.046 ± 0.028 / +0.145 ± 0.111 / +0.152 ± 0.073). Part of it is the refinement (§3.1: +0.07 to +0.09 pp). The ablated ingredients do not account for 25800 (§5). Not resolved.
- **Integration-range scan:** running, output to `$O/intscan/` (no result tables yet when this note was made).
- **Per-run-type spectra for R5nr:** `$O/figures/spectra_r5nr.pdf` (PNGs in `$O/figures/spectra_r5nr/`); not yet reviewed in these notes.
- **S1 residual** at 0.73 GeV and Ta 26138 / 26230 (§3.1): reported, not acted on.

## 9. Sources

`$O/tables/{promotion,r3_r4_r5nr,r5nr_table,ablation_nr,window_gain_check_nr}.md` (with `.tsv` and `.log`), `$O/fit/out_r5.txt` (phase0.py, `$O/fit/phase0.py`), `$O/figures/resolution_fit_r5nr.png`, `output/cluster_profile_20261028/report.md`, `plans/plan_r5_20261028.md` (user decisions of 10-28 and 10-29, contract change of 10-29, S5–S7).
