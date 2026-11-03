# PbWO4 energy resolution: integration range, CB window, Møller over-subtraction, and MF fit quality of small pulses

**Date:** 2026-10-30
**Branch:** `waveform_jz` (HEAD `71e049f`; the `waveform_*.py` scripts of §5 are untracked)
**Previous notes:** `notes/PbWO4_energy_resolution_20261029.md` (R5nr reference, §5 ablation, §8 open items), `notes/PbWO4_energy_resolution_20261028.md` (§10 Møller), `notes/ML_pileup_cnn_20261029.md`

Paths are relative to `analysis/pyscripts/` in `/Users/L00338853/work/PRad/prad2evviewer` unless stated; `$O` = `output/r5_20261028`. R5nr is the working reference of 20261029 (gain refinement off). Δ in pp, paired per split, mean ± SE. "MF" is the matched-filter pulse fit (`mf_pulse.py`); "good" fit = χ²/dof < 2, "bad" fit = χ²/dof > 3. No code, configuration or database change was made in this work.

## 1. Status

- **Nothing adopted.** User decision, verbatim: "Don't adopt anything". R5nr stays the working reference with refinement off, int_tail_ratio 0.1, ped floor off and the current CB window.
- **Integration range (§2):** moot for the mf method. Changing int_tail_ratio rescales each module by a constant that the per-crystal gain absorbs; the scan was stopped as redundant.
- **CB window (§3):** narrowing it fixes χ²/dof at 0.73 GeV (186.67 → 4.84) but moves σ_E/E by ≤ 0.023 pp. A patch exists; it is not applied.
- **Møller rejection at 0.73 GeV (§4):** 51–65 % of the Gaussian gain over R4 is a cut artefact (the rule over-subtracts the elastic peak at the window edge). The Møller rate shows no excess, and the gains are not the cause. The cut widens CB by +0.168 pp; a further +0.273 pp of the CB rise over R4 is unattributed.
- **MF fit quality, 20–50 ADC (§5):** 58.4 % good, 33.9 % bad over the 119 R5 files; bad 43 % at 3.49 GeV and 20 % at 0.73 GeV.
- **Origin of bad fits (§6):** about a quarter is crosstalk that stays outside clusters. About 60 % are normal-shaped in-time cluster members, most likely a noise-model or template-shape mismatch (inferred, untested). Do not cut on MF χ².

## 2. Integration range of the mf template

For the mf method E = A·I_int·κ, and int_tail_ratio r enters only through I_int.

- **Per-module constant.** Each bank holds one template per module and one global f(A) with no amplitude axis, so I_int(r) is a constant per module (`$O/intscan/tables/iint_offline.md:34`). I_int(r)/I_int(0.1) = 1.135 / 1.084 / 1 / 0.885 at r = 0.03 / 0.05 / 0.1 / 0.2, identical across banks to 0.1 % (module-to-module RMS 0.4–0.7 %; `iint_offline.md:45`). Re-derived per-crystal gains absorb the factor. With the fixed r = 0.1 gains, the scan would instead measure an artificial inter-module miscalibration (`$O/intscan/tables/intscan.md:8`).
- **Ped floor on vs off:** |Δσ| ≤ 0.0005 pp at every point (largest Ta 25800 CB, +0.0005 ± 0.0003), none above 2 SE, 41/41 splits rc = 0 (`intscan.md:9,54`; `$O/intscan/compare/pedfloor_vs_r5nr.md`).
- **Stopped.** "the redundant variants were stopped" (coordinator, 2026-10-29; `intscan.md:3`). The r003 and r020 arms finished but are marked not to be quoted. Defaults kept: int_tail_ratio 0.1, tail_break_n 2, ped floor off.

## 3. CB fit window

- **Current window** (`benchmark_hycal_timing_test_20261021.py:311-312`): [μ0 − w, μ0 + w/2] with w = max(400·kS, 12·bw) MeV, where kS = √(E_beam/2239.51) and bw = 20·kS MeV is the histogram bin. μ0 is the fixed-window Gaussian mean.
- **Narrowing** (per-run-type window in units of σ, chosen by the pre-registered criterion in `$O/cbrange/criterion.txt`): at 0.73 GeV the median χ²/ndf goes from 186.67 to 4.84 (a = 5, b = 5), and the CB n parameter is no longer pinned in any of the 16 splits (pinned in all 16 before). σ_E/E moves by ≤ 0.023 pp at every point (largest 2.24 GeV C, −0.023 ± 0.002; 0.73 GeV +0.007). The report concludes "The window is not a lever on the resolution" (`$O/cbrange/tables/cb_range.md` §3, §4, §8).
- **Where the 0.73 GeV χ² comes from:** 94 % from the low-energy shoulder below the Møller dip and 1 % from the dip bins. Masking the dip leaves χ²/ndf at 185.4 (`cb_range.md` §7). The window change does not resolve the CB-vs-Gaussian disagreement (6 of 7 groups).
- **Patch, not applied:** `$O/cbrange/proposed_benchmark.patch` adds `--cb-window-sigma A,B` (default off, current behaviour kept). It was validated by code reading and a patched-copy check, not by a logged run. The benchmark script is unchanged. The aim was fit appearance, not resolution, and the user decided not to adopt it (§1).
- Figures: `$O/cbrange/figures/chi2_heatmaps.png`, `overlays_before_after.png`, `heat_<run type>.png`.

## 4. Møller over-subtraction at 0.73 GeV

Sample: 0.73 GeV H2 inner ring, runs 24590 and 24603, 16 splits, compared with R4 on the same splits. All fits are offline refits of the stored histograms, which reproduce the stored values to ≤ 0.0012 pp. `M` = `output/moller_oversub_20261029`.

**The rule's width is too narrow.** `--moller-reject kin` rejects clusters within 3σ_rule of E_ee(θ) and more than 2σ_rule from E'(θ). σ_rule = E'·√((2.5 %/√E)² + (0.7 %)²) ≈ 3.0 % at 0.73 GeV (`benchmark_hycal_timing_test_20261021.py:253-254`). The measured Gaussian σ/E is 3.45 % with the cut and 3.68 % without. The cut therefore reaches into the elastic peak: it removes 5.51 % of elastic events (1.27–7.21 % across decomposition models), against ≤ 0.71 % at every other run type (`M/tables/kin_summary.md`). Inside the V1 Gaussian window it removes 203 clusters per split, at the low edge (window 676.2–776.0 MeV vs 672.1–778.6 MeV without the cut).

Rule variants, Gaussian and CB σ/E in % (`M/tables/refit_B_variants.md`, `M/tables/refit_V7_summary.tsv`):

| Variant | Gauss | CB (benchmark path) |
|---|---|---|
| R4 | 3.714 ± 0.010 | 2.974 ± 0.020 |
| V0 no rejection | 3.677 ± 0.014 | 3.248 ± 0.013 |
| V1 R5nr rule (n 3, m 2) | 3.450 ± 0.011 | 3.416 ± 0.009 |
| V6 R5nr rule with σ_rule = 3.45 % | 3.552 ± 0.012 | 3.500 ± 0.010 |
| V7 Møller template subtraction, no cluster cut | 3.584–3.622 (nominal 3.622) | 3.249–3.299 (CB1) |

- **51–65 % of the Gaussian gain is a cut artefact.** R4 → V1 gains 0.264 pp. V7 removes Møller events statistically, without cutting clusters, and lands 0.172 pp (nominal) or 0.134 pp (σ_M × 1.2) above V1, so 65 % and 51 % of the gain is over-subtraction. The σ_M × 0.85 configuration (0.194 pp, 73 %) is left out because the fitted Møller peak is 21–25 MeV wide against a template of 18–20 MeV, i.e. wider, not narrower (`M/tables/refit_V7_peakcheck.tsv`). The real Møller-removal gain is V0 − V7 = 0.054 pp (0.054–0.092).
- **No rate excess.** 144,854 rejected clusters against 151,508 Møller events expected from the cross section (ratio 0.956; 0.744–1.031 across fit variants). At 2.24 GeV H2 the ratio is 1.219 (`M/logs/kin_rate.log`).
- **CB:** the cut widens CB by +0.168 ± 0.007 pp (V1 − V0); the stored A1 ablation agrees (−0.156 ± 0.011, 20261029 §5). The remaining +0.273 ± 0.017 pp of the +0.441 ± 0.019 pp CB rise over R4 (V0 − R4) is not attributed. V0 uses the two-peak CB, while R4 is a single CB fitted over the Møller peak.
- **Gains are not the cause.** Per-module ratio of Møller-rejected to unrejected gains: mean 1.0000 / 1.0004, sd 0.0018 / 0.0019 (24590 / 24603). Applied offline, it changes Gauss by +0.008 ± 0.003 and CB by +0.007 ± 0.002 pp (`M/tables/gains_propagate_summary.tsv`).
- **Partner tag:** only 4.5 % of rejected clusters carry a partner tag. The partner tag alone (V5) gives CB 2.840 %, below R4.
- Tools in `M/tools/` (`kin_*.py` decomposition and rate, `refit_*.py` variants, `gains_*.py`); figures in `M/figures/` (`refit_hist_overlay`, `refit_sigma_by_variant`, `refit_gains_ratio`, `kin_q1`–`kin_q4`). No report file was written.

## 5. MF fit quality of small pulses (20–50 ADC)

**Selection and fit** (`waveform_good_bad.py:9-12,33,73,103`; `waveform_examples.py:29-36`). Pulses with 20 ≤ h < 50 ADC and h ≥ 10 × pedestal RMS, one clean peak, RawSum (SSP0) trigger, exactly one HyCal cluster and fewer than 20 fired channels. The fit is the `mf_pulse` matched filter with the run type's per-module spline templates (Ta: n200 bank), free amplitude and t0, and pedestal fixed from samples 0–29. The noise model is σ² = σ_ped² + (0.03·model)².

**Rates** over all 129 evio files, 15 runs (`output/waveform_rates_20261029/tables/rates.md`). The "ALL" rows cover the 119 files of the 14 R5 runs; run 24578 (non-standard trigger, 734 pulses) is excluded.

- **ALL:** 20,979,695 pulses; good 58.40 %, bad 33.85 %, between 7.74 %, no fit 0.01 %.
- **Bad fraction by run type:** 3.49 GeV H2 43.07 %, 2.24 GeV C 36.95 %, 2.24 GeV H2 32.65 %, Ta/X17 31.48 %, 0.73 GeV H2 20.36 %.
- **Weak amplitude dependence:** pooled bad 34.65 / 33.71 / 31.50 % in 20–30 / 30–40 / 40–50 ADC, not monotonic in most run types.
- **Stable over files:** file-to-file RMS of the bad fraction 0.06–0.30 pp for H2 and C runs, 0.73–0.99 pp for Ta. No run shows a trend with |t| > 3 (largest 2.7, run 24503).
- **Per module** (842 modules with ≥ 500 pulses): median 29.4 %, 5–95 % range 16.6–53.7 %. Worst W693 83.6 %, W464 77.4 %, W761 71.7 %, W498 70.5 %, W600 70.3 %. Maps: `output/waveform_rates_20261029/figs/bad_fraction_maps.pdf` (one page per run type).
- Uncertainties in `rates.md` are binomial only. Pulses are correlated within events, so the file-to-file columns give the realistic spread.

**Examples** (`output/waveform_examples_20261029/`): `<run>_<split>.pdf`, one split per run, up to 4 waveforms with the fit per module. `good_bad_by_module.pdf` has 2 good and 2 bad examples per module from 28 splits (`splits.txt`, `splits_extra.txt`); its pooled rates (58.3 / 34.0 / 7.7 %) agree with the full scan.

**Scripts:**
- `waveform_examples.py EVIO out.pdf --mf-templates NPZ`: per-module example pages.
- `waveform_good_bad.py {collect,rates,render}`: classifier, counts and the good/bad PDF.
- `waveform_rates_launch.sh [P]`: resumable launcher over all splits, each run with its own template bank.
- `waveform_rates_aggregate.py`: builds `rates.md` and the maps.
- `waveform_rates_validate.py`: cross-check against `gb_cache`. No recorded run output.

## 6. Origin of the bad fits

Sample: four splits, one per e-p/C run type (24590_00010 0.73 GeV H2, 25228_00050 2.24 GeV H2, 25308_00050 2.24 GeV C, 24386_00010 3.49 GeV H2), pulses of 20–200 ADC, same selection as §5. No Ta split. Tables and figures: `output/crosstalk_20261029/{tables,figs}/` (`X` below). The study's write-up exists only as a draft outside the repo; the numbers here are from the tables.

Crosstalk here means a signal coupled into a channel from another readout channel; it is distinct from external electronic pickup and from shower leakage (real energy in a neighbouring crystal), and the study's tables and scripts in `output/crosstalk_20261029/` call it "pickup".

Decomposition of 20–50 ADC bad fits, pooled (bad 34.0 % of 976,750 pulses; `X/tables/q13_tables.txt`):

| Class | Share of bad |
|---|---|
| Strict bipolar (crosstalk) | 26.3 % |
| Narrow or dip, not strict bipolar | 7.5 % |
| Late (dt ≥ 6 ns), normal shape | 3.3 % |
| Normal-shaped | 62.8 % |

**Bipolar crosstalk (~26 %).**
- **Signature:** bipolar shape arriving at dt_seed +11.6 ns (median; 5–95 % at 6.9–15.2 ns) after a large seed (median 1879 ADC). Bipolar-rich bad non-members have median FWHM 8.7–11.5 ns against 22.7 ns for good hits (`w2_t4.txt`; the tag requires < 15 ns). The bipolar share of all hits is 0.14 % at 50–60 ADC and ≤ 0.04 % above.
- **By run type:** 0 % of bad fits at 0.73 GeV, 24.3 % (2.24 GeV H2), 25.3 % (C), 43.2 % at 3.49 GeV. The 43 % is the crosstalk share of bad fits; the 3.49 GeV bad fraction itself is 43.7 %.
- **Electronics:** enriched on the adjacent channel of the seed's FADC (odds ratio 4.84, `X/tables/q24_b_slot.txt`; 3.4–6.0 across the sub-tables in `q24_d_matched.txt`). The adjacent-pair coupling is k = 1.0–1.4 % of the seed height (median 1.37 / 1.20 / 1.04 % by seed-height bin; `w2_t2c.txt`).
- **Cross-crate:** 74 % of bipolar hits are not on an adjacent channel, and 98.6 % of those are on another crate. They concentrate on specific slot pairs, led by crate 3 slot 8 → crate 4 slot 10 (12,730 hits), then 3/5 → 4/8 and 3/18 → 2/3 (`w2_t1*.txt`). The coupling path is not identified (medium confidence), and the cross-crate k was not measured.
- **Energy impact:** no cluster member lies more than 4 ns from the seed, so the late crosstalk almost never enters a cluster. The 174 bipolar bad members carry 5.1 × 10⁻⁶ of the cluster energy. In-time bad members on the seed's adjacent channel carry 0.10–0.29 % per split, an upper bound on in-gate crosstalk (`w2_t4.txt`).

**Normal-shaped in-time members (~60 %).**
- **Class size:** excluding the in-time crosstalk above, the study estimates about 55–60 %.
- **Behaviour:** bad members have amp/h 0.99–1.03, so the MF amplitude agrees with the pulse height. The pooled bad fraction stays at 31–34 % from 50 to 200 ADC, and the median χ²/dof of good fits rises from 0.80 to 1.19 with amplitude (`q13_tables.txt`, `q56_amp.txt`).
- **Inference, not tested:** the 3 % model term in the noise model is too small, or the template shape mismatches at larger amplitude. A refit with a larger term would test it.

**Consequences.**
- **Do not cut on MF χ².** Bad members carry 0.84–1.49 % of the cluster energy at 20–50 ADC and 3.05–5.07 % at 50–200 ADC per split, and they look like real energy (`q56_energy.txt`).
- **Correction to an earlier statement in this study:** "bad fits are ~40 % low in amplitude" holds only for out-of-cluster crosstalk. Bad non-members have amp/h 0.50–0.60 (1.00 at 0.73 GeV), bad members 0.99–1.03, and all bad hits together a median of 0.93.

## 7. Open items

- **0.73 GeV CB excess:** +0.273 ± 0.017 pp of the CB rise over R4 is not attributed (§4). Needs an ablation at 0.73 GeV.
- **0.73 GeV CB vs Gaussian:** still open (20261029 §8). §4 shows the Gaussian gain there is mostly a cut artefact.
- **Noise model:** refit with a larger model term to test the normal-shaped bad fits (§6).
- **Crosstalk:** cross-crate coupling path and cross-crate k (§6).
- **Coverage:** the §6 study has no Ta split.
- **Data consistency:** 2.6 % of cluster members (32,963 of 1,268,754) have member_e and e_mev differing by more than 10⁻³ MeV, up to about 43 MeV (`q56_amp.txt`, `q56_energy.txt`). Cause not investigated.

## 8. R6: C++ energy-bias correction in place of the Python position map

R6 is the R5nr recipe (refine off) with `--pos-corr` removed and the C++ `energy_bias_correction` enabled per job through `output/energy_bias_20261028/tools/eb_wrap.py --carry --energy-bias`. No global default, `database/` file or C++ file changed. Plan: `plans/plan_r6_20261030.md`; outputs: `output/r6_20261030/` (tables in `tables/`). Evaluation: 41 reference splits, Δ = arm − R5nr in pp (CB core σ headline), paired per split, mean ± SE, Ta per run.

**Arms.** `nopos` (R5nr minus `--pos-corr`), `r6db` (repo `database/energy_bias` tables), `r6rd` (tables re-derived on the R5nr gains, zero mean per module; the R6 candidate), `r6rd_ho` (same, derived without the 41 reference splits). Ta uses the pooled 2.24 GeV H2 grid.

**Regression (S1).** R1: the carried-config path without the bias reproduces `nopos` bit for bit (`result.tsv` and `result_clusters.tsv`, 4 splits). R2: with a uniform −0.01 table the raw cluster energy rises by 1/0.99 for all matched clusters (range 1.010095–1.010107). The selected-cluster count changes by −0.1 to −0.4 % because the factor also enters the in-loop selection on `r.hit.energy`.

**Result.** R6 (`r6rd`) is not promoted. The CB Δ(r6rd − R5nr) is positive at all four e-p/C points (legacy and gated): 3.49 +0.051 ± 0.003, 2.24 H2 +0.026 ± 0.004, C +0.027 ± 0.002, 0.73 +0.004 ± 0.005 pp. Ta 025800 is +0.012 ± 0.006 pp (legacy); 026138 and 026230 are within 1 SE. P4 (harness eps_peak within 0.5 %) holds. The C++ bias reduces σ relative to no correction (r6rd − nopos −0.090 ± 0.016 at 3.49, −0.027 ± 0.004 at 2.24 H2) but does not recover the Python hybrid map (nopos − R5nr +0.141 ± 0.015 at 3.49). The repo tables (`r6db`) are worse than `r6rd` at 3.49 and 2.24 H2 (+0.065, +0.083). The in-sample advantage r6rd − r6rd_ho is −0.018 ± 0.006 at 3.49 and −0.009 ± 0.002 at 2.24 H2. Only 52 (3.49), 78 (2.24 H2) and 32 (C) modules have a cell with n ≥ 200, so most modules use the pooled grid. Full tables: `output/r6_20261030/tables/promotion_r6.md`.

**Tool bug, not fixed here.** `analysis/tools/hycal_energy_bowl_shape.cpp` lines 179-183 and 231-235 pass `atan2(r, z)` in radians to `PhysicsTools::ExpectedEnergy(theta_deg, …)`. The R6 tables do not use that tool.
