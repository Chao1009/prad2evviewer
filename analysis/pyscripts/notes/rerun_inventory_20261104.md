# Rerun inventory under the runinfo calibration rule (2026-11-04)

**Scope.** This note covers every study that produced numbers and is recorded in `notes/*.md` and `plans/*.md`, together with the 2026-11-04 outputs. Each study is checked against two standing rules:

1. Calibration comes from the file that `database/runinfo/general.json` assigns to the run, with nl1/nl2 applied.
2. Results are not compared with old-calibration numbers. They are regenerated as absolute numbers.

This note is read-only and launches nothing. R5nr-ri, R7 and the f(A) ablation are planned separately in `plans/plan_ri_20261104.md`.

## 1. Verified facts that set the frame

- **No study ran with zero calibration.**
  - `grep -rF "(-1 modules"` under `analysis/pyscripts` finds only prose ("ring-1 modules"): `output/offline_20261025/A/study_a.log:286`, `plans/plan_20261025.md:563`, `plans/plan_ml_clustering_20261027.md:19`.
  - Every `[setup] HC calib` loader line in every log reads `(1728 modules)`.
  - Dumps carry nonzero `m_E`. For example, the `module_timing_20261030` dumps have a median of 29-34 MeV.
- **Why no zero calibration occurred.** Commit 165a00c (2026-09-16, which moved files to `old/`) is on `main`. It reached `waveform_jz` only through merge 2daedc9 (2026-10-28 13:41); `git merge-base --is-ancestor 165a00c 2daedc9^1` is false.
  - All jobs up to 10-28 read the top-level `*_fadc.json` files, which still existed. The move was a pure rename, so their content equals the `old/` files.
  - Jobs after the merge used explicit `old/` paths or the runinfo default.
  - The one post-merge missing-file case failed loudly, as a FileNotFoundError: `output/module_timing_20261030/failed_nocalib/run_vs_integral.log:48` (`.../calibration/old/x17_25825_fadc.json`). It produced no data and was moved aside.
- **Basis actually used by the energy-dependent studies of 10-21 to 11-04.** None of these files is in general.json:
  - `calibration_factor_fadc_3p5.json`, `calibration_factor_2p2_fadc.json` and `calibration_factor_fadc_0p7.json` (top level until 10-28, `old/` after).
  - `x17_{25825,26094,26214}_fadc.json` for Ta.
  - Representative headers:
    - `output/r5_20261028/scan/r5nr_splits/*/job.cmd` (`--hc-calib-file .../calibration/old/calibration_factor_fadc_3p5.json`).
    - `output/r5_20261028/gains/024503/gains_024503_all.tsv:3` (`calibration_file=.../old/calibration_factor_fadc_3p5.json sha1=be00ce41`).
    - `output/percrystal_spl/025320/gains_025320_all.tsv:3` (`...calibration_factor_2p2_fadc.json sha1=7ee32263`).
- **Run → general.json mapping** (`calib_resolve.resolve_calib`):

  | Runs | Beam | general.json file |
  |---|---|---|
  | 24386-24512 | 3.49 GeV | `calibration_factor_3p5_June7_nonLinear.json` |
  | 24590, 24603 | 0.73 GeV | same file, by chained inheritance from 24551 |
  | 25228 | 2.24 GeV | `PRad2_2p2_round1_1` |
  | 25257, 25260 | 2.24 GeV | `round1_2` |
  | 25299, 25308 | 2.24 GeV | `round1_3` |
  | 25320 | 2.24 GeV | `round1_4` |
  | 25800 | Ta | `x17_25825.json` |
  | 26138 | Ta | `x17_26094.json` |
  | 26230 | Ta | `x17_26214.json` |

- **nl1/nl2.**
  - Only `calibration_factor_3p5_June7_nonLinear.json` carries nonzero values: nl1 in 1006 modules and nl2 in 986. Its median correction factor is 1.06 at 0.73 GeV and 1.008 at 2.24 GeV (`output/nonlin_20261104/step1/june7_factors.txt`).
  - The `round1_*` and `x17_*` files have nl1 = nl2 = 0, so rule 1's "nl applied" changes numbers only for the 3.49 and 0.73 GeV runs.
  - The old files carry only a legacy `non_linear` key, which the C++ loader does not read (`prad2det/src/HyCalSystem.cpp` reads only nl1/nl2). No recorded study ever applied non-linearity.
- **Size of the change.**
  - Median per-module factor ratio, new over old: 0.954 (SD 5.9 %) for 3.49/0.73 GeV and 0.926 (SD 4.6 %) for 2.24 GeV, round1_3 vs 2p2_fadc.
  - `x17_*.json` over `x17_*_fadc.json`: about 0.94, with SD 3 %. Extremes are 0.53-2.4.
  - Absolute scales and per-module structure both move, so gains, position maps and bias tables fitted on the old basis do not transfer.

## 2. Inventory

Verdicts:

- (i) Rerun on the new basis; needed for a physics conclusion.
- (ii) Regenerate later, low priority.
- (iii) No rerun; the result does not depend on calibration.
- (iv) Invalid as recorded: a wrong calibration even by the rules in force at run time (wrong period, or a private file), or a statement that is false on the new basis.

All calibration-dependent rows not marked (iv) are on the old fadc basis, so rule 2 retires their numbers from comparison either way. "Dep" is Y when the per-module energy scale enters, m when it enters only mildly (a cut or a selection), and N when it does not.

Calibration abbreviations:
- **O** = the three `*_fadc` files: e-p/C runs at top level until 10-28 and under `old/` after.
- **Xf** = `x17_{25825,26094,26214}_fadc.json` per run.
- **RI** = runinfo default at run time.

| Study | Note § / file | Script + calib used | Dep | Verdict | Reason |
|---|---|---|---|---|---|
| Pre-10-21 scans `scan_025260/25299/25308/25320`, bimodality and per-module diagnostics | plan_waveform §4-7; `output/scan_025320/split_00010/benchmark.log:8` | benchmark_hycal_timing.py, scan_module_bimodality.py; RI = June7_nonLinear | Y | (iv) for the 2.24 GeV runs | 3.5 GeV file applied to 2.24 GeV runs; superseded |
| Same scans, 0.73 and 3.49 GeV runs (`scan_024386/…/24590`) | same | June7_nonLinear (RI) | Y | (ii) | File correct under rule 1 but superseded by every later chain |
| Pre-window sweep, 24386 | plan_waveform §3; `output/timing_gate_benchmark_24386.log:5` | sweep_pre_window_24386.py; private `physics_calib_out/.../run24386/calib_factor_iter1.json` | Y | (iv) | Not a database file; superseded by gate Phase 1 |
| LMS and elastic-window validation (8 splits, 25320) | 1021 §3, §5.1; `output/lmstest/r4_on_8s_explicit.tsv:8` | benchmark_test_20261021; O (2p2_fadc) | Y (σE); LMS drift N | (ii); drift (iii) | LMS factors are ratios; on/off deltas small |
| Full rawsum scan, Legacy vs New | 1021 §4 | scan driver `_HC_CALIB_MAP`; O + Xf | Y | (ii) | +0.02-0.04 pp contrast unlikely to flip; the +7.5 % offset is specific to the old file |
| Gate scans: Phase 1, Phase 4, gatescan_20261030, gatescan_tcal, ta/ep_gate_scan | 1022 §2, §7; summary §3.3 | gate_scan_2026102x.py, harness; O + Xf | m | (ii); Ta rows (i) at the R7 or Ta rescan | e-p/C gate effect under 0.004 pp; Ta rests on 9-34 modules |
| mf pulse validation (speed, noise-free, shape, κ, V9) | 1022 §3 | mf_pulse.py, mf_validate.py; none | N | (iii) | Pure waveform |
| mf vs rawsum resolution (Phase 3, full scan) | 1022 §4; `output/scan_20261022_mf/*/benchmark.log:8` | harness; O + Xf | Y | (ii) | Paired −0.04..−0.31 pp gain is robust; absolutes come from R5nr-ri |
| X17 templates | 1022 §5 | fit_pulse_template.py | N | (iii) | Waveform |
| X17 resolution | 1022 §5 | harness; Xf | Y | (i) via R5nr-ri Ta rows | Wrong X17 variant; 2-34 modules |
| mf time offsets, raw timing | 1022 §6 | mf_validate.py | N | (iii) | Timing only |
| Phase 0 decomposition (a, b, c), g_M RMS 2.2 %, decomposition refit | 1023 §2, 1025 §8 | phase0.py; O | Y | (i) | The constant term measures the quality of the old file |
| Per-crystal gains (S_dev, held-out, Round 2 pooled) | 1023 §3-4; `output/percrystal_spl/025320/gains_025320_all.tsv:3` | percrystal_calib.py; O (sha1 7ee32263) | Y | (i) inside R5nr-ri | Gains are relative to the old file; the −0.05..−0.09 pp lever may shrink |
| Seed-time neighbour amplitudes | 1023 §5 | benchmark; O | Y | (ii) | Failed (+0.015 pp); not driven by calibration |
| Free-τ fit study; template integral vs raw sum | 1023 §6, 6.1; 1025 §3.1 | mf_fit_study.py; O | σE arms Y; E/E_raw N | E/E_raw (iii); σE (ii) | Waveform effect |
| Pile-up examples and rates (+12-40 ns artefact), 10-23/10-24 | 1023 §7, §9.3-9.4 | plot_pileup_examples.py; O + Xf | m (E_dom > 50 MeV) | (ii) | Only the threshold is calibrated |
| Pile-up detection grid, synthetic efficiency (B/C/D arms) | 1026 §9.1, 1027 §6.3 | benchmark_pileup_mf.py, pileup_generator.py; none | N | (iii) | Synthetic ADC waveforms |
| pileup_eff_20261024(_fix), benchmark leg | `output/pileup_eff_20261024_fix/*/benchmark.log` | pileup_generator.py; RI: June7_nonLinear ×4 (including the 2.24 GeV H2 and C runs), x17_25825.json | m | (iv) for the 2.24 GeV benchmark leg; synthetic efficiency (iii) | 3.5 GeV file applied to 2.24 GeV runs |
| Module-by-module resolution, centred rectangle | 1023 §8a | plot_module_resolution.py; O (percrystal_rect, sha1 7ee32263) | Y | (i) | Worst-module lists depend on calibration |
| Spline templates: extraction, f(A), production pools | 1023 §9.1-9.2 | spline_template_extract.py; no energy cut | N | (iii) | Waveform |
| Sampled-template calibration A vs B; ±2σ Gaussian definition; f(A) options (i)-(iii) | 1023 §9.5; 1025 §2-3 | harness, regauss_tables.py; O | Y (absolutes) | (ii); f(A) diagnosis (iii) | Contrasts are waveform or method effects; the f(A) ablation is re-planned on R5nr-ri |
| Position correction R1/R2, hybrid map R3, Studies A/B | 1025 §4.1-4.9; 1026 §2 | position_corr.py, offline_eval.py; O | Y | (i) inside R5nr-ri | R(u,v) = E_cal/E′ trained on old-file residuals |
| Best-50 select/build/check | 1025 §4.10; 1026 §8.3-8.5; `output/templates_spline_best50_20261026/logs/select_026230_*.log` | spline_template_extract.py select; O + **x17_25825_fadc for all Ta** | m (E < 1.2·E_beam cap only) | (iii) for e-p/C shapes; **(iv) Ta 26138/26230 cap** | 26230: 2299 candidates rejected per split under the wrong-period file; Ta 29,710 above-cap candidates unexplained |
| p1-shape cut arms A/B | `output/templates_spline_best50_p1shape_20261104/report.md:3-4` | select + build; old/ explicit | m (cap only) | (ii) | The shape-cut effect (fallback-2 fraction, χ²) is a waveform property; the cap moves few candidates except for Ta |
| LMS for all runs, LMS gain variation | 1025 §5; 1028 §7; `output/scan_20261025_p3/lms_2run` | harness (O), lms_gain_variation.py (RI) | drift N; σE Y | (iii); σE (ii) | Null effect |
| Clustering variants cc1, mme | 1025 §6 | benchmark; O | Y | (ii) | Small, method-driven |
| **"non_linear_corr off = exact null"** | 1025 §6.1; summary §3.4; `output/nonlin_20261104/step1/findings_step1.md` §3; walkthrough S7 (:477); plan_mf_stdcal D9 (:76) | statement | Y | **(iv) as a general claim** | True only for files with nl = 0; June7_nonLinear (3.49 and 0.73 GeV) has nl1 ≠ 0 in 1006 modules |
| Software Møller veto (0.73 GeV), Møller diagnostics, Møller kinematic rejection, two-peak CB, CB vs Gaussian | 1025 §7; 1028 §10; summary §3.11-3.12 | benchmark; O (0p7_fadc ×16 at 0.73 GeV) | Y | (i) | 0.73 GeV moves to June7_nonLinear with a median nl factor of 1.06; E_ee/E′ and gap positions shift |
| Møller over-subtraction 0.73 GeV (σ_rule 3.0 % vs 3.45-3.68 %), CB excess | 1030 §4, §7 | offline refits of R4/R5nr; O | Y | (i) | Same reason |
| Module timing and RF (§9.1-9.6, veto alignment, comb), tcal pdf, DAQ memo | 1025 §9; 1026 §4; daq_queries | module_timing_dump.py (RUN_CFG); O, Ta → x17_25825_fadc | m (E > 10/100 MeV cuts) | (iii) for e-p/C | Timing; scale enters only via thresholds |
| RF per-run calibration (e100 tier), σ_t per module (B50 vs T0) | 1025 §9.7; 1026 §3, §9.4; `output/module_timing_20261030/dump/026230_01300.npz` (`calib_file=x17_25825_fadc.json`) | module_timing_dump.py, rf_time_calib.py; Ta 26138/26230 → x17_25825_fadc | m | (ii) for e-p/C; **(iv) Ta 26138/26230** | Wrong period; x17_26214/x17_25825_fadc per-module SD 8.7-10.8 %, max ratio 2.44; >10 % shift in 47 (26138) and 67 (26230) modules |
| σ_t per cluster | 1026 §9.4 | cluster_timing.py; O | Y (membership) | (ii) | Relative B50 vs T0 contrast |
| B50 adoption (118 splits), R4 promotion, best resolution per run type | 1027 §2, §7; 1028 §4; `output/resolution_fit_20261027_best50/r4_splits/*/job.cmd` | benchmark, percrystal_calib.py; O + Xf | Y | (ii) | History, superseded by R5nr(-ri); rule 2 forbids comparing to it |
| CFD-fraction scan | 1026 §7; 1027 §5 | cfd_fraction_scan.py; O, 25800 → x17_25825_fadc | m (E100/E500 bands) | (ii) | Band membership only; no 26138/26230 |
| Wassim single-pulse bank; two-pulse χ² ceiling U0; U5 real-data two-pulse; amplitude-binned shapes; hidden doublets; flat-top | 1026 §8.2; 1027 §6, §8, §11; 1028 §2, §5 | spline_template_extract.py, twopulse_*.py, template_*.py | N (doublet pools: cap only) | (iii) | Waveform/ADC |
| Wassim σE/E three-way | 1027 §6 | harness; O | Y | (ii) | −0.020 pp at 3.49 GeV needs a full scan anyway |
| Cluster regression V0 | 1027 §9 | cluster_regression.py on R3 dumps | Y | (ii) | ML on hold |
| Real-data pile-up rates | 1027 §10; `output/pileup_rates_20261027` (`[setup] run 26230 (sampled) calib x17_25825_fadc.json`) | pileup_rates.py @21f2d76; O, all Ta → x17_25825_fadc | MF pulse counts N; cluster E ≥ 50 MeV, E1 ≥ 100 MeV columns m | pulse rates (iii); E columns (ii); **(iv) Ta 26138/26230 E columns** | Wrong period |
| Baseline noise (noise_diag) | 1028 §6 | noise_diag.py (calib from the pileup_rates dump) | N; σ_t at E ≥ 100 MeV m | (iii); E-band numbers (ii) | ADC |
| Ta degradation, n200 fix | 1028 §3; `resolution_fit_20261028_ta/n200/splits/*/job.cmd` | benchmark; Xf | Y | (ii); superseded by Ta rows of R5nr-ri | — |
| Clustering scan (seed 30), merge check (log_weight 5.2) | 1028 §8-9 | cluster_cfg_wrap.py, merge_driver.py; O + Xf | Y (MeV thresholds) | (ii); log_weight is now R7 | Bit-identity of the merge check (iii) |
| energy_bias, leakage, Geant4 profile; R6 (C++ bias vs pos-corr) | 1028 §11; 1030 §8; summary §3.13; plan_r6 | eb_wrap.py, bias_table.py; O + Xf | Y | (i) | Bias tables are specific to the basis; feeds the R6 decision |
| **R5nr headline** (118 splits), gains, maps, a/√E ⊕ b refit (b² < 0) | 1029 §2-6; summary §1.4, §3.6, §4; walkthrough | benchmark via `r5_jobs.py`; O ×93 + Xf ×25 | Y | (i) = R5nr-ri | Working reference |
| R5 refine-on and ablations (seed 30, lw 5.2, Møller, two-peak, profile) | 1029 §5; plan_r5 | same | Y | (ii); re-test only what R7 or R6 need | — |
| Integration-range scan | 1030 §2 | intscan; O | N in effect | (iii) | Per-module constant absorbed by the gains |
| CB window | 1030 §3 | cbrange offline refit | Y | (ii) | Fit appearance only; nothing adopted |
| MF fit quality 20-50 ADC, waveform rates and examples | 1030 §5; `output/waveform_rates_20261029/logs` | waveform_good_bad.py; RI (correct per run, including 26138 → x17_26094, 26230 → x17_26214) | N | (iii) | ADC selection; already conforming |
| Crosstalk (bad-fit origin) | 1030 §6; crosstalk_figures_20261030.md; `crosstalk_20261029/cache/*.log` | xt_collect.py; RI (June7_nonLinear, round1_1, round1_3: conforming) | N; q56_energy and "1 % cluster energy" m | (iii); energy items (ii) | Seeds in ADC |
| tseed_dt (Ta 26230 width, gated 1.460 / 3.169) | `output/tseed_dt_20261030/tseed_dt.md:36` | benchmark R5 configuration; O + Xf | Y (member E) | dt statistics (ii); energy columns (i) with the Ta timing item | — |
| ML pile-up CNN M1/M2/U0/U5; M3 Ta data | ML_pileup_cnn_20261029.md; `M3/data/real/manifest.json` | twopulse_realdata.py, m2b/m3_decode.py; old/2p2_fadc and x17_25825_fadc (energies recorded only) | N | (iii) | Labels and inputs in ADC |
| mfstd V1 (constants; gate RMS 6.8/10.7/11.7 %) | mfstd_20261102.md V1; `output/mfstd_20261102/tools/run_points.py:18,39` | physics_calib seeded with `-c old/<seed>`; iteration 3 has nl = 0 | Y | (i) if mfstd is pursued | Re-seed from runinfo: round1_4 for 25320, June7_nonLinear with nl for 24590 |
| mfstd V2 (1.843/2.001/3.639 %; wastd) | mfstd_20261102.md V2 | benchmark with mfstd/wastd iteration-3 files, `--lms-corr off` | Y | (i) if mfstd is pursued | The wastd − mfstd contrast is internal to the arm; the R5nr/R6 comparison columns break rule 2 |
| Non-linearity step 2 (seed ratios), step 3 (E/E′ slopes, beam offsets) | `output/nonlin_20261104/step2,3/` | step2_ratios.py, step3_beamoffset.py on R5nr clusters (O) | Y | step 2 (i); step 3 offsets (ii) | Re-derive the offsets on R5nr-ri before applying |
| R7 W0(E) and f(A) ablation (stopped 11-04) | `output/r7_20261104/README.md`, `output/fA_ablation_20261104/README.md` (last line: "Run stopped 2026-11-04 (old/ calibration ...)") | r7_jobs.py and the fA generator inherit `r5_jobs.py:36-37,93-95` (rewrites to old/); 52 and 216 `result.tsv` | Y | (i) = R7 and f(A) on R5nr-ri | Stopped and incomplete, on the old basis. `w0.md` W0 values (3.7073/3.6576/3.4101) are a fit input and independent of calibration |

## 3. Open items and their basis

| Open item | Where | On the new basis? |
|---|---|---|
| 0.73 GeV CB excess (+0.705 ± 0.022 pp CB−G; +0.273 pp unattributed), CB vs Gaussian | summary §1.4, §3.14; 1029 §8; 1030 §4, §7 | Yes. Redo only after R5nr-ri, because nl at 0.73 GeV may change it |
| b² < 0 in a/√E ⊕ b; c² < 0 in the three-term fits | summary §3.6; 1026 §2; 1028 §4 | Yes, from R5nr-ri |
| member_e vs e_mev mismatch (2.6 % of members, up to 43 MeV) | 1030 §7 | Diagnose on R5nr-ri dumps; the gain-file path may be involved (see §4) |
| Ta regression; S1 residual at 0.73 GeV and Ta; Ta gate rescan; Ta gated vs legacy | 1029 §3.1, §8; 1027 §3; summary §5.1 | Yes |
| Ta 26230 timing (+1.99 ns, crate 3 slots 15/8/18) | 1026 §3; daq_queries §3 | The jump is pure timing (no). The E-tiered numbers and tseed_dt need runinfo (x17_26214). LMS slot 15 +110.9 % at 26230 (1028 §7) is consistent with a recalibration |
| Beam offset from the non-linearity study ((0.8-1.0, 1.0) mm at 2.24 GeV, etc.) | `nonlin_20261104/step3/beam_offsets.tsv` | Yes. Re-derive from R5nr-ri, then decide whether to apply |
| In-sample templates, gains and maps (holdout check) | walkthrough :109, :161, :205, :799; 1023 §10; methods §9 | Gains and maps: yes. Template-shape check: no |
| Gains not refit after position correction or per clustering variant; ring dependence; seed-30 promotion; Geant4 profile mechanism; log_weight 5.2 | 1025 §10; 1028 §8-11; 1029 §5 | Yes (log_weight is R7) |
| Energy-bias table placeholders; radians/degrees bug at `hycal_energy_bowl_shape.cpp:179-183,231-235` | plan_r6 D5; 1030 §8 | Tables: yes. Code bug: no, fix it regardless |
| LMS convention at 2.24 GeV (physics_calib gain_corr vs pre-merge blob) | plan_mf_stdcal D13 | Yes |
| Noise-model refit (larger model term, `nps`/`nps_ped`); normal-shaped bad fits | 1030 §6; summary §3.14 | No |
| Cross-crate crosstalk path; no Ta split in crosstalk | crosstalk note :189 | No. Energy-fraction numbers, if extended, use runinfo |
| mfstd V1 gate failure (MF/WA −2.2 % at low amplitude); MF phase bias ±0.4 %; W335 template; RF offsets in database/; 24386 veto; CFD production; U5 threshold; ML U3/U4/U7 | various | No |
| ML U6 (σE/E on 4 splits per type) | plan_ml_pileup :68 | Yes |

## 4. Prerequisites before any rerun

1. **`--gain-file` drops nl and bias_corr.**
   - Member hit energies are `cal_factor*adc` with no nl (`prad2det/include/HyCalSystem.h:153-157`).
   - nl1/nl2 and bias_corr are applied to the cluster energy only (`prad2det/src/HyCalCluster.cpp:458-472`).
   - The benchmark's `--gain-file` path rebuilds the cluster energy as Σ gᵢ·Eᵢ over member energies (`benchmark_hycal_timing_test_20261021.py:3225-3236`), so both corrections are lost.
   - Fix this before R5nr-ri at 0.73 and 3.49 GeV.
2. **`r5_jobs.py` still rewrites paths to `old/`** (`output/r5_20261028/tools/r5_jobs.py:36-37,93-95`). The R7 and f(A) generators inherit it. `percrystal_calib.py` is unchanged.
3. **`x17_*.json` vs `x17_*_fadc.json`.** general.json names the non-fadc files, which differ by about 6 % in median. This note assumes rule 1 means the non-fadc files.
4. **Runinfo tag check.** The benchmark's tag check exits 2 when the runinfo default for the 0.73 GeV runs (a "3p5" file) is used without an explicit file (`benchmark_hycal_timing_test_20261021.py:1426-1452`). Pass the resolved file explicitly or relax the check.
5. **Skipped jobs on rerun.** The harness `job_done` ignores the sha1 of the gain and calibration files (`gate_scan_20261022.py:828-858`). Write reruns to new output directories, or they are skipped silently.

## 5. Prioritized rerun list

Costs assume 118-119 jobs per full chain at about 7 min per job on 1 core, so about 14 core-hours per chain. module_timing_dump is about 111 splits, and its per-split time is not recorded here.

| P | Rerun | Scope | Cost (rough) |
|---|---|---|---|
| 0 | Prerequisites in §4: nl/gain order, `r5_jobs.py`, x17 variant, tag check, new output dirs | code and decisions | none |
| 1 | **R5nr-ri**: headline, per-crystal gains, position maps, a/√E ⊕ b refit | 1 dump chain + gain/map derivation + 1 final chain | 2-3 chains, about 28-42 core-h |
| 2 | **R7** W0(E), then the **f(A) ablation**, on R5nr-ri (replaces the stopped runs) | as planned | R7 1-2 chains; f(A) 216 jobs, about 25 core-h |
| 3 | **0.73 GeV Møller block**: kinematic rejection, two-peak CB, over-subtraction, CB-excess attribution | 16 splits at 0.73 GeV per arm | about 2 core-h per arm |
| 4 | **energy_bias / Geant4 profile / R6** on R5nr-ri gains | as plan_r6 | 1-2 chains, about 14-28 core-h |
| 5 | **Non-linearity steps 2-3**: seed ratios and beam offsets, from R5nr-ri clusters | offline | minutes |
| 6 | **Ta timing**: module_timing_dump + RF e100 tier + tseed_dt + pileup_rates E columns, for 26138/26230 with the runinfo file (the scripts now resolve it) | Ta splits only (17 of about 111), or all 111 | Ta-only is a small fraction of a full dump |
| 7 | Low priority (ii): gate scans, CFD bands, clustering seed 30, Wassim σE three-way, best-50 Ta cap, p1-shape arms, pile-up E columns, mfstd V1/V2 (only if mfstd continues; re-seed from runinfo with nl) | as needed | each up to 1 chain |
| — | No rerun (iii): mf validation, templates and f(A) shapes, synthetic pile-up, CNN, noise, crosstalk, waveform rates, LMS drift, integration range, pure timing | — | — |

## 6. Coverage gaps

- Several output trees are no longer on disk: `scan_20261021_full`, `gate_scan_2026102x`, `scan_20261024`, `scan_20261025` (R2/R3 job.cmd), `module_timing_20261025` and `rf_calib_20261026`. Their calibration is inferred from the unchanged `_HC_CALIB_MAP` (md5-identical from 315e4af to dae28ce) and from the R4 job.cmd.
- The per-split logs of `pileup_20261024` were not kept. Its calibration is inferred from the script at 80e19b6.
- The provenance of `dt_diagnostic_025320`, the `timing_gate_025308` sources, and the starting file of physics_calib for the 25320 investigation (plan_waveform §2) are not recorded.
- It is unchecked whether `physics_calib` overrides its `-c` seed with runinfo.
