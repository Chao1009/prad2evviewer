# PbWO4 energy resolution from waveforms: summary of the 2026-10-21 to 10-24 work

**Author:** Jingyi Zhou (jyzhou), documented 2026-10-24; §3–§6 updated from `plans/plan_20261025.md` (progress log to 10-25 01:42).
**Branch:** `waveform_jz` (HEAD `80e19b6`; everything after it is uncommitted, §5.2)
**Plans:** `plans/plan_benchmark_hycal_timing_test_20261021.md`, `plans/plan_20261021.md`, `plans/plan_20261022.md`, `plans/plan_20261023.md`, `plans/plan_20261024.md`, `plans/plan_20261025.md`
**Previous notes:** `notes/PbWO4_energy_resolution_20261021.md`, `notes/PbWO4_energy_resolution_20261022.md`, `notes/PbWO4_energy_resolution_20261023.md`, `notes/PbWO4_energy_resolution_20261025.md`

Paths are relative to `analysis/pyscripts/` in `/Users/L00338853/work/PRad/prad2evviewer` unless they start with `docs/`, `database/`, `prad2det/` or `prad2dec/` (repo root). "rawsum" is the WaveAnalyzer pulse integral; "mf" is the matched-filter template pulse (20261022 §3). σ/E is the post-refinement σ_E/E, inner ring r < 200 mm, preferred fit (§1.3). Legacy and gated ("New" in the 20261021 note) are the clustering paths of §1.3. Δ is in pp; "paired" means mean ± SE of per-split differences. Daily notes are cited as "20261021 §n".

## 1. Purpose and scope

### 1.1 Question

How far can waveform-level reconstruction and calibration push the elastic e-p σ_E/E of the HyCal PbWO4 inner region?

### 1.2 Data

| Run type | Runs (splits in `output/scan_20261024_pc`) | Beam energy (nominal) | Calibration file (`database/calibration/`) | Preferred fit |
|---|---|---|---|---|
| 0.73 GeV H2 | 24590 (10), 24603 (6) | 0.73 GeV | `calibration_factor_fadc_0p7.json` | Gaussian |
| 2.24 GeV H2 | 25228 (11), 25257 (7), 25260 (8), 25299 (8), 25320 (8) | 2.24 GeV | `calibration_factor_2p2_fadc.json` | Gaussian |
| 2.24 GeV C | 25308 (8) | 2.24 GeV | `calibration_factor_2p2_fadc.json` | CB |
| 2.24 GeV Ta (X17) | 25800 (8), 26138 (4), 26230 (13) | 2.24 GeV | `x17_25825_fadc.json`, `x17_26094_fadc.json`, `x17_26214_fadc.json` (unverified) | CB |
| 3.49 GeV H2 | 24386 (11), 24503 (10), 24512 (7) | 3.49 GeV | `calibration_factor_fadc_3p5.json` | CB |

Run 24578 is excluded (non-standard trigger; `_EXCLUDED_RUNS` in the driver). Split 24386.00060 timed out in the rawsum full scan and is missing there (10 splits, 118 in total); the mf scans have it (119). Source: 20261021 §2.1, §4.1, §6; `plans/plan_benchmark_hycal_timing_test_20261021.md` §3.

### 1.3 Figure of merit

- **σ/E**: post-refinement σ_E/E of the calibrated pool (single-cluster events, seed r < 200 mm, seed E ≥ 200·kE MeV with kE = E_beam/2239.51, seed outside `--exclude-rowcol`), mean over splits, preferred fit of §1.2 (20261021 §4.2; D5 of `plans/plan_20261023.md`).
- **Legacy path**: one hit per module, the largest-integral pulse inside the pre-window, no time cut.
- **Gated path**: every pre-window pulse is a hit; a neighbour joins the cluster only within |t_hit − t_seed| ≤ W (`--window`).

## 2. Road map and script structure

`EVIO → pulse (rawsum | mf) → LMS → time offsets (optional) → clustering (legacy + gated) → per-crystal gains (optional) → seed-indexed refinement → σ/E per split → driver / harness`

| Stage | Script | Purpose | Key flags | Output | Status; git |
|---|---|---|---|---|---|
| Reference | `benchmark_hycal_timing.py` | Original benchmark; bit-identity reference | `--cal-method {median,cb-peak}` | `out_path` | reference; unchanged since 5359d63 |
| Benchmark | `benchmark_hycal_timing_test_20261021.py` | All per-split stages below → σ/E | `evio_path out_path --calibrate --hc-calib-file` | `<out>.tsv`, `_window.tsv`, `_clusters.tsv` | in use; modified |
| Pulse | benchmark (WaveAnalyzer); `mf_pulse.py` (library) | rawsum; mf with h50 two-τ or sampled banks | `--pulse-method {rawsum,mf} --mf-templates --mf-r-model 0.3 --mf-scale 1.010` | `output/pulse_templates_h50.dat` | in use; modified |
| Pulse | `spline_template_extract.py` | Per-module spline templates, f(A) | `select`, `build`, `ampcorr`, `check` | `output/templates_spline/<slug>/` | Phase 3; untracked |
| Pulse | `mf_fit_pulse.py` (library) | Free-τ refit | `--mf-shape {fixed,bounded,prior}` | `output/mf_fit_study_20261023/` | diagnostic; untracked |
| LMS | benchmark; C++ `replay_gainCorr` → `prad_<run>_gain_corr.root` | Per-batch module gain, applied to each pulse integral before clustering | `--lms-corr {auto,off,PATH}` | | in use; null effect |
| Timing | benchmark; offsets from C++ `hycal_module_time_calib` | Per-module offsets (`database/hycal_time_offsets/time_calib.json`) | `--time-calib {off,auto,PATH}` (default off; driver never passes it) | | harness only |
| Clustering | benchmark (`HyCalCluster`) | Legacy and gated, both always run | `--window` (default 8; working point 16), `--pre-window LO HI` | | in use |
| Clustering | benchmark | Constant overrides on both clusterers (item 4) | `--cluster-min-module-energy MEV`, `--cluster-corner-conn {0,1}`, `--cluster-non-linear-corr {0,1}` (default: config value) | `# cluster_cfg` header line when given | P4 variants; modified |
| Calibration | benchmark | Seed-indexed elastic-window gains to E'(θ) | `--cal-method elastic-window --exclude-rowcol --inner-radius 200 --min-per-module 30` | `_window.tsv` | in use |
| Calibration | `percrystal_calib.py`; benchmark | Pooled per-crystal gains | `fit --dumps --splits --fixed-windows --out`; `--cluster-dump auto`, `--gain-file` | `output/percrystal_rect/<run>/` | in use; modified |
| Calibration | benchmark | Position correction E → E/R(u, v), after `--gain-file`, before the window fit, both paths | `--pos-corr FILE` (default off) | `# pos_corr` header line | in use (R2); modified |
| Calibration | `position_corr.py` | Position map R(u, v) per run type; bilinear and poly8 fits; global-vs-type test | `map --scan-dir --splits --path --gauss-window`; `fit --map --form {bilinear,poly8,both}`; `fit --global --maps` | `output/poscorr_20261025/` | in use (R2); untracked |
| Calibration | benchmark | Software Møller veto; Møller-gap re-admission | `--veto-second-cluster EMIN_MEV`, `--veto-dphi DEG` (needs EMIN), `--moller-gap-nsigma X` (default 3; 0 re-admits) | veto count in the result header | add-on study; modified |
| Fit | benchmark; `compare_scan_dirs.py`, `plot_module_resolution.py`, `position_corr.py map` | Reported Gaussian: iterative ±2σ or the old fixed ±150 MeV × kS window (§3.7) | `--gauss-window {2sigma,fixed150}` (`compare_scan_dirs.py` also `stored`) | `# gauss_fit` header line | default 2sigma since 10-24 19:50; modified (`position_corr.py` untracked) |
| Driver | `benchmark_hycal_timing_scan.py` | All runs and splits | `--pulse-method --window --pre-window --gain-file-map auto --gain-dir --mf-templates-map --write-cluster-tsv --output-dir` | `output/scan_<tag>/` | in use; modified |
| Harness | `gate_scan_20261021.py`, `gate_scan_20261022.py` | Gate scans (frozen); rawsum/mf and variant scans | `run --stage --half --configs-file --tags --bench-args`, `select`, `restable` | `output/gate_scan_2026102{1,2}/`, `output/scan_2026102{3,4}/` | 5359d63, unchanged |
| Diagnostics | `plot_module_resolution.py` | Per-module σ_M/μ_M maps | `--scan-dir --group-runs` | `output/module_resolution_20261024/` | untracked |
| Diagnostics | `plot_pileup_examples.py`, `pileup_artefact_study.py` | Pile-up rates, maps; +12–40 ns excess | `--runs --burst-veto --window-rule`; `{scan,…,all}` | `output/pileup_20261024/` | untracked |
| Diagnostics | `pileup_generator.py`, `benchmark_pileup_mf.py` | Synthetic pile-up, mf efficiency | `--template`; `--templates --window-rule` | `output/pileup_eff_20261024/` | modified; untracked |
| Diagnostics | `mf_validate.py`, `mf_fit_study.py`, `compare_scan_dirs.py` | Pulse checks and ported time calibration; free-τ analysis; paired scan comparison | `timecal`; `{diag,res,cost,intraw}`; `--scans` | `output/mf_validate/timecal/`, `output/mf_fit_study_20261023/` | 5359d63; untracked; untracked |
| Diagnostics | `regauss_tables.py` | fixed150 vs ±2σ Gaussian per run type from a scan's stored histograms | `--scan-dir --out` | `output/regauss_20261025/` | untracked |
| Diagnostics | `mf_ampshape_diag.py` | Option (ii) offline refit: shape interpolated in ln A between band templates | `--dumps --templates --out` | `output/templates_spline_20261025/ampshape_diag/` | untracked |

The driver hard-codes `--calibrate --exclude-rowcol 15 20 15 20 --inner-radius 200 --min-per-module 30` and LMS `auto` where a gain_corr file exists. Flags are from argparse (inventory 2026-10-24).

## 3. Findings

Gaussian values in rows sourced from the 20261021–20261024 notes and plans, and from `plans/plan_20261025.md` entries before the 10-24 19:12 definition change, use the old fixed ±150 MeV × kS window (fixed150); this affects the preferred fit at 2.24 and 0.73 GeV H2 only. Rows from `plans/plan_20261025.md` after the change use the iterative ±2σ window unless marked fixed150 (§3.7).

### 3.1 Pulse reconstruction

| Finding | Key numbers | Source |
|---|---|---|
| **mf lowers σ/E in every non-X17 run type** at the gate-matched working point (93 splits) | Δ legacy −0.061 ± 0.006 (3.49 H2), −0.040 ± 0.002 (2.24 H2), −0.073 ± 0.006 (C), −0.306 ± 0.020 (0.73 H2) | 20261022 §4.1 |
| **The gain is largest at 0.73 GeV**, where rawsum loses low-amplitude charge | rawsum captures 88–90 % of the charge at 5–20 ADC; Q_template/Q_raw 1.08–1.13 at 5–20 ADC | 20261022 §3, §8; 20261023 §6.1 |
| **The fixed two-τ template is 1–2 % short of the raw sum** at 100–1000 ADC | template integral / raw sum 0.978–0.994; κ = 1.010 corrects the pooled average only | 20261023 §6.1 |
| **Freeing τ does not fix the shape** | bounded +0.026 ± 0.007 (25320), +0.081 ± 0.003 (24503); prior −0.038 ± 0.006 (25320), −0.002 ± 0.004 (24503); χ²/dof stays 100–730 at high A | 20261023 §6 |
| **Sampled per-module templates cut held-out χ²/dof** to a fraction of the two-τ value | ratio 0.15–0.18 / 0.07 / 0.03 at 100–300 / 300–1000 / 1000–2500 ADC (Checkpoint 1); 0.16–0.20 / 0.06–0.07 / 0.03 after rebuild | `plans/plan_20261024.md` progress log 10-23 18:45, 20:40 |
| **Coverage of the final banks**, modules with fallback 0/1/2 | 0.73 577/543/32; 2.24 H2 1062/89/1; C 742/310/100; Ta 812/297/43; 3.49 858/268/26 | `output/templates_spline/<slug>/summary.tsv`; 20261023 §9.1 |
| **Real amplitude dependence of 2–3 %** in the template integral | I_hi/I_lo 0.966–0.977 (high-band integral 2.7–3.4 % below low band); no time walk (Δt ≈ −0.06 ns) | `plans/plan_20261024.md` log 10-23 18:45 |
| **A three-band f(A) per run type removes most of the deficit** | bands 100–400, 400–1000, 1000–2500 ADC (0.73: no top band); knots f 1.0032–1.0068 at 165–201 ADC, 0.9949–0.9988 at 507–592 ADC, 0.9911–0.9973 at 1365–1592 ADC. Residual: max bin-to-bin spread of template/raw sum on template pulses 0.10–0.53 % with f(A), 0.67–1.14 % without (2-SE flatness rule not met); step between 300–1000 and 1000–4000 ADC on the study pulses 0.6–1.1 % (0.64 % on 25320, 1.09 % on 24503) | 20261023 §9.1, §9.2; `output/templates_spline/ampcorr.json`; `output/mf_fit_study_20261023/int_vs_rawsum/int_vs_rawsum_sampled.md` |
| **Sampled templates close the integral deficit on the study pulses** | median Q_template/Q_raw at 100–300 / 300–1000 / 1000–4000 ADC (25320, 24503, 24590): 0.978–0.995 two-τ fixed → 0.990–1.001 sampled | 20261023 §9.2; same file |
| **Sampled templates bring mf χ²/dof to ≈ 1** at A > 500 ADC | 25320: 233 (two-τ fixed) → 1.08 (sampled), σ² + (0.03·model)² | 20261023 §9.2; `output/mf_fit_study_20261023/uniform/diag_summary_sampled.txt` |
| **Spline banks with refitted gains lower σ/E in every run type** (B − A, centred rectangle, legacy) | −0.114 ± 0.006 (3.49 H2), −0.119 ± 0.002 (2.24 H2), −0.074 ± 0.021 (0.73 H2), −0.108 ± 0.005 (C), −0.082 ± 0.029 (Ta); gated within 0.001 pp of legacy except Ta −0.107 ± 0.032; ε_yield 0.994–1.001, ε_peak 1.000–1.013; template success at 2.24 and 3.49 GeV H2 (`plans/plan_20261024.md` §8.3) | 20261023 §9.5; `output/scan_20261024_spline/compare/arms.md` |
| **Candidate selection agrees with Wassim's macro** on 24386.00010 | 0 mismatches over 227 908 (entry, channel) rows, 798/798 best agree (G1-1) | `plans/plan_20261024.md` log 10-23 18:45 |
| **The reference f(A) leaves a 0.7–1.4 % spread on the study pulses** (sampled arm; acceptance bins 100–250, 250–500, 500–1000, 1000–2000, 2000–4000 ADC) | spread 1.06 % (25320), 1.35 % (24503), 0.72 % (24590); dip at 250–1000 ADC and rise above 1000 ADC at 2.24 and 3.49 GeV; worst pulls 44–90 SE; no fallback-class dependence | `plans/plan_20261025.md` log 10-24 15:30; 20261025 §3.1 |
| **Finer bands (option i) do not flatten the study pulses** | `build --bands 100,250,500,1000,2000,4000`; knots f 1.0044–1.0097 at 100–250 ADC, 0.9819–0.9913 at 2000–4000 ADC (no band above 1000 ADC at 0.73 GeV); spread 1.05 % (25320), 1.39 % (24503), 0.37 % (24590); the 250–1000 dip deepens by 0.2 %; per-module RMS of the G1-7 ratio at A > 300 unchanged (0.22–0.86 %), so the option (iii) trigger fires | `plans/plan_20261025.md` log 10-24 15:45; 20261025 §3.2 |
| **Per-module f(A) (option iii) gives the smallest study-pulse spread** | `ampcorr --per-module-f` on the option (i) banks → `output/templates_spline_20261025_pm/`. Review fix (round 2, blocking): non-own cells first took the run-type knot's absolute f, a step of median 0.9 % (max 10 %) above a module's last own knot; the `fa_percrystal` dump pass (14 of 93 jobs) was stopped and deleted, and the banks rebuilt chaining from the nearest own knot (v4). v4: spread 0.19 % (25320), 0.17 % (24503), 0.36 % (24590), not flat by 2 SE; G1-7 spread 0.11–0.37 %; f range 0.878–1.231 (0.878–1.087 without W335) | `plans/plan_20261025.md` log 10-24 16:15, 16:35, 17:00, 17:10; 20261025 §3.3 |
| **Amplitude-indexed shapes (option ii) flatten the 2.24 and 3.49 GeV study waveforms**, offline only | `mf_ampshape_diag.py` on stored study waveforms, (a) module template + f(A) → (b) shape interpolated in ln A between band templates: spread 0.58 → 0.23 % (25320), 0.58 → 0.16 % (24503), 0.19 → 0.23 % (24590); not flat by 2 SE (worst pull 6.6–8.8); χ²/dof lower in every bin (24503, 2000–4000 ADC: 1.26 → 0.68); t_b − t_a < 0.07 ns. The filter implementation needs a separate go | `plans/plan_20261025.md` log 10-24 16:05; 20261025 §3.4 |
| **Per-module f(A) in production makes σ/E worse; rejected** (`scan_20261025_fa`: pm banks, refitted gains `percrystal_fa`; vs reference, legacy, paired) | +0.043 ± 0.004 (3.49, CB), +0.029 ± 0.002 (2.24 H2; +0.031 ± 0.001 fixed150), +0.044 ± 0.003 (C, CB); 0.73 −0.003 ± 0.007; Ta −0.003 ± 0.009; per-module medians 1.61 → 1.63 % (3.49), 1.83 → 1.90 % (C). The current f(A) stays. Reading (untested): per-module band ratios from ≤ 200 pulses imprint template noise, and Q_raw is not an unbiased truth | `plans/plan_20261025.md` log 10-24 18:53, 20:20; `output/scan_20261025_fa/compare/fa_vs_ref_2sigma.md` |
| **W335 at 0.73 GeV has a broken module template, in the reference bank too** | χ²/ndf 3375 (≤ 119 for every other 0.73 GeV module); I_int 4.52 vs 7.0–7.3 in the other run types and the two-τ value; lohi_iint_ratio 0.65; present in `output/templates_spline/0p73GeV_H2/`, so R1 and R2 carry it; one module of 1103; not rebuilt | `plans/plan_20261025.md` log 10-24 17:10, 20:20 |

### 3.2 LMS gain tracking

| Finding | Key numbers | Source |
|---|---|---|
| **LMS correction leaves post-refinement σ/E unchanged** (025320, 8 splits) | legacy Gauss (fixed150) / CB 2.28 / 2.07 % with and without LMS | 20261021 §3 |
| **In-run drift is small** | ≤ 1.2 % relative to batch 0 (factor 0.988–1.010) | 20261021 §3 |
| **Coverage in the scans** | LMS applied only to 25308 and 25320 (only runs with a gain_corr file) | 20261021 §4.1 |
| **Local gain_corr files exist for two runs** | the local splits are a sparse sample (220–3600 LMS events per run); at `-b 1600` only 24503 (2 batches) and 24512 (1 batch) give files (`output/gain_corr_20261024/`); the tool drops the trailing partial batch | `plans/plan_20261025.md` F12 |
| **Two-run on/off check: no change** (R2 + `--lms-corr output/gain_corr_20261024/prad_{run}_gain_corr.root`, 24503 and 24512, 17 splits) | 3.49 GeV H2 CB legacy 1.619 → 1.619 %, Δ +0.000 ± 0.000 [SD 0.002]; gated +0.000 ± 0.001; ε_yield 0.9999 | `plans/plan_20261025.md` log 10-24 20:58; `output/scan_20261025_p3/lms_2run/compare/lms_vs_r2.md` |
| **24503 carries one real gain step; 24512 is an identity check** | 24503: batch 1 starts at event 8 231 291, 25.7 min after batch 0; batch-1 factors median 1.0005, p05–p95 0.9990–1.0021, min 0.9966; splits 00010–00040 in batch 0 (f ≡ 1), 00050 straddles, 00060–00100 in batch 1. 24512 (one batch, f ≡ 1): `result.tsv` and `result_window.tsv` identical to R2 except `elapsed_s` and the `lms_corr` header in all 7 splits | `plans/plan_20261025.md` log 10-24 20:58 |
| **Drift bound for the runs without a file: ≤ 0.06 %** | `-b 200` diagnostic over the sampled splits of all 12 runs: raw per-module LMS means 0.9996–1.0006; the larger jumps are failed alpha fits on the reference PMTs (W1124 in every run, W781 in the Ta runs); LMS stays off for the 10 runs without a file | `plans/plan_20261025.md` F12, log 10-24 20:58; `output/gain_corr_20261024/b200_diagnostic/` |
| **Full-run files deferred to ifarm** | `replay_gainCorr -b 1600` (or 4000) over every split of each run on ifarm; the 12 files go into `database/gain_factor/gain_correction/` there, never from local partial data | `plans/plan_20261025.md` §6(b) |

### 3.3 Timing and pile-up

| Finding | Key numbers | Source |
|---|---|---|
| **No global W beats legacy on e-p or C** | best gain ≤ 0.004 pp (tcal_w12, half A); at W = 16 Δgate (half B) +0.001 on 3.49 H2, 2.24 H2 and C | `plans/plan_20261021.md` §8 |
| **Working point: raw times, W = 16 ns, pre-window 100–200 ns**, unchanged with mf | offsets cut seed–neighbour Δt σ_MAD 3.10 → 1.05 ns but help only at W ≤ 8 | `plans/plan_20261021.md` §8; 20261022 §7 |
| **Tight gates hurt** | mean of 5 run types, half A: traw_w4 +0.45, tcal_w2 +0.14, tcal_w4 +0.03 pp; traw_w8 fails ε_peak (0.987–0.989) | `plans/plan_20261021.md` §8 |
| **The gate gains only on Ta** | rawsum working point Δgate −0.087 ± 0.046 (half B); mf working point −0.099 ± 0.024 | `plans/plan_20261021.md` §8; 20261022 §7 |
| **mf time offsets fail** | σ_MAD with mf offsets 0.780 / 0.810 ns vs rawsum + database offsets 0.672 / 0.709 ns (24503 / 25257) | 20261022 §6 |
| **Pile-up rates explain why only Ta gains from the gate** | candidates / channels with E_dom > 0.5·E_beam: 0.005–0.013 % (3.49), 0.077–0.097 % (0.73), 0.182–0.241 % (2.24 H2 and C), 1.035–1.250 % (Ta) | 20261023 §7, §9 |
| **The +12–40 ns secondary excess is not a template tail** | mean clean single-pulse residual over [12, 40) ns ≤ 0.36 % of peak (≈ 0.4 %) | 20261023 §9.3; `output/pileup_20261024/artefact/resid_window_means.tsv` |
| **Two populations make the excess** (37–50 % of e-p/C candidates, 65–70 % at 0.73 GeV, pre-fix) | (a) burst events with ≥ 10 candidate channels (34–1149 channels, dominant at 186–275 ns, χ²/dof ~10³): 25 % of 2.24 GeV candidates, 9 % at 3.49, none in the two 0.73 splits studied; (b) genuine isolated second bumps at Δt ≈ 23 ns (16–30), E_ratio ≈ 0.18 (0.14–0.28 for A_dom < 1000 ADC), amplitude-proportional, mostly single-channel, in 0.1–0.3 % of pulses (0.29 % at 0.73, 0.17 % at 2.24, 0.09 % at 3.49); origin not identified | `output/pileup_20261024/artefact/burst_events.tsv`, `window_secondaries_isolated.tsv`, `energy_effect.tsv`; 20261023 §9.3 |
| **Study scope** | 6 H2 splits (24503, 24590, 25320 × 2); C and Ta not studied | `output/pileup_20261024/artefact/` |
| **An accepted secondary takes 6–8 % from the dominant pulse** | dominant energy below the single-pulse refit by a median 6.3 % (0.73), 7.2 % (2.24), 7.7 % (3.49), isolated sample; bursts −20 % | `output/pileup_20261024/artefact/energy_effect.tsv` |
| **Fix: burst veto and window rule 40 ns : 0.16** | events with ≥ 10 candidate channels removed from numerator and denominator; a secondary with 0 ≤ dt < 40 ns and E_ratio < 0.16 is not counted (0.17 lost 2–4 efficiency points at ratio 0.2, dt 24–32 ns, r_model 0.10) | `plans/plan_20261024.md` log 10-24; 20261023 §9.3 |
| **Rates before → after the fix** (50 MeV denominator; upper bound on accidental pile-up, since bumps with E_ratio ≥ 0.16 still count) | 3.49 H2 0.231–0.361 → 0.149–0.234 %; 2.24 H2 0.511–0.588 → 0.397–0.407 %; 0.73 H2 0.415–0.598 → 0.260–0.293 %; C 0.515 → 0.364 %; Ta 10.56–11.24 → 10.32–11.04 %; hot channels change the fixed rates by ≤ 0.011 pp | `output/pileup_20261024/pileup_rate_summary.tsv` → `output/pileup_20261024_fix/pileup_rate_summary.tsv` (rate_fixed) |
| **The window rule leaves the efficiency unchanged** | eff_match changes by 0 points in all 30 cells with ratio ≥ 0.2, dt ≥ 24 ns, every run type, r_model 0.3 and 0.10 | `output/pileup_eff_20261024_fix/<slug>/eff_<slug>.md`; 20261023 §9.4 |

### 3.4 Clustering

| Finding | Key numbers | Source |
|---|---|---|
| **Seed-time neighbour amplitudes give no gain** (S_dev, vs control) | +0.001 ± 0.002 (3.49), +0.015 ± 0.002 (2.24 H2), +0.016 ± 0.002 (C), −0.029 ± 0.040 (0.73) | 20261023 §5 |
| **The premise fails: the FADC readout is zero-suppressed** | PbWO4 channels read out per event 7.1 (0.73), 15–16 (2.24 H2, C), 19.9 (3.49); 99.9 % already give mf candidates | 20261023 §5 |
| **`non_linear_corr` off is an exact null** | every calibration file in use has nl1 = nl2 = 0 (the H2/C `_fadc` files carry only the unread legacy key `non_linear` = 0; the X17 files nl1 = nl2 = 0), so `linear_corr` = 1; with `--gain-file` the cluster energy bypasses it. Identity job (025320.00050, 20k events, `--cluster-non-linear-corr 0`): all TSVs identical except the `# cluster_cfg` line | `plans/plan_20261025.md` F6, log 10-24 19:50; 20261025 §6.1 |
| **`corner_conn` on (`cc1`) misses the winner rule** (full scale, R2 + one flag, paired vs R2, ±2σ Gaussian; rule: Δ ≤ −0.02 pp, > 2 SE at 2.24 and 3.49 GeV H2) | Δ legacy ± SE [SD] / gated: 3.49 GeV H2 CB +0.001 ± 0.003 [0.018] / +0.002 ± 0.003; 2.24 GeV H2 Gauss −0.005 ± 0.001 [0.008] / −0.005 ± 0.001; C CB −0.014 ± 0.005 [0.015] / −0.013 ± 0.005; Ta CB +0.031 ± 0.027 [0.133] / +0.033 ± 0.032; ε_yield 0.9985–1.0001, ε_peak 0.9901–1.0002 | `plans/plan_20261025.md` log 10-24 22:47; `output/scan_20261025/cc1/compare/cc1_vs_r2_pooled.md`; 20261025 §6.2 |
| **`cc1` gains at 0.73 GeV and on C** | 0.73 GeV Gauss 3.718 → 3.640 %, Δ −0.077 ± 0.009 [0.035] (gated −0.071 ± 0.009); CB 2.974 → 2.943 %, Δ −0.031 ± 0.017; 2.24 GeV H2 CB 1.858 → 1.849 %, Δ −0.009 ± 0.003; C −0.014 ± 0.005 | `output/scan_20261025/cc1/compare/cc1_vs_r2_pooled.md`, `cc1_vs_r2_splits.tsv` |
| **`min_module_energy` 2 MeV (`mme2`) is an exact null** | 118/118 splits: Δ 0.000 ± 0.000 in every run type and path; TSVs identical to R2 except the `# cluster_cfg` line and, in 11 splits, the `member_identity maxrel` diagnostic (10⁻⁷). Pulse energies never fall below ≈ 2.2 MeV (2.24, 3.49 GeV) and ≈ 3.7 MeV (0.73 GeV), so a threshold ≤ 2 MeV removes no hit | `plans/plan_20261025.md` log 10-24 23:47; `output/scan_20261025/mme2/compare/mme2_vs_r2_pooled.md` |
| **0.5 MeV (`mme0p5`) is an exact null, as the energy floor predicts** | 118/118 splits at full scale: Δ 0.000 ± 0.000 in every run type and path; `result.tsv` and `result_window.tsv` identical to R2 apart from the `# cluster_cfg` line and the floating-point `member_identity` diagnostic. The threshold is below the current 1 MeV and the ≈ 2.2 MeV floor | `plans/plan_20261025.md` log 10-25 01:42; `output/scan_20261025/mme0p5/compare/mme0p5_vs_r2_pooled.md` |
| **5 MeV (`mme5`) is neutral in every H2/C run type** | 118/118 splits, Δ legacy ± SE: 3.49 GeV H2 CB −0.001 ± 0.001; 2.24 GeV H2 Gauss −0.000 ± 0.000; 0.73 GeV Gauss +0.001 ± 0.001, CB +0.001 ± 0.001; C CB −0.000 ± 0.001; ε_yield 1.0000; Ta −0.027 ± 0.014 (report only) | `plans/plan_20261025.md` log 10-25 00:45; `output/scan_20261025/mme5/compare/mme5_vs_r2_pooled.md` |
| **No P4 variant meets the winner rule; no confirmation pass** | `cc1` gains only at 0.73 GeV and on C; `mme2`, `mme0p5` exact nulls; `mme5` neutral | `plans/plan_20261025.md` log 10-25 00:45; 20261025 §6.2 |

Clustering as configured (20261023 §8; `database/reconstruction_config.json` unless noted). The benchmark clusterers copy 8 of these fields and use SimpleProfile; production uses Geant4 profiles.

| Constant | Value |
|---|---|
| min_module_energy | 1.0 MeV (C++ struct default 0) |
| min_center_energy, min_cluster_energy | 10, 50 MeV |
| corner_conn; split_iter, least_split; log_weight_thres | false (4-connected island); 6, 0.01; 3.6 (seed 3×3) |
| seed_time_window | 4.0 ns in JSON; C++ default −1 (off); benchmark legacy −1, gated `--window` |
| non_linear_corr; leakage_correction | true, factor 1 (`_fadc` files carry no nl1/nl2); false |
| WaveAnalyzer peak (`database/daq_config.json`) | max(5σ_ped, 10 ADC); min_peak_ratio 0.3; int_tail_ratio 0.1 |

Three inconsistencies in `docs/technical_notes/hycal_clustering/hycal_clustering.md` (20261023 §8; not edited):

1. Lines 372–373 call W = 16 ns the production default; line 528 says −1; the JSON has 4.0.
2. Line 58 gives min_module_energy 0 MeV (the struct default); line 520 and the JSON give 1.0 MeV.
3. Step 3 (line 111) sets kSplit on shared hits only; `HyCalCluster.cpp:324` sets it on every cluster the splitter emits.

### 3.5 Calibration and acceptance

| Finding | Key numbers | Source |
|---|---|---|
| **The seed-indexed step removes most of the constant term** | pre → post b 2.48 → 0.77–1.3 % | 20261023 §2 |
| **Per-crystal gains help at 2.24 and 3.49 GeV and C** | held-out (S_dev): −0.066 ± 0.003 (3.49), −0.095 ± 0.003 (2.24 H2), −0.091 ± 0.003 (C); full scan, in-sample: −0.048 ± 0.004, −0.091 ± 0.001, −0.090 ± 0.004 | 20261023 §3.2, §4 |
| **At 0.73 GeV the gain depends on the fit; Ta shows no significant change** | 0.73: −0.066 ± 0.022 (Gauss, fixed150), −0.004 ± 0.014 (CB for all); Ta +0.036 ± 0.039 | 20261023 §3.2, §4 |
| **Per-crystal gains are as stable as seed-indexed gains** | run-to-run SD of per-module run means 0.0019 vs 0.0017–0.0018; real component ≈ 0.0016 at 2.24 and 3.49 GeV for both | 20261023 §3.3 |
| **The gain file alone brings pre close to post** | 25320 CB pre 1.856, post 1.838 (mf full scan CB pre 2.893) | 20261023 §4 |
| **The old exclusion rectangle was off-centre** | 0-based reading: x [−41.5, 83.1], y [−83.0, 41.5] mm; centred: x, y in [−62.3, 62.3] mm | 20261023 §8a.1–8a.2 |
| **Per-module σ_M/μ_M medians, centred** | 1.74 (3.49), 2.05 (2.24 H2), 3.70 (0.73), 1.98 (C) % vs pooled 1.83, 2.07, 4.41, 1.97 | 20261023 §8a.4 |
| **Ring 1 is worse at 3.49 GeV; no row/column banding** | ring 1 excess over ring 2 +0.11 (3.49), +0.01 (2.24 H2), +0.03 (C), +2.2 pp (0.73, 12 modules); SD of row/column medians 0.04–0.11 pp < module SD 0.10–0.19 pp | 20261023 §8a.4 |
| **Hot spot W461, W427, W462 stays** | W461 2.92 → 2.49 % (3.49), 2.57 → 2.36 % (2.24 H2); cause not shown | 20261023 §8a.4 |
| **Gains refitted for the spline banks match the h50 gains** | `output/percrystal_spl/`: 14/14 converged in 2 iterations; vs `percrystal_rect` median ratio 0.998–1.004, per-module RMS 0.3–0.5 % | 20261023 §9.5 |
| **Pooled gains + seed-indexed step (arm B) chosen in every run type** | B − C −0.024 ± 0.004 (3.49), −0.023 ± 0.001 (2.24 H2), −0.025 ± 0.002 (C), −0.099 ± 0.027 (Ta); at 0.73 pooled-only is lower by 0.077 ± 0.041, < 2 SE, so B is kept (rule of `plans/plan_20261024.md` §8.3) | 20261023 §9.5; `output/scan_20261024_spline/compare/arms.md` |
| **The cluster energy depends on the impact position in the seed crystal** (R1, legacy pool, 12 × 12 bins over u, v ∈ [−0.6, 0.6] crystal units) | residual R RMS, statistical part removed: 0.44 % (3.49), 0.17 % (2.24 H2), 0.21 % (C), 0.37 % (0.73); peak-to-peak 1.91 % (3.49), 1.15 % (C); fixed150 maps: 2.24 H2 peak-to-peak 1.31 %, 0.73 RMS 0.24 %, peak-to-peak 3.81 %; Ta 11/144 bins, report only. The gate (RMS < 0.3 % in every H2/C type) fails at 3.49 GeV, so the correction was fitted and scanned | `plans/plan_20261025.md` log 10-24 18:58, 20:10; `output/poscorr_20261025/map_summary.tsv`, `before_2sigma/` |
| **One map per run type** | leave-one-out global test: 31 % (3.49) to 91 % (C) of bins within 2 SE (rule ≥ 95 % in every H2/C type); Ta uses the 2.24 GeV H2 map; even/odd pull RMS 0.94–1.12 (bootstrap errors calibrated); bilinear and poly8 agree within 0.003 pp offline at 2.24, 3.49 GeV and C; bilinear scanned | `plans/plan_20261025.md` log 10-24 18:58, 19:02; 20261025 §4.3, §4.5 |
| **Ring and side dependence is significant and small** | ring maps vs the run-type map, n-weighted RMS of the difference: 3.49 GeV 0.15 % (ring 1), 0.16 % (ring 2), 0.29 % (rings 3–4), ring-1 sides 0.21–0.55 % (bottom worst, max 1.66 %); 2.24 H2 0.11–0.27 % (max 0.81 %); C 0.11–0.32 % (max 1.06 %); max \|ΔR\|/SE 6–17 at 2.24 and 3.49 GeV. One map per type corrects the average only | `plans/plan_20261025.md` log 10-24 19:02 |
| **Held-out halves: −0.098 pp at 3.49 GeV** (map from one half applied to the other, stage R, pooled, paired vs R1) | legacy −0.098 ± 0.005 (3.49, CB), −0.009 ± 0.001 (2.24 H2), −0.021 ± 0.005 (C, CB), −0.006 ± 0.008 (0.73); gated within 0.002 pp; harness ε_yield and ε_peak ≥ 0.9997. Success rule met at 3.49 GeV, no H2/C type worse: adopted as R2 | `plans/plan_20261025.md` log 10-24 20:10; `output/scan_20261025/heldout_pooled.md` |
| **Full scan R2 vs R1** (stages R + RX, 118 splits) | legacy −0.099 ± 0.005 (3.49), −0.010 ± 0.001 (2.24 H2), −0.027 ± 0.003 (C), +0.005 ± 0.008 (0.73), −0.021 ± 0.021 (Ta); gated within 0.010 pp; harness ε_yield 0.9993–0.9999, ε_peak 1.0000–1.0003; peak-yield ratio R2/R1 0.986 (0.73), 0.994 (3.49) | `plans/plan_20261025.md` log 10-24 21:39; `output/scan_20261025/pos_all/compare/pos_vs_r1.md` |
| **The map removes the scale residual** | after-map RMS, statistical part removed: 0.44 → 0.05 % (3.49), 0.17 → 0.00 % (2.24 H2), 0.21 → 0.00 % (C), 0.37 → 0.20 % (0.73); peak-to-peak 1.91 → 0.45 % (3.49); gate passes everywhere | `plans/plan_20261025.md` log 10-24 21:39; `output/poscorr_20261025/after/` |
| **The edge width remains** (σ/μ, edge max(\|u\|, \|v\|) > 0.4 vs centre < 0.2) | after: 1.85 vs 1.55 % (3.49), 2.08 vs 1.81 % (2.24 H2), 1.95 vs 1.76 % (C), 4.37 vs 3.41 % (0.73); before (R1): 1.87 vs 1.54 % (3.49), 1.97 vs 1.74 % (C). The edge excess of 0.2–0.3 pp (≈ 1 pp at 0.73 GeV) is in the width | `plans/plan_20261025.md` log 10-24 18:58, 21:39; `output/poscorr_20261025/after/sigma_*.png` |
| **The software Møller veto does not lower 0.73 GeV σ/E** (R1, 40k events per split, 16 splits, CB, legacy, paired vs no veto) | EMIN 30 / 50 / 100 MeV: 24.4 / 16.7 / 6.8 % vetoed, Δ +0.060 ± 0.031 / +0.019 ± 0.025 / +0.018 ± 0.013; with ±20° back-to-back: 2.6 / 1.6 / 0.6 % vetoed, Δ +0.012 / +0.010 / +0.008 | `plans/plan_20261025.md` log 10-24 20:55; `output/scan_20261025_veto/analysis/veto_table.md` |
| **The veto removes low-side events; the peak width stays** | pool fraction with E/E′ < 0.90: 57.0 % → 42.5 / 46.7 / 52.8 % (EMIN 30 / 50 / 100), 55.8 / 56.2 / 56.7 % (±20°); CB peak μ − 1 +0.69 % → +0.61 % (EMIN 30) | same |
| **Re-admitting the Møller-gap modules is worse** | 12.8 modules per split excluded (union 16); re-admission (`--moller-gap-nsigma 0`) doubles the pool, CB Δ +1.09 ± 0.11 (no veto), +1.11 (EMIN 50) | same |
| **3.49 GeV control** (24503, 10 splits) | EMIN 30: 41.6 % vetoed, ΔCB −0.043 ± 0.009; EMIN 50 ±20°: 22.7 % vetoed, ΔCB −0.002 ± 0.002 | same |
| **40k events per split limit the matrix** | ≈ 27 calibrated modules at 0.73 GeV; the ±2σ Gaussian is unstable there (no-veto σ/E 9.3 % vs 3.71 % at full statistics), so Gaussian Δ are not used | same |

| **Full statistics: no veto lowers the 0.73 GeV CB σ/E** (R1 configuration, all 16 full 0.73 GeV splits, legacy, paired vs R1) | `vf30` (EMIN 30, no φ): 24.4 % vetoed, ε_yield 0.756, Gaussian 3.713 → 3.629 % (Δ −0.083 ± 0.009), CB 2.954 → 2.993 % (Δ +0.039 ± 0.022); `vf50p20` (EMIN 50, ±20°): 1.55 % vetoed, ε_yield 0.984, Gaussian Δ +0.001 ± 0.004, CB 2.954 → 2.993 % (Δ +0.038 ± 0.008) | `plans/plan_20261025.md` log 10-24 22:15; `output/scan_20261025_veto/analysis/cmpfull_vf30.md`, `cmpfull_vf50p20.md` |
| **Recommendation: no veto in the 0.73 GeV reference configuration** | the back-to-back veto fires on 0.6–2.6 % of events and keeps 95–99 % of the E/E′ < 0.90 events; the no-φ veto lowers the Gaussian σ/E by 0.08 pp but discards 24 % of events; no setting lowers the CB σ/E | `plans/plan_20261025.md` log 10-24 22:15 |

### 3.6 Resolution decomposition

Source: `output/resolution_fit_20261024/` (`phase0.py` copied there with only its scan dict changed: spl = `scan_20261024_spline`, pc = `scan_20261024_pc`; `out.txt`); old rows `output/resolution_fit_20261023/` (20261023 §2; `plans/plan_20261023.md` §2.1). E = 0.727, 2.236, 3.484 GeV; C and Ta not fitted. Errors SEM only; χ² for 1 dof. The spline and h50 rows use the fixed150 Gaussian at 0.73 and 2.24 GeV. The R1 and R2 rows (2026-10-25) come from `output/resolution_fit_20261025/out.txt` (`phase0.py --also-fixed150`; scans spl 119, R1 119, R2 118 splits; Gaussian refit from the stored histograms, CB stored; `plans/plan_20261025.md` log 10-24 20:14, 21:41).

| Scan / path / estimator | σ/E 0.73 / 2.24 / 3.49 (%) | 2-term a, b (χ²) | 3-term a, b, c |
|---|---|---|---|
| spline, legacy, preferred | 4.366 / 1.952 / 1.709 | 2.62, 0.89 (575) | a² < 0; b 1.62, c 3.16 |
| spline, legacy, CB for all | 2.954 / 1.875 / 1.709 | 2.20, 1.20 (58) | 1.55, 1.44, 1.33 |
| spline, gated, preferred | 4.364 / 1.950 / 1.710 | 2.61, 0.89 (579) | a² < 0; b 1.62, c 3.17 |
| spline, gated, CB for all | 2.953 / 1.873 / 1.710 | 2.20, 1.20 (62) | 1.52, 1.45, 1.35 |
| h50 centred (`scan_20261024_pc`), legacy, preferred | 4.440 / 2.071 / 1.823 | 2.79, 0.93 (466) | a² < 0; b 1.68, c 3.12 |
| h50 centred, legacy, CB for all | 3.037 / 2.000 / 1.823 | 2.28, 1.32 (25) | 1.81, 1.51, 1.14 |
| R1 (= spline), legacy, preferred (±2σ) | 3.713 / 1.908 / 1.709 | 3.07, −0.58 (1150) | 0.40, 1.55, 2.43 |
| R1, legacy, CB for all | 2.954 / 1.875 / 1.709 | 2.20, 1.20 (58) | 1.55, 1.44, 1.33 |
| R2 (R1 + position correction), legacy, preferred (±2σ) | 3.718 / 1.898 / 1.610 | 3.12, −0.78 (616) | 1.77, 1.16, 2.08 |
| R2, legacy, CB for all | 2.974 / 1.858 / 1.610 | 2.37, 0.98 (3.6) | 2.25, 1.05, 0.64 |
| R2, legacy, preferred (fixed150) | 4.338 / 1.942 / 1.610 | 2.91, 0.11 ± 0.29 (414) | 1.17, 1.24, 2.85 |

| Finding | Key numbers |
|---|---|
| **The estimator still sets b** | spline 2-term b 0.89 (preferred) vs 1.20 % (CB for all); χ² 25–579 for 1 dof with SEM-only errors |
| **The preferred 3-term fit is unphysical on both centred scans** | a² < 0 for spline and h50 centred; the old working figure b ≈ 1.3 % no longer rests on a physical fit |
| **The CB-for-all 3-term fit is physical** | spline legacy a 1.55, b 1.44, c 1.33 |
| **Spline vs h50 centred**, CB for all, 2-term | b 1.32 → 1.20, a 2.28 → 2.20 |
| **The constant term is ≈ 71 % of σ² at 3.49 GeV** | spline CB-for-all 3-term, 1.44²/1.709² (from rounded b); preferred 3-term gives 89 % but with a² < 0 |

2026-10-25 refit (R1, R2; gated within 0.01 of legacy):

| Finding | Key numbers |
|---|---|
| **The preferred estimator is unusable with the ±2σ Gaussian** | 2-term b < 0 (R1 −0.58, R2 −0.78 %; χ² 1150, 616); at 0.73 GeV the ±2σ Gaussian fits the core of a tailed peak (54 MeV half-width), so the mixed estimator is inconsistent across energies |
| **CB for all is the consistent estimator** | independent of the Gaussian definition; R1 equals the spline row |
| **The position correction lowers b to 0.98 %** (CB for all) | R1 → R2 2-term b 1.20 → 0.98 %, a 2.20 → 2.37 %·√GeV, χ² 58 → 3.6; 3-term R2 a 2.25, b 1.05, c 0.64 |
| **The constant term is ≈ 43 % of σ² at 3.49 GeV in R2** | CB-for-all 3-term, b² 1.11 / 1.610² (`out.txt`); 71 % in R1 |

Old scan (`scan_20261022_mf`, old rectangle, no per-crystal):

| Finding | Key numbers |
|---|---|
| **The 2-term fit depends on the estimator** | a = 1.85–3.13 %·√GeV, b = 0.77–1.52 % across rawsum/mf, legacy/gated, preferred/CB-for-all; mf legacy preferred a 3.13, b 0.77 (χ² 178 for 1 dof; real uncertainty ≈ ±0.4 on a, ±0.5 on b) |
| **Working figure b ≈ 1.3 %** | 3-term, mf legacy preferred: a 2.19, b 1.29, c 2.37 %·GeV |
| **The constant term is about half of σ² at 3.49 GeV** | 47 % at b = 1.3 (17 % at b = 0.77, 63 % at b = 1.49) |
| **c/E dominates at 0.73 GeV** | 3-term terms at 0.73 GeV: a/√E 2.57, b 1.29, c/E 3.27 % |

### 3.7 Gaussian definition change (2026-10-24)

Decided by the user 2026-10-24 ≈ 19:12, installed 19:50 (`plans/plan_20261025.md` Decisions, log 10-24 19:30, 19:50; 20261025 §2).

| Item | Definition |
|---|---|
| ±2σ (default) | Start from the old fit, refit in [μ − 2σ, μ + 2σ], iterate until μ and σ change by < 1 %, at most 5 refits; final window and refit count in `# gauss_fit` (`/5*` when capped) |
| fixed150 (old) | Gaussian around the histogram argmax in a fixed ±150 MeV × kS window; `--gauss-window fixed150` |
| Unchanged | CB values, gains, windows, cluster TSVs; the elastic-window μ_global fit, the CB seeds and the peak-yield s_ref stay on the fixed fit |
| Old scans | Refit offline from the inner-ring pre/post histograms stored in each `result.tsv` (the input of the in-job fit); the fixed150 refit reproduces the stored σ within 0.005 MeV on all 952 values of `scan_20261024_spline` |

**Every Gaussian number in this summary from here on uses ±2σ unless marked fixed150.** Older entries keep their values; their Gaussian values are fixed150 (§3 preamble, §4 footnotes).

Post σ/E (%), R1, mean over splits; half-width in MeV (fixed150: nominal 150·kS; ±2σ: mean final); cap = splits at the 5-refit cap. Source: `output/regauss_20261025/r1.{md,tsv}` (`regauss_tables.py`).

| Run type (preferred fit), splits | Legacy fixed150 → ±2σ (Δ) | Gated fixed150 → ±2σ (Δ) | Half-width fixed150 / ±2σ | Refits; cap |
|---|---|---|---|---|
| 3.49 GeV H2 (CB), 28 | 1.783 → 1.745 (−0.038) | 1.783 → 1.745 (−0.039) | 187 / 122 | 2.0; 0 |
| **2.24 GeV H2 (Gauss), 42** | **1.952 → 1.908** (−0.044) | 1.950 → 1.907 (−0.044) | 150 / 85 | 2.0; 0 |
| **0.73 GeV H2 (Gauss), 16** | **4.366 → 3.713** (−0.653) | 4.364 → 3.712 (−0.652) | 86 / 54 | 4.0; 5 |
| 2.24 GeV C (CB), 8 | 1.953 → 1.905 (−0.048) | 1.952 → 1.904 (−0.048) | 150 / 85 | 2.0; 0 |
| 2.24 GeV Ta (X17) (CB), 25 | 2.109 → 2.055 (−0.054) | 2.037 → 1.981 (−0.056) | 150 / 92 | 2.0; 0 |

Bold: the preferred fit is the Gaussian. On the CB run types the preferred σ/E does not change (CB 1.709, 1.858, 2.000 legacy). h50 centred (`scan_20261024_pc`): 2.071 → 2.042 % (2.24 GeV H2), 4.440 → 3.771 % (0.73 GeV H2) (log 10-24 19:30).

- The peak position moves +0.5–0.7 MeV at 2.24 and 3.49 GeV and +2.8 MeV at 0.73 GeV (723.1 → 725.9 MeV, legacy post).
- At 0.73 GeV the refit does not converge in 5 of 16 post spectra per path (8 of 16 pre): it alternates between two bin sets on the tailed peak. The 5-refit cap is kept (no mean-of-last-two estimator); capped spectra are flagged.

## 4. Current best numbers

Reference: **R2 = `output/scan_20261025/pos_all`**, R1 plus the position correction (`--pos-corr output/poscorr_20261025/{runtype}/poscorr_{runtype}_all.json`, bilinear, per run type; Ta with the 2.24 GeV H2 map; harness stages R + RX, 118 splits, no 24386.00060). R1 = `output/scan_20261025_ref`, the driver rerun of `output/scan_20261024_spline` (spline banks, κ = 1.0, centred rectangle, pooled gains `output/percrystal_spl/<run>/gains_<run>_all.tsv`) with the centroid columns; identical to it except `elapsed_s` and the two appended cluster-TSV columns. h50 baseline: `output/scan_20261024_pc` with `output/percrystal_rect/<run>/gains_<run>_all.tsv`. Rows: `output/scan_20261021_full` (rawsum), `output/scan_20261022_mf` (mf), `output/scan_20261023_pc` (mf + per-crystal, old rectangle), `output/scan_20261024_pc` (centred), `output/scan_20261024_spline` = R1 (spline), `output/scan_20261025/pos_all` (R2). σ/E in %, preferred fit, mean over splits. Δ is the step from the row above, legacy / gated, in pp: paired ± SE for the per-crystal, spline and R2 rows (R2 on the splits it shares with R1: 27 at 3.49 GeV), a difference of scan means for the mf row and the definition-change row.

**Gaussian definition (§3.7).** In the 2.24 and 0.73 GeV H2 tables (preferred fit Gaussian), the rows above "R1, ±2σ Gaussian" are fixed150; that row is the same scan refit with the ±2σ window, and the R2 row uses ±2σ. CB tables are unaffected.

**3.49 GeV H2 (CB, 28 splits; rawsum and R2 27)**

| Configuration | Legacy | Gated | Δ legacy / gated |
|---|---|---|---|
| rawsum | 1.936 | 1.971 | |
| mf | 1.873 | 1.873 | −0.063 / −0.098 |
| mf + per-crystal, old rectangle | 1.826 | 1.825 | −0.048 ± 0.004 / −0.048 ± 0.004 |
| mf + per-crystal, centred | 1.823 | 1.823 | −0.003 ± 0.010 / −0.002 ± 0.010 |
| spline + refitted per-crystal, centred | 1.709 | 1.710 | −0.114 ± 0.006 / −0.113 ± 0.006 |
| spline + position correction (R2) | **1.610** | 1.609 | −0.099 ± 0.005 / −0.101 ± 0.005 |

**2.24 GeV H2 (Gaussian, 42 splits)**

| Configuration | Legacy | Gated | Δ legacy / gated |
|---|---|---|---|
| rawsum | 2.263 | 2.293 | |
| mf | 2.223 | 2.220 | −0.040 / −0.073 |
| mf + per-crystal, old rectangle | 2.132 | 2.129 | −0.091 ± 0.001 / −0.090 ± 0.001 |
| mf + per-crystal, centred | 2.071 | 2.069 | −0.061 ± 0.002 / −0.060 ± 0.002 |
| spline + refitted per-crystal, centred | 1.952 | 1.950 | −0.119 ± 0.002 / −0.119 ± 0.002 |
| R1, ±2σ Gaussian | 1.908 | 1.907 | −0.044 / −0.044 (definition change) |
| spline + position correction (R2), ±2σ Gaussian | **1.898** | 1.897 | −0.010 ± 0.001 / −0.010 ± 0.001 |

Footnote: Gaussian fixed150 in the rows above "R1, ±2σ Gaussian" (§3.7).

**0.73 GeV H2 (Gaussian, 16 splits)**

| Configuration | Legacy | Gated | Δ legacy / gated |
|---|---|---|---|
| rawsum | 4.655 | 4.680 | |
| mf | 4.350 | 4.352 | −0.304 / −0.328 |
| mf + per-crystal, old rectangle | 4.285 | 4.285 | −0.065 ± 0.020 / −0.066 ± 0.020 |
| mf + per-crystal, centred | 4.440 | 4.440 | +0.154 ± 0.019 / +0.154 ± 0.019 |
| spline + refitted per-crystal, centred | 4.366 | 4.364 | −0.074 ± 0.021 / −0.075 ± 0.021 |
| R1, ±2σ Gaussian | 3.713 | 3.712 | −0.653 / −0.652 (definition change) |
| spline + position correction (R2), ±2σ Gaussian | **3.718** | 3.714 | +0.005 ± 0.008 / +0.002 ± 0.007 |

Footnote: Gaussian fixed150 in the rows above "R1, ±2σ Gaussian" (§3.7).

**2.24 GeV C (CB, 8 splits)**

| Configuration | Legacy | Gated | Δ legacy / gated |
|---|---|---|---|
| rawsum | 2.166 | 2.19 | |
| mf | 2.093 | 2.092 | −0.073 / −0.095 |
| mf + per-crystal, old rectangle | 2.003 | 2.002 | −0.090 ± 0.004 / −0.090 ± 0.004 |
| mf + per-crystal, centred | 1.966 | 1.965 | −0.037 ± 0.004 / −0.038 ± 0.004 |
| spline + refitted per-crystal, centred | 1.858 | 1.856 | −0.108 ± 0.005 / −0.109 ± 0.007 |
| spline + position correction (R2) | **1.831** | 1.830 | −0.027 ± 0.003 / −0.026 ± 0.003 |

**2.24 GeV Ta, X17 (CB, 25 splits; reported, not interpreted)**

| Configuration | Legacy | Gated | Δ legacy / gated |
|---|---|---|---|
| rawsum | 2.428 | not recorded | |
| mf | 2.359 | not recorded | per run: −0.113, −0.203 (paired, Round 3; unpaired −0.198), +0.000 / −0.148, −0.344, +0.054 |
| mf + per-crystal, old rectangle | 2.395 | 2.274 | +0.036 ± 0.039 / +0.014 ± 0.031 |
| mf + per-crystal, centred | 2.082 | 2.016 | −0.313 ± 0.055 / −0.258 ± 0.054 |
| spline + refitted per-crystal, centred | 2.000 | 1.909 | −0.082 ± 0.029 / −0.107 ± 0.032 |
| spline + position correction (R2) | **1.978** | 1.878 | −0.021 ± 0.021 / −0.031 ± 0.024 |

**Per-module median σ/E post, h50 centred → spline → R1 (±2σ) → R2 (`--group-runs`, N > 1000)**

| Run type | Median σ/E (%) | Modules |
|---|---|---|
| 3.49 GeV H2 | 1.74 → 1.61 → 1.61 → 1.58 | 112 → 111 → 111 → 111 |
| 2.24 GeV H2 | 2.05 → 1.92 → 1.87 → 1.88 | 125 → 124 → 124 → 123 |
| 0.73 GeV H2 | 3.70 → 3.60 → 3.52 → 3.50 | 97 → 99 → 99 → 98 |
| 2.24 GeV C | 1.98 → 1.83 → 1.83 → 1.83 | 78 → 78 → 78 → 78 |
| 2.24 GeV Ta (X17) | 1.93 → 1.99 → 1.99 → 1.91 | 8 → 8 → 8 → 8 |

The 2.24 and 0.73 GeV H2 medians come from per-module Gaussian fits: fixed150 for h50 centred and spline, ±2σ for R1 and R2.

Sources:

- Legacy rawsum, mf and old-rectangle values and their Δ: 20261023 §4.
- Gated rawsum and mf values, rawsum → mf Δ at 3.49, 2.24 and 0.73 GeV: `plans/plan_20261023.md` §2.1.
- C and Ta rawsum → mf Δ (25308; 25800, 26138, 26230): 20261022 §4.2. C rawsum gated: 20261021 §4.1 (two decimals); C mf gated: `plans/plan_20261023.md` Results.
- Old → centred: 20261023 §8a.3 and `plans/plan_20261023.md` "Centred exclusion rectangle" (`output/scan_20261024_pc/compare/rect_old_vs_centred.tsv`).
- Centred → spline: `output/scan_20261024_spline/compare/arms.md` (20261023 §9.5).
- Per-module medians: `output/module_resolution_20261024/summary/types_summary.tsv` → `output/module_resolution_20261024_spline/summary/types_summary.tsv` → `output/module_resolution_20261025_ref/summary/types_summary.tsv` (R1) → `output/module_resolution_20261025_pos/summary/types_summary.tsv` (R2).
- R1 identity with `scan_20261024_spline`: `plans/plan_20261025.md` log 10-24 17:16. Definition change (R1 fixed150 → ±2σ): `output/regauss_20261025/r1.{md,tsv}` (log 10-24 19:30).
- R1 → R2: `output/scan_20261025/pos_all/compare/pos_vs_r1.{md,tsv}` (log 10-24 21:39).

Caveats:

- **In-sample gains.** The `_all` tables (`percrystal_rect` and `percrystal_spl`) are fitted on the splits they are applied to. A held-out test exists only for S_dev (§3.5). R2 applies the all-split position map in-sample; the held-out halves give the same size at 3.49 GeV (−0.098 ± 0.005 vs −0.099 ± 0.005, §3.5).
- **Gains not refit after the position correction.** The pooled gains are fitted on member energies, which the cluster-level factor does not touch; the map has mean 1 (`plans/plan_20261025.md` §5.3).
- **Peak yield with the position correction.** Peak-yield ratio R2/R1 0.986 at 0.73 GeV and 0.994 at 3.49 GeV; harness ε_peak 1.0000–1.0003; not investigated.
- **The centred row measures a sample change.** The centred rectangle removes high-rate top/left modules and halves the calibrated-block events at 2.24 and 3.49 GeV; at 0.73 GeV it adds bottom/right modules (W667, W700, W701), and σ/E rises by 0.15 pp (Møller contamination is an untested hypothesis). The calibration method is the same in both rows (20261023 §8a.3).
- **Gate matching.** The rawsum scan ran at W = 8 with the per-run pre-window (130–200 ns for 24xxx); the mf scans ran at W = 16 and 100–200 ns. The gated rawsum → mf step and the 24xxx legacy step mix gate and pulse-method effects. Gate-matched mf gains are in §3.1 (20261022 §4.1).
- **Ta statistics.** ≈ 16 calibrated modules per split.
- **Uncommitted.** The benchmark changes that produce R1 and R2 (centroid columns, `--pos-corr`, `--gauss-window`) are uncommitted, and `position_corr.py` is untracked (§5.2).

## 5. Status and open items (`plans/plan_20261025.md` log to 10-25 01:42)

### 5.1 Work items

| Item | State | Where |
|---|---|---|
| Item 3 of `plan_20261025`: LMS | Done: two-run check Δ +0.000 ± 0.000 pp; drift bound ≤ 0.06 %; LMS off for the 10 runs without a file; full-run files deferred to ifarm (§3.2) | `plans/plan_20261025.md` log 10-24 20:58; plan §6 |
| Item 4: clustering checks | `non_linear_corr` off: exact null (identity job). P4 done at full scale: cc1, mme2, mme0p5, mme5. `cc1` misses the winner rule (gains at 0.73 GeV and on C); `mme2` and `mme0p5` exact nulls; `mme5` neutral. No winner, no confirmation pass; nothing running | `plans/plan_20261025.md` log 10-24 19:50, 21:39, 22:47; 10-25 00:45, 01:42 |
| Add-on: Gaussian definition | Done: ±2σ default in the benchmark and tools since 10-24 19:50; old scans refit from stored histograms; at 0.73 GeV 5 of 16 post spectra per path end at the 5-refit cap (§3.7) | `plans/plan_20261025.md` Decisions; log 10-24 19:30, 19:50 |
| Add-on: software Møller veto | Done: 40k-event matrix and full-statistics confirmation (`vf30`, `vf50p20`); no setting lowers 0.73 GeV CB σ/E; no veto in the 0.73 GeV reference configuration (§3.5) | `plans/plan_20261025.md` log 10-24 20:55, 22:15 |
| Decomposition refit | Done with R1 and R2; CB for all 2-term b 1.20 → 0.98 % (§3.6) | `plans/plan_20261025.md` log 10-24 20:14, 21:41 |
| Note `notes/PbWO4_energy_resolution_20261025.md` | Complete, including P4 (§6.2); the veto confirmation is in §7 | `plans/plan_20261025.md` log 10-24 20:20; 10-25 01:42 |
| Phase 3 of `plan_20261024` (spline banks in the calibration) | Done: driver scan `output/scan_20261024_spline/` (119/119 splits, no failures); arm B in every run type; B − A −0.074 to −0.119 pp legacy; success at 2.24 and 3.49 GeV H2 (§3.1, §3.5, §4) | `plans/plan_20261024.md` §8, progress log 10-24 11:50; 20261023 §9.5 |
| Pile-up fix rescans | Done: rates `output/pileup_20261024_fix/` (119 splits); efficiency `output/pileup_eff_20261024_fix/` unchanged (0 points in all 30 cells) (§3.3) | 20261023 §9.3, §9.4 |
| f(A) residual (item 2 of `plan_20261025`) | Done, rejected: per-module f(A) flattens the study pulses (spread 0.17–0.19 %) and worsens σ/E (+0.029 to +0.044 pp at 2.24 and 3.49 GeV); the current f(A) stays; option (ii) feasible offline, its filter implementation needs a separate go (§3.1) | 20261023 §9.2, §9.5; `plans/plan_20261025.md` log 10-24 15:30–18:53, 20:20 |
| Origin of the ≈ 23 ns second bump | Not identified; ask DAQ whether it is a reflection (§3.3) | 20261023 §9.3 |
| X17 statistics | Waiting for new X17 evio files (disk space; "next week") | 20261023 §10; `plans/plan_20261023.md` Decisions |
| C++ deconvolver, `PulseTemplate`, `WaveAnalyzer.deconvolve` | Still two-τ; later task | `plans/plan_20261024.md` §12 |
| Energy-dependent W(E) (C++ change) | Approved in principle; not started | `plans/plan_20261022.md` §1, §9 |
| Position-dependent correction (Phase 2 of `plan_20261023`, item 1 of `plan_20261025`) | Done, adopted: R2 = `output/scan_20261025/pos_all` (held-out −0.098 ± 0.005, full −0.099 ± 0.005 pp at 3.49 GeV); ring dependence and edge width not corrected (§3.5, §4). The Phase 4 quick checks became item 4 | `plans/plan_20261025.md` log 10-24 18:58–21:39 |
| mf time offsets | Failed acceptance (§3.3); open physics item | 20261022 §6 |
| Hot channels | Pile-up: W335 (≈ 90 % at 0.73 GeV), W1 (37–56 %), W1019/W1020 (20–64 %), W934 (+25 ns echo). Resolution: W461, W427, W462. Template: W335 at 0.73 GeV broken in the reference bank, not rebuilt (§3.1) | `plans/plan_20261024.md` log 10-23 21:50; 20261023 §8a.4; `plans/plan_20261025.md` log 10-24 17:10, 20:20 |

### 5.2 Git state

| State | Files |
|---|---|
| Committed in `315e4af` (10-21) and `5359d63` (10-23) | benchmark, driver, `gate_scan_2026102{1,2}.py`, `mf_pulse.py`, `mf_validate.py`, `percrystal_calib.py`, 20261021 and 20261022 notes, plans 10-21 to 10-23, `docs/research-lms-gain-tracking.md` |
| Committed in `80e19b6` (10-24) | benchmark, driver, `mf_pulse.py`, `percrystal_calib.py`, `pileup_generator.py`, `spline_template_extract.py`, `mf_fit_pulse.py`, `mf_fit_study.py`, `plot_module_resolution.py`, `plot_pileup_examples.py`, `benchmark_pileup_mf.py`, `pileup_artefact_study.py`, `compare_scan_dirs.py`, 20261021–20261023 notes, this note, `plans/plan_20261024.md`, the move of the plans into `plans/`, `docs/research-spline-template-pipeline.md`, `dev.md` |
| Reference, unchanged since `5359d63` | `benchmark_hycal_timing.py` |
| Modified since `80e19b6`, uncommitted | `benchmark_hycal_timing_test_20261021.py`, `spline_template_extract.py`, `mf_pulse.py`, `mf_fit_study.py`, `plot_module_resolution.py`, `compare_scan_dirs.py`, this note (this update) |
| Untracked | `position_corr.py`, `mf_ampshape_diag.py`, `regauss_tables.py`, `plans/plan_20261025.md`, `notes/PbWO4_energy_resolution_20261025.md` |

`output/` is gitignored. Checked with `git status --short` and `git show --stat 80e19b6` (read-only) on 10-24 21:42.

### 5.3 Source disagreements

| Topic | Disagreement |
|---|---|
| Exclusion rectangle | Scans up to `scan_20261023_pc` and `output/percrystal/` read `--exclude-rowcol 15 20 15 20` 0-based (off-centre, `16 21 16 21` in the new reading); the current benchmark reads it 1-based. The 10-21/10-22 commands no longer reproduce their tables. The driver hard-codes `15 20 15 20` |
| Split counts | rawsum full scan 118, mf scans 119 (24386: 10 vs 11); 3.49 GeV 27 in harness tables vs 28 in `scan_20261023_pc`; `scan_20261021_full` holds 129 `result.tsv` (119 plus 10 of the excluded 24578; the summary TSV keeps 118, without 24386.00060) |
| 25800 rawsum σ/E | 2.34 % (20261021 §4.1) vs 2.33 % (20261022 §4.2), same scan |
| rawsum low-amplitude loss at 5–20 ADC | 10–12 % (rawsum/mf 0.88–0.90, 20261022 §3) vs 8–13 % (template/raw-sum ratio, 20261023 §6.1) |
| X17 templates | X17 under two template sets, checked 2026-10-24: not identical. The gate-rescan jobs (output/gate_scan_20261022/mf_traw_w16_p100-200/) used the 26138-only templates (column 26138, 170 modules); the mf full scan (output/scan_20261022_mf/) used the combined templates (column 27000, 377 modules). Post CB σ_E/E: 25800 2.21 vs 2.22, 26230 2.56 vs 2.48; 26138 agrees at 2.23 only by rounding (per-split values differ; means 2.2300 vs 2.2325, n = 4). The −0.203 Δ quoted for 26138 in resolution_mf_vs_rawsum_fullscan.md is the Round-3 paired value; the full scan's own unpaired definition gives −0.198. A direct one-split test (25800.00100, 20k events) confirms `--mf-templates auto` resolves to the combined file and that the two files give different pulse counts (56 416 vs 56 570). |
| Legacy pre-window | `plan_benchmark_hycal_timing_test_20261021.md` caption 130–200 ns; note and driver 100–200 ns for 25xxx/26xxx |
| Ring numbering | 20261023 §8a.2 counts rings around the beam hole (ring 3 at r ≈ 73 mm); §8a.4 counts from the rectangle (ring 1 = r 73–103 mm, used here) |
| 0.73 GeV fit model | Preferred Gaussian 4.37–4.55 % per run (fixed150; ±2σ 3.771 % mean, §3.7) vs CB 3.05 / 3.00 % (24590 / 24603, `scan_20261024_pc`); only 20261021 §4.3 comments (radiative tail) |
| Pile-up rate denominator | 20261023 §7: E_dom > 0.5·E_beam (pre-fix); `plan_20261024`: 50 MeV, pre-fix 0.231–0.598 % e-p/C, 10.56–11.24 % Ta, after the fix 0.149–0.407 %, 10.32–11.04 %. Not comparable |
| Window-rule ratio | `plan_20261024` B-result entry 0.17; its later entry, 20261023 §9.3 and the fixed rescan 0.16 (final) |
| Isolated-bump fraction | `plan_20261024` 0.2–0.3 %; `energy_effect.tsv` 0.09–0.29 %; 20261023 §9.3 0.1–0.3 % |
| Bursts at 0.73 GeV | None in the 6-split study; present in the full rescan (24590 7.1 %, 24603 29.0 % of candidates) |
| f(A) residual, 300–1000 vs 1000–4000 ADC | 20261023 §9.5 0.7–1.1 %; `int_vs_rawsum_sampled.md` 0.64 % (25320), 1.09 % (24503) |
| Git state in the `plan_20261025` header | Header: driver and `percrystal_calib.py` uncommitted; `spline_template_extract.py`, `compare_scan_dirs.py`, `plot_module_resolution.py`, `mf_fit_study.py` untracked. `git` on 10-24 21:42: all committed in `80e19b6`; the last four modified since (§5.2) |
| LMS drift bound | `plan_20261025` Decisions ±0.05 %; F12 and log 10-24 20:58 ≤ 0.06 % (per-module means 0.9996–1.0006) |
| Template-tail residual per ns | 20261023 §9.3 +2.5 / +4 % (A < 150 ADC) and ≤ 0.9 % at 400–1000 ADC; `output/pileup_20261024/artefact/summary.json` +2.24 / +3.58 % and 0.15–0.59 % |

## 6. Reproduction

Environment and the current working point (20261023 §11):

```bash
cd /Users/L00338853/work/PRad/prad2evviewer/analysis/pyscripts && source setup_env.sh
E=/Users/L00338853/work/PRad/data/evio; DB=$PRAD2_DATABASE_DIR/calibration
# --exclude-rowcol is 1-based inclusive: 15 20 15 20 = centred 6 x 6 rectangle
WP="--calibrate --exclude-rowcol 15 20 15 20 --pre-window 100 200 --inner-radius 200 --min-per-module 30 --window 16 --pulse-method mf"
```

(a) One split, mf with pooled per-crystal gains (reference settings):

```bash
python3 benchmark_hycal_timing_test_20261021.py $E/prad_025320.evio.00400 output/<tag>/025320_00400 $WP \
  --lms-corr auto --hc-calib-file $DB/calibration_factor_2p2_fadc.json \
  --gain-file output/percrystal_rect/025320/gains_025320_all.tsv
```

(b) Per-crystal gains: dump one split's calibration pool, then the pooled fit per run (the reference used the Round-2 harness dumps):

```bash
python3 benchmark_hycal_timing_test_20261021.py $E/prad_025320.evio.00050 output/<tag>/025320_00050 $WP \
  --lms-corr auto --hc-calib-file $DB/calibration_factor_2p2_fadc.json --cal-method percrystal --cluster-dump auto
python3 percrystal_calib.py fit --dumps "output/scan_20261023/mf_percrystal/scan_025320/split_*/result_dump.npz" \
  --splits all --fixed-windows --out output/percrystal_rect/025320/
```

(c) Full scan as run for `output/scan_20261024_pc` (20261023 §8a.3, §11(e)):

```bash
python3 benchmark_hycal_timing_scan.py --pulse-method mf --window 16 --pre-window 100 200 --gain-file-map auto \
  --gain-dir output/percrystal_rect --write-cluster-tsv --parallel 6 --timeout-seconds 1800 --output-dir output/scan_<tag>
python3 compare_scan_dirs.py --scans old=output/scan_20261023_pc centred=output/scan_<tag> --out output/scan_<tag>/compare/rect_old_vs_centred
```

(d) Module maps:

```bash
MPLBACKEND=Agg python3 plot_module_resolution.py --scan-dir output/scan_<tag> --out output/module_resolution_<tag> --group-runs
#   --min-n 1000 default; --path gated; --exclude-rowcol 16 21 16 21 for scans before 20261024
```

(e) Spline templates and the reference scan (`output/scan_20261024_spline`; 20261023 §9.1, §9.5):

```bash
T=output/templates_spline
# templates: select per split (Ta: --trigger physics-any), build per run type, ampcorr over the five banks (rewrites the npz), check
python3 spline_template_extract.py select $E/prad_025320.evio.00050 --out $T/2p24GeV_H2/pulses/025320_00050.npz
python3 spline_template_extract.py select $E/prad_025800.evio.00100 --trigger physics-any --out $T/2p24GeV_Ta_X17/pulses/025800_00100.npz
python3 spline_template_extract.py build --pulses "$T/2p24GeV_H2/pulses/*.npz" --run-type 2p24GeV_H2 --out $T/2p24GeV_H2/templates_2p24GeV_H2
python3 spline_template_extract.py ampcorr --out $T/ampcorr.json --templates $T/{0p73GeV_H2,2p24GeV_H2,2p24GeV_C,2p24GeV_Ta_X17,3p49GeV_H2}/templates_*.npz
python3 spline_template_extract.py check --pulses "$T/2p24GeV_H2/pulses/*.npz" --templates $T/2p24GeV_H2/templates_2p24GeV_H2.npz --out $T/2p24GeV_H2/check/
# gains: harness dump pass, stages R and RX (configs entry bench_args "--mf-templates $T/{runtype}/templates_{runtype}.npz --mf-scale 1.0 --cal-method percrystal --cluster-dump auto"), then the pooled fit per run
python3 gate_scan_20261022.py run --stage R --configs-file cfgs_spl.json --tags spl_percrystal --out output/scan_20261024 --parallel 12 --timeout-seconds 1800
python3 percrystal_calib.py fit --dumps "output/scan_20261024/spl_percrystal/scan_025320/split_*/result_dump.npz" --splits all --fixed-windows --out output/percrystal_spl/025320/
#   024386 also takes output/percrystal_spl/dumps_extra/024386_00060_dump.npz
# driver, compare and plot argv below are reconstructed from --help (the logs do not record argv)
python3 benchmark_hycal_timing_scan.py --pulse-method mf --window 16 --pre-window 100 200 \
  --mf-templates-map "$T/{runtype}/templates_{runtype}.npz" --mf-scale 1.0 \
  --gain-file-map auto --gain-dir output/percrystal_spl --write-cluster-tsv \
  --parallel 11 --timeout-seconds 1800 --output-dir output/scan_20261024_spline
python3 compare_scan_dirs.py --scans A=output/scan_20261024_pc B=output/scan_20261024_spline --out output/scan_20261024_spline/compare/arms
MPLBACKEND=Agg python3 plot_module_resolution.py --scan-dir output/scan_20261024_spline --out output/module_resolution_20261024_spline --group-runs
# pile-up with the fix (rates, maps, examples; then efficiency per run type)
MPLBACKEND=Agg python3 plot_pileup_examples.py --jobs 12 --splits all --mf-templates "$T/{runtype}/templates_{runtype}.npz" --dominant-min-mev 50 --mf-scale 1.0 --burst-veto 10 --window-rule 40:0.16 --out output/pileup_20261024_fix
python3 benchmark_pileup_mf.py output/pileup_eff_20261024_fix/2p24GeV_H2/synth --templates $T/2p24GeV_H2/templates_2p24GeV_H2.npz --out output/pileup_eff_20261024_fix/2p24GeV_H2 --tolerance-ns 4 --window-rule 40:0.16
```

(f) Package of `plans/plan_20261025.md`: R1, position maps and fit, held-out and full harness runs (R2), veto matrix, LMS two-run check, P4 variants, Gaussian tables, decomposition (20261025 §10). Configs (`cfgs_20261025_p2.json`, `cfgs_veto.json`, `cfgs_p34_withpos.json`) and `run_*.sh` are in `$S`. The tools do not log argv; (r) marks commands reconstructed from the plan and `--help`.

```bash
S=output/scan_20261025/configs
T0=output/templates_spline; R1=output/scan_20261025_ref; PC=output/poscorr_20261025; OUT=output/scan_20261025
# one split with the R2 settings (the harness prepends --mf-r-model 0.3 for mf) (r)
python3 benchmark_hycal_timing_test_20261021.py $E/prad_025320.evio.00400 output/<tag>/025320_00400 $WP --lms-corr auto \
  --hc-calib-file $DB/calibration_factor_2p2_fadc.json --mf-r-model 0.3 --mf-templates $T0/2p24GeV_H2/templates_2p24GeV_H2.npz --mf-scale 1.0 \
  --gain-file output/percrystal_spl/025320/gains_025320_all.tsv --pos-corr $PC/2p24GeV_H2/poscorr_2p24GeV_H2_all.json
# R1: driver rerun of the reference configuration with the centroid columns (r)
python3 benchmark_hycal_timing_scan.py --pulse-method mf --window 16 --pre-window 100 200 --mf-templates-map "$T0/{runtype}/templates_{runtype}.npz" \
  --mf-scale 1.0 --gain-file-map auto --gain-dir output/percrystal_spl --write-cluster-tsv --parallel 8 --timeout-seconds 1800 --output-dir $R1
# position maps (the fitted maps are fixed150; the ±2σ view is before_2sigma/), fits per run type, global test, Ta symlinks
for SET in all even odd; do MPLBACKEND=Agg python3 position_corr.py map --scan-dir $R1 --splits $SET --gauss-window fixed150 --out $PC/; done
MPLBACKEND=Agg python3 position_corr.py map --scan-dir $R1 --splits all --out $PC/before_2sigma/
for SL in 0p73GeV_H2 2p24GeV_H2 2p24GeV_C 3p49GeV_H2; do for SET in all even odd; do
  python3 position_corr.py fit --map $PC/map_${SL}_$SET.tsv --form both --out $PC/$SL/poscorr_${SL}_$SET.json; done; done
for SET in all even odd; do python3 position_corr.py fit --global --maps "$PC/map_*_$SET.tsv" --out $PC/global_$SET/; done
mkdir -p $PC/2p24GeV_Ta_X17 && for SET in all even odd; do ln -sf ../2p24GeV_H2/poscorr_2p24GeV_H2_$SET.json $PC/2p24GeV_Ta_X17/poscorr_2p24GeV_Ta_X17_$SET.json; done
# harness: pos_even on odd splits, pos_odd on even splits (stage R), pos_all (R + RX, + --write-cluster-tsv); bench_args per tag (r):
#   "--mf-templates $T0/{runtype}/templates_{runtype}.npz --mf-scale 1.0 --gain-file output/percrystal_spl/{run}/gains_{run}_all.tsv --pos-corr $PC/{runtype}/poscorr_{runtype}_<set>.json"
zsh $S/run_p2.sh
python3 compare_scan_dirs.py --scans ref=$R1 pos=$OUT/pos_all --gauss-window 2sigma --out $OUT/pos_all/compare/pos_vs_r1   # stored: tables made before 10-24 19:50
MPLBACKEND=Agg python3 plot_module_resolution.py --scan-dir $OUT/pos_all --out output/module_resolution_20261025_pos --group-runs
MPLBACKEND=Agg python3 position_corr.py map --scan-dir $OUT/pos_all --out $PC/after/
# veto matrix (40k events; tags vb, v30, v50, v100, v30p20, v50p20, v100p20, vb_rm, v50_rm, v50p20_rm), v30 rerun + P3 lms_2run, full-statistics veto
#   e.g. bench_args "<R1 args> --max-events 40000 --write-cluster-tsv --veto-second-cluster 50 --veto-dphi 20 --moller-gap-nsigma 0" (r)
zsh $S/run_veto.sh; zsh $S/run_v30_p3.sh; zsh $S/run_vetofull.sh   # -> output/scan_20261025_veto/, output/scan_20261025_p3/lms_2run/
#   veto analysis: output/scan_20261025_veto/analysis/veto_table.py; P3 compare vs $OUT/pos_all -> lms_2run/compare/lms_vs_r2
# P4 variants, all four run at full scale, compare vs R2 after each
zsh $S/run_p4.sh cc1 mme2 mme0p5 mme5                              # -> $OUT/<tag>/compare/<tag>_vs_r2
# Gaussian tables (fixed150 vs ±2σ from the stored histograms) and the decomposition refit
for X in spline:scan_20261024_spline pc:scan_20261024_pc r1:scan_20261025_ref; do
  python3 regauss_tables.py --scan-dir output/${X#*:} --out output/regauss_20261025/${X%%:*}; done
(cd output/resolution_fit_20261025 && python3 phase0.py --also-fixed150 > out.txt)   # (r) cwd and redirect
```

**Plan §13 caveat.** The commands in `plans/plan_20261024.md` §13 are stale (old `_spl` path, `A=scan_20261023_pc`, no `ampcorr` step, no Ta `--trigger physics-any`, no fix flags); use (e). The `pileup_generator.py` synth step and `mf_fit_study.py --arm sampled` stay in §13 (κ = 1.0, `_fix` output dirs).

**Rectangle caveat.** For scans up to `scan_20261023_pc`, pass `--exclude-rowcol 16 21 16 21` to the benchmark, `percrystal_calib.py fit` and `plot_module_resolution.py`. The driver has no rectangle flag, so older driver scans do not reproduce with the current benchmark.

Other commands: harness variant scans in `plans/plan_20261023.md` §13 (configs file in a scratch directory); gate scans in `plans/plan_20261021.md` §3.7 and 20261022 §10; spline templates, spline gains, reference scan and pile-up fix in (e); the `plan_20261025` package in (f) and 20261025 §10.

### See also

- Plans: `plans/plan_benchmark_hycal_timing_test_20261021.md` (benchmark design), `plans/plan_20261021.md` §8 (Phase 1 gate scan), `plans/plan_20261022.md` (mf, Rounds 3–4), `plans/plan_20261023.md` (Phases 0–4, Round 2, free τ, centred rectangle), `plans/plan_20261024.md` (spline templates, pile-up, Phase 3; progress log).
- `docs/technical_notes/waveform_analysis/fit_pulse_template.md`: two-τ template extraction and Table 4.
- `docs/technical_notes/hycal_clustering/hycal_clustering.md`: clustering and the earlier calibration-equalised resolution.
