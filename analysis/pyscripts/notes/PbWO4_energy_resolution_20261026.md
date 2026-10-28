# PbWO4 energy resolution: R3 adopted, per-run RF time calibration, best-50 spline templates

**Date:** 2026-10-26
**Branch:** `waveform_jz` (HEAD `3566071`, 10-25 14:14; everything since is uncommitted, §1.2)
**Plan:** `analysis/pyscripts/plans/plan_20261025.md` (Todo 28, 34–40; Progress log 10-26), `analysis/pyscripts/plans/plan_rf_calib_20261026.md`
**Previous notes:** `analysis/pyscripts/notes/PbWO4_energy_resolution_20261025.md` (R1, R2, hybrid map, module timing, RF calibration), `analysis/pyscripts/notes/PbWO4_energy_resolution_summary.md`

Paths are relative to `analysis/pyscripts/` in `/Users/L00338853/work/PRad/prad2evviewer` unless they start with `database/`, `docs/`, `prad2dec/`, `test/` or `analysis/src/` (repo root). σ/E, preferred fit, legacy/gated, Δ in pp and the run-type names are as in the summary note; legacy path, ±2σ Gaussian. R2 = `output/scan_20261025/pos_all`; R3 = R2 with the hybrid per-module map at 3.49 GeV H2, 2.24 GeV H2 and C (§2).

## 1. Status

Version 1, written 10-26; best-50 pipeline (reduced scope) finished 18:35 (benchmark jobs on `output/percrystal_best50_20261026/` and `output/scan_20261026_best50/`). Finished today: R3 adopted as the reference (§2); per-run RF time calibration (§3); `module_timing_tcal.pdf` regenerated (§4); DAQ memo (§5); RF period fixed in three repo-root comments (§6); CFD-fraction question answered (§7); flat-top diagnostic and the best-50 template bank (§8). §9.1 (pile-up) and §9.4 (σ_t per module and per cluster) filled from the reduced pipeline of 10-26 18:35; §9.2/§9.3 not run (reduced scope).

### 1.1 Code changes

| Script | Change (defaults unchanged unless stated) |
|---|---|
| `spline_template_extract.py` | `select`: `--mf-bank`, `--per-module-best K`, `--emax-frac F`, `--calib-file`, `--mf-strict`, `--max-modules-per-event N`, `--module-count-height-min H`, `--exclude-trigger-bits MASK`. `build`: `--pool-order {random,energy,energy-phase}`, `--iterations`, `--converge-tol`, `--converge-df`, `--max-iterations`, `--fit-points {bins,direct}`, `--lambda-select {fixed,heldout}`, `--lambda-min`, `--lambda-fallback`, `--lambda-chi2-sigma`, `--shape-sigma`, `--osc-escalate` (§8.3, §8.4). Default `select` (025320.00050, 025800.00100) and the Ta build reproduce the unmodified script and production (plan log 10-26 12:40). File last modified 15:56. |
| New: `rf_time_calib.py` | Per-run RF offsets (20261025 §9.7). |
| `plot_module_timing_rf.py` | `RF:prep` exposes the event mask `sel` (S0; outputs unchanged on `--runs 26138`). |
| Repo root, comments only | `docs/REPLAYED_DATA.md:201`, `prad2dec/include/TdcData.h:38` and `:95`, `test/evio_dump.cpp:685`: RF reference period ≈ 131.3 ns → 128.25 ns (§6). |
| Unchanged | Benchmark, driver, harness, `position_corr.py` (the hybrid-map edits are from 10-25), C++ code, `database/`. |

### 1.2 Git

| Item | State |
|---|---|
| HEAD | `3566071` (10-25 14:14, "energy_cali_benchmarking updates"). Nothing committed since. |
| Modified | `benchmark_hycal_timing_test_20261021.py`, `position_corr.py`, `spline_template_extract.py`, `notes/PbWO4_energy_resolution_20261022.md`, `notes/PbWO4_energy_resolution_20261025.md`, `notes/PbWO4_energy_resolution_summary.md`, `plans/plan_20261024.md`, `plans/plan_20261025.md`; repo root `docs/REPLAYED_DATA.md`, `prad2dec/include/TdcData.h`, `test/evio_dump.cpp`. |
| Untracked | `module_timing_dump.py`, `plot_module_energy.py`, `plot_module_timing.py`, `plot_module_timing_rf.py`, `plot_module_timing_vs_integral.py`, `plot_template_overlay.py`, `rf_time_calib.py`, `rf_time_dump.py`, `notes/daq_queries_20261026.md`, `plans/plan_rf_calib_20261026.md`, this note. `output/` is gitignored. |

## 2. R3 adopted as the reference

**Decision (user, 10-26).** R3 uses the hybrid per-module position map (`poscorr_hybrid_v1`, splits `output/scan_20261025/pos_hybrid`) at 3.49 GeV H2, 2.24 GeV H2 and C, and keeps the bilinear R2 map (`pos_all`) at 0.73 GeV H2 and Ta. At 0.73 GeV the hybrid map worsens σ/E (+0.191 ± 0.040 pp Gaussian, +0.127 ± 0.014 pp CB), and Ta has 0 of 40 own fits. No rescan: every job uses only its own run type's map. Details in 20261025 §4.9.

| Run type (fit) | R2 | R3 | Δ in-sample (paired) | Δ held-out | Module median R2 → R3 |
|---|---|---|---|---|---|
| 3.49 GeV H2 (CB) | 1.610 | **1.492** | −0.118 ± 0.007 | −0.097 ± 0.007 | 1.582 → 1.519 |
| 2.24 GeV H2 (Gauss) | 1.898 | **1.854** | −0.044 ± 0.002 | −0.035 ± 0.002 | 1.879 → 1.842 |
| 2.24 GeV C (CB) | 1.831 | **1.797** | −0.034 ± 0.005 | −0.016 ± 0.004 | 1.827 → 1.810 |
| 0.73 GeV H2 (Gauss) | 3.718 | **3.718** | 0 (R2 map) | – | 3.497 (unchanged) |
| Ta X17 (CB) | 1.978 | **1.978** | 0 (R2 map) | – | 1.908 (unchanged) |

Module medians: `--group-runs`, N > 1000 (`output/module_resolution_20261025_hybrid/summary/types_summary.tsv`; the R2 file at 0.73 GeV and Ta).

**Decomposition refit.** Source `output/resolution_fit_20261026/out_r3.txt`. It was made with a copy of `phase0.py` that adds the scan key r3 (per-run symlinks to the `pos_hybrid` / `pos_all` splits), and it reproduces the R2 rows exactly. σ² = a²/E + b² (+ c²/E²), with E = 0.73 / 2.24 / 3.49 GeV, H2 only, SEM-only errors.

| Scan / estimator | σ/E % at 0.73 / 2.24 / 3.49 | 2-term a, b (χ², 1 dof) | 3-term a, b, c |
|---|---|---|---|
| R2, CB for all | 2.974 / 1.858 / 1.610 | 2.37 ± 0.02, 0.98 ± 0.02 (3.56) | 2.25, 1.05, 0.64 |
| **R3, CB for all** | 2.974 / 1.806 / 1.492 | **2.49 ± 0.02, 0.68 ± 0.03 (3.42)** | c² = −0.40: unphysical |
| R2, preferred (±2σ) | 3.718 / 1.898 / 1.610 | 3.12, −0.78 (616) | 1.77, 1.16, 2.08 |
| R3, preferred (±2σ) | 3.718 / 1.854 / 1.492 | 3.18, −0.98 (374) | 2.26, 0.69, 1.82 |

- CB for all stays the consistent estimator (20261025 §8). R2 → R3 lowers the 2-term constant term from 0.98 to 0.68 % and raises a from 2.37 to 2.49 %·√GeV.
- The R3 3-term CB fit has c² < 0, so the constant-term share of σ² at 3.49 GeV stays the R2 value (43 %).
- The 2.24 and 3.49 GeV inputs carry the in-sample hybrid gain; the held-out gains are smaller by 0.01–0.02 pp.

## 3. Per-run RF time calibration

Executed from `plans/plan_rf_calib_20261026.md` (S0–S7). Script `rf_time_calib.py`, output `output/rf_calib_20261026/` (13 JSON, `rf_calib_modules.tsv`, `rf_calib_runs.tsv`, `rf_calib_validation.txt`, 13 PDFs, `run.log`; 80 s with 12 jobs). Method and full tables: 20261025 §9.7.

- **Gates.** (ii) PASS: held-out median |dhalf| 0.0167 ns, 0.9985 below 0.30 ns over the 13 011 of 13 254 measured (run, module) pairs with dhalf defined. (v) PASS: all 33 central + halo X17 modules firm. (vi) PASS: circ_mad 0.020–0.040 ns (pass B); the pooled fit reproduces `rf_summary.tsv` to 0.00000 ns.
- **Usable runs.** 12 of 13. 24386 is unusable: 78.2 % of its measured modules are bimodal (every other run ≤ 5.5 %), and its JSON has `usable: false` and empty `modules`. 24503 is flagged `drift` (split span 0.133 ns; the other 11 runs ≤ 0.044 ns).
- **Global phase jump.** 25800 → 26230: +1.990 ns over 1142 modules. 80 of 83 boards are within ±0.1 ns; the exceptions are crate 3 slots 15 / 8 / 18 (+1.578 / +0.322 / +0.143 ns from the global value). These are the boards where the shifted central modules are about half or more of the measured modules. The phase change is global, not per board.
- **Central block.** The 26 central X17 modules are shifted in 26230 only (D² −2.10 … −0.85 ns; in 25800 and 26138 the other 29 of the 33 agree with the prior). W525, W528, W597 and W697 are early in all three X17 runs (−0.89 … −1.35 ns at 25800 and 26138, −1.10 … −2.00 ns at 26230).
- **Open.** (a) Gate (ii) denominator: 0.9985 over pairs with dhalf defined, 0.9802 over all measured pairs (below 0.99); which reading applies is for the user. (b) Adoption into `database/hycal_rf_offsets/` is blocked by the time definition. The offsets centre the Python matched-filter `p_time`, while `analysis/src/Replay.cpp:1093` applies `hycal_rf_offsets` to the WaveAnalyzer CFD time of the cluster hit. The format is compatible, but the values do not transfer without a constant per module type (plan Open Question 1).

## 4. `module_timing_tcal.pdf` regenerated

`plot_module_timing.py --ref veto --time-calib` on the counter-aligned code (381 pages; companions `module_timing_veto_tcal.tsv`, `module_timing_veto_tcal_summary.txt`, `veto_summary_tcal.tsv`, log `run_veto_tcal.log`; 20261025 §9.2).
- X17 (t_module − t_seed, 374 modules): median σ 2.640 → 1.036 ns, and modules with |mean − median| > 3 ns 74 → 0. The seed-referenced offset blocks of 20261025 §9.5 are the database offsets themselves. 19 modules remain outliers in mean (|mean − median| > 5 × MAD-σ, 0.60–1.26 ns).
- Veto-referenced H2: median σ 0.630 → 0.623 ns (0.73 GeV, 93 modules) and 0.622 → 0.618 ns (2.24 GeV, 39 modules). The module-referenced alignment already absorbs the offsets.

## 5. DAQ memo

`notes/daq_queries_20261026.md`, three topics, 10 questions:
1. **Run 24386:** RF phase drift between splits (0.88 ns span), two phases 1.33 ns apart within splits, veto peak at +18.49 ns, 78 % bimodal modules; asks for restarts, resyncs or log entries.
2. **Veto V1–V4:** cable delays (V1 − V4 = −2.13 … −2.36 ns fitted), the two V1 timing populations at 0.73 GeV, no coincidence peak in 25228, 25257, 25260; asks for delays, changes and any offsets outside `time_calib.json`.
3. **FADC clock vs RF:** bunch period 4.00755 ± 0.00013 ns (59.5σ from 4.000 ns: the FADC clock is not RF-locked), the +1.99 ns global jump 25800 → 26230 and the crate 3 slot 15 / 8 / 18 exceptions; asks for the clock topology and changes at 26230.

## 6. RF period in the docs

Comments only, no build: ≈ 131.3 ns ≈ 7.61 MHz → 128.25 ns (32 × 4.0077 ns RF buckets) ≈ 7.80 MHz in `docs/REPLAYED_DATA.md:201`, `prad2dec/include/TdcData.h:38` and `:95`, and `test/evio_dump.cpp:685`.

## 7. CFD fraction (question answered)

- The CFD fraction 0.5 is hard-coded in `prad2dec/src/WaveAnalyzer.cpp:523` (`const float cfd_fraction = 0.5f;`).
- It enters only the pulse time, never the energy, in both the raw-sum and matched-filter paths.
- The matched-filter time (fitted onset + template t50) has no fraction parameter.
- A fraction scan is proposed; it has not been run.

## 8. Spline template extraction

### 8.1 Algorithm (`spline_template_extract.py`: cut-flow note `CUTFLOW_SEQ`, `--help`)

- **select** (per split): PbWO4 candidates pass the cut flow c1 (clip) → c5 (A ≥ 50 ADC) → c4 (contained) → c3a (first threshold) → c8 (pedestal) → c9 (tail) → c2 (no secondary matched-filter pulse) → c3b (matched-filter onset). By default the first `--per-module-cap` per module and split are stored.
- **build** (per run type): each module's pool is aligned on its fitted t0, and a spline with roughness penalty λ ∫ f″² dt is fitted to the 0.5 ns bin medians (default) or, with `--fit-points direct` as in the best-50 bank, to the pooled points themselves. The worst 20 % of the pool are dropped in iterations 1 and 2 (pool 50 → 40), the pool is frozen from iteration 3; `summary.tsv` records the fallback class 0/1/2 per module. `build` stores a per-run-type f(A); `ampcorr` then writes one global f(A) if the run types agree (within 2 combined SE), otherwise each keeps its own.

### 8.2 Flat-top diagnostic

Source `output/templates_spline_diag_20261026/report.md` (read-only; 2.24 GeV H2 and Ta production banks; the 48 rebuilt pools reproduce production exactly).
- **Alignment is not the cause.** Residual misalignment is ≤ about 0.1 ns (final ±0.5 ns pass moves pulses by 0.07–0.08 ns rms); 0.15 ns of jitter would give about 18 times the fitted rise excess. The closure test recovers +0.370 ns for δ = +0.37 ns.
- **λ = 1e4 on small pools.** Of the modules with w99 > 2.5 ns (17 H2, 23 Ta), 82 % / 91 % have N0 < 100. In 8 of the 12 flattest H2 modules (median N0 29) and 11 of 12 Ta (N0 27), w99 grows from 0.90 to 2.90 ns (H2) and from 0.98 to 3.03 ns (Ta) between λ = 1e2 and 1e4.
- **Flat in the data.** W18, W704, W885, W1148 are flat at every λ (w99 2.6–3.6 ns; bin-median curvature −0.002 to −0.019 per ns² against a typical −0.031). Ta W325 has two humps 1.8 ns apart.
- Median w99 is 1.66 ns (production) against 2.12 ns for the h50 two-τ shape.
- **LMS (bit 24) in the production banks:** 0 candidates in H2 and C (exact `trigger_bits == 0x100`), 5 in Ta (3.0e-6; `bits & 0xf00`).

### 8.3 Best-50 selection

The options add a per-module selection by energy:
- `--per-module-best 50` keeps, per module and split, the 50 candidates with the highest E = mf_amp · I_int · cal, with E < 1.2 E_beam (`--emax-frac 1.2`; calibration from `module_timing_dump.RUN_CFG` unless `--calib-file`).
- `--mf-strict` requires exactly one matched-filter pulse.
- `--max-modules-per-event 20` rejects events with more than 20 PbWO4 channels at A ≥ 50 ADC.
- `--exclude-trigger-bits 0x1000000` drops LMS events.
- `build --pool-order energy` fits the top 50 by energy and holds out the next 13 per module for `check` (`holdout_keys`).

All 119 splits were selected into `output/templates_spline_best50_20261026/<runtype>/pulses/` with the production spline bank as `--mf-bank` (`logs/select_*.log`).

| Item | 0.73 GeV H2 | 2.24 GeV H2 | 3.49 GeV H2 | C | Ta 25800 / 26138 / 26230 |
|---|---|---|---|---|---|
| `--mf-strict` rejection (of c3b survivors) | 15.6 % | 12.9 % | 9.9 % | 11.9 % | 21.3 / 17.8 / 21.8 % |
| `--max-modules-per-event 20` (events) | 0.0008 % | 0.0010 % | 0.0044 % | 0.0013 % | 0.27 / 0.20 / 2.52 % |
| Bit-24 events rejected | 0 | 0 | 0 | 0 | 2 / 0 / 0 |
| Modules with 50 stored, per split (mean of 1152) | 226 | 610 | 404 | 366 | 305 / 314 / 288 |

- **Strict needs `--mf-bank`.** With the default h50 two-τ bank nearly every A ≥ 100 ADC pulse has two or more matched-filter pulses (template mismatch; `--help`). The 10-26 smoke (10-26 session diagnostic) rejected 87–91 % of c3b survivors with h50, against the 10–22 % above with the spline bank.
- **Multiplicity.** The cut is negligible in H2 and C and in 25800 / 26138; in 26230 (per-event count mean 11.7, 99th percentile 22–23) it removes 2.3–2.7 % of events.

### 8.4 Pass-1 oscillation and the escalation

- **Pass 1** (`pass1_osc/`, λ grid {0, 1e1, 1e2, 1e3, 1e4} selected per module by held-out χ²): 33–46 % of the module templates per run type had spurious wiggles (10-26 session diagnostic).
- **S2′ grid** (2.24 GeV H2 and Ta, `grid/`; npz `meta.argv`). Common (d overrides `--pool-order`): `--pool-order energy --n-max 50 --lambda-select heldout --lambda-chi2-sigma full --lambda-fallback 1e3 --shape-sigma 0.012 --converge-tol 0.01 --converge-df 5e-3 --max-iterations 10`. The full tables go into 20261025 §4.10, which is being written while this note is drafted. Modules at λ = 1e2 / 1e3 / 1e4, 2.24 GeV H2 (`grid/<v>/2p24GeV_H2/build.log`):
  - a: `--fit-points direct --lambda-min 1e2`; 179 / 542 / 430.
  - b: as a with `--lambda-min 1e3`; 0 / 664 / 487.
  - c: as a with `--fit-points bins`; 295 / 520 / 336 (converged 471 of 1151, against 778 for a).
  - d: as a with `--pool-order energy-phase`; 173 / 567 / 411.
- **Choice.** The bank takes variant a plus `--osc-escalate` (`a_noesc/` is a without escalation for all five run types). Escalation multiplies λ by 10 and rebuilds the pool while a module has extrema of prominence ≥ 0.2 % of the peak outside ±3 ns (2–150 ns); a module still wiggly at λ = 1e4 takes the run-type average (`osc_fallback`).
- **Result.** 0.0 % of the module templates oscillate at 0.2 % prominence (10-26 session diagnostic). Per run type (`<runtype>/build.log`):

| Run type | Flagged before | osc_fallback | Converged (≤ 10 iterations) | Fallback 0 / 1 / 2 |
|---|---|---|---|---|
| 0.73 GeV H2 | 175 | 27 | 811 / 1088 | 607 / 481 / 64 |
| 2.24 GeV H2 | 96 | 24 | 773 / 1127 | 1059 / 68 / 25 |
| 2.24 GeV C | 132 | 19 | 729 / 1034 | 753 / 281 / 118 |
| Ta X17 | 156 | 20 | 767 / 1089 | 840 / 249 / 63 |
| 3.49 GeV H2 | 184 | 42 | 718 / 1065 | 869 / 196 / 87 |

At 3.49 GeV, 42 of the 87 fallback-2 modules are `osc_fallback`; the pass-1 bank had 45 fallback-2 modules.

### 8.5 Best-50 bank vs production

Bank `output/templates_spline_best50_20261026/<runtype>/templates_<runtype>.npz`, moved from `esc_stage/` at 15:36 (the build logs, `ampcorr.log` and the npz argv still name `esc_stage/`). `ampcorr` scope is per run type (`ampcorr.log`).

| Gate (`check/gates.tsv`) | Production (`output/templates_spline/`) | Best-50 |
|---|---|---|
| G1-4 held-out χ²/dof ≤ 0.5 × two-τ | PASS (5/5) | PASS (5/5) |
| G1-5 sub-sample rise | PASS | PASS |
| G1-6 I_int vs the C++ integral | **FAIL** (5/5: 3–5 modules with phase-mean > 0.5 %, 33–73 no-peak) | **PASS** (5/5: 0 over 0.5 %, 0 no-peak) |
| G1-7 deficit flatness | FAIL | FAIL (0.73 GeV: N/A, only the 300–1000 bin is populated) |
| G1-8 coverage (fallback ≤ 1 ≥ 1100) | PASS / REPORT | 3.49 GeV **FAIL** (1065; production 1126) |
| G1-9 amplitude dependence | FAIL | FAIL |

- 10-26 session diagnostic: w99 within +0.11 ns of production, rise 0.15–0.2 ns faster, held-out χ² lower than production on the same pulses.
- `pass2/` holds a re-selection of 2.24 GeV H2 and Ta with the variant-a bank without escalation (now `a_noesc/`; sha1 in `pass2/logs/`) as `--mf-bank`, built 15:40; not evaluated here.

## 9. Best-50 pipeline results (reduced scope)

Gains `output/percrystal_best50_20261026/`, scans `output/scan_20261026_best50/` (runs 24503, 24590, 25308, 25320, 25800; `control/`).

Reduced scope of `plans/plan_best50_pipeline_20261026.md` (user decision "do reduced"): S0 code, S1 pile-up, S5 module timing, S6 cluster timing; S2–S4 (gains for nine runs, 118-split scan, module resolution) not run. Every output is under `output/best50_20261026/`; tables in `compare/`. Bank B50 = `output/templates_spline_best50_20261026/`, production bank T0 = `output/templates_spline/`. Code added for this section (defaults unchanged; this supersedes the "Unchanged: Benchmark" entry of §1.1): `benchmark_hycal_timing_test_20261021.py --cluster-times` (with `--write-cluster-tsv`: four columns `evnum seed_time_ns n_members members` appended to `result_clusters.tsv`; flag off byte-identical on 025308 split 00100, `scan/identity/`); `module_timing_dump.py --templates-root`; `rf_time_calib.py --dump-dir --rf-dir --xcheck` (defaults reproduce `output/rf_calib_20261026/` byte for byte except `created`); new `cluster_timing.py`, `best50_compare.py`.

### 9.1 Pile-up rejection table

Detection-efficiency grid of `benchmark_pileup_mf.py` (the table of 20261023 §9.4; no table is literally a rejection efficiency): `eff_match` = % of 50 events per Δt × amplitude-ratio cell with an mf pulse within ±4 ns of both onsets, MIN_SEP_NS 12. Arms: A = production synth × T0 (`output/pileup_eff_20261024_fix/`); **B = production synth × B50** (same pseudo-data, only the reconstruction bank changes); C = synth regenerated with B50 × B50; D = that synth × T0. The C/D synth have A's 3000 base events per run type (base lists array-identical). Pooled Δ = mean over cells of (arm − A) in points, headline block r_model 0.3:

| Run type | B all 60 cells | B cells with A < 100 % | C all | D all | B cells with Δ < 0 (A → B, %) |
|---|---:|---:|---:|---:|---|
| 0.73 GeV H2 | +0.37 | +0.59 (41) | +0.30 | +0.17 | Δt 32 ratio 0.3: 2 → 0; Δt 32 ratio 0.5: 100 → 98; Δt 20 ratio 0.5: 2 → 0 |
| 2.24 GeV C | +0.63 | +0.95 (40) | +0.47 | +0.13 | none |
| 2.24 GeV H2 | +0.67 | +1.00 (40) | +0.67 | +0.23 | none |
| 2.24 GeV Ta (X17) | +0.07 | +0.07 (56) | +0.73 | +1.13 | Δt 24 ratio 0.5: 4 → 0 |
| 3.49 GeV H2 | +0.13 | +0.25 (40) | +0.63 | +0.47 | Δt 32 ratio 0.5: 100 → 98; Δt 20 ratio 1: 42 → 40 |

With the best-50 bank the grid is unchanged or better: the gains are at Δt 20–24 ns, ratio 1 (2.24 GeV H2 Δt 20 ns: 26 → 52 %, Δt 24 ns: 90 → 100 %; +24 at 2.24 C, +12 at 0.73), the losses are single cells of one or two events out of 50. `fake_rate` changes by ≤ 0.24 per event in any cell; `n2_wa` of B equals A cell by cell. Full grids (all four blocks, `fake_rate`, `n2_mf`): `compare/pileup.md`, cells `compare/pileup_cells.tsv`.

### 9.2 Full-scan σ_E/E
Not run — reduced scope, 4-split result in 20261025 §4.10.

### 9.3 Per-module σ_E/E
Not run — reduced scope, 4-split result in 20261025 §4.10.

### 9.4 σ_t per module and per cluster

**Per module.** Module timing dump of the 13 RF runs (111 splits) with B50 (`module_timing/dump/`, 0 failed, `templates` metadata under B50), `plot_module_timing.py --ref veto`, `plot_module_timing_rf.py`, and `rf_time_calib.py` on the new dump (cross-check = the new `rf_summary.tsv`). Gates (ii), (v), (vi) PASS; 24386 bimodal (unusable) and 24503 drift as with T0. Medians over modules (ns), B50 / T0 / median per-module Δ:

| Run type | RF calibration per run, E ≥ 100 MeV (tier e100 in both) | RF, runs aligned, all energies (`module_timing_rf.tsv`) | Veto (`module_timing_veto.tsv`) |
|---|---|---|---|
| 0.73 GeV H2 | 0.173 / 0.187 / −0.014 (569) | 0.517 / 0.512 / −0.002 (1150) | 0.620 / 0.630 / −0.003 (93) |
| 2.24 GeV H2 | 0.172 / 0.189 / −0.020 (4773) | 0.461 / 0.461 / −0.006 (1151) | 0.618 / 0.622 / −0.005 (39) |
| 3.49 GeV H2 | 0.157 / 0.182 / −0.026 (1646) | 0.542 / 0.545 / −0.012 (1132) | no fit |
| 2.24 GeV Ta (X17) | 0.185 / 0.203 / −0.020 (944) | 0.650 / 0.647 / −0.001 (1150) | 2.650 / 2.640 / +0.030 (374) |

(Number of modules or (run, module) pairs in parentheses.) The per-run module offsets (tier e100 in both) move by a median +0.04 to +0.07 ns (folded); the tier changes for ≤ 0.3 % of modules per run. Tables: `compare/rfcalib.md`, `compare/timing_rf.md`, `compare/timing_veto.md`.

**Per cluster.** For each single-cluster event of `result_clusters.tsv` (legacy path, `--cluster-times`): seed time t_seed − o_seed, and energy-weighted time Σ E_i (t_i − o_i) / Σ E_i over the members, with o the per-run module offsets of the same bank (B50: `rf_calib/`; T0: `output/rf_calib_20261026/`). Each member's t_i − o_i is first unwrapped to the RF bucket of the seed's corrected time (51 % of members are in another bucket); E_i are the member energies before the per-crystal gains. φ = fold(t − t_RF) with t_RF the nearest RF channel-A hit of the same TI event (the matching of `RF.prep`), σ_t = `fit_unwrapped` on the 80-bin circular histogram (the module-level estimator). One split per RF run (13 benchmark jobs per bank: R3 `job.cmd` of the first split with an RF dump, `--cluster-times` appended, B50 arm with only `--mf-templates` changed; both arms with the `percrystal_spl` gains and the R3 maps); 24386 skipped (unusable RF). Members with an offset: 100.0 % in both arms. Pooled per run type, tier all (ns):

| Run type | N (B50 / T0) | seed B50 | seed T0 | Δ | E-weighted B50 | E-weighted T0 | Δ |
|---|---:|---:|---:|---:|---:|---:|---:|
| 0.73 GeV H2 | 360 779 / 360 593 | 0.167 ± 0.002 | 0.162 ± 0.002 | +0.005 | 0.143 ± 0.001 | 0.139 ± 0.001 | +0.003 |
| 2.24 GeV H2 | 657 895 / 657 027 | 0.107 ± 0.002 | 0.135 ± 0.002 | −0.028 | 0.096 ± 0.001 | 0.110 ± 0.001 | −0.014 |
| 3.49 GeV H2 | 238 455 / 238 402 | 0.097 ± 0.002 | 0.143 ± 0.003 | −0.046 | 0.089 ± 0.002 | 0.113 ± 0.001 | −0.024 |
| 2.24 GeV Ta (X17) | 19 739 / 19 729 | 0.127 ± 0.005 | 0.146 ± 0.003 | −0.020 | 0.114 ± 0.003 | 0.118 ± 0.002 | −0.004 |

Member-count tiers: ≥ 5 members hold 72 % (0.73 GeV) to > 98 % (other types) of the clusters and carry the values above. Tier 2–4: 0.73 GeV seed 0.178 / 0.176, E-weighted 0.157 / 0.153; 2.24 GeV H2 seed 0.180 / 0.204, E-weighted 0.156 / 0.177 (B50 / T0). Single-member clusters fit only at 0.73 GeV (0.213 / 0.221, N 2243); the fit fails at 2.24 GeV H2 (N 1168), and at 3.49 GeV and X17 tiers 1 and 2–4 have < 300 clusters. Per-run tables and tiers: `compare/cluster_timing.md`; `cluster_timing/{new,r3_ref}/cluster_timing.{tsv,pdf}`.

## 10. Plan state (`plans/plan_20261025.md`)

| Todo | State |
|---|---|
| 28 R3 | Done, checked (log 10-26). |
| 34 RF calibration per run | Done, checked (log 10-26 10:30). |
| 35, 36 run 24386, veto hardware | Memo written (§5); still open in the plan. |
| 37 doc fix 131 → 128.25 ns | Done in the working tree (§6); unchecked in the plan. |
| 38 `module_timing_tcal.pdf` | Done (§4); unchecked in the plan. |
| 39 Ta-only gate re-scan | Open; waits for new X17 data. |
| 40 best-K templates | Options added (log 10-26 12:40), bank built (§8); reduced pipeline done (§9.1, §9.4); full rescan not run. |

## 11. Caveats and open items

| Item | Detail |
|---|---|
| R3 in-sample | The hybrid maps are fitted on the splits they are applied to; the held-out gains are smaller by 0.01–0.02 pp (§2). |
| R3 decomposition | 3-term CB fit unphysical (c² = −0.40); the 43 % constant-term share is R2-based. |
| RF offsets | Python matched-filter time definition, not `Replay.cpp` CFD time; not in `database/`. Gate (ii) denominator open. 24386 unusable, 24503 drift. |
| Best-50 coverage | 3.49 GeV fails G1-8: 87 fallback-2 modules, against 26 in production. |
| Held-out sets differ | Production `check` uses the hash hold-out, best-50 the next 13 by energy, so the G1-4 χ² values in the two `gates.tsv` are not on the same pulses. |
| Session numbers | Pass-1 oscillation (33–46 %), 0.0 % after escalation, h50 strict rejection (87–91 %), w99, rise and held-out χ² vs production come from the 10-26 session diagnostic and have no project file. |
| Ta multiplicity | 26230 loses 2.5 % of events to `--max-modules-per-event 20`; 25800 / 26138 0.20–0.27 %. |
| Best-50 σ_t per cluster (§9.4) | One split per RF run, both banks with the `percrystal_spl` gains (best-50 gains not derived for these runs); the E-weighted time uses pre-gain member energies; σ ≈ 0.1 ns is two bins of the 80-bin (0.05 ns) estimator. The plan's sanity rule (seed σ_t within 0.05 ns of the `module_timing_rf.tsv` N-weighted median σ, 0.46–0.67 ns) does not hold and does not apply: that σ pools all energies with one offset per run; the per-run E ≥ 100 MeV medians are 0.16–0.20 ns. |

## 12. How to run

```bash
cd /Users/L00338853/work/PRad/prad2evviewer/analysis/pyscripts && source setup_env.sh
E=/Users/L00338853/work/PRad/data/evio; T0=output/templates_spline; BK=output/templates_spline_best50_20261026

# §2 R3 decomposition refit (copy of phase0.py with scan key r3)
(cd output/resolution_fit_20261026 && python3 phase0.py --scans spl,r2,r3 --also-fixed150 > out_r3.txt)

# §3 RF time calibration
MPLBACKEND=Agg python3 rf_time_calib.py --jobs 12 2>&1 | tee output/rf_calib_20261026/run.log

# §4 module_timing_tcal.pdf
MPLBACKEND=Agg python3 plot_module_timing.py --ref veto --time-calib

# §8.2 flat-top diagnostic
MPLBACKEND=Agg python3 output/templates_spline_diag_20261026/diag.py --banks 2p24GeV_H2,2p24GeV_Ta_X17 --nflat 12 --nrand 12 --seed 20261026

# §8.3 best-50 select, per split (H2/C: --trigger physics; Ta: --trigger physics-any), RT = run type
python3 spline_template_extract.py select $E/prad_025320.evio.00050 --out $BK/2p24GeV_H2/pulses/025320_00050.npz --trigger physics \
  --per-module-best 50 --emax-frac 1.2 --mf-strict --mf-bank $T0/2p24GeV_H2/templates_2p24GeV_H2.npz \
  --max-modules-per-event 20 --exclude-trigger-bits 0x1000000
# §8.4 build per run type: options from the npz meta.argv; --out at the post-move path (argv: esc_stage/)
python3 spline_template_extract.py build --pulses "$BK/$RT/pulses/*.npz" --run-type $RT --pool-order energy --n-max 50 --fit-points direct \
  --shape-sigma 0.012 --lambda-select heldout --lambda-chi2-sigma full --lambda-min 1e2 --lambda-fallback 1e3 --converge-tol 0.01 \
  --converge-df 5e-3 --max-iterations 10 --osc-escalate --out $BK/$RT/templates_$RT
python3 spline_template_extract.py ampcorr --out $BK/ampcorr.json --templates $BK/*/templates_*.npz
python3 spline_template_extract.py check --pulses "$BK/$RT/pulses/*.npz" --templates $BK/$RT/templates_$RT.npz --out $BK/$RT/check/

# §9 best-50 pipeline, reduced scope (O=output/best50_20261026; MPLBACKEND=Agg); SL = run-type slug
python3 benchmark_pileup_mf.py output/pileup_eff_20261024_fix/$SL/synth --templates $BK/$SL/templates_$SL.npz --out $O/pileup/B_prodsynth_newbank/$SL --tolerance-ns 4 --window-rule 40:0.16
python3 pileup_generator.py <base split of output/pileup_eff_20261024_fix/$SL/synth/manifest.json> --template $BK/$SL/templates_$SL.npz --out-dir $O/pileup/synth_newbank/$SL/synth \
  --material PbWO4 --n-per-config 50 --height-min 500 --t0-min 22 --residual-veto-sigma 999 --noise-model none --post-pad 20 --seed 12345 --force
python3 benchmark_pileup_mf.py $O/pileup/synth_newbank/$SL/synth --templates $BK/$SL/templates_$SL.npz --out $O/pileup/C_newsynth_newbank/$SL --tolerance-ns 4 --window-rule 40:0.16
python3 benchmark_pileup_mf.py $O/pileup/synth_newbank/$SL/synth --templates $T0/$SL/templates_$SL.npz --out $O/pileup/D_newsynth_prodbank/$SL --tolerance-ns 4 --window-rule 40:0.16
python3 best50_compare.py pileup --ref output/pileup_eff_20261024_fix --new $O/pileup/B_prodsynth_newbank --arms C=$O/pileup/C_newsynth_newbank D=$O/pileup/D_newsynth_prodbank --out $O/compare/pileup
python3 module_timing_dump.py --all --parallel 10 --templates-root $BK --out $O/module_timing/dump
python3 plot_module_timing.py --ref veto --dump-dir $O/module_timing/dump --out-dir $O/module_timing
python3 plot_module_timing_rf.py --dump-dir $O/module_timing/dump --rf-dir output/module_timing_20261025/rf --out-dir $O/module_timing
python3 rf_time_calib.py --dump-dir $O/module_timing/dump --rf-dir output/module_timing_20261025/rf --xcheck $O/module_timing/rf_summary.tsv --out-dir $O/rf_calib --jobs 12
python3 best50_compare.py timing --new $O/module_timing/module_timing_rf.tsv --ref output/module_timing_20261025/module_timing_rf.tsv --out $O/compare/timing_rf
python3 best50_compare.py timing --new $O/module_timing/module_timing_veto.tsv --ref output/module_timing_20261025/module_timing_veto.tsv --out $O/compare/timing_veto
python3 best50_compare.py rfcalib --new $O/rf_calib/rf_calib_modules.tsv --ref output/rf_calib_20261026/rf_calib_modules.tsv --out $O/compare/rfcalib
# cluster jobs (13 per arm, 6 parallel): output/resolution_fit_20261026/r3_splits/scan_<run>/split_<n>/job.cmd with argv[3] -> $O/cluster_timing/<arm>/scan_<run>/split_<n>/result,
#   --cluster-times appended, arm new: --mf-templates $T0/... -> $BK/...; the commands as run are in $O/cluster_timing/<arm>/scan_*/split_*/job.cmd
python3 cluster_timing.py --scan-dir $O/cluster_timing/new --rf-dir output/module_timing_20261025/rf --offsets-dir $O/rf_calib --out-dir $O/cluster_timing/new --jobs 6
python3 cluster_timing.py --scan-dir $O/cluster_timing/r3_ref --rf-dir output/module_timing_20261025/rf --offsets-dir output/rf_calib_20261026 --out-dir $O/cluster_timing/r3_ref --jobs 6
python3 best50_compare.py cluster --new $O/cluster_timing/new --ref $O/cluster_timing/r3_ref --out $O/compare/cluster_timing
```
