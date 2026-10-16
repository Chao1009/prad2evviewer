# Waveform Analysis — Calibration Notes

Analysis notes specific to calibration measurements arising from the
waveform-analysis arc. For the overall workstream plan and phase
structure, see [`docs/plans/plan_waveform.md`](../plans/plan_waveform.md).

---

## Calibration Investigation: Run 24386 vs 25320 (Oct 2026)

**Status:** COMPLETE

### Findings

- Run 24386 (11 EVIO files): σ_CB/E = 2.0% (physics_calib iter1), clean single-peak spectrum
- Run 25320 (88 EVIO files): σ_CB/E = 3.75% (calibration_factor_2p2_fadc.json baseline), double-peak contamination present
- Physics_calib methodology is mathematically sound (convergence ratio_stdev → 0.0018 on 25320 full data)
- Per-module statistics with 11–12 EVIO files are underconstrained; physics_calib iter 1–7 converge but to suboptimal local minimum
- 25320 double-peak structure is data-specific (trigger configuration or beam instability), not an algorithm issue

### Technical Details

- 25320 iter1 convergence: 652 modules attempted, 474 good fits (72.7%), ratio_stdev = 0.0148
- 25320 iter7 convergence: 651 modules calibrated, 474 good fits, ratio_stdev = 0.0018
- 25320 iter7 benchmark result: σ_CB/E = 5.89% (worse than 3.75% baseline), indicating calibration divergence due to limited statistics
- Seed calibration mismatch identified: old calibration_factor_3p5_June7_nonLinear.json had zeroed non-linearity coefficients; corrected to new_calibration_3p5_June7_nonLinear.json, but did not improve iter7 result on 25320

### Decision

- **Use 24386 as reference calibration:** Achieves 2.0% σ_CB/E and validates hycal_clustering.md baseline
- **Use 25320 baseline:** calibration_factor_2p2_fadc.json gives 3.75%, acceptable for analysis despite contamination flag
- **25320 analysis caveat:** Flag for potential pile-up or trigger contamination when used in downstream physics
- **Future calibration work:** Full-run physics_calib (50+ EVIO files per run) required for stable convergence

### Reference Files

- Validated 24386 iter1 calibration: `/Users/L00338853/work/PRad/data/physics_calib_out/Physics_calib/run24386/calib_factor_iter1.json`
- 25320 working baseline: `database/calibration/calibration_factor_2p2_fadc.json`
- Pulse template comparison: `analysis/pyscripts/output/pulse_templates_h50.dat`

---

## Timing-Gate Benchmark Study on Run 24386 (Oct 2026)

**Status:** COMPLETE

### Objective

Validate optimal pre-window width `[LO, HI]` for run 24386 using the validated iter1 physics_calib calibration. Sweep the lower edge `LO` across 7 points while fixing `HI = 200` ns.

### Implementation

- Script: `analysis/pyscripts/sweep_pre_window_24386.py`
- Plan: `docs/plans/plan-timing-gate-benchmark-24386.md`
- Sweep configuration: `LO ∈ {50, 75, 100, 125, 130, 150, 175}` ns, fixed `HI = 200` ns
- EVIO file: `prad_024386.evio.00050`, 50k events
- Fixed parameters: `--window 16 --calibrate --cal-method cb-peak --target-peak 3488.5 --exclude-rowcol 15 20 15 20`
- Calibration seed: `calib_factor_iter1.json` (24386 iter1)

### Deliverables

- Aggregate TSV: `analysis/pyscripts/output/timing_gate_benchmark_24386_summary.tsv` (14 rows: 7 LOs × 2 paths)
- Trend plots:
  - `analysis/pyscripts/output/timing_gate_benchmark_24386_curves_sigma_E.png`
  - `analysis/pyscripts/output/timing_gate_benchmark_24386_curves_acceptance.png`
  - `analysis/pyscripts/output/timing_gate_benchmark_24386_curves_mu.png`
- Per-window benchmark outputs: `analysis/pyscripts/output/timing_gate_benchmark_24386/`

### Key Findings

- 24386 achieves σ_CB/E = 2.0% at all tested pre-window configurations (variations < 0.1 pp)
- Pre-window width has minimal impact on energy resolution in the 25–150 ns range
- Narrowest window (LO=175, width=25 ns) begins to degrade statistics
- Runinfo baseline (LO=130, width=70 ns) is a reasonable production default

---

## Multi-Run EVIO Data Quality Scan (Oct 2026)

**Status:** COMPLETE

### Objective

Scan all available EVIO files across 7 runs to characterize data quality, identify secondary peaks in energy spectra, and compare σ_E/E across runs.

### Implementation

- Script: `analysis/pyscripts/benchmark_hycal_timing_scan.py`
- Input: All EVIO files in `/Users/L00338853/work/PRad/data/evio/`
- Processing: Full-file (no max_events cap), 8-way parallelism, per-file calibration
- Secondary peak detection: `scipy.signal.find_peaks` with prominence > 40, distance > 100 MeV from primary
- Secondary peak definition: any peak distinct from the primary (elastic) peak, regardless of position

### Runs Scanned

| Run | Beam Energy | Target | Files | Events/file (avg) |
|-----|-------------|--------|-------|-------------------|
| 24386 | 3.49 GeV | H2 | 11 | ~30k |
| 24503 | 3.49 GeV | H2 | 10 | ~105k |
| 24512 | 3.49 GeV | H2 | 7 | ~95k |
| 24590 | 0.73 GeV | H2 | 10 | ~140k |
| 25308 | 2.24 GeV | 1μm C | 8 | ~151k |
| 25320 | 2.24 GeV | H2 | 88 | ~131k |
| 26138 | 2.24 GeV | 1μm Ta (X17) | 4 | ~5k |

### Deliverables

- Aggregate TSV: `analysis/pyscripts/output/benchmark_hycal_timing_scan_summary.tsv` (138 rows)
- Trend plot: `analysis/pyscripts/output/benchmark_hycal_timing_scan_sigma_E_vs_file.png`
- Combined PDF (138 pages): `analysis/pyscripts/output/benchmark_hycal_timing_scan_all_calibrated_plots.pdf`
- Findings markdown: `analysis/pyscripts/output/benchmark_hycal_timing_scan_findings.md`
- Per-file outputs: `analysis/pyscripts/output/scan_NNNNN/split_NNNNN/`

### Secondary Peak Findings (Position-Agnostic Detection)

| Run | Beam E | σ_E/E (CB) mean | Secondary peaks | Secondary location | Interpretation |
|-----|--------|-----------------|-----------------|---------------------|-----------------|
| 24386 | 3.49 GeV | 1.98% | 4/11 (36%) | ~3660 MeV | Sporadic pile-up |
| 24503 | 3.49 GeV | 2.02% | 10/10 (100%) | ~1500 or ~3660 MeV | Inelastic + pile-up |
| 24512 | 3.49 GeV | 2.01% | 7/7 (100%) | ~1500 or ~3660 MeV | Inelastic + pile-up |
| 24590 | 0.73 GeV | 4.91% | 1/10 (10%) | ~680 MeV | Mostly clean |
| 25308 | 2.24 GeV | 2.10% | 8/8 (100%) | ~900 MeV | Møller / multi-cluster |
| 25320 | 2.24 GeV | 2.05% | 88/88 (100%) | ~1951 MeV | **Δ-excitation (physics)** |
| 26138 | 2.24 GeV | 2.39% | 0/4 (0%) | — | Clean |

### Key Interpretations

**Run 25320 "double peak" is physics, not contamination.** The consistent secondary at 1951 MeV (±3 MeV across 88 files, prominence 1600–1850 — 40× above threshold) is ~290 MeV below the elastic peak at 2241 MeV. This is the **Δ(1232) resonance excitation**: ep → eΔ → epπ⁰. Expected for H2 target at 2.24 GeV beam energy.

**Run 25308 secondary at 900 MeV** is far below the elastic peak and consistent with Møller scattering or multi-cluster artifacts in the carbon target data.

**Runs 24503 and 24512** show the ONLY above-primary peaks (at 3660 MeV, 130 MeV above the 3530 MeV elastic). This is the signature of **pile-up contamination** — two interactions summing into a single cluster. Both runs have 100% file-level contamination.

**Run 24386 is the cleanest reference at 3.49 GeV.** Only 4/11 files show pile-up, and the mean σ_E/E = 1.98% is the best resolution of any run.

### Recommended for Precision Physics Analysis

- 24386 (3.49 GeV, H2) — cleanest 3.49 GeV reference
- 25320 (2.24 GeV, H2) — clean + highest statistics; Δ-excitation is physics
- 26138 (2.24 GeV, Ta) — clean X17 run
- 24590 (0.73 GeV, H2) — mostly clean, low-energy run

### Flagged for Pile-up Contamination

- 24503 (3.49 GeV, H2) — pile-up in all files
- 24512 (3.49 GeV, H2) — pile-up in all files

### Technical Notes

- Known issue: `benchmark_hycal_timing.py` can hang on specific EVIO files (observed: 024386/split_00040 hung for 1h47m with 99% idle CPU). Workaround: kill the stuck child; orchestrator continues gracefully.
- The sweep script's `--skip-benchmark` flag allows re-aggregation of existing TSVs with revised peak-detection logic without rerunning the ~2-hour benchmark phase.

---

## Pulse Template Comparison (Oct 2026)

**Status:** COMPLETE

### Objective

Compare per-module tau_r and tau_f pulse template parameters between run 25308 (Carbon, 2.24 GeV) and run 25320 (H2, 2.24 GeV) PbWO4 modules.

### Implementation

- Script: `analysis/pyscripts/dump_pbwo4_table.py`
- Input templates:
  - `output/pulse_templates_025308_h50_RawSum_max20modules.json`
  - `output/pulse_templates_025320_SSP_RawSum_h50_max20modules.json`

### Deliverables

- Comparison table: `analysis/pyscripts/output/pulse_templates_h50.dat`
  - 1151 PbWO4 modules
  - 11 columns: module name, per-run tau_r/tau_f mean/sigma, relative differences

### Columns

- `module` — module name (W1 to W1156)
- `tau_r_mean_25308`, `tau_r_sigma_25308`, `tau_f_mean_25308`, `tau_f_sigma_25308`
- `tau_r_mean_25320`, `tau_r_sigma_25320`, `tau_f_mean_25320`, `tau_f_sigma_25320`
- `d_tau_r/(tau_r_25308)` — relative rise-time difference
- `d_tau_f/(tau_f_25308)` — relative fall-time difference

### Key Findings

- Mean τ_r relative difference: +1.18% (25308 ~1% higher than 25320)
- Mean τ_f relative difference: -0.24% (nearly identical)
- τ_r stdev across modules: ±4.49%
- τ_f stdev across modules: ±0.76% (very tight)
- One cut parameter differs: `t0_min = None` (25308) vs `22.0` (25320)

### Interpretation

Pulse shapes between the two runs are highly consistent at the module level. The small τ_r shift (~1%) is likely due to the different trigger configurations (25308 uses tighter time cut).
