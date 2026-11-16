# P15: calibrated-time gate scan

Date 2026-11-13. Package P15 of `notes/plan_energy_resolution_improvement_20261112.md`. Question: does subtracting the per-module time offsets before the time gate move the optimum gate or the resolution? Decisions: the offset correction is applied before the gate scan (user, 11-13); the best gate per run type is shown, not a fixed 12 ns.

## Method

- New switch `--set TIME_OFFSETS=True` in `energy_resolution_benchmark.py` (default False): at report time the stored raw hit time is replaced by `t - offset[module]` before every `add_hit` (legacy and gated paths). Table `$PRAD2_DATABASE_DIR/hycal_time_offsets/time_calib.json`, resolved rows only (calibrated = raw - offset_ns), indexed by `Module.index`; lead glass and unresolved rows 0. No rescan (the scan stores raw `h_t`); C++ `HyCalCluster` unchanged.
- Arm as the P0 baseline: gains + position correction, Moller rejection and non-linearity off. Gates 1, 2, 3, 4, 5, 6, 8, 10, 12, 16 ns plus legacy for five run types, extended to 20 and 24 ns for the four H2 run types (Ta X17 stops at 16 ns; its minimum, 8 ns, lies inside). Output `output/energy_resolution/variants/P15/<run type>/` with `logs/` and `run.sh`; tables by `improvement/tools/p15_table.py`.
- Best gate by the `summary()` rule: smallest gate within one SE of the out-of-sample (oos) minimum, with yield >= 99 % of legacy. "Scan min" is the plain minimum. Raw reference: the full 12 Nov scan `energy_resolution/<run type>/summary.tsv` (gates 2 to 24 ns in 0.5 steps, no 1 ns point; its 8 and 12 ns rows equal P0 `baseline_20261112/pre/`, checked on oos res_cb for all five run types). The calibrated scan covers 1 to 24 ns (H2; Ta to 16 ns), the raw one 2 to 24 ns.
- Pool continuity: the legacy row reproduces baseline `pre/` exactly (res_cb 3.389 / 1.876 / 2.308 / 1.585 / 1.559; oos and n_used also identical). Legacy takes no time coincidence, so the offsets do not enter it.
- Columns: res_cb is sigma/E (%) from the two-peak Crystal Ball, mean over splits, with SE; tail, pile-up flag and yield as in `summary.tsv`; d = calibrated minus legacy (pp); the last two columns compare with the raw scan at the same gate (oos).

## Result

### 0p73GeV_H2

Legacy continuity (P15 vs baseline pre/): res_cb 3.3888 vs 3.3888, res_cb_oos 3.42 vs 3.42, n_used 4.4301e+05 vs 4.4301e+05

| gate (ns) | res_cb in | SE | res_cb oos | SE | tail | pile-up | yield % | d in vs legacy | d oos vs legacy | raw oos (12 Nov) | cal - raw oos |
|---|---|---|---|---|---|---|---|---|---|---|---|
| legacy | 3.389 | 0.014 | 3.420 | 0.014 | 44.77 | 9.78 | 100.00 | +0.000 | +0.000 | 3.420 | +0.000 |
| 1 | 3.717 | 0.031 | 3.768 | 0.034 | 43.51 | 9.65 | 93.39 | +0.328 | +0.348 | nan | +nan |
| 2 | 3.461 | 0.016 | 3.496 | 0.012 | 44.22 | 9.76 | 98.75 | +0.072 | +0.076 | 4.367 | -0.871 |
| 3 | 3.397 | 0.016 | 3.437 | 0.013 | 44.61 | 9.76 | 99.62 | +0.008 | +0.017 | 3.641 | -0.205 |
| 4 | 3.392 | 0.014 | 3.417 | 0.013 | 44.71 | 9.76 | 99.86 | +0.003 | -0.003 | 3.543 | -0.126 |
| 5 | 3.384 | 0.015 | 3.414 | 0.013 | 44.80 | 9.76 | 99.94 | -0.005 | -0.006 | 3.521 | -0.108 |
| 6 | 3.385 | 0.015 | 3.415 | 0.014 | 44.81 | 9.76 | 99.97 | -0.004 | -0.005 | 3.459 | -0.044 |
| 8 | 3.387 | 0.014 | 3.419 | 0.014 | 44.77 | 9.77 | 100.03 | -0.002 | -0.001 | 3.428 | -0.009 |
| 10 | 3.377 | 0.015 | 3.414 | 0.013 | 44.82 | 9.77 | 100.11 | -0.012 | -0.006 | 3.413 | +0.000 |
| 12 | 3.373 | 0.012 | 3.406 | 0.011 | 44.83 | 9.75 | 100.22 | -0.016 | -0.014 | 3.401 | +0.005 |
| 16 | 3.364 | 0.015 | 3.397 | 0.013 | 44.82 | 9.73 | 100.45 | -0.025 | -0.023 | 3.400 | -0.003 |
| 20 | 3.365 | 0.013 | 3.397 | 0.014 | 44.77 | 9.71 | 100.64 | -0.024 | -0.023 | 3.396 | +0.000 |
| 24 | 3.363 | 0.013 | 3.398 | 0.014 | 44.75 | 9.70 | 100.80 | -0.026 | -0.022 | 3.397 | +0.001 |

### 2p24GeV_H2

Legacy continuity (P15 vs baseline pre/): res_cb 1.8759 vs 1.8759, res_cb_oos 1.8888 vs 1.8888, n_used 1.2236e+06 vs 1.2236e+06

| gate (ns) | res_cb in | SE | res_cb oos | SE | tail | pile-up | yield % | d in vs legacy | d oos vs legacy | raw oos (12 Nov) | cal - raw oos |
|---|---|---|---|---|---|---|---|---|---|---|---|
| legacy | 1.876 | 0.005 | 1.889 | 0.005 | 28.17 | 24.29 | 100.00 | +0.000 | +0.000 | 1.889 | +0.000 |
| 1 | 2.369 | 0.006 | 2.379 | 0.006 | 28.58 | 23.75 | 78.21 | +0.493 | +0.490 | nan | +nan |
| 2 | 2.023 | 0.005 | 2.033 | 0.006 | 27.92 | 24.09 | 95.21 | +0.147 | +0.144 | 2.751 | -0.718 |
| 3 | 1.941 | 0.005 | 1.953 | 0.005 | 27.90 | 24.17 | 97.89 | +0.065 | +0.064 | 2.402 | -0.448 |
| 4 | 1.912 | 0.005 | 1.921 | 0.005 | 27.95 | 24.21 | 98.79 | +0.036 | +0.032 | 2.146 | -0.225 |
| 5 | 1.900 | 0.005 | 1.909 | 0.005 | 27.97 | 24.22 | 99.24 | +0.024 | +0.021 | 2.018 | -0.109 |
| 6 | 1.889 | 0.005 | 1.898 | 0.005 | 28.03 | 24.22 | 99.54 | +0.013 | +0.009 | 1.955 | -0.057 |
| 8 | 1.876 | 0.005 | 1.887 | 0.005 | 28.05 | 24.23 | 99.94 | +0.000 | -0.001 | 1.895 | -0.007 |
| 10 | 1.866 | 0.004 | 1.878 | 0.005 | 28.04 | 24.23 | 100.19 | -0.010 | -0.011 | 1.886 | -0.008 |
| 12 | 1.862 | 0.005 | 1.873 | 0.005 | 27.98 | 24.22 | 100.36 | -0.014 | -0.016 | 1.874 | -0.001 |
| 16 | 1.856 | 0.005 | 1.867 | 0.005 | 27.83 | 24.20 | 100.62 | -0.020 | -0.021 | 1.868 | -0.000 |
| 20 | 1.856 | 0.005 | 1.870 | 0.005 | 27.69 | 24.18 | 100.80 | -0.020 | -0.019 | 1.869 | +0.001 |
| 24 | 1.865 | 0.005 | 1.876 | 0.004 | 27.58 | 24.17 | 100.90 | -0.011 | -0.013 | 1.875 | +0.000 |

### 2p24GeV_Ta_X17

Legacy continuity (P15 vs baseline pre/): res_cb 2.3081 vs 2.3081, res_cb_oos 2.332 vs 2.332, n_used 1.7585e+05 vs 1.7585e+05

| gate (ns) | res_cb in | SE | res_cb oos | SE | tail | pile-up | yield % | d in vs legacy | d oos vs legacy | raw oos (12 Nov) | cal - raw oos |
|---|---|---|---|---|---|---|---|---|---|---|---|
| legacy | 2.308 | 0.035 | 2.332 | 0.036 | 49.11 | 36.45 | 100.00 | +0.000 | +0.000 | 2.332 | +0.000 |
| 1 | 2.870 | 0.154 | 2.880 | 0.139 | 41.67 | 22.35 | 75.71 | +0.562 | +0.548 | nan | +nan |
| 2 | 2.232 | 0.048 | 2.270 | 0.039 | 44.32 | 31.36 | 97.64 | -0.076 | -0.062 | 2.770 | -0.501 |
| 3 | 2.040 | 0.026 | 2.049 | 0.025 | 45.38 | 33.73 | 102.50 | -0.269 | -0.283 | 2.409 | -0.360 |
| 4 | 2.023 | 0.024 | 2.030 | 0.023 | 45.80 | 34.24 | 103.04 | -0.285 | -0.302 | 2.336 | -0.306 |
| 5 | 1.992 | 0.030 | 2.027 | 0.026 | 46.13 | 34.51 | 103.15 | -0.316 | -0.305 | 2.187 | -0.161 |
| 6 | 2.008 | 0.023 | 2.022 | 0.025 | 46.17 | 34.70 | 103.27 | -0.300 | -0.310 | 2.097 | -0.076 |
| 8 | 2.004 | 0.024 | 2.019 | 0.027 | 46.39 | 34.90 | 103.31 | -0.304 | -0.313 | 2.054 | -0.035 |
| 10 | 2.016 | 0.024 | 2.026 | 0.027 | 46.44 | 34.96 | 103.25 | -0.292 | -0.306 | 2.041 | -0.015 |
| 12 | 2.020 | 0.022 | 2.031 | 0.027 | 46.56 | 34.97 | 103.20 | -0.288 | -0.301 | 2.050 | -0.019 |
| 16 | 2.047 | 0.021 | 2.062 | 0.025 | 46.67 | 34.86 | 102.98 | -0.261 | -0.270 | 2.066 | -0.004 |

### 3p49GeV_H2

Legacy continuity (P15 vs baseline pre/): res_cb 1.5854 vs 1.5854, res_cb_oos 1.6136 vs 1.6136, n_used 3.8156e+05 vs 3.8156e+05

| gate (ns) | res_cb in | SE | res_cb oos | SE | tail | pile-up | yield % | d in vs legacy | d oos vs legacy | raw oos (12 Nov) | cal - raw oos |
|---|---|---|---|---|---|---|---|---|---|---|---|
| legacy | 1.585 | 0.006 | 1.614 | 0.006 | 11.31 | 33.74 | 100.00 | +0.000 | +0.000 | 1.614 | +0.000 |
| 1 | 2.071 | 0.006 | 2.100 | 0.008 | 9.92 | 32.89 | 70.47 | +0.486 | +0.487 | nan | +nan |
| 2 | 1.762 | 0.009 | 1.788 | 0.011 | 10.59 | 33.46 | 92.40 | +0.176 | +0.174 | 2.592 | -0.804 |
| 3 | 1.671 | 0.007 | 1.701 | 0.008 | 10.85 | 33.62 | 96.44 | +0.086 | +0.087 | 2.191 | -0.490 |
| 4 | 1.636 | 0.006 | 1.670 | 0.007 | 10.97 | 33.66 | 97.83 | +0.051 | +0.056 | 1.873 | -0.203 |
| 5 | 1.629 | 0.008 | 1.660 | 0.008 | 10.99 | 33.69 | 98.51 | +0.044 | +0.046 | 1.784 | -0.124 |
| 6 | 1.612 | 0.008 | 1.641 | 0.008 | 11.12 | 33.70 | 99.00 | +0.027 | +0.028 | 1.705 | -0.064 |
| 8 | 1.589 | 0.007 | 1.616 | 0.006 | 11.23 | 33.70 | 99.72 | +0.004 | +0.002 | 1.637 | -0.021 |
| 10 | 1.577 | 0.006 | 1.604 | 0.007 | 11.24 | 33.70 | 100.10 | -0.009 | -0.009 | 1.610 | -0.006 |
| 12 | 1.571 | 0.007 | 1.600 | 0.007 | 11.21 | 33.71 | 100.33 | -0.015 | -0.014 | 1.602 | -0.002 |
| 16 | 1.565 | 0.006 | 1.595 | 0.007 | 11.13 | 33.70 | 100.60 | -0.020 | -0.019 | 1.595 | +0.000 |
| 20 | 1.569 | 0.006 | 1.597 | 0.007 | 10.98 | 33.68 | 100.77 | -0.017 | -0.017 | 1.597 | +0.000 |
| 24 | 1.574 | 0.006 | 1.602 | 0.006 | 10.84 | 33.67 | 100.88 | -0.012 | -0.012 | 1.601 | +0.002 |

### 3p49GeV_H2_old

Legacy continuity (P15 vs baseline pre/): res_cb 1.5591 vs 1.5591, res_cb_oos 1.579 vs 1.579, n_used 5.7891e+05 vs 5.7891e+05

| gate (ns) | res_cb in | SE | res_cb oos | SE | tail | pile-up | yield % | d in vs legacy | d oos vs legacy | raw oos (12 Nov) | cal - raw oos |
|---|---|---|---|---|---|---|---|---|---|---|---|
| legacy | 1.559 | 0.006 | 1.579 | 0.006 | 10.92 | 25.37 | 100.00 | +0.000 | +0.000 | 1.579 | +0.000 |
| 1 | 2.051 | 0.011 | 2.067 | 0.010 | 9.51 | 24.53 | 70.55 | +0.492 | +0.488 | nan | +nan |
| 2 | 1.734 | 0.007 | 1.752 | 0.008 | 10.23 | 25.09 | 92.39 | +0.175 | +0.173 | 2.639 | -0.887 |
| 3 | 1.651 | 0.006 | 1.668 | 0.006 | 10.43 | 25.24 | 96.44 | +0.092 | +0.089 | 2.175 | -0.507 |
| 4 | 1.622 | 0.006 | 1.638 | 0.006 | 10.51 | 25.31 | 97.82 | +0.063 | +0.059 | 1.845 | -0.207 |
| 5 | 1.607 | 0.005 | 1.623 | 0.007 | 10.60 | 25.33 | 98.51 | +0.048 | +0.044 | 1.747 | -0.125 |
| 6 | 1.600 | 0.006 | 1.615 | 0.007 | 10.65 | 25.36 | 99.00 | +0.041 | +0.036 | 1.676 | -0.061 |
| 8 | 1.566 | 0.006 | 1.584 | 0.007 | 10.83 | 25.35 | 99.69 | +0.007 | +0.005 | 1.602 | -0.018 |
| 10 | 1.555 | 0.006 | 1.573 | 0.006 | 10.88 | 25.36 | 100.06 | -0.004 | -0.006 | 1.584 | -0.011 |
| 12 | 1.547 | 0.006 | 1.566 | 0.006 | 10.89 | 25.36 | 100.25 | -0.012 | -0.013 | 1.570 | -0.004 |
| 16 | 1.546 | 0.006 | 1.564 | 0.006 | 10.80 | 25.35 | 100.47 | -0.013 | -0.015 | 1.564 | +0.000 |
| 20 | 1.547 | 0.006 | 1.565 | 0.006 | 10.69 | 25.34 | 100.61 | -0.012 | -0.014 | 1.564 | +0.002 |
| 24 | 1.550 | 0.006 | 1.568 | 0.006 | 10.60 | 25.34 | 100.72 | -0.009 | -0.011 | 1.569 | -0.001 |

### Best gate per run type (out-of-sample res_cb, %)

| run type | cal rule-best (ns) | res_cb oos | SE | cal scan min (ns) | res_cb oos | raw rule-best (ns) | res_cb oos | raw scan min (ns) | res_cb oos | legacy oos | cal best - legacy | cal best - raw best |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0p73GeV_H2 | 12 | 3.406 | 0.011 | 20 | 3.397 | 11.5 | 3.404 | 20.5 | 3.394 | 3.420 | -0.014 | +0.002 |
| 2p24GeV_H2 | 16 | 1.867 | 0.005 | 16 | 1.867 | 13 | 1.871 | 15.5 | 1.867 | 1.889 | -0.021 | -0.004 |
| 2p24GeV_Ta_X17 | 4 | 2.030 | 0.023 | 8 | 2.019 | 7.5 | 2.061 | 10.5 | 2.040 | 2.332 | -0.302 | -0.031 |
| 3p49GeV_H2 | 12 | 1.600 | 0.007 | 16 | 1.595 | 13 | 1.601 | 16 | 1.595 | 1.614 | -0.014 | -0.001 |
| 3p49GeV_H2_old | 12 | 1.566 | 0.006 | 16 | 1.564 | 12.5 | 1.568 | 18 | 1.562 | 1.579 | -0.013 | -0.002 |

## What the calibrated gate buys (out-of-sample res_cb, pp; negative = better)

| run type | (i) best calibrated gate vs legacy | (ii) best calibrated gate vs best raw gate | gate that matches legacy: calibrated vs raw |
|---|---|---|---|
| 0p73GeV_H2 | -0.014 (12 ns; SE 0.011) | +0.002 (same within SE) | about 4 ns vs about 8 ns |
| 2p24GeV_H2 | -0.021 (16 ns; SE 0.005) | -0.004 (same within SE) | about 8 ns vs about 8 to 10 ns |
| 2p24GeV_Ta_X17 | -0.302 (4 ns) | -0.031 (about 1 SE) | already better than legacy at 2 ns (-0.06) vs about 4 ns |
| 3p49GeV_H2 | -0.014 (12 ns; SE 0.007) | -0.001 (same) | about 8 ns vs about 8 to 10 ns |
| 3p49GeV_H2_old | -0.013 (12 ns; SE 0.006) | -0.002 (same) | about 8 ns vs about 8 to 10 ns |

- At equal gate the calibrated times are far better below 8 ns (for example 2.24 GeV at 4 ns: 1.921 against raw 2.146; 3.49 GeV at 2 ns: 1.788 against 2.592). The offset spread was the steep part of the old scan, as the research note inferred.
- Above about 8 to 10 ns the calibrated and raw curves coincide (differences at or below 0.01 pp at 10 ns and beyond), so the offsets do not move the plateau.
- The best calibrated gate is no better than the best raw gate at the four H2 run types: the difference is below one SE everywhere, and the sign flips. The calibrated gate lets a tighter gate (4 to 8 ns) reach legacy resolution, which matters only if a tight gate is wanted for another reason (pile-up or tail rejection). At the H2 run types the pile-up flag and tail change by at most about 0.1 and 0.4 pp between gates of 8 ns and beyond, and tight calibrated gates show no gain in either (the Ta tail falls 49 to 46 % with any gate of 3 ns or more).
- Ta X17 is the exception: the calibrated gate buys 0.30 pp over legacy and 0.03 pp (about 1 SE) over the best raw gate, with a flat minimum from 3 to 12 ns (2.02 to 2.05 %). Yield exceeds legacy by 3 % at 3 to 16 ns because the gated path keeps secondary pulses legacy drops.
- With 20 and 24 ns added, the H2 scan minima lie inside the range: 16 ns at 2.24 and 3.49 GeV (and 3.49 old), 20 ns at 0.73 GeV (3.397, against 3.406 at 12 ns). The curve is flat from 12 ns to 24 ns within about 0.01 pp and rises slightly at 24 ns at 2.24 GeV (1.876). The `summary()` rule returns 12 ns at 0.73, 3.49 and 3.49 old GeV and 16 ns at 2.24 GeV.
- Seed-member time difference after correction: not recomputed here. The research note's recount on the crosstalk cache gives a calibrated robust sigma of 0.7, 1.05 and 1.4 ns at 0.73, 2.24 and 3.49 GeV (raw 2.3, 2.8 and 3.1 ns).

## Caveats

- The offset table was derived on WaveAnalyzer times (neighbour-pair network, `hycal_module_time_calib.cpp`); applying it to matched-filter times is a cross-basis use. Self-referenced MF offsets agree on pairs within RMS 0.41 to 0.46 ns.
- The table has no lead-glass rows, so lead glass keeps offset 0 and PbWO4 shifts by about +1.4 ns relative to it. Transition-ring (T) clusters and the lead-glass region rows are not comparable between the raw and calibrated scans. The pool (inner PbWO4) is unaffected.
- Ta X17 coverage rests on inference: no raw-versus-calibrated Ta width was recounted (calibrated core |dt| < 2 ns holds 69 to 83 % of hits there against 78 to 95 % for H2).
- No time-walk correction (out of scope); small late members remain untested.
- `energy_resolution_benchmark.py` now differs from the sha in `baseline_20261112/p0_code.sha1`; with `TIME_OFFSETS` off the code path is unchanged.

## Files

`improvement/P15_tgate_calibrated.tsv` (every run type and gate); `output/energy_resolution/variants/P15/{<run type>,logs,run.sh,run.out}`; `improvement/tools/p15_table.py`. Overrides are recorded in each `regions_splits.tsv` (`# overrides: TIME_OFFSETS=True`).

## Figure

`figures/P15_seed_member_dt.png` (and `.pdf`), made by `tools/p15_dt_hist.py`: histogram of dt = t_seed - t_member (0.5 ns bins, log y) for non-seed members of the legacy (ungated) clusters, raw and offset-corrected (calibrated = raw - offset_ns, resolved rows), with the P15 best calibrated gate marked. Source: `crosstalk_study/cache/<run type>/*_legacy.npz`, first 3 splits per cache dir (Ta: the three run dirs pooled), no e-p cut. Robust sigma = half the 16-84 % range.

| run type | members | robust sigma raw (ns) | robust sigma calibrated (ns) | best gate (ns) |
|---|---|---|---|---|
| 0.73 GeV H2 | 2.53 M | 2.32 | 0.70 | 12 |
| 2.24 GeV H2 | 4.12 M | 2.81 | 1.06 | 16 |
| 2.24 GeV Ta X17 | 0.67 M | 3.20 | 1.43 | 4 |
| 3.49 GeV H2 | 2.70 M | 3.13 | 1.45 | 12 |

The late tail (member later than seed, negative dt) is the larger one and the offsets do not remove it; the flat floor beyond about 10 ns is unchanged by the calibration.
