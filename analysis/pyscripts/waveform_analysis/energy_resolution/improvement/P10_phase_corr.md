# P10: onset-phase correction I(phi) (D5 multi-run verdict + residual-phase check)

Control v2 (`control_phase2_v2_20261113/<RT>/`), pool ep, control gate path (16 ns; 0.73 GeV 20 ns); Delta = variant - control, paired per (run, split), mean +- SE (D3); res_cb in pp, leak_cut / eff_cut in %. Full per-run table with all guards: `P10_phase_corr.tsv` (tool `tools/p10_table.py`; override string and run in every row). `bank_sha1`: 2.24 GeV control scans 15dfd867, rescans 7f125b75 (hits identical, P5 section 6; same-bank rule relaxed as in P5), 3.49 GeV 294184b0 on both arms; refinement reports use the control scans. Code: `energy_resolution_benchmark.py` sha1 843782c7 at the end of the work (concurrent editors); new switches `PHASE_STUDY` (dump) and `PHASE_RESID` (table), both default off. Variants: `output/energy_resolution/variants/P10/` (`phase_off/`, `G2/ M1/ M2/ G2x/ M1x/`, `tables/`, `phase_ctrl/`, scripts, logs); rescans `output/energy_resolution/scan_variants/phase_off/`.

## 1. D5 verdict: PHASE_CORR=False vs control v2

| run type | run | pairs | Delta res_cb | Delta res_cb oos | Delta tail_hi_f | Delta leak_cut | Delta eff_cut | yield |
|---|---|---|---|---|---|---|---|---|
| 2.24 | 25228 | 11 | +0.040 +- 0.007 | +0.043 +- 0.009 | +0.004 +- 0.002 | +0.43 +- 0.03 | -0.08 +- 0.02 | +0.15 % |
| 2.24 | 25257 | 7 | +0.060 +- 0.005 | +0.062 +- 0.007 | +0.004 +- 0.003 | +0.45 +- 0.06 | -0.08 +- 0.03 | +0.31 % |
| 2.24 | 25260 | 8 | +0.067 +- 0.006 | +0.062 +- 0.007 | +0.002 +- 0.002 | +0.47 +- 0.03 | -0.09 +- 0.01 | +0.05 % |
| 2.24 | 25299 | 8 | +0.051 +- 0.006 | +0.058 +- 0.008 | +0.003 +- 0.002 | +0.40 +- 0.03 | -0.08 +- 0.03 | +0.19 % |
| 2.24 | 25320 | 8 | +0.063 +- 0.006 | +0.058 +- 0.006 | +0.008 +- 0.003 | +0.51 +- 0.07 | -0.03 +- 0.03 | +0.19 % |
| 2.24 | pooled | 42 | +0.055 +- 0.003 | +0.055 +- 0.004 | +0.004 +- 0.001 | +0.45 +- 0.02 | -0.07 +- 0.01 | +0.17 % |
| 3.49 | 24698 | 8 | +0.038 +- 0.006 | +0.038 +- 0.005 | +0.039 +- 0.012 | +0.61 +- 0.12 | -0.13 +- 0.02 | +0.26 % |
| 3.49 | 24911 | 8 | +0.062 +- 0.005 | +0.060 +- 0.005 | +0.047 +- 0.008 | +0.95 +- 0.10 | -0.11 +- 0.04 | +0.40 % |
| 3.49 | pooled | 16 | +0.050 +- 0.005 | +0.049 +- 0.004 | +0.043 +- 0.007 | +0.78 +- 0.09 | -0.12 +- 0.02 | +0.33 % |

Verdict: the phase correction is real and worth 0.055 pp (2.24 GeV) and 0.050 pp (3.49 GeV) of res_cb, plus 0.45 / 0.78 % of leak_cut, on every one of the 7 runs with the same sign (13 SE pooled at 2.24; replicated across all runs; Delta mu < 0.01 %). The run 25320 single-run number (+0.063) is within the run spread. Keep `PHASE_CORR=True` (production default; the "-0.04 pp production only" expectation was too small: the correction is a gain in the benchmark itself). Not run: 3.49 old, 0.73 GeV, Ta (no candidate; the scans of those types were not rescanned).

## 2. Residual phase dependence after the correction

`PHASE_STUDY` dumps, for the control-gate pool clusters (control v2 code and defaults), E/E' and for the three largest members the energy fraction f_j, module and onset phase x_j (table position 0-256, same formula as `scan`). Model (clusters with 0.95 < E/E' < 1.05): E/E' = 1 + sum_j f_j delta(x_j). First-harmonic delta (per unit f, permille):

| sample | cos | sin | amplitude |
|---|---|---|---|
| 2.24 GeV, correction on (5 runs; per run -3.9..-4.2 / +1.5..+2.0, SE 0.1) | -4.0 | +1.9 | 4.5 |
| 3.49 GeV, on (24698 / 24911) | -3.4 / -3.6 | +1.1 / +1.2 | 3.6 |
| run 25320, correction off (rescan) | -8.8 | -2.9 | 9.3 |
| mean applied correction (phase_table, all modules) | +2.2 | +4.4 | 4.9 |

The correction removes only about half of the real dependence: it is about 0.75 of the needed size and rotated by about 25 degrees. Binned median E/E' of the seed phase spans -5.3 .. +3.0 permille (16 bins) after the correction. The residual is proportional to the seed fraction, equal in all five 2.24 GeV runs and in both 3.49 GeV runs (stable), and module-specific: per-module amplitude scatter about 2.4 permille per component against a median SE of 0.7 (chi2/dof 27, 125 modules with >= 800 clusters; by quadrant the cos term runs -2.2 .. -4.9). Cross-checked on held-out clusters: the window std of E/E' falls 18.295 -> 18.166 (global k=1), 18.139 (k=1,2), 18.059 (per-module k=1), 18.033 permille (per-module k=1,2).

## 3. Refined forms under D3

Residual factor applied to each PbWO4 matched-filter hit energy, E -> E / (1 + delta(x)), at report level (`PHASE_RESID`); the table is fitted on one half of each run's splits (ordinal parity) and applied to the other, so the reports are cross-fitted (gains, position correction and shape factor refit in the report as usual). G2: global k=1,2 per run type (4 numbers). M1 / M2: G2 + per-module deviations, ridge lam 100 (module k=1 / k=1,2; about 500 module rows; lam chosen by held-out std, flat from 10 to 300). Pooled Delta vs control v2:

| form | run type | pairs | Delta res_cb | Delta res_cb oos | Delta tail_hi_f | Delta tail_lo_f | Delta leak_cut | Delta eff_cut | yield |
|---|---|---|---|---|---|---|---|---|---|
| G2 | 2.24 | 42 | -0.021 +- 0.002 | -0.022 +- 0.002 | -0.001 +- 0.001 | +0.007 +- 0.006 | -0.12 +- 0.01 | +0.07 +- 0.01 | -0.12 % |
| G2 | 3.49 | 16 | -0.021 +- 0.003 | -0.019 +- 0.003 | -0.013 +- 0.004 | -0.023 +- 0.009 | -0.22 +- 0.06 | +0.12 +- 0.02 | -0.16 % |
| G2 | 3.49 old (own table) | 28 | -0.015 +- 0.003 | -0.022 +- 0.002 | -0.007 +- 0.004 | -0.051 +- 0.013 | -0.10 +- 0.06 | +0.13 +- 0.02 | -0.15 % |
| G2x | 3.49 old (3.49 table, other runs) | 28 | -0.016 +- 0.003 | -0.022 +- 0.003 | -0.005 +- 0.004 | -0.045 +- 0.014 | -0.12 +- 0.07 | +0.13 +- 0.02 | -0.15 % |
| G2 | 0.73 | 16 | -0.009 +- 0.008 | -0.025 +- 0.010 | 0 | -0.022 +- 0.023 | 0.00 | +0.03 +- 0.03 | +0.89 % |
| M1 | 2.24 | 42 | -0.028 +- 0.003 | -0.031 +- 0.002 | -0.001 +- 0.001 | +0.014 +- 0.007 | -0.21 +- 0.01 | +0.09 +- 0.01 | -0.16 % |
| M1 | 3.49 | 16 | -0.040 +- 0.004 | -0.037 +- 0.003 | -0.024 +- 0.004 | -0.051 +- 0.012 | -0.44 +- 0.07 | +0.18 +- 0.02 | -0.30 % |
| M1 | 3.49 old (own table) | 28 | -0.034 +- 0.004 | -0.037 +- 0.004 | -0.021 +- 0.004 | -0.071 +- 0.014 | -0.41 +- 0.08 | +0.19 +- 0.02 | -0.23 % |
| M1x | 3.49 old (3.49 table, other runs) | 28 | -0.035 +- 0.004 | -0.037 +- 0.003 | -0.019 +- 0.004 | -0.088 +- 0.014 | -0.35 +- 0.08 | +0.19 +- 0.02 | -0.26 % |
| M1 | 0.73 | 16 | -0.019 +- 0.007 | -0.021 +- 0.009 | 0 | -0.048 +- 0.026 | 0.00 | -0.005 +- 0.023 | +0.61 % |
| M2 | 2.24 | 42 | -0.033 +- 0.003 | -0.036 +- 0.002 | -0.001 +- 0.001 | +0.014 +- 0.006 | -0.22 +- 0.02 | +0.10 +- 0.01 | -0.21 % |
| M2 | 3.49 | 16 | -0.042 +- 0.003 | -0.039 +- 0.003 | -0.030 +- 0.005 | -0.077 +- 0.012 | -0.45 +- 0.06 | +0.21 +- 0.01 | -0.21 % |

Reading. Every form moves res_cb in the right direction in all runs with every guard favourable (tail_hi_f lower or null, eff_cut up, yield loss <= 0.3 %, leak_cut lower); the low tail is null at 2.24 (+0.014 +- 0.007 at 1.9 SE, tolerated by the user's tail_lo_f rule only if it is the fixed-offset tail, which it is) and lower at 3.49. M1 meets the D3 adoption rule (Delta <= -2 SE at 3.49 with the same sign on 3.49 old, and at 2.24); 0.73 GeV res_cb in-sample is null for G2 and -0.019 +- 0.007 for M1, yield +0.6 / +0.9 %. M1x (table fitted on 24698 + 24911, applied to the earlier 3.49 GeV runs) reproduces the own-table result, so the per-module table is a stable detector property, not a fit to the evaluated clusters. M2 adds only -0.005 / -0.002 pp over M1 (plateau). The size is comparable with the P9 shape factor. Caveats: one fit per run type; 2.24 vs 3.49 amplitudes differ by about 15 % (G2 fitted per energy; not tested across energies); the table must be refitted whenever the bank or the t0 definition changes (it assumes the bank's t50 of the scan); no Ta test.

## Recommendation

1. Keep PHASE_CORR as is (verdict, section 1; the switch stays).
2. Refine: adopt form M1 (G2 global k=1,2 plus per-module k=1 residual, lam 100) as a calibration table next to the gains, gain about -0.03 pp (2.24) / -0.04 pp (3.49, 3.49 old) in res_cb and -0.2 / -0.4 % in leak_cut, cross-fitted; G2 (four numbers per energy) takes 60 % of it at 2.24 and 50 % at 3.49 without a per-module table. A production implementation needs a fitter (tools/p10_phase_fit.py is the prototype) and a decision on where the table lives (U3-type default, off until the user confirms). 0.73 GeV and Ta are not covered.
3. Not done: refit across energies, Ta, 24937; the root cause (why the template-derived I(phi) is 0.75 x and rotated; likely matched-filter amplitude / t0 coupling beyond the integral estimate) is not investigated.

Review: critic (code, light, PASS); SHOULD C1 (bank_sha1 check of t50 at phase_study / phase_resid; hits unchanged by the 2.24 rebuild, P5) and C2 (the fit window cuts on E/E'; mild bias toward zero amplitudes, so the table is conservative; wider window not rerun) noted, not applied.
