# P9: edge-aware position correction

Control: `control_phase2_20261113` defaults (TIME_OFFSETS, control gate 20/16/16/16 ns, R7ri, 24937 excluded), re-run with `--set EVAL_REGIONS=True` into `output/energy_resolution/variants/P9/control/` (pool rows identical to the phase-2 control, max diff 0). Variants: `output/energy_resolution/variants/P9/<variant>/<RT>/` (each with `delta.tsv` from `compare_reports.py`, `poscorr_extra.tsv`, logs in `logs/`). Full per-variant table: `P9_edge_poscorr.tsv` (D12; tool `tools/p9_table.py`, runner `tools/p9_run.sh`). Override string: `--set EVAL_REGIONS=True --set 'PC_EXTRA="<variant>"'`; `bank_sha1` as the scans of the control. Values: mean ± SE of paired per-split differences (variant − control), **out of sample** (even/odd cross fit, the extra factor is refitted inside each half); units % (res_cb) and pp (others). Ta not run (it borrows the 2.24 GeV map; the switch is not applied to it, a warning is printed).

## Method

BM switch `PC_EXTRA` (default `""`, off; byte-identical outputs when off): after the existing gains + position correction, an extra multiplicative factor `E / R_extra`, trained per run type on the same pool as the position correction (target y = E/(R E'), window PC_R, iterated 2.5 MAD window, polynomial least squares, scaled so that the median of y does not move). Variables: `big`/`sm` = larger/smaller of |u|,|v| (distance to the module boundary), `f1,f2,f3` = calibrated energy fractions of the three largest members, `nm` = member count, `rad` = seed radius. Edge class = |u| or |v| > 0.35 (the `pool_edge` definition), centre = rest. Variants: `pos_edge` (edge only, position variables only, degree 4, pooled over modules = edge-distance parametrisation), `shape_edge` (edge only, position + shape, degree 3), `shape_all` / `shape_all2` / `shape_all4` (edge and centre fitted separately, degree 3/2/4), `shape_all_ring` (+ rad). `shape_all2`: 28 coefficients per class (56 total) from 1.7e5 to 2.8e5 training clusters per class and half. A gradient-boosted tree on the same variables (offline, random half split) gave no more than the degree-3 polynomial on the shape variables, so capacity is not the limit.

## Region table, shape_all2 (recommended), Δ vs control, out of sample

| run type | variant | region | control res_cb oos | Δ res_cb oos | Δ eff_cut oos | Δ tail_hi3 oos | Δ tail_hi_f oos | Δ tail_lo_f oos | Δ leak_cut oos | Δ res_cb in | Δ gap | yield loss % |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0p73GeV_H2 | shape_all2 | pool | 3.397 | -0.005 ± 0.013 | +0.839 ± 0.052 | +0.029 ± 0.004 | -0.000 ± 0.001 | -1.855 ± 0.070 | +0.003 ± 0.001 | +0.002 ± 0.011 | -0.006 ± 0.010 | 2.02 |
| 0p73GeV_H2 | shape_all2 | pool_edge | 3.441 | -0.016 ± 0.025 | +0.786 ± 0.085 | +0.029 ± 0.007 | -0.000 ± 0.001 | -1.627 ± 0.084 | +0.004 ± 0.002 | -0.038 ± 0.029 | +0.022 ± 0.015 | 2.02 |
| 0p73GeV_H2 | shape_all2 | pool_centre | 3.309 | +0.024 ± 0.023 | +0.862 ± 0.064 | +0.032 ± 0.011 | +0.001 ± 0.001 | -2.130 ± 0.085 | +0.003 ± 0.001 | +0.060 ± 0.022 | -0.036 ± 0.016 | 2.02 |
| 2p24GeV_H2 | shape_all2 | pool | 1.858 | -0.057 ± 0.003 | +0.343 ± 0.015 | +0.007 ± 0.003 | -0.005 ± 0.002 | -0.047 ± 0.009 | -0.391 ± 0.022 | -0.057 ± 0.003 | +0.000 ± 0.002 | -0.07 |
| 2p24GeV_H2 | shape_all2 | pool_edge | 1.902 | -0.051 ± 0.006 | +0.482 ± 0.023 | +0.010 ± 0.004 | -0.003 ± 0.002 | -0.115 ± 0.014 | -0.293 ± 0.038 | -0.043 ± 0.005 | -0.008 ± 0.004 | -0.07 |
| 2p24GeV_H2 | shape_all2 | pool_centre | 1.811 | -0.059 ± 0.005 | +0.265 ± 0.019 | -0.000 ± 0.004 | -0.006 ± 0.002 | +0.005 ± 0.011 | -0.447 ± 0.029 | -0.060 ± 0.004 | +0.000 ± 0.002 | -0.07 |
| 3p49GeV_H2 | shape_all2 | pool | 1.561 | -0.035 ± 0.005 | +0.412 ± 0.033 | +0.019 ± 0.006 | -0.004 ± 0.012 | -0.284 ± 0.025 | -0.451 ± 0.104 | -0.037 ± 0.005 | +0.002 ± 0.003 | 0.03 |
| 3p49GeV_H2 | shape_all2 | pool_edge | 1.615 | -0.017 ± 0.006 | +0.651 ± 0.054 | +0.025 ± 0.006 | +0.026 ± 0.014 | -0.500 ± 0.044 | -0.127 ± 0.130 | -0.010 ± 0.006 | -0.007 ± 0.005 | 0.03 |
| 3p49GeV_H2 | shape_all2 | pool_centre | 1.506 | -0.044 ± 0.006 | +0.233 ± 0.031 | +0.005 ± 0.007 | -0.017 ± 0.011 | -0.151 ± 0.024 | -0.716 ± 0.113 | -0.050 ± 0.007 | +0.006 ± 0.005 | 0.03 |
| 3p49GeV_H2_old | shape_all2 | pool | 1.532 | -0.032 ± 0.003 | +0.475 ± 0.026 | +0.015 ± 0.005 | -0.012 ± 0.007 | -0.276 ± 0.016 | -0.353 ± 0.073 | -0.030 ± 0.004 | -0.002 ± 0.003 | -0.05 |
| 3p49GeV_H2_old | shape_all2 | pool_edge | 1.590 | -0.011 ± 0.007 | +0.724 ± 0.043 | +0.019 ± 0.009 | +0.025 ± 0.009 | -0.475 ± 0.034 | -0.014 ± 0.104 | -0.020 ± 0.005 | +0.009 ± 0.005 | -0.05 |
| 3p49GeV_H2_old | shape_all2 | pool_centre | 1.467 | -0.043 ± 0.005 | +0.310 ± 0.027 | -0.003 ± 0.007 | -0.043 ± 0.010 | -0.133 ± 0.019 | -0.540 ± 0.096 | -0.037 ± 0.004 | -0.006 ± 0.006 | -0.05 |

Control values (oos):

| run type | region | res_cb oos | eff_cut | tail_hi3 | tail_hi_f | tail_lo_f | leak_cut |
|---|---|---|---|---|---|---|---|
| 0p73GeV_H2 | pool | 3.397 | 49.370 | 0.041 | 0.002 | 21.033 | 0.005 |
| 0p73GeV_H2 | pool_edge | 3.441 | 47.561 | 0.043 | 0.003 | 22.043 | 0.007 |
| 0p73GeV_H2 | pool_centre | 3.309 | 51.343 | 0.045 | 0.001 | 19.869 | 0.004 |
| 2p24GeV_H2 | pool | 1.858 | 69.108 | 0.125 | 0.075 | 26.928 | 4.557 |
| 2p24GeV_H2 | pool_edge | 1.902 | 66.328 | 0.078 | 0.035 | 29.099 | 4.660 |
| 2p24GeV_H2 | pool_centre | 1.811 | 71.390 | 0.159 | 0.107 | 25.106 | 4.349 |
| 3p49GeV_H2 | pool | 1.561 | 84.626 | 0.113 | 0.277 | 12.176 | 22.080 |
| 3p49GeV_H2 | pool_edge | 1.615 | 81.867 | 0.084 | 0.299 | 14.451 | 22.007 |
| 3p49GeV_H2 | pool_centre | 1.506 | 86.858 | 0.128 | 0.217 | 10.340 | 21.891 |
| 3p49GeV_H2_old | pool | 1.532 | 85.052 | 0.103 | 0.238 | 11.703 | 21.697 |
| 3p49GeV_H2_old | pool_edge | 1.590 | 82.518 | 0.087 | 0.277 | 13.739 | 21.794 |
| 3p49GeV_H2_old | pool_centre | 1.467 | 87.036 | 0.093 | 0.166 | 10.085 | 21.150 |

## Variant family at 3.49 GeV (pool and edge)

| run type | variant | region | control res_cb oos | Δ res_cb oos | Δ eff_cut oos | Δ tail_hi3 oos | Δ tail_hi_f oos | Δ tail_lo_f oos | Δ leak_cut oos | Δ res_cb in | Δ gap | yield loss % |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 3p49GeV_H2 | shape_edge | pool | 1.561 | -0.011 ± 0.003 | +0.261 ± 0.020 | +0.020 ± 0.005 | +0.009 ± 0.009 | -0.199 ± 0.015 | -0.085 ± 0.066 | -0.012 ± 0.002 | +0.001 ± 0.002 | -0.07 |
| 3p49GeV_H2 | shape_edge | pool_edge | 1.615 | -0.016 ± 0.006 | +0.662 ± 0.053 | +0.041 ± 0.007 | +0.053 ± 0.015 | -0.512 ± 0.038 | -0.039 ± 0.134 | -0.016 ± 0.005 | -0.000 ± 0.005 | -0.07 |
| 3p49GeV_H2 | shape_all | pool | 1.561 | -0.036 ± 0.005 | +0.418 ± 0.033 | +0.027 ± 0.008 | +0.007 ± 0.012 | -0.292 ± 0.026 | -0.422 ± 0.110 | -0.041 ± 0.004 | +0.005 ± 0.002 | 0.06 |
| 3p49GeV_H2 | shape_all | pool_edge | 1.615 | -0.016 ± 0.006 | +0.672 ± 0.053 | +0.040 ± 0.007 | +0.051 ± 0.015 | -0.520 ± 0.039 | -0.030 ± 0.133 | -0.017 ± 0.005 | +0.000 ± 0.005 | 0.06 |
| 3p49GeV_H2 | shape_all2 | pool | 1.561 | -0.035 ± 0.005 | +0.412 ± 0.033 | +0.019 ± 0.006 | -0.004 ± 0.012 | -0.284 ± 0.025 | -0.451 ± 0.104 | -0.037 ± 0.005 | +0.002 ± 0.003 | 0.03 |
| 3p49GeV_H2 | shape_all2 | pool_edge | 1.615 | -0.017 ± 0.006 | +0.651 ± 0.054 | +0.025 ± 0.006 | +0.026 ± 0.014 | -0.500 ± 0.044 | -0.127 ± 0.130 | -0.010 ± 0.006 | -0.007 ± 0.005 | 0.03 |
| 3p49GeV_H2 | shape_all4 | pool | 1.561 | -0.036 ± 0.005 | +0.398 ± 0.036 | +0.044 ± 0.006 | +0.034 ± 0.011 | -0.277 ± 0.027 | -0.403 ± 0.111 | -0.039 ± 0.005 | +0.003 ± 0.003 | 0.01 |
| 3p49GeV_H2 | shape_all4 | pool_edge | 1.615 | -0.015 ± 0.006 | +0.669 ± 0.049 | +0.058 ± 0.007 | +0.086 ± 0.015 | -0.508 ± 0.043 | +0.024 ± 0.127 | -0.013 ± 0.005 | -0.002 ± 0.005 | 0.01 |
| 3p49GeV_H2 | shape_all_ring | pool | 1.561 | -0.046 ± 0.007 | +0.432 ± 0.040 | +0.036 ± 0.006 | +0.003 ± 0.012 | -0.300 ± 0.032 | -0.598 ± 0.126 | -0.045 ± 0.006 | -0.001 ± 0.003 | 0.05 |
| 3p49GeV_H2 | shape_all_ring | pool_edge | 1.615 | -0.025 ± 0.008 | +0.814 ± 0.053 | +0.054 ± 0.008 | +0.062 ± 0.015 | -0.610 ± 0.042 | -0.107 ± 0.146 | -0.017 ± 0.006 | -0.008 ± 0.007 | 0.05 |
| 3p49GeV_H2 | pos_edge | pool | 1.561 | +0.007 ± 0.002 | +0.105 ± 0.018 | +0.015 ± 0.003 | +0.033 ± 0.004 | -0.088 ± 0.010 | +0.176 ± 0.032 | +0.005 ± 0.001 | +0.002 ± 0.001 | -0.07 |
| 3p49GeV_H2 | pos_edge | pool_edge | 1.615 | +0.017 ± 0.004 | +0.245 ± 0.031 | +0.023 ± 0.006 | +0.073 ± 0.008 | -0.203 ± 0.017 | +0.373 ± 0.076 | +0.018 ± 0.003 | -0.001 ± 0.004 | -0.07 |

## Edge-only variant (the literal P9 target)

| run type | variant | region | control res_cb oos | Δ res_cb oos | Δ eff_cut oos | Δ tail_hi3 oos | Δ tail_hi_f oos | Δ tail_lo_f oos | Δ leak_cut oos | Δ res_cb in | Δ gap | yield loss % |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0p73GeV_H2 | shape_edge | pool | 3.397 | +0.017 ± 0.009 | +0.356 ± 0.049 | +0.015 ± 0.003 | +0.000 ± 0.000 | -0.814 ± 0.053 | +0.002 ± 0.001 | +0.016 ± 0.012 | +0.002 ± 0.006 | 1.05 |
| 0p73GeV_H2 | shape_edge | pool_edge | 3.441 | +0.008 ± 0.025 | +0.760 ± 0.071 | +0.040 ± 0.006 | +0.000 ± 0.001 | -1.782 ± 0.089 | +0.008 ± 0.003 | -0.019 ± 0.028 | +0.027 ± 0.020 | 1.05 |
| 0p73GeV_H2 | shape_edge | pool_centre | 3.309 | -0.006 ± 0.005 | +0.036 ± 0.050 | +0.000 ± 0.000 | +0.000 ± 0.000 | +0.116 ± 0.057 | -0.000 ± 0.001 | +0.005 ± 0.010 | -0.011 ± 0.008 | 1.05 |
| 2p24GeV_H2 | shape_edge | pool | 1.858 | -0.025 ± 0.003 | +0.185 ± 0.010 | +0.004 ± 0.002 | -0.002 ± 0.001 | -0.035 ± 0.007 | -0.165 ± 0.017 | -0.025 ± 0.002 | +0.000 ± 0.002 | -0.11 |
| 2p24GeV_H2 | shape_edge | pool_edge | 1.902 | -0.049 ± 0.006 | +0.475 ± 0.021 | +0.014 ± 0.005 | -0.001 ± 0.002 | -0.098 ± 0.015 | -0.271 ± 0.039 | -0.039 ± 0.005 | -0.010 ± 0.004 | -0.11 |
| 2p24GeV_H2 | shape_edge | pool_centre | 1.811 | +0.000 ± 0.000 | +0.000 ± 0.004 | +0.000 ± 0.000 | +0.000 ± 0.000 | +0.002 ± 0.004 | +0.000 ± 0.002 | -0.000 ± 0.000 | +0.001 ± 0.000 | -0.11 |
| 3p49GeV_H2 | shape_edge | pool | 1.561 | -0.011 ± 0.003 | +0.261 ± 0.020 | +0.020 ± 0.005 | +0.009 ± 0.009 | -0.199 ± 0.015 | -0.085 ± 0.066 | -0.012 ± 0.002 | +0.001 ± 0.002 | -0.07 |
| 3p49GeV_H2 | shape_edge | pool_edge | 1.615 | -0.016 ± 0.006 | +0.662 ± 0.053 | +0.041 ± 0.007 | +0.053 ± 0.015 | -0.512 ± 0.038 | -0.039 ± 0.134 | -0.016 ± 0.005 | -0.000 ± 0.005 | -0.07 |
| 3p49GeV_H2 | shape_edge | pool_centre | 1.506 | -0.001 ± 0.000 | -0.004 ± 0.006 | +0.002 ± 0.001 | +0.001 ± 0.001 | +0.003 ± 0.004 | -0.014 ± 0.008 | -0.000 ± 0.000 | -0.000 ± 0.000 | -0.07 |
| 3p49GeV_H2_old | shape_edge | pool | 1.532 | -0.010 ± 0.004 | +0.320 ± 0.012 | +0.015 ± 0.005 | +0.008 ± 0.004 | -0.217 ± 0.011 | -0.085 ± 0.066 | -0.011 ± 0.004 | +0.000 ± 0.003 | -0.17 |
| 3p49GeV_H2_old | shape_edge | pool_edge | 1.590 | -0.013 ± 0.008 | +0.776 ± 0.040 | +0.029 ± 0.009 | +0.041 ± 0.010 | -0.492 ± 0.029 | -0.001 ± 0.108 | -0.024 ± 0.007 | +0.011 ± 0.005 | -0.17 |
| 3p49GeV_H2_old | shape_edge | pool_centre | 1.467 | -0.000 ± 0.001 | +0.008 ± 0.006 | -0.000 ± 0.000 | -0.001 ± 0.001 | -0.006 ± 0.005 | -0.006 ± 0.014 | +0.001 ± 0.001 | -0.001 ± 0.001 | -0.17 |

Edge − centre `res_cb` oos (pp):

| run type | control | shape_edge | shape_all2 |
|---|---|---|---|
| 0p73GeV_H2 | 0.132 | 0.146 | 0.092 |
| 2p24GeV_H2 | 0.090 | 0.042 | 0.099 |
| 3p49GeV_H2 | 0.109 | 0.093 | 0.137 |
| 3p49GeV_H2_old | 0.123 | 0.109 | 0.154 |

## Reading

- Position-only edge form (`pos_edge`): null to slightly worse in the core (+0.008 / +0.007 pp pool at 2.24 / 3.49), because a pooled polynomial replaces nothing the per-module hybrid map does not already fit. The mean response is already flat across the module; the edge excess is fluctuation, not bias.
- Shower-shape factor: the core gain is mostly in the centre, so the edge − centre core gap does not close with `shape_all2` (0.109 → 0.136 at 3.49, 0.090 → 0.100 at 2.24). `shape_edge` narrows it at 2.24 (0.090 → 0.043) and 3.49 (0.109 → 0.092) but has a 0.73 GeV core cost (+0.017 ± 0.009).
- Edge tail: `eff_cut` +0.65 / +0.48 pp and `tail_lo_f` −0.50 / −0.12 pp at the edge (3.49 / 2.24). The edge `eff_cut` deficit against the centre (−5.0 pp at 3.49) shrinks by about 0.4 pp, i.e. under 10 %. The P5 tail clause (target `eff_cut` −4.4 / −5.7 pp) is mostly not closed; what remains is event-by-event loss that these variables do not predict.
- Pool (`shape_all2`): `res_cb` oos −0.057 ± 0.003 (2.24), −0.035 ± 0.005 (3.49), −0.032 ± 0.003 (3.49 old); `leak_cut` −0.39 ± 0.02 / −0.45 ± 0.10 / −0.35 ± 0.07 pp; `eff_cut` +0.34 / +0.41 / +0.48 pp; `tail_lo_f` lower at all energies; U2 gap change +0.000 ± 0.002 / +0.002 ± 0.003 / −0.002 ± 0.003 pp (no overfit). At 0.73 GeV: `res_cb` null (−0.005 ± 0.013), `tail_lo_f` −1.86 pp, `eff_cut` +0.84 pp, `leak_cut` +0.003 pp (negligible).
- `tail_hi3` is up by +0.019 ± 0.006 (3.49), +0.007 ± 0.003 (2.24), +0.015 ± 0.005 (old), +0.029 (0.73): a 3σ_CB line that moves inward as σ narrows. The fixed-offset `tail_hi_f` is −0.004 ± 0.012 / −0.005 ± 0.002 / −0.012 ± 0.007 / 0.000. D3 amended `tail_lo3` → `tail_lo_f` for σ-narrowing variants but did not state the high side; applying `tail_hi_f` there is a ruling needed (open question).
- Guards: yield loss ≤ 0.1 % at 2.24 / 3.49 / old; **2.0 % at 0.73 GeV** (> 1 %): the fitted peak moves −0.44 pp against E' there (−0.06 pp at 3.49, −0.07 at 2.24); the scale would be re-anchored by the gain solve in production, and 0.73 GeV alone never decides.
- First pass used a mean normalisation (`P9_old_norm`, deleted): same widths, peak shift −0.07 pp; replaced by the median normalisation.

## Recommendation (U7)

Adopt `shape_all2` as a switched, default-off option, then default-on once the `tail_hi_f` ruling is made: D3 rule (1) met at 3.49 (−7 SE) with the same sign in the old runs and at 2.24; `leak_cut` and `eff_cut` improve together at 2.24 / 3.49 (separation gain about 2.4 % in σ at 3.49); no overfit by the gap. It is a cluster-shape correction (below-threshold and sharing loss, near P12), not an edge fix: the edge-specific P9 goal (edge tails to the centre values) is not met, and the 0.73 GeV yield guard fails. `shape_edge` alone is smaller (−0.011 pp at 3.49) and costs 0.73 GeV core. Refit after P6 / P16 change the gains (the factor is trained on calibrated member energies); BM assumes member-indexed gains (not `GAIN_MODE seed`).

Log: P9 shape factor `shape_all2` found (core −0.057/−0.035 pp at 2.24/3.49, leak_cut −0.4 pp, edge−centre gap not closed, edge eff_cut deficit closed <10 %); position-only edge forms null.
