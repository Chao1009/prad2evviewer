# Phase-2 control v3 (2026-11-14)

v3 = control v2 (`control_phase2_v2.md`) + three levers the user adopted on 2026-11-14:
- P16, the below-threshold (KTeV-style) correction `BELOW_THR`, at the four H2 run types. The user waived the 2.24 GeV `tail_lo_f` guard.
- P10 form M1, the refined residual phase correction `PHASE_RESID` (global k = 1, 2 plus a per-module k = 1 table, ridge λ 100), with its own table per run type, Ta included.
- P7, the ring-scoped calibration. The M ring is calibrated at 0.73, 2.24 and 3.49 GeV old, and the H ring at 2.24, 3.49 and 3.49 GeV old. The user waived the U2 gap clause.

The P11 fiducial cut is not applied. `PHASE_CORR` stays True.

Reports: `output/energy_resolution/control_phase2_v3_20261114/`.
- `v3/<RT>/` is the headline, run with the `--set` chain in `v3_overrides.txt` and gates around the control gate.
- `v3_defaults/` is the same chain, run from the new benchmark defaults.
- `s1_p16/`, `s1_ta_nop16/`, `s2_m1/` and `s2_ta_nop16/` are the stack steps.
- `repro_v2*` and `repro_legacy*` are the defaults checks.
- `prep/` holds the Ta threshold edges and the rejected Ta tables.
- `stage*.sh`, `lib.sh` and `check_identical.sh` are the drivers, and `stage*_code.sha1` the code hashes per stage.
- `logs/walltimes.txt` and `resfit.py` / `resfit.out` complete the directory.

Tables: `control_phase2_v3.tsv` (v3 control), `control_phase2_v3_delta.tsv` (v3 − v2, paired), `control_phase2_v3_attr.tsv` (per step and region), `control_phase2_v3_gates.tsv` and `control_phase2_v3_final.tsv`. All are written by `tools/p2v3_table.py`.

Units and conventions:
- `res_cb` is the Crystal Ball σ/E of the ep pool in %, the mean ± SE over splits.
- Shares are in %.
- Δ is v3 (or a step) minus the previous chain, paired per (run, split) (D3).
- "oos" means out of sample. Gains, position map and shape factor come from the even/odd cross fit. Each split is scored with the P16 and M1 tables built on the other half of its run's splits.
- The scans are those of control v2. Bank sha1: 0.73 5d588cce, 2.24 15dfd867, 3.49 and 3.49 old 294184b0, Ta 8faf5d6f.

## (a) Final table: σ/E (%), in sample / out of sample, ep pool

| run type | gate (v3) | legacy | control v1 | control v2 | **control v3** |
|---|---|---|---|---|---|
| 0.73 GeV H2 | 20 ns | 3.389 / 3.420 | 3.365 / 3.397 | 3.365 / 3.397 | **3.155 ± 0.014 / 3.160 ± 0.018** |
| 2.24 GeV H2 | 16 ns | 1.876 / 1.889 | 1.846 / 1.858 | 1.790 / 1.802 | **1.734 ± 0.004 / 1.746 ± 0.004** |
| 3.49 GeV H2 (24698 + 24911) | 16 ns | 1.585 / 1.614 | 1.538 / 1.561 | 1.501 / 1.526 | **1.445 ± 0.006 / 1.470 ± 0.007** |
| 3.49 GeV H2, old runs | 16 ns | 1.559 / 1.579 | 1.509 / 1.532 | 1.479 / 1.500 | **1.422 ± 0.007 / 1.443 ± 0.007** |
| 2.24 GeV Ta X17 | 4 ns | 2.308 / 2.332 | 2.004 / 2.042 | 2.008 / 2.028 | **1.954 ± 0.028 / 1.987 ± 0.025** |

Gates: legacy is the ungated path of `baseline_20261112/pre` (3.49 GeV with run 24937). v1 and v2 use their control gates (v1 Ta 8 ns, v2 Ta 4 ns). SE of the earlier columns: `control_phase2_v3_final.tsv`.

The v3 control table at the control gate, out of sample:

| run type | gate | splits | res_cb in | res_cb oos | sep_pi | tail_hi_f | tail_lo_f | tail_hi3 | leak_cut | eff_cut | pile-up flag | yield vs legacy |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.73 | 20 ns | 16 | 3.155 ± 0.014 | 3.160 ± 0.018 | 6.29 | 0.001 | 19.58 | 0.056 | 0.01 | 49.4 | 9.7 | 100.66 |
| 2.24 | 16 ns | 42 | 1.734 ± 0.004 | 1.746 ± 0.004 | 3.70 | 0.068 | 26.94 | 0.136 | 3.78 | 69.6 | 24.2 | 100.55 |
| 3.49 | 16 ns | 16 | 1.445 ± 0.006 | 1.470 ± 0.007 | 2.83 | 0.226 | 11.82 | 0.148 | 20.96 | 85.3 | 33.7 | 100.52 |
| 3.49 old | 16 ns | 28 | 1.422 ± 0.007 | 1.443 ± 0.007 | 2.88 | 0.189 | 11.30 | 0.134 | 20.66 | 85.8 | 25.3 | 100.41 |
| Ta X17 | 4 ns | 25 | 1.954 ± 0.028 | 1.987 ± 0.025 | 3.27 | 0.287 | 45.10 | 0.349 | 4.49 | 51.2 | 34.2 | 103.18 |

**Best gate on the v3 chain** (`control_phase2_v3_gates.tsv`). The table lists oos `res_cb` per gate, with the paired Δ to the control gate:

| run type | gates (oos res_cb) | Δ of the neighbours vs the control gate |
|---|---|---|
| 0.73 | 16: 3.164, **20: 3.160**, 24: 3.163 | +0.004 ± 0.004, +0.003 ± 0.003 (null) |
| 2.24 | 12: 1.747, **16: 1.746**, 20: 1.749 | +0.001 ± 0.002 (null), +0.003 ± 0.001 (worse) |
| 3.49 | 12: 1.473, **16: 1.470**, 20: 1.470 | +0.003 ± 0.002, +0.001 ± 0.001 (null) |
| 3.49 old | 12: 1.441, **16: 1.443**, 20: 1.443 | −0.002 ± 0.002, +0.000 ± 0.001 (null) |
| Ta | **4: 1.987**, 5: 2.000, 6: 1.990, 8: 1.985 | +0.013 ± 0.008, +0.003 ± 0.011, −0.003 ± 0.020 (null) |

No neighbouring gate beats the control gate by 2 SE, so the control gates 20 / 16 / 16 / 16 / 4 ns stay.

## (b) Attribution per lever, out of sample, paired, pool at the control gate (`control_phase2_v3_attr.tsv`)

| run type | v2 → +P16 | +P16 → +M1 | +M1 → +P7 (pool) | total v3 − v2 | sum of isolated measurements (P16 cross / P10 M1) |
|---|---|---|---|---|---|
| 0.73 | −0.220 ± 0.013 | −0.017 ± 0.013 (null) | 0 (by construction) | **−0.237 ± 0.013** | −0.215 − 0.021 = −0.236 |
| 2.24 | −0.026 ± 0.002 | −0.030 ± 0.002 | 0 | **−0.056 ± 0.003** | −0.026 − 0.031 = −0.057 |
| 3.49 | −0.019 ± 0.003 | −0.038 ± 0.002 | 0 | **−0.057 ± 0.003** | −0.019 − 0.037 = −0.056 |
| 3.49 old | −0.017 ± 0.002 | −0.040 ± 0.004 | 0 | **−0.057 ± 0.004** | −0.018 − 0.037 = −0.055 |
| Ta | −0.001 ± 0.006 (null; only via the 2.24 map, no Ta table) | −0.039 ± 0.015 | 0 (no ring) | **−0.041 ± 0.014** | n/a (new) |

The same pairs enter every step, so the steps add up exactly to the total. On the pool, the combined gain equals the sum of the isolated P16 and P10 measurements to ≤ 0.002 pp: the two levers act independently.
- P16 dominates at 0.73 GeV and M1 at 2.24 and 3.49 GeV.
- At 0.73 GeV the M1 step is null in the stack (−0.017 ± 0.013 oos, −0.011 ± 0.006 in sample). In isolation P10 measured −0.021 ± 0.009. M1 stays on at 0.73 GeV by the user's per-energy decision, and no guard moves.
- In sample, the totals are −0.211 / −0.055 / −0.056 / −0.056 / −0.053.

`leak_cut` falls at both the P16 and the M1 step at 2.24, 3.49 and 3.49 GeV old (totals −0.39 / −0.67 / −0.68) and at the Ta M1 step (−0.24); it is null at 0.73 GeV. Per-step values are in `control_phase2_v3_attr.tsv`; `sep_pi` rises by +0.44 / +0.11 / +0.10 / +0.11 / +0.07.

**P7 acts on the rings only.** The pool rows of `v3` are byte-identical to those of `s2_m1` at every shared path and split (`regions_splits.tsv` pool / I1-I3 / edge / centre / O / T, and the whole `splits.tsv`), at all four H2 run types. Only the trained rings differ. Ring rows, out of sample:

| run type | ring | P16 + M1 on the ring (from v2) | P7 step in the stack | P7 alone on v2 (P7 note) | v2 → v3 total | ring oos res_cb v3 (pool 1.746 / 1.470 / 1.443 / 3.160) | gap Δ (U2, waived) |
|---|---|---|---|---|---|---|---|
| 0.73 | M | −0.224 ± 0.031 | −0.461 ± 0.046 | −0.459 | −0.685 ± 0.049 | 3.444 | +0.30 |
| 2.24 | H | −0.076 ± 0.003 | −0.262 ± 0.008 | −0.302 | −0.338 ± 0.008 | 1.822 | +0.12 |
| 2.24 | M | −0.065 ± 0.009 | −0.031 ± 0.009 | −0.059 | −0.096 ± 0.010 | 1.871 | +0.04 |
| 3.49 | H | −0.076 ± 0.006 | −0.648 ± 0.013 | −0.688 | −0.725 ± 0.009 | 1.663 | +0.08 |
| 3.49 old | H | −0.089 ± 0.007 | −0.804 ± 0.052 | −0.854 | −0.893 ± 0.048 | 2.187 | −0.13 |
| 3.49 old | M | −0.036 ± 0.031 (null) | −0.183 ± 0.036 | −0.211 | −0.219 ± 0.030 | 1.788 | +0.13 |

The "from v2" and "total" columns use `variants/P7/ctrl`, which is v2 with `EVAL_REGIONS` on and whose pool equals v2; `control_phase2_v3_attr.tsv` holds them as the ring rows from `v2 (P7 ctrl)`. The P7 ring guards move favourably:
- `tail_lo_f`: H −0.4 / −4.0 / −4.3, M −0.1 / −0.0 / −1.8.
- `eff_cut`: H +1.3 / +5.6 / +5.7, M +0.8 / +0.2 / +2.6.
- `leak_cut`: H −0.6 / −1.4 / −1.5, M −0.5 / −0.5 / −0.9.
- `sep_pi`: +0.06 to +0.70.

One exception is the 2.24 GeV H ring, whose `tail_hi_f` rises by +0.005 ± 0.001.

**Where the stack is less than the sum.** The ring steps are 5-15 % smaller than P7 measured alone: H −0.26 vs −0.30, −0.65 vs −0.69, −0.80 vs −0.85. The 2.24 GeV M ring step is half (−0.031 vs −0.059 pp). Two changes from P7's setting explain this, but the two are not separated:
- P16 and M1 already narrow the rings by 0.04-0.09 pp from v2, so part of what the ring gains absorbed on v2 is gone. The below-threshold loss acts like a per-module scale, which a ring gain solve absorbs.
- The adopted training sample is pool + adopted rings, not P7's pool + H + M + O + T. The M solve therefore lacks the O and T clusters whose members sit in M.

On the v2 → v3 total the rings still gain −0.10 to −0.89 pp.

## (c) Guards, v3 − v2, paired, out of sample (`control_phase2_v3_delta.tsv`)

| run type | pairs | Δ res_cb oos | Δ tail_hi_f | Δ tail_lo_f | Δ tail_hi3 | Δ leak_cut | Δ eff_cut | Δ pile-up flag | Δ peak μ | yield change | Δ gap (U2) | per run Δ res_cb oos |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.73 | 16 | −0.237 ± 0.013 | −0.001 ± 0.001 (null) | −1.45 ± 0.09 | +0.015 ± 0.008 (null) | +0.000 ± 0.001 (null) | +0.02 ± 0.04 (null) | +0.01 ± 0.01 (null) | −0.003 ± 0.015 (null) | +2.84 % | −0.026 ± 0.013 | 24590 −0.218, 24603 −0.268 |
| 2.24 | 42 | −0.056 ± 0.003 | −0.001 ± 0.001 (null) | **+0.056 ± 0.007** | +0.005 ± 0.002 | −0.39 ± 0.01 | +0.14 ± 0.01 | −0.005 ± 0.003 (null) | +0.015 ± 0.002 | −0.10 % | −0.001 ± 0.002 | −0.050 to −0.064 (5 runs) |
| 3.49 | 16 | −0.057 ± 0.003 | −0.047 ± 0.006 | −0.070 ± 0.015 | +0.016 ± 0.004 | −0.67 ± 0.06 | +0.24 ± 0.02 | +0.00 ± 0.01 (null) | +0.016 ± 0.004 | −0.32 % | −0.000 ± 0.003 | 24698 −0.060, 24911 −0.053 |
| 3.49 old | 28 | −0.057 ± 0.004 | −0.037 ± 0.005 | −0.125 ± 0.021 | +0.016 ± 0.004 | −0.68 ± 0.08 | +0.27 ± 0.03 | −0.02 ± 0.01 (null) | +0.017 ± 0.003 | −0.23 % | −0.001 ± 0.003 | −0.052 to −0.063 (3 runs) |
| Ta | 25 | −0.041 ± 0.014 | +0.007 ± 0.013 (null) | −0.00 ± 0.05 (null) | +0.004 ± 0.015 (null) | −0.24 ± 0.07 | +0.05 ± 0.07 (null) | +0.00 ± 0.01 (null) | +0.016 ± 0.012 (null) | +0.05 % | +0.012 ± 0.012 (null) | 25800 −0.029, 26138 −0.053, 26230 −0.044 |

D3 holds at every run type. `res_cb` improves by 2.9 SE (Ta) to 22 SE (2.24), on every run with the same sign. `leak_cut` is lower or null, and `tail_hi_f` is lower or null. `eff_cut` is up or null. The yield loss is ≤ 0.32 %, and the U2 gap does not grow.

Two items need reading:
- At 2.24 GeV `tail_lo_f` is +0.056 ± 0.007. This is P16's known +0.052, which the user waived; the M1 step adds a null +0.005.
- `tail_hi3` is +0.015 to +0.016 at the σ-narrowing H2 types. This is the 3σ line moving in. The D3 guard for σ-narrowing variants is `tail_hi_f`, which is lower or null.

The pool peak moves by +0.015 % (P16 adds energy, +0.018 at its step).

**Ta decisions** (D3 out of sample, Ta's own tables):
- **P16 for Ta: not adopted.** Ta got its own table, built with Ta's own thresholds (TET × MeV/ADC from the Ta scans, median 8.7 MeV against 9.5 MeV for 2.24 H2; `prep/edge_2p24GeV_Ta_X17.npz`). It was built on 1.4 × 10⁴ peak clusters per half, cross-fitted. Against Ta without it (`s1_p16` vs `s1_ta_nop16`), `res_cb` oos is −0.077 ± 0.015, `leak_cut` −0.50, `tail_hi_f` −0.05. However, `tail_lo_f` is **+0.86 ± 0.08** and `eff_cut` **−0.64 ± 0.09**, so both guards fail at 8-11 SE; the user's waiver covered 2.24 H2 only. The tables are kept in `prep/below_thr_Ta/`.
- Ta borrows the 2.24 GeV correction (P16 flagged this). With the 2.24 GeV report under P16 and no Ta table, Ta is null against v2 on every column (`res_cb` oos −0.001 ± 0.006, `tail_lo_f` −0.045 ± 0.015). The borrowed map change does not hurt Ta.
- **M1 for Ta: adopted.** The table was fitted on the Ta chain without P16 (`s1_ta_nop16` dump). The step gives `res_cb` oos −0.039 ± 0.015 (2.7 SE) and `leak_cut` −0.21 ± 0.07, with all guards null. The alternative, fitted on the own-P16 chain, is kept in `prep/M1_2p24GeV_Ta_X17_ownp16.tsv` and is not used. The `s2_m1/2p24GeV_Ta_X17` report is that rejected arm: own P16 plus the own-P16 M1 table, which then sat at the stable path.
- P7 is not used for Ta, which was not adopted there.

## (d) Resolution function on v3

The fit uses the oos `res_cb` at the control gate, the user's "best gate", and the model and beam energies of `energy_resolution_fit.py` (`resfit.py`, `resfit.out`). `resolution_fit.tsv` used the one-SE rule gates of the older reports.

| fit | a (%) | b (% √GeV) | c (% GeV) | χ²/dof |
|---|---|---|---|---|
| resolution_fit.tsv (older reports) | 1.10 ± 0.04 | 1.99 ± 0.07 | 1.62 ± 0.06 | 0/0 |
| v2, same method | 1.03 ± 0.04 | 1.89 ± 0.08 | 1.72 ± 0.05 | 0/0 |
| **v3**, 0.73 / 2.24 / 3.49 | **0.89 ± 0.05** | **2.06 ± 0.07** | **1.34 ± 0.07** | 0/0 |
| v3, with 3.49 old as a fourth point | 0.81 ± 0.04 | 2.16 ± 0.06 | 1.25 ± 0.07 | 6.9/1 |

With zero degrees of freedom the terms trade against each other. The changes have sources:
- The fall of c (1.62 → 1.34 % GeV) is the P16 gain at 0.73 GeV. P16 removes threshold-loss fluctuations, which scale like a noise term.
- The fall of a (1.10 → 0.89 %) is M1, P16 and the v2 shape factor at 2.24 and 3.49 GeV.

The two-term form a ⊕ b/√E does not describe the data (χ² 84/1).

## Defaults decision

All three levers are now benchmark defaults, as P9's `PC_EXTRA_BY_RT` is. New dicts in `energy_resolution_benchmark.py`, with the tables under `waveform_analysis/energy_resolution/tables/` (`TABLE_DIR`):
- `BELOW_THR_BY_RT` holds the four H2 run types → `below_thr/` (`<RT>_even.npz` / `_odd.npz`, the P16 cross-fit tables built on the v2 control dumps). Ta is absent.
- `PHASE_RESID_BY_RT` holds all five run types → `phase_resid/M1_<RT>.tsv`, refitted on the v3 chain: the H2 tables on the `s1_p16` dumps, Ta on `s1_ta_nop16`.
- `RING_CALIB_BY_RT` holds 0.73 ("M",), 2.24 ("H", "M"), 3.49 ("H",) and 3.49 old ("H", "M"). It sets `TRAIN_REGIONS`, `EVAL_REGIONS`, `GAIN_REGION_SCOPE` and `GAIN_COMBINE_RUNS` when `TRAIN_REGIONS` is (). A report with only one run skips it (a one-run ring solve overfits, P7).

`BELOW_THR` and `PHASE_RESID` default to `"auto"`, resolved at the start of `report` from the dicts; `""` still means off and an explicit path still works. The levers are robust:
- Every table is cross-fitted, and a split always uses the half built on the other half of its run.
- The P16 single-report cross-fit reproduces P16's two-report `kt_cross` to ≤ 0.005 pp: 0.73 −0.220 vs −0.215, and 2.24 / 3.49 / 3.49 old −0.026 / −0.019 / −0.017 vs −0.026 / −0.019 / −0.018.
- P7 cannot move the pool.

Costs and caveats:
- The combined ring solves cost time: the 2.24 GeV report took 547 s with three gates in v3, against 200 s with one gate for the same chain without rings (`s2_m1`).
- The M1 tables assume the scans' bank (t50). They must be refitted when a bank or the t0 definition changes (P10 caveat), as must the P16 tables when the thresholds change.

Reproducing older states:
- Control v2: `--set "BELOW_THR_BY_RT={}" --set "PHASE_RESID_BY_RT={}" --set "RING_CALIB_BY_RT={}"`.
- Legacy: those three plus `--set TIME_OFFSETS=False --set "LOG_WEIGHT_BY_RT={}" --set "EXCLUDE_RUNS=()" --set "PC_EXTRA_BY_RT={}"`.

## Defaults checks

- Before the defaults were switched on (code sha1 in `stage0_code.sha1`), the benchmark with all new switches at their defaults (empty dicts) reproduced control v2 to the byte (`repro_v2/` vs `control_phase2_v2_20261113/`). The comparison covered `summary`, `splits`, `regions_splits`, `gains`, `poscorr_*` and `runs`, comment lines dropped, at all five run types with v2's gate lists. This also covers P10's "defaults unchanged" claim and P16's and P7's.
- In the same state, the legacy `--set` string of v2 reproduced `baseline_20261112/pre/` to the byte (`repro_legacy/`, gates 8 12, 30 of 30 files).
- After the switch (`stage4_code.sha1`):
  - `v3_defaults` (no `--set`) equals `v3` to the byte on all 32 compared files. The one exception is Ta `regions_splits.tsv`: v3 ran Ta with `EVAL_REGIONS` on, and the default Ta report has only the pool rows, which equal v3's 125 rows exactly.
  - `repro_v2_off` (the v2 off string) equals control v2 (33 of 33 files).
  - `repro_legacy_off` (the extended legacy string) equals `baseline_20261112/pre/` (30 of 30 files).
  - With the levers on by default, the brief's four-term legacy string alone no longer reproduces `baseline_20261112/pre/`; it needs the three `*_BY_RT={}` terms (as v2 added `PC_EXTRA_BY_RT={}`). Before the switch it did reproduce it (`repro_legacy/`).
- Ta in `v3` equals `s2_ta_nop16` on every shared row (the 2.24 GeV report it borrows from differs only in the ring rows).

## Benchmark edits

Diff summary of `energy_resolution_benchmark.py` (the file was at sha1 1546db13 before this work):
- `split_parity()` is factored out of `phase_resid`.
- `below_thr(..., par)` prefers the cross-fitted pair `<RT>_even.npz` / `<RT>_odd.npz`; a split of half p uses the other half's table. The table cache is keyed by file.
- `TABLE_DIR` and the three `*_BY_RT` dicts are new. `BELOW_THR` and `PHASE_RESID` default to `"auto"`. The dicts are in `OVERRIDABLE`.
- The `report` start resolves the `"auto"` values and the ring dict before the workers fork, and prints a `[note]` line with the active levers. The ring is skipped when the report holds only one run.

Tools:
- `p16_table.py`: E′ uses the target mass of the run type (`bm.pool_mass(rt)`). The earlier code used H2 for every run type; for the H2 tables nothing changes.
- `p16_common.thresholds`: also looks in `control_phase2_v3_20261114/prep/` for the Ta edge table.
- New: `tools/p2v3_table.py`.

## (e) Open

- 2.24 GeV `tail_lo_f` +0.056: the mechanism behind P16's low-tail rise is not measured. It is waived, not explained.
- Ta P16: the own table fails the low-side guards. A Ta-specific mechanism (heavier radiative tail, different threshold classes; the T class edges 9.0 / 9.8 MeV were set for H2 and most Ta modules fall below 9.0) is not investigated.
- 0.73 GeV M1 is null in the stack (kept by decision). The 0.73 GeV P16 gain partly acts as the missing position correction (P16 note). A 0.73 GeV position-correction arm on v3 is untested.
- The 2.24 GeV M-ring step is half of P7 alone. Its cause is not separated: P16 absorbing module-scale effects, or the training set without O and T.
- The U2 gap rises on the rings (+0.04 to +0.30 pp, waived). The O and T rings stay on the prior (P14 and Møller constants are deferred).
- The resolution function has zero degrees of freedom without the 3.49 GeV old runs. The four-point fit has χ² 6.9/1, because 3.49 and 3.49 old differ by 0.027 ± 0.010 pp.
- `poscorr_comparison.py` uses `cluster_events` directly and so sees neither P16 nor M1. It is out of the chain, unchanged.
- The P7 ring constants still need a production consumer (a per-cluster seed-ring lookup).
- Second-order paths of the cross fit, as in P10 and P16 alone: the M1 targets and the P16 peak selection use the in-sample gains, and the oos gains of one half are fitted on clusters of the other half corrected with tables built on the first.
