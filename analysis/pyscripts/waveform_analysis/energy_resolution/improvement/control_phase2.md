# Phase-2 control (step 0, 2026-11-13)

Control = benchmark defaults after this step: calibrated per-module time offsets on, best calibrated gate per run type, log-weight
threshold R7ri at 2.24 and 3.49 GeV, run 24937 dropped; Møller rejection and non-linearity off. Reports: `output/energy_resolution/control_phase2_20261113/<RT>/`
(`run.sh`, `logs/`, `p2_code.sha1`; gates 5 6 8 10 12 16 20 ns, as P5). Table: `control_phase2.tsv` (`tools/p2_control.py`). Units as P5: `res_cb` = Crystal Ball σ/E of the ep pool (%),
mean ± SE over splits; `tail_hi3` = share (%) of window clusters above μ+3σ; `tail_lo_f` = share below μ−Δ_π (fixed offset, D3 guard); `leak_cut`, `eff_cut` in % (elastic cut `SEP_CUT_MEV`, frozen); pool region, control gate.

## Control table

| run type | gate | splits | res_cb in | res_cb oos | tail_hi3 | tail_lo_f | leak_cut | eff_cut | pile-up flag (%) | yield vs legacy (%) |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.73 GeV H2 | 20 ns | 16 | 3.365 ± 0.013 | 3.397 ± 0.014 | 0.049 ± 0.005 | 21.06 ± 0.25 | 0.01 ± 0.00 | 49.4 ± 0.6 | 9.7 ± 0.1 | 100.64 |
| 2.24 GeV H2 | 16 ns | 42 | 1.846 ± 0.005 | 1.858 ± 0.005 | 0.125 ± 0.004 | 26.92 ± 0.05 | 4.48 ± 0.03 | 69.2 ± 0.1 | 24.2 ± 0.1 | 100.59 |
| 3.49 GeV H2 (24698 + 24911) | 16 ns | 16 | 1.538 ± 0.008 | 1.561 ± 0.008 | 0.108 ± 0.011 | 12.08 ± 0.09 | 21.76 ± 0.12 | 84.8 ± 0.1 | 33.7 ± 0.2 | 100.61 |
| 3.49 GeV H2 old runs | 16 ns | 28 | 1.509 ± 0.006 | 1.532 ± 0.007 | 0.105 ± 0.009 | 11.64 ± 0.14 | 21.32 ± 0.10 | 85.2 ± 0.1 | 25.3 ± 0.4 | 100.46 |
| 2.24 GeV Ta X17 | 8 ns | 25 | 2.004 ± 0.023 | 2.042 ± 0.026 | 0.607 ± 0.038 | 45.4 ± 0.9 | 4.98 ± 0.15 | 50.9 ± 0.9 | 34.9 ± 1.1 | 103.36 |

Pile-up flag: the `pu` column of `summary.tsv`, the share (%) of used clusters flagged as pile-up (a measurement, not a switch); the pile-up handling itself is the scan's
(matched filter, `MIN_SEP_NS` 8, `R_MODEL` 0.02; unchanged). Log weight per row: 0.73 GeV 5.2 (null at 0.73, P4b), 2.24 GeV 3.6576, 3.49 GeV and 3.49 old 3.7073, Ta 5.2 (not measured, P4/P5).

**Against P5** (§1 "best gate + log weight": 3.365 / 1.847 / 1.538 / 1.509 / 2.004, oos 3.397 / 1.858 / 1.561 / 1.532 / 2.019): in-sample values agree to ≤ 0.001 pp (all within 0.1 SE) on every run type,
including 3.49 GeV (1.538, P5's two-run column; see below). Out of sample the four H2 rows agree to ≤ 0.001 pp; Ta is 2.042 vs 2.019 (+0.023, 0.9 SE of 0.026). Cause: Ta borrows the H2 2.24 GeV position correction
from the report root, which is now the log-weighted 2.24 GeV control (the P15 Ta report borrowed the unweighted one); consistent with the Ta legacy row moving 2.308 → 2.311 (it has no gate or log weight of its own).
The Ta control therefore depends on the 2.24 GeV control being rebuilt first.

## Gate choice (both recorded)

| run type | P5 choice (scan minimum, used) | P15 one-SE rule | oos scan in this control (5 / 6 / 8 / 10 / 12 / 16 / 20) |
|---|---|---|---|
| 0.73 | 20 ns | 12 ns | 3.414 / 3.415 / 3.419 / 3.414 / 3.406 / 3.397 / 3.397 |
| 2.24 | 16 ns | 16 ns | 1.899 / 1.887 / 1.877 / 1.869 / 1.864 / 1.858 / 1.863 |
| 3.49 | 16 ns | 12 ns | 1.623 / 1.610 / 1.585 / 1.574 / 1.565 / 1.561 / 1.563 |
| 3.49 old | 16 ns | 12 ns | 1.588 / 1.583 / 1.558 / 1.542 / 1.533 / 1.532 / 1.534 |
| Ta | 8 ns | 4 ns | 2.026 / 2.028 / 2.042 / 2.052 / 2.049 / 2.082 / 2.108 |

The oos minimum is at the P5 gate on four types (0.73: 20 and 16 tie within 0.001). On Ta the minimum of this coarse grid is 5 ns (2.026, 0.6 SE below the 8 ns gate); P15's 1-16 ns scan put it at 8 ns, kept as briefed.
(The P15 one-SE entries are quoted from P15, not recomputed on this control; the Ta one-SE gate, 4 ns, is not on this control's grid.)

## Diff summary of the default change

`waveform_analysis/energy_resolution_benchmark.py`:
- `TIME_OFFSETS = True` (was False; the constant itself was already in the tree, uncommitted, from P15).
- New `CONTROL_GATE` (ns per run type: 20 / 16 / 16 / 16 / 8). `report <scan dir>` without gates now uses it (gates argument optional; explicit gates unchanged).
- New `LOG_WEIGHT_BY_RT` (2p24GeV_H2 3.6576, 3p49GeV_H2 and _old 3.7073): applied in `clusterers.make` to every path before `CLUSTER_SET`, which still wins. Absent types keep the builder's 5.2.
- New `EXCLUDE_RUNS = (24937,)`: `report` skips those runs' scan files and prints a note. `RUN_TYPE_OF` keeps the 24937 entry (scan stays possible); `run_energy_resolution.sh` drops 24937 from its 3.49 GeV list.
- `OVERRIDABLE` gains `CONTROL_GATE`, `LOG_WEIGHT_BY_RT`, `EXCLUDE_RUNS`. The old defaults are recorded in a comment above `TIME_OFFSETS`. `overrides:` lines record `--set` only, so defaults do not appear there; `run.sh` and `p2_code.sha1` carry the provenance.
- Old control: `--set TIME_OFFSETS=False --set "LOG_WEIGHT_BY_RT={}" --set "EXCLUDE_RUNS=()"` plus the gate list.

`waveform_analysis/compare_reports.py`: `tail_lo_f` (the D3 guard) added to the paired quantities, the verdict columns and the printed summary (it was not emitted before; `tail_lo3` kept). Output gains the columns `d_tail_lo_f`, `d_tail_lo_f_se`, `verdict_tail_lo_f`.

Legacy reproduction (`repro_legacy/`, the three `--set` above, gates 8 12): `summary.tsv` and the `res_cb`, `tail_hi3`, `tail_lo_f`, `leak_cut`, `eff_cut` rows of `regions_splits.tsv` match `baseline_20261112/pre/<RT>` exactly (max |diff| 0) for 0.73, 2.24, Ta, 3.49 old and 3.49 (3.49 with `EXCLUDE_RUNS=()`, 20 splits). The no-gate default was checked on 0.73 GeV (gate 20 row identical to the control's).

## Run 24937 effect (3.49 GeV)

24937's 4 splits are NaN in every report, so dropping it removes 4 of 20 split rows and leaves the 16 valid splits; the only other channel is its 389 clusters in the gain solve.
Same defaults with `EXCLUDE_RUNS=()` (`control_phase2_20261113/with_24937/`), paired over the 16 valid splits, with − without:

| path | Δ res_cb in | Δ res_cb oos | Δ tail_lo_f | Δ leak_cut | Δ eff_cut |
|---|---|---|---|---|---|
| legacy | +0.0005 ± 0.0004 | −0.0001 ± 0.0003 | +0.000 ± 0.002 | +0.009 ± 0.007 | −0.001 ± 0.003 |
| g16 | +0.0004 ± 0.0003 | −0.0004 ± 0.0003 | −0.001 ± 0.001 | +0.014 ± 0.005 | −0.001 ± 0.002 |

Null on the resolution (≤ 0.0005 pp, ≤ 1.8 SE); `leak_cut` +0.014 ± 0.005 pp is 2.8 SE but 0.07 % of its value. The unweighted pre/ legacy row moves 1.585 → 1.557 between old and new control only through the log weight (P4b), not through 24937.
The reported split count at 3.49 GeV goes from 20 to 16 (n_pairs already 16 in P5).

(The first Ta legacy reproduction failed with rc=1 because the 2.24 GeV repro report it borrows was not yet built; rerun in `with_24937.sh`, rc=0.) The `p2_code.sha1` hashes are those at the run; afterwards only a comment in `CONTROL_GATE` changed (P15 one-SE values corrected).

Open: Ta rebuild order dependence (above); Ta gate 5 vs 8 ns.
