# Phase-2 control v2 (2026-11-13)

v2 = v1 (`control_phase2.md`) + Ta control gate 4 ns (P17) + shower-shape factor `shape_all2` (P9) at 2.24 GeV, 3.49 GeV and 3.49 GeV old runs (D13 amendment). Reports: `output/energy_resolution/control_phase2_v2_20261113/<RT>/` (`run.sh`, `logs/`, `p2v2_code.sha1`; gates 5 6 8 10 12 16 20 ns, Ta also 4). Tables: `control_phase2_v2.tsv` (control), `control_phase2_v2_delta.tsv` (v2 − v1, paired); tool `tools/p2v2_table.py`. Units as v1: `res_cb` = Crystal Ball σ/E of the ep pool (%), mean ± SE over splits; shares in %; pool region at the control gate. `tail_hi_f` (fixed-offset high tail) added to the table; the D3 high-side guard for σ-narrowing variants (supervisor ruling 11-13).

## v2 control table

| run type | gate | splits | res_cb in | res_cb oos | tail_hi_f | tail_lo_f | tail_hi3 | leak_cut | eff_cut | pile-up flag (%) | yield vs legacy (%) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.73 GeV H2 | 20 ns | 16 | 3.365 ± 0.013 | 3.397 ± 0.014 | 0.002 ± 0.001 | 21.06 ± 0.25 | 0.049 ± 0.005 | 0.01 ± 0.00 | 49.4 ± 0.6 | 9.7 ± 0.1 | 100.64 |
| 2.24 GeV H2 | 16 ns | 42 | 1.790 ± 0.004 | 1.802 ± 0.004 | 0.068 ± 0.003 | 26.86 ± 0.05 | 0.128 ± 0.004 | 4.08 ± 0.02 | 69.5 ± 0.1 | 24.2 ± 0.1 | 100.51 |
| 3.49 GeV H2 (24698 + 24911) | 16 ns | 16 | 1.501 ± 0.007 | 1.526 ± 0.007 | 0.235 ± 0.015 | 11.82 ± 0.09 | 0.126 ± 0.010 | 21.30 ± 0.11 | 85.2 ± 0.1 | 33.7 ± 0.2 | 100.47 |
| 3.49 GeV H2 old runs | 16 ns | 28 | 1.479 ± 0.006 | 1.500 ± 0.006 | 0.201 ± 0.011 | 11.35 ± 0.15 | 0.115 ± 0.006 | 21.08 ± 0.10 | 85.7 ± 0.1 | 25.4 ± 0.4 | 100.37 |
| 2.24 GeV Ta X17 | 4 ns | 25 | 2.008 ± 0.027 | 2.028 ± 0.022 | 0.282 ± 0.019 | 45.06 ± 0.91 | 0.358 ± 0.025 | 4.56 ± 0.17 | 51.2 ± 0.9 | 34.2 ± 1.1 | 103.08 |

oos gate scan (5 / 6 / 8 / 10 / 12 / 16 / 20 ns): 2.24: 1.842 / 1.834 / 1.819 / 1.811 / 1.808 / 1.802 / 1.802; 3.49: 1.582 / 1.566 / 1.546 / 1.534 / 1.529 / 1.526 / 1.527; 3.49 old: 1.548 / 1.540 / 1.516 / 1.505 / 1.502 / 1.500 / 1.501 (the minimum stays at 16 ns on all three); Ta (4 / 5 / 6 / 8 / 10 / 12 / 16 / 20): 2.028 / 2.026 / 2.028 / 2.042 / 2.052 / 2.049 / 2.082 / 2.108.

## v2 − v1 (paired per run and split, pool ep, control gate path of each; "(null)" = |Δ| < 2 SE; pp and % units)

| run type | pairs | Δ res_cb in | Δ res_cb oos | Δ tail_hi_f | Δ tail_lo_f | Δ tail_hi3 | Δ leak_cut | Δ eff_cut | Δ pile-up flag | Δ peak mu | yield change (%) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.73 GeV H2 | 16 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.00 |
| 2.24 GeV H2 | 42 | -0.057 ± 0.003 | -0.057 ± 0.003 | -0.007 ± 0.002 | -0.054 ± 0.009 | +0.002 ± 0.003 (null) | -0.40 ± 0.02 | +0.35 ± 0.02 | +0.00 ± 0.00 (null) | -0.057 ± 0.003 | +0.07 |
| 3.49 GeV H2 | 16 | -0.037 ± 0.005 | -0.035 ± 0.005 | -0.003 ± 0.009 (null) | -0.264 ± 0.023 | +0.018 ± 0.005 | -0.46 ± 0.08 | +0.40 ± 0.03 | +0.00 ± 0.01 (null) | -0.057 ± 0.006 | -0.03 |
| 3.49 GeV H2 old | 28 | -0.030 ± 0.004 | -0.032 ± 0.003 | -0.016 ± 0.006 | -0.287 ± 0.021 | +0.011 ± 0.005 (null) | -0.24 ± 0.09 | +0.48 ± 0.03 | +0.02 ± 0.01 (null) | -0.056 ± 0.005 | +0.05 |
| 2.24 GeV Ta X17 | 25 | +0.004 ± 0.024 (null) | -0.015 ± 0.022 (null) | -0.236 ± 0.026 | -0.366 ± 0.051 | -0.249 ± 0.025 | -0.42 ± 0.10 | +0.31 ± 0.08 | -0.67 ± 0.03 | -0.406 ± 0.023 | +0.06 |

Reading. The three H2 high-energy rows reproduce the P9 `shape_all2` deltas (−0.057 / −0.035 / −0.032 oos; P9: −0.057 / −0.035 / −0.032): every D3 guard holds (`tail_hi_f` lower or null, `tail_lo_f` lower, `leak_cut` lower, `eff_cut` up, pile-up flag null, yield change within ±0.1 %). `tail_hi3` is up at 3.49 GeV (+0.018 ± 0.005, the 3σ line moves in as σ narrows; the fixed-offset `tail_hi_f` is null/lower). 0.73 GeV is unchanged to the last digit (gate 20 ns, no factor, no log weight: v2 = v1 there). Ta: res null (−0.015 ± 0.022 oos), tail_hi3 −0.25, leak −0.42, pile-up flag −0.67, as P17 predicted; the peak shifts −0.41 pp (the tighter gate drops late clusters; the elastic window re-centres the fit, yield change +0.06 %).

## Ta decision: no shape factor on Ta

BM now supports a borrowed factor: `PC_EXTRA_BY_RT["2p24GeV_Ta_X17"] = "shape_all2"` makes Ta read the 2.24 GeV factor (`poscorr_extra.tsv` of the report root) on top of the borrowed map (`read_extra`); Ta then uses its own gains for the shape variables. Measured (`output/energy_resolution/variants/v2_ta_extra/`, same code, gate 4 ns; Ta with − Ta without, 25 pairs, pool, oos where marked):

| gate | Δ res_cb oos | Δ tail_hi_f | Δ tail_lo_f | Δ tail_hi3 | Δ leak_cut | Δ eff_cut | Δ peak mu | yield change |
|---|---|---|---|---|---|---|---|---|
| 4 ns | -0.032 ± 0.017 (null) | -0.032 ± 0.012 | **+1.54 ± 0.11** | -0.018 ± 0.017 | -0.13 ± 0.09 | **-1.16 ± 0.14** | +0.48 ± 0.05 | +0.94 % |
| 8 ns | -0.076 ± 0.020 | -0.105 ± 0.023 | **+1.32 ± 0.12** | -0.085 ± 0.024 | -0.50 ± 0.09 | **-0.92 ± 0.14** | +0.26 ± 0.04 | +0.86 % |

The factor is fitted on the 2.24 GeV H2 pool (not on Ta data) and moves the Ta peak by +0.5 pp, which drops `eff_cut` by 1.2 pp and raises the low tail by 1.5 pp (both guards fail by more than 8 SE, D3). The resolution gain at 4 ns is 1.8 SE (null). Not adopted: Ta stays on the borrowed 2.24 map without the factor. Refitting the factor on Ta's own pool is a possible follow-up, not done.

## 0.73 GeV decision: stays off

At 0.73 GeV (position correction "none", so the factor acts alone) `shape_all2` has a core effect of −0.005 ± 0.013 (P9) and shifts the fitted peak by −0.44 pp with 2.0 % yield loss. Likely cause (consistent with the anchor test below): the factor is normalised so that the median of E/(R E') over the fit window keeps its place, but the 0.73 GeV spectrum is strongly skewed (21 % low tail), the factor sharpens the core, and the median moves towards the tail while the peak does not. Tried a median-preserving → core-peak anchor (new switch `PC_EXTRA_ANCHOR = "core"`: iterated median inside ±`PC_EXTRA_CORE_K` × 1.4826 MAD, 6 passes; default `"median"`, which is what 2.24 / 3.49 use), 0.73 GeV, paired with the P9 control (16 pairs, pool, `variants/v2_anchor/`):

| anchor | Δ res_cb oos | Δ res_cb in | Δ tail_hi_f | Δ tail_lo_f | Δ leak_cut | Δ eff_cut | Δ peak mu | yield change |
|---|---|---|---|---|---|---|---|---|
| median (P9) | -0.005 ± 0.013 | +0.002 ± 0.011 | -0.000 ± 0.000 | -1.82 ± 0.07 | +0.0024 ± 0.0011 | +0.84 ± 0.06 | -0.441 ± 0.015 | -2.02 % |
| core, K = 1 | -0.016 ± 0.009 | -0.002 ± 0.008 | +0.000 ± 0.000 | -1.58 ± 0.08 | +0.0027 ± 0.0011 | +0.69 ± 0.06 | -0.133 ± 0.017 | -0.51 % |
| core, K = 0.5 | -0.009 ± 0.008 | +0.000 ± 0.007 | -0.000 ± 0.001 | -1.56 ± 0.08 | +0.0027 ± 0.0012 | +0.69 ± 0.05 | -0.142 ± 0.014 | +0.43 % |

The core anchor removes two thirds of the peak shift and brings the yield change inside ±1 %, but the resolution stays null (best −0.016 ± 0.009, 1.9 SE) and `leak_cut` is +0.003 pp (2.4 SE; the control value is 0.005 %). D3 adoption needs a 2 SE gain in a figure of merit: none at 0.73 GeV, and a residual peak shift of −0.13 pp is left. The only favourable moves are guards (`tail_lo_f` −1.6, `eff_cut` +0.7). 0.73 GeV remains off in v2. At 3.49 GeV the core anchor is equal to the median anchor (Δ res_cb oos −0.033 vs −0.035, peak shift −0.049 vs −0.057), so there is no reason to change the anchor of the adopted types.

Re-anchoring the gains after the factor was not tried: the gain solve takes the member energies before the factor (the factor's training input) and a gain re-solve after it would also move the position correction; with a null resolution effect it buys nothing at 0.73 GeV.

## Fixes

- `GAIN_SEED_STEP="run"` branch of `corrections` called `pc_eval(maps[m, n], p, geo)` without gains and exited when a shape factor is on: now `pc_eval(maps[m, n], p, geo, g["g"])` (member gains; the seed gains are appended afterwards). Check: `report --set 'GAIN_SEED_STEP="run"'` on `run_24698.npz` alone (8 splits of 3.49 GeV, shape factor on by default), gate 16: rc 0, `poscorr_extra.tsv` written.
- `poscorr_comparison.py`: `if run in bm.EXCLUDE_RUNS: continue` where it globs the scan npz (run 24937 no longer enters the cache and training).
- Both py_compile clean.

## Diff summary of the default change

`waveform_analysis/energy_resolution_benchmark.py` (the file also carries the uncommitted P6/P8/P9/P17 additions; only these are new here):
- `CONTROL_GATE["2p24GeV_Ta_X17"]` 8.0 → 4.0 (comment: P17).
- `PC_EXTRA` default `""` → `"auto"` (= `PC_EXTRA_BY_RT[run type]`, `""` off, a form name forces it for every run type; resolved at the start of `report`); new `PC_EXTRA_BY_RT = {2.24 H2, 3.49 H2, 3.49 H2 old: "shape_all2"}`, absent (0.73 GeV, Ta) = off. Old defaults (v1): `PC_EXTRA` off for all, Ta gate 8; recorded in the comment block above `TIME_OFFSETS` (reproduce with `--set "PC_EXTRA_BY_RT={}"` and an explicit gate list).
- A run type that borrows a position correction (Ta) reads the borrowed factor when it is listed in `PC_EXTRA_BY_RT` (`read_extra`; the old "not applied to Ta" warning is gone). Ta is not listed (see above). The borrowed report must be built with the same form first.
- New switches `PC_EXTRA_ANCHOR` (`"median"` default, `"core"`) and `PC_EXTRA_CORE_K` (1.0), via `anchor()` in `fit_extra`; defaults are byte-identical to P9. All three new constants plus `PC_EXTRA_BY_RT` are in `OVERRIDABLE`.
- Fix (a) above. `poscorr_comparison.py`: fix (b).
- New tool `improvement/tools/p2v2_table.py`.

Legacy reproduction (`control_phase2_v2_20261113/repro_legacy/`; `--set TIME_OFFSETS=False --set "LOG_WEIGHT_BY_RT={}" --set "EXCLUDE_RUNS=()" --set "PC_EXTRA_BY_RT={}"`, gates 8 12): `summary.tsv` identical and the `res_cb`, `tail_hi3`, `tail_lo_f`, `leak_cut`, `eff_cut` rows of `regions_splits.tsv` match `baseline_20261112/pre/<RT>` with max |diff| 0 on all five run types (3.49 GeV with 24937 included, 20 splits).

`p2v2_code.sha1` holds the hashes of the code at the run (benchmark, compare_reports, poscorr_comparison); `p2v2_code.sha1.at_run` is the run-time hash; the benchmark changed afterwards by one review fix (a stale `poscorr_extra.tsv` is removed when a report is built without the factor; no effect on any number) and `p2v2_code.sha1` is the current hash. `p9_run.sh` control now passes `PC_EXTRA=""`.

Open: Ta shape factor refit on Ta's own pool; 0.73 GeV factor with a residual −0.13 pp peak shift (kept off); `tail_hi3` +0.018 at 3.49 GeV is a 3σ-line effect.
