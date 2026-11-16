# P2 results note (2026-11-13)

Tables: `P2_regions.md` / `P2_regions.tsv` (from `tools/p2_tables.py`); raw: `output/energy_resolution/variants/P2/<run type>/{regions,regions_splits}.tsv`.
Run: `report --set EVAL_REGIONS=True --set MOLLER_PAIRS=True --set MOLLER_BEAM_OFFSET=...` gates 8 12 on five run types; offsets = mean of beam_offsets.tsv rows (0.73: (0.1546, 0.6638); 2.24: (0.9264, 0.9909); `3p49GeV_H2_old` = the 3488 rows (0.7452, 0.5010)); `3p49GeV_H2` and Ta have no rows: (0, 0).
Continuity: pool legacy `res_cb` equals the baseline `pre` value to 4 digits at all five run types (3.3888 / 1.8759 / 1.5854 / 1.5591 / 2.3081).

Legacy `res_cb` (%), 0.73 / 2.24 / 3.49 GeV (3.49 old in brackets):
- pool 3.389 / 1.876 / 1.585 [1.559]; edge 3.461 / 1.936 / 1.662 [1.630]; centre 3.256 / 1.844 / 1.540 [1.510]. Edge worse than centre by 0.20 / 0.09 / 0.12 pp.
- I1 3.418 / 1.875 / 1.575; I2 3.337 / 1.877 / 1.658; I3 ep-starved (330 / 945 / 176 in-range clusters, 1264 old).
- M 4.149 / 2.043 / 1.983; O 4.081 / 2.043 / 1.951; T 4.821 / 2.326 / 2.378; H 7.131 / 2.365 / 2.411 (H is excluded from the pool; prior fraction 0.47-0.81).
- sep_pi: pool 5.86 / 3.45 / 2.62; M 4.79 / 3.17 / 2.10; T 4.21 / 2.95 / 1.80. leak_cut (2.24 / 3.49): pool 4.7 / 22.3 %, M 8.6 / 25.0 %, T 10.3 / 21.6 %. tail_hi3 at 2.24 GeV: pool 0.16 %, M 1.23, O 1.46, T 3.07 (outer rings carry a larger high tail).
- Prior fraction (gains): I1 0.02-0.03, I2 0.07, I3 0.09-0.12, M/O/T 0.98-1.0 (outside the pool the gains are the prior; I3 is in the pool but starved).
- Ta (pool only; no I2/I3, O and T ep-starved): pool 2.308, g8 2.035; H is not meaningful (res_cb 19 %).

Moller:
- moller2 (E1+E2)/E_beam, Gaussian within 2 sigma, sigma/E: 4.10 / 2.15 / 2.21 [2.27] % at 0.73 / 2.24 / 3.49 GeV, Ta 3.53; peak at -3.2 / -2.5 / -1.2 [-2.2] % (3.49 SE 0.36 from one outlier split). `res_cb` is NaN for moller2 (not a CB row).
- moller E/E_ee, pool minus ep peak: -0.84 / +0.41 / +0.16 [+0.55] pp at 0.73 / 2.24 / 3.49 GeV (I1 alone: -0.83 / +0.45 / +0.27 [+0.63]; Ta -1.64). Not a flat offset across energies, so no common scale error; I2, I3, M have no usable moller fit at 2.24 and 3.49 GeV (too few); 0.73 M fit is unreliable (mu -94 %). Width is angle dominated (g2 4.5 / 2.9 / 4.7 % pool): upper bound only.

T split by lead-glass coincidence (legacy; unfitted, pooled):
| run type | T clusters with / without | share with | median E/E' with | without |
|---|---|---|---|---|
| 0.73 | 288 / 3510 | 7.6 % | 0.804 | 1.001 |
| 2.24 | 2191 / 12650 | 14.8 % | 0.877 | 0.979 |
| 3.49 | 1326 / 3497 | 27.5 % | 0.923 | 0.981 |
| 3.49 old | 1342 / 3143 | 29.9 % | 0.921 | 0.998 |
| Ta 2.24 | 85 / 405 | 17.3 % | 0.814 | 0.882 |
D7 gate input (P11): clusters seeded in the outermost two layers share energy with lead glass in 8-30 % of the cases and those lose a median 7.7-19.6 % of E/E' against 0.0-2.1 % for those without (H2). Gate condition is met (a coincidence split exists and carries a deficit well beyond the 2 % of the rest); the minimum-lead-glass-gain / fiducial-cut question for P11 is therefore open for go. The decision is the user's.

Drift (pool mu vs split order, per run, legacy): slopes are within 0.001 pp/split (|slope| <= 2 SE except 0.73 run 24590, -0.0013 +- 0.0005; Ta 25800 -0.0003 +- 0.0001, 26230 -0.0002 +- 0.0001; 2.24 run 25228 +0.0001 +- 0.00003); first-last up to -0.17 pp (24590) and -0.28 pp (Ta 25800). No drift beyond 0.05 pp at 2.24 / 3.49 GeV (U5 holds).

Not done: the Assumptions check "E/E_ee flat against azimuth" (not in the P2 block; no azimuth output). Run 24937 has no finite split (n_used 0, as in the baseline).
