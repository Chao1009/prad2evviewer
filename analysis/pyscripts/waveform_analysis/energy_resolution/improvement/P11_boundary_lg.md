# P11: boundary clusters and lead glass (U1)

Paths relative to `analysis/pyscripts/`. Control v2 (`control_phase2_v2_20261113`, gate g20 / g16 / g16 / g16 / g4 ns, 0.73 / 2.24 / 3.49 / 3.49 old / Ta). Variant run:
`output/energy_resolution/variants/P11/run.sh` (`--set EVAL_REGIONS=True --set T_LG_STUDY=True`, benchmark copy `energy_resolution_benchmark_P11.py`, `p11_code.sha1`, `logs/`).
Offline tables: `tools/p11_table.py` (-> `P11_boundary_lg.tsv`), `tools/p11_layers.py` (-> `P11_layers.tsv`), `tools/p11_null.py`. Lead glass stays out of the acceptance and uncalibrated (U1, D4).

## Result

**Keep treatment (a) for the transition ring; do not add a lead-glass drop or a correction. The coincidence is an outermost-layer problem, and no lead-glass treatment fixes it.**
- (b) Dropping coincident clusters leaves `res_cb` and `sep_pi` unchanged within 2 SE at 0.73, 2.24 and 3.49 old (2.24 `sep_pi` 3.19 vs 3.26), and at 3.49 lowers `res_cb` by 0.144 ± 0.062 pp (2.3 SE) at a cost of 30 % of the in-range T clusters. `tail_hi3` rises (+0.04 to +0.49 pp). The clusters passing the cut fall to 98.9 / 94.5 / 80.2 / 84.2 % of (a).
- (c), (d) cannot be done as specified. The scans store no lead-glass energy: every lead-glass hit has `h_e` = 0, `h_a` = NaN (run 24698, 16564 hits), because the scan's wa_hit branch (BM:836) keeps no integral. A rescan of one run takes about 2290 s (`wall_s`), excluded by the brief. Substitute: a held-out scale factor on the coincident clusters' energy from what the scans hold: the flag (c0), the number of lead-glass modules with a hit within 60 mm (c1), c1 weighted by the nearest one's distance (d1). k is fitted on one half of the splits (even / odd ordinal in the run) so that the median of the coincident clusters' E / E' (unwindowed) equals that of the others (0.5 to 1.3), and applied to the other half.
  All three make T worse: `res_cb` +0.23 to +0.52 pp at 2.24 and 3.49 (up to 6 SE; null at 0.73), `tail_hi3` 2.8 to 7.6 %, `sep_pi` down 0.1 to 0.4 at 2.24 and above. They raise the clusters that pass the cut by 5 to 17 % while `leak_cut` rises 1 to 9 pp. The deficit of a coincident cluster is broad and a hit count does not track it. This is a null for the proxy, not for a lead-glass-energy correction.
- Seed layer: coincidences sit in the outermost PbWO4 layer (layer 0): 17 / 28 / 56 / 62 / 40 % of its clusters (0.73 / 2.24 / 3.49 / 3.49 old / Ta) against 0.2 / 1.0 / 4.3 / 3.4 / 1.4 % in layer 1 (`P11_layers.tsv`). Layer 0 is wide with or without a coincidence (IQR / 1.349 of E / E' in 0.5 to 1.3: 13.7 / 25.4 / 17.9 / 15.8 % against 5.1 / 3.4 / 4.3 / 5.0 % in layer 1).
  The seed-layer cut improves T `res_cb` by 0.43 ± 0.12 (0.73), 0.20 ± 0.05 (2.24), 0.43 ± 0.08 (3.49) pp, is null at 3.49 old, and `sep_pi` rises 4.37 -> 4.80, 3.26 -> 3.64, 1.71 -> 2.12 (1.86 -> 1.83 old). It removes 43 to 47 % of the in-range T clusters and leaves 66 to 72 % of the clusters that pass the cut.

## T ring at the control gate

Mean ± SE over splits; d = variant - (a) per split; `res_cb`, d, `tail_hi3`, `eff_cut`, `leak_cut` in %; kept = in-range clusters, % of (a); passing = kept x `eff_cut` / `eff_cut`(a). Offline fits reproduce the report's T rows (3.49: T 2.485, T_nolg 2.341).

| RT | variant | T in range kept (%) | res_cb (%) | d res_cb ± SE | tail_hi3 | sep_pi | eff_cut | leak_cut | passing cut vs (a) (%) | k even/odd |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.73 | (a) current | 100.0 | 4.640 ± 0.170 | 0.000 ± 0.000 | 0.85 | 4.37 | 75.3 | 1.0 | 100.0 |  |
| 0.73 | (b) drop lg | 93.8 | 4.637 ± 0.167 | -0.004 ± 0.037 | 0.91 | 4.37 | 79.4 | 1.0 | 98.9 |  |
| 0.73 | (c0) flag factor | 101.6 | 4.634 ± 0.220 | -0.007 ± 0.079 | 2.80 | 4.50 | 77.5 | 2.9 | 104.6 | 0.31/0.37 |
| 0.73 | (c1) n_lg factor | 101.7 | 4.628 ± 0.218 | -0.012 ± 0.070 | 2.87 | 4.49 | 77.5 | 2.8 | 104.7 | 0.29/0.36 |
| 0.73 | (d1) n_lg, distance | 101.8 | 4.722 ± 0.138 | 0.082 ± 0.091 | 2.77 | 4.26 | 77.8 | 3.1 | 105.1 | 0.33/0.40 |
| 0.73 | seed-layer cut | 57.1 | 4.213 ± 0.138 | -0.428 ± 0.116 | 0.63 | 4.80 | 86.6 | 0.7 | 65.7 |  |
| 2.24 | (a) current | 100.0 | 2.092 ± 0.071 | 0.000 ± 0.000 | 2.28 | 3.26 | 59.3 | 8.7 | 100.0 |  |
| 2.24 | (b) drop lg | 81.9 | 2.105 ± 0.064 | 0.013 ± 0.033 | 2.63 | 3.19 | 68.4 | 10.0 | 94.5 |  |
| 2.24 | (c0) flag factor | 101.5 | 2.317 ± 0.071 | 0.225 ± 0.039 | 7.07 | 2.91 | 66.4 | 16.0 | 113.7 | 0.16/0.13 |
| 2.24 | (c1) n_lg factor | 101.7 | 2.330 ± 0.072 | 0.238 ± 0.048 | 6.03 | 2.89 | 66.4 | 15.1 | 113.9 | 0.12/0.11 |
| 2.24 | (d1) n_lg, distance | 101.8 | 2.354 ± 0.061 | 0.261 ± 0.048 | 5.52 | 2.82 | 66.3 | 14.8 | 114.0 | 0.13/0.13 |
| 2.24 | seed-layer cut | 53.1 | 1.891 ± 0.063 | -0.201 ± 0.052 | 2.15 | 3.64 | 79.0 | 9.3 | 70.8 |  |
| 3.49 | (a) current | 100.0 | 2.485 ± 0.109 | 0.000 ± 0.000 | 1.38 | 1.71 | 58.5 | 23.4 | 100.0 |  |
| 3.49 | (b) drop lg | 70.0 | 2.341 ± 0.076 | -0.144 ± 0.062 | 1.87 | 1.80 | 67.0 | 25.7 | 80.2 |  |
| 3.49 | (c0) flag factor | 101.0 | 2.891 ± 0.086 | 0.405 ± 0.083 | 5.23 | 1.45 | 64.4 | 31.2 | 111.2 | 0.08/0.07 |
| 3.49 | (c1) n_lg factor | 101.2 | 2.865 ± 0.093 | 0.380 ± 0.077 | 4.26 | 1.47 | 64.2 | 30.0 | 110.9 | 0.06/0.05 |
| 3.49 | (d1) n_lg, distance | 101.4 | 2.814 ± 0.084 | 0.329 ± 0.100 | 5.06 | 1.50 | 64.5 | 30.3 | 111.7 | 0.08/0.06 |
| 3.49 | seed-layer cut | 57.0 | 2.055 ± 0.111 | -0.430 ± 0.078 | 0.86 | 2.12 | 70.0 | 23.2 | 68.1 |  |
| 3.49 old | (a) current | 100.0 | 2.336 ± 0.109 | 0.000 ± 0.000 | 0.91 | 1.86 | 54.0 | 21.1 | 100.0 |  |
| 3.49 old | (b) drop lg | 72.0 | 2.339 ± 0.122 | 0.003 ± 0.168 | 0.95 | 1.94 | 63.1 | 24.4 | 84.2 |  |
| 3.49 old | (c0) flag factor | 101.6 | 2.860 ± 0.089 | 0.524 ± 0.130 | 7.56 | 1.48 | 62.5 | 31.7 | 117.5 | 0.10/0.11 |
| 3.49 old | (c1) n_lg factor | 101.7 | 2.854 ± 0.102 | 0.518 ± 0.138 | 5.51 | 1.49 | 61.9 | 30.5 | 116.6 | 0.08/0.08 |
| 3.49 old | (d1) n_lg, distance | 101.8 | 2.841 ± 0.143 | 0.505 ± 0.190 | 5.41 | 1.62 | 61.8 | 30.3 | 116.5 | 0.09/0.09 |
| 3.49 old | seed-layer cut | 57.4 | 2.337 ± 0.094 | 0.001 ± 0.123 | 0.57 | 1.83 | 68.1 | 25.8 | 72.4 |  |
| Ta | all | - | not fitted: T has 502 clusters, lg 17.3 %, below REGION_MIN_FIT per split | | | | | | | |

## Pool null check

`tools/p11_null.py`: the pool rows (kind ep, all columns, shared run / split / path) of the five variant reports equal control v2: 0 differences over 32 / 84 / 32 / 56 / 50 rows. The treatment acts on T only.

## Recommendation under U7 (separation and tails) and fiducial consequence

1. Transition ring: keep (a). Neither the drop (b) nor a correction from lead-glass presence improves `sep_pi` or the tails.
2. Fiducial: the lever is the seed layer. Dropping the outermost PbWO4 layer raises T `sep_pi` by about 0.4 at 0.73, 2.24 and 3.49 GeV (not at 3.49 old; Ta not testable) at the price of 43 to 47 % of T clusters; T is about 1 % of the pool's in-range clusters (P2 counts), so the pool-level cost is about 0.5 %. An acceptance decision (U7), not taken here.
3. A true lead-glass correction (c) needs the raw lead-glass integral stored by a rescan (one run per energy, about 38 min each). Its ceiling is low: only layer 0 coincides and layer 0 is wide without a coincidence. Not recommended unless item 2 is refused.
4. Ta: T holds 502 clusters, below the per-split fit floor; no number.

## Caveats

- Timing: lead-glass hits carry offset 0 in calibrated time (no lead-glass rows in time_calib.json); the flag ignores time (any time, 60 mm), so it is unaffected.
- `T_lg` alone is ep-starved at 0.73 and 3.49 in `regions.tsv`; read the offline table.
- Benchmark edit: `T_LG_STUDY` (default False; records gain nlg / dlg, `T_lg` / `T_nolg` fitted like T, `t_clusters.npz`). Concurrent P7 / P10 edits are not part of it. Critic round 1: PASS (2 NITs, kept as documented).
