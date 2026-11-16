# P14: snake-scan constants through LMS ratios (2026-11-14), results of slices 1-5

Plan: `notes/plan_p14_snake_scan_20261114.md` (D1-D10). SS = `waveform_analysis/snake_scan/`, V3 = `output/energy_resolution/control_phase2_v3_20261114/`. Stages: execute, write (P14.1-P14.5 results and the P14.4 gate; P14.6 is not part of this note). Study table `P14_snake_scan.tsv` (written by `python SS/snake_scan.py transfer --rt 2p24GeV_H2 --study ER/improvement/P14_snake_scan.tsv`), per-module tables `SS/out/transfer_<rt>.tsv` (primary), `transfer_<rt>_merged.tsv` (sensitivity), `prior_<rt>.tsv`. L_scan: `SS/out/lms_scan.tsv` is the primary (23701 alone, D4 "else 23701 alone"), `lms_scan_merged.tsv` the 23701 + 23744 per-module merge kept as the sensitivity variant. The per-split LMS cache (64 MB) is outside the source tree, in `output/snake_scan/lms_npz/` (git-ignored). Tests: `SS/test_p14_1..5.py`.

## Result in one paragraph

The scan reaches 9 PbWO4 modules with a clean dwell, of which 4 are elastic-pool (W552, W557, W590, W599), 2 H (W564, W594), 1 M (W585), 1 O (W547) and 1 T (W580). The primary L_scan is the 23701-only value (D4 "else 23701 alone"). The pool spread of r = g_pool / g_T (1.48 MAD, 4 modules) is 2.29 % at 2.24 GeV and 3.03 % at 3.49 GeV against an expected 0.87 % and 0.88 %. Against the 2 % gate both are `fail`. The earlier per-module 2-SE merge with run 23744 (kept as the sensitivity variant) gives 1.57 % (`pass`, marginal) and 3.02 %; that 2.24 GeV pass is an artefact of the merge: 23744 sits about 0.65 % below 23701 (median L23744 / L23701 = 0.9935, 459 of 1151 modules disagree), and merging only the modules that agree within 2 SE turns a common factor into module-to-module scatter (for the 4 pool modules it alone moves the 2.24 GeV spread from 2.29 % to 1.57 %). (Offset recomputed as 0.9935; this argument must be re-assessed after the farm pass.) The 4-pool-module spread is fragile either way: leave-one-out spread 1.2-4.6 % (primary, 2.24 GeV), the outlier W599 sits at r = 1.09, adding the H modules gives 5.6 % (4.55 % without W594, whose 3x3 has only 7 members). The spread of g_pool itself over the same four modules is 1.1-1.3 %, so the transfer adds scatter to what a constant would give on the pool. The observed scatter exceeds the expected one by a factor 2.6-3.5 (primary). Decision for P14.6 is left to the user; the numbers are below.

## Dwell coverage (P14.1, `SS/out/dwell_23701.tsv`)

Run 23701, 9 local files (00000, 00010, ..., 00080), map `SS/hycal_map_til_run023849.json`. One dwell per file, 9 dwell modules; ramp partners W589, W586, W563 carry no edge. The EPICS `ptrans` to module-centre affine map has x: a = 1, b = 126.75 mm and y: a = -1 (slope fixed), b = 10.11 mm, residual 0 mm on the 9 dwells (`epics_affine_23701.txt`; zero residual because the module pitch quantises ptrans and a is fixed, so it is a consistency check, not a precision measurement). Dwell shares 0.83-0.98 (W564 weak, 0.83).

Region membership (benchmark regions as replicated in `SS/cosmic_mip.py:region_groups`; solved gain = row in `ring/gains.tsv` with `prior_frac < 0.5`):

| module | row, col | radius mm | region | solved g, 2.24 GeV (runs) | solved g, 3.49 GeV (runs) |
|---|---|---|---|---|---|
| W547 | 17, 3 | 301 | O | no (prior_frac 1.0) | no |
| W552 | 17, 8 | 198 | pool (I3) | 5 of 5 | 2 of 2 |
| W557 | 17, 13 | 94 | pool (I1) | 5 of 5 | 2 of 2 |
| W564 | 17, 20 | 53 | H | yes (region solve, all runs combined) | yes |
| W580 | 18, 2 | 322 | T | no | no |
| W585 | 18, 7 | 218 | M | yes (prior_frac 0.10, same in all 5 runs) | no (prior_frac 0.94-0.95) |
| W590 | 18, 12 | 115 | pool (I1) | 5 of 5 | 2 of 2 |
| W594 | 18, 16 | 33 | H | yes (g = 1.063 and 1.098, 34-466k members, prior_frac 0.02-0.05) | yes |
| W599 | 18, 21 | 73 | pool (I1) | 5 of 5 | 2 of 2 |

Counts: pool 4, H 2, M 1, O 1, T 1. Fewer than 3 would have meant `cannot be assessed`; 4 pool and 6 pool+H allow a number, but it is thin.

Gain file used: `V3/v3_defaults/<rt>/ring/gains.tsv` (the control v3 report of D7; byte-identical to `V3/v3/<rt>/ring/gains.tsv`, both checked with `cmp`). It is the ring (pool + region) solve, which is the only file that carries solved H and M gains; the pool-only `gains.tsv` carries the same pool rows within sigma but prior-dominated H rows (W564 0.999, W594 0.964 at prior_frac 0.12 and 0.96). The directories `output/energy_resolution/scan/*` hold only the per-run npz, no report. Per module g_pool is the inverse-variance mean of the run rows with prior_frac < 0.5 (H rows are one combined solve, identical in every run, so their sigma is not divided by the number of runs).

## Edge table (P14.2, `SS/out/edge_23701.tsv`)

s = E_edge / 2108.85 of the 3x3 sum in MeV of the production calibration W; s_err = hypot(fit error, half the window spread x s).

| module | region | s | s_err | window spread | n_single | n_pileup | s_central |
|---|---|---|---|---|---|---|---|
| W547 | O | 0.9874 | 0.0031 | 0.0061 | 283475 | 268880 | 0.803 |
| W552 | pool | 1.0970 | 0.0037 | 0.0065 | 280455 | 331145 | 0.909 |
| W557 | pool | 1.0237 | 0.0042 | 0.0079 | 216823 | 171041 | 0.829 |
| W564 | H | 1.1081 | 0.0097 | 0.0155 | 122366 | 112942 | 0.900 |
| W580 | T | 0.9594 | 0.0031 | 0.0063 | 341605 | 280293 | 0.785 |
| W585 | M | 1.0013 | 0.0028 | 0.0052 | 187778 | 179535 | 0.818 |
| W590 | pool | 1.0236 | 0.0032 | 0.0057 | 103248 | 85656 | 0.838 |
| W594 | H | 0.9984 | 0.0033 | 0.0065 | 247028 | 103330 | 0.853 |
| W599 | pool | 1.0289 | 0.0040 | 0.0075 | 352064 | 309568 | 0.824 |

Caveat on the plan's `edge_consistency` test. The test expected s of the 4 probed modules to agree within 1 %. It fails by design: s is measured in MeV of the production calibration W, and W was calibrated for the production HV, not for the scan; s spans 0.959-1.108 over the 9 modules (0.96-1.10 over W547, W552, W580, W599). The ADC edges agree within 4.3 % (the equalised scan HV). So s_i is not a noise measurement around one constant; it is the module's W-gain error at the scan HV, which is the input D5 uses as the gain error. s_central is 0.79-0.91 (central module only, spot not centred or containment 80-90 %), unused.

W594 flag. The 3x3 of W594 has 7 members: two of its neighbours are beam-hole positions, so its containment differs from the 9-member modules (all other eight dwell modules have 9). Its s (0.9984) and g_T are therefore not directly comparable with the others; it is flagged in `P14_snake_scan.tsv` (`flags` column) and the gate table gives the pool+H spread also without W594.

Pile-up fraction by module (n_pileup / (n_single + n_pileup)): W547 48.7 %, W552 54.1 %, W557 44.1 %, W564 48.0 %, W580 45.1 %, W585 48.9 %, W590 45.3 %, W594 29.5 %, W599 46.8 %, i.e. 29-54 % against about 10 % expected by the plan. Open systematic: the pile-up rejection may depend on pulse energy (the edge sits at the top of the amplitude distribution, so a rejection that is more efficient for large pulses would bias the single-pulse edge), and the single-pulse edge could differ from the all-event edge. A cross-check fit on all events (no pile-up rejection) was not done; s_err does not include this.

## LMS statistics (P14.3)

- L_scan (`lms_scan.tsv`, primary): 23701 alone, 25 LMS events in the 9 local files, 1778 alpha events, no file used the epoch alpha median; median statistical error 0.81 % (25 events). Variant `lms_scan_merged.tsv` adds 202 LMS events of run 23744 with the merge rule (agree within 2 SE per module): 692 of 1151 modules merged, 459 disagree and keep 23701 alone, median error 0.31 %. The disagreements are systematic: median L23744 / L23701 = 0.9935 (z up to 6.6), so the 0.3 % of the merged table is not the real uncertainty and the merge is not used for the verdict. Of the 9 dwell modules the merge touched W557, W580, W585, W594 (the other five are identical in both tables).
- L_prod (`lms_prod_<rt>.tsv`): 2.24 GeV 42 of 42 local splits (runs 25228, 25257, 25260, 25299, 25320), 2910 LMS events, 5187 alpha events, 1151 modules, median error 0.08 %; 3.49 GeV 20 of 20 splits (24698, 24911, 24937), 2424 LMS events.
- Segments (`segments_23701.tsv`): no step > 3 sigma in LMS2/LMS3 reference or pedestal across the 9 files (sigma_R2 4.3 %, sigma_R3 2.4 %, sigma_ped 0.14 ADC), so no epoch split.
- The L_scan / L_prod ratio q itself spans 0.918 (W599) to 1.142 (W594) at 2.24 GeV, far outside the 0.1-0.5 % statistical errors; this is the HV change between scan and production plus the LMS-to-shower response difference (Risks, plan).

## Transfer (P14.4), 2.24 GeV (primary L_scan = 23701 alone, N_rt = 0.98150 from 4 pool modules)

g_T = N (L_scan / L_prod) / s, sigma_T = g_T hypot(s_err/s, L_scan err, L_prod err) (N error not included). Full tables in `SS/out/transfer_2p24GeV_H2.tsv` (columns module, region, s, s_err, L_scan, L_prod, g_T, sigma_T, g_pool, r, then sigma_pool, group, exp_r, L errors, L_source).

| module | region | s | L_scan | L_prod | L_scan / L_prod | g_T | sigma_T | g_pool | r |
|---|---|---|---|---|---|---|---|---|---|
| W547 | O | 0.9874 | 8297.8 | 8436.3 | 0.9836 | 0.9777 | 0.0099 | - | - |
| W552 | pool | 1.0970 | 5406.3 | 5052.2 | 1.0701 | 0.9574 | 0.0058 | 0.9612 | 1.0040 |
| W557 | pool | 1.0237 | 2989.8 | 2869.8 | 1.0418 | 0.9989 | 0.0124 | 0.9720 | 0.9731 |
| W564 | H | 1.1081 | 1566.0 | 1502.1 | 1.0425 | 0.9234 | 0.0150 | 0.9690 | 1.0494 |
| W580 | T | 0.9594 | 4909.2 | 5066.5 | 0.9689 | 0.9913 | 0.0123 | - | - |
| W585 | M | 1.0013 | 5505.8 | 5562.6 | 0.9898 | 0.9702 | 0.0111 | 0.9668 | 0.9965 |
| W590 | pool | 1.0236 | 5633.9 | 5512.8 | 1.0220 | 0.9799 | 0.0071 | 0.9760 | 0.9960 |
| W594 (7-member 3x3) | H | 0.9984 | 4461.6 | 3925.8 | 1.1365 | 1.1172 | 0.0075 | 1.0628 | 0.9513 |
| W599 | pool | 1.0289 | 7161.0 | 7799.6 | 0.9181 | 0.8758 | 0.0066 | 0.9557 | 1.0912 |

3.49 GeV (N_rt = 0.96295; pool W552, W557, W590, W599 over runs 24698 and 24911; H W564, W594; W585 has prior_frac 0.94, not solved): r = W552 1.0200, W557 0.9800, W564 1.0022, W590 0.9791, W594 1.0209, W599 1.0709; g_T of the unsolved modules W547 0.9588, W580 0.9727, W585 0.9462 (table `transfer_3p49GeV_H2.tsv`). Sensitivity with the merged L_scan (`transfer_<rt>_merged.tsv`): 2.24 GeV N_rt unchanged (0.98150), W557 L_scan / L_prod 1.0316 instead of 1.0418 and g_T 0.9891 instead of 0.9989, W580 g_T 0.9745 instead of 0.9913, W585 0.9638, W594 r 0.9467; 3.49 GeV N_rt = 0.96760 instead of 0.96295, r W552 1.0151, W557 0.9849, W590 0.9744, W599 1.0657.

## The gate

Figure of merit: 1.48 MAD of r over the pool modules; expected = rms over modules of hypot(sigma_g/g, s_err/s, L_scan err, L_prod err); gate spread <= 2 % with at least 3 pool modules. The primary verdict uses the 23701-only L_scan; the merged L_scan is the sensitivity.

| run type | L_scan | modules | spread | sd | expected | verdict |
|---|---|---|---|---|---|---|
| 2p24GeV_H2 | 23701 only (primary) | pool, 4 | 2.29 % | 5.2 % | 0.87 % | fail |
| 2p24GeV_H2 | 23701 only | pool+H, 6 | 5.60 % | 5.1 % | 1.01 % | fail |
| 2p24GeV_H2 | 23701 only | pool+H without W594, 5 | 4.55 % | 4.7 % | 1.07 % | fail |
| 2p24GeV_H2 | 23701 only | pool+H+M, 7 | 3.48 % | 4.7 % | 1.03 % | fail |
| 2p24GeV_H2 | merged (sensitivity) | pool, 4 | 1.57 % | 4.9 % | 0.67 % | pass (marginal) |
| 2p24GeV_H2 | merged | pool+H, 6 | 4.94 % | 5.1 % | 0.88 % | fail |
| 2p24GeV_H2 | merged | pool+H without W594, 5 | 3.14 % | 4.5 % | 0.94 % | fail |
| 2p24GeV_H2 | merged | pool+H+M, 7 | 3.02 % | 4.7 % | 0.83 % | fail |
| 3p49GeV_H2 | 23701 only (primary) | pool, 4 | 3.03 % | 4.3 % | 0.88 % | fail |
| 3p49GeV_H2 | 23701 only | pool+H, 6 | 3.00 % | 3.4 % | 1.02 % | fail |
| 3p49GeV_H2 | 23701 only | pool+H without W594, 5 | 3.29 % | 3.8 % | 1.07 % | fail |
| 3p49GeV_H2 | merged (sensitivity) | pool, 4 | 3.02 % | 4.1 % | 0.67 % | fail |
| 3p49GeV_H2 | merged | pool+H, 6 | 2.23 % | 3.2 % | 0.88 % | fail |
| 3p49GeV_H2 | merged | pool+H without W594, 5 | 2.63 % | 3.6 % | 0.94 % | fail |

(`P14_snake_scan.tsv` has all 16 rows, including pool+H+M at 3.49 GeV, identical to pool+H because W585 is unsolved there, and a `flags` column marking W594.) The primary verdict is the first row of each run type (pool, 23701-only L_scan): 2.24 GeV `fail` (2.29 %), 3.49 GeV `fail` (3.03 %). The merged variant gives `pass` (1.57 %) and `fail` (3.02 %); the 2.24 GeV `pass` does not survive the choice of L_scan, because the merge converts the 0.65 % common offset between 23744 and 23701 into scatter. (Offset recomputed as 0.9935; this argument must be re-assessed after the farm pass.) The 4-pool-module spread is fragile either way: leave-one-out spread 1.2-4.6 % (primary) and 1.2-3.1 % (merged) at 2.24 GeV, 0.1-6.0 % (primary) and 1.6-6.0 % (merged) at 3.49 GeV; W599 (r = 1.09 at 2.24 GeV, 1.07 at 3.49 GeV; L_scan / L_prod = 0.918, the lowest of the nine) is the outlier at both energies; adding the H modules moves the 2.24 GeV spread to 4.5-5.6 % (primary). Without W594 (7-member 3x3) the pool+H spread is 4.55 % and 3.29 % (primary), 3.14 % and 2.63 % (merged); with W594 it is 5.60 % and 3.00 % (primary). Excess over expectation (primary): sqrt(2.29^2 - 0.87^2) = 2.1 % at 2.24 GeV, sqrt(3.03^2 - 0.88^2) = 2.9 % at 3.49 GeV, i.e. an unquantified systematic of the transfer (the LMS-to-shower response change, the edge-versus-W systematic of the s_i) of 2-3 % per module; the observed spread is 2.6 (2.24 GeV) and 3.5 (3.49 GeV) times the expected one.

Attribution (spread of x/median(x) over the 4 pool modules, 2.24 GeV | 3.49 GeV, primary L_scan): g alone 1.14 | 1.32 %; g s 1.18 | 2.36 %; g q 4.46 | 2.14 %; s q 3.11 | 4.60 %; g s q (the r spread) 2.29 | 3.03 % (q = L_prod / L_scan). Merged variant: g q 4.44 | 2.85 %; s q 2.41 | 4.58 %; g s q 1.57 | 3.02 %. Reading: g_pool is already within 1.3 % across the four pool modules; the edge factor s changes the spread little at 2.24 GeV (the pool s are 1.02-1.10), the LMS ratio q is the larger contributor (g q 4.5 % at 2.24 GeV, s q 4.6 % at 3.49 GeV), and the merge shifts the 2.24 GeV pool spread by 0.7 pp through W557 alone (r 0.9731 versus 0.9827). Edge versus LMS cannot be separated further with 4 modules: s and q partly cancel. The transfer predicts g_pool to 2.3-3.0 % per module, worse than a constant on the pool (g spread 1.1-1.3 %); it is not tested on modules where g varies more than that, because the pool g varies little.

## Prior tables (`SS/out/prior_<rt>.tsv`, O and T dwell modules only)

Only one O (W547) and one T (W580) module have an edge. sigma_g = max(pool spread x g_T, sigma_T); computed from the primary (23701-only) transfer.

| run type | module | g | sigma_g |
|---|---|---|---|
| 2p24GeV_H2 | W547 (O) | 0.9777 | 0.0223 |
| 2p24GeV_H2 | W580 (T) | 0.9913 | 0.0227 |
| 3p49GeV_H2 | W547 (O) | 0.9588 | 0.0291 |
| 3p49GeV_H2 | W580 (T) | 0.9727 | 0.0295 |

Caveat on sigma_g: it uses 1.48 MAD of only 4 values (2.29 % and 3.03 %), while the sd of r is 5.2 % (2.24 GeV) and 4.3 % (3.49 GeV), dominated by the W599 outlier; sigma_g is therefore likely overconfident by a factor 1.4-2.3 (an sd-based value would be about 0.05 at 2.24 GeV). Recommendation only: if P14.6 is run with fewer than 7 pool modules, use max(MAD, sd) of r. The 3.49 GeV table is written although its gate failed; two modules of the 252 O and 256 T modules is too few to test any region metric in P14.6 whatever the gate, and the scalar prior of 0.965 sits within 0.03 of these four values. Not used further here.

## MIP groups (P14.5, `SS/out/mip_23744.tsv`) versus the dwell groups (D9)

23744 files 00000, 00030, 00060; interior hits of vertical tracks; MIP MPV by group in MeV of the production W (`mev_w`) and nominal 0.114 MeV/ADC (`mev_nominal`; 0.114 is the plan's value from an 18.5 k ADC edge, P14.2 measured 15.4 k ADC, about 0.137 MeV/ADC, so the nominal-MeV column is about 17 % low; ratios are unaffected, noted in the `mip_23744.tsv` header). Both the MIP and the edge measure the group's response at the scan HV in production-W MeV, so the group ratio to the pool is compared directly at the scan HV (no LMS ratio enters; the LMS-transfer ratios are in the r column above).

| group | dwell modules | edge s ratio to pool (median s) | MIP mev_w ratio to pool (+- err) | MIP mev_nominal ratio | median r (2.24 / 3.49 GeV) |
|---|---|---|---|---|---|
| T | 1 | 0.935 | 0.987 +- 0.006 | 0.976 | - |
| O | 1 | 0.962 | 0.973 +- 0.006 | 1.031 | - |
| M | 1 | 0.976 | 0.997 +- 0.006 | 1.007 | 0.997 / - |
| H | 2 | 1.026 | 1.087 +- 0.019 | 1.022 | 1.000 / 1.012 |
| pool | 4 | 1.000 | 1.000 +- 0.006 | 1.000 | 1.000 / 1.000 |

The MIP group means (groups with 187-348 modules) differ from the single-module edge ratios by 1-6 %: T 0.935 vs 0.987, O 0.962 vs 0.973, M 0.976 vs 0.997, H 1.026 vs 1.087. A group of one dwell module cannot be compared with a group mean to better than the module-to-module scatter (about 3 % in s), so the check neither confirms nor refutes the edge ratios. The MIP group means among T, O, M and pool agree to 2.7 %. The H group MIP (28 of 32 modules, 580 hits) has 1.9 % error.

## Slice tests (`SS/test_p14_1..5.py`, plan criteria unmodified)

| slice | test | result | measured, diagnosis |
|---|---|---|---|
| 1 | dwell_00, dwell_30, epics_module, map_diff | pass | EPICS-to-module residual 0.00 mm over 9 dwells (quantised, see above) |
| 2 | edge_w547 | FAIL | edge 2082 MeV against 2503 MeV expected (-16.8 %): the plan's 18.5 k ADC uses another integral convention; the edge is at 15.4 k ADC in the WaveAnalyzer integral (sigma 3.74 %, below the 5 % limit) |
| 2 | edge_pileup | FAIL | 48.7 % against 2-20 %; the pile-up fraction is 29-54 % in all modules, an open systematic (edge section) |
| 2 | edge_windows | FAIL | window spread of W547 0.61 % against 0.5 %; every module 0.52-1.55 %, plateau curvature; included in s_err |
| 2 | edge_consistency | FAIL | s of W547, W552, W580, W599 spans 0.959-1.097 against 1 %; fails by design (s is the W-gain error at the scan HV, caveat above) |
| 3 | lms_count_23701, lms_ref_by_name, lms_prod_25320, segments_none | pass | 9 LMS events, 1567 channels, peak 370 ADC; LMS2 by name in both maps; 72 LMS events and L_prod for 1151 modules; no step > 3 sigma |
| 4 | transfer_units, transfer_pool_median, transfer_spread | pass | spread verdict agrees with `P14_snake_scan.tsv` (the gate itself fails, 2.29 % and 3.03 %) |
| 5 | tracks_rate, mip_groups | pass | 211 tracks in 30k events (150-300); groups of 187-348 modules (H 28, 580 hits) |
| 5 | mip_mpv | FAIL | MPV of files 00000 and 00030 differ by 2.3 ADC = 2.08 SE (200.7 +- 0.75 against 198.4 +- 0.81) against the 2 SE limit; peak median and pooled MPV (200.2 ADC) pass |

Slice 2 and the 5-mip_mpv failures are criteria of the plan that the data do not meet (a wrong ADC reference, the 2.08 SE is marginal); none changes the transfer arithmetic, but the edge failures bear on s_err (pile-up, plateau curvature).

## Caveats

- Four pool modules (six with H, seven with M) define the figure of merit; one module (W599, r = 1.09) is the outlier and the leave-one-out spread ranges 1.2-4.6 %. W594 has a 7-member 3x3 (containment differs); the pool+H spread is given with and without it. With the planned full scan (P14.7) the number of pool modules rises to about 160; nothing here replaces that.
- s_i is not a noise measurement around a constant: W is the production calibration, so s_i is the production-W error at the scan HV (D5 input); edge model systematic (window spread, 0.5-1.5 %) is in s_err, the 3x3 containment and the spot position are common only if the spot is centred; s_central (0.79-0.91) shows that module-to-module containment can differ. The pile-up fraction of 29-54 % is an open systematic of the edge (cross-check fit on all events not done).
- L_scan systematic: the 23701-only and merged values differ by up to 2.5 % per module (median L23744 / L23701 = 0.9935), larger than the 0.3 % statistical error of the merged table; the primary 23701-only values rest on 25 events (median error 0.81 %). Which of the two is closer to the truth cannot be decided here; the merge rule (per-module selective, 2 SE) turns the common offset into scatter, which is why it is the sensitivity and not the primary. N_rt is unchanged at 2.24 GeV (0.0 %) and moves by 0.5 % at 3.49 GeV; the 2.24 GeV pool spread moves by 0.7 pp.
- g_pool at 3.49 GeV is from two runs (24698, 24911), the H gains are one solve over all runs.
- MIP nominal MeV (`mev_nominal`) is about 17 % low (MEV_PER_ADC from the plan's 18.5 k ADC edge, 15.4 k measured); ratios are unaffected.
- The gate is on r, the quantity D5 defines; it tests the transfer, not the downstream effect on O/T resolution.
- Not measured: absolute scale (stays with the pool), position dependence, per-module non-linearity, spot size (Q3).
