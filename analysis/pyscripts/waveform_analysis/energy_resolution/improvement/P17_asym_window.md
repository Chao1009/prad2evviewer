# P17: asymmetric member time window

Date 2026-11-13. Package P17 of `notes/plan_energy_resolution_improvement_20261112.md`. Question: does a tight early / loose late member window beat the symmetric control gate (phase-2 control, `control_phase2.md`)? Decisions bound: D3 (paired, out of sample, 2 SE), D3 amendment (low-tail guard `tail_lo_f`), D4 (3.49 GeV = 24698 + 24911), D13.

## Implementation

- `--set ASYM_WINDOW="{tag: (early, late), ...}"` (new constant in `energy_resolution_benchmark.py`, default `{}`, in `OVERRIDABLE`; default behaviour byte-unchanged, the code is inert when empty). A gate G given on the command line and named in the dict (path `gG`) takes members with `-early <= t_member - t_seed <= +late`; `t` is the offset-corrected time, the seed is the event's largest hit (the first seed of the C++ energy-ordered BFS). Gate tags are only path labels: tag = late + early/100 (12.02 = early 2, late 12).
- C++ `HyCalCluster` takes one symmetric bound (`seed_time_window`; `grow_island` and `calculate_energy_square` test `|t - t_seed|`), so the late bound goes to the C++ window and the early bound is a python pre-filter in `cluster_events`, at the place where the offsets are applied: hits outside the window are moved by +1e4 ns before `add_hit`. They cannot join the seed's cluster and stay in the pool for another seed (equal to what the symmetric gate does with an out-of-window pulse). Weakness: the anchor is the event's largest hit, not each cluster's own seed; every hit of the event is filtered against the first seed, so an in-time second cluster can lose members earlier than the anchor by more than `early` (it may fragment, which can change the single/pair selection; single-cluster pool rows dominate, so the results are not affected beyond the check below); early <= late is assumed (the C++ bound `late` also cuts the early side otherwise). Check: (8, 8) against symmetric gate 8 on Ta: d res_cb_oos -0.003 ± 0.002, d tail_hi3 +0.013 ± 0.006, d pile-up -0.01 ± 0.02, yield +0.23 %, so the anchor change is a 0.01 pp effect, below every difference below.
- Sign: the brief wrote `-late <= t_member - t_seed <= +early`; with the figure's convention (dt = t_seed - t_member; late members, dt < 0, are the tail) "early" bounds the members that are earlier than the seed, so the implemented window is `-early <= t_member - t_seed <= +late` (the brief's inequality has the bounds swapped; the intent, tight on the early side, is what is implemented).
- Corrections (gains, position map) do not depend on the gate, so every variant uses the control's; Ta borrows the 2.24 GeV map (copied from the control into `variants/P17/2p24GeV_H2/` before the Ta scan; the later 2.24 GeV report rewrote it identically).
- Table tool `tools/p17_table.py` (pairs per (run, split) the pool ep rows of the variant path with the control gate path; compare_reports.py cannot pair different path names).

## Result

Pool, ep, mean ± SE over paired splits, d = variant - control, out of sample (oos) unless "in". Control gates: 0.73 GeV 20, 2.24 GeV 16, 3.49 GeV 16, 3.49 old 16, Ta 8 ns. tail_hi3, tail_lo_f, leak_cut in %, pile-up flag in %.

### Ta X17: grid (31 windows, 25 splits), d res_cb oos (pp; control 2.042 ± 0.026)

| early \ late (ns) | 6 | 8 | 10 | 12 | 16 | 20 |
|---|---|---|---|---|---|---|
| 1 | +0.158 | +0.228 | +0.255 | +0.253 | +0.299 | +0.293 |
| 2 | +0.094 | +0.091 | +0.092 | +0.089 | +0.113 | +0.111 |
| 3 | -0.005 | -0.006 | -0.008 | -0.009 | +0.008 | +0.016 |
| 4 | -0.007 | -0.013 | -0.008 | -0.010 | +0.000 | +0.021 |
| 6 | -0.014 | -0.011 | -0.003 | -0.010 | -0.000 | +0.017 |

SE 0.004 to 0.10 per cell (early 3 to 6 ns: 0.004 to 0.017; early 1 and 2: 0.03 to 0.1). d tail_hi3 (pp; control 0.607):

| early \ late | 6 | 8 | 10 | 12 | 16 | 20 |
|---|---|---|---|---|---|---|
| 1 | -0.27 | -0.23 | -0.25 | -0.07 | -0.10 | -0.00 |
| 2 | -0.29 | -0.19 | -0.17 | -0.11 | -0.01 | +0.08 |
| 3 | -0.26 | -0.14 | -0.12 | -0.08 | +0.03 | +0.12 |
| 4 | -0.17 | -0.07 | -0.06 | -0.02 | +0.10 | +0.20 |
| 6 | -0.15 | -0.03 | -0.03 | +0.02 | +0.13 | +0.21 |

d pile-up flag (pp; control 34.9): early 1: -7.0 to -5.4; early 2: about -2.5; early 3: -0.5 to -0.8; early 4: -0.3 to -0.6; early 6: -0.1 to -0.4. Yield against the control: early 1 -27 %, early 2 -2.3 to -3.0 %, early 3 to 6 +0.04 to +0.6 %.

Reading: an early bound of 1 or 2 ns cuts the core (robust sigma of the corrected dt on Ta is 1.4 ns, `P15`): resolution +0.09 to +0.30 pp worse, yield -3 to -27 %. From 3 ns on the early bound costs nothing, and late bounds of 6 to 12 ns are as good as the control 8 ns; a late bound of 16 to 20 ns is worse than the control on Ta (tail_hi3 up to +0.21, leak up) because the accidental floor is 4 times higher against the peak.

### Ta X17: symmetric gate against the asymmetric windows (settles the Ta gate)

The control grid has no 3 or 4 ns; `variants/P17/sym_extra/` adds symmetric 3 and 4 ns on the same configuration (62 s).

| window | d res_cb oos | d res_cb in | d tail_hi3 | d tail_lo_f | d leak_cut | d eff_cut | d pile-up | yield vs ctl |
|---|---|---|---|---|---|---|---|---|
| symmetric 3 ns | +0.006 ± 0.027 | +0.031 ± 0.024 | -0.344 ± 0.032 | -0.64 ± 0.10 | -0.45 ± 0.10 | +0.38 ± 0.11 | -1.17 ± 0.07 | -0.25 % |
| symmetric 4 ns | -0.015 ± 0.022 | +0.004 ± 0.024 | -0.249 ± 0.025 | -0.37 ± 0.05 | -0.42 ± 0.10 | +0.31 ± 0.08 | -0.67 ± 0.03 | +0.06 % |
| asym (4, 8) | -0.013 ± 0.007 | -0.004 ± 0.008 | -0.070 ± 0.025 | -0.15 ± 0.04 | -0.16 ± 0.05 | +0.13 ± 0.03 | -0.29 ± 0.03 | +0.58 % |
| asym (6, 8) | -0.011 ± 0.004 | -0.007 ± 0.004 | -0.031 ± 0.014 | -0.05 ± 0.01 | -0.11 ± 0.03 | +0.02 ± 0.02 | -0.10 ± 0.02 | +0.27 % |

Asymmetric windows against the symmetric 4 ns gate (paired): (4, 8) d res_cb oos +0.001 ± 0.019, d tail_hi3 +0.178 ± 0.022, d leak +0.26 ± 0.09, d pile-up +0.39 ± 0.02; (3, 6) +0.009 ± 0.017, -0.007 ± 0.021, +0.10 ± 0.08, -0.02 ± 0.05; (6, 12) +0.004 ± 0.023, +0.267 ± 0.034, +0.46 ± 0.11, +0.55 ± 0.04. No asymmetric window beats the symmetric 4 ns gate in resolution; the ones with a long late side are worse in tail_hi3, leak and pile-up, the ones with a short late side equal it (a 3 to 4 ns window is nearly symmetric already).

### H2 (best few points from the Ta scan, late side extended to 12 / 16 / 20 ns)

| run type | control | window (early, late) | res_cb oos (control) | d res_cb oos | d res_cb in | d tail_hi3 | d tail_lo_f | d leak_cut | d pile-up | yield vs ctl |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.73 GeV | 20 ns, 3.397 | (4, 20) | 3.3972 | +0.0004 ± 0.0008 | -0.0003 ± 0.0004 | -0.0017 ± 0.0012 | +0.003 ± 0.002 | -0.001 ± 0.001 | -0.000 ± 0.000 | 0.00 % |
| 0.73 GeV | | (6, 16) | 3.3974 | +0.0006 ± 0.0036 | -0.001 ± 0.005 | -0.0031 ± 0.0016 | -0.016 ± 0.037 | -0.002 ± 0.001 | +0.012 ± 0.005 | -0.16 % |
| 2.24 GeV | 16 ns, 1.8584 | (4, 16) | 1.8583 | -0.0001 ± 0.0002 | -0.0003 ± 0.0001 | -0.0075 ± 0.0012 | +0.007 ± 0.001 | -0.013 ± 0.002 | +0.009 ± 0.001 | +0.01 % |
| 2.24 GeV | | (3, 16) | 1.8584 | -0.0000 ± 0.0002 | -0.0002 ± 0.0002 | -0.0093 ± 0.0013 | +0.010 ± 0.002 | -0.013 ± 0.002 | +0.009 ± 0.001 | +0.01 % |
| 3.49 GeV | 16 ns, 1.5611 | (6, 16) | 1.5609 | -0.0002 ± 0.0001 | -0.0003 ± 0.0002 | -0.0066 ± 0.0011 | +0.003 ± 0.001 | -0.009 ± 0.003 | +0.002 ± 0.001 | +0.01 % |
| 3.49 GeV | | (3, 16) | 1.5614 | +0.0003 ± 0.0004 | +0.0002 ± 0.0003 | -0.0100 ± 0.0012 | +0.006 ± 0.003 | -0.006 ± 0.007 | +0.001 ± 0.001 | +0.00 % |
| 3.49 GeV old | 16 ns, 1.5323 | (3, 16) | 1.5310 | -0.0013 ± 0.0009 | -0.0004 ± 0.0006 | -0.0067 ± 0.0014 | +0.008 ± 0.003 | -0.015 ± 0.012 | -0.003 ± 0.002 | 0.00 % |

Other H2 windows: late 12 ns costs +0.004 to +0.010 pp (real at 0.73, 2.24 and 3.49 GeV: +0.010 ± 0.004, +0.005 ± 0.002, +0.004 ± 0.002; null at 3.49 old; yield -0.2 to -0.4 %); late 20 ns costs +0.004 ± 0.002 at 2.24 GeV (real, 2.3 SE), +0.001 to +0.002 at 3.49 GeV (null) and equals the control at 0.73 GeV (control = 20 ns there). Early 3, 4 and 6 ns are nearly indistinguishable on H2 (spread <= 0.001 pp in res_cb oos, <= 0.0035 pp in tail_hi3): the early side holds nothing beyond 3 ns once the offsets are applied.

## Best (early, late) per run type

Lowest oos d among the asymmetric points (early != late); ties within 1 SE named.

| run type | best window | d res_cb oos vs control | what it buys (tail_hi3, leak, pile-up) | verdict D3 |
|---|---|---|---|---|
| 0.73 GeV H2 | (4, 20) = control gate | +0.0004 ± 0.0008 | nothing | null |
| 2.24 GeV H2 | (4, 16); (3, 16), (6, 16) tied | -0.0001 ± 0.0002 | tail_hi3 -0.0075 ± 0.0012 (6 % of 0.125), leak -0.013, pile-up +0.009 ± 0.001 (0.04 % of 24.2) | res null; tail_hi3 real, tiny |
| 3.49 GeV H2 | (6, 16) | -0.0002 ± 0.0001 | tail_hi3 -0.0066 ± 0.0011 (6 % of 0.108), pile-up +0.002 | res null; tail_hi3 real, tiny |
| 3.49 GeV H2 old | (3, 16) | -0.0013 ± 0.0009 | tail_hi3 -0.0067 ± 0.0014 (6 % of 0.105) | null |
| 2.24 GeV Ta X17 | (4, 8); (6, 8) tied (SE 0.004) | -0.013 ± 0.007 ((6, 8): -0.011 ± 0.004) | tail_hi3 -0.070 ± 0.025 (12 % of 0.607), leak -0.16 ± 0.05, pile-up -0.29 ± 0.03 | (4, 8) null (1.8 SE); (6, 8) real but 0.5 % relative |

All guards hold for the Ta and H2 rows (pile-up flag lower or within +0.01 pp, tail_hi3 lower, tail_lo_f lower on Ta, +0.001 to +0.010 pp on H2 where it is 12 to 27 %).

## Recommendation

- Asymmetric window: do not adopt. At H2 the late bound is bounded below by the control gate (12 ns costs resolution, 16 ns is the control) and the early side is inert from 3 ns, so the best asymmetric window equals the symmetric control gate in resolution (|d| <= 0.0013 pp, within 1.5 SE) and gains 6 % of an already tiny tail_hi3 (0.007 to 0.010 pp). On Ta the early side must stay >= 3 ns (2 ns costs +0.09 pp) and nothing beats a symmetric 4 ns gate. The code stays behind `ASYM_WINDOW` (default empty) for P6/P9 users; if a tail trim is wanted at H2, (3, 16) is free.
- Ta gate (open item from step 0): symmetric 8 ns -> 4 ns. Resolution is the same (oos -0.015 ± 0.022, 5 / 6 ns at 2.026 / 2.028 in the control grid agree within errors), and 4 ns improves what D13 ranks first: tail_hi3 -0.249 ± 0.025 (41 % of 0.607), leak_cut -0.42 ± 0.10, pile-up flag -0.67 ± 0.03, tail_lo_f -0.37 ± 0.05 and eff_cut +0.31 ± 0.08 (both favourable). 3 ns gains more tail (-0.344) and costs yield (-0.25 %) with res_cb +0.006 ± 0.027; 4 ns is the compromise. Applying it means `CONTROL_GATE["2p24GeV_Ta_X17"] = 4.0` (a one-line change and a rebuilt Ta control; not done here).
- 2.24 GeV Ta X17 in the budget: the Ta resolution stays 2.00 to 2.04 % at 4 to 8 ns.

## Files

`improvement/P17_asym_window.tsv` (every run type and window, per quantity: control, variant, d, SE; Ta rows also against symmetric 4 ns `s4_d_*`; symmetric 3 and 4 ns Ta rows have `symmetric = 1`); `improvement/tools/p17_table.py`; `output/energy_resolution/variants/P17/{2p24GeV_Ta_X17 (31 windows), 0p73GeV_H2, 2p24GeV_H2, 3p49GeV_H2, 3p49GeV_H2_old (8 windows each), sym_extra/2p24GeV_Ta_X17 (symmetric 3, 4 ns), logs, run_Ta.sh, run_H2.sh, tags.txt}`. Overrides are recorded in each `regions_splits.tsv` (`# overrides: ASYM_WINDOW=...`).
