# P8: secondary-pulse acceptance (2026-11-13)

Paths relative to `analysis/pyscripts/`. Control: phase-2 control `output/energy_resolution/control_phase2_20261113/<RT>` (D13: calibrated times,
CONTROL_GATE 16 ns H2 / 8 ns Ta, log weight R7ri at 2.24 and 3.49 GeV, 24937 excluded). Variants: `output/energy_resolution/variants/P8/<variant>/<RT>`
(`run.sh`, `logs/`, `p8_code.sha1`, benchmark copy `energy_resolution_benchmark_P8.py`). Table: `P8_secondaries_acceptance.tsv`
(`tools/p8_table.py`): paired Δ = variant − control per (run, split), mean ± SE, pool ep rows; Ta also for run 26138 alone (4 splits).
Decision rule: D3 with the 2026-11-13 guard (`tail_lo_f`); P8 kill (D7): Ta gated `res_cb` worse by > 2 SE.

## Result

**No candidate passes. No production change is recommended.** At the control gate no rule improves `res_cb` or `sep_pi` by 2 SE at any H2
energy. Both merge rules lower the low tail (Δ`eff_cut` +0.3 to +0.45 pp) but raise `tail_hi3` by +0.008 to +0.030 pp (≥ 2 SE in all 6 H2 cells)
and fail the Ta guard on 26138. The drop rule is worse on every H2 figure of merit. The pre-pulse drop is null on H2.

## Mechanism: the gated clusterer already merges in-gate extras

The scan stores every in-window pulse of a module, and the report's legacy path keeps only each module's largest. The gated path hands all pulses to
`HyCalCluster`, whose BFS (`grow_island`, `prad2det/src/HyCalCluster.cpp:184-218`) adds, for every neighbour visit, the largest unconsumed
pulse of that module within ±gate of the seed. A second in-gate pulse on a module therefore joins the cluster as an extra hit when another adjacent cluster
member is expanded (always for the seed module once a neighbour joins); a peripheral member can keep its in-gate extra out. Direct check (run 25320 geometry, `make` clusterers): seed 1000 MeV at t 0, extra 30 MeV on the seed
module, neighbour 100 MeV: g16 gives 1130 MeV at extra dt 10 ns and 1100 MeV at dt 30 ns; legacy gives 1100 MeV at both.
So the energy that P3 found missing from the legacy path (median −(e_extra/E), dt 8-30 ns) is already recovered at the control gate up to dt ≈ 16 ns.

## Candidates (cuts fixed a priori from P3 cells.tsv; `--set SECONDARY_RULE=...`, report level, no rescan)

Rule: an extra pulse (not its module's largest; matched-filter hits only) with lo < t − t_main ≤ hi and A/A_main < ar_max is removed, and its
energy is either added to the module's largest pulse (merge) or discarded (drop). On the first 3 M hits of the 25320 scan (2.24 GeV), C1 selects 80 % of the extras.
On Ta 26138 it selects 13 %, because Ta extras sit mostly in the pile-up locus.

| variant | override (`--set`) | rationale |
|---|---|---|
| C1 merge late low | `SECONDARY_RULE={"mode": "merge", "dt": (0.0, 48.0), "ar_max": 0.12}` | P3/P5 tail locus (dt 8-48, A2/A1 ≈ 0.05) carries ~1 % energy dropped by the main fit: add it back |
| C2 merge late wide | `SECONDARY_RULE={"mode": "merge", "dt": (0.0, 1.0e9), "ar_max": 0.25}` | every late pulse below the P3 pile-up line A2/A1 0.25 |
| C3 drop late low | `SECONDARY_RULE={"mode": "drop", "dt": (0.0, 48.0), "ar_max": 0.12}` | the brief's "drop late low-amplitude extras" (R_MODEL-0.3-like removal without the main-fit change) |
| C4 drop pre-pulse (added after C3) | `SECONDARY_RULE={"mode": "drop", "dt": (-1.0e9, 0.0), "ar_max": 1.0e9}` | P3 hi-tail locus dt < 0: test whether in-gate pre-pulses feed `tail_hi3` |
| S1 merge core (sensitivity) | `SECONDARY_RULE={"mode": "merge", "dt": (0.0, 30.0), "ar_max": 0.07}` | narrowest excess cells only: does the `tail_hi3` rise shrink with the merged set? |

## D3 / U7 table (control gate: g16 H2, g8 Ta; Δ ± SE in pp; yield loss in %, negative = gain)

| variant | RT | n | Δ res_cb | Δ res_cb oos | Δ sep_pi | Δ tail_hi3 | Δ leak_cut | Δ tail_lo_f | Δ eff_cut | yield loss |
|---|---|---|---|---|---|---|---|---|---|---|
| C1 | 2.24 | 42 | −0.002 ± 0.003 | −0.003 ± 0.003 | +0.003 ± 0.006 | **+0.013 ± 0.002** | **+0.029 ± 0.012** | −0.245 ± 0.006 | +0.378 ± 0.014 | −0.06 |
| C1 | 3.49 | 16 | −0.006 ± 0.004 | −0.002 ± 0.004 | +0.009 ± 0.006 | **+0.030 ± 0.004** | +0.064 ± 0.062 | −0.347 ± 0.013 | +0.424 ± 0.022 | −0.09 |
| C1 | 3.49 old | 28 | +0.001 ± 0.003 | −0.003 ± 0.004 | −0.001 ± 0.006 | **+0.015 ± 0.005** | +0.086 ± 0.044 | −0.248 ± 0.018 | +0.334 ± 0.026 | −0.14 |
| C1 | Ta | 25 | +0.015 ± 0.011 | +0.009 ± 0.015 | −0.025 ± 0.017 | +0.009 ± 0.019 | **+0.149 ± 0.055** | −0.238 ± 0.056 | +0.488 ± 0.073 | −0.87 |
| C1 | Ta 26138 | 4 | **+0.055 ± 0.024** | +0.057 ± 0.066 | **−0.095 ± 0.041** | −0.001 ± 0.079 | **+0.298 ± 0.048** | −0.326 ± 0.196 | +0.393 ± 0.090 | −1.43 |
| C2 | 2.24 | 42 | −0.004 ± 0.003 | −0.005 ± 0.003 | +0.007 ± 0.006 | **+0.015 ± 0.002** | +0.023 ± 0.013 | −0.263 ± 0.006 | +0.387 ± 0.015 | −0.07 |
| C2 | 3.49 | 16 | −0.004 ± 0.004 | −0.001 ± 0.004 | +0.006 ± 0.007 | **+0.028 ± 0.004** | +0.105 ± 0.060 | −0.362 ± 0.015 | +0.448 ± 0.019 | −0.10 |
| C2 | 3.49 old | 28 | +0.002 ± 0.004 | −0.002 ± 0.003 | −0.004 ± 0.007 | **+0.017 ± 0.005** | +0.083 ± 0.066 | −0.257 ± 0.018 | +0.338 ± 0.025 | −0.11 |
| C2 | Ta | 25 | **+0.031 ± 0.014** | +0.007 ± 0.021 | **−0.051 ± 0.022** | +0.029 ± 0.022 | **+0.184 ± 0.064** | −0.182 ± 0.065 | +0.455 ± 0.083 | −0.70 |
| C2 | Ta 26138 | 4 | **+0.070 ± 0.017** | +0.072 ± 0.074 | **−0.124 ± 0.036** | +0.025 ± 0.064 | **+0.302 ± 0.068** | −0.245 ± 0.251 | +0.336 ± 0.141 | −1.10 |
| C3 | 2.24 | 42 | **+0.014 ± 0.002** | **+0.015 ± 0.002** | **−0.026 ± 0.005** | −0.005 ± 0.002 | +0.020 ± 0.013 | **+0.252 ± 0.008** | **−0.434 ± 0.014** | +0.43 |
| C3 | 3.49 | 16 | **+0.014 ± 0.003** | **+0.013 ± 0.002** | **−0.024 ± 0.005** | −0.008 ± 0.004 | +0.082 ± 0.048 | **+0.234 ± 0.013** | **−0.329 ± 0.017** | +0.45 |
| C3 | 3.49 old | 28 | **+0.008 ± 0.003** | **+0.008 ± 0.002** | **−0.015 ± 0.005** | −0.005 ± 0.002 | +0.023 ± 0.040 | **+0.215 ± 0.015** | **−0.264 ± 0.018** | +0.33 |
| C3 | Ta | 25 | **+0.008 ± 0.004** | +0.006 ± 0.004 | **−0.013 ± 0.006** | +0.005 ± 0.006 | +0.016 ± 0.020 | +0.015 ± 0.024 | +0.010 ± 0.025 | +0.02 |
| C3 | Ta 26138 | 4 | **+0.027 ± 0.004** | **+0.019 ± 0.005** | **−0.048 ± 0.007** | −0.013 ± 0.013 | +0.049 ± 0.029 | −0.021 ± 0.019 | +0.005 ± 0.046 | 0.00 |
| C4 | 2.24 | 42 | −0.000 ± 0.000 | −0.000 ± 0.000 | +0.000 ± 0.000 | **+0.004 ± 0.001** | +0.002 ± 0.001 | **+0.007 ± 0.001** | **−0.008 ± 0.001** | 0.00 |
| C4 | 3.49 | 16 | **+0.001 ± 0.000** | +0.000 ± 0.000 | **−0.001 ± 0.000** | −0.001 ± 0.001 | +0.001 ± 0.005 | +0.004 ± 0.003 | **−0.007 ± 0.003** | 0.00 |
| C4 | 3.49 old | 28 | +0.000 ± 0.000 | +0.000 ± 0.000 | −0.000 ± 0.001 | −0.002 ± 0.001 | +0.001 ± 0.006 | +0.003 ± 0.002 | −0.002 ± 0.002 | −0.01 |
| C4 | Ta | 25 | **+0.020 ± 0.009** | +0.008 ± 0.012 | **−0.032 ± 0.016** | −0.048 ± 0.017 | +0.016 ± 0.038 | −0.226 ± 0.056 | +0.160 ± 0.058 | −1.08 |
| C4 | Ta 26138 | 4 | +0.011 ± 0.017 | −0.043 ± 0.057 | −0.024 ± 0.031 | −0.078 ± 0.030 | +0.029 ± 0.093 | +0.017 ± 0.272 | −0.126 ± 0.191 | −0.95 |
| S1 | 2.24 | 42 | −0.002 ± 0.003 | −0.003 ± 0.003 | +0.003 ± 0.006 | **+0.008 ± 0.002** | +0.015 ± 0.015 | −0.163 ± 0.006 | +0.284 ± 0.013 | −0.13 |
| S1 | 3.49 | 16 | −0.005 ± 0.003 | −0.000 ± 0.004 | +0.008 ± 0.005 | **+0.020 ± 0.004** | +0.061 ± 0.051 | −0.300 ± 0.014 | +0.352 ± 0.022 | −0.20 |
| S1 | 3.49 old | 28 | +0.002 ± 0.003 | −0.002 ± 0.003 | −0.004 ± 0.005 | +0.008 ± 0.005 | +0.057 ± 0.041 | −0.203 ± 0.015 | +0.274 ± 0.022 | −0.20 |
| S1 | Ta | 25 | +0.004 ± 0.011 | −0.004 ± 0.015 | −0.009 ± 0.018 | −0.003 ± 0.013 | +0.098 ± 0.058 | −0.202 ± 0.033 | +0.395 ± 0.057 | −0.48 |
| S1 | Ta 26138 | 4 | +0.024 ± 0.037 | +0.023 ± 0.076 | −0.044 ± 0.072 | −0.015 ± 0.037 | +0.030 ± 0.165 | −0.259 ± 0.071 | +0.332 ± 0.097 | −0.30 |

Bold: a quantity moved the wrong way by ≥ 2 SE. Legacy-path rows are in the `.tsv`: merges gain −0.01 to −0.02 pp there with `tail_hi3` up; drops
are exactly 0 because legacy already drops every extra; the legacy path is not the control. Bank sha1 (npz `bank_sha1`): 2.24 GeV 15dfd867, 3.49 and 3.49 old
294184b0, Ta 8faf5d6f. Control and variants use the same scans.

## Verdicts

- C1, C2, S1 (merge): **not adopted (D3 null on figure of merit (1), figure of merit (2) worse).** `res_cb` and `sep_pi` are null at all three H2 sets.
  `tail_hi3` rises at 2.24 and 3.49, replicated on 3.49 old (S1 old 1.7 SE). `leak_cut` rises +0.03 to +0.10 (≤ 2.5 SE). The Ta kill fires on 26138:
  C1 +0.055 ± 0.024 (2.3 SE), C2 +0.070 ± 0.017 (4.1 SE); over all 25 Ta splits C2 also fires (+0.031 ± 0.014). On the U7 metric the merges trade
  efficiency (Δ`eff_cut` +0.3 to +0.45 pp, `tail_lo_f` −0.16 to −0.36 pp) for a larger high tail and no separation gain. That is the wrong direction for
  inelastic rejection at 3.49 GeV, where `tail_hi3` and `leak_cut` set the pion-threshold leak. Shrinking the merged set (S1) shrinks both the gain and
  the `tail_hi3` cost, so no merge set was found that keeps one without the other.
- C3 (drop late low): **killed.** It is worse in `res_cb` (+0.008 to +0.014 pp, ≥ 3 SE, replicated), `tail_lo_f` (+0.22 to +0.25) and yield
  (0.33-0.45 % lost). It confirms that late low-amplitude extras are signal that the gate already collects, not spurious pulses.
- C4 (drop pre-pulses): **not adopted** (|Δ res_cb| ≤ 0.001 pp; `tail_hi3` +0.004 / −0.001 / −0.002, mixed sign, not replicated; small guard failures at 2.24). The dt < 0 hi-tail of P3's cells does not come
  from pre-pulses joining the cluster. Ta: `tail_hi3` −0.048 ± 0.017 and `tail_lo_f` −0.23 ± 0.06, but `res_cb` +0.020 ± 0.009; not interpreted (Ta is guard only).
- P8 kill (D7): no acceptance rule gives an H2 gain ≥ 2 SE in either figure of merit at the control gate. The secondary-pulse lever is exhausted at
  report level once the calibrated gate is on.

## Production change

None. The control's gated clusterer, with all in-window pulses stored, is the acceptance to keep. The `SECONDARY_RULE` switch stays in the benchmark
(default `{}`, outputs unchanged) for later re-tests, for example if the gate or the pulse fitter changes.

## Benchmark edit (flagged; additive, default unchanged)

`waveform_analysis/energy_resolution_benchmark.py`: constant `SECONDARY_RULE = {}` (after `SECONDARY_DUMP`), its name in `OVERRIDABLE`, function
`secondary_rule(h)` (before `cluster_events`), and in `split_job` `h_a` is loaded and the rule applied when the switch is set. Default regression:
`variants/P8/ctl/3p49GeV_H2` (no `--set`) reproduces the control's `splits.tsv`, `summary.tsv` and `gains.tsv` rows byte for byte, and
`compare_reports` gives every Δ = 0. Other packages (P6, P17) edited the same file while these reports ran. The drop variants' legacy rows are exactly 0
on all four run types, run from start to finish, which shows that no default moved under these runs.

## Caveats and open items

- Cuts were fixed from P3 cells (pre/ basis, legacy path) and not tuned on these reports. All five rules were judged on the same data, so a single ≥ 2 SE
  move is not by itself evidence; the conclusions rest on replicated signs.
- The Ta position correction is borrowed from the same variant's 2.24 GeV report, as a production change would do. So the Ta Δ includes the 2.24 GeV
  map change.
- Merge-then-gate (C1, C2) loses a member whose largest pulse lies outside the gate when the merged extra was inside it. This is a candidate cause of
  part of the `tail_hi3` and Ta cost, not measured. A seed-time-aware merge would need a C++ change to `grow_island`, outside this package.
- 0.73 GeV was not run (no candidate to confirm). The scan-level option (R_MODEL) stays killed (P5 iii).
