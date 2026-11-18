# Prior sub-ns timing results vs today's t_cluster − t_veto (2026-11-07)

**Scope.** This note compares the earlier studies that reported sub-ns HyCal timing with `time_study/tcluster_minus_tveto.py`, which gives σ = 1.2–4.2 ns (`time_study/tables/tcluster_minus_tveto_fits.md`). Sources: `notes/*.md`, `output/module_timing_20261030/`, `output/tseed_dt_20261030/`, `plot_module_timing.py` (PMT), `module_timing_dump.py` (MTD), the benchmark `benchmark_hycal_timing_test_20261021.py` (B, read only), and today's benchmark logs in the session scratchpad (`ct_*_e100/<split>/benchmark.log`). I ran no analysis and changed no code. **[doc]** marks a statement taken from a note, an output file or code. **[inf]** marks my inference.

## 1. Prior sub-ns results

**Veto-referenced results.** Every one of these is **per module**, not per cluster.
Quantity: t_module − t_veto.
- t_module: MF `p_time` (fitted onset + template t50) of the module's dominant pulse, raw (MTD:12-16).
- t_veto: WaveAnalyzer CFD `pk.time` (MTD:219-226), shifted by the per-counter alignment.

The procedure:
- **Trigger.** The event passes (trig & 0xf00) ≠ 0 and (trig & 0x3000000) = 0 (PMT:115, 257).
- **Veto pulse.** The pulse has v_int > 2000 and is in time: |v_t − t_event − D0| < 10 ns. If several pulses of V1–V4 qualify, the one with the largest integral wins (PMT:196-213).
- **Definitions.** t_event is the event median of the PbWO4 times with E > 10 MeV, from at least 3 pulses (PMT:130-134). D0 is the fitted peak of v_t − t_event per run type, on a 0.5 ns histogram (PMT:454-493).
- **Run selection.** Runs without a significant D0 peak are excluded: 25228, 25257, 25260 and 24386 (`output/module_timing_20261030/run_veto.log`).
- **Counter alignment.** Per-counter shifts come from the fit median[m,k] = a_m + b_k (default `--veto-align counter`, PMT:314, 530-573).
- **Module cut.** No module energy cut (`--min-e 0`).
- **Histogram and fit.** 0.2 ns bins (PMT:101). Iterative ±2σ Gaussian core with no constant term, at most 5 refits (PMT:372-435). A module is fitted only if N > 10 000.
- **Statistic.** The median σ over the fitted modules.

| Result (median σ over modules) | 0.73 GeV H2 | 2.24 GeV H2 | Runs | Source |
|---|---|---|---|---|
| counter-aligned, raw | **0.630** (93 mod.) | **0.622** (39) | 24590+24603; 25299+25320 | `notes/PbWO4_energy_resolution_20261025.md:660` |
| same, best-50 bank (10-26 dump = `module_timing_20261030`) | 0.620 | 0.618 | same | `output/module_timing_20261030/module_timing_veto_summary.txt:5,10`; `notes/..._20261026.md:199-201` |
| counter-aligned + `--time-calib` | 0.623 / 0.631 | 0.618 | same | `notes/..._20261026.md:73`; `..._veto_tcal_summary.txt:5,10` |
| no counter alignment (`type`) | 0.730 / 0.709 | 1.016 / 1.009 | same | `notes/..._20261025.md:693`; `..._veto_typealign_summary.txt:5,12` |
| accidental-comb tooth σ (module ⊕ veto ⊕ bunch; no gate). Method: tooth centroids k = ±2…±24 over ±100 ns, fitted as kP + c | 0.569 | 0.553 | as above; 3.49: 0.558 | `notes/..._20261025.md:703-709` |

- **No sub-ns veto result exists at 3.49 GeV.** No module reached N > 10⁴ there; the best module, W630, has N = 2960 (`notes/..._20261025.md:667, 680`). **No veto result at all exists for a cluster or seed time.** The only cluster-vs-veto computation is today's script.
- **Non-veto sub-ns results, for context:**
  - t_module − t_event, calibrated: 0.53–0.58 ns (`..._20261025.md:662`).
  - Per-module RF phase: 0.46–0.58 ns over all energies, 0.18–0.30 ns above 100 MeV (`:664, 753`).
  - Cluster seed time vs RF: 0.097–0.167 ns (`..._20261026.md:206-212`).
    - Events: single-cluster, legacy path (`--cluster-times`), one split per run.
    - Quantity: t_seed − o_seed, with per-run module offsets subtracted, folded against the nearest RF hit.
    - Fit: `fit_unwrapped` on an 80-bin circular histogram.
  - Member − seed, calibrated: robust σ 0.64 / 0.84 ns at 0.73 / 3.49 GeV, against 2.28 / 2.91 ns raw (`output/tseed_dt_20261030/tables/tseed_dt.md:13-14, 21-22`).

## 2. Ingredients side by side

| Ingredient | Prior (PMT `--ref veto`, 10-25/10-30) | Today (`tcluster_minus_tveto.py`) |
|---|---|---|
| HyCal time | MF `p_time` of one module's dominant pulse | MF seed time `seed_time_ns` = `cluster.center.time` (B:2801), from MF `p_time` (B:1851) |
| HyCal calibration | raw; `--time-calib` makes no difference **per module** (0.630 → 0.623) | **raw**: the benchmark log says `time calib : off (0 modules with offset)` for 024590_00050, 025228_00100 and 024503_00050, so `toff = 0` (B:2611-2614, 1616-1625). The README and the docstring (:5) say "time_calib.json applied", which is wrong. The `HC t calib` line of the C++ setup sets `Module::time_offset`, which the clusterer and the benchmark do not apply. |
| Pooling | one module per histogram; median σ over 93 / 39 modules | all clusters' seed modules pooled into one histogram per counter |
| Veto time | CFD `pk.time` | same |
| Veto pulse selection | v_int > 2000, in time (\|Δ − D0\| < 10 ns), largest integral among V1–V4 | largest integral per counter; **no integral cut, no in-time gate** |
| Counter offsets | b_k subtracted (V1 −0.91 … V4 +1.44 ns at 0.73 GeV) | none; each counter fitted separately |
| Run selection | runs without a coincidence peak excluded | **025228 used** for 2.24 GeV |
| Trigger | 0xf00 and not LMS/alpha | none |
| Energy cut | none on the module; t_event uses E > 10 MeV | single-cluster events, E > 100 MeV, legacy path (benchmark log `signal cut`); none in the script |
| Binning | 0.2 ns | 2 ns |
| Fit | iterative ±2σ Gaussian, no constant | Gauss + constant in ±10 ns of the mode, σ start 4 ns |
| Sample | all splits of 2 runs per type; N > 10⁴ per module | one split per type; N = 1.5k–17k per counter |

## 3. Likely causes, ranked

1. **t_cluster is uncalibrated, and pooling seed modules adds the raw module offsets.**
   - [doc] The benchmark logs show time calibration off.
   - [doc] Across the 93 veto-coincident modules at 0.73 GeV, the module means of t_module − t_veto spread over a 16–84 % range of [31.41, 36.50] ns raw. With calibration the range is [34.78, 34.99] ns. At 2.24 GeV the ranges are [31.25, 36.17] and [34.28, 34.44] ns (`module_timing_veto_summary.txt:6,11`; `_tcal_summary.txt:6,11`).
   - [doc] The same switch takes member − seed from 2.28 to 0.64 ns at 0.73 GeV (`tseed_dt.md:13-14`). It takes Ta t_module − t_seed from 2.64 to 1.04 ns (`notes/..._20261026.md:72`).
   - [inf] A ±2.5 ns spread of module means, convolved with ≈ 0.6 ns, gives σ ≈ 1.5–3 ns. That alone accounts for today's 0.73 GeV widths. The prior per-module σ was immune because each histogram holds one module.
2. **The 2.24 GeV sample, run 25228, has no veto coincidence peak.**
   - [doc] Its in-time fraction is 0.39 %, against 3.0 % in 25299 and 25320. Its D0 fit is not significant (−67.69 ns), and PMT excludes the run (`output/module_timing_20261030/veto_summary.tsv` row 25228; `run_veto.log`; `notes/daq_queries_20261026.md:15`).
   - [inf] The expected coincidence is at t_cluster − t_veto ≈ −D0 = +33.4 ns, and the 0.73 and 3.49 GeV peaks sit at +32.3 to +35.4 ns. The 2.24 GeV peaks at −34 to −40 ns are not at that position, so they cannot be compared with the prior gated result. They may be a shifted correlation: the sister runs 25257 and 25260 fit at v_t − t_event ≈ +42.4 and +40.2 ns, not significant (`run_veto.log:12-13`). The cause is open (`notes/daq_queries_20261026.md:15`).
3. **There is no in-time gate and no v_int > 2000 cut, and the fit includes accidentals.**
   - [doc] At 3.49 GeV, 89 % of the veto pulses above 2000 are accidental (`..._20261025.md:680, 714`). The no-gate signal-to-background ratio in ±10 ns is 1.4 at 3.49 GeV and 8.1 at 2.24 GeV (`:723-724`). Neighbouring-bunch satellites at ±4 and ±8 ns are 10–20 % of the main peak (`:699`).
   - [inf] Taking the largest-integral pulse per counter regardless of time selects accidentals. A ±10 ns window then holds 5 bunch teeth, and Gauss + constant absorbs them into σ: 3.49 GeV V2 4.18 ns, V4 3.76 ns.
4. **V1 has two timing populations at 0.73 GeV.**
   - [doc] The populations sit at ≈ −37.0 and −33.8 ns against the event, and the origin is not established (`..._20261025.md:694`).
   - [inf] This is likely part of V1's 2.67 ns. The prior method also blended counters, but it aligned them.
5. **Binning and fit model.**
   - [inf] 2 ns bins add w²/12 = 0.33 ns² in quadrature, which turns 0.62 ns into ≈ 0.85 ns. The ±10 ns window with a free constant also widens σ against a ±2σ core fit. This is a secondary effect and cannot account for widths above about 1 ns by itself.
6. **Statistics and selection.**
   - [inf] One split per type, and no trigger mask; the prior method's mask removed the LMS/alpha spike. These are minor.

## 4. Minimal changes to reproduce the prior method

1. **Calibrated t_cluster.** Rerun the `ct_*_e100` benchmark jobs with `--time-calib database/hycal_time_offsets/time_calib.json`. Alternatively, subtract the seed module's offset in the script, mapping `seed_module_id` to time_calib.json. I have not verified that mapping.
2. **2.24 GeV sample.** Replace 025228 with a split of 25299 or 25320.
3. **Veto pulse.** Require v_int > 2000 and \|t_veto − t_ref − D0\| < 10 ns before choosing the largest integral. t_ref can be the calibrated t_cluster with D0 ≈ −33.4 ns, or PMT's `pick_veto_idx` with t_event. Apply the trigger mask from the dump `trig`.
4. **Fit.** Use 0.2 ns bins (0.5 ns at most) and PMT's `fit_core`, an iterative ±2σ Gaussian without a constant.
5. **Optional, to compare like with like.** Add a per-seed-module version: histogram per `seed_module_id` with N > 10⁴, then the median σ. Also apply the PMT counter shifts b_k if the counters are pooled.

## Gaps

- The rank between causes 1 and 3 at 3.49 GeV is a judgement call.
