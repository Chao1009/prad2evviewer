# PbWO4 energy resolution with template (matched-filter) pulses, and the timing-gate scans

**Date:** 2026-10-22
**Branch:** `waveform_jz` (baseline commit `315e4af`; everything below is uncommitted)
**Plan:** `analysis/pyscripts/plan_20261022.md` (mf method, Rounds 3–4); Phase 1 gate scan in `analysis/pyscripts/plan_20261021.md` §8
**Previous note:** `analysis/pyscripts/notes/PbWO4_energy_resolution_20261021.md` (LMS, elastic-window refinement, `RUN_TARGET`, calibration map, Legacy vs New)

Paths are relative to `analysis/pyscripts/` in `/Users/L00338853/work/PRad/prad2evviewer` unless they start with `analysis/` or `docs/`. "rawsum" is the WaveAnalyzer pulse (integral and tiered time, as before); "mf" is the template pulse of §3. Legacy and New (gated) paths are as in the 20261021 note §5. σ/E is post-refinement σ_E/E, inner ring r < 200 mm; Δ in pp.

## 1. Status

- `benchmark_hycal_timing_test_20261021.py` is still the benchmark. New flags: `--time-calib`, `--pulse-method`, `--mf-*` (templates, template-column, nsigma, min-amp, kmax, r-model, scale, int-tail-ratio, ped-floor), `--pulse-dump`, `--single-pulse`.
- `rawsum` is the default and byte-identical to the previous output (V1). mf defaults are the decided values: `--mf-r-model 0.3`, `--mf-scale 1.010`.
- `benchmark_hycal_timing_scan.py` (driver) passes `--pulse-method`, `--window` and `--pre-window` through. Defaults are unchanged and rawsum commands are byte-identical. mf outputs carry the `_mf` stem.
- Git: nothing committed since `315e4af`. Modified: the benchmark, the driver, `plan_20261021.md`. Untracked: `mf_pulse.py`, `mf_validate.py`, `gate_scan_20261021.py`, `gate_scan_20261022.py`, `plan_20261022.md`, this note. `output/` is gitignored. Other uncommitted changes (`benchmark_hycal_timing.py`, `plan_waveform.md`, bimodality and plotting scripts) are not covered here.

## 2. Timing-gate scan, rawsum (Phase 1)

Output `output/gate_scan_20261021/`; 634 jobs, 4 h 17 min at `--parallel 12`. Runs 24503, 25320, 25308, 24590 and Ta 25800+26138+26230, split into a tuning half A and an evaluation half B.

| Item | Behaviour |
|---|---|
| Time reference | Raw `pk.time` (`traw`) or with per-module offsets from `time_calib.json` (`tcal`, `--time-calib auto`). Offsets cut the seed–neighbour Δt σ_MAD from 3.10 to 1.05 ns. |
| Stage A (half A) | W ∈ {4, 8, 16, 100} raw, {2, 3, 4, 6, 8, 12, 16} corrected, per-run pre-window. W = 100 is the gate-off reference. |
| Stages B1, B2 (half A) | Pre-window LO ∈ {100, 120, 130, 140}, then HI ∈ {180, 190, 200, 220}. |
| Stage C (half B) | C1 `traw_w8_pmap`, C3 `traw_w16_p100-200`, C4 `tcal_w4_p130-200`. Floor: ε_yield and ε_peak ≥ 0.99 in every run type. |

Results:

- Tight gates make σ/E worse on e-p and C, also on corrected times: mean of 5 run types `traw_w4` +0.45, `tcal_w2` +0.14, `tcal_w4` +0.03 pp. Corrected times help only at W ≤ 8; at W = 16 raw and corrected agree within 0.002 pp.
- `traw_w8` (the old setting) fails the ε_peak floor on 3.49 and 2.24 GeV H2 and C (0.987–0.989), and on half B is +0.017 to +0.033 pp worse on e-p and C. C4 (production replay settings) is +0.024 to +0.056 pp worse there.
- LO = 100–140 is flat (0.007 pp). HI = 220 fails ε_yield on Ta (0.968).
- **Working point: raw times, W = 16 ns, pre-window 100–200 ns (`traw_w16_p100-200`).** The gate is neutral on e-p and C. On Ta it gives Δgate −0.087 ± 0.046 pp.
- Half A and half B agree for the working point in every run type, so no over-tuning flag was raised.

20 of 634 jobs (3.2%) timed out during PNG output after writing a complete `result.tsv`; their numbers are used. Both harnesses set `MPLBACKEND=Agg` to avoid plotting hangs, and the σ/E-vs-W plots use a linear W axis.

## 3. Template (matched-filter) pulse method

| Item | Behaviour |
|---|---|
| Template file | `output/pulse_templates_h50.dat`: per-module median (and MAD) of two-τ fits, τ_r and τ_f, for runs 25308 and 25320. Physics trigger `SSP_RawSum` (type 41), height ≥ 50 ADC, χ² ≤ 3, `min_pulses` 50, t0 ≥ 22 ns for 25320. 1151 rows. Column 25320, else 25308; validity τ_r 1–6, τ_f 15–35 ns; W1019, W1020, W1124 and any W module without a valid row take the PbWO4 aggregate τ_r 2.499, τ_f 24.41 ns. PbGlass and Veto keep WaveAnalyzer peaks. |
| Model | p(t) = (1 − e^{−t/τ_r}) e^{−t/τ_f} / T_max for t > 0, unit peak. |
| Search | Own pedestal (clipped mean of samples 0–29). Kernel bank per module at onset phases 0–3 ns (1 ns grid), 24 samples. Candidate channel if the 3-sample moving-average maximum exceeds max(3σ, 3 ADC). |
| Accept and subtract | Iterative. Accept if A ≥ max(5σ_A, 5 ADC) (`--mf-nsigma`, `--mf-min-amp`). A secondary under an accepted pulse's kernel span also needs A ≥ r_model·A_parent (`--mf-r-model`, 0.3). Up to `--mf-kmax` 4 pulses. Onsets closer than `MIN_SEP_NS` = 12 ns are one pulse. |
| Time | Parabola on the 1 ns score grid, then 3 Gauss–Newton steps on (a, t0). Multi-pulse channels: joint Gauss–Newton on all (a_k, t0_k). Reported t = t0 + t50 (template leading-edge 50%), the `pk.time` convention. |
| Energy | E_units = A · I_int[m](φ) · κ. I_int is the WaveAnalyzer walk-out integral of the noiseless unit template at the fitted onset, with `int_tail_ratio` (0.1) and `tail_break_n` (2) from `wave_cfg` (`--mf-int-tail-ratio` overrides). κ = `--mf-scale` 1.010. Then × f_LMS, `energize()` and the elastic-window refinement, exactly as rawsum (20261021 note §2.1). |
| Paths | Pre-window cut unchanged. Legacy takes the largest-energy mf pulse per module; New takes all. |

Validation (Phase 2). Rows marked ¹ are quoted from the Phase 2 summary of 2026-10-22; they are not in `plan_20261022.md`, and their source files were deleted (§9).

| Check | Result |
|---|---|
| Speed¹ | 1.14× the rawsum runtime (target ≤ 1.5×). |
| Noise-free fit¹ | Amplitude within 0.1%, time within 0.001 ns. |
| Template shape¹ | h50 fits better than the background shape in every module. χ²/dof 1.76: residual shape remains. |
| Integration range¹ | The walk-out rule on the template keeps 0.948 of the full template area; matches the C++ integral to 0.004%. |
| Energy scale (V6) | Pooled E_mf/E_rawsum 0.9894 (r_model 0.05), 0.9898 (r_model 0.3, module spread 0.64%) → κ = 1.010. |
| Low amplitude¹ | rawsum captures only 88–90% of the charge at 5–20 ADC; mf does not lose it. |
| Timing, raw times¹ | mf narrower than rawsum in every neighbour-energy bin, e.g. 3.03 → 2.73 ns at low neighbour energy. |
| Extra pulses | Without `MIN_SEP_NS`, ~15% of large pulses split into a main pulse and a −6 ns "pre-pulse" (ratio ~0.25), which tripped the V11 stop rule (+0.30 pp on 025320 legacy Gauss). |
| Pile-up truth (V9) | Fails the ≤ 3% amplitude-bias target (+4%)¹; both pulses are found after a pedestal-window artefact in the test was fixed. 5-sample-pedestal diagnostic: both found 0% at ΔT = 16 ns, 80–100% at 20 ns, so the two-pulse resolution is ~20 ns (`MIN_SEP_NS` is 12 ns). Merged 12–20 ns pairs: 67 of 383,734 mf pulses on 025320. |

**r_model 0.3** replaces 0.05. Multi-peak channels drop from 13.4–14.0% to 3.6–4.0% (0.26% for rawsum¹). Preferred-fit post σ/E changes by −0.002 ± 0.002 pp (r0.3 − r0.05), and single-cluster yields agree to ≤ 0.02%. V6 then gives κ = 1.010. Every later mf job uses r_model 0.3 and κ 1.010; the r_model 0.05 Phase 3 run is archived in `output/gate_scan_20261022/phase3_r005/`.

## 4. Energy resolution with mf

### 4.1 Gate-matched comparison (Phase 3)

`output/gate_scan_20261022/resolution_table.md`. Both methods run fresh at the Phase 1 working point (raw times, W = 16, pre-window 100–200), 11 non-X17 runs, 93 splits. Δ = mf − rawsum of the preferred fit, paired over splits, mean ± SE.

| Run type (fit) | Splits | Legacy rawsum → mf (%) | Δ Legacy (pp) | New rawsum → mf (%) | Δ New (pp) |
|---|---|---|---|---|---|
| 3.49 GeV H2 (CB) | 27 | 1.94 → 1.87 | −0.061 ± 0.006 | 1.94 → 1.87 | −0.062 ± 0.006 |
| 2.24 GeV H2 (Gauss) | 42 | 2.26 → 2.22 | −0.040 ± 0.002 | 2.26 → 2.22 | −0.042 ± 0.002 |
| 2.24 GeV C (CB) | 8 | 2.17 → 2.09 | −0.073 ± 0.006 | 2.17 → 2.09 | −0.075 ± 0.006 |
| 0.73 GeV H2 (Gauss) | 16 | 4.66 → 4.35 | −0.306 ± 0.020 | 4.65 → 4.35 | −0.303 ± 0.020 |

mf is better in every run (per run: −0.032 to −0.077 pp at 2.24 and 3.49 GeV, −0.294 and −0.313 pp at 0.73 GeV). Ta is in §5.

In the Phase 3 tables of 24590 (108 → 113 modules), 24603 (93 → 96), 25308 (173 → 174) and 26230 (10 → 13), mf uses extra modules that rawsum excluded under the Møller-gap rule (20261021 note §2.1), released by its narrower σ_glob, so Δ there is not a same-module comparison.

### 4.2 mf full scan through the driver

`output/scan_20261022_mf/resolution_mf_vs_rawsum_fullscan.md`: 119 jobs, 14 runs (24578 excluded), 20:20–21:18, 0 failures. X17 runs use the combined X17 templates (§5).

Preferred-fit (★) legacy σ/E (%), rawsum (`output/scan_20261021_full/`) → mf, mean ± split SD. Δ = mf − rawsum, unpaired difference of split means (pp). CB values for every run are in the source file.

| Run | Type | Splits raw/mf | ★ fit | Legacy σ/E | Δ legacy | Δ gated |
|---|---|---|---|---|---|---|
| 24386 | 3.49 GeV H2 | 10/11 | CB | 1.93 ± 0.02 → 1.85 ± 0.04 | −0.078 | −0.117 |
| 24503 | 3.49 GeV H2 | 10/10 | CB | 1.94 ± 0.02 → 1.89 ± 0.02 | −0.057 | −0.091 |
| 24512 | 3.49 GeV H2 | 7/7 | CB | 1.94 ± 0.01 → 1.89 ± 0.01 | −0.046 | −0.076 |
| 24590 | 0.73 GeV H2 | 10/10 | Gauss | 4.60 ± 0.12 → 4.28 ± 0.13 | −0.312 | −0.332 |
| 24603 | 0.73 GeV H2 | 6/6 | Gauss | 4.75 ± 0.22 → 4.46 ± 0.21 | −0.291 | −0.322 |
| 25228 | 2.24 GeV H2 | 11/11 | Gauss | 2.26 ± 0.01 → 2.22 ± 0.03 | −0.035 | −0.070 |
| 25257 | 2.24 GeV H2 | 7/7 | Gauss | 2.27 ± 0.02 → 2.21 ± 0.02 | −0.053 | −0.082 |
| 25260 | 2.24 GeV H2 | 8/8 | Gauss | 2.27 ± 0.02 → 2.22 ± 0.02 | −0.042 | −0.076 |
| 25299 | 2.24 GeV H2 | 8/8 | Gauss | 2.26 ± 0.02 → 2.23 ± 0.01 | −0.032 | −0.064 |
| 25308 | 2.24 GeV C | 8/8 | CB | 2.17 ± 0.01 → 2.09 ± 0.01 | −0.073 | −0.095 |
| 25320 | 2.24 GeV H2 | 8/8 | Gauss | 2.27 ± 0.01 → 2.23 ± 0.02 | −0.042 | −0.075 |
| 25800 | 2.24 GeV Ta | 8/8 | CB | 2.33 ± 0.06 → 2.22 ± 0.07 | −0.113 | −0.148 |
| 26138 | 2.24 GeV Ta | 4/4 | CB | 2.43 ± 0.11 → 2.23 ± 0.06 | −0.203 | −0.344 |
| 26230 | 2.24 GeV Ta | 13/13 | CB | 2.48 ± 0.24 → 2.48 ± 0.31 | +0.000 | +0.054 |

The rawsum full scan used W = 8 and the per-run pre-window (130–200 ns for 24xxx), so only the legacy Δ of the 25xxx and X17 runs is gate-matched; gated Δ and the 24xxx legacy Δ mix gate and pulse-method effects. Gate-matched numbers: §4.1.

### 4.3 Single-pulse events

`--single-pulse cluster` keeps a single-cluster event only if every member module has exactly one pulse in its full readout; `seed` requires this of the seed only (seed ran on 2 splits per run). Preferred-fit legacy σ/E, rawsum → mf (%), Δ (pp); `output/gate_scan_20261022/resolution_compare.md`.

| Run type (fit) | Cluster | Seed |
|---|---|---|
| 3.49 GeV H2 (CB) | 1.93 → 1.92, −0.010 ± 0.008 | 1.94 → 1.91, −0.023 ± 0.010 |
| 2.24 GeV H2 (Gauss) | 2.26 → 2.24, −0.016 ± 0.003 | 2.26 → 2.24, −0.015 ± 0.005 |
| 2.24 GeV C (CB) | 2.17 → 2.11, −0.053 ± 0.008 | 2.17 → 2.11, −0.053 ± 0.001 |
| 0.73 GeV H2 (Gauss) | 4.65 → 4.49, −0.160 ± 0.024 | 4.63 → 4.50, −0.134 ± 0.061 |
| Kept (legacy, rawsum / mf) | 98.1–99.7% / 48.0–73.1% | 99.6–99.9% / 80.0–87.0% |

mf reports more secondary pulses (4% of channels at r_model 0.3), so the cut keeps different events and modules for the two methods (25320 cluster: 129 vs 82 modules). This table therefore does not isolate the multi-pulse contribution to the all-pulse gain (§4.1).

## 5. X17

| Item | 26138 only | Combined 25800 + 26138 + 26230 |
|---|---|---|
| File (column) | `output/pulse_templates_x17_26138_h50.dat` (26138) | `output/pulse_templates_x17_combined_h50.dat` (27000, a label for the combined sample) |
| Input | 4 splits, 21,243 RawSum events | 25 splits, 173,299 RawSum events |
| Modules with ≥ 50 pulses | 170 (median 208 pulses) | 377 of 1151 (median 304); 20 with 0 pulses |
| Aggregate τ_r / τ_f (ns), χ²/dof | 2.759 / 23.085, 1.13 | 2.653 / 23.62, 1.17 |
| Per-module median vs 25320, τ_r / τ_f | +4.3% / −2.3% | +3.5% / −1.4% |
| Used in | Phase 3 RX and Phase 4 (§7) | `--mf-templates auto` now; mf full scan (§4.2) |

Cuts are identical to 25320. Modules below 50 pulses take the X17 aggregate.

**Statistics (rough estimate).** Pulses per module per split: median 0.60 (10/25% quantiles 0.12/0.20). Reaching ≥ 50 pulses needs about 83 splits for 50% of modules, 250 for 75% and 420 for 90%. Local data has 25 X17 splits of about 101k events each.

**Resolution.** Legacy CB σ/E, rawsum → mf (%), Δ (pp). Phase 3 RX columns are gate-matched (raw W 16, 100–200) with the 26138 templates. The full-scan column uses the combined templates (§4.2).

| Run | Splits | Modules raw/mf | All pulses (26138 tmpl.) | Cluster cut (kept raw/mf %; modules) | Full scan (combined tmpl.) |
|---|---|---|---|---|---|
| 25800 | 8 | 34/33 | 2.33 → 2.21, −0.125 ± 0.024 | 2.28 → 2.15, −0.138 ± 0.055 (64/37; 28/16) | 2.33 → 2.22, −0.113 |
| 26138 | 4 | 18/18 | 2.43 → 2.23, −0.203 ± 0.065 | 2.36 → 1.93, −0.428 ± 0.046 (64/34; 14/4) | 2.43 → 2.23, −0.203 |
| 26230 | 13 | 10/13 | 2.48 → 2.56, +0.075 ± 0.060 | rawsum 2.70; mf no fit (45/22; 2/0) | 2.48 → 2.48, +0.000 |
| Ta pooled | 25 (cluster 12) | | −0.034 ± 0.040 (gated −0.098 ± 0.039) | −0.235 ± 0.056 (gated −0.222 ± 0.056) | |

These numbers rest on 2–34 modules per run, with split SDs up to 0.36 pp. More X17 evio files are expected next week, once disk space allows; the templates and X17 resolution will be redone then.

## 6. mf time offsets (failed)

`time_calib.json` comes from the C++ tool `analysis/tools/hycal_module_time_calib.cpp`. It is ported as `mf_validate.py timecal`. Input: 8 × 125k events of 025320 (r_model 0.3, κ 1.010), `output/mf_validate/timecal/`. Only 463 (rawsum) and 377 (mf) modules connect at > 100 entries per pair, against 1150 in the database.

Rawsum control vs the database: median −0.05 ns, RMS 0.11 ns, 92% within 0.2 ns (463 modules). The mf table `mf_time_calib_hybrid.json` has the 377 solved modules; the other 773 take the database offset + 0.266 ns.

Acceptance on independent runs (20k events each), seed–neighbour σ_MAD (ns):

| Run | rawsum + DB offsets | rawsum + ported offsets | mf + mf offsets | mf raw | rawsum raw |
|---|---|---|---|---|---|
| 24503 | 0.672 | 0.651 | 0.780 | 2.22 | 2.43 |
| 25257 | 0.709 | 0.691 | 0.810 | 2.26 | 2.44 |

mf is narrower on raw times and wider after per-module offsets, also with equal-statistics tables. Phase 4 therefore runs on raw times only, and the mf time-calib plumbing was removed from the harness. Why per-module offsets help rawsum more than mf is an open physics item.

## 7. Gate rescan with mf (Phase 4)

`output/gate_scan_20261022/`, 17:35–18:52. mf, all pulses, raw times, r_model 0.3, κ 1.010, X17 with the 26138 templates. Runs and halves as Phase 1. Reference: rawsum legacy `traw_w8_pmap` from `--baseline output/gate_scan_20261021`. 0 failures; determinism check passed (58 comparisons).

- **A** (W ∈ {8, 16, 100}): W16 −0.140 ± 0.016 pp wins; W8 −0.125 is in band; W100 (gate off) fails ε_peak.
- **B1** (LO ∈ {100, 120, 130}): all in band (−0.139 to −0.142); the tie-break (widest pre-window) keeps LO 100.
- **B2** (HI ∈ {190, 200, 220}): HI 200 wins; 220 fails ε_yield (0.967, Ta).
- **Working point `mf_traw_w16_p100-200`**: the same gate as the Phase 1 rawsum working point.

Stage C (half B, preferred fit; σ/E in %, Δ in pp ± SE, paired) against the rawsum working point `traw_w16_p100-200`:

| Run type (fit) | rawsum WP (§2 C3) | mf WP | Δ vs rawsum WP | Δ vs ref `traw_w8_pmap` | Δgate (mf) |
|---|---|---|---|---|---|
| 3.49 GeV H2 (CB) | 1.947 | 1.891 | −0.056 ± 0.013 | −0.055 ± 0.013 | −0.001 ± 0.000 |
| 2.24 GeV H2 (Gauss) | 2.266 | 2.213 | −0.053 ± 0.007 | −0.052 ± 0.008 | −0.003 ± 0.001 |
| 2.24 GeV C (CB) | 2.161 | 2.090 | −0.071 ± 0.012 | −0.070 ± 0.012 | −0.001 ± 0.001 |
| 0.73 GeV H2 (Gauss) | 4.568 | 4.258 | −0.311 ± 0.039 | −0.314 ± 0.040 | +0.004 ± 0.005 |
| 2.24 GeV Ta X17 (CB) | 2.345 | 2.287 | −0.057 ± 0.058 | −0.145 ± 0.054 | −0.099 ± 0.024 |
| Mean of 5 | | | | −0.127 ± 0.014 | −0.020 ± 0.005 |

- On H2 and C the gain is the pulse method: within mf, W = 16 moves σ/E by ≤ 0.004 pp.
- On Ta, mf vs the rawsum WP is −0.057 ± 0.058 (not significant on half B); the gate supplies most of the gain (Δgate −0.099 ± 0.024).
- ε on half B: ε_yield 0.9932–0.9990, ε_peak 1.0002–1.0151.
- The over-tuning guard (half A vs half B) is ok in every run type.

## 8. Interpretation

**Resolution model.** A rough three-point estimate read off the mf CB post values at 0.73, 2.24 and 3.49 GeV is σ/E ≈ 2.5%/√E ⊕ ~1.2%; the fit below checks it. It uses one fit model across run types, although 2.24 and 0.73 GeV H2 prefer Gauss. An unweighted least-squares fit of a/√E ⊕ b to the per-run mf CB legacy values in `resolution_mf_vs_rawsum_fullscan.md` (10 H2 runs) gives a ≈ 2.3%, b ≈ 1.27%.

**Stochastic vs constant.** At 3.49 GeV a/√E is 1.23–1.34%, comparable to b; the constant term is about half of σ² (45–52%) and its share grows with E.

**Implications.**

- The pulse method helps most at low energy: −0.31 pp at 0.73 GeV against −0.04 to −0.07 pp at 2.24 and 3.49 GeV (§4.1). This is where the stochastic and noise part dominates and where mf recovers the low-amplitude charge (§3).
- Timing is not a resolution lever on e-p or C, with either pulse method (§2, §7). Ta is the only run type where the gate matters.
- The next levers target the constant term: per-crystal calibration, a position-dependent correction, and seed-time-constrained neighbour amplitudes. `plan_20261023.md` is in preparation.

## 9. Caveats and open items

| Item | Detail |
|---|---|
| Lost Phase 2 outputs | `output/mf_validate/` was deleted, except `timecal/`; see ¹ in §3. |
| X17 statistics | See §5. |
| mf time offsets | Failed acceptance (§6); open. |
| Harness | Cosmetic bug in the interim `summarise` of `gate_scan_20261022.py`; the final outputs are correct. |
| Driver plot titles | mf titles end in nested parentheses, "(…, mf (template))": `title_suffix` appends `_PULSE_LABEL` (`benchmark_hycal_timing_scan.py:1139`). |
| Git | Nothing committed since `315e4af` (§1). |

## 10. How to run

```bash
cd /Users/L00338853/work/PRad/prad2evviewer/analysis/pyscripts && source setup_env.sh
E=/Users/L00338853/work/PRad/data/evio
# (a) single run, 025320 split 00400, mf at the working point (rawsum: drop --pulse-method mf)
python3 benchmark_hycal_timing_test_20261021.py $E/prad_025320.evio.00400 output/<tag>/example \
  --calibrate --exclude-rowcol 15 20 15 20 --pre-window 100 200 --window 16 --inner-radius 200 --min-per-module 30 \
  --hc-calib-file $PRAD2_DATABASE_DIR/calibration/calibration_factor_2p2_fadc.json --lms-corr auto --pulse-method mf
#     variants: --single-pulse {cluster,seed}; --mf-templates PATH [--mf-template-column RUN]; --pulse-dump F.npz --max-events N
# (b) mf full scan through the driver (as run for output/scan_20261022_mf)
python3 benchmark_hycal_timing_scan.py --pulse-method mf --window 16 --pre-window 100 200 \
  --output-dir output/scan_<tag> --parallel 12 --timeout-seconds 1800
# (c) harness (sets MPLBACKEND=Agg itself)
OUT=output/gate_scan_<tag>; mkdir -p $OUT/logs; BL="--baseline output/gate_scan_20261021"
G="python3 gate_scan_20261022.py"; P="--out $OUT $BL --parallel 12 --timeout-seconds 1800"
$G run --stage R $P && $G restable --out $OUT                                 # Phase 3, 11 runs
# §5 RX and §7 ran when `auto` picked the 26138 file; `auto` now picks the combined file, so those X17 numbers do not reproduce exactly. For both RX lines, add --bench-args "--mf-templates output/pulse_templates_x17_26138_h50.dat" (--bench-args applies to every run of a stage, so do not use it for A–C).
$G run --stage RX $P && $G restable --out $OUT --stage RX                     # X17 runs
$G run --stage RX --sp cluster $P && $G restable --out $OUT --stage RX --sp cluster
$G run --stage R --sp cluster $P && $G restable --out $OUT --sp cluster
$G run --stage R --sp seed --max-splits-per-run 2 $P && $G restable --out $OUT --sp seed
for S in A B1 B2; do $G run --stage $S --half A $P && $G select --stage $S --out $OUT $BL; done   # Phase 4
$G run --stage C --half B $P && $G summarise --out $OUT $BL
# (d) X17 templates (cuts as 25320), then the per-module table
python3 fit_pulse_template.py "$E/prad_025800.evio.*" "$E/prad_026138.evio.*" "$E/prad_026230.evio.*" \
  -o output/pulse_templates_x17_combined_SSP_RawSum_h50_max20modules.json \
  --model two_tau --chi2-max 3 --height-min 50 --height-rms-mult 10 --min-pulses 50 --t0-min 22 \
  --model-err-floor 0.03 --trigger-event-type SSP_RawSum --max-modules-per-event 20 --max-pulses-per-channel 0
python3 dump_pbwo4_table.py --file-a output/pulse_templates_x17_combined_SSP_RawSum_h50_max20modules.json \
  --file-b output/pulse_templates_025320_SSP_RawSum_h50_max20modules.json \
  --label-a '#27000' --label-b '#25320' --out output/pulse_templates_x17_combined_h50.dat
```

`--dry-run` on any `run` call prints the job count. `--mf-templates auto` picks the combined X17 file for 25800, 26138 and 26230.
