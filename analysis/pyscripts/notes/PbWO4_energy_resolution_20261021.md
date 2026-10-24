# PbWO4 energy resolution with LMS correction and elastic-window refinement

**Date:** 2026-10-21
**Branch:** `waveform_jz`
**Plan:** `analysis/pyscripts/plans/plan_benchmark_hycal_timing_test_20261021.md`
**LMS research:** `docs/research-lms-gain-tracking.md`

All paths are relative to the repo root `/Users/L00338853/work/PRad/prad2evviewer` unless absolute.

## 1. Status

- `benchmark_hycal_timing_test_20261021.py` is the benchmark now in use.
- `benchmark_hycal_timing_scan.py` calls it.
- `benchmark_hycal_timing.py` is kept unchanged as the reference.

## 2. What changed

### 2.1 Benchmark: `benchmark_hycal_timing_test_20261021.py`

| Item | Behaviour |
|---|---|
| LMS gain tracking | Reads the C++ `replay_gainCorr` output `prad_NNNNNN_gain_corr.root` via uproot (lazy import). Factor = alpha-normalised (LMS2+LMS3)/2, renormalised to the run's first batch (`--lms-ref first-batch`, default). `ref-run` keeps the stored `ref_run` scale, as the C++ replay does. |
| Where it applies | Multiplies every pulse integral before `energize()`, in both paths. Batch looked up by TI event number. PbWO4 only; PbGlass and Veto keep 1.0. |
| Invalid entries | Failed-fit sentinel 1.0 and missing reference 0 are handled. Uncorrected module-batches are counted. |
| `--lms-corr {auto,off,PATH}` | `auto` = `$PRAD2_DATABASE_DIR/gain_factor/gain_correction/prad_<run>_gain_corr.root`. A missing file exits 2; there is no silent fallback. |
| Elastic-window refinement | `--cal-method elastic-window`, the new default. `cb-peak` and `median` are kept unchanged. |
| Kinematics | E'(θ) is a port of `PhysicsTools::ExpectedEnergy` `"ep"`, including `EnergyLoss` (2.056·secθ MeV). θ = atan(r/hycal_z) from the seed-module centre, as in `physics_calib.cpp`. |
| `--target {H2,C,Ta}` | Default from the per-run map `RUN_TARGET` (below). C (12 u) and Ta (180.948 u) use plain elastic recoil E·M/(M+E(1−cosθ)) plus the same `EnergyLoss`, because the C++ has no nucleus mode. An unmapped run with elastic-window exits 2. |
| Møller exclusion | A module is not used if the elastic–Møller gap (mu_global/Ē')·(E'_el − E_ee) < 3σ_M (`MOLLER_GAP_NSIGMA = 3.0`). |
| Beam-energy scaling | kE = E_beam/2239.51, kS = sqrt(kE). Defaults: `--signal-min` 500·kE, `--e-range` 500–3500·kE, `--min-seed-energy` 200·kE, `--bin-width` 20·kS, `--target-peak` E_beam. Explicit CLI values override. Resolved values are printed and written to the `.tsv` header. |
| Calibration-mismatch check | The energy tag in the calibration file name (`0p7`, `2p1`, `2p2`, `3p5`, `x17`) must lie within 3% of E_beam. A mismatching runinfo default exits 2; a mismatching explicit `--hc-calib-file` only warns. |

Elastic-window refinement, per path:

| Step | Rule |
|---|---|
| 1. Global fit | CB fit of the pool gives mu_global and sigma_glob. |
| 2. Module pre-fit | CB on mu_global·[1−0.15/kS, 1+0.10/kS], lower edge raised to max(…, (1+E_ee/E')/2) (Møller-aware). |
| 3. Window | Fixed μ_M ± 2σ_M (`--window-nsigma 2`). |
| 4. Gain | g_M = E'(θ_M) / median(in-window). |
| 5. Passes | Single pass; a second pass returns g = 1 exactly. |
| Fallback | μ_M = median of the pre-fit range, σ_M = sigma_glob, when the module fit fails or n_prefit < 100. |
| Statistics cut | `--min-per-module` (30) counts in-window events. |

`RUN_TARGET`:

| Target | Runs |
|---|---|
| H2 | 24386, 24503, 24512, 24578, 24590, 24603, 25228, 25257, 25260, 25299, 25320 |
| C | 25308 |
| Ta | 25800, 26138, 26230 |

New and changed outputs:

| Output | Content |
|---|---|
| `<out>_lms.tsv` | Per-batch factor statistics. |
| `<out>_lms.png` | Per-batch factor plot. Not written with `--lms-corr off`. |
| `<out>_window.tsv` | Per path and module: θ, E', E_ee, prefit range, μ_M, σ_M, `sigma_source`, window, `n_in_window`, gain, `moller_gap_nsigma`, `exclude_reason`, `used`. |
| `<out>.tsv` | Comment header lines added only. Schema unchanged; `parse_tsv` still works. |
| `<out>_calibrated.png` | Also overlays the post-refinement Gaussian fit as a dashed curve. |

### 2.2 Scan driver: `benchmark_hycal_timing_scan.py`

| Item | Change |
|---|---|
| `_BENCHMARK_SCRIPT` | Points to `benchmark_hycal_timing_test_20261021.py`. |
| Flags dropped | `--cal-method cb-peak`, `--target-peak`, `--signal-min`, `--e-range`, `--bin-width`, `--min-seed-energy` are no longer passed; the benchmark derives them from E_beam. Per-run `--pre-window` is still passed. |
| `--lms-corr` | `auto` if the gain_corr file exists, else `off` (logged as a warning). |
| `_EXCLUDED_RUNS` | `{24578: "non-standard trigger"}` |
| Trend plots | Marker shape by run type (0.73 GeV H2, 2.24 GeV H2, 2.24 GeV C, 2.24 GeV Ta (X17), 3.49 GeV H2). Filled marker + dashed line = CB; open marker + solid line = Gaussian. One colour per run. |

`_HC_CALIB_MAP` (per-energy `_fadc` calibration, under `$PRAD2_DATABASE_DIR/calibration/`):

| Runs | File |
|---|---|
| 3.49 GeV | `calibration_factor_fadc_3p5.json` |
| 0.73 GeV | `calibration_factor_fadc_0p7.json` |
| 2.24 GeV H2 and C | `calibration_factor_2p2_fadc.json` |
| 25800 | `x17_25825_fadc.json` |
| 26138 | `x17_26094_fadc.json` |
| 26230 | `x17_26214_fadc.json` |

## 3. Validation

Outputs in `analysis/pyscripts/output/lmstest/`.

| Check | Result |
|---|---|
| V1: E' port vs compiled C++ `PhysicsTools` | max \|Δ\| 2.5e-4 MeV (`"ep"`). |
| V3: regression | `--lms-corr off --cal-method cb-peak` reproduces `benchmark_hycal_timing.py` bit-identically on 50k events of 025320 split 00400. Non-`#` lines of `r4b_old_00400_50k.tsv` and `r4b_newoff_00400_50k.tsv` are identical. |
| V4: batch-0 identity | Split 00050 sits in batch 0 (`batches_touched=0`). On 50k events, LMS on and off give identical non-`#` lines (`r4_v4_on` vs `r4_v4_off`). |
| Auto defaults | `r4_on_8s_auto` reproduces the explicit-flag run exactly. |

8-split 025320 (splits 00050–00400, batches 0, 1, 2, 4, 5, 6, 7, 8 touched; 264 modules). Post σ_E/E (%):

| Run | Legacy Gauss | Legacy CB | New Gauss | New CB |
|---|---|---|---|---|
| old script, cb-peak (`old_cbpk_8s`) | 2.28 | 2.07 | 2.32 | 2.08 |
| new script, elastic-window, LMS on (`r4_on_8s_explicit`) | 2.28 | 2.07 | 2.31 | 2.09 |
| new script, elastic-window, LMS off (`r4_off_8s_explicit`) | 2.28 | 2.07 | 2.31 | 2.09 |

- Elastic-window and cb-peak agree within 0.01 pp after refinement.
- Post σ_E/E is identical on and off. CB pre σ_E/E is 2.97 vs 2.96% (Legacy) and 3.02 vs 3.01% (New), on vs off.
- The New-path Gaussian pre σ is 84.62 MeV on vs 81.97 MeV off (Legacy 80.95 vs 80.91 MeV).
- In-run LMS drift is ≤ 1.2% relative to batch 0 (factor range 0.988–1.010).

## 4. Findings

### 4.1 Per-run σ_E/E, full scan `output/scan_20261021_full/`

Post-refinement σ_E/E (%), elastic-window method, inner ring r < 200 mm. Values are means over splits; ± is the split-to-split standard deviation (legacy only). **Legacy** is the largest in-window peak per module (per-run pre-window; see §5). **New** is the gated clustering path. Δ = New − Legacy, in percentage points (pp), using the preferred fit model (marked ★; see the Fit-model note). Run 24578 is excluded.

| Run | E (GeV) | Target | Splits | Modules | CB post: Legacy | CB post: New | Gauss post: Legacy | Gauss post: New | Δ pref. (pp) |
|---|---|---|---|---|---|---|---|---|---|
| 24386 | 3.49 | H2 | 10 | 66 | ★1.93 ± 0.02 | ★1.97 | 1.96 ± 0.02 | 2.00 | +0.04 |
| 24503 | 3.49 | H2 | 10 | 151 | ★1.94 ± 0.02 | ★1.98 | 1.98 ± 0.02 | 2.02 | +0.04 |
| 24512 | 3.49 | H2 | 7 | 143 | ★1.94 ± 0.01 | ★1.97 | 1.98 ± 0.01 | 2.01 | +0.03 |
| 24590 | 0.73 | H2 | 10 | 108 | 3.39 ± 0.07 | 3.40 | ★4.60 ± 0.12 | ★4.62 | +0.02 |
| 24603¹ | 0.73 | H2 | 6 | 93/94 | 3.39 ± 0.06 | 3.42 | ★4.75 ± 0.20 | ★4.79 | +0.04 |
| 25228 | 2.24 | H2 | 11 | 126 | 2.03 ± 0.02 | 2.04 | ★2.26 ± 0.02 | ★2.29 | +0.03 |
| 25257 | 2.24 | H2 | 7 | 130 | 2.04 ± 0.02 | 2.05 | ★2.27 ± 0.02 | ★2.29 | +0.02 |
| 25260 | 2.24 | H2 | 8 | 130 | 2.02 ± 0.02 | 2.03 | ★2.27 ± 0.02 | ★2.30 | +0.03 |
| 25299 | 2.24 | H2 | 8 | 130 | 2.02 ± 0.03 | 2.03 | ★2.26 ± 0.02 | ★2.29 | +0.03 |
| 25320² | 2.24 | H2 | 8 | 129 | 2.04 ± 0.03 | 2.05 | ★2.27 ± 0.01 | ★2.30 | +0.03 |
| 25308² | 2.24 | C | 8 | 173 | ★2.17 ± 0.01 | ★2.19 | 2.24 ± 0.01 | 2.27 | +0.02 |
| 25800 | 2.24 | Ta (X17) | 8 | 34 | ★2.34 ± 0.06 | ★2.32 | 2.39 ± 0.04 | 2.38 | −0.02 |
| 26138 | 2.24 | Ta (X17) | 4 | 18 | ★2.43 ± 0.10 | ★2.34 | 2.49 ± 0.09 | 2.48 | −0.09 |
| 26230 | 2.24 | Ta (X17) | 13 | 10/9 | ★2.48 ± 0.23 | ★2.34 | 2.54 ± 0.18 | 2.42 | −0.14 |

¹ 24603 re-run on 2026-10-21 with the 130–200 ns pre-window (as for 24590); the change from the default-window scan is ≤ 0.01 pp.
² LMS gain correction was applied only to 25308 and 25320 (the only runs with a `gain_corr` file). All other runs ran with LMS off.

### 4.2 Fit model

| Run type | Preferred fit |
|---|---|
| 3.49 GeV e-p H2 (24386, 24503, 24512) | CB |
| 2.24 GeV carbon (25308) | CB |
| 2.24 GeV X17 Ta (25800, 26138, 26230) | CB |
| 2.24 GeV e-p H2 (25228, 25257, 25260, 25299, 25320) | Gaussian |
| 0.73 GeV e-p H2 (24590, 24603) | Gaussian |

Quote σ_E/E from the preferred fit model for each run type.

### 4.3 Observations

- **Calibration offset.** With `calibration_factor_2p2_fadc.json`, μ_global sits ≈ +7.5% above Ē' for H2 (025320 8-split: +7.69% Legacy, +7.41% New, LMS on) and C (25308 splits: Legacy +7.36 to +7.62%, New +7.08 to +7.34%). The refinement absorbs it.
- **3.49 GeV pre vs post.** Pre-refinement CB σ_E/E (split means from `benchmark_hycal_timing_scan_summary.tsv`, 15:27 version) is ~2.3% for 24386 and ~4.2% for 24503 and 24512. All three converge to ~1.94% post.
- **0.73 GeV.** The refinement barely changes σ_E/E. The Gaussian value (4.6–4.8%) is driven by the radiative tail.
- **X17 Ta.** New is better than Legacy, on weak statistics; see §5.

## 5. Legacy vs New clustering

Legacy takes the largest-integral in-window pulse per module (per-run `--pre-window`: 130–200 ns for 24xxx, 100–200 ns for 25xxx and 26xxx; 24603 see ¹), with `seed_time_window` ≤ 0 and one hit per module. New is gated clustering: every in-window pulse becomes a separate hit, and `HyCalCluster` applies a seed-anchored time-coincidence cut W (`seed_time_window`, `--window`, default 8 ns) during clustering.

| Sample | Δ = New − Legacy (pp) | Reading |
|---|---|---|
| e-p (all energies) and C | +0.02 to +0.04 | New worse, consistently across splits. |
| X17 Ta | −0.02 to −0.14 | New better, on 9–34 modules with large split spreads. Suggestive. |

### 5.1 Old vs new gain refinement

The old method is described in `docs/technical_notes/hycal_clustering/hycal_clustering.md`, section "Calibration-equalised resolution" (line ~375).

| | Old (`cb-peak`) | New (`elastic-window`) |
|---|---|---|
| Input spectra | No LMS | LMS-corrected |
| Per-module reference | CB μ | Median of events in μ_M ± 2σ_M |
| Target | Fixed `--target-peak` (3488.5 MeV in the technical note) | E'(θ_M) of the run's target |
| Passes | 2 iterations | Single pass |
| Møller handling | None | Møller-aware pre-fit edge; exclusion if gap < 3σ_M |

## 6. Caveats and open items

| Item | Detail |
|---|---|
| LMS coverage | See ² in §4.1. The X17 runs need `replay_gainCorr`. |
| 3-minute LMS batches | The `-b 1600` files are not yet made. The existing 025320 file uses `-b 4000` (about 7.5 min, 9 batches). Regeneration needs all 400+ splits of the run (plan §8). |
| Re-runs | 24603, 25228 and 25257 were re-run in full (24 splits, 2026-10-21); 25228 and 25257 are unchanged, and every other run is byte-identical. |
| 24386 split 00060 | Timed out (1800 s limit); excluded. |
| 24578 | Excluded: non-standard trigger, ~2000 physics-trigger events per split, no per-module refinement. |
| X17 calibration | Taken from the X17 `_fadc` files; unverified. |
| Git | Nothing committed. `benchmark_hycal_timing_test_20261021.py` is untracked; V3 used the working-tree `benchmark_hycal_timing.py`. |

## 7. How to run

```bash
cd /Users/L00338853/work/PRad/prad2evviewer/analysis/pyscripts && source setup_env.sh
# full scan (settings as run for output/scan_20261021_full; use a new output dir)
python3 benchmark_hycal_timing_scan.py --output-dir output/scan_<tag> --parallel 12 --timeout-seconds 1800
# single run, 025320 split 00400, LMS on, elastic-window, energy defaults from E_beam
python3 benchmark_hycal_timing_test_20261021.py /Users/L00338853/work/PRad/data/evio/prad_025320.evio.00400 output/lmstest/example \
  --calibrate --exclude-rowcol 15 20 15 20 --pre-window 100 200 --inner-radius 200 --min-per-module 30 \
  --hc-calib-file $PRAD2_DATABASE_DIR/calibration/calibration_factor_2p2_fadc.json --lms-corr auto
# regression path (old method, no LMS): add  --cal-method cb-peak --lms-corr off
```

A quoted glob (`".../prad_025320.evio.*"`) in place of the single file processes several splits in one job.
