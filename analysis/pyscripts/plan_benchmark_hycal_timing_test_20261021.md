# Plan: `benchmark_hycal_timing_test_20261021.py`, LMS-corrected benchmark with elastic-window gain refinement

**Status:** draft, awaiting user approval
**Date:** 2026-10-21
**Branch:** `waveform_jz`
**Research input:** `docs/research-lms-gain-tracking.md` (spot-checked against the code and the 025320 `gain_corr.root` on 2026-10-21)

All paths are relative to the repo root `/Users/L00338853/work/PRad/prad2evviewer` unless absolute.

---

## 1. Goal

Create `analysis/pyscripts/benchmark_hycal_timing_test_20261021.py`, a new version of `benchmark_hycal_timing.py` (the working-tree copy, including its uncommitted pre-calibration CB changes). It adds two things:

1. **LMS gain-drift correction.** It reads the per-batch, per-PbWO4-module LMS factors that the C++ `replay_gainCorr` writes to `prad_NNNNNN_gain_corr.root`. It renormalises them to the run's first batch and multiplies every HyCal pulse integral by its module's factor before `energize()`. The batch is looked up by TI event number.
2. **Elastic-window gain refinement.** For each seed module and each path (legacy and gated), it fits a Crystal-Ball to the module's LMS-corrected, pre-refinement cluster-energy spectrum. It keeps the events within ±2σ of the fitted peak and sets the module gain so that the median of those events lands on the expected elastic e-p energy E'(θ_M).

`benchmark_hycal_timing.py` stays untouched.

## 2. Approach and rationale

| Choice | Rationale | Rejected alternative (one line each) |
|---|---|---|
| Read the existing C++ `replay_gainCorr` output via uproot (user decision 2) | It already implements alpha-normalised reference-PMT LMS gains per batch with the bit-24 AND `nch>1000` selection. | Reimplementing LMS in Python: duplicates validated C++, and is out of scope. |
| Run `replay_gainCorr -b 1600` (about 3 min at the measured 8.9 Hz LMS rate) | Matches the user's "every 3 mins". | Fixed-time binning: the C++ bins by LMS count. 1600 events correspond to 3.0 min at 8.9 Hz. |
| Factor relative to the run's first batch, `A_b / A_0` with `A = (c[LMS2]+c[LMS3])/2` (decision 2) | Pure in-run drift. It leaves the absolute scale of `calibration_factor_2p2_fadc.json` unchanged. | The file's native `ref_run 24827` scale: it moves the 2p2 calibration by about 2% on average and up to 12% per module. It stays available as `--lms-ref ref-run`. |
| Apply to `pk.integral` before `energize()`, per module, for both paths (mirrors `analysis/src/Replay.cpp:492,527`) | Clustering then sees drift-corrected hits. This is the same insertion point as the replay. | Scaling the cluster energy afterwards: mixes modules inside a cluster incorrectly. |
| CB pre-fit seeded at E'(θ_M) on `[E'·(1−0.15), E'·(1+0.10)]`, **with a Moller-aware floor** (see §5.4) | Keeps the fit on the elastic peak. | The plain −15% edge: for modules at r ≈ 53 mm the Moller peak (−13.6%) sits inside it. |
| Window `[μ_M − 2σ_M, μ_M + 2σ_M]`, fixed from the pre-fit (decision 3); peak estimator is the median | The CB tail is unconstrained inside a ±2σ cut, and the median is robust there. | Re-centring the window each iteration: rejected by the user. |
| **Single pass, no iteration** | With the window fixed in pre-refinement energy, the selected event set is fixed. `g_M = T_M / median` then reproduces itself exactly on a second pass, so iterating is a no-op. | `n_iter = 2` as in the old `fit_per_module_gains`: redundant here. |
| Keep the old `--cal-method median|cb-peak` paths unchanged | They give an exact regression test against the old script, and the later scan-script swap passes `--cal-method cb-peak`. | Removing them: loses the bit-exact regression check. |

## 3. Verified facts the plan relies on

| Fact | Evidence |
|---|---|
| `replay_gainCorr` binary is built: `build/bin/prad2ana_replay_gainCorr`. Usage: `<evio_file_or_dir> [...] -o dir [-b batch] [-r ref_run] [-j threads] ...`. `-b` defaults to 4000. | `analysis/tools/replay_gainCorr.cpp:~215,~247` |
| The output path is **hardcoded** to `<PRAD2_DATABASE_DIR>/gain_factor/gain_correction/prad_<run>_gain_corr.root` (RECREATE). There is no output-path flag. | `replay_gainCorr.cpp:~312-326` |
| A directory argument globs every file whose name contains `.evio`, so mixed runs are possible. Pass the explicit `prad_025320.evio.*` list instead. | `replay_gainCorr.cpp:55-70` |
| Tree `gain_corr` branches: `batch_id, event_num_start, event_num_end, n_lms_events, n_alpha_events, ref_run` (int32), `unix_time` (uint32), `gain_corr_W[1156][3]` (float), plus fit means. Row `i` is module id `1001+i` (W(i+1)). Columns are LMS1, LMS2, LMS3. | uproot inspection of `database/gain_factor/gain_correction/prad_025320_gain_corr.root`; `analysis/include/gain_factor.h` |
| Failed fits are encoded as `gain_corr_W = 1.0` exactly (45 entries per column in 025320). A missing reference is encoded as `0` (W1141 in every batch). The C++ loader maps 0 to 1.0 per column. | `GainCorrCompute.cpp` `FlushGainBatch`; `gain_factor.h` `LoadGainCorrTimeSeries` |
| Batch lookup: the last batch with `event_num_start <= ev`. Events before batch 0 use batch 0, events after the last batch use the last batch. The trailing partial batch is dropped. | `gain_factor.h:206-226` (`GetCorr`); `GainCorrCompute.cpp` (no flush after the loop) |
| The batch key is the CODA/TI event number, `info.event_number`. It increases monotonically across the whole run: 1222 → 70,146,636 in the 025320 file. | `Replay.cpp:1925`; `prad2dec/src/EvChannel.cpp:442-449` |
| Python exposes `fadc_evt.info.event_number` (int32, same source). | `python/bind_dec.cpp:146`; `prad2dec/include/Fadc250Data.h:85` |
| Existing file for 025320: 9 batches of 4000 LMS events (about 7.5 min each), `ref_run 24827`. Relative to batch 0, `A_b/A_0` over 1155 modules spans 0.988–1.010, with per-batch median 0.999–1.000. Five rows have both LMS2 and LMS3 at the 1.0 sentinel. Four of them (W561, W562, W595, W596) are beam-hole ids with no module. Only W1124 is a real failed fit. | uproot inspection (scratch); `database/hycal_map.json` |
| One split is about 160k events, about 7.7 s of beam at about 20.8 kHz. A 1600-LMS batch covers about 3.75M event numbers, about 23 splits. | `output/baseline_025320_split00010_cbpeak.tsv` header (`n_phys=160765`); batch event ranges |
| Locally present: only `prad_025320.evio.{00050,00100,…,00400}` (8 splits, 2.1 GB each) in `/Users/L00338853/work/PRad/data/evio/`. The run has 400 or more splits. | `ls` |
| All 296 module centres with r < 200 mm are PbWO4. The innermost PbGlass centre is at r = 371.9 mm. W1141 (id 2141, index 1712) is at r = 343.8 mm, outside the inner ring. | prad2py geometry scan (scratch) |
| `row`/`column` are 0-based. The beam hole is (16,16), (16,17), (17,16), (17,17). `--exclude-rowcol 15 20 15 20` (the scan default) is asymmetric: modules at r ≈ 53 mm in row 14 and col 14 remain. | prad2py geometry scan |
| For run 25320, `p.geo.Ebeam = 2239.51` MeV, `p.geo.hycal_z = 6267.0` mm, and `hycal_x/y = 0`. Runinfo selects the **3.5 GeV** calibration, so `--hc-calib-file …/calibration_factor_2p2_fadc.json` must be passed. | `setup_pipeline(run_num=25320)`; `Commands.md:167`; `database/runinfo/general.json:76-77,26-28` |
| `PhysicsTools::ExpectedEnergy` `"ep"`: `E·M_p/(M_p+E(1−cosθ)) − EnergyLoss`, with `EnergyLoss = 2.056·secθ` MeV. `"ee"` (Moller): `m_e·((γ+1)+(γ−1)cos²θ)/((γ+1)−(γ−1)cos²θ) − EnergyLoss`, with `γ = E/m_e`. Constants: `M_PROTON=938.272`, `M_ELECTRON=0.51099895`, `DEG2RAD=3.14159265/180` (float). | `analysis/src/PhysicsTools.cpp:23-25,284-332` |
| `physics_calib.cpp` computes θ as `atan(hypot(mod.x, mod.y)/hycal_z)` on the module centre, with no shower depth. | `analysis/tools/physics_calib.cpp:~571-576` |
| Module bindings: `id, index, name, x, y, row, column, type, is_pwo4(), is_glass(), is_hycal(), energize()`, where `energize(adc) = cal_factor·adc`. `hycal.module_count()` and `hycal.module(i)` exist. `r.cluster.hits` is a list of `ModuleHit` (`index, energy, time`). | `python/bind_det.cpp:531-556,597-614,736-771` |
| uproot is **not installed** in the miniforge Python 3.13 that runs prad2py. PyROOT is also absent. | `python3 -c 'import uproot'` fails; `build/python/prad2py.cpython-313-darwin.so` |
| `analysis/pyscripts/output/` is gitignored. `database/gain_factor/gain_correction/*.root` is **git-tracked**. | `.gitignore:38`; `git ls-files` |

### 3.1 Consequence for what the LMS correction can show

- A single split (about 8 s) sits inside one LMS batch, so on one split the factor is a per-module constant. The per-module refinement absorbs a constant per-module scale for the seed module. On one split, the LMS correction can therefore change only the **pre-refinement** spectrum and the cross-module sharing inside clusters. The post-refinement σ_E/E should not move.
- The LMS correction becomes visible only on samples that span several batches, such as the 8 local splits 00050–00400. They touch batches {0,1,2,4,5,6,7,8} of the existing file; no local split falls in batch 3.
- The in-run drift is small: at most 1.2% per module, with a median near 0. The expected improvement in pre σ_E/E on the 8-split sample is therefore small, of order 0.1 pp or less.
- Validation (§9) is built around this. A null result on one split is the expected outcome, not a failure.

## 4. Files

| File | Action |
|---|---|
| `analysis/pyscripts/benchmark_hycal_timing_test_20261021.py` | **Create.** Start as a copy of the working-tree `benchmark_hycal_timing.py`, then apply §5. |
| `analysis/pyscripts/output/lms/prad_025320_gain_corr_b1600.root` | Produced by the `replay_gainCorr` step (§8). Gitignored location. |
| `analysis/pyscripts/output/lmstest/*` | Validation outputs (gitignored). |
| `benchmark_hycal_timing.py`, `benchmark_hycal_timing_scan.py`, C++ sources, technical note, `database/**` | **Not modified.** |

## 5. Design

### 5.1 CLI (new flags; all existing flags kept with unchanged defaults and meaning; energy-scale defaults superseded by §11)

| Flag | Default | Meaning |
|---|---|---|
| `--lms-corr {auto,off,PATH}` | `auto` | `auto`: `$PRAD2_DATABASE_DIR/gain_factor/gain_correction/prad_<run:06d>_gain_corr.root`. Without the env var, use `<repo>/database/…`. The run number comes from `--run-num` if it is > 0, else `C.extract_run_number(evio_path)`. `off`: all factors are 1.0 and uproot is not imported. `PATH`: an explicit file. |
| `--lms-ref {first-batch,ref-run}` | `first-batch` | `first-batch`: `A_b/A_0` with the validity rules of §5.2. `ref-run`: `A_b` as stored, with C++ loader semantics (0 → 1.0 per column, mean of LMS2 and LMS3). This matches the C++ replay. |
| `--cal-method {elastic-window,cb-peak,median}` | **`elastic-window`** | The new method (§5.4), with target `T_M = E'(θ_M)` always. `cb-peak` and `median` are the unchanged old refinement, which uses `--target-peak` and 2 iterations. `--target-peak` is ignored by elastic-window. |
| `--window-nsigma` | `2.0` | Half-width of the refinement window, in units of σ_M. |

- The E' method is not a flag. Only the PhysicsTools `"ep"` port was decided, and choosing `elastic-window` selects it.
- Two tunables are module-level constants, not flags:
  - `PREFIT_RANGE = (-0.15, 0.10)`: the per-module CB pre-fit range as fractions of E'(θ_M), before the Moller floor.
  - `MIN_FIT_EVENTS = 100`: the minimum number of events in the pre-fit range for a per-module CB fit.

Behaviour on errors:

- If `--lms-corr auto` or a PATH points to a missing file, print the expected path and the §8 command, suggest `--lms-corr off`, and `sys.exit(2)`. **There is no silent fallback**, because an "LMS-corrected" output that is secretly uncorrected would be indistinguishable from a real one.
- If uproot is missing and `--lms-corr` is not `off`, `sys.exit(2)` with the install hint.
- The existing check that `--write-cluster-tsv` requires `--calibrate` stays.

### 5.2 LMS factor table

New function `load_lms_gain_corr(path, hycal, ref_mode) -> LmsTable`, using a small dataclass or dict:

- Read with `uproot.open(path)["gain_corr"]` and `.array(library="np")` for `batch_id, event_num_start, event_num_end, n_lms_events, unix_time, ref_run, gain_corr_W`. The template is `scripts/online_gain_monitor.py:213-256`.
- Sort by `event_num_start` and assert that it is strictly increasing. The reference batch is the first row after sorting.
- Build `factors`, of shape `(n_batches, hycal.module_count())`, initialised to 1.0. It is indexed by `mod.index`, because that is what the event loop has. For each row `i` in 0..1155:
  - Take `mod = hycal.module_by_id(1001+i)`, and skip `None` (the hole positions).
  - Let `c_bj = gain_corr_W[b, i, j]` for j = 1 (LMS2) and j = 2 (LMS3).
  - `valid(c) := c > 0 and c != 1.0`. 1.0 is the failed-fit sentinel. A genuine value of exactly 1.0f from a ratio of fitted floats is not expected.
  - **first-batch mode**:
    - If both columns are valid in batch b and in batch 0: `f = (c_b1 + c_b2) / (c_01 + c_02)`. This is the literal `A_b/A_0`.
    - Else, if exactly one column j is valid in both batches: `f = c_bj / c_0j`.
    - Else: `f = 1.0`, and the module-batch is counted as "uncorrected".
  - **ref-run mode**: map 0 to 1.0 per column, then `f = (c_b1 + c_b2)/2`. This is exactly `Replay.cpp:492,527`.
- PbGlass, Veto and LMS channels keep 1.0. They have no rows. W1141 (all zeros) gets 1.0 in first-batch mode by the rule above, and 1.0 in ref-run mode by the 0 → 1 mapping.
- `batch_of(ev) = clip(searchsorted(event_num_start, ev, side="right") − 1, 0, n_batches−1)`. This mirrors `GetCorr`.
- Keep `factors_list = [row.tolist() for row in factors]`. Scalar indexing into a Python list is faster than into numpy in the per-channel loop.
- `--lms-corr off`: a 1-batch table of ones, with `path = "off"`.

### 5.3 Event-loop changes

These go in the `main()` loop, at `benchmark_hycal_timing.py:367-445` in the copy.

1. After the trigger check:
   - Set `ev_num = int(fadc_evt.info.event_number)`, `b = lms.batch_of(ev_num)` and `fvec = lms.factors_list[b]`.
   - Increment `n_ev_per_batch[b]`, counting triggered events that reach clustering.
2. In the channel loop, compute `f = fvec[mod.index]` once per channel. Use `mod.energize(pk.integral * f)` in **both** the `cl_new.add_hit` call and the `cl_legacy.add_hit` call.
   - The selection of `best_pk` stays on the raw `pk.integral`. f is a positive per-module constant within an event, so the order is unchanged.
3. For each accepted single cluster (legacy and gated), compute `f_eff = Σ_h h.energy / Σ_h (h.energy / fvec[h.index])` over `r.cluster.hits`.
   - Compute it before the next `clear()`, because the cluster is invalidated then.
   - Use 1.0 if the sum is 0.
   - `E_c / f_eff` is the approximate no-LMS cluster energy. It is used only for the §5.5 diagnostic.
4. Extend the per-event tuple to `(seed_id, x, y, row, col, seed_E, cluster_E, f_eff)`. `ev_num` is used only for the batch lookup and count; it is not stored.
   - **Trap:** the old code unpacks with `for sid, *_, ce in pool`, which would bind `ce = f_eff` after the extension.
   - Update every unpacking to the explicit form `sid, x, y, row, col, se, ce, *_`. There are such unpackings in `fit_per_module_gains`, in its `raw`/`corr` builders (`e[6]` stays correct), and in the cluster-TSV writer.

### 5.4 Elastic-window refinement

There are four new functions. The kinematics ports go next to `gauss` and `crystal_ball`; the refinement functions go next to `fit_per_module_gains`.

**Kinematics ports.** These are double-precision ports of `PhysicsTools.cpp:284-332`.

```python
M_PROTON, M_ELECTRON, DEG2RAD = 938.272, 0.51099895, 3.14159265 / 180.0
def energy_loss(theta_deg):           # PhysicsTools::EnergyLoss
    c = math.cos(theta_deg * DEG2RAD)
    sec = 1.0 / c if c > 0.01 else 100.0
    return (0.5*1.6 + 0.12*1.6 + 0.1*2.0 + 0.48*1.8) * sec
def expected_energy_ep(theta_deg, ebeam):   # ExpectedEnergy(..., "ep")
def expected_energy_ee(theta_deg, ebeam):   # ExpectedEnergy(..., "ee"); 0 if den <= 0
def module_theta_deg(x, y, hycal_z):        # atan(hypot(x,y)/hycal_z), as physics_calib.cpp
```

The kinematics use `p.geo.Ebeam` and `p.geo.hycal_z`. θ comes from the seed-module centre, with no shower depth and no beam offset, matching `physics_calib.cpp`.

**`fit_module_cb(vals, lo, hi, seed, bin_width, sigma_glob)`** returns `(mu, sigma, alpha, n, ok)`.

- It histograms `vals` in `[lo, hi]` with `bin_width` (the `--bin-width`, 20 MeV).
- It runs `curve_fit(crystal_ball, …)` with `p0 = [max, seed, 50, 1.5, 2.0]` and bounds `([0, lo, 5, 0.1, 1.01], [inf, hi, hi−lo, 20, 50])`. This is the same pattern as the CB block of `per_module_peak`.
- `ok` requires:
  - no exception,
  - finite parameters,
  - `lo < mu < hi`,
  - `0.25·sigma_glob ≤ sigma ≤ 2·sigma_glob`.

**`fit_elastic_window_gains(events, *, ebeam, hycal_z, inner_r, min_seed_E, min_per_mod, exclude_rowcol, nsigma, bin_width)`** returns `(gains, raw, corr, n_modules, rows, sigma_glob, eprime_mean)`. It uses the module constants `PREFIT_RANGE` and `MIN_FIT_EVENTS`.

1. **Pool.** Apply the same `keep()` filter as `fit_per_module_gains`: seed r < `inner_r`, seed E ≥ `min_seed_E`, and `exclude_rowcol`. The single-cluster > `--signal-min` cut is already applied in the loop. Factor `keep()` into a shared `_select_pool()` that both functions use.
2. **Global σ.** Fit the whole pool with `fit_peak` and then `fit_peak_cb`, on `[Ē' − 800, Ē' + 800]` with the Gaussian μ as `mu0_hint`. Ē' is the mean of E'_M over the pool's modules. This mirrors the existing pre-CB call with `target_peak` replaced by Ē'.
   - `sigma_glob` is the CB σ, falling back to the Gaussian σ. If both are NaN, print an error and skip the refinement for this path; the script still writes the headline outputs.
   - `mu_global` is the CB μ of the same fit, falling back to the Gaussian μ, computed per path (user decision, 2026-10-21; replaces the original E'-centring). The global offset `mu_global/Ē' − 1` is reported per path. On 025320 with the 2p2 calibration it is about +7.7% (legacy) and +7.4% (gated), not the ~2% assumed when this plan was written.
3. **Per module** (seeds with any pool events):
   - Compute `θ_M`, `E'_M = expected_energy_ep` (superseded by §11: `expected_energy_elastic` for the run's target), and `E_ee,M = expected_energy_ee`.
   - **Pre-fit range** (in pre-refinement, LMS-corrected energy), centred on `mu_global`: `lo = mu_global·max(1+LO, (1 + E_ee,M/E'_M)/2)` and `hi = mu_global·(1+HI)`, with `(LO, HI) = PREFIT_RANGE`. The relative shape and the Moller-aware floor are those of the original E'-centred design, now measured from `mu_global`.
     - The Moller floor binds wherever E_ee > 0.70·E', which is r ≲ 88 mm. The Moller peak lies inside the plain −15% edge only for r ≲ 56 mm.
     - The MeV values below are for the original E'-centred range; with `mu_global` centring they scale by `mu_global/Ē'`. At r = 53 mm the floor moves `lo` from 1902 to 2086 MeV. At r = 75 mm it moves `lo` from 1902 to 1970, which is harmless (about −6σ from the elastic peak). At r = 33 mm (where E_ee = 2109) it moves `lo` to 2173.
   - If `n_prefit ≥ MIN_FIT_EVENTS`, run `fit_module_cb(vals, lo, hi, seed=mu_global, …)` on the full module sample; the histogram edges restrict the range (critic C2).
     - If `ok`: `μ_M, σ_M` come from the fit, with `sigma_source = "module"`.
     - Otherwise (or when `n_prefit` is too small): `μ_M = median(vals in [lo, hi])`, `σ_M = sigma_glob`, `sigma_source = "global"`.
     - The centre fallback completes the supervisor's σ-fallback default, because a failed fit leaves no μ either.
     - If `n_prefit == 0`, the module is not used.
   - **Window:** `[w_lo, w_hi] = [μ_M − nsigma·σ_M, μ_M + nsigma·σ_M]`, fixed. It is not clipped (decision 3). If `w_lo < lo`, set `flag_moller_overlap` (see Open Question 4). (Superseded by §11: replaced by the Moller-gap exclusion.)
   - Set `in_win = events with w_lo ≤ E_c ≤ w_hi`. The module is used if `len(in_win) ≥ min_per_mod`: `--min-per-module` applies to the in-window sample, which is the refinement sample.
   - Set `T_M = E'_M` and `g_M = T_M / median(in_win)`.
4. **Single pass.**
   - There is no loop. A second pass on `g_M·E_c` over the same fixed event set returns `g' = 1` exactly, so iterating cannot change anything.
   - This is stated in the docstring and asserted in V5 (§9).
5. **Return values.**
   - `raw` and `corr` hold **all** pool events of used modules, uncorrected and multiplied by `g_M`. This is the same semantics as the old function, so the pre and post σ_E/E stay comparable with old outputs.
   - `rows` holds the per-module table (§5.5).

**Wiring in the `--calibrate` branch.**

- Dispatch on `--cal-method`. `elastic-window` calls the new function once per path (`events_legacy`, `events_new`). The others call the unchanged `fit_per_module_gains`.
- For elastic-window, the calibrated-spectrum fit range is `cal_lo, cal_hi = Ē' ± 800` and the `_calibrated.png` reference line is at `Ē'`.
- For the old methods these stay at `--target-peak ± 800`, so they are unchanged.

### 5.5 Outputs

The existing outputs are kept, and their schemas are unchanged so that `benchmark_hycal_timing_scan.parse_tsv` keeps working:

- `<out>.tsv`: add **comment lines only** to the header:
  - `# lms_corr=<path|off> lms_ref=<mode> batches_touched=<ids> n_uncorrected_module_batches=<n>`
  - `# cal_method=<m> target_mean_MeV=<Ē' or --target-peak> window_nsigma=<n>`
  - The calibration-block comment line shows the same target value.
- `<out>.png`: unchanged, apart from an `LMS: on/off` tag in the title.
- `<out>_calibrated.png`: unchanged, apart from the target line and the title tag.
- `<out>_clusters.tsv`: same columns. `cluster_E_cal_MeV` uses the chosen method's gains. No new columns.

New outputs:

| Output | Content |
|---|---|
| `<out>_lms.tsv` | Header: file path, `ref_mode`, `ref_run`, `n_batches`, `n_lms_events` of batch 0, event range, uncorrected counts. One row per batch: `batch_id, event_num_start, event_num_end, unix_time, minutes_since_batch0, n_lms_events, n_events_this_job, f_median, f_p05, f_p95, f_min, f_max`. The `f_*` columns are taken over PbWO4 modules with at least one valid column in that batch. The file is written even with `--lms-corr off`, as one row of 1.0. |
| `<out>_lms.png` | The per-batch LMS factor against minutes since batch 0: the median line, the 5–95% band, and dashed min/max. Batches with `n_events_this_job > 0` are shaded. The title gives the file and the mode. It is skipped when `--lms-corr off`. |
| `<out>_window.tsv` (only with `--calibrate --cal-method elastic-window`) | One row per (path, seed module): `path, module, module_id, row, col, r_mm, theta_deg, eprime_MeV, emoller_MeV, prefit_lo_MeV, prefit_hi_MeV, n_total, n_prefit, mu_M_MeV, sigma_M_MeV, sigma_source, win_lo_MeV, win_hi_MeV, n_in_window, frac_in_window_lms, frac_in_window_nolms, median_in_window_MeV, target_MeV, gain, flag_mu_off, flag_moller_overlap, used`. Superseded by §11: `flag_moller_overlap` is replaced by `moller_gap_nsigma` and `exclude_reason`. |
| stdout summary | An LMS block: the file, mode, batches touched with event counts, the factor range over touched batches, and the uncorrected modules (first 10 names). A window block per path: modules used, module fits against global fallback, `flag_mu_off` and `flag_moller_overlap` counts (superseded by §11: Moller-excluded count), and the pooled in-window fraction with and without LMS. |

Diagnostic definitions:

- `frac_in_window_lms = n_in_window / n_total`.
- `frac_in_window_nolms` counts `E_c / f_eff` against the **same** window.
  - This comparison is biased in favour of "with LMS", because the window was derived from the LMS-corrected data.
  - It is a sanity indicator only. The unbiased measure is the pre σ_E/E with `--lms-corr` on against off (§9 V6).
- `flag_mu_off = |μ_M − mu_global| > 2·σ_M` (user decision, 2026-10-21; originally compared with E'_M).
  - It marks modules whose peak is far from the path's global peak, or fits locked onto the wrong peak. The global scale offset itself is reported separately as `mu_global/Ē' − 1`.

## 6. Interfaces and data changes

- New functions, all internal to the new script:
  - `energy_loss`, `expected_energy_ep`, `expected_energy_ee`, `module_theta_deg`
  - `resolve_lms_path`, `load_lms_gain_corr`, `LmsTable.batch_of`
  - `cluster_lms_feff`
  - `fit_module_cb`, `_select_pool`, `fit_elastic_window_gains`
  - `write_lms_tsv`, `plot_lms_factors`, `write_window_tsv`
- Modified in the copy:
  - `fit_per_module_gains`: tuple unpacking only, and it now uses `_select_pool`.
  - `main`: the CLI, the loop, the calibrate dispatch, and the output additions.
- Reused unchanged: `gauss`, `crystal_ball`, `fit_peak_cb`, `per_module_peak`, `fit_peak`, all existing plot and TSV code.
- No database, C++ or binding changes.
- New runtime dependency: `uproot` (only when `--lms-corr` is not `off`).
- **Interface the later scan-script step will need** (not done now):
  - The scan calls the script as a subprocess with `<evio> <out_stem>` plus flags. It would change `_BENCHMARK_SCRIPT` and add `--hc-calib-file` per run.
  - It **must** pass `--lms-corr off` (or an explicit path) for every run without a gain_corr file. Under the default `auto`, such runs exit with code 2, and only 21 runs have files today.
  - It would replace `--cal-method cb-peak` with `elastic-window` or keep it.
  - `parse_tsv` needs no change.
  - Per-split scans see one batch per split (§3.1).

## 7. Risks and compatibility

| Risk | Mitigation |
|---|---|
| A single split sits inside one batch, so the LMS effect is absorbed by the refinement and looks like a null result. | Primary LMS validation on 8 splits spanning the run (§9 V6). The expectation is stated in §3.1. |
| Failed-fit sentinel 1.0 gives spurious relative jumps of about 2%. | The validity rules in §5.2, and an uncorrected count reported in the outputs. |
| Moller electrons inside the plain pre-fit range (r ≲ 56 mm) or the window. | The Moller floor in §5.4 (binds for r ≲ 88 mm), plus `flag_moller_overlap`. The 15–20 exclusion removes r ≤ 103 mm only partially, because of the asymmetric block (§3). |
| Low per-module statistics: on one split most modules have fewer than 100 events; on eight splits most have 100 or more. | The global-σ fallback. One-split runs will mostly use it, so per-module σ is meaningful only on multi-split samples. |
| The gain_corr file must cover the whole run from split 0. The first batch is the run's first 1600 LMS events. | Always run `replay_gainCorr` on **all** splits of the run (§8). A partial-run file would move t0. |
| Python and C++ integral definitions may differ (unverified). | Only ratios to batch 0 are used, so a constant definition offset cancels. |
| uproot is not installed. | Prerequisite step, plus a hard error with an install hint. `--lms-corr off` still works. |
| Overwriting the git-tracked `database/.../prad_025320_gain_corr.root` while regenerating. | Use the `PRAD2_DATABASE_DIR` copy in §8. Never run in place. |
| Calibration: runinfo picks 3.5 GeV for 25320. | Every validation command passes `--hc-calib-file …2p2_fadc.json`. |

## 8. Producing `gain_corr.root` at 3-min batches

- **Required?** Yes, to meet the user's 3-minute spec. The existing `database/gain_factor/gain_correction/prad_025320_gain_corr.root` uses `-b 4000` (about 7.5 min, 9 batches).
- It is **not** required for developing or validating the code. The -b 4000 file exercises the identical path, and §9 V2–V6 run on it first.
- The regeneration needs **all 400+ splits** of 025320, because LMS events are spread over the run (about 70 per split). Only 8 splits are local, which would give about 550 LMS events and zero full 1600-event batches.
- So it must run on a host that has the full run (see Open Question 1).

```bash
# On a host with the full run.  REPO = prad2evviewer checkout with build/ built; RAW = dir holding prad_025320.evio.*
DB=/scratch/$USER/db_b1600
cp -R "$REPO/database" "$DB"     # 102 MB; the tool writes $DB/gain_factor/gain_correction/ (hardcoded, RECREATE)
mkdir -p /scratch/$USER/lms_025320  # the tool does not create -o; Phase 1 fails without it
PRAD2_DATABASE_DIR="$DB" "$REPO/build/bin/prad2ana_replay_gainCorr" \
    -o /scratch/$USER/lms_025320 -b 1600 -r 24827 -j 4 \
    "$RAW"/prad_025320.evio.*       # options first: BSD getopt stops at the first positional
# Expect about 22 full batches (about 36k LMS events / 1600; the trailing partial is dropped).
cp "$DB/gain_factor/gain_correction/prad_025320_gain_corr.root" \
   "$REPO/analysis/pyscripts/output/lms/prad_025320_gain_corr_b1600.root"
```

- Pass the file glob, not the directory, because directory mode picks up every `.evio` file in it.
- `-r 24827` is the runinfo default, written out for reproducibility. In first-batch mode the reference nearly cancels.

## 9. Validation (run 025320; all from `analysis/pyscripts`)

Common setup:

```bash
cd /Users/L00338853/work/PRad/prad2evviewer/analysis/pyscripts && source setup_env.sh
EVIO=/Users/L00338853/work/PRad/data/evio
CAL=$PRAD2_DATABASE_DIR/calibration/calibration_factor_2p2_fadc.json
GC4000=$PRAD2_DATABASE_DIR/gain_factor/gain_correction/prad_025320_gain_corr.root
OUT=output/lmstest; mkdir -p $OUT
COMMON="--calibrate --target-peak 2239.5 --exclude-rowcol 15 20 15 20 --signal-min 500 \
  --pre-window 100 200 --e-range 500 3500 --bin-width 20 --inner-radius 200 \
  --min-seed-energy 200 --min-per-module 30 --hc-calib-file $CAL \
  --daq-config $DAQ_CONFIG --hc-map-file $HC_MAP"      # = scan-script settings for 25320 + 2p2 calib
```

**V1. E' port matches C++.** The reference values below come from the real `PhysicsTools.cpp`, compiled and run on 2026-10-21. The Python port in double precision differs by at most 2.5e-4 MeV for "ep".

```bash
python3 - <<'EOF'
import benchmark_hycal_timing_test_20261021 as B
ep = {0:2237.454102, 0.3:2237.380859, 0.5:2237.250244, 1.0:2236.639893,
      1.5:2235.623291, 1.83:2234.729980, 2.5:2232.376221, 3.5:2227.524170}
ee = {0.3:2110.535400, 0.5:1917.009644, 1.0:1340.913086, 1.83:690.024841}
dep = max(abs(B.expected_energy_ep(t, 2239.51) - v) for t, v in ep.items())
dee = max(abs(B.expected_energy_ee(t, 2239.51) - v) for t, v in ee.items())
print(dep, dee); assert dep < 1e-3 and dee < 0.5    # ee: C++ float32 cancellation near cos^2 ~ 1
EOF
```

The reference values were produced by a harness that links the real `analysis/src/PhysicsTools.cpp` against the `build/` static libraries with `root-config` (root 6.40.02).

**V2. LMS reader on the existing file.** Call `load_lms_gain_corr($GC4000, p.hycal, "first-batch")` after `C.setup_pipeline(evio_path=$EVIO/prad_025320.evio.00400, hc_calib_file=$CAL, …)`, then assert:

- 9 batches, and every factor in batch 0 is 1.0.
- For modules corrected in batches 1–8, all factors are in [0.985, 1.015]. The observed range is 0.988–1.010, which satisfies "within 1.2% in-run".
- W1141 (`module_by_id(2141).index`) is 1.0 in every batch, and every `is_glass()` module is 1.0.
- The uncorrected count is exactly 2 real modules in every batch: W1124 (sentinel) and W1141 (zero). The beam-hole rows have no module and are not counted.
- `batch_of(1221) == 0`, `batch_of(9387322) == 1`, and `batch_of(10**8) == 8`.

Also run the same load in `ref-run` mode and check that factor 0 for a sample module equals `(c[1]+c[2])/2` from uproot.

**V3. Regression against the old script.** With LMS off and the old method, the output must be bit-identical:

```bash
python3 benchmark_hycal_timing.py              $EVIO/prad_025320.evio.00400 $OUT/old_00400_50k     $COMMON --cal-method cb-peak --max-events 50000
python3 benchmark_hycal_timing_test_20261021.py $EVIO/prad_025320.evio.00400 $OUT/newoff_00400_50k $COMMON --cal-method cb-peak --lms-corr off --max-events 50000
diff <(grep -v '^#' $OUT/old_00400_50k.tsv) <(grep -v '^#' $OUT/newoff_00400_50k.tsv) && echo IDENTICAL
```

**V4. Identity in batch 0.**

- Split 00050 (event number about 8M) lies in batch 0 of the -b 4000 file. With `--lms-corr $GC4000`, every factor is 1.0.
- So the non-`#` lines of `<out>.tsv` must equal the `--lms-corr off` run, both using `--cal-method elastic-window` and `--max-events 50000`.
- `<out>_lms.tsv` must show only batch 0 touched.
- If 00050 straddles into batch 1, take the first 50k events, which are well inside batch 0. Check `batches_touched` in the header.

**V5. Single split with non-trivial factors** (00400, batch 8, full split, elastic-window, LMS on and off). Check:

- Exit code 0. With LMS on, all 6 files are written: `.tsv`, `.png`, `_calibrated.png`, `_lms.tsv`, `_lms.png`, `_window.tsv`. With LMS off, 5 are written (no `_lms.png`).
- `_lms.tsv` shows batch 8 touched.
- `_window.tsv` has rows for both paths. The `used` modules are ≥ 60% of pool modules with ≥ 30 events. `n_in_window ≥ 30` for used modules. There are no NaN gains.
- The post-refinement in-window median equals `T_M` to 1e-6 relative. This is the single-pass assertion; recompute it from `_clusters.tsv` or in-script.
- The CB post μ (legacy and gated) is within 10 MeV of Ē' (about 2236 MeV). The median estimator leaves a few-MeV bias from the CB tail inside the window; the technical note reports 6 MeV for the all-event median.
- `|σ_E/E post (LMS on) − σ_E/E post (LMS off)| ≤ 0.05 pp`. This is the absorption expectation of §3.1.
- Report the pre σ_E/E on and off, and the median `frac_in_window_lms` and `frac_in_window_nolms`.
- The sanity range for the median `frac_in_window_lms` is [0.4, 0.95]. Modules at r ≈ 53–75 mm (the bimodal ones) are expected lower.

**V6. Multi-split: the real LMS test.** Run on all 8 local splits, which touch batches {0,1,2,4,5,6,7,8} (about 45 min of run time; about 12 min of CPU each):

```bash
G="$EVIO/prad_025320.evio.*"
python3 benchmark_hycal_timing_test_20261021.py "$G" $OUT/new_on_8s   $COMMON --lms-corr $GC4000 &
python3 benchmark_hycal_timing_test_20261021.py "$G" $OUT/new_off_8s  $COMMON --lms-corr off &
python3 benchmark_hycal_timing.py              "$G" $OUT/old_cbpk_8s $COMMON --cal-method cb-peak &
wait
python3 -c "import benchmark_hycal_timing_scan as S; from pathlib import Path; print(S.parse_tsv(Path('$OUT/new_on_8s.tsv')))"
```

Report a table for legacy and gated. Its rows are `new_on_8s`, `new_off_8s` and `old_cbpk_8s`. Its columns are:

- n modules, n inner events
- σ_E/E: Gaussian pre, CB pre, Gaussian post, CB post
- post CB μ against Ē'
- the number of `sigma_source=module` and `flag_mu_off` modules

Expectations (report these; do not hard-fail on them):

- The LMS-on pre σ_E/E is less than or equal to LMS off. The gain should be at most about 0.1 pp. A larger gain, or a degradation, is a red flag to investigate.
- `parse_tsv` returns without error and has finite σ fields.
- The elastic-window and cb-peak post σ_E/E are reported side by side, with no pass threshold. That comparison is the user's physics question.

**V7. After the -b 1600 file exists.** Repeat V2 (with batch count about 22 and the factor range within 1.2%), V5 and V6 with `--lms-corr output/lms/prad_025320_gain_corr_b1600.root`. The 8 splits then touch about 8 distinct batches.

## 10. Acceptance criteria

1. The new script exists. `benchmark_hycal_timing.py`, the scan script, C++ code, the technical note and `database/` are unmodified. The Todo P0 snapshot of their `git diff --stat` is unchanged at the end. (The tree already has unrelated modifications, so a bare `git status` check does not work.)
2. V1–V4 pass as written, including the V3 identical diff.
3. V5 and V6 complete. All outputs are written, the V5 checks hold, and the V6 table is reported.
4. With `--lms-corr off` the script runs without uproot installed.
5. A missing gain_corr file under `auto` exits with code 2 and prints the §8 command.

## 11. Addendum (2026-10-21): all beam energies, target kinematics, Moller exclusion

User request: support e-p runs at every beam energy, and on the X17 runs do gain refinement only. User decisions: (a) the X17 runs keep the LMS correction and anchor the refinement to the e-Ta elastic E'; (b) a module is excluded from refinement when its elastic-Moller gap is below 3σ_M, which replaces `flag_moller_overlap`, while the Moller-aware pre-fit floor stays.

- **Energy scaling.** `E_REF = 2239.51` MeV, `kE = round(E_beam/E_REF, 6)`, `kS = round(sqrt(E_beam/E_REF), 6)`. Both are exactly 1 at 2.24 GeV, so outputs there do not change.
  - Defaults when the flag is not given: `--signal-min` 500·kE, `--e-range` (500, 3500)·kE, `--min-seed-energy` 200·kE, `--bin-width` 20·kS (σ_exp = 3.3%/√E), `--target-peak` E_beam.
  - The peak-shape fit windows scale with kS: 150, 400, ±600, the p0 widths and the σ bound. The broad ±800 MeV ranges scale with kE.
  - `PREFIT_RANGE` is multiplied by `1/kS`, so it is a fixed multiple of σ_exp.
  - Explicit CLI values override the defaults. The resolved values are printed on stdout and in the `.tsv` `#` header.
- **Calibration check.** The effective file is `--hc-calib-file` if given, else the runinfo file. Its energy is inferred from the file name (`0p7`, `2p2`, `2p1`, `3p5`, `x17`) and matches if it is within 3% of E_beam. That separates `2p1` (2108.8) from `2p2` (2239.5).
  - A runinfo default that mismatches, or has no recognisable tag, exits 2 and lists the matching files.
  - An explicit file that mismatches prints a warning banner and the run continues, because an explicit CLI value wins.
  - Rationale: runinfo is known to give the 3.5 GeV file for 25320 and for the 0.73 GeV runs.
- **Kinematics.** `--target {H2,C,Ta}`, defaulting to a per-run map:
  - H2: 24386, 24503, 24512, 24578, 24590, 24603, 25228, 25257, 25260, 25299, 25320. The 24603, 25228 and 25257 assignments come from `data/evio/README.md`.
  - C: 25308.
  - Ta: 25800, 26138, 26230.
  - An unmapped run with `--cal-method elastic-window` exits 2.
  - H2 uses the `"ep"` port. C (12 u) and Ta (180.948 u, u = 931.494 MeV) use plain elastic recoil `E·M/(M+E(1−cosθ))` with the same `EnergyLoss`, because the C++ `ExpectedEnergy` has no nucleus mode.
  - T_M, Ē', the Moller floor and the gap all use this E'_elastic.
- **Moller exclusion.** `gap = (mu_global/Ē')·(E'_el − E_ee)`, which puts the gap in the same measured-energy units as σ_M. The module is not used if `gap < 3·σ_M`.
  - `_window.tsv` gains `moller_gap_nsigma` and `exclude_reason` (`moller_gap`, `few_in_window`, `no_prefit_events`, `bad_median`), and loses `flag_moller_overlap`.
- **Validation.**
  - V1, V3 and V4 must still hold.
  - The 025320 8-split result must reproduce the previous V6 to within 0.01 pp, both with explicit flags and with auto defaults.
  - Smoke runs on the local 0.73 GeV, 3.49 GeV, carbon and Ta splits.
  - The runs other than the 25320 family and 25308 have no `gain_corr` file, so those smoke runs use `--lms-corr off`.

## Open Questions

1. **Where do we regenerate at -b 1600?** It needs all 400+ splits of 025320, and only 8 are local. Options are the JLab farm or wherever the full run lives, or copying the run (about 900 GB). *Default:* develop and validate on the existing -b 4000 file, then run V7 once you produce the -b 1600 file.
2. **May we install `uproot` into the miniforge env** (`python3 -m pip install uproot`)? It is missing today and needed to read the ROOT file. *Default:* yes. Without it only `--lms-corr off` works.
3. **Moller-aware pre-fit floor.** `lo = max(E'(1−15%), (E' + E_Moller)/2)` is an addition to the agreed −15%/+10% range.
   - It binds for r ≲ 88 mm.
   - It matters for r ≲ 56 mm, where the Moller peak lies 6–14% below E' and would otherwise sit inside the fit range.
   - *Default:* on. Say if you want the plain −15% edge.
4. **A window that reaches below the Moller floor** (`w_lo < lo`, possible at r ≈ 33 mm with the global-σ fallback). *Default:* keep the window literal and flag the module (`flag_moller_overlap`). The alternative is to clip `w_lo` to `lo`.
5. **`--min-per-module` counts in-window events** (the refinement sample), not all events of the module. *Default:* in-window. This uses fewer modules than the old method on the same data.

## Todo

- [x] P0. Snapshot the protected paths, run from the repo root: `git diff --stat -- analysis/pyscripts/benchmark_hycal_timing.py analysis/pyscripts/benchmark_hycal_timing_scan.py database analysis/src analysis/tools analysis/include docs/technical_notes > <scratch>/protected_before.txt` (any location outside the repo).
- [x] P1. After Open Question 2: `python3 -m pip install uproot` in `/Users/L00338853/Software/opt/miniforge3`. Verify with `python3 -c "import uproot, prad2py"` after `source setup_env.sh`.
- [x] 1. Copy the working-tree `benchmark_hycal_timing.py` to `benchmark_hycal_timing_test_20261021.py`. Update the module docstring (methods, new flags, new outputs, the single-pass note).
- [x] 2. Add `energy_loss`, `expected_energy_ep`, `expected_energy_ee` and `module_theta_deg` (§5.4). Run V1.
- [x] 3. Add the CLI flags from §5.1, and the missing-file and missing-uproot exits.
- [x] 4. Add `resolve_lms_path` and `load_lms_gain_corr` with the §5.2 validity rules, `batch_of` and `factors_list`. The uproot import is lazy. Run V2.
- [x] 5. Event loop: `ev_num`, the batch lookup, `n_ev_per_batch`, and the factor applied in both `add_hit` paths (§5.3 items 1–2).
- [x] 6. Add `cluster_lms_feff`. Extend the event tuples to 9 fields and fix every unpacking (`fit_per_module_gains`, the cluster-TSV writer) (§5.3 items 3–4).
- [x] 7. Factor out `_select_pool()` and use it in `fit_per_module_gains`. Run V3; it must be identical.
- [x] 8. Add `fit_module_cb` and `fit_elastic_window_gains`: global σ, the Moller floor, the fallback, the fixed window, the in-window median, single pass (§5.4).
- [x] 9. Wire the `--cal-method` dispatch, with the Ē'-centred cal range and plot line for elastic-window.
- [x] 10. Outputs: the TSV header comment lines, `_lms.tsv`, `_lms.png`, `_window.tsv`, the stdout LMS and window blocks, and the title tags (§5.5).
- [x] 11. Run V4 and V5. Record the numbers.
- [x] 12. Run V6. Record the comparison table in the hand-off report.
- [ ] 13. (Host with the full run) Produce `prad_025320_gain_corr_b1600.root` per §8, then run V7.
- [x] 14. Addendum §11: energy scaling, calibration check, `--target` kinematics, Moller exclusion.
- [x] 15. Addendum validation: V1/V3/V4, 025320 8-split explicit and auto, smoke runs at the other energies and targets.

## Fit-model note (2026-10-21, user observation from the full scan `output/scan_20261021_full/`)
- **The Crystal Ball (CB) fit describes the elastic peak well** for 3.49 GeV e-p (H2: 24386, 24503, 24512), 2.24 GeV carbon (25308) and the 2.24 GeV X17 Ta runs (25800, 26138, 26230).
- **The Gaussian fit works better** for the other runs: 2.24 GeV e-p H2 (25228, 25257, 25260, 25299, 25320) and 0.73 GeV e-p H2 (24590, 24603).
- Quote σ_E/E from the preferred fit model for each run type.
- **Run 24578 is excluded** from the scan (2026-10-21): its trigger was not the standard one, leaving only ~2000 physics-trigger events per split and no per-module refinement. It is excluded via `_EXCLUDED_RUNS` in `benchmark_hycal_timing_scan.py`.

---

## Energy-resolution summary: full scan `output/scan_20261021_full/` (2026-10-21)

Post-refinement σ_E/E (%), elastic-window method, inner ring r < 200 mm. Values are means over splits; ± is the split-to-split standard deviation (legacy only). **Legacy** is the largest in-window peak per module (130–200 ns). **New** is the gated clustering path. Δ = New − Legacy, in percentage points (pp), using the preferred fit model (marked ★; see the Fit-model note). Run 24578 is excluded.

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

**Observations**
- **e-p and carbon runs:** New is 0.02–0.04 pp worse than Legacy on every run, consistently across splits.
- **X17 Ta runs:** New is better (−0.02 to −0.14 pp). These runs have very few modules (9–34) and large split-to-split spreads, so the gain is suggestive, not established.
- **0.73 GeV:** the refinement barely changes σ_E/E. The Gaussian value (4.6–4.8 %) is driven by the radiative tail.
- 24603, 25228 and 25257 were re-run in full (24 splits, 2026-10-21); 25228 and 25257 are unchanged, and every other run is byte-identical.
- **24386:** split 00060 timed out and is excluded.
