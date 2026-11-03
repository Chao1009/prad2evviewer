# e-p energy resolution, step by step: from EVIO to σ_E/E and a/√E ⊕ b (R5nr, 2026-11-02)

## Scope

This note traces how the current working-reference numbers for the e-p (H2 target) points at 0.73, 2.24 and 3.49 GeV are produced. The reference is **R5nr**: R5 with the online gain refinement switched off, commit 71e049f. The executed benchmark is `benchmark_hycal_timing_test_20261021.py` with sha1 08b19869, which is enforced by `output/r5_20261028/tools/run1_sha.sh:5`. HEAD's working copy has the same sha.

- The headline is the **legacy path, Crystal-Ball (CB) core σ/μ**. Where the gated path differs, the step says so.
- 2.24 GeV C and Ta are out of scope. They are mentioned only where they share a step.
- Every flag value comes from the R5nr job lists, `output/r5_20261028/jobs/abl_norefine.txt` (41 reference splits) and `final_nr.txt` (all 118 splits), not from argparse defaults.

**Path and label conventions**

- Paths are relative to `analysis/pyscripts/` unless they start with `prad2det/`, `prad2dec/`, `python/` or `database/`; those are relative to the repository root.
- `B` = `benchmark_hycal_timing_test_20261021.py`.
- `HC.cpp` = `prad2det/src/HyCalCluster.cpp`; `HC.h` = `prad2det/include/HyCalCluster.h`.
- `WA.cpp` = `prad2dec/src/WaveAnalyzer.cpp`.
- `$O` = `output/r5_20261028`; `$NR` = `output/resolution_fit_20261029_r5`.
- Line numbers are at HEAD 71e049f.
- "Verified" means read in the code or an output file. "Inference" marks anything not checked directly.

Example split used throughout: run 24503, split 00020 (3.49 GeV H2), with outputs in `$O/scan/r5nr_splits/scan_024503/split_00020/`.

## Flow diagram

```
 PREPARATION (offline, once per run type / run)                         PER-SPLIT JOB (one EVIO split)
 ------------------------------------------------                      ------------------------------
 P1 spline template bank best-50 + ampcorr f(A), I_int                 S0  r5_jobs.py -> chain_r5nr.sh -> run1_sha.sh -> run1.sh
    spline_template_extract.py select/build/ampcorr ---------------+       -> merge_cfg_wrap.py (cluster cfg r5) -> B main
 P4 fixed inputs: hc-calib json, runinfo general.json,             |   S1  EVIO decode, Physics records, trigger == 0x100
    daq_config analyzer, Geant4 profiles, reconstruction_config_r5 |   S2  pedestal / WaveAnalyzer (glass, fallback; defines I_int)
                                                                   +-> S3  matched-filter fit -> E_u = A * I_int * f(A) * kappa(=1)
 P2 per-crystal gains (S2 of plan_r5)                                  S4  time: pre-window (100,200) ns raw; no offsets; LMS off*
    dump pass (B --cal-method percrystal --cluster-dump auto)          S5  module MeV = cal_factor * E_u ; legacy (1 hit/module) vs gated (all hits)
    -> percrystal_calib.py fit pass-0 -> pass-1 Moller-rejected        S6  island clustering (C++ HyCalCluster, Geant4 profile, seed 30 MeV)
    -> $O/gains/<run>/gains_<run>_all.tsv ----------------------+      S7  cluster-E corrections: non-linear (inert), leakage OFF, bias OFF
                                                                |      S8  log-weighted centroid (W0 = 5.2) -> x, y ; theta = atan(r/hycal_z)
 P3 position maps (S3 of plan_r5)                               |      S9  exactly one cluster > 500 kE ; Ebeam, hycal_z from runinfo
    map scan (B with P2 gains, no pos-corr)                     +----> S10 per-crystal gains: E <- sum_i g_i E_i
    -> norefine_scan.sh view -> position_corr.py map/fit  -----------> S11 position correction: E <- E / R(u, v)
    -> $O/poscorr/maps_norefine/<slug>/poscorr_r5nr_*.json             S12 Moller kinematic reject (n=3, m=2) using E'(theta), E_ee(theta)
                                                                       S13 E'(theta) e-p elastic (used in S12, S14; not a rescale)
                                                                       S14 inner pool + elastic-window pass (--calibrate); refine OFF (gain := 1)
                                                                       S15 histogram of absolute E [MeV], target_ref +- 800 kE, 20 kS bins
                                                                       S16 fits: fixed Gauss -> CB core [mu0-w, mu0+w/2]; +-2 sigma Gauss (secondary)
                                                                       S17 result.tsv: sigma_cb/mu_cb (legacy row = headline)
 AGGREGATION                                                           
 S18 compare_scan_dirs.py (pair splits) -> nr_table.py mean +- SE over splits -> tables/r5nr_table.md
 S19 phase0.py / nr_table.py fit2: sigma^2 = a^2/E + b|b| on the three H2 points (1 dof)
 * LMS on for runs 25308 (C) and 25320 (2.24 H2) only.
```

Order inside one job (verified, B:3225-3400): S10 gain file → S11 pos-corr → S12 Møller tag → S14 window pass → S15/S16 spectrum and fits.

---

## Part A. Preparation (offline ingredients)

### P1. Spline template bank (best-50) and I_int / f(A)

**What.** For each module and run type there is one sampled unit-peak pulse template. It is a penalised cubic B-spline fitted to the 50 highest-energy single-pulse candidates per module, pooled over the run type's splits and iteratively aligned.

- I_int is the phase-averaged C++ walk-out integral (see S2) of the noiseless unit template. There is no stored κ.
- `ampcorr` adds a per-run-type amplitude correction f(A). It makes the MF energy agree with the C++ integral across amplitude.

**Where.** `spline_template_extract.py`:

| Function | Line |
|---|---|
| `cmd_select` | 1113 |
| `Selector` | 501 |
| `cmd_build` | 2826 |
| `build_shape` | 1964 |
| `build_lambda_select` | 2147 |
| `cmd_ampcorr` | 4051 (`ampcorr_rule` 3503) |
| `cmd_check` | 5556 |
| CLI | 5705 |

The I_int definition is `mf_pulse.build_iint` (`mf_pulse.py:600-634`; `iint` :636-640).

**Provenance** (from the npz `meta.argv` and the notes; git 3566071, 2026-10-26):
```
spline_template_extract.py select <evio split> --trigger physics \
  --mf-bank output/templates_spline/3p49GeV_H2/templates_3p49GeV_H2.npz \
  --per-module-best 50 --emax-frac 1.2 --mf-strict --max-modules-per-event 20 \
  --exclude-trigger-bits 0x1000000 --out .../3p49GeV_H2/pulses/<run>_<split>.npz
spline_template_extract.py build --pulses ".../3p49GeV_H2/pulses/*.npz" --run-type 3p49GeV_H2 \
  --pool-order energy --n-max 50 --fit-points direct --shape-sigma 0.012 \
  --lambda-select heldout --lambda-chi2-sigma full --lambda-min 1e2 --lambda-fallback 1e3 \
  --converge-tol 0.01 --converge-df 5e-3 --max-iterations 10 --osc-escalate --out .../templates_3p49GeV_H2
spline_template_extract.py ampcorr --out $B/ampcorr.json --templates $B/{0p73GeV_H2,2p24GeV_H2,2p24GeV_C,2p24GeV_Ta_X17,3p49GeV_H2}/templates_*.npz
```
Recipe sources: `notes/PbWO4_energy_resolution_20261025.md:985-1002` and `notes/PbWO4_energy_resolution_20261026.md:110-160`. Logs: `output/templates_spline_best50_20261026/{ampcorr.log, 3p49GeV_H2/build.log, logs/select_*.log}`.

**Banks used in R5nr** (counted in `final_nr.txt`):

| Point | Bank | Splits | Template runs |
|---|---|---|---|
| 3.49 GeV H2 | `output/templates_spline_best50_20261026/3p49GeV_H2/templates_3p49GeV_H2.npz` | 27 | 24386, 24503, 24512 (28 selection splits) |
| 2.24 GeV H2 | `.../2p24GeV_H2/templates_2p24GeV_H2.npz` | 42 | 25228, 25257, 25260, 25299, 25320 |
| 0.73 GeV H2 | `.../0p73GeV_H2/templates_0p73GeV_H2.npz` | 16 | 24590, 24603 |
| 2.24 GeV C | `.../2p24GeV_C/...` | 8 | — |
| Ta | `output/templates_spline_ta_fix_20261028/n200/templates_2p24GeV_Ta_X17.npz` | 25 | — |

At 3.49 GeV, 1065 modules have their own shape and 87 fall back to the run-type average. ampcorr at 3.49 GeV is a 17-point f(A) over A = 218-1932 ADC, with f between 0.9940 and 0.9994.

**Caveats**

- The template pool contains the evaluated splits, for example 24503.00020. Treating this as in-sample is an inference; confirming it needs `holdout_keys`.
- The npz `argv` and the logs still name `esc_stage/`, but the files were moved one level up (`notes/..._20261026.md` §8.5).
- The select step was seeded with the earlier production bank `output/templates_spline/...`.

### P2. Per-crystal gains (plan_r5 S2)

**What.** One gain table per run, fitted offline from a "dump" of the calibration pool.

- `percrystal_calib.py` does **not** fit a peak per crystal. It solves one ridge-regularised joint least-squares problem, Σ_i g_i E_ki = E'_k.
  - Events k: pool clusters whose energy falls in the elastic window of their seed module.
  - Unknowns: g_i for every member module of those clusters.
  - Prior: a pull toward g0 = E'_mean/μ_global.
- Pass 0 has no rejection. Pass 1 repeats the fit after removing clusters that pass-0 gains tag as Møller (kinematic rule n=3, m=2; see S12).
- There is no further iteration beyond pass 0 → pass 1. Inside each fit, the `--fixed-windows` loop converges in 2 passes.

**Where**

| Item | Location |
|---|---|
| Dump jobs | `$O/tools/r5_jobs.py:124-127` (`variant == "dump"`), listed in `$O/jobs/dump.txt` (119 = 118 splits + 24386.00060) |
| Solver `solve()` | `percrystal_calib.py:114-279`; normal equations :90-111 |
| `cmd_fit` | :589-694 |
| `moller_kin_mask` | :560-574 |
| `write_gain_table` | :298-328 |
| `read_gain_table` | :358-379 |
| Two-pass driver | `$O/tools/fit1.sh:6-8` |
```
# r5_jobs.py:124-127  (dump jobs: no gains, no map, no Moller flags)
if variant == "dump":
    t = base(tok, run, tail(f"{O}/gains/dumps"))
    drop(t, "--gain-file"); drop(t, "--pos-corr"); drop(t, "--write-cluster-tsv", 0)
    return t + ["--cal-method", "percrystal", "--cluster-dump", "auto"] + PROF + gr
# percrystal_calib.py:219-224  (one joint ridge solve)
Nm, b = _normal_equations(ptr_s, inv_m, E_s, tgt_s, M)
Nm *= w
b *= w
D = np.diag(Nm).copy()
A = Nm + lam * np.eye(M)
gsol = np.linalg.solve(A, b + lam * g0)
```

**Settings.**
- Dump jobs: the R5nr flags of S0, without gain file, map and Møller flags, plus `--cal-method percrystal --cluster-dump auto --cluster-profile geant4`. The executed `job.cmd` files (e.g. `$O/gains/dumps/scan_024503/split_00020/job.cmd`, Oct 29 07:19) carry no `--gain-refine` flag, because benchmark d91af9b7 did not have it. `jobs/dump.txt` was regenerated afterwards (13:47) with `--gain-refine off` appended. The physics is unaffected, because refinement is skipped in percrystal mode.
- `fit1.sh`, pass 0: `percrystal_calib.py fit --dumps "$O/gains/dumps/scan_$R/split_*/result_dump.npz" --splits all --fixed-windows --out $O/gains_pass0/$R/`.
- `fit1.sh`, pass 1: the same arguments plus `--moller-reject kin --moller-kin-n 3 --moller-kin-m 2 --moller-gains $O/gains_pass0/$R/gains_${R}_all.tsv --out $O/gains/$R/`.

**Outputs.** `$O/gains/<run>/gains_<run>_all.{tsv,json}`.
- The header carries `path=legacy g0=…` and the calibration-file sha1.
- Columns: `path module index row col g sigma_g n_seed n_member D_ii prior_frac g_seedindexed used`.
- Run 24503: `g0 = 0.929768`, 568 modules listed. Unlisted modules get g0 (B:3228).

**Caveats**
- **In-sample.** `--splits all` pools every split of the run, and the same splits are evaluated later. The only exception is 24386.00060, which is never evaluated. This is stated in plan_r5:244.
- The JSON field `n_rejected` (for example 36679 for 24386) counts all tagged dump clusters. The pool itself shrinks by only 1.1 %.
- The dumps were produced with benchmark d91af9b7, while evaluation uses 08b19869. plan_r5:52 says S2 is unaffected; this was not re-verified here.

### P3. Position (bowl) correction maps (plan_r5 S3)

**What.** The maps give R(u,v) = E/E'(seed) as a function of the cluster centroid relative to the seed-crystal centre, (u,v) = (x_cl − x_seed, y_cl − y_seed)/module size. They are built from a map scan of B run with the P2 gains and no `--pos-corr`.

- `position_corr.py map` bins (u,v) on 12×12 bins over [−0.6, 0.6]².
- `fit --form bilinear` produces a `poscorr_v1` grid.
- `fit --per-module` produces a `poscorr_hybrid_v1` map. A seed module with at least 3000 clusters gets its own poly8; otherwise the embedded bilinear fallback is used. Clusters with `moller_rej == 1` are dropped.

**Where**

| Item | Location |
|---|---|
| Map jobs | `$O/tools/r5_jobs.py:131-135` → `$O/jobs/map.txt` (118 jobs). The executed `job.cmd` (e.g. `$O/poscorr/scan/scan_024503/split_00020/job.cmd`) has no `--gain-refine` flag (d91af9b7 always refines). `map.txt` was regenerated later with `--gain-refine on`. The refine-off view (next row) removes the refinement. |
| Refine-off view | `$O/tools/norefine_scan.sh:12-13`: sets `cluster_E_cal := cluster_E_raw` where not NaN → `$O/poscorr/scan_norefine/` |
| Fitter | `position_corr.py`: docstring 1-60, `load_pool` 114-157, `cmd_fit_hybrid` 930-998 |
| Logs | `$O/poscorr/logs_norefine/` |
```
# norefine_scan.sh:12-13
awk 'BEGIN{FS=OFS="\t"} NR==1{if($9!="cluster_E_raw_MeV"||$10!="cluster_E_cal_MeV"){print "bad header" > "/dev/stderr"; exit 1}; print; next}
     {if($10!="nan") $10=$9; print}' "$d/result_clusters.tsv" > "$DST/$rel/result_clusters.tsv"
```

**Commands.** These are not logged verbatim. They are reconstructed from plan_r5:140-152, the JSON `source` field and the logs:
```
position_corr.py map --scan-dir $O/poscorr/scan_norefine --out $O/poscorr/maps_norefine --splits all
position_corr.py fit --map map_<slug>_all.tsv --form bilinear --out .../poscorr_r5nr_bilinear_<slug>_all.json
position_corr.py fit --per-module --scan-dir ... --splits all --slug <slug> --fallback <bilinear json> --out .../poscorr_r5nr_<slug>_all.json
```

**Outputs.** `$O/poscorr/maps_norefine/<slug>/poscorr_r5nr_<slug>_all.json`.

| Point | Map form | Built from |
|---|---|---|
| 3.49 H2 | hybrid: 52 own modules, bilinear fallback normalised to 1.000794 | 27 splits |
| 2.24 H2 | hybrid: 76 own modules | 42 splits |
| 0.73 H2 | bilinear only (symlink) | 16 splits |
| 2.24 C | hybrid: 28 own modules | — |
| Ta | symlink to the 2.24 H2 bilinear map | — |

**Caveats**
- The maps are in-sample: each is built on exactly the splits later evaluated, and with the in-sample P2 gains.
- The `r5_jobs.py:132` comment says refine-off maps never use the synthesized view, but the JSON `source.scan_dir` is `scan_norefine`. The comment is stale.
- The map estimator at 0.73 GeV is Gaussian, while the headline is CB.

### P4. Fixed inputs (not re-derived for R5nr)

**Per-point fixed inputs**

| Ingredient | Path | Used in |
|---|---|---|
| ADC → MeV factor per module (`--hc-calib-file`) | `database/calibration/old/calibration_factor_fadc_3p5.json` (3.49, 27 jobs), `calibration_factor_2p2_fadc.json` (2.24 H2 and C, 50 jobs), `calibration_factor_fadc_0p7.json` (0.73, 16 jobs); Ta `database/calibration/x17_*_fadc.json` | S5 |
| Ebeam, hycal_z, target position | `database/runinfo/general.json`: the entry with the largest `from_run` ≤ run. Run 24503 → 3488.43 MeV (:50); hycal_z 6265.887 mm (:27); target at (0,0,0) | S9, S13 |
| WaveAnalyzer configuration | `database/daq_config.json` `fadc250_waveform.analyzer` (:100-113) | S2, P1 |
| Geant4 shower profile | `database/cluster_profiles/prof_pwo.dat`, `prof_lg.dat`, loaded at `prad2det/src/PipelineBuilder.cpp:274-284` | S6 |

**Cluster configuration.** `$O/configs/reconstruction_config_r5.json` differs from `database/reconstruction_config.json` only in `min_center_energy` (30 vs 10).

The C++ per-module time offsets `database/hycal_time_offsets/time_calib.json` are loaded into `Module::time_offset` because the configuration has `time_calib: true`. B never reads that field, and `--time-calib off` means raw times (S4).

---

## Part B. Per-split job (inside one benchmark run)

### S0. Job orchestration and the R5nr command

**What.**
1. `r5_jobs.py` rewrites the R4 `job.cmd` files. It changes the wrapper, the output stem, the calibration path, the LMS blob, `--gain-file` and `--pos-corr`, then appends the R5 flags.
2. `chain_r5nr.sh` runs the jobs with `xargs -P 4` through `run1_sha.sh`. That script refuses to start (rc=97) unless B's sha1 is 08b19869.
3. `run1.sh` writes `job.cmd`, runs the job with a 3600 s limit, and appends `rc=`.
4. `merge_cfg_wrap.py` patches `_common.setup_pipeline`. It copies every scalar `hycal` key of the r5 configuration onto `p.hc_clusterer`, then calls `B.main`.

**Where**

| File | Lines |
|---|---|
| `$O/tools/r5_jobs.py` | flags 38-40, 157; pos-corr choice 110-112, 155-156 |
| `$O/tools/chain_r5nr.sh` | 1-25 (xargs at :21) |
| `$O/tools/run1_sha.sh` | 5-11 |
| `$O/tools/run1.sh` | 1-10 |
| `output/merge_check_20261028/tools/merge_cfg_wrap.py` | 54-76 (`patched` to `return B.main`) |

**Exact R5nr command for one 3.49 GeV H2 split.** This is line 1 of `$O/jobs/abl_norefine.txt` (the same command with the output in `$O/scan/r5nr_splits/` is `final_nr.txt:12`), run from `analysis/pyscripts`:
```
bash output/r5_20261028/tools/run1_sha.sh output/resolution_fit_20261029_r5/r5_norefine_splits/scan_024503/split_00020 \
 /Users/L00338853/Software/opt/miniforge3/bin/python3 /Users/L00338853/work/PRad/prad2evviewer/analysis/pyscripts/output/merge_check_20261028/tools/merge_cfg_wrap.py \
 --hycal-cfg output/r5_20261028/configs/reconstruction_config_r5.json /Users/L00338853/work/PRad/data/evio/prad_024503.evio.00020 \
 output/resolution_fit_20261029_r5/r5_norefine_splits/scan_024503/split_00020/result --calibrate --exclude-rowcol 15 20 15 20 \
 --pre-window 100.0 200.0 --inner-radius 200 --min-per-module 30 --lms-corr off \
 --hc-calib-file /Users/L00338853/work/PRad/prad2evviewer/database/calibration/old/calibration_factor_fadc_3p5.json \
 --window 16.0 --time-calib off --pulse-method mf --mf-scale 1.01 --mf-r-model 0.3 \
 --mf-templates output/templates_spline_best50_20261026/3p49GeV_H2/templates_3p49GeV_H2.npz --mf-scale 1.0 \
 --gain-file output/r5_20261028/gains/024503/gains_024503_all.tsv \
 --pos-corr output/r5_20261028/poscorr/maps_norefine/3p49GeV_H2/poscorr_r5nr_3p49GeV_H2_all.json \
 --write-cluster-tsv --moller-reject kin --moller-kin-n 3 --moller-kin-m 2 --moller-reject-gains on \
 --fit-select chi2 --cb-two-peak --cluster-profile geant4 --gain-refine off
```

**Where the outputs land**
- The job writes to the directory given as its first argument: `$NR/r5_norefine_splits/scan_024503/split_00020/` (`result.tsv`, `result_window.tsv`, `result_clusters.tsv`, `result_lms.tsv`, `result.png`, `result_calibrated.png`, `benchmark.log`, `job.cmd`).
- `chain_r5nr.sh:9-13` then copies the 41 reference dirs to `$O/scan/r5nr_splits/scan_024503/split_00020/`, which is the tree the tables read.
- The other 77 splits are written directly under `$O/scan/r5nr_splits/` (`final_nr.txt`). For a split that is not one of the 41, use its line from `final_nr.txt`.
- The copied `job.cmd` therefore still shows the `$NR/...` output path.

**Caveat.** `--mf-scale` appears twice. argparse keeps the last value, so κ = 1.0; the log shows `kappa 1`. See S3.

### S1. Event decode and trigger selection

**What.** B opens the EVIO split with the pybind `dec.EvChannel` and keeps only Physics records. It decodes FADC data without SSP and keeps events whose trigger word is exactly 0x100 (SSP0). Each kept event has 100 samples at 4 ns per channel. Sync, EPICS and control records drop out at the Physics-type test.

**Where.** B:2942-2975 is B's own loop. It does not use `_common.iter_physics_events`, which uses the 0xf00 mask (`_common.py:498-548`). The gate constant is `_common.py:478`.
```
2942  if ch.open_auto(fpath) != dec.Status.success:
2954  if ch.get_event_type() != dec.EventType.Physics: continue
2958  decoded = ch.decode_event(i, with_ssp=False)
2963  if fadc_evt.info.trigger_bits != C.PHYSICS_TRIGGER_BITS:   # == 0x100
```

**R5nr settings.**
- There is no `--trigger` flag; the gate is hard-coded.
- `--max-events` defaults to 0 (all events).
- The pipeline is set up by `_common.setup_pipeline` (`_common.py:202-255`), which loads runinfo, the DAQ configuration and calibration. `--hc-calib-file` is passed to `PipelineBuilder.set_hycal_calib` (`_common.py:220`).
- Example split: 158820 physics events, of which 157234 are triggered.

**Legacy vs gated.** Identical.

### S2. Pedestal and C++ WaveAnalyzer

**What.** `WaveAnalyzer::Analyze` proceeds as follows:
1. Smooths the samples.
2. Estimates the pedestal from the first 30 samples, using a median/MAD seed followed by a sigma clip.
3. Finds local maxima above max(5·ped_rms, 10 ADC).
4. Integrates each peak outward until the signal stays below 10 % of the peak height for 2 samples.
5. Times the peak by CFD.

In R5nr **mf mode**, PbWO4 channels do not use this integral for energy. WaveAnalyzer is run only for lead-glass channels, odd-length waveforms and clipped-channel fallback (B:1508, 1551-1600). For PbWO4 it supplies the *definition* of I_int (P1) and the `ped_nsamples` configuration.

**Where**

| Function | Location |
|---|---|
| `smooth` | `WA.cpp:244-260` |
| `findPedestal` | :272-346 |
| `findPeaks` | :358-585 (integral walk :411-441, CFD :511-547) |
| `Analyze` | :587-665 (threshold :642) |
| Configuration struct | `prad2dec/include/WaveAnalyzer.h:53-69` |
| Binding | `python/bind_dec.cpp:59-66` |
```
// WA.cpp:429-436 (integral walk)
float integral = buf[peak_pos] - ped_mean;
const float tail_cut = ped_height * cfg.int_tail_ratio;
auto is_below = [&](float v) { return v < tail_cut || v < ped_rms || v * ped_height < 0; };
```

**Settings** (`database/daq_config.json:100-113`): `smooth_order 2`, `peak_nsigma 5.0`, `min_peak_height 10.0`, `int_tail_ratio 0.1`, `tail_break_n 2`, `ped_nsamples 30`, `overflow 4095`, `clk_mhz 250`. `nnls_deconv.enabled` is false.

**Caveat.** Fallback pulses take the C++ integral with no κ and no f(A).

### S3. Matched-filter pulse fit: E_u = A · I_int · f(A) · κ

**What.** For each PbWO4 channel the fit runs these steps:
1. Subtract a robust pedestal: clipped mean of samples 0-29, with σ = max(1.4826·MAD, 0.3).
2. Scan the module template over onset t0 on a 1 ns grid.
3. Accept a pulse if A ≥ max(5σ_A, 5 ADC).
4. Refine (A, t0) with 3 Gauss-Newton steps.
5. Subtract the pulse and repeat, up to 4 pulses. A secondary pulse under a parent's span needs A ≥ 0.3·A_parent (`--mf-r-model 0.3`).

The energy in ADC-integral units is E_u = A·I_int(m)·f(A)·κ. I_int is the continuous phase-mean of the walk-out integral of the unit template.

**Where**

| Item | Location |
|---|---|
| Bank and reconstructor setup | B:2624-2636; `mf_pulse.build_bank` :773-823; `SampledTemplateBank.from_npz` :708-770 |
| `MFReconstructor.run` | `mf_pulse.py:1165-1242` |
| `_pedestal` | :845-853 |
| `_group` | :1026-1163 |
| `_gn` | :903-934 |
| Energy (`pulse_energy_units`) | :1245-1258, called in `MFBlock.flush` (B:1846-1847) |
```
# mf_pulse.py:1245-1258
def pulse_energy_units(bank, midx_of_pulse, t0, amp, kappa=1.0, ped_rms=None):
    if ped_rms is None: I = bank.iint(midx_of_pulse, t0)
    else:               I = bank.iint_floor(midx_of_pulse, t0, amp, ped_rms)
    if getattr(bank, "ampcorr", None) is not None:
        I = I * (bank.f_amp(amp) if bank.ampcorr[1].ndim == 1 else bank.f_amp(amp, midx_of_pulse))
    return amp * I * kappa
```

**R5nr settings.**
- Pulse fit: `--pulse-method mf`, `--mf-templates` as in P1, `--mf-r-model 0.3`, `--mf-nsigma 5.0` (default), `--mf-min-amp 5.0`, `--mf-kmax 4`.
- `--mf-ped-floor` off, `--mf-shape fixed`, `--single-pulse` off.
- **κ.** `--mf-scale 1.01 … --mf-scale 1.0` resolves to κ = 1.0 because argparse keeps the last value. The 1.010 default (B:2301) is the 1/V6 rescale that belonged to the older two-τ bank (`notes/PbWO4_energy_resolution_20261022.md:58,64`). The spline banks carry their own I_int and f(A).
- **Fallbacks.** No accepted pulse gives no hit. A ≤ 0 is dropped. Clipped channels with fewer than 6 unclipped samples fall back to the C++ peaks (`mf_pulse.py:1197-1214`, B:1851-1860). The example split has 6 saturated channels and 1 fallback.

**Legacy vs gated.** The pulse list is the same; the paths diverge in S5.

### S4. Time: pre-window, no offsets, LMS

**What.** Pulse times are raw, in ns from the start of the readout window.

- The only front-end time cut is the open interval 100 < t < 200 ns, applied to every pulse before either clusterer sees it.
- `--window 16.0` is **not** a pulse-time window. It is only the seed-time gate of the gated clusterer (B:2720; see S6). The legacy path has no ±16 ns cut.
- No per-module time offsets are subtracted.
- The LMS gain factor is 1 with `--lms-corr off`.

**Where.**
- Offsets: `resolve_time_calib_path` B:1113-1122 and `load_time_offsets` B:1125-1140.
- LMS: `resolve_lms_path` B:1100-1110; `off` gives one batch of ones (B:1194-1215).
- The cut is in `_event_hits`, B:1616-1621.
```
toff = t_off[mod.index] if t_off else 0.0          # B:1616
...
t = tm - toff
if t <= pre_lo or t >= pre_hi:                      # B:1620-1621
    continue
```

**R5nr settings.** `--pre-window 100.0 200.0`, `--time-calib off` (the log says "off (0 modules with offset)"), `--lms-corr off` for 3.49 and 0.73. 2.24 H2 run 25320 and C run 25308 instead use `output/lms_premerge/prad_0253{08,20}_gain_corr.root` (16 jobs in `final_nr.txt`).

**Decisions.**
- MF time offsets were rejected: σ_MAD 0.78/0.81 ns with MF offsets versus 0.67/0.71 ns with rawsum offsets (`notes/PbWO4_energy_resolution_20261022.md` §6).
- LMS on versus off is a null result at 3.49 GeV: CB 1.619 → 1.619 % (`notes/PbWO4_energy_resolution_summary.md:127`).
- The configuration key `time_calib: true` has no effect on B (P4).

### S5. Per-module energy and the legacy/gated hit lists

**What.** Each in-window pulse becomes a module hit with E_MeV = cal_factor·E_u·f_LMS, via `Module::energize` (zero if negative). cal_factor comes from `--hc-calib-file`.

- **Legacy:** one hit per module, the in-window pulse with the largest E_u. The variable is named `best_int`, but it holds E_u.
- **Gated:** every in-window pulse is pushed as its own hit.

The per-crystal gains are **not** applied here (see S10), so all clustering thresholds act on energies before the gains.

**Where.** `MFBlock._event_hits` is B:1602-1627; the rawsum equivalent is B:3008-3027. `Module::energize` is `prad2det/include/HyCalSystem.h:154-158`. The calibration is loaded by `HyCalSystem::LoadCalibration` (`prad2det/src/HyCalSystem.cpp:504-529`).
```
for tm, eu in plist:                                       # B:1618-1626
    t = tm - toff
    if t <= pre_lo or t >= pre_hi:
        continue
    hn.append((mod.index, mod.energize(eu * f), float(t)))   # gated: every pulse
    if eu > best_int:
        best_int, best_t, found = eu, t, True
if found:
    hl.append((mod.index, mod.energize(best_int * f), float(best_t)))   # legacy: one per module
```

**Settings and inputs.** `--hc-calib-file` per point (P4). At 3.49 GeV, 1152 PbWO4 modules have a factor > 0.

### S6. Island clustering (C++ HyCalCluster)

**What.**
1. Hits above `min_module_energy` from live modules enter the clusterer.
2. Seeds are taken in descending energy and need E ≥ `min_center_energy`.
3. An island grows by BFS over edge neighbours (`corner_conn` false). In the gated path, neighbour pulses must lie within ±`seed_time_window` of the seed.
4. Islands with several local maxima are split. Starting fractions come from the shower profile; then `split_iter` iterations recompute centroids and profile fractions, and fractions below `least_split` are dropped.
5. Clusters must pass `min_cluster_energy` and `min_cluster_size`.

**Where.**

| Function | `HC.cpp` lines |
|---|---|
| `AddHit` | 94-100 |
| `FormClusters` | 102-117 |
| `ReconstructMatched` | 134-146 |
| `group_hits` | 148-182 |
| `grow_island` | 184-219 (time gate 187-206) |
| `split_cluster` | 221-244 |
| `find_maxima` | 246-279 |
| `split_hits` | 281-328 |
| `eval_fraction` | 348-409 |

Configuration defaults are `HC.h:35-79`; the Geant4 profile is `HC.h:257-291`; the binding is `python/bind_det.cpp:766-790`. The Python builder `_make_cluster` is B:2696-2723:
```
base_cfg = p.hc_clusterer.get_config()
def _make_cluster(seed_window):
    cfg = det.HyCalClusterConfig()
    for f in ("min_module_energy", "min_center_energy", "min_cluster_energy",
              "min_cluster_size", "corner_conn", "split_iter",
              "least_split", "log_weight_thres"):
        setattr(cfg, f, getattr(base_cfg, f))
    cfg.seed_time_window = seed_window
    ...
    if args.cluster_profile == 'geant4':
        cl.set_config(base_cfg)    # installs base_cfg.profile
    cl.set_config(cfg)
# cl_legacy = _make_cluster(-1.0) ; cl_new = _make_cluster(W)   (B:2719-2720)
```

**Effective R5nr cluster configuration.** Verified from the `[merge_cfg_wrap]` log line (L30) of the example split and from the code:

| Field | Value | Source |
|---|---|---|
| min_module_energy | 1.0 MeV | configuration file, copied |
| min_center_energy (seed) | **30.0 MeV** | r5 configuration via the wrapper (database: 10) |
| min_cluster_energy / min_cluster_size | 50.0 MeV / 1 | copied |
| corner_conn | false | copied |
| split_iter / least_split | 6 / 0.01 | copied |
| log_weight_thres | 5.2 | copied |
| seed_time_window | −1 (legacy) / 16 ns (gated, `--window`) | B:2719-2720; the configuration's 4.0 is overwritten |
| profile | Geant4 tables, 200-2100 MeV (clamped above) | `--cluster-profile geant4` |
| non_linear_corr / leakage_correction / energy_bias_correction | true (inert) / **false** / false | C++ struct defaults; not copied |

**Caveats.**
- Log line L15 prints `min_ctr_E=10` because it is written before the wrapper patch takes effect.
- L30 prints `leakage_correction=True`. That value is on `p.hc_clusterer`, which B uses only as the configuration source, so the benchmark clusterers have leakage **off**.
- B:2724-2751 self-checks that the Geant4 and simple profiles split a test row differently (L42).
- Applying the 2100 MeV profile row to 3.49 GeV clusters is an inference from the clamp; it was not checked numerically.

### S7. Cluster-energy corrections in C++ (all effectively off)

**What.** `reconstruct_pos` multiplies the cluster energy by `bias_corr` and `linear_corr`, where 1/linear_corr = 1 + nl1·ΔE + nl2·ΔE² and ΔE is in GeV relative to the seed's `base_energy`.

- **Non-linear correction is inert.** The calibration files carry only the legacy key `non_linear` (= 0), while `LoadCalibration` reads `nl1`/`nl2` (`HyCalSystem.cpp:524-525`). linear_corr is therefore 1. The log confirms "member identity … 0 violations", with maximum relative deviation 6.7e-7.
- Leakage correction (`HC.cpp:485-570`) and energy bias are off (S6).
- `shower_depth` (`HC.cpp:59-65`) is never used by B.

**Where.** `HC.cpp:412-476` (`reconstruct_pos`), formula at :458-473.
```
const float delta_gev = (cl.energy - base_energy) / 1000.f;
float non_linear_factor = 1.f / (1.f + nl1 * delta_gev + nl2 * delta_gev * delta_gev);
if (cl.energy > 3800.f || non_linear_factor < 0.7f || non_linear_factor > 1.3f) non_linear_factor = 1.f;
result.energy = cl.energy * result.bias_corr * result.linear_corr;
```

**Evidence for the decision.** C++ energy-bias correction was tested as R6 arm "r6rd". It was worse than R5nr: +0.051 ± 0.003 pp at 3.49 GeV and +0.026 ± 0.004 pp at 2.24 GeV (`output/r6_20261030/tables/promotion_r6.md`).

### S8. Position reconstruction

**What.** The cluster position is a log-weighted centroid over the seed and its 3×3 neighbours, with w_i = max(0, W0 + ln(E_i/E_cl)) and W0 = 5.2. Offsets are in module units, scaled by the seed's size. There is no depth correction and no lab transform; x and y are HyCal-local mm.

Downstream, B uses the position in three ways:
- θ = atan(r/hycal_z), with the target at z = 0 (`module_theta_deg`, B:290-293).
- The window pass uses the seed-module centre for θ (B:874-875).
- The Møller tag and the position correction use the cluster centroid.

**Where.** `HC.cpp:22-55` (`LogCentroid`) and :478-483 (`get_weight`). The centroid is stored at tuple positions 10 and 11 of each event (B:2795, 2891-2916).
```
float w = config_.log_weight_thres + std::log(E / E_total);   // HC.cpp get_weight
return (w > 0.f) ? w : 0.f;
x = center.x + (wx / wtot) * center.size_x;                   // HC.cpp:48
```

**Decision.** `lw3p6` (W0 = 3.6 with the R3 maps) as an ablation: ΔCB = −0.041 ± 0.005 pp at 3.49 GeV (`$O/tables/ablation_nr.md`).

### S9. Event selection and beam energy

**What.**
1. An event enters the per-path event list if the clusterer returns **exactly one cluster above `signal_min` = 500·kE MeV**, where kE = Ebeam/2239.51 and kS = √kE.
2. Ebeam and hycal_z come from runinfo (P4).
3. The headline pool, `_select_pool`, requires:
   - seed-module centre at r < 200 mm,
   - seed energy ≥ 200·kE,
   - seed **outside** the central 6×6 block, rows and columns 15-20 (1-based, inclusive) around the beam hole.
4. Finally, only seeds that pass the elastic-window selection enter the spectrum (S14).

No second-cluster veto is applied (`--veto-second-cluster` is None).

**Where.** Gate in `finish_event`, B:2859-2923 (cut at :2872, append at :2881-2899). `_select_pool` is B:615-632. The scaling defaults are B:2564-2579.
```
big_l = [r for r in rec_l if r.hit.energy > args.signal_min]          # B:2872
if (len(big_l) == 1 and not (veto_on and veto_other(...))):
    events_legacy.append((r.hit.center_id, sm.x, sm.y, sm.row, sm.column,
                          r.cluster.center.energy, r.hit.energy, ...) + member_rec(...))
# _select_pool: x^2+y^2 < inner_r^2, seed_E >= min_seed_E,
#   not (excl[0] <= row+1 <= excl[1] and excl[2] <= col+1 <= excl[3])
```

**R5nr settings.** `--inner-radius 200`, `--exclude-rowcol 15 20 15 20`, `--min-per-module 30` (used in S14).

At 3.49 GeV: kE = 1.55768, signal_min = 778.84 MeV, min_seed = 311.54 MeV. Example split: 122692 legacy single-cluster events; 32867 in the pool.

**Caveat.** `notes/PbWO4_energy_resolution_summary.md` §1.4 calls this cut a "2-ring edge exclusion". It is the central-hole block, applied to the seed only.

### S10. Per-crystal gains applied to the cluster energy (`--gain-file`)

**What.** Each cluster energy in both paths' event lists is replaced by E_cl = Σ_i g_i E_i over the cluster's members, using their split-shared energies. Modules missing from the table get the table's g0.

- This happens after clustering. The seed threshold, `signal_min` and the pool cuts still act on energies before the gains.
- It bypasses `linear_corr`, which is inert anyway.
- **The headline CB spectrum includes these gains**, because S15 histograms `e[6]` of these same lists.

**Where.** B:3225-3239.
```
gtab, gtab_g0 = PCC.read_gain_table(args.gain_file)
garr = np.full(max(p.hycal.module_count(), max(gtab) + 1), gtab_g0)
for i_, g_ in gtab.items():
    garr[i_] = g_
for evs_ in (events_legacy, events_new):
    for k_, t_ in enumerate(evs_):
        evs_[k_] = (t_[:6] + (float(np.dot(garr[t_[8]], t_[9])), t_[7]) + t_[8:])
```

**Input.** `$O/gains/<run>/gains_<run>_all.tsv` (P2).

**Caveat.** The `--gain-file` help text (B:2231-2239) says "the headline spectrum [is] unchanged". "Headline spectrum" there means B's whole-spectrum table (all radii, the first block of `result.tsv`), not the inner-ring CB σ that this note calls the headline.

### S11. Position (bowl) correction (`--pos-corr`)

**What.** Each cluster energy is divided by R(u,v), with u = (x_cl − x_seed)/size_x and v defined the same way. The correction is applied after S10 and before the Møller tag.

- **Hybrid map:** a seed with its own fit uses the poly8 R = c0(1 + c1u² + c2v² + c3u²v² + c4u⁴ + c5v⁴ + c6u + c7v), with u and v clipped to ±0.6.
- **Other seeds:** bilinear interpolation on the 12×12 fallback grid.
- Events with non-finite u, v or R are left unchanged.

**Where.** `poscorr_load` is B:653-680, `poscorr_eval` B:689-744, and the application B:3241-3280.
```
u_ = (np.array([e_[10] - e_[1] for e_ in evs_], dtype=np.float64) / sz[:, 0])     # B:3255
...
R_, clp_ = poscorr_eval(pcd, u_, v_, seed_)                                         # hybrid
...
evs_[k_] = t_[:6] + (t_[6] / float(R_[k_]),) + t_[7:]                               # E <- E / R
```

**Input.** `$O/poscorr/maps_norefine/3p49GeV_H2/poscorr_r5nr_3p49GeV_H2_all.json` (P3). Example split: applied to 245380 events (both paths); 172 clamped; 46868 with their module's own fit.

**Value.** See the budget list (R6 nopos arm).

### S12. Møller kinematic rejection

**What.** A single-arm e-e (Møller) electron can sit under the e-p peak. A cluster is tagged and rejected if its energy is near the Møller line and away from the e-p line:

  |E − E_ee(θ)| < n·σ  and  |E − E'(θ)| > m·σ,

with n = 3, m = 2, σ = res·E′ and res = √((0.025/√E′_GeV)² + 0.007²). The definitions:
- E is the corrected cluster energy after S10 and S11.
- θ comes from the cluster centroid.
- E′ is defined in S13.

The filtered lists `evf_*` feed the headline pool. `--moller-reject-gains on` also feeds them to the window pass. The `pair` tag (two-cluster sum and acoplanarity) is written out but not used.

**Where.** `moller_kin` is B:3290-3302; constants B:253-254; filtering B:3308-3316.
```
th_ = math.degrees(math.atan(math.hypot(x_cl, y_cl) / hz_mr))
ep_ = expected_energy_elastic(th_, ebeam_mr, m_target)
ee_ = expected_energy_ee(th_, ebeam_mr)
res_ = (args.moller_kin_res if args.moller_kin_res is not None
        else math.sqrt((MOLLER_RES_A / math.sqrt(ep_ / 1000.0)) ** 2 + MOLLER_RES_B ** 2))
sig_ = res_ * ep_
return int(abs(e_cl - ee_) < args.moller_kin_n * sig_
           and abs(e_cl - ep_) > args.moller_kin_m * sig_)
```

**R5nr settings.** `--moller-reject kin --moller-kin-n 3 --moller-kin-m 2 --moller-reject-gains on`. Example split: 18495 of 122692 legacy events are kin-tagged, but only 563 of the 32867 pool events are rejected.

**Caveat.** The rule over-subtracts at 0.73 GeV; see the budget list.

### S13. Expected energy E'(θ)

**What.** E'(θ) is the e-p elastic scattered energy for a proton at rest, minus a constant material energy-loss term:

  E′ = E·M/(M + E(1 − cosθ)) − 2.056 MeV / cosθ, with M = 938.272 MeV.

θ = atan(r/hycal_z), with no shower depth and no beam offset. In R5nr, E' is used for the Møller tag (S12), the window-pass targets and the histogram centre `target_ref` (S14-S15). **It does not rescale the fitted energies**: the fitted variable is the absolute energy (S15). Across the pool radii E′ varies by only 0.17 % at 3.49 GeV and 0.03 % at 0.73 GeV (inference, by hand arithmetic).

**Where.** `energy_loss` B:257-261; `expected_energy_elastic` B:264-270; `expected_energy_ee` B:278-287; `module_theta_deg` B:290-293; `M_PROTON` B:225.
```
e = ebeam * m_target / (m_target + ebeam * (1.0 - cos_t)); return e - energy_loss(theta_deg)
energy_loss = (0.5*1.6 + 0.12*1.6 + 0.1*2.0 + 0.48*1.8) * sec     # = 2.056 MeV / cos(theta)
```

### S14. `--calibrate`: elastic-window pass with refinement off

**What.** `--calibrate` fills the event lists and, after the event loop, runs one elastic-window pass per seed module (`fit_elastic_window_gains`). The pass has these steps:
1. A global fit gives μ_global and σ_glob.
2. Each module gets a CB pre-fit, then a window μ_M ± 2σ_M.
3. A module is dropped if its e-p/Møller gap is under 3σ_M, if it has fewer than 30 in-window events, or if it has no pre-fit events.
4. gain = E′(θ_M)/median(E in window).

With `--gain-refine off` these gains are measured and written to `result_window.tsv`, but every used module's gain is set to 1.0, so the "post" spectrum equals the "pre" spectrum. **No step inside a job is iterative**, apart from the Gaussian refit loop of S16 and `split_iter` in clustering.

**Where.** `fit_elastic_window_gains` B:823-955 (gain at :950-951); the call is B:3380-3386; the headline pool with Møller rejection is B:3387-3400; refine off is B:3463-3470.
```
if args.gain_refine == 'off' and not pcm:
    gains_l = {sid: 1.0 for sid in gains_l}
    gains_n = {sid: 1.0 for sid in gains_n}
    corr_l, corr_n = np.array(raw_l, copy=True), np.array(raw_n, copy=True)
# headline sample, B:3395-3397:
raw_l = np.array([e[6] for e in pool_l if e[0] in gains_l])
corr_l = np.array([e[6] * gains_l[e[0]] for e in pool_l if e[0] in gains_l])
```

**R5nr settings.** `--calibrate`, `--cal-method elastic-window` (default), `--window-nsigma 2` (default), `--min-per-module 30`, `--gain-refine off`.

The calibration actually applied in R5nr is therefore cal_factor (S5) × per-crystal g_i (S10) ÷ R(u,v) (S11).

**Caveats.**
- "post = pre" holds exactly for the energies, but the window pass still decides **which seed modules** enter the spectrum: 138 modules and 29712 legacy events in the example split. That selection includes the 3σ Møller-gap exclusion.
- The decision to turn refinement off is recorded in `plans/plan_r5_20261028.md:48-64` (contract change 2026-10-29). Refined R5 was 0.02-0.04 pp better on e-p (`$O/tables/r5nr_table.md`), but refinement is in-sample.

### S15. The fitted spectrum

**What.** The fitted spectrum is a histogram of the **absolute cluster energy in MeV**, not E/E′ or E − E′. It covers [target_ref − 800·kE, target_ref + 800·kE], where target_ref is the mean E′ of the pool, with bin width 20·kS. There are four spectra (legacy/gated × pre/post); with refinement off, post equals pre. σ_E/E is σ/μ of the fit (S16).

At 3.49 GeV the range is 3482.6 ± 1246 MeV in 24.96 MeV bins.

**Where.** B:3472-3479.

**Caveat.** The first table in `result.tsv` (the whole-spectrum Gaussian, all radii, 500-3500·kE) is a different quantity. Do not mix it with the headline.

### S16. Fits: fixed Gaussian → CB core (headline); ±2σ Gaussian (secondary)

**What.**
1. **Fixed-window Gaussian** (`fit_peak`, B:983-992): ±max(150·kS, 4·bw) around the argmax. It supplies μ0 for the other fits.
2. **Reported Gaussian** (`fit_gauss_2sigma_hist`, B:1016-1067): iterative ±2σ refit, at most 5 refits, 1 % tolerance. This is the "G" column.
3. **CB core.** A Crystal Ball fitted on [μ0 − w, μ0 + w/2] with **w = max(400·kS, 12·bw)**, where kS = √(Ebeam/2239.51) and bw is the bin width. With the default bw, w = 400·kS: 499 MeV at 3.49 GeV and 228 MeV at 0.73 GeV. Bounds: α 0.1-20, n 1.01-50, σ 1-w. The fit needs at least 50 entries. In R5nr the CB used for the headline comes from `fit_cal_cb` (B:348-508), which with `--cb-two-peak` given evaluates mode `two-peak`.
4. **Two-peak CB.** With `--moller-reject` on it immediately falls back to the plain single CB (`single_reason=moller_reject`, B:414-417). Two-peak CB therefore **never runs in R5nr** at any energy. The `result.tsv` header reads `cb_fit mode=two-peak … legacy_post_cb_mode=single legacy_post_single_reason=moller_reject`.
5. **`--fit-select chi2`** (`select_fit_chi2`, B:511-562) compares Gaussian and CB χ²/ndf over μ_G ± 3σ_G and writes `preferred_fit`. It is informational only; the headline is always CB.

**Where.** `fit_peak_cb` B:296-328 (window at :311-312); `fit_cal_cb` B:348-508 (window at :409); CB shape B:211-222.
```
win  = max(400.0 * escale, 12.0 * bin_width)                       # B:311
mask = (centres > mu0 - win) & (centres < mu0 + 0.5 * win)         # B:312
...
reason = ('moller_reject' if moller_reject else ...)                # B:414
if reason:
    out['single_reason'] = reason
    _plain()
```

**Caveat.** At 3.49 GeV the single-CB χ²/ndf is 7.97. For the CB-window study, see the budget list.

### S17. Per-split outputs

**`result.tsv`** (writer B:3787-4025)
- Header lines: resolved Ebeam/kE/kS, `gain_file=`, `pos_corr=`, `moller_reject` counts, `cb_fit` mode.
- The whole-spectrum table.
- The **`# --calibrate inner ring`** table, with columns `path n_modules n_inner_events mu_pre_MeV sigma_pre_MeV mu_post_MeV sigma_post_MeV sigma_over_E_post mu_cb_MeV sigma_cb_MeV alpha_cb n_cb sigma_over_E_cb …`.
- The headline is row `legacy`, column `sigma_over_E_cb` (a fraction). Example split (`result.tsv:222`):
```
legacy  138  29712  3493.23 54.66 3493.23 54.66 0.0156 3494.00 53.32 1.258 2.718 0.0153 ...
```
  This is 1.53 % (legacy); the gated row 223 is 1.52 %.

**`result_window.tsv`**: the per-module window table, with measured gains not applied.

**`result_clusters.tsv`**: one row per accepted event. `cluster_E_raw` already includes S10 and S11; `cluster_E_cal` is raw × 1.0, or NaN for unused seeds.

---

## Part C. Aggregation

### S18. Paired split aggregation (mean ± SE over splits)

**What.** `compare_scan_dirs.py` reads every `scan_*/split_*/result.tsv` of each scan and keeps only the (run, split) keys present in every scan ("paired").

- CB σ/E is recomputed as `sigma_cb_MeV/mu_cb_MeV` from the inner-ring row, with no refit.
- G is refit from the stored histogram with `--gauss-window 2sigma`.

`nr_table.py` then takes the mean and SD/√n (ddof = 1) over splits, per energy point (per run for Ta). It flags "CB≠G" when |mean(Δ_CB − Δ_G)| > 0.05 pp or > 2 SE.

**Where.**
- `compare_scan_dirs.py`: `load_scan` 73-93, `split_rows` 376-396, Gaussian refit 109-153.
- `gate_scan_20261022.py`: `sigma` 535-543, `_SIGCOLS` 529-532.
- `$NR/tools/nr_table.py`: `ms` 39-44, flag 102-104.
- The driver is `$O/tools/s6_nr.sh:13-15, 21`.
```
def ms(v):                                       # nr_table.py:39-44
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    ...
    return v.mean(), (v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else np.nan), len(v)
```

**Commands** (from `s6_nr.sh`):
- `compare_scan_dirs.py --scans R5=$O/scan/r5_splits R5nr=$O/scan/r5nr_splits --gauss-window 2sigma --out $O/compare/r5nr_vs_r5`
- `nr_table.py --ep … --ta … --ref $O/compare/r5nr_vs_r5_splits.tsv --out $O/tables/r5nr_table`

**Splits per point** (`final_nr.txt`):

| Point | Splits | Runs (splits) |
|---|---|---|
| 3.49 GeV H2 | 27 | 24386 (10), 24503 (10), 24512 (7) |
| 2.24 GeV H2 | 42 | 25228 (11), 25257 (7), 25260 (8), 25299 (8), 25320 (8) |
| 0.73 GeV H2 | 16 | 24590 (10), 24603 (6) |

The 41 reference splits (4 each of 24503, 24590, 25308, 25320, plus all 25 Ta) carry the ablations, so the ablation table has n = 4 per e-p point.

### S19. a/√E ⊕ b fit

**What.** The fit is weighted least squares of σ² (in %²) against 1/E on the three H2 means, σ² = p0/E + p1. The weights are 1/(2σ·SE)², with nominal E = 0.727, 2.236 and 3.484 GeV, giving 1 degree of freedom. The results are reported as a = √p0 and **b = sign(p1)·√|p1|**. A negative b therefore means a *negative constant term* p1 = −b², which is unphysical; it is not a fitted negative b. A 3-term variant (adding c/E) is printed alongside it.

**Where.** `$O/fit/phase0.py:125-151` (`fit2`, `do_fit`; b at :140). It is duplicated in `$NR/tools/nr_table.py:47-59`. The driver is `s6_nr.sh:19, 21`.
```
s2 = ys ** 2; err = 2 * ys * ses
A = np.c_[1 / Es, np.ones_like(Es)]
W = 1 / err ** 2
cov = np.linalg.inv(A.T @ (A * W[:, None])); p = cov @ (A.T @ (W * s2))
chi2 = float(np.sum(W * (s2 - A @ p) ** 2))
a = np.sqrt(p[0]); b = np.sqrt(abs(p[1])) * np.sign(p[1])
```

**Why b < 0 and χ² = 203.** The 0.73 GeV point (3.416 %) lies above every a/√E ⊕ b curve through the 2.24 and 3.49 GeV points. A recomputation from the table inputs gives p = (8.604, −0.5006), a = 2.933, b = −0.7075, χ² = 203.7, which matches the table to input rounding.

---

## Current R5nr numbers

Verified from `$O/tables/r5nr_table.md`, CB core σ_E/E in %, mean ± SE over paired splits, with [G] the ±2σ Gaussian:

| Point | n | Legacy CB (headline) | Gated CB | Legacy G |
|---|---|---|---|---|
| 3.49 GeV H2 | 27 | **1.494 ± 0.008** | 1.493 ± 0.008 | 1.541 |
| 2.24 GeV H2 | 42 | **1.800 ± 0.004** | 1.798 ± 0.004 | 1.845 |
| 0.73 GeV H2 | 16 | **3.416 ± 0.009** | 3.415 ± 0.009 | 3.450 |

The 2.24 GeV C point (out of scope) is 1.799 ± 0.009.

**a/√E ⊕ b, three H2 points:**

| Path | Series | a (%·√GeV) | b (%) | χ² (1 dof) |
|---|---|---|---|---|
| legacy | CB | **2.934 ± 0.011** | **−0.709 ± 0.021** | 202.99 |
| legacy | G | 2.858 ± 0.011 | −0.347 ± 0.039 | 331.67 |
| gated | CB | 2.935 ± 0.011 | −0.713 ± 0.021 | 205.30 |

The same fit appears in `$O/tables/r3_r4_r5nr.md`, where R5nr is relabelled "R5" by `s6_nr.sh:29`, and in `$O/fit/out_r5.txt`. No conflicting a/b values were found.

## Where the resolution budget goes (existing studies)

- **Integration range: moot.** `int_tail_ratio` is a per-module constant absorbed by the gains; the pedestal floor moves σ by ≤ 0.0005 pp. Sources: `notes/PbWO4_energy_resolution_20261030.md` §2, `$O/intscan/tables/intscan.md`.
- **CB window: ≤ 0.023 pp on σ** (largest at C), although χ²/ndf at 0.73 GeV drops from 186.67 to 4.84. The patch was not applied. Sources: `$O/cbrange/tables/cb_range.md`, notes 20261030 §3.
- **Møller over-subtraction at 0.73 GeV:** the rule's σ (≈3.0 %) is below the measured 3.45-3.68 %, so the cut removes 5.51 % of elastic events and widens CB σ by +0.168 ± 0.007 pp (`nomoller` − R5nr: ΔCB −0.156 ± 0.011, ΔG +0.234 ± 0.017). Sources: notes 20261030 §4, `output/moller_oversub_20261029/`, `$O/tables/ablation_nr.md`.
- **Position correction is worth** (nopos − R5nr, legacy CB) +0.141 ± 0.015 pp at 3.49, +0.053 ± 0.005 at 2.24 H2 and +0.015 ± 0.002 at 0.73. Sources: `output/r6_20261030/tables/promotion_r6.md:118`, `plans/plan_r6_20261030.md`.
- **Gate scans:** no global gate W beats legacy on e-p by more than 0.004 pp, and tight gates hurt (for example +0.45 pp). The 16 ns gate matters only for Ta. Sources: `notes/PbWO4_energy_resolution_summary.md:136-138`, `plans/plan_20261021.md` §8, `gate_scan_20261021.py`, `gate_scan_20261022.py`.
- **Other ablations** (variant − R5nr, CB pp, n = 4, at 3.49 / 2.24 / 0.73): seed 10 MeV +0.009 / +0.012 / +0.016; lw 3.6 −0.041 / −0.015 / +0.008; simple profile 0.000 / −0.004 / +0.003; R4 gains +0.003 / −0.008 / +0.025; legacy fit 0.000. Source: `$O/tables/ablation_nr.md`.

## Is anything iterative?

- **Within a split job: no.** The only loops are `split_iter` (6) in clustering and the ±2σ Gaussian refit (at most 5). The window pass is a single pass, and its gains are not applied.
- **Offline, in sequence:**
  1. P1 templates (independent).
  2. P2 gains, pass 0 → pass 1 (Møller-rejected), each with an internal `--fixed-windows` loop that converges in 2 passes.
  3. P3 maps, built *with* the P2 gains.
  4. The final R5nr scan, with gains and maps.

  Nothing feeds back from the final scan to the gains or maps. Gains and maps (and probably the templates) are **in-sample** with respect to the evaluated splits.

## Discrepancies found (notes vs code at HEAD)

1. **Two-peak CB never runs in R5nr.** `--moller-reject` forces the single CB (B:414; `result.tsv` header `legacy_post_cb_mode=single legacy_post_single_reason=moller_reject`). The summary §3.12 statement "two-peak = single CB, Δ +0.000" is true by construction, which the note does not say.
2. **`--window 16` is not a "±16 ns raw time window".** It is the seed-time gate of the *gated* clusterer only. Both paths cut pulses only on the raw pre-window (100, 200) ns.
3. **The `--gain-file` help text** (B:2231-2239) says the headline spectrum is unchanged. That refers to the whole-spectrum table; the inner-ring CB headline includes the gains (S10) and the position map (S11).
4. **The fitted variable is absolute E in MeV, with σ/μ_CB**, not an E/E′ or E − E′ distribution. E′ is used only for selection and the histogram centre.
5. **Stale constants.** `notes/PbWO4_methods_summary_20261027.md` §4 and `notes/PbWO4_energy_resolution_summary.md` §3.4 list log-weight 3.6, seed 10 MeV and SimpleProfile. Those are R3/R4 values; R5nr uses 5.2, 30 MeV and Geant4.
6. **Calibration file.** methods_summary §3 names only `calibration_factor_2p2_fadc.json`. 3.49 GeV uses `calibration_factor_fadc_3p5.json` and 0.73 GeV uses `calibration_factor_fadc_0p7.json`.
7. **Stale line numbers.** The `B:`, `HC:`, `WA:` and `MF:` line numbers in methods_summary §2, §4 and §6 are stale at HEAD, as is the `mf_pulse.py` docstring reference "WaveAnalyzer.cpp:391-462". The `STE:` references are still correct.
8. **"2-ring edge exclusion"** (summary §1.4) is the central 6×6 hole block, applied to the seed only.
9. **Legacy hit choice.** The legacy path takes the pulse with the largest E_u (MF energy), not the largest raw integral (B:1622-1626).

Code- and configuration-internal inconsistencies are flagged as caveats in their steps: the wrapper log line (S6), `time_calib: true` (P4), the duplicate `--mf-scale` (S0, S3), LMS on for 25308/25320 (S4), the stale `r5_jobs.py:132` comment (P3), the meaning of `n_rejected` and the d91af9b7/08b19869 split (P2), and the copied `job.cmd` paths (S0).

## Open questions and coverage gaps

- See also the inference-labelled caveats in P1, P3 and S6.
- The count of pool modules that could contribute, versus the 138 used, was not derived.
- The +0.0078 ± 0.0047 pp result for plan_r5 decision 2 (R4 gains vs R5 gains) was not re-verified.
