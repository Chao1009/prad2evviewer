# e-p resolution walkthrough: scripts per step (2026-11-02)

Companion to `notes/ep_resolution_walkthrough_20261102.md`. Lists the scripts and functions involved in each step of the R5nr e-p energy-resolution chain (R5nr: commit 71e049f, benchmark sha 08b19869).

Conventions:
- Paths are relative to `analysis/pyscripts/` unless noted; `prad2det/`, `prad2dec/`, `python/` and `database/` are relative to the repository root.
- **B** = `benchmark_hycal_timing_test_20261021.py`.
- Line numbers are as of 2026-11-02 and drift with edits; function names are the stable reference.

## Preparation (inputs made before the evaluation pass)

| Step | What it does | Scripts and key functions |
|---|---|---|
| P1 | Spline templates, I_int, amplitude correction f(A) | `spline_template_extract.py` (`Selector` :501, `cmd_select` :1113, `build_shape` :1964, `cmd_build` :2826, `ampcorr_rule` :3503, `cmd_ampcorr` :4051); `mf_pulse.py` (`SampledTemplateBank.build_iint` :600, `iint` :636) |
| P2 | Per-crystal gains (joint ridge least squares, Møller-rejected) | `output/r5_20261028/tools/r5_jobs.py` (dump jobs, :124-127), `fit1.sh`; `percrystal_calib.py` (`_normal_equations` :90, `solve` :114, `write_gain_table` :298, `read_gain_table` :358, `moller_kin_mask` :560, `cmd_fit` :589) |
| P3 | Position (bowl) correction maps, hybrid poly8/bilinear | `r5_jobs.py` (map jobs, :131-135), `norefine_scan.sh`; `position_corr.py` (`load_pool` :114, `cmd_map` :461, `cmd_fit_hybrid` :930, `cmd_fit` :1001) |
| P4 | Fixed inputs | `database/calibration/old/calibration_factor_{fadc_3p5,2p2_fadc,fadc_0p7}.json`, `x17_*_fadc.json`; `runinfo/general.json`; `daq_config.json`; `cluster_profiles/prof_pwo.dat`, `prof_lg.dat` (loaded in `PipelineBuilder.cpp` :275-278); `reconstruction_config_r5.json`, `reconstruction_config.json`; `database/hycal_time_offsets/time_calib.json` |

Note: the starting bank `output/templates_spline/3p49GeV_H2/templates_3p49GeV_H2.npz` for P1 is missing; the best-50 bank used by R5nr exists.

## Evaluation pass

| Step | What it does | Scripts and key functions |
|---|---|---|
| S0 | Orchestration | `r5_jobs.py`, `chain_r5nr.sh`, `run1_sha.sh`, `run1.sh`; `output/merge_check_20261028/tools/merge_cfg_wrap.py`; B `main` :2021 |
| S1 | Decode, trigger selection | B :2942-2975; `_common.py` (`setup_pipeline` :202) |
| S2 | Pedestal, WaveAnalyzer peaks | `prad2dec/src/WaveAnalyzer.cpp` (`smooth` :244, `findPedestal` :272, `findPeaks` :358, `Analyze` :587); `python/bind_dec.cpp`; `daq_config.json` :100-113 |
| S3 | Matched-filter fit (kmax 4) | `mf_pulse.py` (`MFReconstructor` :833, `run` :1165, `pulse_energy_units` :1245); B `MFBlock.flush` :1829 |
| S4 | Hit time, LMS gain correction | B `_event_hits` :1602, `resolve_lms_path` :1100, `resolve_time_calib_path` :1113, `load_time_offsets` :1125; `output/lms_premerge/*.root` |
| S5 | Hit energy in MeV | B `_event_hits`; `prad2det/include/HyCalSystem.h` (`energize` :153); `prad2det/src/HyCalSystem.cpp` (`LoadCalibration` :504) |
| S6 | Island clustering | `prad2det/src/HyCalCluster.cpp` (`FormClusters` :102, `grow_island` :184, `split_cluster` :221); `python/bind_det.cpp`; B `_make_cluster` :2697 |
| S7 | Cluster corrections (leakage, energy bias; off in R5nr) | `HyCalCluster.cpp` (leakage block); `HyCalSystem.cpp` :524-525 |
| S8 | Cluster position, angle | `HyCalCluster.cpp` (`LogCentroid` :31); B `module_theta_deg` :290 |
| S9 | Event selection, cluster pool | B `finish_event` :2859, `_select_pool` :615 |
| S10 | Per-crystal gains applied | B :3225-3239; `percrystal_calib.py` (`read_gain_table`) |
| S11 | Position correction applied | B `poscorr_load` :653, `poscorr_eval` :689, :3241-3280 |
| S12 | Møller kinematic rejection | B `moller_kin` :3290; `percrystal_calib.py` (`moller_kin_mask`) |
| S13 | Expected elastic energy E′(θ) | B `energy_loss` :257, `expected_energy_elastic` :264 |
| S14 | Elastic-window calibration pass | B `fit_elastic_window_gains` :823 |
| S15 | Histogram fill | B :3472-3479 |
| S16 | Peak fits (CB core headline, ±2σ Gaussian secondary) | B `fit_peak` :983, `fit_gauss_2sigma_hist` :1016, `fit_peak_cb` :296, `fit_cal_cb` :348, `select_fit_chi2` :511 |
| S17 | Result writer (`result.tsv`) | B :3787-4025 |
| S18 | Aggregation over paired splits | `compare_scan_dirs.py`; `output/resolution_fit_20261029_r5/tools/nr_table.py`; `s6_nr.sh` |
| S19 | a/√E ⊕ b fit | `output/r5_20261028/fit/phase0.py` (`fit2` :125, `do_fit` :136); `nr_table.py` |
