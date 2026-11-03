# Plan: spline-template MF pulse extraction under the standard HyCal energy calibration (2026-11-02)

**Status:** draft for user approval. Drafted by planner 2026-11-02. **Branch:** `waveform_jz` (HEAD 71e049f).
**User request (verbatim):** "docs/analysis_notes/replay_functions_and_tools.md, if we want to adapt the current spline fit and Matched filter extraction to the standard energy calibration introduced in the note, how do we do it?"
**Inputs read:** `docs/analysis_notes/replay_functions_and_tools.md` (RF below), `docs/analysis_notes/physics_calibration.md` (PC), `docs/analysis_notes/gain_correction.md` (GC), `notes/ep_resolution_walkthrough_20261102.md` (WT), `plans/plan_r6_20261030.md`, and the code cited below.

Paths are relative to `analysis/pyscripts/` unless they start with `docs/`, `database/`, `prad2det/`, `prad2dec/`, `python/` or `analysis/`. `B` = `benchmark_hycal_timing_test_20261021.py`. `WA` = `prad2dec/src/WaveAnalyzer.cpp`. `PCc` = `analysis/tools/physics_calib.cpp`.

```
E=/Users/L00338853/work/PRad/data/evio
DB=/Users/L00338853/work/PRad/prad2evviewer/database
BIN=/Users/L00338853/work/PRad/prad2evviewer/build/bin        # rebuilt at HEAD in S1 (Open Question 7)
O=output/mfstd_20261102                                    # every output of this plan
TB=output/templates_spline_best50_20261026                     # R5nr banks (WT P1)
RAW=<raw-replay dir chosen in Open Question 5>
```

## Goal

Answer the question with a route that lets the standard calibration chain consume MF pulse energies, and validate it.

**The standard chain** (RF and the notes it builds on):

| # | Step | Tool, doc section | Consumes | Produces | Code |
|---|---|---|---|---|---|
| 1 | EVIO → raw ROOT; `--Ecalib` gives the compact calibration sample (sum trigger, peaks only) | `prad2ana_replay_rawdata`, RF Step 1, §1.2, §3.2, §6 | FADC samples | `events` tree: `hycal.samples` and/or `hycal.npeaks/peak_{height,time,integral}`, `hycal.gain_factor` | `prad2det/include/EventData_io.h:135-171` |
| 2 | Software pulse analysis | WaveAnalyzer, RF §1.2-1.5 | samples | per peak: integral = Σ(smoothed − ped) over an outward walk until 2 samples < max(0.1·height, ped_rms), in ADC·samples; time = CFD50 in raw ns | `WA:244` smooth, `:272` pedestal, `:358` peaks, `:429-455` integral; config `database/daq_config.json:85-125` |
| 3 | LMS gain time series | `prad2ana_replay_gainCorr`, RF §1.5, §5; GC §1-4 | first-peak integral of LMS events (`analysis/src/GainCorrCompute.cpp:406`) | `$DB/gain_factor/gain_correction/prad_<run>_gain_corr.root`; per batch c = g_ref/g_cur | `analysis/include/gain_factor.h:176` (`ModuleGain`), `:221-229` (`GetCorr` by event_num) |
| 4 | Per-module calibration factor | `prad2ana_physics_calib -i N` (PC; not linked from RF) | raw ROOT `peak_integral × c` of sum-trigger, single-cluster, 5×5-central e-p events | `Physics_calib/run<N>/calib_factor_iterN.json` = list of `{name, factor, base_energy, nl1, nl2}` | branches `PCc:535-543`; stored peaks preferred, WA only if no peak branches `PCc:559`; gain `PCc:562,568`; fit `PCc:380`; damped update ratio = 1 + 0.85(E_exp/peak − 1), clamp [0.5, 2] `PCc:394-398`; seed/iteration files `PCc:256-272` |
| 5 | Nonlinearity | `prad2ana_nonlinearity` on recon files | cluster E vs E_exp at 3 beam energies | nl1 (nl2 = 0) per module | applied per cluster: `prad2det/src/HyCalCluster.cpp:458-470` |
| 6 | Constants installed per run period | `database/runinfo/general.json` `calibration.file` (:44, :73-97, :119) | — | — | `prad2det/src/PipelineBuilder.cpp:299-300`; loader `prad2det/src/HyCalSystem.cpp:504-529` |
| 7 | Reconstruction | `prad2ana_replay_recon`, RF §1.3-1.4, §4 | samples (WA re-run) or, if no samples, stored peaks | E_module = factor · integral · c (`prad2det/include/HyCalSystem.h:153-157`); cluster E × bias × nonlinear | `analysis/src/Replay.cpp:993-1020` (branch choice), `:1168-1172`, `:1178-1237` |

**The calibrated quantity** is the per-module factor in MeV per (LMS-corrected WA integral in ADC·samples). An MF substitute must hand the chain a number with the same meaning in the same slot, `hycal.peak_integral`.

**Our MF chain today** (WT P1, S2-S5): E_u = A · I_int(m) · f(A) · κ with κ = 1 (`mf_pulse.py:1245-1258`). I_int is the WA walk-out integral of the noiseless unit template, phase-averaged (`mf_pulse.py:600-640`). f(A) maps MF to the WA integral across amplitude (`$TB/ampcorr.json`, 0.9940-0.9994 at 3.49 GeV). B already feeds E_u into the C++ `Module::energize` with a standard-format factor file (B:1618-1626). So the units are drop-in by construction. What is not standard is how the constants were derived: R5nr uses the old WA-derived files in `database/calibration/old/`, then per-crystal ridge gains (`percrystal_calib.py`), a position map, and an elastic-window pass (WT S5, S10, S11, S14).

**Goal.** Run the standard step 4 (`physics_calib`) on MF energies without changing any C++ or `database/` file. Then compare the MF-derived constants with WA-derived constants obtained by the same tool on the same events, and with the R5nr constants, and measure the resolution on the R5 reference splits.

**Non-goals.** No C++, binding or `database/` edits (R5 decision 5). No commit of constants to `runinfo/general.json`. No MF re-derivation of nonlinearity (step 5) or of LMS gains (step 3). No full `replay_recon` from MF (needs option (a), or the full raw schema; see Risks). No Ta or C points (Open Question 3).

**Expected change size.** One new tracked script, `mf_raw_writer.py`, about 300 lines. About 200 lines of helpers under `$O/tools/` (git-ignored, sha1 recorded as in R5/R6). No edits to tracked files. Machine time is in Open Question 1.

## Approach

Option (b): a Python writer reads a standard raw ROOT file that still has waveforms, runs the R5nr MF on its PbWO4 channels, and writes a peaks-only raw ROOT file with `hycal.peak_integral := E_u`. The unchanged `prad2ana_physics_calib` then calibrates on it. The same tool run on the original raw file (no peak branches, so it runs WA itself, `PCc:559`) gives the WA arm on identical events. B then measures resolution with each constant file through `--hc-calib-file`.

Options considered (effort in Open Question 1):
- **(a) C++ port.** A spline-MF fitter in `prad2dec` beside WaveAnalyzer, wired like the existing template deconvolution (`nnls_deconv`, `WA:667`; `PulseTemplateStore`, `prad2dec/src/PulseTemplateStore.cpp:115-122`). Every standard tool then gets MF. It is the production route, but the cost is large and it touches shared code; deferred until (b) shows a gain.
- **(b) Python peaks writer.** Chosen: no C++ edit, and the standard tool and its output format are used unchanged. It covers steps 4 and 6 only. Steps 3 and 5 are reused; step 7 from MF needs the full raw schema (Risks).
- **(c) Re-implement `physics_calib` in Python inside B.** Rejected: it duplicates the 5×5 and transition-module selection (PC §3-4), which would drift from C++. R5's `percrystal_calib.py` already is our own Python calibration.

## Decisions

1. **The observable and the substitution.** The calibration consumes `hycal.peak_integral` (ADC·samples) × c (LMS) and returns factor = MeV per unit (step 4). For PbWO4 channels the writer sets:
   - `peak_integral` = E_u = A·I_int·f(A)·κ, κ = 1, as `pulse_energy_units` (`mf_pulse.py:1245-1258`) is called at B:1846-1847;
   - `peak_height` = A;
   - `peak_time` = t0 + t50 (`mf_pulse.py:1236`), the CFD50-like time in raw ns. WA also stores raw times with no module offset (RF §1.2); `physics_calib` subtracts `time_offset` and applies `hc_time_cuts` itself (`PCc:572-590`). With `seed_time_window` 4.0 (`database/reconstruction_config.json:36`) it runs in multi-pulse mode: every in-window peak becomes a hit, and the clusterer applies the seed-time coincidence. MF times therefore enter the selection; `time_gate` (S1) and `paired_events` (S2) check this.
2. **Units: existing constants are drop-in to about 1 %.** I_int applies the WA integral rule to the noiseless template, and f(A) removes the amplitude dependence against WA. The median template/raw ratio was 0.978-0.994 for 100-1000 ADC (`notes/PbWO4_energy_resolution_20261023.md:193-215`), but that was measured on the earlier fixed-template bank; for the spline banks, gate G1-6 (I_int vs C++ integral) passed 5/5 (`notes/PbWO4_energy_resolution_20261026.md:159`). There is no per-module κ. Slice 1 re-measures on the R5nr bank the per-module ratio r_m = median(E_u / WA integral) over matched pulses with A ≥ 100 ADC (low-amplitude pulses, where WA is biased by 8-13 %, are excluded); if r_m is near 1 everywhere, the existing files remain usable with MF and a re-derivation is a refinement, not a fix.
3. **Route (b)**, per Approach. (a) is decided after Slice 3 on its result (Open Question 1).
4. **Input is the default raw replay (with waveforms), not `--Ecalib`.** `--Ecalib` has no samples, so MF cannot run, and it drops `event_num` (`EventData_io.h:137-140`), which `physics_calib` needs for the LMS lookup (`PCc:562`). The default raw file also gives the paired WA arm for free (`PCc:559`).
5. **Output schema.** One `events` tree with exactly the branches `physics_calib` enables (`PCc:535-543`), copied unchanged except the peak branches:
   - copied: `event_num/I`, `trigger_type/b`, `trigger_bits/i`, `hycal.nch/I`, `hycal.module_id[hycal.nch]/s`, `hycal.module_type[hycal.nch]/b`, `hycal.gain_factor[hycal.nch]/F`;
   - written: `hycal.npeaks[hycal.nch]/b`, `hycal.peak_{height,time,integral}[hycal.nch][8]/F`, `hycal.peak_quality[hycal.nch][8]/b`;
   - no `hycal.samples`. If samples were present, `replay_recon` would re-run WA and ignore the peaks (`Replay.cpp:1171`).
   - Writer backend: uproot 5.7.7 (installed). If its leaf layout is not read correctly by `prad2::SetRawReadBranches` (Slice 1, first test), the fallback is a ROOT macro under `$O/tools/` that fills the same leaf lists from a flat uproot tree (`root-config` is installed; PyROOT is not).
   - Only sum-trigger events are written (`physics_calib` skips the rest, `PCc:552`). `hycal.nch` keeps every channel, so the `nch > 100` cut (`PCc:554`) behaves as on the source file.
   - The output name keeps the source basename (`prad_<run>.evio.<split>_raw.root`), because `physics_calib` takes the run number from the first file name (PC §2.1).
6. **Fallbacks follow R5nr.**
   - Non-PbWO4 channels (lead glass, veto, LMS), odd-length waveforms, and clipped PbWO4 channels flagged `wa_fb` (`mf_pulse.py:1197-1214`) get the WA peaks from the bound analyzer (`python/bind_dec.cpp:388-394`, built as in `_common.py:237`), with no f(A). These are the same rules as B:1508, B:1551-1600 and B:1851-1860.
   - PbWO4 modules with no template of their own take the run-type average spline (`mf_pulse.py:740-748`).
   - Pulses with A ≤ 0 are dropped. At most 8 peaks are written (kmax 4 plus fallbacks; `MAX_PEAKS` 8).
   - `peak_quality` bit 7 marks an MF peak, so fallbacks can be counted afterwards.
7. **MF settings are the R5nr ones** (WT S3): `--mf-r-model 0.3`, `--mf-min-amp 5`, `--mf-kmax 4`, nsigma 5, shape fixed, ped-floor off, κ = 1. The bank is the R5nr bank of the point (`$TB/<runtype>/templates_<runtype>.npz`, with `$TB/ampcorr.json` folded in as B does), and its sha1 goes in the output header.
8. **LMS (step 3) is reused, not redone with MF.** c = g_ref/g_cur is a ratio of the same pulse type, so a scale common to both extractors cancels. LMS and alpha pulses are not e-p shower shapes, and the bank has no templates for them. `physics_calib` applies the `gain_corr` file if it exists for the run (none for 024xxx, so c = 1 at 3.49 and 0.73 GeV; 25320 has one), so both of its arms see the same c. B's LMS convention differs (stored (LMS2+LMS3)/2 normalised to the first batch, B:1194-1201; R5nr passes the pre-merge blob for 25320), so at 2.24 GeV only the `mfstd − wastd` contrast is clean (Decision 13).
9. **Seed and nonlinearity.** Iteration 1 is seeded (`-c`) with the R5nr file of the point: `database/calibration/old/calibration_factor_fadc_3p5.json`, `old/calibration_factor_2p2_fadc.json` or `old/calibration_factor_fadc_0p7.json`. These carry no nl1/nl2, so nl = 0 is inherited, as in R5nr (WT S7, "non-linear correction is inert"). Both arms use the same seed and the same N.
10. **Iterations.** Run `-i 1..N` until both hold, with a cap of N = 8:
    - median |ratio − 1| < 0.3 %;
    - ≤ 2 % of fitted modules at the 0.5 or 2.0 clamp (PC §9.3, here with numbers).
    Each arm stops on its own criterion. The report gives the iteration count of each arm.
11. **Calibration runs.** One `physics_calib` call per run, with all local splits of the run: 24503 (3.49 GeV, 10 splits), 25320 (2.24 GeV H2, 8) and 24590 (0.73 GeV, 10). These are the runs of the reference splits in `output/r5_20261028/configs/ref41.txt`, so constants and evaluation share the calibration period. One run per call avoids applying the first file's RunConfig (Ebeam, seed, output directory) to other runs (`PCc:244`). Modules with fewer than 40 entries are not fitted and keep the seed (`PCc:363`; PC §5.1 says 100).
12. **Template banks stay where they are for (b).** Production location, format and run selection are an (a) question (Open Question 1): spline knots and coefficients per module (`spline_t`, `spline_c`; much smaller than the 4081-point `shape` grid) as JSON under `database/waveform/`, selected like `analyzer.nnls_deconv.template_file` (`database/daq_config.json:118`), with a per-period override in `runinfo/general.json`.
13. **Validation**, all pre-declared:
    - **V1, constants.** Per module: C_MF/C_WA against 1/r_m. Gate: RMS of (C_MF/C_WA)·r_m − 1 < 0.5 % over modules with ≥ 100 entries. Also report C_MF / C_seed, C_MF / C_runinfo and C_WA / C_runinfo (C_runinfo = the file `runinfo/general.json` installs for the run, e.g. `calibration_factor_3p5_June7_nonLinear.json`), and C_MF / (C_seed·g_i) with g_i the R5 per-crystal gains (`output/r5_20261028/gains/<run>/gains_<run>_all.tsv`): the median (a global offset; g0 ≈ 0.93 at 3.49 GeV, WT P2) and the spread after removing it.
    - **V2, resolution.** B on the reference splits of the three points (from `output/r5_20261028/configs/ref41.txt`), with these arms:
      - `mfstd`: R5nr command, `--hc-calib-file <MF iterN>`, no `--gain-file`, no `--pos-corr`;
      - `wastd`: the same with `--pulse-method rawsum` and `<WA iterN>`;
      - baselines: R6 `nopos` (`output/r6_20261030/nopos`) and R5nr.
      Statistic as R6 decision 8: paired per-split Δ, mean ± SE, CB headline, legacy and gated, with `compare_scan_dirs.py`. Contrasts: `mfstd − wastd` (MF under the standard calibration) and `mfstd − nopos` (standard calibration vs per-crystal ridge gains, both MF). At 2.24 GeV the headline is `mfstd − wastd` only; `mfstd − nopos` is reported with the LMS caveat of Decision 8.
    - No promotion is decided here. The results feed Open Question 1.

## Outline

**`mf_raw_writer.py`** (new, tracked):
```
mf_raw_writer.py IN_RAW.root --mf-templates NPZ --out OUT_DIR
    [--ampcorr $TB/ampcorr.json] [--daq-config $DB/daq_config.json]
    [--mf-r-model 0.3 --mf-min-amp 5 --mf-kmax 4]          # κ = 1 fixed
    [--max-events N] [--diag]              # --diag also writes <stem>_ratio.tsv
-> OUT_DIR/<IN basename>  (Decision 5 schema)
   OUT_DIR/<stem>_writer.json   {in_sha1, bank_sha1, argv, n_events, n_mf, n_wa_fb, n_odd, n_lg, n_dropped}
   <stem>_ratio.tsv  module  n_match  r_med  r_iqr  dt_med_ns   (--diag; pulses matched |dt| < 8 ns, A ≥ 100 ADC)
```
- Read with uproot in chunks (`hycal.samples[nch][MAX_SAMPLES]`, `hycal.nsamples`). Build the bank and reconstructor as B does (B:2624-2636; `mf_pulse.build_bank` :773-823; `SampledTemplateBank.from_npz` :708-770; `MFReconstructor.run` :1165-1242). Map module_id to bank index the way MFBlock does (B:1466 onward).
- Gotcha: `hycal.module_type` decides PbWO4 vs other. Use the same `is_pwo4` test as C++ (`Module.is_pwo4`, bound in `python/bind_det.cpp`), not a module-name prefix.
- Gotcha: the WA binding returns times with `time_offset = 0` only when called that way. Keep raw times on both arms.
- Invariant: for a channel with no accepted MF pulse and no fallback, `npeaks = 0` (no hit), as in R5nr.

**`$O/tools/`** (git-ignored; sha1 in `$O/git_state.txt`):
- `run_points.sh`: per point, `replay_rawdata` (default mode) → `mf_raw_writer.py --diag` on every split → `physics_calib` loops for the `wa` and `mf` arms. Records the sha1 of every `$BIN` binary used in `$O/git_state.txt`.
- `compare_constants.py`: reads both arms' `calib_factor_iterN.json` and `calib_result_iterN.json`, the seed, the R5 gains and the `--diag` ratio TSVs. Writes `$O/tables/constants_<point>.tsv` and the V1 summary `$O/tables/v1.md`.
- `v2_jobs.py`: rewrites the reference-split lines of `output/r5_20261028/jobs/final_nr.txt` into the `mfstd` and `wastd` arms (Decision 13), in the R5/R6 job format.

**Output tree**: `$O/raw_mf/<run>/`, `$O/calib/{wa,mf}/Physics_calib/run<run>/`, `$O/scan/{mfstd,wastd}/`, `$O/tables/`, `$O/compare/`, `$O/git_state.txt`.

## Slices

### S1. Writer on one split, schema and unit gates (≈ 0.5 day)
Delivers the writer and the per-module MF/WA ratio on 24503.00020.
- Units: `mf_raw_writer.py`; `$RAW/prad_024503.evio.00020_raw.root` from `$BIN/prad2ana_replay_rawdata $E/prad_024503.evio.00020 -o $RAW`.
- Tests:
  - `schema_roundtrip`: a `--max-events 2000` output → `physics_calib -i 1` prints "waveform branches: no, peak branches: yes" and selects > 0 events.
  - `energy_identity`: for 100 random PbWO4 pulses, `peak_integral` equals `pulse_energy_units(bank, m, t0, A, 1.0)` within 1e-6 relative.
  - `fallback_counts`: in the writer JSON, n_lg > 0, `n_wa_fb` equals `stats["n_sat_fallback"]` of `MFReconstructor` (`mf_pulse.py:1238-1241`), and the channel counts n_mf + n_wa_fb + n_lg + n_odd equal the number of channels with peaks.
  - `unit_gate`: in `_ratio.tsv`, |r_med − 1| < 0.02 for ≥ 95 % of modules with n_match ≥ 50.
  - `time_gate`: |dt_med| < 1 ns for the same modules.
- Checkpoint: `python3 mf_raw_writer.py $RAW/prad_024503.evio.00020_raw.root --mf-templates $TB/3p49GeV_H2/templates_3p49GeV_H2.npz --out $O/raw_mf/024503 --diag && $BIN/prad2ana_physics_calib $O/raw_mf/024503/prad_024503.evio.00020_raw.root -o $O/calib/mf -i 1 -c $DB/calibration/old/calibration_factor_fadc_3p5.json -j 4`.

### S2. Both calibration arms at 3.49 GeV (≈ 0.5 day, mostly machine)
Delivers converged C_WA and C_MF for run 24503 and the V1 table.
- Units: `run_points.sh`; `compare_constants.py`.
- Tests:
  - `convergence`: each arm meets Decision 10 by N ≤ 8.
  - `v1_gate`: the RMS of (C_MF/C_WA)·r_m − 1 is < 0.5 %.
  - `paired_events`: both arms report the same number of selected events within 1 %. They differ only through the energies entering the cluster cuts; a larger gap is reported, not fixed.
- Checkpoint: `bash $O/tools/run_points.sh 24503 && python3 $O/tools/compare_constants.py --point 3p49GeV_H2 --out $O/tables`.

### S3. Resolution at 3.49 GeV (≈ 2 h machine)
Delivers V2 for the 3.49 GeV reference splits.
- Units: `v2_jobs.py`, then the jobs.
- Tests:
  - `jobs_diff`: each `mfstd` command differs from its R5nr command only in `--hc-calib-file`, the dropped `--gain-file` and `--pos-corr`, and the stem; `wastd` also differs in `--pulse-method rawsum` (token diff saved).
  - `all_ran`: every job log ends `rc=0`.
- Checkpoint: `python3 compare_scan_dirs.py --scans mfstd=$O/scan/mfstd wastd=$O/scan/wastd nopos=output/r6_20261030/nopos r5nr=output/r5_20261028/scan/r5nr_splits --gauss-window 2sigma --out $O/compare/3p49/cmp` (`--out` is a prefix).
- Stop point: report to the user before S4 if `mfstd` is worse than `wastd` by more than 1 SE.

### S4. 2.24 GeV H2 and 0.73 GeV (≈ 0.5 day, mostly machine)
Delivers S2 and S3 for runs 25320 and 24590, with the same tests and checkpoints per point.

### S5. Report (≈ 2 h)
Delivers `notes/mfstd_20261102.md`: the step table of Goal, the V1 and V2 tables, the discrepancy list, and the recommendation for Open Question 1.
- Test `report_complete`: the note holds V1 and V2 tables for 3 points, the iteration count per arm, and the discrepancy list.
- Checkpoint: `grep -c '^| ' notes/mfstd_20261102.md && ls $O/tables/constants_{3p49GeV_H2,2p24GeV_H2,0p73GeV_H2}.tsv`.

## Verification

- `$O/tables/v1.md` covers three points, each with both arms converged.
- `$O/compare/{3p49,2p24,0p73}/` hold paired Δ tables for `mfstd − wastd` and `mfstd − nopos`.
- `git status --short` shows only `mf_raw_writer.py`, the note and this plan as new tracked files, and nothing under `database/`, `prad2det/`, `prad2dec/`, `python/` or `analysis/src|tools`.

## Risks

- **Selection differs from R5.** `physics_calib` uses a 5×5 sum and per-module damped fits (PC §3-5), with the database clustering config (`min_center_energy` 10 at `database/reconstruction_config.json:16` vs 30 at `output/r5_20261028/configs/reconstruction_config_r5.json:14`). R5 uses cluster sums and a joint ridge solve. A module-to-module spread in V1's C_MF/(C_seed·g_i) is therefore expected; it is reported, not tuned.
- **In-sample.** The banks contain the evaluated splits (WT P1 caveat), and run 24503's splits include reference splits. R5 is also in-sample, so the comparison is like for like; neither number is a held-out estimate.
- **LMS convention at 2.24 GeV.** For 25320, `physics_calib` uses the database `gain_corr` (ref-run ratio), while B (R5nr) uses its own pre-merge LMS blob (WT S4). Constants derived under one convention and applied under the other carry a per-module scale error, so `mfstd − nopos` at 2.24 GeV is caveated (Decisions 8, 13; Open Question 6).
- **Full recon from MF is not covered.** `replay_recon` on MF output would need GEM, VTP and TDC branches too. The writer copies only the `physics_calib` set. Copying the full raw schema through uproot is untested. This is the main argument for (a) if MF is promoted.
- **Disk.** Default raw replay keeps waveforms and GEM strips for about 31 splits (Open Question 5).

## Scope Notes

- In: steps 4 and 6 of the standard chain on MF energies; V1, V2; the report.
- Out: C++, binding and `database/` edits (R5 decision 5, `plans/plan_r5_20261028.md`); LMS and nonlinearity re-derivation; Ta, C and X17; full recon from MF; option (a) implementation.
- Do not reopen: κ = 1 and the R5nr MF settings (WT S3); refine-off for B arms (`plans/plan_r5_20261028.md:48-64`).

## Assumptions

- `hycal_time_cuts` from `runinfo/general.json` accept MF times as they do WA times (checked by `time_gate` in S1 and by `paired_events` in S2).
- Each calibration run's local splits give ≥ 40 entries for the inner modules (checked by the count of fitted modules in S2/S4).

## Open Questions

1. **Route and effort.** Recommended: (b) now; decide (a) after S4.

   | Option | What | Effort |
   |---|---|---|
   | (a) C++ spline-MF fitter in `prad2dec`, bank converter, config keys, bindings, parity tests | all standard tools use MF; production route | ≈ 1500-2000 lines; 3-5 weeks including parity with Python at 1e-4 and maintainer review |
   | (b) this plan | `physics_calib` on MF; constants in the standard format | ≈ 500 lines; 3-4 days plus ≈ 8 h machine |
   | (c) Python `physics_calib` | duplicate of the C++ selection | ≈ 350 lines; 1-1.5 weeks including parity with C++ |

2. **What the constants are for.** Recommended: analysis-only. They are a candidate replacement for R5's seed file × per-crystal gains in B, not a `database/` submission. A submission would also need nonlinearity and the C++ route.
3. **Points.** Recommended: 3.49, 2.24 H2 and 0.73 GeV only. C shares the 2.24 calibration period. Ta uses the X17 files and its own bank.
4. **Seed.** Recommended: the R5nr `old/` files (Decision 9), which give nl = 0 and parity with R5nr. The alternative is the run-period default (`calibration_factor_3p5_June7_nonLinear.json`, with nl1/nl2 on 1006 W modules), which is closer to the collaboration's current state.
5. **Raw replay location.** Recommended: outside the repo, next to `$E` (for example `/Users/L00338853/work/PRad/data/replay_raw_mfstd`), deleted after S4.
6. **LMS at 2.24 GeV (run 25320).** Recommended: keep the standard behaviour in `physics_calib` (database `gain_corr`) and the R5nr LMS in B, headline only `mfstd − wastd` there (Decision 13). Alternative: LMS off on both sides (a shadow `PRAD2_DATABASE_DIR` without the 25320 `gain_corr`, as R6 did for `energy_bias/`, plus `--lms-corr off`), which makes `mfstd − nopos` comparable but departs from both the standard chain and R5nr.
7. **Rebuild the tools before S1.** The built `$BIN/prad2ana_physics_calib` and `prad2ana_replay_rawdata` date from Sep 14 and predate the source changes the plan relies on: the CB fit, `FillPeaksFromWaveforms` (`PCc:559`) and the LMS loader (commits 2711d35, 165a00c, d606375, b7f9dbc, 6e638c6). Recommended: yes, `cmake --build build --target physics_calib replay_rawdata -j4` at HEAD (build directory only, no source edit), with the binary sha1s recorded.

## User Decisions (2026-11-02)

- OQ1: route (b) now; decide (a) C++ port after the results.
- OQ2: analysis-only (candidate replacement for seed x per-crystal gains in B); no database submission.
- OQ3 (default, no objection): points 3.49, 2.24 H2 (run 25320), 0.73 GeV only; no C, Ta, X17.
- OQ4: iteration-1 seed = R5nr `database/calibration/old/` files (nl = 0).
- OQ5 (default): raw replays at `/Users/L00338853/work/PRad/data/replay_raw_mfstd`; not deleted at the end (size reported; deletion left to the user).
- OQ6: keep both LMS conventions at 2.24 GeV; the 2.24 headline is `mfstd - wastd` only.
- OQ7 (default): tools rebuilt at HEAD, no source edits, binary sha1s recorded. Executed in a separate build dir `/Users/L00338853/work/PRad/data/build_mfstd` (BUILD_ANALYSIS=ON, python off; the shared `build/` has BUILD_ANALYSIS=OFF and is left untouched). `-include stdlib.h -include unistd.h` was added to CMAKE_CXX_FLAGS because `GainCorrCompute.cpp:143` uses `mkdtemp` without a header on this macOS toolchain (finding; no source edit).
- Output root is `output/mfstd_20261102` (`$O`), replacing `output/mfstd_20261102`.
- Finding to report: `--Ecalib` drops `event_num` (`EventData_io.h:137-140`) but `physics_calib` uses it for the LMS lookup (`physics_calib.cpp:562`); our runs use the default raw replay, which keeps `event_num`.

## Log

- 2026-11-02 executor: S1-S5 done, all three points; V1 gate failed (median -2.2 %, see notes/mfstd_20261102.md); V2 mfstd < wastd by 0.23-0.50 pp, nopos/r5nr better than mfstd; deferred: none

## Todo

- [x] S1: rebuild `physics_calib` and `replay_rawdata` (Open Question 7); raw replay of 24503.00020 into `$RAW`
- [x] S1: write `mf_raw_writer.py` (Decisions 1, 5-7) and run the five S1 tests
- [x] S1: if `schema_roundtrip` fails, write the macro fallback and rerun
- [x] S2: `run_points.sh`; raw replay of the 10 splits of 24503
- [x] S2: both arms to convergence; `compare_constants.py`; V1 gate
- [x] S3: `v2_jobs.py`; run `mfstd` and `wastd` on the 3.49 reference splits; compare; stop-point check
- [x] S4: repeat S2-S3 for 25320 and 24590
- [x] S5: write `notes/mfstd_20261102.md`; record tool sha1s in `$O/git_state.txt`
