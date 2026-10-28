# Plan: best-50 bank pipeline — pile-up efficiency, σ_E/E and σ_t per module and per cluster (2026-10-26)

**Status:** draft for user approval (planner, 2026-10-26 17:05; critic round 1 applied). **Branch:** `waveform_jz`.
**User request (verbatim):** "With the new extracted templates, go through pile-up deconvolution studies and I want the table of the pile-up rejection efficiency, then run everything again to get sigma_E/E and sigma_t per module and per cluster. Put all of these in a new folder. Make notes for today."
**Inputs read:** `notes/PbWO4_energy_resolution_2026102{2,3,5}.md` (§7, §9.4, §9, §9.7, §11), `notes/PbWO4_energy_resolution_summary.md` §3.8, §4–§6, `plans/plan_2026102{4,5}.md`, `plans/plan_rf_calib_20261026.md`, `output/scan_20261025/configs/{cfgs_hybrid.json,run_hybrid.sh}`, `output/scan_20261025/manifest.json`, the running job's `contract.md` and `run_s3b_s4.sh` (session scratchpad `best50/`), and the scripts named in §2.

Paths are relative to `/Users/L00338853/work/PRad/prad2evviewer/analysis/pyscripts/` unless they start with `database/` or `prad2det/` (repo root). Shell variables used throughout: `E=/Users/L00338853/work/PRad/data/evio`, `DB=/Users/L00338853/work/PRad/prad2evviewer/database/calibration`, `T0=output/templates_spline` (production bank), `TN=output/templates_spline_best50_20261026` (new bank), `O=output/best50_20261026` (every output of this plan), `S=<session scratchpad>/best50_pipeline` (scratch only), `G="python3 gate_scan_20261022.py"`. Every command runs after `cd /Users/L00338853/work/PRad/prad2evviewer/analysis/pyscripts && source setup_env.sh && export MPLBACKEND=Agg`.

---

## 1. Goal and approach

**Goal.** With the escalated best-50 spline bank (`$TN/<slug>/templates_<slug>.npz`, five slugs), (i) rerun the pile-up detection-efficiency study and report its table against the production bank; (ii) rerun the R3 energy-resolution pipeline (per-crystal gains, the 118-split scan, module-by-module resolution) with the new bank, position maps fixed to the R3 maps, and report σ_E/E per run type paired against R3; (iii) rerun the module timing chain (dump, plots, RF calibration) with the new bank and report σ_t per module against the 10-25/10-26 results; (iv) define and measure σ_t per cluster against the RF; (v) consolidate the comparison tables and hand them to a new daily note. Everything lands under `output/best50_20261026/`.

**Approach.** Every stage is the recorded recipe of the previous run with only the template, gain, position-map and output paths changed (the harness `bench_args` placeholders `{run}`/`{runtype}` do the substitution; the timing scripts get the three missing path flags). The template effect is isolated: position maps are the R3 maps (hybrid for 3.49/2.24 H2/C, bilinear for 0.73/Ta, through a per-run-type symlink tree so one harness tag covers all 118 splits); the pile-up headline arm reuses the production synthetic pile-up so only the reconstruction bank changes. The running best-50 job already produced the gain dumps and fits of five runs with the final bank and the percrystal_spl recipe; they are reused, and only the nine remaining runs are dumped and fitted. Per-cluster times do not exist in any output today; the benchmark gets one default-off flag that appends per-member times to `result_clusters.tsv`, and a new offline script folds seed and energy-weighted cluster times against the RF with the per-run offsets of stage S5. Stages run serially at 12 jobs (the setting of every previous full run), after the running job releases the machine.

**Rejected alternatives.** A single new-bank synth arm (confounds the injected shape with the bank; kept as arm C with D as control); re-deriving all 14 gains (the five finished fits use the identical recipe); two harness tags for the map mix (one posmap symlink tree is simpler); refitting the maps (option, §9); an in-benchmark offset-corrected cluster time (couples S3 to S5); rebuilding clusters from the timing dump (it has no clusters); overlapping the dump with the scan (option, §7).

## 2. Facts the plan relies on (verified 2026-10-26, 15:50–16:30)

| # | Fact | Where |
|---|---|---|
| F1 | New bank at its final path: `$TN/<slug>/templates_<slug>.npz` for 0p73GeV_H2, 2p24GeV_H2, 2p24GeV_C, 2p24GeV_Ta_X17, 3p49GeV_H2 (npz 15:28, dirs moved 15:36; `esc_stage/` no longer exists). Every `summary.tsv` has the `osc_escalated` column; escalated modules 172 / 83 / 127 / 152 / 176, fallback-0 modules 607 / 1059 / 753 / 840 / 869 (0.73 / 2.24 H2 / C / Ta / 3.49). The npz is a superset of the production keys (adds `converged, holdout_keys, iter_hist, lam_module, n_iterations`); `mf_pulse.build_bank` loads by key, and the running job's benchmark dumps run on it (rc = 0). Build settings: `--pool-order energy --n-max 50 --fit-points direct --shape-sigma 0.012 --lambda-select heldout --lambda-chi2-sigma full --lambda-min 1e2 --lambda-fallback 1e3 --converge-tol 0.01 --converge-df 5e-3 --max-iterations 10 --osc-escalate` | `ls`, `awk` on summary.tsv, `zipfile` listing, `best50/run_esc.sh` |
| F2 | Pile-up study (10/23–10/24): no table is literally named "rejection efficiency"; the reported table is the **detection-efficiency grid** of `benchmark_pileup_mf.py` on `pileup_generator.py` pseudo-data, current version `output/pileup_eff_20261024_fix/<slug>/eff_<slug>.{md,json,png}` (Oct 24 08:43, production bank sha1 `d6c49edeb6bc` for 2p24GeV_H2 = current `$T0`). Grid: Δt ∈ {8, 12, 16, 20, 24, 32, 48, 64, 96, 128} ns (generator Δt = injected peak − base 50 % CFD; onset separation ≈ Δt − 6.2 ns) × amplitude ratio ∈ {0.1, 0.2, 0.3, 0.5, 1, 2}, 50 events per cell, per run type (= the energy axis). Quantities per cell: `eff_match` (% events with an mf pulse within ±4 ns of the injected onset **and** one within ±4 ns of the base onset), `fake_rate` (mean mf pulses per event matched to neither onset), `n2_mf` (% events with ≥ 2 mf pulses), `n2_wa` (% events with ≥ 2 WaveAnalyzer peaks; bank-independent), `*_clean` rows (base onset ≥ 60 ns). Blocks: onset-convention check; r_model 0.3; r_model 0.3 + window rule 40:0.16; r_model 0.1; r_model 0.1 + window rule (json keys `dt_ns, ratios, results{headline, second}, results_window{headline, second}, window_rule_check, settings`; headline = r_model 0.3, second = 0.1). `MIN_SEP_NS = 12.0` is a constant in `mf_pulse.py:62`, printed in the md header. The synth depends on the bank (the injected pulse is the per-module template shape; `template_sha1` in the json); the base-event selection does not. Generator 6–7 s and benchmark 7.2–7.4 s per slug. Headline of note 20261023 §9.4: r_model 0.3, ratio ≥ 0.5 found 100 % from Δt = 32 ns (Ta 96 %), ratio 2 from 20 ns, ratio ≤ 0.3 only at 128 ns; fakes 0 except 0.1–0.8 per event at Δt 16–24 ns, ratio 1 | `benchmark_pileup_mf.py:255-289, 333, 534-684`; `pileup_generator.py:1732-1810`; `eff_2p24GeV_H2.md:6-40`; summary §3.8, §6 (lines 654-667); note 20261023 §9.4 |
| F3 | Pile-up **rates** (`plot_pileup_examples.py --jobs 12 --splits all …`, 119 splits, 59 min) are a separate product (`output/pileup_20261024_fix/pileup_rate_summary.tsv`); not the requested table | summary §6 line 661; plan_20261024 log |
| F4 | R3 = harness tag `pos_hybrid` for 24386, 24503, 24512, 25228, 25257, 25260, 25299, 25308, 25320 and tag `pos_all` for 24590, 24603, 25800, 26138, 26230, both with `bench_args "--mf-templates $T0/{runtype}/templates_{runtype}.npz --mf-scale 1.0 --gain-file output/percrystal_spl/{run}/gains_{run}_all.tsv --pos-corr <map> --write-cluster-tsv"`, `time_calib off, W 16, LO 100, HI 200, pulse mf`; `<map>` = `output/poscorr_20261025/hybrid/{runtype}/poscorr_hybrid_{runtype}_all.json` (hybrid) or `output/poscorr_20261025/{runtype}/poscorr_{runtype}_all.json` (bilinear; the Ta file is a symlink to the 2p24GeV_H2 bilinear map). The R3 mix exists as the symlink tree `output/resolution_fit_20261026/r3_splits/scan_<run>` (118 splits, driver layout). `build_cmd` adds `--calibrate --exclude-rowcol 15 20 15 20 --pre-window 100.0 200.0 --inner-radius 200 --min-per-module 30 --lms-corr auto|off --hc-calib-file … --window 16.0 --time-calib off --pulse-method mf --mf-scale 1.01 --mf-r-model 0.3` (a later `--mf-scale 1.0` in bench_args wins). Stage R = 93 (run, split) pairs (no 24578, no X17, no 24386.00060), RX = 25 X17 pairs; `--runs R1,R2` (6-digit or plain), `--tags`, `--dry-run`, `--force`; layout `<out>/<tag>/scan_<run>/split_<n>/{job.cmd, benchmark.log, result.tsv, result_clusters.tsv, …}` | note 20261025 §4.9; `cfgs_hybrid.json`; `run_hybrid.sh`; `gate_scan_20261022.py:144-150, 349-377, 395, 2150-2179` |
| F5 | R3 scan cost: stage R 93 jobs at `--parallel 12` 22:30:24→23:23:22 (53 min), RX 25 jobs 7.9 min; per-split `elapsed_s` 38–564 s (median 382, mean 348; 11.4 core-hours); slowest runs 24590/24603 (617 k / 364 k dump events). Peak memory not logged; the running job's benchmark processes use 2.5–2.8 GB RSS each. Machine: 18 cores, 48 GB, 78 GB free on the volume | `output/scan_20261025/manifest.json`; result.tsv headers; `ps`, `sysctl`, `df` |
| F6 | percrystal_spl recipe: harness dump pass (`bench_args "--mf-templates $T0/{runtype}/templates_{runtype}.npz --mf-scale 1.0 --cal-method percrystal --cluster-dump auto"`, stages R + RX, 118 jobs ≈ 90 min) + the 24386.00060 extra dump (one direct benchmark call, `--lms-corr off`, `--cluster-dump <file>`), then `percrystal_calib.py fit --dumps "<tag>/scan_<run>/split_*/result_dump.npz" [extra] --splits all --fixed-windows --out output/percrystal_spl/<run>/` × 14 → `gains_<run>_all.{tsv,json}` (json records `argv, splits, n_events, n_free, n_prior, g0, converged`). **Running job:** `output/percrystal_best50_20261026/dumps/scan_<run>/split_<n>/` for 024503 ×10, 024590 ×10, 025308 ×8, 025320 ×8, 025800 ×8 (all splits of those runs), job.cmd = the harness dump command with `--mf-templates $TN/…` (started 15:36, after the final bank), and `output/percrystal_best50_20261026/<run>/gains_<run>_all.{tsv,json}` with argv `fit --dumps … --splits all --fixed-windows`, all converged; n_free 283 / 267 / 284 / 284 / 138 (spl: 283 / 267 / 284 / 284 / 139); g0 e.g. 025308 0.9199 vs spl 0.9123. The 4-split R3 benchmark of the running job writes `output/scan_20261026_best50/scan_<run>/split_*/` (20 jobs, 10 done at 16:28) and `control/` (5 R3 reruns with the production bank, one checked diff-identical to scan_20261025) | summary §6 lines 643-646; plan_20261025 §13 lines 226-249; json `argv`; `fit_*.log` (`real 3.95` s) |
| F7 | Aggregation: `compare_scan_dirs.py --scans LABEL=DIR … --out PREFIX [--gauss-window 2sigma] [--exclude-runs 24578]` pairs the common (run, split) keys and writes PREFIX.md/.tsv (`type slug path fit n_splits ref cmp ref_mean cmp_mean delta delta_sd delta_se nmod_ref nmod_cmp nev_ref nev_cmp eps_yield eps_peak`; fit ∈ {G, CB}, path ∈ {legacy, gated}, pre and post arms) and PREFIX_splits.tsv; `output/scan_20261025/pooled_delta.py LABEL SPLITS_TSV` prints the pooled paired-Δ markdown table (mean ± SE over splits, pp). Preferred fit: CB for 3.49 / C / Ta, Gaussian for 2.24 H2 / 0.73 | `compare_scan_dirs.py:447-458`; `output/scan_20261025/pos_hybrid/compare/`; note 20261025 §4.9 |
| F8 | Module-by-module resolution: `plot_module_resolution.py --scan-dir DIR --out DIR --group-runs` (defaults `--path legacy --min-n 1000 --fit auto --gauss-window 2sigma --inner-radius 200 --exclude-rowcol 15 20 15 20`), 5–7 s per run; per module σ_M/μ_M (`post_g_s_over_mu`, `post_cb_s_over_mu`, pre values, `post_resid`), low_stat for 0 < N ≤ 1000, `summary/types_summary.tsv` (`post_median post_p16 post_p84 …`), per-type dirs `3p49_H2 2p24_H2 0p73_H2 2p24_C Ta_X17`, `<run>/<run>_modules.tsv`. `output/module_resolution_20261025_hybrid/` was run on `pos_hybrid` for **all** 14 runs (not the R3 mix) | `plot_module_resolution.py:521-541`; `run_hybrid.sh` |
| F9 | Module timing dump: `python3 module_timing_dump.py --all --parallel 12 --out output/module_timing_20261025/dump`; **no template flag**: `_T = "output/templates_spline/{0}/templates_{0}.npz"` and `RUN_CFG` (lines 67-80; types 3p49 → 24386, 24503, 24512; 2p24 → 25228, 25257, 25260, 25299, 25320; 0p73 → 24590, 24603; x17 → 25800, 26138, 26230 with the `2p24GeV_Ta_X17` bank; no 25308); mf settings hard-coded `nsigma 5, min_amp 5, kmax 4, r_model 0.3`. 111 splits, 13 runs; per-split `elapsed` 111–510 s, sum 43 610 s → 61 min at 12 parallel; dump 4.9 GB (8–60 MB per npz). Template-dependent keys: `m_t, m_E, m_npk, m_inwin, m_fb, m_ev, m_idx, templates`; independent: `ev_num, trig, v_*` | `module_timing_dump.py:60-82, 304-314`; `dump_run.log`; note 20261025 §11(h) line 820 |
| F10 | Plots (all take `--dump-dir`, `--out-dir`): `plot_module_timing.py --ref veto` (83 s; `module_timing.pdf`, `module_timing_veto.tsv` with per-type `N_<t> mean_<t> sigma_<t> fit_ok_<t>`, `veto_summary.tsv`; reference = V1–V4 WaveAnalyzer CFD time, H2 runs; X17 uses `t_seed`), `plot_module_timing_vs_integral.py` (107 s), `plot_module_timing_rf.py [--rf-dir]` (312 s; `module_timing_rf.tsv` with `module id type N N_excl n_runs offset_raw offset_tcal sigma sigma_err fit_ok …`, `rf_summary.tsv`, `module_timing_rf{,_tcal}.pdf`). RF dumps `output/module_timing_20261025/rf/` (111 npz, 582 MB, keys `ev_num trig ts rf_off rf_t rf_ch`; `rf_time_dump.py` reads no templates; 38 s) are reused | `plot_module_timing.py:1113-1134`; `plot_module_timing_vs_integral.py:290-304`; `plot_module_timing_rf.py:483-494`; `run_veto.log`, `run_rf.log`, `rf_dump_run.log` |
| F11 | RF calibration: `MPLBACKEND=Agg python3 rf_time_calib.py --jobs 12` (80 s); flags only `--out-dir --jobs --runs --no-pdf --bimodal-ratio`; `DUMP_DIR`, `RF_DIR`, `XCHECK` are module constants (lines 73-77) pointing at `module_timing_20261025`; gate (vi) reproduces `rf_summary.tsv:run_offset_tcal` of the **same** dump, so the new `rf_summary.tsv` (from `plot_module_timing_rf.py` on the new dump) must be the xcheck. Outputs `<run>.json` ×13 (5-digit run, e.g. `24386.json`; the RF npz use the 6-digit run, `024386_00010.npz`; `database/hycal_rf_offsets/24840.json` schema: `default, modules[{name, offset_ns, …}]`, plus `usable, quality, time_definition`), `rf_calib_modules.tsv` (`run type module id crate slot n_e100 n_allE tier offset_ns offset_err_ns sigma_ns tau_ns rho_ns rho_run_ns … dhalf_ns dhalf_err_ns flags`), `rf_calib_runs.tsv` (`run type n_splits P_ns usable quality … median_sigma_e100 frac_sigma_lt_0p3`), `rf_calib_validation.txt`, PDFs; exit 1 on a failed gate. σ_t per module = `RF.fit_unwrapped` (80-bin circular histogram of fold(t_m − t_RF), unwrapped around the smoothed argmax, ±2σ Gaussian core) on E ≥ 100 MeV, N ≥ 300, else all energies. All 13 runs have RF (match rate 1); 24386 is bimodal (unusable), 24503 drifts (0.13 ns) | `rf_time_calib.py:70-80, 452-456`; note 20261025 §9.7; plan_rf_calib §3–§4 |
| F12 | `result_clusters.tsv` (`--write-cluster-tsv`, requires `--calibrate`): header `event_num config_name seed_module_id seed_x_mm seed_y_mm seed_row seed_col seed_energy_MeV cluster_E_raw_MeV cluster_E_cal_MeV cluster_x_mm cluster_y_mm`; one row per single-cluster event per path (`legacy`, `gated`); `event_num` is the running index; the TI event number is `ex_[5]` of the event tuple (`evnum`, kept by `member_rec`, benchmark lines 2317-2329, which returns `(mi, me, r.hit.x, r.hit.y, evnum)` appended after `f_eff`). Hit times fed to the clusterer are the mf `p_time` (minus the `--time-calib` offset, off in R3) in ns (`_event_hits`, lines 1289-1314); `r.cluster.center.time` (seed time) and `r.cluster.hits[i].time` are exposed by the binding (`python/bind_det.cpp:734, 750`; `HyCalCluster.cpp:462`). 23 MB per split today. No time columns anywhere; `--cluster-dump` npz has membership (`ptr, midx, E`) but no times | benchmark lines 1289-1314, 2317-2329, 3055-3082; `bind_det.cpp` |
| F13 | Benchmark edits are unsafe while the running job's `xargs … runjob.sh` launches `benchmark_hycal_timing_test_20261021.py` processes (it runs the working-tree script); `module_timing_dump.py`, `rf_time_calib.py` and new files are not used by it. The running job also edits `notes/PbWO4_energy_resolution_20261025.md`, `notes/PbWO4_energy_resolution_summary.md`, `plans/plan_20261025.md` (its contract S5) | `best50/contract.md`, `run_s3b_s4.sh`, `runjob.sh` |

## 3. Definitions

| Quantity | Definition used in this plan |
|---|---|
| Pile-up efficiency table | Per run type and per block (r_model 0.3; 0.3 + window rule 40:0.16; 0.1; 0.1 + window rule): the 10 × 6 grids of `eff_match` (%), `fake_rate` (per event) and `n2_mf` (%) of `benchmark_pileup_mf.py`, exactly as in `eff_<slug>.md` (F2). Arms: **A** = production synth × production bank (exists, `output/pileup_eff_20261024_fix/`); **B** = production synth × new bank (headline: isolates the reconstruction bank); **C** = new synth (generator with the new bank, same base splits, seed 12345) × new bank (the study's own pattern); **D** = new synth × production bank (control for the injected-shape effect). The reported table is B with Δ(B − A) per cell in points, plus the pooled Δ over the 60 cells and over the cells with A < 100 %, and C/D summarised as pooled Δ only |
| σ_E/E per run type | `compare_scan_dirs.py` output: per type × path (legacy, gated) × fit (G = iterative ±2σ Gaussian, CB) the post means of ref (R3) and cmp (new) and the paired per-split Δ (cmp − ref) mean ± SE in pp, on the 118 common splits; preferred fit per type as in F7 |
| σ_E/E per module | `plot_module_resolution.py --group-runs` σ_M/μ_M per module (legacy path, N > 1000 pooled over the runs of the type), new bank vs the R3 reference tree; comparison = per-module Δ with median / p16 / p84 per type and the number of modules fitted in both |
| σ_t per module | S5: `module_timing_rf.tsv:sigma` per module and type (RF-referenced, runs aligned) and `rf_calib_modules.tsv:sigma_ns` per run × module (tier e100 / allE), new dump vs `output/module_timing_20261025/` and `output/rf_calib_20261026/`; plus the veto-referenced `module_timing_veto.tsv:sigma_<type>` for the H2 runs |
| σ_t per cluster | For each single-cluster event row of `result_clusters.tsv` (legacy path by default) of the 13 RF runs with usable RF calibration: (a) **seed time** t_seed = seed_time_ns − o_{seed,r}; (b) **energy-weighted cluster time** t_ewt = Σ_i E_i (t_i − o_{i,r}) / Σ_i E_i over the cluster members i (from the new `members` column) that have an offset o_{i,r} in `$O/rf_calib/<run>.json` and a finite t_i (the seed must be among them; members without an offset are dropped and counted). φ = fold(t − t_RF) with t_RF the nearest channel-A hit of the same TI event in `output/module_timing_20261025/rf/<split>.npz` and P = `rf_calib_runs.tsv:P_ns` of the run; σ_t = `plot_module_timing_rf.fit_unwrapped` on the 80-bin circular histogram (the estimator of F11). Reported per run, per run type (pooled) and per member-count tier (1, 2–4, ≥ 5), for both definitions, with N, σ ± err, μ and the fraction of members with an offset. Clusters whose seed lies in the centred 6 × 6 exclusion are absent from the TSV (benchmark rectangle cut), so the X17 central block is not covered here (the module-level S5 covers it) |

## 4. Stages

Gate for every heavy stage (S2 onwards): the running job has released the machine —
```
pgrep -fl "benchmark_hycal_timing_test_20261021.py|spline_template_extract.py" | grep -v pgrep | wc -l    # 0
ls output/scan_20261026_best50/scan_*/split_*/result.tsv | wc -l                                          # 20
```
If the count is not 0 after the 20 results exist (e.g. a new stage of that job started), stop and report; do not start S2–S6 alongside it.

### S0. Preparation and code (executor, ≈ 1 h; the benchmark edit only after the gate)

S0.1 Bank check (A5) and output tree:
```
for s in 0p73GeV_H2 2p24GeV_H2 2p24GeV_C 2p24GeV_Ta_X17 3p49GeV_H2; do test -f $TN/$s/templates_$s.npz && head -1 $TN/$s/summary.tsv | tr '\t' '\n' | grep -c '^osc_escalated$'; done   # five lines "1"
test ! -e $TN/esc_stage && echo moved
mkdir -p $O/{pileup/logs,gains/{dumps/b50_percrystal,dumps_extra,logs},scan/{configs/posmaps,logs,compare,identity},module_resolution,module_timing,rf_calib,cluster_timing,compare} $S
```
S0.2 `module_timing_dump.py`: add `--templates-root DIR` (default `output/templates_spline`); build `_T`/`RUN_CFG` from it at run time (the per-split subprocess must receive the same value). Default behaviour unchanged.
S0.3 `rf_time_calib.py`: add `--dump-dir`, `--rf-dir`, `--xcheck` with the current constants as defaults; thread them to `prep`/`XCHECK` reads. Default behaviour unchanged.
S0.4 `benchmark_hycal_timing_test_20261021.py` (after the gate): flag `--cluster-times` (store_true, default off; error exit 2 if given without `--write-cluster-tsv`). `member_rec` additionally collects `mt = np.fromiter((h.time for h in hits), np.float64, len(hits))` and `t_seed = float(r.cluster.center.time)` and returns `(mi, me, r.hit.x, r.hit.y, evnum, mt, t_seed)` (appended, so `ex_[0..5]` keep their meaning; grep every consumer of the event tuples — the cluster-dump writer at ≈ 2650-2730 and the TSV writer — for fixed-arity unpacking and keep them working). With the flag on, the TSV header gains `\tevnum\tseed_time_ns\tn_members\tmembers` and each row `\t{ex_[5]}\t{ex_[7]:.3f}\t{len(ex_[1])}\t` + `";".join(f"{idx_name[i]}:{e:.1f}:{t:.3f}" for i, e, t in zip(ex_[1], ex_[2], ex_[6]))`, where `idx_name` is the module-index → name list built once as the cluster-dump writer does (`mod_name = np.array([m_.name for m_ in mods_all])`, ≈ lines 2706-2709); **not** `_mod_name`, which maps PrimEx IDs (1–900 → G, 1001–2152 → W) while `h.index` is the module index. With the flag off, byte-identical output.
S0.5 Identity check of S0.4 (two jobs, ≈ 7 min, parallel): take `output/resolution_fit_20261026/r3_splits/scan_025308/split_00100/job.cmd`, replace the output stem (argv[3]) by `$O/scan/identity/off/result` and run; the same with `--cluster-times` appended into `$O/scan/identity/on/`. Acceptance: `off/result.tsv` equals `…/r3_splits/scan_025308/split_00100/result.tsv` except the `elapsed_s` header line and `off/result_clusters.tsv` is byte-identical; `on/result.tsv` identical to `off/`, `on/result_clusters.tsv` has 16 columns whose first 12 equal `off/`, `n_members ≥ 1`, `100 < seed_time_ns < 200`, and the seed module appears in `members` with `t == seed_time_ns` in every row (a 10-line check in `$S`).
S0.6 New scripts (no machine): `cluster_timing.py` (§6) and `best50_compare.py` (§6). Self-check each on existing files (`best50_compare.py pileup` on A vs A must give Δ = 0 everywhere; `modres`/`timing`/`rfcalib` on a file against itself likewise).

### S1. Pile-up detection efficiency, new bank vs production (≈ 2 min machine; may run before the gate)

```
for SL in 0p73GeV_H2 2p24GeV_C 2p24GeV_H2 2p24GeV_Ta_X17 3p49GeV_H2; do        # arm B
  python3 benchmark_pileup_mf.py output/pileup_eff_20261024_fix/$SL/synth --templates $TN/$SL/templates_$SL.npz --out $O/pileup/B_prodsynth_newbank/$SL --tolerance-ns 4 --window-rule 40:0.16 > $O/pileup/logs/B_$SL.log 2>&1; done
# arm C/D synth: the recorded generator command per slug with --template $TN/… and --out-dir $O/pileup/synth_newbank/$SL/synth
#   bases: 0p73GeV_H2 $E/prad_024590.evio.00010; 2p24GeV_C $E/prad_025308.evio.00050; 2p24GeV_H2 $E/prad_025320.evio.00050; 3p49GeV_H2 $E/prad_024503.evio.00010; 2p24GeV_Ta_X17 $E/prad_025800.evio.* (all 8 splits)
python3 pileup_generator.py $E/prad_025320.evio.00050 --template $TN/2p24GeV_H2/templates_2p24GeV_H2.npz --out-dir $O/pileup/synth_newbank/2p24GeV_H2/synth --material PbWO4 --n-per-config 50 --height-min 500 --t0-min 22 --residual-veto-sigma 999 --noise-model none --post-pad 20 --seed 12345 --force > $O/pileup/logs/gen_2p24GeV_H2.log 2>&1
for SL in …; do                                                                 # arms C and D
  python3 benchmark_pileup_mf.py $O/pileup/synth_newbank/$SL/synth --templates $TN/$SL/templates_$SL.npz --out $O/pileup/C_newsynth_newbank/$SL --tolerance-ns 4 --window-rule 40:0.16
  python3 benchmark_pileup_mf.py $O/pileup/synth_newbank/$SL/synth --templates $T0/$SL/templates_$SL.npz --out $O/pileup/D_newsynth_prodbank/$SL --tolerance-ns 4 --window-rule 40:0.16; done
python3 best50_compare.py pileup --ref output/pileup_eff_20261024_fix --new $O/pileup/B_prodsynth_newbank --arms C=$O/pileup/C_newsynth_newbank D=$O/pileup/D_newsynth_prodbank --out $O/compare/pileup
```
Before arm B confirm `benchmark_pileup_mf.py` writes only under `--out` (`grep -n "synth_dir\|open(.*\"w\|savez" benchmark_pileup_mf.py`); if it writes into the synth dir, copy the production synth to `$O/pileup/synth_prod/<slug>/synth` first. The generator's base-event selection is bank-independent (F2): C/D synth must have the same number of events per cell as A (3000 per slug) and the same base event list (`manifest.json`); report if not.
**Acceptance:** five `eff_<slug>.{md,json,png}` per arm; `$O/compare/pileup.md` with, per slug and block, the B grid of `eff_match` with Δ(B − A) in points, the B and A `fake_rate` grids, `n2_mf`, pooled Δ over all 60 cells and over the cells with A < 100, the C and D pooled Δ, and the md header lines (sha1, MIN_SEP_NS = 12, r_model, tolerance). `n2_wa` of B equals A cell by cell (same synth). No expectation on the sign of Δ; report.

### S2. Per-crystal gains with the new bank (reuse 5 runs; dump and fit 9 runs; ≈ 50 min machine)

```
for R in 024503 024590 025308 025320 025800; do                                 # reuse (F6)
  ln -s $PWD/output/percrystal_best50_20261026/dumps/scan_$R $O/gains/dumps/b50_percrystal/scan_$R
  ln -s $PWD/output/percrystal_best50_20261026/$R $O/gains/$R; done
cat > $O/scan/configs/cfgs_best50.json <<EOF
[{"tag":"b50_percrystal","time_calib":"off","W":16,"LO":100,"HI":200,"pulse":"mf","sp":"off","bench_args":"--mf-templates $TN/{runtype}/templates_{runtype}.npz --mf-scale 1.0 --cal-method percrystal --cluster-dump auto"},
 {"tag":"r3n","time_calib":"off","W":16,"LO":100,"HI":200,"pulse":"mf","sp":"off","bench_args":"--mf-templates $TN/{runtype}/templates_{runtype}.npz --mf-scale 1.0 --gain-file $O/gains/{run}/gains_{run}_all.tsv --pos-corr $O/scan/configs/posmaps/{runtype}/poscorr_r3_{runtype}_all.json --write-cluster-tsv --cluster-times"}]
EOF
P="--out $O/gains/dumps --parallel 12 --timeout-seconds 1800"
$G run --stage R  --configs-file $O/scan/configs/cfgs_best50.json --tags b50_percrystal --runs 024386,024512,024603,025228,025257,025260,025299 $P --dry-run | tail -3   # 57 jobs
nohup $G run --stage R  --configs-file $O/scan/configs/cfgs_best50.json --tags b50_percrystal --runs 024386,024512,024603,025228,025257,025260,025299 $P > $O/gains/logs/b50_R.out 2>&1
nohup $G run --stage RX --configs-file $O/scan/configs/cfgs_best50.json --tags b50_percrystal --runs 026138,026230 $P > $O/gains/logs/b50_RX.out 2>&1                     # 17 jobs
python3 benchmark_hycal_timing_test_20261021.py $E/prad_024386.evio.00060 $S/x24386_60 --calibrate --exclude-rowcol 15 20 15 20 --pre-window 100.0 200.0 --inner-radius 200 --min-per-module 30 --lms-corr off --hc-calib-file $DB/calibration_factor_fadc_3p5.json --window 16.0 --time-calib off --pulse-method mf --mf-scale 1.01 --mf-r-model 0.3 --mf-templates $TN/3p49GeV_H2/templates_3p49GeV_H2.npz --mf-scale 1.0 --cal-method percrystal --cluster-dump $O/gains/dumps_extra/024386_00060_dump.npz > $O/gains/logs/x24386_60.log 2>&1
for R in 024386 024512 024603 025228 025257 025260 025299 026138 026230; do
  X=""; [ $R = 024386 ] && X="$O/gains/dumps_extra/024386_00060_dump.npz"
  python3 percrystal_calib.py fit --dumps "$O/gains/dumps/b50_percrystal/scan_$R/split_*/result_dump.npz" $X --splits all --fixed-windows --out $O/gains/$R/ > $O/gains/logs/fit_$R.log 2>&1; done
python3 best50_compare.py gains --new $O/gains --ref output/percrystal_spl --out $O/compare/gains_vs_spl
python3 compare_scan_dirs.py --scans spl=output/scan_20261024/spl_percrystal new=$O/gains/dumps/b50_percrystal --gauss-window 2sigma --out $O/compare/pregain_peak
```
Run the extra dump as the 13th job alongside the harness (the harness holds 12). Stage R and RX one after the other (as R3).
**Acceptance:** 74 harness jobs + 1 extra, 0 failures (`b50_*.out`, `manifest.json`); 14 `gains_<run>_all.{tsv,json}` under `$O/gains/<run>/` (5 symlinked), all `converged`; **gains consistency** (`$O/compare/gains_vs_spl.md`): per run n_free within ± 5 of percrystal_spl, median ratio g_new/g_spl over modules free in both within 0.98–1.02, per-module RMS of the ratio < 1.5 % (10-24 reference for a bank change, plan_20261025 log 17:58: 0.23–0.63 %), the five reused runs' json `argv` contain `--splits all --fixed-windows` with the `--dumps` glob under `output/percrystal_best50_20261026/dumps/scan_<run>/`, and `grep -L templates_spline_best50_20261026 $O/gains/dumps/b50_percrystal/scan_<run>/split_*/job.cmd` prints nothing for those five runs; **pre-gain elastic-peak shift** (`$O/compare/pregain_peak.tsv`, arms A′ vs C = the raw template effect before any pooled gain): `eps_peak` per type reported with the pre σ/E paired Δ. Failures are reported, not loosened.

### S3. Full scan, R3 configuration with new templates + new gains, R3 position maps (118 splits, ≈ 70 min machine)

```
for s in 3p49GeV_H2 2p24GeV_H2 2p24GeV_C; do mkdir -p $O/scan/configs/posmaps/$s; ln -s $PWD/output/poscorr_20261025/hybrid/$s/poscorr_hybrid_${s}_all.json $O/scan/configs/posmaps/$s/poscorr_r3_${s}_all.json; done
for s in 0p73GeV_H2 2p24GeV_Ta_X17; do mkdir -p $O/scan/configs/posmaps/$s; ln -s $PWD/output/poscorr_20261025/$s/poscorr_${s}_all.json $O/scan/configs/posmaps/$s/poscorr_r3_${s}_all.json; done
cmp $O/scan/configs/posmaps/2p24GeV_Ta_X17/poscorr_r3_2p24GeV_Ta_X17_all.json output/poscorr_20261025/2p24GeV_H2/poscorr_2p24GeV_H2_all.json   # Ta = H2 bilinear (F4)
P="--out $O/scan --parallel 12 --timeout-seconds 1800"
$G run --stage R  --configs-file $O/scan/configs/cfgs_best50.json --tags r3n $P --dry-run | tail -3            # 93 jobs; check one printed command: $TN bank, $O gains, posmaps path, --cluster-times
nohup $G run --stage R  --configs-file $O/scan/configs/cfgs_best50.json --tags r3n $P > $O/scan/logs/r3n_R.out 2>&1
nohup $G run --stage RX --configs-file $O/scan/configs/cfgs_best50.json --tags r3n $P > $O/scan/logs/r3n_RX.out 2>&1   # 25 jobs
python3 compare_scan_dirs.py --scans r3=output/resolution_fit_20261026/r3_splits new=$O/scan/r3n --gauss-window 2sigma --out $O/scan/compare/new_vs_r3
python3 output/scan_20261025/pooled_delta.py new $O/scan/compare/new_vs_r3_splits.tsv > $O/scan/compare/new_vs_r3_pooled.md
```
**Acceptance:** 118 `result.tsv` + `result_clusters.tsv` (16 columns), 0 failures; `new_vs_r3.tsv` has n_splits = 118 common keys with `eps_yield` within 0.98–1.02 per type; the table of §3 (legacy and gated × G and CB × pre/post, Δ mean ± SE per type) in `$O/scan/compare/new_vs_r3.md` and the pooled md. Sanity: the post-gain `eps_peak` per type within ± 0.5 % (the pooled gains absorb the template-integral shift seen in S2's g0); otherwise report.

### S4. Module-by-module resolution (≈ 4 min)

```
python3 plot_module_resolution.py --scan-dir $O/scan/r3n --out $O/module_resolution/new --group-runs > $O/module_resolution/new.log 2>&1
python3 plot_module_resolution.py --scan-dir output/resolution_fit_20261026/r3_splits --out $O/module_resolution/r3_ref --group-runs > $O/module_resolution/r3_ref.log 2>&1
python3 best50_compare.py modres --new $O/module_resolution/new --ref $O/module_resolution/r3_ref --out $O/compare/modres
```
The reference is recomputed on the R3 mix because `output/module_resolution_20261025_hybrid/` used the hybrid maps for all 14 runs (F8).
**Acceptance:** both `summary/types_summary.tsv` present; `n_fitted` per type within ± 3 of each other; `$O/compare/modres.md` with per type the median / p16 / p84 of σ_M/μ_M (new, ref) and of the per-module Δ, the number of modules fitted in both, and the ten largest |Δ| with their N; `modres.tsv` per module.

### S5. Module timing with the new bank (dump ≈ 61 min at 12 parallel; plots + RF calibration ≈ 10 min)

```
python3 module_timing_dump.py --all --parallel 12 --templates-root $TN --out $O/module_timing/dump > $O/module_timing/dump_run.log 2>&1
python3 plot_module_timing.py --ref veto --dump-dir $O/module_timing/dump --out-dir $O/module_timing > $O/module_timing/run_veto.log 2>&1
python3 plot_module_timing_vs_integral.py --dump-dir $O/module_timing/dump --out-dir $O/module_timing > $O/module_timing/run_vs_integral.log 2>&1
python3 plot_module_timing_rf.py --dump-dir $O/module_timing/dump --rf-dir output/module_timing_20261025/rf --out-dir $O/module_timing > $O/module_timing/run_rf.log 2>&1
python3 rf_time_calib.py --dump-dir $O/module_timing/dump --rf-dir output/module_timing_20261025/rf --xcheck $O/module_timing/rf_summary.tsv --out-dir $O/rf_calib --jobs 12 > $O/rf_calib/run.log 2>&1; echo "rf_time_calib rc=$?" | tee -a $O/rf_calib/run.log
python3 best50_compare.py timing --new $O/module_timing/module_timing_rf.tsv --ref output/module_timing_20261025/module_timing_rf.tsv --out $O/compare/timing_rf
python3 best50_compare.py timing --new $O/module_timing/module_timing_veto.tsv --ref output/module_timing_20261025/module_timing_veto.tsv --out $O/compare/timing_veto
python3 best50_compare.py rfcalib --new $O/rf_calib/rf_calib_modules.tsv --ref output/rf_calib_20261026/rf_calib_modules.tsv --out $O/compare/rfcalib
```
**Acceptance:** 111 dump npz, `dump_run.log` "finished, 0 failed", every npz `templates` metadata under `$TN`; the four plot/TSV sets of F10 under `$O/module_timing/`; `rf_time_calib.py` exit 0 with 24386 flagged bimodal and 24503 `drift` as before (if a gate fails, keep the outputs and report the failing gate; do not change `--bimodal-ratio`); `$O/compare/timing_rf.md`: per type median σ (new, ref) and median Δσ with p16/p84 over modules with `fit_ok` in both, offset shift median; `rfcalib.md`: per run median Δoffset_ns (folded) and median Δsigma_ns over measured (run, module) pairs, and the fraction of modules whose tier changed.

### S6. σ_t per cluster (≈ 15 min: 13 reference jobs ≈ 7 min, offline ≈ 5 min)

S6.1 Production-bank reference on one split per RF run (13 jobs, parallel 12): for each run of F9 take the first split of its R3 pair list, copy `output/resolution_fit_20261026/r3_splits/scan_<run>/split_<n>/job.cmd` with argv[3] → `$O/cluster_timing/r3_ref/scan_<run>/split_<n>/result` and `--cluster-times` appended (the S0.5 helper), run with `/usr/bin/time -p`.
S6.2 
```
python3 cluster_timing.py --scan-dir $O/scan/r3n --rf-dir output/module_timing_20261025/rf --offsets-dir $O/rf_calib --out-dir $O/cluster_timing/new --jobs 12
python3 cluster_timing.py --scan-dir $O/cluster_timing/r3_ref --rf-dir output/module_timing_20261025/rf --offsets-dir output/rf_calib_20261026 --out-dir $O/cluster_timing/r3_ref --jobs 12
python3 best50_compare.py cluster --new $O/cluster_timing/new --ref $O/cluster_timing/r3_ref --out $O/compare/cluster_timing
```
**Acceptance:** `cluster_timing.tsv` rows for 12 usable runs (24386 skipped as unusable, reported) × {seed, ewt} × tiers, `fit_ok` for every per-type row, N ≥ 1000 per type and definition; per-split table joinable on the 13 common splits; `$O/compare/cluster_timing.md` with σ_t (seed, ewt) per type for new and ref on the common splits and for new on all splits, the fraction of members with an offset (expected > 0.95 at the H2 runs), and the member-count tiers. Sanity: the per-type seed-definition σ_t on the common splits is within 0.05 ns of the S5 `module_timing_rf.tsv` N-weighted median σ of the same type (same estimator, same reference); otherwise report.

### S7. Consolidated tables and the daily note (≈ 30 min executor)

```
python3 best50_compare.py summary --root $O --out $O/compare/summary.md
```
`summary.md` concatenates, in this order: bank facts (F1 counts), S1 table (B with Δ vs A; C/D pooled), S2 gains consistency and pre-gain peak shift, S3 σ_E/E table and pooled Δ, S4 module medians, S5 σ_t per module (RF, veto) and RF-calibration deltas, S6 σ_t per cluster, run times per stage (from the logs). Then write `notes/PbWO4_energy_resolution_20261026.md` (new file; the 10-25 note and the summary note belong to the running job) with sections: 1 Status (what ran, what was reused, wall times); 2 Bank (build settings, F1 counts, pointer to the running job's §4.10); 3 Pile-up detection efficiency (the table of §3 and the headline sentence in the form of note 20261023 §9.4); 4 Gains; 5 σ_E/E per run type; 6 Module-by-module resolution; 7 σ_t per module and RF calibration; 8 σ_t per cluster (definition of §3 verbatim); 9 Caveats and open items; 10 How to run (the commands of §4 as executed, with `(r)` on none — every command here is logged). No numbers in the note that are not in `$O/compare/`.

## 5. Files to create or modify

| File | Change |
|---|---|
| `benchmark_hycal_timing_test_20261021.py` | `--cluster-times` flag (S0.4); `member_rec` appends member times and seed time; TSV writer appends four columns when the flag is on. Defaults byte-identical (S0.5) |
| `module_timing_dump.py` | `--templates-root DIR` (default `output/templates_spline`) |
| `rf_time_calib.py` | `--dump-dir`, `--rf-dir`, `--xcheck` (defaults = current constants) |
| `cluster_timing.py` (new) | §6 |
| `best50_compare.py` (new) | §6 |
| `$O/scan/configs/cfgs_best50.json`, `$O/scan/configs/posmaps/<slug>/poscorr_r3_<slug>_all.json` (symlinks), `$O/scan/configs/run_*.sh` (the command blocks of §4 as run) | S2/S3 |
| `plans/plan_best50_pipeline_20261026.md` | this plan (progress log, todo) |
| `notes/PbWO4_energy_resolution_20261026.md` (new) | S7 |

Nothing is written under `output/templates_spline*/`, `output/percrystal_*/`, `output/scan_2026102*/`, `output/module_timing_20261025/`, `output/rf_calib_20261026/`, `output/pileup_*/`, `output/resolution_fit_*/` or `database/`.

## 6. Interfaces and data changes

- `result_clusters.tsv` with `--cluster-times`: columns 13–16 `evnum` (TI event number, int), `seed_time_ns` (the seed hit time as fed to the clusterer, ns, 3 decimals), `n_members` (int), `members` (`name:E_MeV:t_ns;…`, E 1 decimal, t 3 decimals). Size ≈ 2.5 × today (≈ 60 MB per split, ≈ 7 GB for 118 splits). Without the flag: unchanged.
- `cluster_timing.py --scan-dir DIR --rf-dir DIR --offsets-dir DIR --out-dir DIR [--path legacy|gated] [--jobs N] [--min-n 300]`: reads `DIR/scan_*/split_*/result_clusters.tsv` (requires the 16-column form), `rf/<run6>_<split>.npz` (6-digit run, the naming of `output/module_timing_20261025/rf/`), `<offsets-dir>/<run>.json` (5-digit run) and `rf_calib_runs.tsv` (`P_ns`, `usable`); imports `fold`, `_bin`, `NB`, `fit_unwrapped` from `plot_module_timing_rf.py` and `TYPE_OF_RUN` from `plot_module_timing.py`. Writes `cluster_timing.tsv` (`run type path definition tier N P_ns mu_ns sigma_ns sigma_err_ns fit_ok median_n_members frac_members_with_offset`), `cluster_timing_splits.tsv` (the same per split, for pairing), `cluster_timing.pdf` (one page per run type: folded histograms with fits for both definitions), `run.log`. Rows with a missing RF event, a seed without an offset, or a non-finite time are counted in the log and skipped.
- `best50_compare.py` subcommands, each writing `<out>.md` and `<out>.tsv`: `pileup --ref DIR --new DIR [--arms L=DIR …]` (joins `<slug>/eff_<slug>.json` grids per block); `gains --new DIR --ref DIR` (joins `gains_<run>_all.tsv` — columns `path module index row col g sigma_g n_seed n_member D_ii prior_frac g_seedindexed used`, `#` header lines with g0/n_free/n_prior — on (path, module), restricted to modules free in both; the executor takes the free/prior definition from `percrystal_calib.py`, where `prior_frac`/`used` are written); `modres --new DIR --ref DIR` (joins `<run>_modules.tsv` per type / `types_summary.tsv`); `timing --new TSV --ref TSV` (joins on `module`, σ columns per type, `fit_ok` in both); `rfcalib --new TSV --ref TSV` (joins on (run, module), tiers e100/allE, folded Δoffset); `cluster --new DIR --ref DIR` (joins `cluster_timing_splits.tsv` on (run, split)); `summary --root $O --out FILE` (concatenates the md files in the §4 S7 order with stage wall times parsed from the logs). Numpy/pandas only; no plots except where stated.
- `module_timing_dump.py --templates-root`, `rf_time_calib.py --dump-dir/--rf-dir/--xcheck`: path flags; the dump npz `templates` metadata and the calibration JSON `time_definition` record the new paths automatically.

## 7. Runtime, parallelism and critical path

Estimates from the logs of the previous runs (F2, F5, F6, F9–F11); 12 jobs = the setting of every full run so far (18 cores, 48 GB; 12 × 2.8 GB = 34 GB), which leaves the machine usable. Disk: scan ≈ 7.5 GB (with `members`), dump 4.9 GB, rest < 1 GB; 78 GB free.

| Stage | Jobs | Per job | Parallel | Wall | Depends on |
|---|---|---|---|---|---|
| Running job (16:28: 10 of 20 benchmark jobs done, gains fitted, builds done) | — | ≈ 400 s | 10 | ends ≈ 16:40, then note edits | — |
| S0 edits (dump, calib, new scripts) | — | — | — | 1 h executor, 0 machine; S0.4 after the gate | — |
| S0.5 identity | 2 | 6 min | 2 | 7 min | gate |
| S1 pile-up (B, synth, C, D) | 5 + 5 + 10 | 6–8 s | 1 | 2 min | bank only (can run now) |
| S2 dumps 74 + 1, fits 9 | 75 + 9 | 400 s; fit 4 s | 12 + 1 | 45 min + 1 min; checks 3 min | gate, S0.4 |
| S3 scan | 118 | 348 s mean (38–564) | 12 | 60–75 min (R3: 61) + compare 2 min | S2, S0.5 |
| S4 module resolution ×2 | 2 | 1.5 min | 1 | 4 min | S3 |
| S5 dump; plots; calib | 111; 3; 1 | 393 s mean; 83/107/312 s; 80 s | 12 | 61 min + 10 min | S0.2, S0.3 (bank only) |
| S6 reference jobs; offline | 13; 2 | 400 s; ≈ 2 min | 12 | 7 + 5 min | S3, S5 |
| S7 tables, note | — | — | — | 30 min executor | all |

**Critical path (serial, 12 jobs):** gate ≈ 16:40 → S0.5 + S1 (16:50) → S2 (17:40) → S3 (18:55) → S4 (19:00) → S5 (20:10) → S6 (20:25) → S7 (21:00). ≈ 4 h 20 min wall from the gate, of which ≈ 3 h 10 min machine. Option (not default): start the S5 dump at `--parallel 6` together with S3 at `--parallel 8` (14 jobs, ≈ 40 GB) and finish it after S3 at 12 → saves ≈ 45 min. Stop points for a shorter day: after S3 (σ_E/E delivered), after S5 (σ_t per module delivered).

## 8. Risks and compatibility

| Risk | Mitigation |
|---|---|
| Editing the benchmark while the running job still launches it (F13) | The gate; S0.4 strictly after it; S0.5 identity before S2 |
| `member_rec` tuple consumers with fixed arity (cluster-dump writer) | Append only; grep every unpack site (critic check: the cluster-dump writer is index-based, the `--gain-file` rebuild keeps `t_[8:]`, other sites unpack with `*_`); S0.5 covers the TSV path |
| 25308 (2.24 C) has no timing dump, RF npz or offsets (F9) | S5/S6 cover the 13 RF runs as before; 25308 is in S2–S4 only; recorded in the note's caveats |
| `rf_time_calib.py` gates (ii), (v), (vi) may fail on the new bank | Outputs are kept; the failing gate is reported with its numbers; no threshold change |
| Energy-weighted time drops members without an offset (tier `pred` has an offset; `unresolved` W1123/W1124 do not) | Fraction reported per run; seed must have an offset |

## 9. Out of scope

- Adopting the best-50 bank into any default (`_T` of `module_timing_dump.py`, `--mf-templates auto`, the harness configs) or into `database/`.
- C++ changes (`prad2det`, `python/bind_det.cpp`); the hit and centre times are already exposed.
- New X17 data; the pile-up rates scan (`plot_pileup_examples.py --splits all`, 59 min; Open Question 1).
- Refitting the position maps on the new bank (`position_corr.py map/fit` on an S3-like scan without `--pos-corr`, ≈ 1.3 h machine + held-out): recorded as an option for a later day, after the template effect is known.
- The `--ref evt`, `--time-calib` and `nogate` variants of `plot_module_timing.py` (note 20261025 §9 ran four; only the veto-referenced one is in the brief; Open Question 3).
- Commits (agents never commit); edits to the 10-25 note, the summary note or plan_20261025 (owned by the running job).

## 10. Acceptance criteria and validation commands

- A1 Stages S1–S7 as in §4, each with its acceptance block; S2 reuse decision stated (five runs reused, nine derived) with the json `argv` check.
- A2 `find $O -maxdepth 1 -type d` shows exactly `pileup gains scan module_resolution module_timing rf_calib cluster_timing compare` (configs and position-map symlinks live under `scan/configs/`); `git status --porcelain output/` lists nothing outside `output/best50_20261026/` that this plan touched (the running job's own files excepted).
- A3 §7 numbers come from the cited logs; the executor records the actual wall time per stage in the progress log.
- A4 Sanity checks: S2 gains consistency and pre-gain peak shift; S3 `eps_yield`, post `eps_peak`; S5 24386 bimodal / 24503 drift reproduced; S6 seed-definition σ_t vs S5 module σ_t.
- A5 S0.1 prints five "1" and "moved".
- A6/A7 §9 and the preamble; `MPLBACKEND=Agg` exported in every shell; scratch only under `$S`.
- A8 Every command above is the recorded recipe with paths changed; the three script changes are path/flag additions with unchanged defaults, verified by S0.5 and the self-checks of S0.6.

Validation commands (run by the executor at the end):
```
ls $O/pileup/B_prodsynth_newbank/*/eff_*.json | wc -l                       # 5
ls $O/gains/*/gains_*_all.tsv | wc -l                                        # 14
ls $O/scan/r3n/scan_*/split_*/result.tsv | wc -l                             # 118
head -1 $O/scan/r3n/scan_025320/split_00050/result_clusters.tsv | tr '\t' '\n' | wc -l   # 16
grep -c . $O/scan/compare/new_vs_r3.tsv                                      # > 1; n_splits column = 118 in every row
ls $O/module_timing/dump/*.npz | wc -l                                       # 111
ls $O/rf_calib/*.json | wc -l                                                # 13
wc -l $O/cluster_timing/new/cluster_timing.tsv                               # > 1
test -s $O/compare/summary.md && test -s notes/PbWO4_energy_resolution_20261026.md
```

## Open Questions

1. Pile-up **rates** with the new bank (`plot_pileup_examples.py --jobs 12 --splits all --mf-templates "$TN/{runtype}/templates_{runtype}.npz" --dominant-min-mev 50 --mf-scale 1.0 --burst-veto 10 --window-rule 40:0.16 --out $O/pileup/rates`, 59 min at 12 jobs): run it after S7 or skip? Default in this plan: skip (the requested table is the efficiency grid, F2/F3).
2. Per-cluster σ_t reference: 13 production-bank splits (default, 7 min) or the full 118 (≈ 65 min, gives a paired per-split comparison on everything)?
3. The three extra `plot_module_timing.py` variants (`--ref evt`, `--time-calib`, `--veto-gate none`; ≈ 4 min): include for completeness of the §9-style PDFs? Default: only `--ref veto`.
4. Overlap option of §7 (S5 dump alongside S3, 14 jobs): allowed, or keep 12 jobs serial (default)?

## Todo

- [x] S0.1 bank check printed five "1" and "moved"; `$O` tree created
- [x] S0.2 `module_timing_dump.py --templates-root`
- [x] S0.3 `rf_time_calib.py --dump-dir/--rf-dir/--xcheck`
- [x] S0.6 `cluster_timing.py`, `best50_compare.py` written and self-checked
- [x] S1 arms B, C, D run; `$O/compare/pileup.md` written
- [x] Gate: running job's benchmark processes 0, 20 results present
- [x] S0.4 `--cluster-times` in the benchmark
- [x] S0.5 identity check (off byte-identical; on 16 columns, seed in members)
- [x] S2 symlinks for 5 runs; json `argv` check recorded (outputs: `output/percrystal_best50_20261026/`, 10-27)
- [x] S2 harness R (57) + RX (17) + extra dump; 0 failures (outputs: `output/percrystal_best50_20261026/dumps/`, 10-27)
- [x] S2 nine fits converged; `gains_vs_spl.md`; `pregain_peak.tsv` (outputs: `output/percrystal_best50_20261026/`, consistency in `logs/gains_vs_spl.txt`, 10-27; no `gains_vs_spl.md` or `pregain_peak.tsv` written; Ta gains redone with the n200 bank in `output/percrystal_ta_fix_20261028/n200/`, 10-28)
- [x] S3 posmaps symlinks; `cmp` Ta = H2 bilinear; dry-run command inspected (outputs: `output/resolution_fit_20261027_best50/logs/scan_jobs.txt`, 10-27)
- [x] S3 harness R (93) + RX (25); 118 results, 16-column TSVs (outputs: `output/resolution_fit_20261027_best50/r4_splits`, 10-27; Ta rescanned with the n200 bank in `output/resolution_fit_20261028_ta/n200/splits`, 10-28)
- [x] S3 `new_vs_r3.{md,tsv}` (118 common), pooled md; `eps_yield`, `eps_peak` recorded (outputs: `output/resolution_fit_20261027_best50/compare/r4_vs_r3.{md,tsv}`, `output/resolution_fit_20261027_best50/compare_vs_r3.md`, 10-27; Ta `output/resolution_fit_20261028_ta/compare/n200_vs_r3.{md,tsv}`, 10-28)
- [x] S4 new + r3_ref module resolution; `modres.md` (outputs: `output/resolution_fit_20261027_best50/module_resolution/{r4,r3_ref}`, table in `output/resolution_fit_20261027_best50/compare_vs_r3.md` instead of `modres.md`, 10-27; Ta `output/resolution_fit_20261028_ta/module_resolution/n200`, 10-28)
- [x] S5 dump 111 npz, 0 failed; `templates` metadata under `$TN`
- [x] S5 veto, vs_integral, rf plots (vs_integral not run, reduced scope); `rf_time_calib.py` exit code and gates recorded
- [x] S5 `timing_rf.md`, `timing_veto.md`, `rfcalib.md`
- [x] S6.1 13 reference jobs; S6.2 `cluster_timing` new + ref; `cluster_timing.md`; seed-vs-module sanity recorded
- [x] S7 (reduced: no `summary.md`; note §9.1/§9.4 filled from `$O/compare/`) `summary.md`; `notes/PbWO4_energy_resolution_20261026.md` written from it
- [ ] Progress log: actual wall times per stage; validation commands of §10 pasted with their output

## Progress log

- 2026-10-26 16:30 — plan drafted (planner). Running job state at 16:28: gains fitted for 5 runs, 10/20 benchmark jobs done, template builds finished.
- 2026-10-26 17:05 — critic round 1 (REVISE: C1 blocking — member names by module index, not `_mod_name`; C2–C8 applied: configs under `scan/configs/`, `b50_percrystal` mkdir, argv/job.cmd check, rf_time_calib exit code captured, 38–564 s, 5/6-digit run naming, `--slug` dropped). Scope-critic skipped on the coordinator's instruction; risks and alternatives compacted.
- 2026-10-26 18:35 — reduced scope executed (user: "do reduced"): S0 (no S0.4 gate wait needed), S1 3 min, S5 dump 73.6 min at 10 parallel + plots 3.7 min + rf_time_calib 85 s (gates PASS), S6 26 benchmark jobs (one split per RF run × 2 banks, production gains) 39 min at 6 parallel alongside the dump, cluster_timing 11 s each. S2–S4 not run. Note 20261026 §9.
