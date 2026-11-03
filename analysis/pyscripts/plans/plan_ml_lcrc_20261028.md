# Plan: ML pile-up work (M2+) on ANL LCRC Swing, with the Mac keeping the benchmark jobs
**Status (2026-10-29): ON HOLD** — user decision: the LCRC/Swing work waits until the resource situation (allocation, project name, access) is clear. Local results since the plan was written: M2b passes all e-p criteria (`output/ml_pileup_20261027/M2b/`); M3 Ta per run done (`output/ml_pileup_20261027/M3/`) — Ta MF singles are 6–18 % contaminated by unresolved pile-up, so arm (c) must use cleaned negatives and a per-run operating point set by close-pair excess (< 0.5 %/channel); arm (d) sweep is not needed. S0 questions remain open.


**Date:** 2026-10-28. **Status:** draft for user approval; nothing executed on LCRC. No LCRC access was used to write this plan (web docs only); every LCRC-specific name below carries a docs URL or "confirm on login".
**Reads first:** `plans/plan_ml_pileup_20261027.md` (user decisions, U0–U7, M1–M3 naming); `notes/PbWO4_energy_resolution_20261027.md` §8, §10; `notes/PbWO4_energy_resolution_20261028.md` §5; `output/ml_pileup_20261027/U0/report.md`, `U5/report.md`; `output/ml_pileup_20261027/M1/tools/m1_{data,model,train,eval}.py`; root `README.md` "Building", root `CMakeLists.txt` (BUILD_PYTHON), `python/CMakeLists.txt`.
Paths are relative to `analysis/pyscripts/` unless absolute. Standing constraints carried over: Python only; no commits or pushes by agents; `database/` never written; Ta runs 25800/26138/26230 per run, never pooled; outputs under `output/ml_pileup_20261027/`; no package installs on the Mac without approval.

**Correction to the brief.** Swing is scheduled by **PBS Pro, not Slurm** (https://docs.lcrc.anl.gov/swing/running-jobs-swing/ ; https://lcrc.anl.gov/systems/swing). S3 therefore gives PBS templates (`qsub`, `-A`, `-q`, `-l select=1:ngpus=N`); there is no partition/QoS/`--gres`/`srun`. Everything else in the brief stands.

## 1. Goal and approach

Goal: move the GPU- and scale-bound parts of the ML pile-up programme (M2+: real-waveform training, Ta per run, hyper-parameter sweeps; later the ML Ta/three-cluster work) to LCRC Swing, where the user runs Claude Code agents directly, while the Mac keeps everything that needs the full `output/` tree (123 GB), the EVIO splits (257 GB, 130 files, 2.1 GB each) or `prad2py`.

Approach: ship a self-contained data bundle of npz files plus the git branch; train on Swing; sync results back into the mirrored `output/ml_pileup_20261027/M2/` layout. Rejected alternatives, one line each: (i) building `prad2py` on Swing and shipping EVIO for raw-waveform generation — 2.1 GB per split, needs the C++ tree, FetchContent from GitHub and the read-only `database/`; held in reserve (S2.3) for the case that the npz bundle is too small; (ii) running the benchmark/R5 σ_E/E jobs on Swing — they are CPU-bound, need all 118 splits and `prad2py`, and give nothing from a GPU.

Why M2 is not information-limited at 8 ns: the U0 χ² fit reaches 100 % at Δt 8 ns ratio 1 (target > 80 %) on synthetic data (20261027 §8), so the open problems that justify GPU scale are the real-data ones from U5 (20261028 §5): the U0 thresholds give a fake rate 70–170× the genuine close-pair rate on e-p, the operating points that cure it cost the ratio-0.3 efficiency (22–46 % at 8 ns), and the Ta runs have 12–18 % of channels with ≥ 2 pulses (§10) with a fake-subtracted close-pair fraction of 1.21 / 1.31 / 1.87 % per run. M2 on Swing targets those, not the synthetic ceiling.

## 2. Facts this plan rests on

Local (measured 2026-10-28):
- Git: branch `waveform_jz` at 287af16 = `origin/waveform_jz` (pushed); remote `git@github.com:Chao1009/prad2evviewer.git` (SSH). Uncommitted: 7 modified, 1 deleted, 1 untracked (notes, benchmark scripts, `position_corr.py`) — not on origin.
- `database/` is git-tracked (185 files, 55 MB) and arrives with the clone; `output/` is almost entirely untracked (70 tracked files, none under `ml_pileup_*`), so the M1 tools in `output/ml_pileup_20261027/M1/tools/` travel by rsync, not git.
- Mac ML venv `~/.venvs/prad-ml-cpu`: Python 3.13.13, torch 2.15.0 (CPU), scikit-learn 1.9.1, numpy 2.5.4, scipy 1.18.1, matplotlib 3.11.2.
- Bundle candidates and sizes: best-50 banks `output/templates_spline_best50_20261026/<slug>/templates_<slug>.npz` 397 MB (5 slugs); best-50 clean pulse pools `<slug>/pulses/` 45 + 36 + 257 + 81 + 112 = 531 MB (0.73 H2 / C / 2.24 H2 / Ta / 3.49 H2); `doublet_audit/` 35 MB (M1 removes its rows); Ta n200 bank `output/templates_spline_ta_fix_20261028/n200/` 81 MB; frozen 60-cell test set `output/pileup_eff_20261024_fix/` 17 MB; `output/pileup_rates_20261027/` 823 MB (`dump/` 25 npz, MF pulse lists of 25 splits, no samples) plus TSVs; U0 740 KB, U5 528 KB; M1 352 MB on disk (data generated 2026-10-28: train 163 / real_2p24 160 / test_synth 32 / val 14 MB, plus tools, model, tables). The U5 decode caches (3.00 M real waveforms of the 25 `pileup_rates` splits, incl. 9 Ta splits, `uint16[N,100]` + MF pulse-1 metadata) are 1.4 GB and live in the session scratchpad (`twopulse_realdata.py` `SCRATCH`), not in `output/`; they are the only real-waveform set that exists without EVIO.
- Hard-coded scratch paths: `m1_data.py` `SCR_M1` and `twopulse_realdata.py` `SCRATCH` point at the Mac session scratchpad; both must become overridable before anything runs on Swing (S1.4).
- `prad2py` build (root `README.md`, `CMakeLists.txt`, `python/CMakeLists.txt`): CMake ≥ 3.14, C++17, Python ≥ 3.7 with headers, ZLIB, Threads; pybind11 v2.13.6, evio, nlohmann_json, websocketpp, asio are fetched from GitHub by CMake (needs outbound network at configure time); module lands in `build/python/`.

Swing (researcher-worker, web, 2026-10-28; base https://docs.lcrc.anl.gov/):
- Nodes: 6 DGX nodes, 8× A100 each (gpu1–4, gpu6: 40 GB; gpu5: 80 GB), 2× EPYC 7742 (128 cores), 1 TB RAM (2 TB on gpu5), InfiniBand HDR; 2 login nodes `gpulogin1/2.lcrc.anl.gov` without GPUs (`swing/hardware-overview-swing/`).
- Scheduler PBS Pro: queues `gpu` (default, 5 nodes), `gpu-large` (gpu5), `backfill` (max 4 h, for exhausted allocations); `-A <project>` and walltime mandatory; default walltime 1 h, max 24 h; GPUs by `-l select=1:ngpus=N`, N ∈ {1,2,4,8}; each GPU = 1/8 of a node (16 cores, ≈ 128 GB), nodes shared; per-user limits 4 running / 10 queued jobs, 16 GPUs; interactive `qsub -I -l select=1:ngpus=1 -l walltime=01:00:00 -q gpu -A <project>`; batch scripts `#!/bin/bash -l` with `#PBS` lines and `cd $PBS_O_WORKDIR`; `qstat`, `qdel`, `sbank-list-allocations -p <project>` (`swing/running-jobs-swing/`, `running-jobs-at-lcrc/pbs-pro/`). Charging in node-hours, 1 node-hour = 8 GPU-hours (`allocation-management/allocations/`); confirm on login that a 1-GPU job bills 1/8 node-hour per hour.
- Software: Lmod (`module avail/spider/load`); Anaconda is no longer supported, use `module load miniforge3/<version>` (`using-software/software-specific-guides/conda/`); conda envs in the home directory per the docs; no PyTorch recipe documented, the TensorFlow guide uses `module load singularity` + NGC container with `singularity exec --nv` (`using-software/software-specific-guides/tensorflow/`). CUDA module names and the miniforge3 version: confirm on login (`module spider cuda`, `module spider miniforge3`).
- Storage: `/home/<user>` 100 GB soft / 1 TB hard; `/lcrc/project/<project>` 1 TB default, backed up; `/lcrc/globalscratch` GPFS, not backed up, cleaned on maintenance days without a stated age; node-local `/scratch` 15–20 GB RAM-backed tmpfs cleared at job end (the hardware table also lists a 14 TB SSD per node — conflict; `df -h /scratch` on login) (`data-management/filesystem-and-storage/disk-quota/`, `.../file-systems/`). `lcrc-quota` shows usage.
- Access and transfer: no direct SSH; ProxyJump through `logins.lcrc.anl.gov` to `swing.lcrc.anl.gov`, Duo MFA (`swing/getting-started-swing/`); Globus collection "LCRC Improv DTN" (`data-management/data-transfer/using-globus/`), SCP/SFTP "not recommended" for large transfers. Outbound internet from login and compute nodes (git clone, pip): not documented — confirm on login.
- Rules: no compute on login nodes (cgroup-throttled; `best-practices-and-policies/shared-resource-policy/`); maintenance the second Monday of each month from ≈ 08:00 Central, jobs that would overlap are held (`.../monthly-maintenance-day/`).

## 3. Steps

Owner key: **user** = only the user does it (logins, credentials, git pushes, allocation); **local agent** = Claude Code on the Mac (executor); **LCRC agent** = Claude Code on a Swing login node, submitting PBS jobs.

### S0 Prerequisites (owner: user)

| # | Item | Verification |
|---|---|---|
| S0.1 | LCRC project name for `-A` (allocation with node-hours on Swing) | `sbank-list-allocations -p <project>` on login shows a Swing allocation with hours left (`sbank-list-projects` to list your projects; the per-user lookup form is "confirm on login") |
| S0.2 | Swing login works: `~/.ssh/config` with `ProxyJump logins.lcrc.anl.gov` → `swing.lcrc.anl.gov`, Duo | `ssh swing hostname` prints `gpulogin1` or `gpulogin2` |
| S0.3 | Claude Code installed on Swing in `$HOME` (login node only; it submits jobs, it never trains) and outbound HTTPS works from the login node | `claude --version`; `curl -sI https://api.anthropic.com \| head -1`; `curl -sI https://github.com \| head -1` |
| S0.4 | GitHub reachable from Swing: an SSH key for `git@github.com:Chao1009/prad2evviewer.git` or an HTTPS token; decide what to push first (the 10 uncommitted local changes are not on origin; `waveform_jz` itself is pushed at 287af16) | `git ls-remote origin waveform_jz` from Swing returns the commit the user intends |
| S0.5 | Transfer path chosen: Globus (needs Globus Connect Personal on the Mac + "LCRC Improv DTN"; recommended for the 2–3 GB bundle and mandatory if EVIO splits ever move) or `rsync -e ssh` through the ProxyJump (fine for ≤ 3 GB, the docs discourage it for large sets) | a 100 MB test file arrives with matching `shasum` |
| S0.6 | Landing directories: code `/lcrc/project/<project>/<user>/prad2evviewer` (backed up, 1 TB); data bundle and results `/lcrc/project/<project>/<user>/prad_ml/`; intermediates `/lcrc/globalscratch/<user>/prad_ml/` (not backed up, may be purged) | `lcrc-quota` shows headroom ≥ 50 GB |

### S1 Division of labour and the interface (owner: local agent prepares; user transfers)

S1.1 Stays on the Mac: the benchmark and R5/σ_E/E jobs (`benchmark_hycal_timing_test_20261021.py`, 118 splits), U5b `--pulse-list` hook, U6, U7 paper material, every EVIO decode (`prad2py`), `pileup_rates.py dump`, `twopulse_realdata.py decode`, template extraction, the M1 baseline (already running, CPU). Goes to Swing: M2 (a)–(d) of S4, Ta per run, HP sweeps, later `output/ml_cluster_*` V1/V1b at scale and the X17 three-cluster work (same bundle mechanism; the R3 per-cluster dumps are not sized here — measure before that step).

S1.2 Bundle `output/ml_lcrc_bundle_20261028/` (local agent builds it with `rsync -a --relative` so the paths below are preserved relative to `analysis/pyscripts/`; a `MANIFEST.tsv` with path, bytes, sha1 of every file; never `database/`, never EVIO by default):

| Item | Source path | Size | Needed by |
|---|---|---|---|
| best-50 banks | `output/templates_spline_best50_20261026/<slug>/templates_<slug>.npz` (5) | 397 MB | all |
| clean pulse pools | `output/templates_spline_best50_20261026/<slug>/pulses/*.npz` (5 slugs, 119 files) | 531 MB | (b), (c), (d) |
| doublet audit | `output/templates_spline_best50_20261026/doublet_audit/` | 35 MB | M1 data stage |
| Ta n200 bank | `output/templates_spline_ta_fix_20261028/n200/templates_2p24GeV_Ta_X17.npz` (+ `summary.tsv`, `check/`) | 81 MB | (c) |
| frozen test set | `output/pileup_eff_20261024_fix/` | 17 MB | all |
| pile-up candidate lists | `output/pileup_rates_20261027/` (`dump/` 26 npz + TSVs + `report.md`) | 823 MB | (a), (c) |
| real waveforms (U5 decode caches) | promote scratch `U5/` to `output/ml_pileup_20261027/U5/decode/` first (S1.4) | 1.4 GB | (a), (c) |
| U0, U5 results (thresholds, OPs) | `output/ml_pileup_20261027/U0/`, `U5/` | 1.3 MB | (a), (c) |
| M1 everything | `output/ml_pileup_20261027/M1/` (tools, data, model, tables once done) | ≈ 370 MB | all |
| **Total** | | **≈ 3.7 GB** | |

Optional, only after the S2.3 decision: EVIO splits, 2.1 GB each — the 25 `pileup_rates` splits are 53 GB, the 9 Ta splits 19 GB (`rates_by_type.tsv` lists them); via Globus to `/lcrc/project/<project>/<user>/evio/`.

S1.3 Comes back (LCRC agent → `/lcrc/project/.../prad_ml/output/ml_pileup_20261027/M2/<exp>/`, user rsyncs → Mac `output/ml_pileup_20261027/M2/<exp>/`, same relative layout): `model/*.pt` + `config.json` (M1-size CNNs are < 1 MB), `tables/*.tsv`, `figures/*.png`, `logs/*.log`, `report.md`, `manifest.json` (input sha1s, git commit, PBS job ids, GPU type, wall time). Never `database/`, never git-tracked files, never raw data copies, size per experiment ≤ 100 MB (checkpoint chains stay in globalscratch).

S1.4 Code changes before shipping (local agent, Python only, byte-identical defaults on the Mac): (i) `twopulse_realdata.py` `SCRATCH` and `m1_data.py` `SCR_M1` read `os.environ.get("PRAD_ML_SCRATCH", <current default>)`; (ii) `twopulse_realdata.py` gains a `promote` sub-command (or a one-off script) copying the decode caches into `output/ml_pileup_20261027/U5/decode/` with a manifest; (iii) `m1_train.py` and `m1_eval.py` pick `device = cuda if available else cpu` and take `--threads N` for `torch.set_num_threads` (currently fixed at 4), with CPU results unchanged on the Mac (one smoke run, same loss curve to 1e-6). Verification: `OMP_NUM_THREADS=4 ... m1_train.py --smoke` and `m1_eval.py --smoke` before and after the edit give identical TSVs.

### S2 Environment on Swing (owner: LCRC agent, inside an interactive job for the GPU test; user approves the env creation)

S2.1 Conda env (recommended):
```
module spider miniforge3                      # pick the version shown; confirm on login
module load miniforge3/<version>
conda create -y -n prad-ml-gpu python=3.13    # home dir, as the LCRC conda page recommends; ~6–8 GB
source activate prad-ml-gpu                   # the docs use `source activate`
pip install torch==2.15.0 --index-url https://download.pytorch.org/whl/<cuXYZ>   # pick the CUDA tag matching the driver shown by `nvidia-smi` on a compute node; confirm on login
pip install scikit-learn==1.9.1 numpy==2.5.4 scipy==1.18.1 matplotlib==3.11.2
```
Versions pinned to the Mac venv so Mac and Swing results differ only by device. Test in an interactive job (S3.1): `python -c "import torch; print(torch.__version__, torch.cuda.is_available(), torch.cuda.get_device_name(0))"` prints `2.15.0 True NVIDIA A100...`. Fallback if pip wheels cannot be downloaded from compute nodes: install on the login node (allowed for small installs; compiling is not) or use the NGC PyTorch container through `module load singularity` as the LCRC TensorFlow guide does (image tag: confirm on login).

S2.2 Code: `git clone` the repo into `/lcrc/project/<project>/<user>/prad2evviewer`, `git checkout waveform_jz`; make the clone's `analysis/pyscripts/output` a symlink to `/lcrc/project/<project>/<user>/prad_ml/output` so there is one `output/` tree on Swing (the bundle lands there once, M2 tools/pbs/results are written there and come back with S3.3); rsync the bundle into it (preserving relative paths); `export PRAD_ML_SCRATCH=/lcrc/globalscratch/<user>/prad_ml`, `export MPLBACKEND=Agg`; always `qsub` from `analysis/pyscripts`. Verification: `python -c "import twopulse_fit, mf_pulse, twopulse_realdata"` from `analysis/pyscripts/` succeeds without `prad2py` (`pileup_rates.py` imports `prad2py` only inside `process_split`; `_common.py` is not needed on Swing); `sha1sum -c MANIFEST.tsv` passes.

S2.3 `prad2py` on Swing: **not for M2** (recommendation). Build only if the real-waveform set must grow beyond the 3.00 M channels of the U5 caches (e.g. more than 120 k channels per split or more splits for Ta). If so, in an interactive job (login nodes are throttled):
```
module load cmake gcc            # names: confirm on login (`module spider cmake`, `module spider gcc`)
cd /lcrc/project/<project>/<user>/prad2evviewer
cmake -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_PYTHON=ON -DBUILD_ANALYSIS=OFF -DBUILD_GUI=OFF -DWITH_ET=OFF -DEVIO_SOURCE=fetch -DPython_EXECUTABLE=$(which python)
cmake --build build -j 16
export PYTHONPATH=$PWD/build/python:$PYTHONPATH
export PRAD2_DATABASE_DIR=$PWD/database        # read-only use; the clone's tracked copy
python -c "import prad2py; print(prad2py.__version__)"
```
Needs outbound GitHub at configure time (FetchContent: pybind11 v2.13.6, evio, nlohmann_json, websocketpp, asio) and ZLIB headers; if compute nodes have no internet, configure on the login node (network only, no build) and build in the job. Then `twopulse_realdata.py decode` on the shipped EVIO splits. This path is an Open Question, not a default.

### S3 PBS templates (owner: LCRC agent writes them under `output/ml_pileup_20261027/M2/pbs/`; they are results-side files and come back with the sync)

S3.1 Interactive debug (1 GPU, 1 h; `-q gpu`):
```
qsub -I -A <project> -q gpu -l select=1:ngpus=1 -l walltime=01:00:00
```
Use it for S2.1's CUDA test, a `--smoke` run of `m1_train.py`, and `nvidia-smi` / `df -h /scratch` / `nproc` to confirm the 1/8-node share (expect 16 cores).

S3.2 Batch training template `train_m2.pbs` (one GPU; the models are < 10⁶ parameters, so one A100 is the unit; sweeps use many 1-GPU jobs, not multi-GPU):
```
#!/bin/bash -l
#PBS -A <project>
#PBS -q gpu
#PBS -l select=1:ngpus=1
#PBS -l walltime=04:00:00
#PBS -N m2_<exp>
#PBS -j oe
#PBS -o /lcrc/project/<project>/<user>/prad_ml/output/ml_pileup_20261027/M2/<exp>/logs/
cd $PBS_O_WORKDIR
module load miniforge3/<version>; source activate prad-ml-gpu
export OMP_NUM_THREADS=16 MKL_NUM_THREADS=16 MPLBACKEND=Agg
export PRAD_ML_SCRATCH=/lcrc/globalscratch/<user>/prad_ml
RUN=/lcrc/project/<project>/<user>/prad_ml/output/ml_pileup_20261027/M2/<exp>
python output/ml_pileup_20261027/M2/tools/m2_train.py --exp <exp> --out $RUN --threads 16 \
       --checkpoint-every 1 --resume --max-minutes 220 --seed <seed>
python output/ml_pileup_20261027/M2/tools/m2_eval.py --exp <exp> --out $RUN
```
`mkdir -p $RUN/logs` before `qsub` (PBS cannot deliver the job log otherwise). `--max-minutes` stays under walltime so the final checkpoint and `report.md` are written; `--resume` restarts from the last checkpoint if the job is requeued. Sweeps: PBS job arrays (`#PBS -J 0-35`, index in `$PBS_ARRAY_INDEX`; array support on Swing and whether subjobs count toward the 4 running / 10 queued limits: confirm on login); with 1-GPU jobs the concurrency is min(4 jobs, 16 GPUs) = 4, so a 36-run sweep takes 9 waves (or 4 configs per 4-GPU job, 3 waves). `backfill` only when the allocation is exhausted (4 h cap). `gpu-large` is not needed (datasets ≤ 2 GB fit a 40 GB card many times over). Verification: a `--smoke` submission of `train_m2.pbs` exits 0 and writes `logs/` and `report.md`.

S3.3 Results sync (user, from the Mac): `rsync -av --exclude 'ckpt_*' swing:/lcrc/project/<project>/<user>/prad_ml/output/ml_pileup_20261027/M2/ output/ml_pileup_20261027/M2/` then `python3 -c` sha1 check against each `manifest.json`.

### S4 M2 experiments on Swing (owner: LCRC agent; only the ones that need GPU or scale)

Common protocol: the M1 data protocol (module and pulse hold-out, audit rows removed, even/odd split files, Ta never pooled), the M1 eval metrics (`eff_dt4`, `eff4_u0`, fake per channel at the synthetic WP = 2 % fire on calibration singles and the real WP = 2 % fire on real half-0 singles, Wilson 68 %), 3 seeds, a model card per run. Target cells: Δt ∈ {4, 8, 12, 16, 20, 24, 32} ns (onset separation), ratio {0.3, 0.5, 1, 2}. Headline acceptance for every arm: **eff ≥ 80 % at Δt 8 ns ratio 1** on M1's `data/test_synth.npz` cells (onset Δt; the `pileup_eff_20261024_fix` set is in the pileup_3way label convention and has no 8 ns onset cell), at a synthetic fake ≤ 0.03/event; the arm-specific criteria below address the real-data transfer problem.

| Arm | What | Data | Accept | Depends on |
|---|---|---|---|---|
| (a) real singles as negatives / domain adaptation | M1 architecture trained with the real U5-cache channels that have no MF secondary and no fit detection above the U5 operating point as (noisy) negatives, weighted 1:1 with synthetic singles; second variant: a domain-adversarial head (real vs synthetic) on the CNN trunk | M1 train/val + U5 decode caches for 2.24 H2 (4 splits, 480 k channels; half 0 for the WP, half 1 for the rate) | real-WP fake detections with |dt| < 16 ns on half 1 < 2× the MF-expected genuine rate (0.042–0.045 %/channel, 20261028 §5), with ratio-0.3 efficiency at 8 ns > the U5 OP's 45.6 %; synthetic 8 ns r1 ≥ 80 % kept; pre-pulse symmetry (N(−16<dt<0) ≈ fake) reported | M1's `real_fake.tsv`: if M1 already meets the real-WP criterion, (a) is a one-seed confirmation; if not, (a) is first |
| (b) learned deconvolution / Δt–ratio regression at scale | plan_ml_pileup §3 "M2": 5–10 unrolled ISTA iterations with a learned per-module kernel bank initialised from the best-50 spline, 1 ns output grid; plus the M1 regression heads trained on all five run types jointly (run-type input), 36 → 119 pulse files | all best-50 pulse pools + banks, overlays built on Swing with the M1 generator at 10× M1's rows (≈ 6.5 M waveforms, ≈ 1.5 GB) | Δt RMS ≤ 0.3 ns and A1 bias < 1 % at 8 ns r1 (U0 fit: 0.2 ns, < 1 %); ratio-0.3 efficiency at 8 ns ≥ 90 % at the synthetic WP; no run type below 80 % at 8 ns r1 | M1 synthetic headline (if M1 < 80 % at 8 ns r1, run (d) first) |
| (c) Ta per run | model from the best of (a)/(b) fine-tuned per run on the Ta pulse pool + Ta real singles, evaluated per run on the U5 Ta caches (360 k channels per run); Ta n200 bank as the kernel init | Ta pools (81 MB), n200 bank, U5 caches for 25800 / 26138 / 26230 (3 splits each) | per run: fake-subtracted close-pair fraction (|dt| < 16 ns) with its error, to compare with U5's 1.21 / 1.31 / 1.87 %; real-WP fake below the U5 OP residual 0.69 / 0.77 / 1.14 %; RF on-comb fraction at 16–80 ns reported against the MF-complete 36.0 % reference; three separate tables, never pooled | (a) or (b) accepted |
| (d) hyper-parameter sweep | 12 configs (conv depth 4–6, width 32–64, kernel 5–7, dilation on/off, lr {1e-3, 2e-3}, w-reg {0.5, 1, 2}, dropout {0.1, 0.2}) × 3 seeds = 36 one-GPU jobs ≤ 2 h each | M1 data (or (b)'s set) | one table; the chosen config has val eff within 2 points of train and the best 8 ns r1 at the synthetic WP; sweep cost ≤ 72 GPU-h = 9 node-hours | none; runs first if M1 misses the 80 % headline |

Compute (estimate, unsourced): M1 on the Mac is ≈ 1–3 h per 20-epoch run on ≈ 0.6 M waveforms at 4 threads; on one A100 with the data resident on the GPU the same run is minutes, so (a)+(b)+(c)+(d) ≈ 50–100 GPU-h ≈ 6–13 node-hours in all.

### S5 Agent protocol on LCRC (owner: LCRC agent; the user pastes this into the Swing session's CLAUDE.md or the first prompt)

- Same rules as on the Mac: Python only; no git commits or pushes (the user pushes); never write under `database/` (the clone's copy is read-only input); outputs only under `.../prad_ml/output/ml_pileup_20261027/M2/`; Ta per run; bounded reports (`report.md` ≤ 150 lines, tables as TSV); no package installs beyond S2.1 without the user's word; no EVIO transfers initiated by the agent.
- PBS-specific: never run training, data generation or `prad2py` builds on `gpulogin1/2` (cgroup-throttled, no GPU); every compute step goes through `qsub` (interactive for ≤ 1 h debugging, batch otherwise); always `-A <project>`, `-q gpu`, explicit walltime ≤ 24 h, `--max-minutes` under it, checkpoint every epoch, `--resume`; stay under 4 running / 10 queued / 16 GPUs; request 1 GPU per job unless a measured need says otherwise; intermediates in `/lcrc/globalscratch/<user>/prad_ml` (purgeable, so nothing only there is a result), results on `/lcrc/project`; check `sbank-list-allocations -p <project>` before a sweep and stop when the remaining hours are < 2× the sweep's estimate; do not queue jobs that would run into the second-Monday maintenance window; `manifest.json` per run with input sha1s, git commit, PBS job id, `nvidia-smi` name, wall time.
- Reporting: one `report.md` per arm in the S1.3 format; numbers quoted with intervals; "confirm on login" items of this plan resolved in the first session and written to `M2/swing_facts.md` (actual module names, CUDA tag, `/scratch` size, array support, internet from compute nodes).

### S6 Verification and sync-back (owner: user with the local agent)

- [ ] `swing_facts.md` exists and every "confirm on login" item of §2 has a value.
- [ ] S2.1 CUDA test line and the `--smoke` train/eval on Swing reproduce the Mac smoke TSVs (CPU vs GPU: losses within 1e-4, same ROC AUC to 3 decimals).
- [ ] For each arm: `report.md`, `tables/`, `model/`, `manifest.json` present; sha1s of inputs match `MANIFEST.tsv`; acceptance row of S4 answered yes/no with the number.
- [ ] Sync-back rsync excludes checkpoints; `git status` on the Mac shows no change outside `output/ml_pileup_20261027/M2/` (nothing tracked touched); `database/` untouched on both sides (`git status database/` clean on Swing).
- [ ] Notes: a `notes/PbWO4_energy_resolution_<date>.md` section summarising M2 per arm with the Ta runs separate (local agent, after user review).
- [ ] Allocation: node-hours used per arm recorded from `sbank-list-allocations` in `M2/report.md`.

## 4. Files to create or modify

- Create: `plans/plan_ml_lcrc_20261028.md` (this file); `output/ml_lcrc_bundle_20261028/MANIFEST.tsv` (+ a `make_bundle.sh` of `rsync --relative` lines, local agent); `output/ml_pileup_20261027/M2/tools/m2_train.py`, `m2_eval.py`, `m2_data.py` (LCRC agent, derived from the M1 tools); `output/ml_pileup_20261027/M2/pbs/{interactive.txt,train_m2.pbs,sweep_m2.pbs}`; `output/ml_pileup_20261027/M2/swing_facts.md`; `output/ml_pileup_20261027/M2/<exp>/{report.md,manifest.json,...}`.
- Modify (local, Python, defaults byte-identical): `twopulse_realdata.py` (env-overridable `SCRATCH`, `promote`), `output/ml_pileup_20261027/M1/tools/m1_data.py` (env-overridable `SCR_M1`), `m1_train.py`, `m1_eval.py` (device selection, `--threads`).
- Not touched: `database/`, any C++ file, `mf_pulse.py`, `twopulse_fit.py`, the benchmark scripts, other plans and notes.

## 5. Interfaces and data changes

- Env var `PRAD_ML_SCRATCH` (default = current Mac scratchpad path) read by `twopulse_realdata.py` and `m1_data.py`.
- New result directory `output/ml_pileup_20261027/U5/decode/` (promoted U5 caches, 1.4 GB, with `manifest.json`).
- Bundle layout = mirrored `output/` relative paths; results layout `output/ml_pileup_20261027/M2/<exp>/` identical on Swing and Mac; `manifest.json` keys: `inputs{path: sha1}`, `git_commit`, `pbs_job_id`, `gpu`, `wall_s`, `seed`, `config`.
- No change to any TSV schema of U0/U5/M1; M2 tables reuse the M1 `cells.tsv`/`real_fake.tsv` columns.

## 6. Risks and compatibility

- Outbound internet from compute nodes unknown: pip/conda and CMake FetchContent may only work on login nodes; S2.1 and S2.3 give the fallback (install/configure on login, run/build in the job).
- The U5 caches live only in the session scratchpad: if the scratchpad is cleaned before S1.4 `promote`, they must be regenerated locally (`twopulse_realdata.py decode`, 63 s per the 20261028 note plus EVIO reads) — do S1.4 first.
- `/lcrc/globalscratch` may be purged without notice: nothing that is a result lives only there.
- Version drift: Swing torch 2.15.0+CUDA vs Mac 2.15.0 CPU; nondeterministic cuDNN kernels can change losses at 1e-5; the S6 smoke comparison sets the tolerance (1e-4), not bit-identity.
- Ta: the n200 bank (`templates_spline_ta_fix_20261028`) and the best-50 Ta bank differ; (c) must state which bank each number uses, as the 20261028 note §3 does.

## 7. Validation commands (the acceptance checklist is S6; the per-arm criteria are the S4 table)

- `shasum -a1 -c output/ml_lcrc_bundle_20261028/MANIFEST.tsv` passes on both sides; bundle ≤ 4 GB.
- Mac regression after S1.4: `OMP_NUM_THREADS=4 MPLBACKEND=Agg ~/.venvs/prad-ml-cpu/bin/python output/ml_pileup_20261027/M1/tools/m1_train.py --smoke <file> --out <scratch>` and `m1_eval.py --stage all --smoke` give TSVs identical to the pre-edit run (diff empty).
- Swing: `qsub -I ... ; python -c "import torch; print(torch.cuda.is_available())"` → `True`.

## Open Questions

1. Should the 9 uncommitted local changes be committed and pushed to `origin/waveform_jz` before Swing clones, or does Swing work from 287af16? (User's git operation either way.)
2. Transfer: Globus (install Globus Connect Personal on the Mac) or rsync over the ProxyJump for the 3.7 GB bundle?
3. Is building `prad2py` on Swing and shipping EVIO (53 GB for the 25 `pileup_rates` splits, 19 GB Ta-only) wanted at all, or is the 3.00 M-channel U5 cache the ceiling for real-waveform training? (Recommendation: npz only; revisit after arm (a).)
4. May the LCRC agent create the `prad-ml-gpu` conda env in `$HOME` on Swing with the pinned versions of S2.1 (the Mac rule "no installs without approval" is carried over)?
5. Order of arms given M1's outcome: (a) first if M1's real-WP fake is above expectation (the U5 pattern), else (b) first — confirm, or fix an order now.
6. Does the user want the Swing CLAUDE.md written by the local agent (a copy of S5) and shipped with the bundle, or typed by the user?

## Todo

- [ ] S0.1 project name and allocation checked (`sbank-list-allocations`)
- [ ] S0.2 ProxyJump login works; S0.3 Claude Code on Swing, outbound HTTPS confirmed
- [ ] S0.4 GitHub key/token on Swing; decision on pushing the uncommitted changes
- [ ] S0.5 transfer path chosen and a 100 MB test file verified; S0.6 landing dirs created
- [ ] S1.4 code edits: `PRAD_ML_SCRATCH`, `promote` of U5 caches to `output/ml_pileup_20261027/U5/decode/`, device/threads flags; Mac smoke TSVs unchanged
- [ ] S1.2 bundle built with `MANIFEST.tsv`, size ≤ 4 GB
- [ ] bundle transferred and `shasum -c` passes on Swing
- [ ] S2.1 conda env created and CUDA test passes in an interactive job; `swing_facts.md` started
- [ ] S2.2 clone + bundle in place; `import twopulse_fit, mf_pulse` works without `prad2py`
- [ ] S3 PBS templates written; `m1_train.py --smoke` runs on one A100 and matches the Mac smoke within 1e-4
- [ ] S4 (d) sweep if M1 misses the 80 % headline, else skip to the arm chosen in Open Question 5
- [ ] S4 (a) real-negatives / domain adaptation, report with the real-WP fake rate
- [ ] S4 (b) learned deconvolution and joint five-type regression at 10× M1 rows
- [ ] S4 (c) Ta 25800 / 26138 / 26230 per run, three tables
- [ ] S6 sync-back, manifests verified, node-hours recorded, notes section drafted for user review
- [ ] (later) `output/ml_cluster_*` V1/V1b and X17 three-cluster work through the same bundle mechanism; dumps sized first
