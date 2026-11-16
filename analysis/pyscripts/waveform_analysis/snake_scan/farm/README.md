# P14.7 farm package: snake-scan constants for the whole scan (runs 23699-23723) and the cosmic run 23744

Everything here is run by the user on JLab ifarm (account jyzhou). No step needs a decision from me; each step prints
what it did. All paths below are relative to `<PRAD2_ROOT>/analysis/pyscripts/waveform_analysis/snake_scan/farm`
(`FARM`), on the Mac `/Users/L00338853/work/PRad/prad2evviewer/...`. Tests and measured numbers: `TEST_RESULTS.md`.

## What a job does

One job = one EVIO split, one core, writes one small file `<OUTDIR>/<run>/<file>.npz` and nothing else.

- Scan runs (23699-23723): ONE Python decode loop (`farm_pass.py`, scan mode) records, per 5000-event block, what the
  dwell classifier needs (`snake_scan._scan_blocks` rules, EPICS `ptrans`), and per event the 3x3 WaveAnalyzer integrals
  and peak counts around the running beam module (edge inputs for all pile-up settings), plus exact per-block counters.
  The LMS/alpha/pedestal series come from `lms_pass.lms_pass` (unchanged) in a second pass over the same file that costs
  2-5 s (the C++ decode is about 2 s per file; the cost of the first loop is the per-event Python channel loop).
  So: 2 decodes per file, 1 expensive. Wave 1 needed 3 passes (dwell, edge after the dwell table, LMS), each expensive but
  the LMS one.
- Run 23744: `cosmic_mip.scan_file` unchanged (vertical-track MIP hits, LMS, alpha events).
- Dwell selection, the edge fit (with the neighbour pass), the pooled LMS table, the transfer and the MIP table are made
  later by `farm_merge.py` from the products, with `snake_scan.blocks_to_rows`, `edge_pass.make_hist/fit_edge/second_pass`,
  `lms_pass.reduce_group`, `transfer.build`, `cosmic_mip.group_table` called unchanged (no existing module is modified).

Size and cost (measured, one core, Mac under load 7; details in TEST_RESULTS.md):

| item | scan file | 23744 file |
|---|---|---|
| wall time | 320-590 s | 92-95 s |
| max RSS | 280-676 MB | 460-540 MB |
| product, default `--keep-every 2` | up to about 10 MB (measured 9.8 MB; about 1.7 MB for a short last file) | 0.4-0.5 MB |
| product, `--keep-every 1` | 3.4-19.7 MB (the 20 MB target) | same |
| requested (swif2/sbatch) | 2 h, 2 GB RAM, 2 GB disk | 1 h, 1.5 GB RAM, 1 GB disk |

Job count: every 3rd split of 25 runs is about 1.02e9 events / 0.7 M per split / 3 = roughly 500 jobs (the brief said about 190:
that is the figure for every 9th split; the exact number comes out of step 2), plus all 23744 splits. At 6-10 min per scan job that is
50-90 CPU-h.

## Steps (in this order)

### 0. Put the code on ifarm (once)

The `snake_scan/` directory is untracked in git, so copy it (from the Mac checkout, directory `analysis/pyscripts`):

```
rsync -avP -e "ssh -J jyzhou@login.jlab.org" --exclude out --exclude __pycache__ --exclude '*.npz' \
    waveform_analysis/snake_scan jyzhou@ifarm:<PRAD2_ROOT>/analysis/pyscripts/waveform_analysis/
rsync -avP -e "ssh -J jyzhou@login.jlab.org" _common.py jyzhou@ifarm:<PRAD2_ROOT>/analysis/pyscripts/
```

`<PRAD2_ROOT>` is a checkout of prad2evviewer on ifarm (any branch with `database/`; `database/calibration/
calibration_factor_3p5_June7_nonLinear.json`, `database/runinfo/general.json` and `database/daq_config.json` are read). The
bindings must be built from the same `python/bind_dec.cpp` as on the Mac (the Mac working tree has an uncommitted change in
that file: copy it over or commit it first, otherwise step 1b fails or, worse, gives different integrals).

Settings (all optional; the defaults derive from the location of the scripts): `PRAD2_ROOT`, `PRAD2_DATABASE_DIR` (default
`$PRAD2_ROOT/database`), `PRAD2_BINDINGS` (directory with `prad2py*.so`, default `$PRAD2_ROOT/build/python`, as `setup_env.sh` does on the
Mac), `PYTHON` (needs numpy and scipy) and `PY_SETUP` (one line evaluated first, e.g. `source /work/.../venv/bin/activate` or
`module load python`). swif2 jobs inherit no environment, so write them once into `FARM/farm_env.sh` (template
`farm_env.sh.example`); `farm_job.sh` and `check_build.sh` source it first.

### 1. Build check

a. Bindings and packages (about 5 s):

```
cd <FARM>; ./check_build.sh
```
It must end with `OK` (it fails if the map copy differs from `snake_scan/`'s; python, numpy, scipy, `from prad2py import dec, det`, the study modules, the map copy equals
`snake_scan/hycal_map_til_run023849.json`). If `prad2py` is missing, build it: from `<PRAD2_ROOT>`
`cmake -S . -B build -DBUILD_PYTHON=ON -DPython_EXECUTABLE=$(which python3) && cmake --build build -j8` (README of the repo, RHEL 9
caveat about the interpreter), and set `PRAD2_BINDINGS` to the directory with the `.so`.

b. Same numbers as the Mac (needs `prad_023701.evio.00000` readable, i.e. after step 2 or from `/mss` via jcache; 10-20 s):

```
./check_build.sh /cache/clas12/rg-o/data/prad_023701/prad_023701.evio.00000
```
`check_build.sh` with an EVIO argument runs 1a and then `check_golden.py` on that file in the same environment (`farm_env.sh`,
`PY_SETUP`, `PRAD2_BINDINGS`), so it needs no manual setup. Expected last line `OK: identical` (counters, integral sums of the first 10000
events against `golden_23701_00000.json`; exit code 1 and `MISMATCH` otherwise).

### 2. File list and jcache

The file count per run is not known locally, so the list is built on ifarm from the tape stub listing:

```
./make_filelist.sh jcache_list.txt        # scan runs: splits NNNNN % 3 == 0 (00000, 00003, ...); 23744: all splits
./jcache.sh jcache_list.txt               # jcache get in batches of 100
./jcache.sh jcache_list.txt --check       # repeat until "missing: 0" (jcache get is asynchronous)
```
`make_filelist.sh` writes `jcache_list.txt` (the `/mss/clas12/rg-o/data/prad_0<run>/prad_0<run>.evio.<NNNNN>` paths) and
`jcache_list_cache.txt` (the same files under `/cache/clas12/rg-o/data/...`, the job input) and prints, per run, the splits on tape
and the selected ones. There is deliberately no shipped `jcache_list.txt`: its content depends on what is on tape. Runs 23721-23723 (lead
glass) are in the list as requested. Change the every-3rd rule with `STEP=` / `OFFSET=` (e.g. `STEP=9`) before running it. Local
naming check: splits are `prad_0<5-digit run>.evio.<5-digit split>` (e.g. `prad_023701.evio.00040`). jcache edits `/cache` by its
own rules (scratch with cleanup): submit soon after `--check` is clean.

Staging volume: each scan split is about 2.2 GB, so about 500 files x 2.2 GB is roughly 1.1 TB on `/cache` (plus the 23744 splits). If your
cache quota is smaller, stage and submit per run: make one list per run (`grep prad_0<run>. jcache_list.txt`), `jcache.sh` and `--check` that list,
submit it with the `submit_*` script, and continue with the next run when its products exist.

### 3. Submit

```
OUT=/volatile/<group>/jyzhou/p14_scan        # same form as the repo examples; /work also works; products about 3-5 GB, logs small
./submit_swif2.sh jcache_list_cache.txt $OUT            # writes $OUT/swif2_jobs.sh; review it
./submit_swif2.sh jcache_list_cache.txt $OUT --submit   # creates workflow p14_snake_scan, adds the jobs, runs it
```
or, with Slurm directly:
```
./submit_sbatch.sh jcache_list_cache.txt $OUT            # prints the two sbatch commands (scan array, 23744 array)
./submit_sbatch.sh jcache_list_cache.txt $OUT --submit
```
`submit_sbatch.sh` names its lists and command file with a timestamp (`slurm_scan_<tag>.list`, `slurm_cosmic_<tag>.list`,
`sbatch_cmds_<tag>.sh`; `TAG=` overrides), so a later resubmission of `missing.txt` does not overwrite the lists that running array tasks read.
Account and partition are the variables `ACCOUNT` and `PARTITION`, default `hallb` and `production` (the Slurm defaults of
`scripts/shell/replay_common.sh`); `submit_swif2.sh` passes them as `swif2 add-job -account -partition`. The repo has no other swif2 usage and I
could not run it on the farm: before `--submit` check the generated `swif2_jobs.sh` against `swif2 add-job -help` (option names, `-time/-ram/-disk` syntax) and
the account against `sacctmgr`/`sinfo`. Resource requests are the
variables `SCAN_TIME SCAN_RAM SCAN_DISK COS_TIME COS_RAM COS_DISK` (sbatch: `SCAN_MEM COS_MEM`). Smoke test first: put 2 lines of the
cache list in a file and submit that.

### 4. Monitor, rerun failures

```
swif2 status -workflow p14_snake_scan          # or: squeue -u jyzhou
find $OUT -name '*.npz' | wc -l                # products so far
./farm_missing.sh jcache_list_cache.txt $OUT > missing.txt     # lines without a product (also removes stale partial files)
./submit_swif2.sh missing.txt $OUT --submit    # (or submit_sbatch.sh) only the missing files
```
A job is idempotent: a finished file is skipped, a failed one leaves no product (the npz is written to a hidden `.partial.npz` and
renamed at the end), `./farm_job.sh <evio> $OUT --force` redoes one file by hand. Logs: `$OUT/logs/`. A file that is not on `/cache`
fails at once with the message "not on /cache yet".

### 5. Merge (on ifarm; light, a few minutes, under 1 GB RAM)

```
./farm_merge.sh --in $OUT --out $OUT/merged                          # all steps; the transfer step is skipped there
```
`farm_merge.sh` is `farm_merge.py` run in the environment of `farm_env.sh` (bindings path, `PY_SETUP`); it prints a warning when a product was made with
`--max-events` (partial) or products have different `VERSION`s. Without further options: wave-1 pile-up rejection (`--pileup w1`), pooled L_scan over all scan runs, edges from all scan runs except the
lead-glass runs 23721-23723. Output in `$OUT/merged` with the wave-1 names where one exists:
`dwell_<run>.tsv`, `dwell_scan.tsv`, `edge_scan.tsv` (columns as `edge_23701.tsv`) and `edge_hist/`, `lms_scan.tsv`,
`lms_scan_per_run.tsv` (+ `.npz`), `mip_23744.tsv`, `lms_23744.npz`. `edge_scan.tsv` replaces `edge_23701.tsv` because
`transfer.py` reads a fixed name: `farm_merge.py` stages a copy under that name for `transfer.build` and writes
`transfer_<rt>.tsv`, `prior_<rt>.tsv` with an explanatory header line.

Check `lms_scan_per_run.tsv` first: the pooling assumes the same PMT HV over the whole scan; a step (|z| > 3 in `L_rel_med`,
`R_LMS2/3` or `ped_med`, listed in the header) means the runs must be split (`--lms-runs RUN ...`, `--tag NAME` writes
`lms_scan_NAME.tsv`). `--lms-runs 23701` reproduces the wave-1 primary L_scan.

### 6. Bring the results back and run the transfer

```
rsync -avP -e "ssh -J jyzhou@login.jlab.org" jyzhou@ifarm:/volatile/<group>/jyzhou/p14_scan/merged ./
```
(about 5-20 MB). The transfer needs the production tables and the benchmark gains, which are on the Mac: on the Mac,
from `FARM`, either run step 5 with the products rsynced too (`rsync ... jyzhou@ifarm:/volatile/<group>/jyzhou/p14_scan ./p14_products`, 3-5 GB) or just
```
./farm_merge.sh --out ./merged --steps transfer --rt 2p24GeV_H2 3p49GeV_H2 [--lms-runs ... --tag ...]
```
(reads `./merged/edge_scan.tsv`, `lms_scan[_TAG].tsv`, `SS/out/lms_prod_<rt>.tsv`, V3 `ring/gains.tsv`).
Then: copy the tables you want to keep into `SS/out/` and tell me; P14.7 writes its results into `P14_snake_scan.md`.

## Switching the pile-up rejection

The products hold the peak count of every 3x3 member of every kept event, so the rejection is chosen at merge time, no re-decode:
`--pileup w1` (wave 1, default: no member with >= 2 peaks), `--pileup central` (only the central module must have exactly one
peak), `--pileup none`. Rerun step 5 with `--steps dwell,edge` (and `--tag` is not used by the edge table: use another `--out`).
`n_single`/`n_pileup` in the edge table are the exact counts for the chosen setting, for the events after the warm-up of the running module (see below). A recommendation that needs more than the peak
counts (for example peak times or amplitudes) needs a new job pass: add the quantity in `farm_pass.py`, bump `VERSION` and rerun.

## Assumptions and limits

- Beam energy 2108.85 MeV (`edge_pass.E_BEAM`) and the calibration file `calibration_factor_3p5_June7_nonLinear.json` (chain of
  `runinfo/general.json` for each run; the job prints it, the merge warns if the products disagree) are valid for the whole scan.
- The running beam module of a file is the mode of the argmax of the last 400 large-peak events (stable only, share >= 0.6). Inside a
  dwell it equals the wave-1 dwell module; the first 100 large-peak events of each file are not used and the module is stable only after the window (400 large-peak events) has settled, so the lag is up to about 240 large-peak events at the start of a file and of each dwell, and a
  module that is not the argmax of the large-peak events does not get edge events. The dwell classifier itself is the wave-1 code on the
  same block quantities, so dwell rows are identical to wave 1.
- Per-event arrays are stored only for blocks with a large-peak argmax share >= 0.6, for events with a 3x3 sum of 650-4100 MeV (the
  fit never uses less than about 700 MeV) and, by default, for every 2nd event (`--keep-every 2`; the statistical error of the
  edge grows by about 1.4 and stays below the window systematic). `n_single`/`n_pileup` are exact only for the events after the warm-up of the running module (the earlier ones have no module and are not counted).
- Modules with fewer than 20000 exact single events in their dwell(s) (`--min-n`) or a failed fit get no edge; the merge prints them.
- The pooled L_scan is a median over all scan-epoch LMS events (about 3 per file, so about 1500 over the scan) of the alpha-normalised
  single-peak integral (`lms_pass` conventions, wave-1 `MIN_ALPHA` rule); the per-run table gives each run's median ratio to the pool.
  Run 23744 (cosmic epoch) is shown in the per-run table (`epoch cosmic`, its own median ratio to the pool) and is not pooled. With `lms_pass.reduce_group`
  the median ratio 23744/23701 measured in the review of this package is 0.9935 (a 0.7 % offset); the wave-1 note quotes 0.975, which does not reproduce.
- Lead-glass runs 23721-23723: dwell rows are written (PbWO4 argmax; not meaningful for a lead-glass scan), their LMS events are pooled,
  their dwells are not used for edges (`--edge-runs` overrides).
- The C++ pipeline setup prints a long configuration block into every job log; harmless.
- The derived map is a copy (`FARM/hycal_map_til_run023849.json`, identical to `SS/hycal_map_til_run023849.json`; `check_build.sh`
  verifies it) so the farm package is self-contained.

## Files

`farm_job.sh` (job), `farm_pass.py` (job body), `farm_merge.py` + `farm_merge.sh` (wrapper), `make_filelist.sh`, `jcache.sh`, `submit_swif2.sh`, `submit_sbatch.sh`
+ `farm_array.sbatch`, `farm_missing.sh`, `check_build.sh`, `farm_env.sh.example`, `check_golden.py` + `golden_23701_00000.json`, `test_farm.py` (local
comparison with wave 1), `hycal_map_til_run023849.json`, `TEST_RESULTS.md`.
