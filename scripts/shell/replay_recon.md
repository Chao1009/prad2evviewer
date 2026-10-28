# replay_recon.sh

A one-shot shell script that automates the full replay pipeline for a single PRad-II run on JLab ifarm:
replay recon → filter → live charge → quick check.

`prad2ana_replay_recon` performs its own grouped `hadd` merge by default
(62 split ROOT files per merged recon ROOT).  The script then filters all
merged recon ROOTs in one invocation, runs live charge over all filtered ROOTs
together, and writes every product directly under one run output directory.

---

## Usage

```bash
cd /path/to/working/dir
chmod +x /path/to/prad2evviewer/scripts/shell/replay_recon.sh
bash /path/to/prad2evviewer/scripts/shell/replay_recon.sh
```

> **Important**: Run with `bash` directly. Do **not** use `source` — ifarm's default shell is `tcsh` and the bash-specific syntax (`${VAR:-default}`) will cause errors.

Run the script by its path inside the prad2evviewer checkout; do not copy it
elsewhere.  It sources `replay_common.sh` from the same directory, which holds
the defaults, prompts, jcache staging and pipeline shared with
`submit_replay_recon.sh`.  The working directory can be anywhere.

---

## Prerequisites

- ROOT environment is set up (`root` command available):
  ```bash
  source /path/to/root/bin/thisroot.sh
  ```
- `hadd` and `jcache` are available (standard on JLab ifarm)
- `prad2evviewer` has been compiled (`build/bin/` contains the executables)

---

## Interactive Parameters

The script prompts for the following inputs at startup. Press Enter to accept the default value shown in brackets.


| Prompt                                           | Description                                                                                 | Default                                          |
| ------------------------------------------------ | ------------------------------------------------------------------------------------------- | ------------------------------------------------ |
| `Enter run number`                               | 6-digit run number, e.g.`024650`                                                            | *(required)*                                     |
| `Enter replay mode`                              | `prad2`, `x17`, `x17_full`, `prad1`, or `random`                                            | `prad2`                                          |
| `Include GEM hit-level branches in output?`      | `yes` passes `-gem_hit` to replay                                                           | `no`                                             |
| `Enter prad2evviewer directory`                  | Source/installation directory                                                               | the checkout containing the script, or `$PRAD2_SOFT` |
| `Enter executable directory`                     | Directory containing the `prad2ana_*` executables                                           | `<PRAD2_SOFT>/build/bin`                         |
| `Enter EVIO cache base directory`                | Base directory for EVIO input data                                                          | `/cache/clas12/rg-o/data`                        |
| `Enter output base directory`                    | Root directory for all output files                                                         | `./`                                             |
| `Enter number of parallel jobs (-j)`             | Shared CPU count for`replay_recon -j`, `replay_filter -t`, and `quick_check -j`             | `15`                                             |
| `Enter GEM zero suppression (-z)`                | GEM zero-suppression sigma threshold                                                        | `5`                                              |
| `Enter max number of files to process (-f)`      | Maximum number of EVIO sub-files to process                                                 | `10000`                                          |
| `Enter replay merge group size (-m, 0 disables)` | Number of split recon ROOT files per merged output                                          | `62`                                             |
| `Cut JSON [default]`                             | Path to the cut config file for`replay_filter`; type `default` to use the built-in template | `<PRAD2_SOFT>/analysis/cuts/prad2_default.json`  |

All defaults can also be overridden via environment variables before running:

```bash
PRAD2_SOFT=/data/soft/prad2evviewer OUTPUT_BASE=/data/recon bash replay_recon.sh
```


| Variable               | Description                                                             |
| ---------------------- | ----------------------------------------------------------------------- |
| `PRAD2_SOFT`           | Root directory of the prad2evviewer source/installation (default: the checkout containing the script) |
| `PRAD2_BIN`            | Executable directory (default:`$PRAD2_SOFT/build/bin`)                  |
| `CACHE_BASE`           | Base directory for EVIO input data                                      |
| `MSS_BASE`             | Tape directory for jcache staging (default: `/mss/clas12/rg-o/data`)    |
| `OUTPUT_BASE`          | Base directory for all output files                                     |
| `REPLAY_CORES`         | Shared CPU count for replay, filter, and quick check                    |
| `REPLAY_ZERO_SUPPRESS` | GEM zero-suppression sigma threshold                                    |
| `REPLAY_MAX_FILES`     | Maximum number of EVIO sub-files to process                             |
| `REPLAY_MERGE_FILES`   | Number of split recon ROOT files per merged output;`0` disables merging |
| `DEFAULT_CUTS`         | Path to the default cut JSON file                                       |

If `PRAD2_BIN` or `DEFAULT_CUTS` are not set explicitly, they follow the
prad2evviewer directory entered at the prompt.

---

## Input Files


| Content           | Path                                                               |
| ----------------- | ------------------------------------------------------------------ |
| EVIO raw data     | `<CACHE_BASE>/prad_<RUN>/`, default `/cache/clas12/rg-o/data/prad_<RUN>/` |
| Cut configuration | User-specified, or `<PRAD2_SOFT>/analysis/cuts/prad2_default.json` |

---

## Output Files

All output files are written to `<OUTPUT_BASE>/prad_<RUN>/`:


| File                                                                        | Description                                                                  |
| --------------------------------------------------------------------------- | ---------------------------------------------------------------------------- |
| `prad_<RUN>.XXXXX_recon.root`                                               | Per-split recon ROOT output from replay                                      |
| `prad_<RUN>_recon_000.root`, `prad_<RUN>_recon_001.root`, ...               | Merged recon ROOT files, grouped by`REPLAY_MERGE_FILES`                      |
| `prad_<RUN>_filter_000.root`, `prad_<RUN>_filter_001.root`, ... | Filtered ROOT files, one per recon input                                     |
| `prad_<RUN>_epics.root`                                                     | Run-level slow-control ROOT containing only`scalers`, `epics`, and `runinfo`; written only when there is more than one recon input |
| `prad_<RUN>_filter_report.json`                                             | Run-level per-checkpoint filter verdict from`replay_filter`                  |
| `prad_<RUN>_live_charge.json`                                               | Live charge result in nC from all filtered ROOT files                        |
| `prad_<RUN>_quick_check.root`                                               | Quick-check histograms from all recon ROOT files                             |

---

## Pipeline Steps

```
0. Check ROOT environment
1. Collect parameters interactively
2. Check for regular files in the cache input directory
   └─ If missing or empty → submit jcache tape-staging request and exit
3. Create output directory
4. prad2ana_replay_recon    (multi-threaded replay; grouped merge via -m)
5. Select merged recon ROOTs, or per-split recon ROOTs when -m 0
6. prad2ana_replay_filter   (apply slow-control cuts to all recon ROOTs)
7. prad2ana_live_charge     (compute live charge from all filtered ROOTs)
8. prad2ana_quick_check     (fast quality-check histograms from recon ROOTs)
```

Each step prints its full command line.  The script stops at the first step
that fails and exits with that step's status.  If the filter executable or the
cut JSON is missing, the filter is skipped with a warning; live charge, which
reads the filter outputs, is skipped with it, and the quick check still runs.

---

## Data Not Yet in Cache (jcache Tape Staging)

If `<CACHE_BASE>/prad_<RUN>/` is missing or holds no regular files, the data is
likely still on tape.  The script will prompt for your email address and submit
a staging request:

```bash
jcache get <MSS_BASE>/prad_<RUN>/* -e <your@email>
```

`MSS_BASE` defaults to `/mss/clas12/rg-o/data`.  The script exits with status 1
if the email is empty or the jcache request fails.  You will receive an email
notification when the files have been moved to cache.  Re-run the script after
receiving the notification.

---

## Example Session

```
$ bash /path/to/prad2evviewer/scripts/shell/replay_recon.sh
ROOT is available: 6.32.02
Enter run number (e.g. 024650): 024650
Enter replay mode (prad2, x17, x17_full, prad1, or random) [prad2]:
Replay mode: PRad2
Include GEM hit-level branches in output? [no]:
Enter prad2evviewer directory [/path/to/prad2evviewer]:
Enter executable directory [/path/to/prad2evviewer/build/bin]:
Enter EVIO cache base directory [/cache/clas12/rg-o/data]:
Enter output base directory [./]: /home/liyuan/PRad2Analysis/data/recon
Enter number of parallel jobs (-j) [15]:
Enter GEM zero suppression (-z) [5]:
Enter max number of files to process (-f) [10000]:
Enter replay merge group size (-m, 0 disables) [62]:
Enter cut JSON file for replay_filter (...):
Cut JSON [default]:

Checking input directory: /cache/clas12/rg-o/data/prad_024650
Output directory: /home/liyuan/PRad2Analysis/data/recon/prad_024650

Starting replay...
Command: /path/to/prad2evviewer/build/bin/prad2ana_replay_recon ...
...
Replay finished.
Downstream input ROOT file(s): 2

Running replay filter...
...
Running live charge calculation...
...
Live charge JSON: .../prad_024650/prad_024650_live_charge.json

Running quick check...
...
Quick check ROOT: .../prad_024650/prad_024650_quick_check.root

Replay pipeline finished.
Recon inputs: 2
Output dir: /home/liyuan/PRad2Analysis/data/recon/prad_024650
```
