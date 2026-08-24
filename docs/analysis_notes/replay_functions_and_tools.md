# Replay Functions and Command-Line Tools

Primary sources:

- [`analysis/src/Replay.cpp`](../../analysis/src/Replay.cpp)
- [`analysis/include/Replay.h`](../../analysis/include/Replay.h)
- [`analysis/tools/replay_rawdata.cpp`](../../analysis/tools/replay_rawdata.cpp)
- [`analysis/tools/replay_recon.cpp`](../../analysis/tools/replay_recon.cpp)
- [`analysis/tools/replay_gainCorr.cpp`](../../analysis/tools/replay_gainCorr.cpp)
- [`prad2det/include/EventData_io.h`](../../prad2det/include/EventData_io.h)

Related notes:

- [Replay and slow-control quick start](replay_cuts.md)
- [HyCal gain correction](gain_correction.md)
- [HyCal time calibration](hycal_time_calib.md)
- [RF-time extraction](extract_rf_time.md)

This note is a compact reference for the four main `Replay` processing
functions and the three replay executables built on top of them.

---

## Quick Overview of Replay Procedure for Ordinary Analysis

The recommended reusable workflow is to decode EVIO once, preserve raw-replay
ROOT files, and run reconstruction from those files. Raw-to-recon avoids
repeating the expensive EVIO decode and is typically much faster than direct
EVIO reconstruction; the exact speedup depends on storage, compression, and
whether waveform samples are retained.

```text
                              EVIO files
                                  |
                                  v
                    prad2ana_replay_rawdata
                                  |
                                  v
                         *_raw.root files
                           |            |
                           |            +--> replay_filter
                           |                 + live_charge
                           |                 (run-quality checks)
                           v
                     prad2ana_replay_recon
                     (raw ROOT -> recon)
                                  |
                                  v
                        *_recon.root files

Alternative shortcut when reusable raw files are not needed:

EVIO files ----------------> prad2ana_replay_recon ----------------> *_recon.root
```

### Step 0: define paths and create output directories

```bash
RUN=024327
EVIO_DIR=/data/evio/prad_${RUN}
RAW_DIR=/data/replay_raw/prad_${RUN}
RECON_DIR=/data/replay_recon/prad_${RUN}
QC_DIR=/data/replay_qc/prad_${RUN}
NCPU=8

mkdir -p "$RAW_DIR" "$RECON_DIR" "$QC_DIR"
```

If the database is not in the installed default location:

```bash
export PRAD2_DATABASE_DIR=/path/to/prad2evviewer/database
```

### Step 1: EVIO to raw-replay ROOT

Recommended command when the raw files may be reconstructed again with changed
waveform, timing, or clustering settings:

```bash
prad2ana_replay_rawdata "$EVIO_DIR" \
    -o "$RAW_DIR" \
    -j "$NCPU"
```

This keeps waveform samples and GEM strip-level information. Use `-p` as well
when software/firmware peak branches are needed for raw-level diagnostics:

```bash
prad2ana_replay_rawdata "$EVIO_DIR" \
    -o "$RAW_DIR" \
    -j "$NCPU" \
    -p
```

For a smaller raw format when the waveform peak algorithm is already fixed,
use `--noWaveform`; it stores peaks instead of samples and remains usable by
the current raw-to-recon path:

```bash
prad2ana_replay_rawdata "$EVIO_DIR" \
    -o "$RAW_DIR" \
    -j "$NCPU" \
    --noWaveform
```

### Step 2: optional filtering and live-charge calculation

Run-quality filtering can be performed as soon as raw replay has produced the
`scalers`, `epics`, and `runinfo` side trees:

```bash
prad2ana_replay_filter "$RAW_DIR"/*_raw.root \
    -o "$RAW_DIR" \
    -c /path/to/prad2_default.json \
    -t "$NCPU"
```

Calculate live charge from the filtered outputs together:

```bash
prad2ana_live_charge "$QC_DIR"/*_filter.root \
    -j "$QC_DIR/prad_${RUN}_live_charge.json"
```

This branch of the workflow is useful for establishing good intervals and run
normalization before or alongside reconstruction. See
[replay_cuts.md](replay_cuts.md) for the cut-file format and report products.

**After filtering, the original `*_raw.root` files can be deleted to save disk space.**

`replay_recon` recognizes a raw replay input when its filename contains both
`_raw` and `.root`. This includes standard `*_raw.root` files and filtered
variants such as `*_raw_filter.root`, so either can be used in Step 3. A
standard `sample_raw.root` becomes `sample_recon.root`; a variant such as
`sample_raw_filter.root` becomes `sample_raw_filter_recon.root`.

### Step 3: raw-replay ROOT to reconstructed ROOT

Recommended reusable production path:

```bash
prad2ana_replay_recon "$RAW_DIR" \
    -o "$RECON_DIR" \
    -j "$NCPU"
```

`-m 0` keeps one reconstructed ROOT file per raw split. Omit it to merge
successful outputs in the default groups of 62:

```bash
prad2ana_replay_recon "$RAW_DIR" \
    -o "$RECON_DIR" \
    -j "$NCPU"
```

Add `-gem_hit` only when individual GEM-hit branches are needed; standard
HyCal–GEM match branches are written without it. Add `-X17` for X17 runs.

### Shortcut: EVIO directly to reconstructed ROOT

Use this when raw-level output will not be inspected or reused:

```bash
prad2ana_replay_recon "$EVIO_DIR" \
    -o "$RECON_DIR" \
    -j "$NCPU"
```

Both Step 1 and direct-EVIO reconstruction automatically create the run's gain
correction first if it is missing. For production, inspect that subprocess log
and verify that `<db>/gain_factor/gain_correction/prad_${RUN}_gain_corr.root`
was produced successfully.

---

## 1. Replay Processing Functions

### 1.1 Summary

| Function | Input | Main output tree | Decodes | Reconstructs | Applies corrections |
|---|---|---|---|---|---|
| `Replay::Process` | EVIO | `events` | FADC, GEM SSP strips, raw SSP/VTP/TDC snapshots, DSC2, EPICS, run control | Waveform peaks when requested; GEM strip processing | HyCal gain factor is looked up and stored, but not multiplied into waveform samples; GEM pedestal/common-mode/zero suppression |
| `Replay::ProcessWithRecon` | EVIO | `recon` | EVIO FADC, SSP, VTP clusters, RF TDC, DSC2, EPICS, run control | HyCal clusters, GEM hits, HyCal–GEM matching | Gain, energy calibration, module time offset/window, cluster linear correction, detector alignment/transforms, RF module offset, GEM pedestal/common-mode/zero suppression |
| `Replay::ProcessRaw2Recon` | `_raw.root` | `recon` | Re-decodes saved VTP/TDC words; reads saved FADC and GEM-strip data | Same HyCal/GEM/matching products as direct recon | Same reconstruction corrections as `ProcessWithRecon`; reuses saved GEM strips rather than decoding EVIO |
| `Replay::Process_LMSgainFactor` | EVIO | `lms_gain` | FADC waveforms for LMS/alpha events, DSC2, EPICS, run control | Soft waveform peaks only | No physics energy calibration; produces the measurements later used to calculate gain corrections |

All four paths also write or preserve the side trees `scalers`, `epics`, and
`runinfo` when those records are available.

### 1.2 `Replay::Process`: EVIO to raw replay ROOT

**Provides**

- Decoded event number, trigger information, and timestamp.
- One unified FADC channel array for Pb-glass, PbWO4, Veto, and LMS channels.
- Module ID/type lookup from `hycal_map.json`.
- Raw waveform samples by default.
- With peak output enabled: software waveform peaks, pedestal diagnostics, and
  emulated FADC250 firmware-mode peaks.
- Processed GEM strip-level hits after pedestal subtraction, common-mode
  treatment, zero suppression, mapping, and cross-talk tagging.
- Raw SSP trigger, VTP, and TDC bank snapshots for later decoding.
- A per-channel `hycal.gain_factor` selected by event number:
  PbWO4 uses the average of LMS channels 2 and 3, Pb-glass uses its stored
  average, and Veto/LMS use 1.

**Does not provide**

- No HyCal energy conversion or HyCal clustering.
- No GEM X/Y hit reconstruction.
- No HyCal–GEM matching.
- No detector-coordinate alignment or shower-depth projection.
- No module time-offset correction applied to the stored raw peak times.
- No RF-hit decoding or per-cluster `cl_dt_rf`.
- The gain factor is stored beside the raw data; waveform samples and peak
  integrals remain unscaled.

**Output**

- `events`: raw channel/waveform/optional peak and GEM-strip data.
- `scalers`: DSC2 scaler snapshots.
- `epics`: EPICS slow-control snapshots.
- `runinfo`: PRESTART/GO/END records and embedded DAQ configuration text.

Reference: `Replay::Process` and `prad2::SetRawWriteBranches`.

### 1.3 `Replay::ProcessWithRecon`: EVIO directly to reconstructed ROOT

**Provides**

- Direct EVIO decode followed by full reconstruction.
- Software waveform analysis for HyCal, Veto, and LMS channels.
- VTP `PRAD_CLUSTER` decode and RF TDC channel-A/channel-B decode.
- HyCal hit energies and clusters.
- GEM strip processing, X/Y clustering, Cartesian GEM hits, and HyCal–GEM
  matching for PRad-II.
- Optional individual GEM-hit branches; match branches are written regardless.
- PRad-II, PRad-1, X17 blind/full, and random-trigger paths.

**Corrections and transformations**

- Time-dependent HyCal gain correction is multiplied into pulse integral.
- HyCal energy calibration converts corrected ADC integral to MeV.
- Per-module HyCal time offset and run-dependent time window are applied.
- Seed-time coincidence, cluster energy correction, and cluster linear
  correction come from the reconstruction configuration/clusterer.
- Shower depth is assigned and HyCal clusters are transformed to the lab frame.
- GEM pedestal, common-mode, zero-suppression, mapping, reconstruction, and lab
  transform are applied.
- Per-cluster RF difference is folded, corrected by the center-module RF
  offset, and stored as `cl_dt_rf`.

**Does not provide**

- No raw HyCal waveform branches in the `recon` tree.
- The PRad-II reconstruction loop builds HyCal clusters from PbWO4 modules;
  Pb-glass handling is part of the separate PRad-1 path.
- Individual GEM-hit branches are omitted unless `-gem_hit` is passed.
- PRad-1 has no GEM reconstruction and no RF readout; `cl_dt_rf` remains NaN.
- This is reconstruction, not the derivation of gain, energy, time, alignment,
  or RF calibration constants.

**Output**

- `recon`: clusters, matches, trigger summaries, RF information, Veto/LMS
  summaries, and optionally GEM hits.
- `scalers`, `epics`, and `runinfo` side trees.

Reference: `Replay::ProcessWithRecon`, `PipelineBuilder`, `HyCalCluster`,
`GemSystem`, `MatchingTools`, and `SetReconWriteBranches`.

### 1.4 `Replay::ProcessRaw2Recon`: raw replay ROOT to reconstructed ROOT

**Provides**

- The same principal `recon` output as direct EVIO reconstruction.
- Reuses saved raw waveform samples when available; otherwise uses saved
  software peak branches.
- Re-decodes saved flat VTP and TDC bank snapshots.
- Reconstructs GEM hits from saved strip-level data, then performs matching.
- Clones `scalers`, `epics`, and `runinfo` from the raw file.

**Corrections**

- Loads current gain corrections again by event number. If a valid current
  correction is unavailable, it falls back to the raw file's stored
  `hycal.gain_factor`.
- Applies energy calibration, module time offset/window, clustering corrections,
  detector transforms, RF offset, and GEM reconstruction exactly at the
  raw-to-recon stage.

**Does not provide / restrictions**

- Does not support PRad-1 mode.
- `-z` cannot redo GEM pedestal/common-mode/zero suppression: raw input already
  contains strips after those operations.
- Requires the current GEM strip schema; legacy raw ROOT must be regenerated.
- If a raw file has neither waveform branches nor peak branches usable for
  reconstruction, HyCal reconstruction cannot recover the missing information.

Reference: `Replay::ProcessRaw2Recon`.

### 1.5 `Replay::Process_LMSgainFactor`: EVIO to LMS/alpha samples

**Provides**

- Selects events carrying LMS or alpha trigger bits.
- Runs software waveform analysis on mapped FADC channels.
- Classifies an LMS event as LMS-triggered with `nch > 1000`.
- Classifies an alpha event as alpha-triggered with `nch < 50`.
- Stores module IDs/types and software peak height, time, and integral in the
  `lms_gain` tree.
- Writes `scalers`, `epics`, and `runinfo` for batch time association.

**Does not provide**

- No HyCal energy calibration or cluster reconstruction.
- No GEM reconstruction or matching.
- Does not itself calculate the final correction ratios; that is phase 2 of
  `replay_gainCorr`, implemented by `ComputeGainCorrections`.
- The current phase-2 correction calculation produces the 1156 PbWO4 W-module
  correction arrays; it does not derive a new Pb-glass gain time series.

Reference: `Replay::Process_LMSgainFactor` and `GainCorrCompute.cpp`.

### 1.6 Supporting `Replay` methods

| Method | Purpose |
|---|---|
| `LoadDaqConfig` | Loads event tags, ADC/waveform settings, detector map paths, and DAQ decoder configuration. |
| `LoadHyCalMap` | Builds `(crate,slot,channel) -> name`, module type, module ID, and reverse-location maps. |
| `moduleName` / `moduleType` / `moduleID` | Convert DAQ addresses to detector identity. |
| `moduleLocation` | Reverse lookup from global module ID to crate/slot/channel. |
| `clear*Event` | Reset event buffers and default RF values before filling. |
| `setup*Branches` | Delegate the exact ROOT schema to `EventData_io.h`. |

---

## 2. Common Setup

Build the tools with ROOT analysis enabled:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_ANALYSIS=ON
cmake --build build --target replay_rawdata replay_recon replay_gainCorr -j4
```

Build-tree executables:

```text
build/bin/prad2ana_replay_rawdata
build/bin/prad2ana_replay_recon
build/bin/prad2ana_replay_gainCorr
```

For nonstandard layouts, set the database explicitly:

```bash
export PRAD2_DATABASE_DIR=/path/to/prad2evviewer/database
```

All three tools accept one or more input files/directories. Directory contents
are sorted. Raw replay and gain-correction tools collect EVIO files; recon also
accepts files ending in `_raw.root`.

The tools do not create the user-supplied `-o` directory. Create it before
running. `replay_gainCorr` separately creates its final database correction
directory when possible.

Both `replay_rawdata` and direct-EVIO `replay_recon` check for a run gain-file
under `<db>/gain_factor/gain_correction`. If missing, they automatically launch
`prad2ana_replay_gainCorr`. Failure is nonfatal: replay continues with the
available fallback/identity correction and prints a warning.

---

## 3. `prad2ana_replay_rawdata`

### 3.1 Purpose and syntax

Decode EVIO to per-file `_raw.root` files without physics cluster reconstruction.

```text
prad2ana_replay_rawdata <evio_file_or_dir> [more inputs...] -o <output_dir>
    [-f max_files] [-n max_events] [-j threads]
    [-c daq_config.json] [-d hycal_map.json]
    [-z zerosup_sigma] [-p] [-x] [--Ecalib] [--noWaveform]
```

### 3.2 Arguments

| Argument | Required? | Default | Meaning |
|---|---:|---|---|
| EVIO file/directory | Yes | — | One or more EVIO inputs. |
| `-o <dir>` | Yes | — | Output directory; it must already exist. |
| `-f <N>` | No | all | Process only the first `N` collected files. |
| `-n <N>` | No | all | Process at most `N` physics events per file. |
| `-j <N>` | No | 4 | Worker threads, capped by file count. |
| `-c <file>` | No | `<db>/daq_config.json` | DAQ/event-decoder configuration. |
| `-d <file>` | No | `<db>/hycal_map.json` | Unified detector/DAQ map. |
| `-z <sigma>` | No | 5 | Override GEM zero-suppression threshold. |
| `-p` | No | off | Add software peaks, pedestal diagnostics, and emulated firmware peaks. |
| `--noWaveform` | No | off | Omit waveform samples but force peak branches on. Useful for smaller diagnostic raw files. |
| `--Ecalib` | No | off | Keep only sum-trigger events, omit event number/timestamp, waveforms, GEM, SSP/VTP/TDC snapshots, and other extra branches; peak branches remain enabled. |
| `-x` | No | off | Select `reconstruction_config_x17.json` while building the pipeline. It does **not** reconstruct X17 clusters in the raw output. |

`-p`, `--noWaveform`, and `--Ecalib` differ as follows:

| Mode | Waveforms | Peak branches | GEM strips/raw trigger snapshots | Event filter |
|---|---:|---:|---:|---|
| default | yes | no | yes | all physics events |
| `-p` | yes | yes | yes | all physics events |
| `--noWaveform` | no | yes | yes | all physics events |
| `--Ecalib` | no | yes | no | sum-trigger only |

### 3.3 Examples

Standard raw replay with waveforms:

```bash
mkdir -p /data/replay_raw/prad_024327
prad2ana_replay_rawdata /data/evio/prad_024327 \
    -o /data/replay_raw/prad_024327 -j 8
```

Waveform and peak-algorithm diagnostics:

```bash
prad2ana_replay_rawdata prad_024327.evio.00000 \
    -o ./raw_check -n 100000 -p
```

Smaller raw output that keeps peaks but drops samples:

```bash
prad2ana_replay_rawdata /data/evio/prad_024327 \
    -o ./raw_peaks -j 8 --noWaveform
```

Compact energy-calibration sample:

```bash
prad2ana_replay_rawdata /data/evio/prad_024327 \
    -o ./raw_ecalib -j 8 --Ecalib
```

Test custom detector settings on a limited sample:

```bash
prad2ana_replay_rawdata /data/evio/prad_024327 \
    -o ./raw_test -f 2 -n 50000 -j 2 \
    -c ./daq_config_test.json -d ./hycal_map_test.json -z 6
```

---

## 4. `prad2ana_replay_recon`

### 4.1 Purpose and syntax

Reconstruct directly from EVIO or from current `_raw.root` files.

```text
prad2ana_replay_recon <evio_or_raw_file_or_dir> [more inputs...] -o <output_dir>
    [-f max_files] [-j threads] [-m merge_size]
    [-r reconstruction.json] [-c daq_config.json]
    [-d hycal_map.json] [-g gem_pedestal.json] [-z zerosup_sigma]
    [-prad1 | -x17 | -x17_full] [-random] [-gem_hit]
```

### 4.2 Arguments

| Argument | Required? | Default | Meaning |
|---|---:|---|---|
| EVIO/`_raw.root` input | Yes | — | Files, directories, or a mixture. |
| `-o <dir>` | Yes | — | Output directory; it must already exist. |
| `-f <N>` | No | all | Process only the first `N` collected files. |
| `-j <N>` | No | 4 | Reconstruction workers and maximum concurrent merge jobs. |
| `-m <N>` | No | 62 | Merge each group of `N` successful split outputs with `hadd`, then delete those inputs. `0` disables merging. A one-file remainder is left unmerged. |
| `-r <file>` | No | PRad-II or X17 database default | Reconstruction/cluster configuration. |
| `-c <file>` | No | `<db>/daq_config.json` | DAQ configuration. |
| `-d <file>` | No | `<db>/hycal_map.json` | HyCal/DAQ map. |
| `-g <file>` | No | run-config default | GEM pedestal file. |
| `-z <sigma>` | No | 5 | GEM zero-suppression override for direct EVIO only; unavailable for raw-to-recon. |
| `-gem_hit` | No | off | Write individual GEM-hit branches. Matching branches are present without it. |
| `-prad1` | No | off | PRad-1 map, ADC1881M pedestal, and no GEM/RF path. EVIO input only. |
| `-x17` | No | off | X17 blind mode; uses X17 config and retains only 3-cluster events whose event number ends in 8, while retaining raw-sum events. |
| `-x17_full` | No | off | X17 full mode without blind event-number thinning. |
| `-random` | No | off | Treat events as sum-triggered and widen the HyCal time window to 0–400 ns. |

Mode notes:

- With no mode flag, PRad-II mode is used.
- `-prad1` and either X17 mode are rejected together.
- `_raw.root` input does not support `-prad1`.
- `-random` with PRad-1/X17 prints a warning and should not be treated as a
  supported combined mode.
- The source header mentions `-n max_events`, but the current parser has no
  `-n` option and neither reconstruction function accepts a maximum-event
  argument. Use `-f`, a smaller input list, or first make a limited raw file.

### 4.3 Examples

Standard PRad-II reconstruction from EVIO, keeping per-split files:

```bash
mkdir -p /data/replay_recon/prad_024327
prad2ana_replay_recon /data/evio/prad_024327 \
    -o /data/replay_recon/prad_024327 -j 8 -m 0
```

Standard reconstruction with default merge groups of 62:

```bash
prad2ana_replay_recon /data/evio/prad_024327 \
    -o /data/replay_recon/prad_024327 -j 8
```

Reconstruct existing raw ROOT and include all individual GEM hits before matching:

```bash
prad2ana_replay_recon /data/replay_raw/prad_024327 \
    -o ./recon_from_raw -j 8 -gem_hit
```

Use custom reconstruction, map, pedestal, and zero suppression on EVIO:

```bash
prad2ana_replay_recon /data/evio/prad_024327 \
    -o ./recon_test -f 3 -j 3 \
    -r ./reconstruction_test.json \
    -c ./daq_config_test.json \
    -d ./hycal_map_test.json \
    -g ./gem_pedestal_test.json -z 6 -gem_hit
```

X17 blind and full samples:

```bash
prad2ana_replay_recon /data/evio/x17_run \
    -o ./x17_blind -x17 -j 8

prad2ana_replay_recon /data/evio/x17_run \
    -o ./x17_full -x17_full -j 8
```

Random-trigger reconstruction:

```bash
prad2ana_replay_recon /data/evio/random_run \
    -o ./random_recon -random -j 4
```

PRad-1 reconstruction:

```bash
prad2ana_replay_recon /data/evio/prad1_run \
    -o ./prad1_recon -prad1 -j 4
```

---

## 5. `prad2ana_replay_gainCorr`

### 5.1 Purpose and phases

The tool derives time-dependent HyCal gain corrections from EVIO data:

1. EVIO files are processed in parallel with `Process_LMSgainFactor`, producing
   per-file `_lms.root` trees.
2. All `lms_gain` trees are chained. Histograms are accumulated in batches of
   LMS events, fitted, and compared with the reference gain table.
3. Optionally, diagnostic plots are written.
4. Intermediate files are deleted by default or merged when `-s` is used.

Final output is always written to

```text
<db>/gain_factor/gain_correction/prad_<run:06d>_gain_corr.root
```

The `-o` directory is for intermediate LMS files, not the final correction
location.

### 5.2 Syntax and arguments

```text
prad2ana_replay_gainCorr <evio_file_or_dir> [more inputs...] -o <work_dir>
    [-f max_files] [-j threads]
    [-c daq_config.json] [-d hycal_map.json]
    [-b lms_events_per_batch] [-r reference_run]
    [-s] [-p] [-w id1,id2,...]
```

| Argument | Required? | Default | Meaning |
|---|---:|---|---|
| EVIO input | Yes | — | One or more EVIO files/directories. `_raw.root` is not accepted. |
| `-o <dir>` | Yes | — | Working directory for intermediate `_lms.root` files. |
| `-f <N>` | No | all | Process only the first `N` files. |
| `-j <N>` | No | 4 | Phase-1 EVIO replay workers. |
| `-c <file>` | No | `<db>/daq_config.json` | DAQ configuration. |
| `-d <file>` | No | `<db>/hycal_map.json` | HyCal/DAQ map. |
| `-b <N>` | No | 4000 | LMS events per gain-correction batch. Must be positive. |
| `-r <run>` | No | `gain_ref_run` from run config | Reference run for the baseline gain table. |
| `-s` | No | off | Preserve intermediate content by merging per-file LMS ROOT files into `<work_dir>/prad_<run>_lms.root`; otherwise delete them. |
| `-p` | No | off | Produce `<db>/gain_factor/gain_correction/prad_<run>_gain_corr_plots.pdf`. |
| `-w <list>` | No | none | Comma-separated W IDs to include in diagnostic plots; meaningful with `-p`. |

The correction ROOT `gain_corr` tree contains batch boundaries/counts,
reference-PMT ratios, current/reference W gains, correction factors, fitted
means, reference run, and batch time. See [gain_correction.md](gain_correction.md)
for formulas and branch-level interpretation.

### 5.3 Examples

Normal explicit gain-correction production:

```bash
mkdir -p ./lms_work
prad2ana_replay_gainCorr /data/evio/prad_024327 \
    -o ./lms_work -j 8
```

Use smaller time batches for finer time resolution:

```bash
prad2ana_replay_gainCorr /data/evio/prad_024327 \
    -o ./lms_work -j 8 -b 2000
```

Force a specific reference run:

```bash
prad2ana_replay_gainCorr /data/evio/prad_024327 \
    -o ./lms_work -j 8 -r 24185
```

Keep merged intermediate LMS data for debugging:

```bash
prad2ana_replay_gainCorr /data/evio/prad_024327 \
    -o ./lms_debug -j 8 -s
```

Produce diagnostic plots for selected crystals:

```bash
prad2ana_replay_gainCorr /data/evio/prad_024327 \
    -o ./lms_plots -j 8 -p -w 565,735,892
```

Quick limited-file configuration test:

```bash
prad2ana_replay_gainCorr /data/evio/prad_024327 \
    -o ./lms_test -f 2 -j 2 -b 500 \
    -c ./daq_config_test.json -d ./hycal_map_test.json -p
```

---

## 6. Which Tool Should Be Used?

| Goal | Recommended command |
|---|---|
| Preserve waveforms for detector/DAQ studies | `replay_rawdata` default |
| Compare software and firmware peak finding | `replay_rawdata -p` |
| Reduce raw output size but retain peaks | `replay_rawdata --noWaveform` |
| Produce a compact energy-calibration raw sample | `replay_rawdata --Ecalib` |
| Run normal physics reconstruction directly from EVIO | `replay_recon` |
| Change reconstruction without decoding EVIO again | `replay_recon <*_raw.root>` |
| Save individual GEM hits as well as matches | `replay_recon -gem_hit` |
| Derive/replace the run gain-correction time series | `replay_gainCorr` |

For ordinary physics production, direct `replay_recon` is sufficient; running
`replay_rawdata` first is necessary only when raw-level information must be
preserved or reconstruction will be repeated with different settings.
