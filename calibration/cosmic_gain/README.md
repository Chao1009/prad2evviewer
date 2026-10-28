# Cosmic Gain

ROOT macros for iteratively tuning HyCal PMT high voltages from
cosmic-ray peak-height measurements.

## Files

| File | Role |
|---|---|
| `read_cosmic_json.C` | Loads a `cosmic_modules_run<N>.json` summary (per-module peak height, integral, σ, count from a cosmic run) and histograms the per-module means for both PWO (W1–W1156) and LG (G-IDs) channels. |
| `new_vset.C` | Reads the current `vset_iter<N>.json` together with the latest cosmic summary and produces the next `vset_iter<N+1>.json`.  V0Set is nudged by ±5/±10/±20 V per module to drive the peak-height mean into the **33–37 ADC** target band, clamped at **1270 V** (PWO) and **1800 V** (LG). |
| `cosmic_common.h` | Readers shared by the two macros: the per-module lines of a cosmic summary (first, i.e. oldest, run of each module) and the channel/V0Set lines of a vset file. |
| `vset_iter0.json` | Seed voltage set used at the start of the iteration. |

## Procedure

1. Take a cosmic-ray run with the current `vset_iter<N>.json` loaded.
2. `prad2ana_cosmic_test -r <run> -n <file_count>` (run from the directory holding `prad_023<run>.000NN_raw.root`, with `./fit_canvas/` present) produces `cosmic_modules_run<run>.json` (per-module peak summary).
3. `root -l read_cosmic_json.C` to inspect the per-module distributions and confirm the run has enough statistics.
4. Edit the input/output filenames at the top of `new_vset.C`, then `root -l new_vset.C` to write `vset_iter<N+1>.json`.
5. Load the new HV set on the next cosmic run and repeat until every module sits within the target band.
