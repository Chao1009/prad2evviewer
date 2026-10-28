# HyCal Calibration Tools

PRad-II, Jefferson Lab Hall B.

The directory groups three families of tools:

| Family | Entry point | Role |
|---|---|---|
| Gain equaliser | `hycal_gain_equalizer.py` | Closed-loop HV adjustment driving each module's edge ADC into a target band. |
| Snake scan | `hycal_snake_scan.py` | Position-driven scan along a predefined or auto-generated path; the operator-facing tool that parks the beam on each module long enough to measure it. |
| Cosmic-ray HV iteration | [`cosmic_gain/`](cosmic_gain/README.md) | Offline ROOT macros that propose a new `vset_iter{N+1}.json` from the latest cosmic peak summary. |

Each scan tool runs in three modes selected by command-line flag:
**expert** (full control over real EPICS), **observer** (real EPICS
reads, no writes), and **simulation** (fake motor / scaler so the
GUI can be exercised offline).

---

## Gain Equalizer — Operator Manual

### Quick start

On a counting-house machine, logged in as **clasrun**:

```bash
cd /home/clasrun/prad2_daq/prad2evviewer
python calibration/hycal_gain_equalizer.py --expert
```

The defaults (server, HV address, password, target ADC, minimum counts,
…) are designed to work out of the box.  **Discuss with the run
coordinator (RC) before changing them.**

### Running a gain scan

1. Select **Path**: `snake-all-pwo-r`.
2. Set **Start** module.  The plan is to finish the bottom part first, then return to the top.  Ask the RC if you do not know which module to start from.
3. Click **Start**.

### During the scan

Pay attention to:

- **Red marker** (expected beam position) vs. **hot region** on the scaler map — they should be at the same location.
- **Histogram** in the right panel — a clear spectrum should be building up.
- **Event log** — watch for `WARN` and `ERROR` messages.

### On `ERROR` (scan stops)

An `ERROR` stops the scan automatically.  The **Start** button changes
to **Resume** and the **Stop** button changes to **Reset**.

1. Note the module name and error message from the log.
2. Click **Resume** to retry the failed module.
3. To start over from a different module, click **Reset** to clear the scan, pick a new starting module, and click **Start**.
4. If the same error repeats, **contact the RC.**

### After each row

1. Find the screenshots in `/home/clasrun/prad2_daq/prad2evviewer/calibration/logs/`:

   ```
   GE_20260406_143025_W100_success.png
   GE_20260406_143530_W101_failure.png
   ```

   Files are named `GE_{timestamp}_{module}_{status}.png` and sort by time.

2. Post a log entry to **PRADLOG** with all screenshots from the completed row.
3. Failed modules need a redo or manual gain equalisation — discuss with the RC if uncertain.

### Event log files

Full event logs are saved to `calibration/logs/gain_eq_YYYYMMDD.log`
(one file per day, appended to across sessions).  Upload these with
the PRADLOG entry.

---

## Snake Scan — Operator Manual

### Quick start

```bash
cd /home/clasrun/prad2_daq/prad2evviewer
python calibration/hycal_snake_scan.py --expert
```

### Running a scan

1. Select **Path** profile (`(autogen)` or a predefined entry from `paths.json`).
2. For autogen: set **LG layers** (0 = PbWO4 only; 1–2 to include PbGlass).
3. Talk to the RC if you are unsure about the scan path.
4. Set **Start** module and **Count** (0 = scan all from start to end).
5. Set **Dwell time** and **Position threshold**.
6. Click **Start Scan**.

### During the scan

Verify that:

1. The **red marker** position matches the **scaler hot spot** — the current module should have the highest scaler reading.
2. **WARN** / **ERROR** messages in the event log are addressed promptly.

### Resume after interruption

1. Find the last completed module (the **Done** colour on the map, or check the event log).
2. Select the next module as **Start**.
3. Click **Start Scan**.

---

## Tools

| Script | Purpose |
|---|---|
| `hycal_gain_equalizer.py` | Closed-loop gain equalisation GUI (expert / observer / simulation). |
| `hycal_snake_scan.py` | Snake-scan GUI with dwell control (expert / observer / simulation). |
| `scan_path_editor.py` | Manual path-builder GUI. |
| `gain_scanner.py` | Gain-scan engine, spectrum analyser, and HTTP / HV clients. |
| `scan_geoview.py` | HyCal map widget with scaler overlay; subclasses `scripts/hycal_geoview.py`'s `HyCalMapWidget`, so `scripts/` must sit next to `calibration/` (as in the source checkout). |
| `scan_epics.py` | EPICS PV utilities (motor, scaler). |
| `scan_engine.py` | Scan-path engine and motion executor. |
| `scan_utils.py` | Shared types, constants, coordinate transforms, path-profile loading, and theme. |
| `scan_gui_common.py` | Shared GUI code for both scan GUIs: the `ScanWindowBase` main-window shell (top bar, map with scaler controls, path profiles, event log, beam / encoder displays), the `run_scan_gui` command-line entry point, focus-guarded input widgets, the position-check panel, and EPICS bring-up. |
| `pmt_response.py` | Power-law PMT gain model (`edge = A · V^k`); proposes ΔV for the next iteration when the fit is trustworthy and falls back to a static lookup otherwise. Pure Python, no third-party dependencies — testable in isolation. |
| `paths.json` | Predefined scan-path profiles. |

### Command-line modes

```bash
python calibration/hycal_gain_equalizer.py             # simulation (read-only)
python calibration/hycal_gain_equalizer.py --expert    # expert (full control)
python calibration/hycal_gain_equalizer.py --observer  # observer (real reads, no writes)

python calibration/hycal_snake_scan.py                 # simulation
python calibration/hycal_snake_scan.py --expert        # expert
python calibration/hycal_snake_scan.py --observer      # observer

python calibration/scan_path_editor.py                 # path editor
```

## Coordinate system

| Beam at HyCal (0,0) | `ptrans_x = -126.75` | `ptrans_y = 10.11` |
|---|---|---|

```
ptrans_x = -126.75 + module_x       (x same direction)
ptrans_y =   10.11 - module_y       (y inverted)
```

Travel limits: `ptrans_x` ∈ [−582.65, 329.15] mm,
`ptrans_y` ∈ [−672.50, 692.72] mm.

## Safety

- **Expert mode** writes to four motor PVs only: `ptrans_{x,y}.VAL` and `ptrans_{x,y}.SPMG`.
- **Gain equaliser** additionally sends HTTP commands to `prad2hvd` (HV set) and `prad2_server` (histogram clear).  HV limits are enforced server-side by `prad2hvd`.
- **Observer mode** uses real EPICS reads but blocks all writes.
- **Simulation mode** uses a fake motor (no real EPICS writes) and blocks all HV / server writes.
