# Scripts

Python utilities for detector visualization and monitoring.

## hycal_scaler_map.py

PyQt6 live colour-coded HyCal FADC scaler map. Polls `B_DET_HYCAL_FADC_<name>` EPICS channels every 10 s. Requires `pyepics` for real EPICS; `--sim` works without it.

```bash
python scripts/hycal_scaler_map.py          # real EPICS (default)
python scripts/hycal_scaler_map.py --sim    # simulation (random values)
```

## hycal_pedestal_monitor.py

PyQt6 GUI for measuring and monitoring FADC250 pedestals on all 7 HyCal crates. Reads pedestal means from `.cnf` files, parses per-channel RMS from `faV3peds` stdout, and flags irregular channels. No numpy/matplotlib -- only PyQt6.

```bash
python scripts/hycal_pedestal_monitor.py          # view existing data
python scripts/hycal_pedestal_monitor.py --sim     # test with simulated data
```

### Shift pedestal check (operator procedure)

Pedestals must be measured **before the first DAQ run of each shift** while DAQ is idle.

1. Make sure the DAQ is **stopped**.
2. Launch: `python scripts/hycal_pedestal_monitor.py`
3. Click **Measure Pedestals** and confirm. The tool SSHs to `adchycal1`--`adchycal7` and runs `faV3peds` (takes a few minutes).
4. Inspect the two maps (left: current mean, right: difference from configured) and the report panel for flagged channels:
   - `DEAD` -- avg < 1, rms < 0.1
   - `OUT OF RANGE` -- mean outside 50--300
   - `HIGH RMS` -- sigma > 1.5
   - `DRIFT` -- shifted > 3 counts from configured
5. Click **Save Report** to save for the shift log.
6. If new issues appear, notify the run coordinator before starting data taking.

Thresholds are defined at the top of the script and can be adjusted.

## trigger_mask_editor.py

PyQt6 visual editor for FAV3 trigger masks. Displays a HyCal geo view (with LMS1-3, LMSP, V1-V4 below) and lets you click or drag modules to toggle channels off/on. Generates trigger mask `.cnf` files -- only slots with disabled channels are written. Unmapped DAQ channels (slot positions with no module) are always masked off.

```bash
python scripts/trigger_mask_editor.py                     # start fresh
python scripts/trigger_mask_editor.py -i existing.cnf     # load existing mask
python scripts/trigger_mask_editor.py -o output.cnf       # set default save path
```

## gem_layout.py

Visualize GEM strip layout from `gem_map.json`.

```bash
python scripts/gem_layout.py [gem_map.json]
```

## gem_cluster_view.py

Visualize GEM clustering from `gem_dump -m evdump` JSON output.

```bash
python scripts/gem_cluster_view.py <event.json> [gem_map.json] [--det N] [-o file.png]
```

## tdc_viewer.py

PyQt6 viewer for V1190 TDC hits from the tagger crate (ROC `0x008E`, bank
`0xE107`). Shows a per-slot bar chart of hits/channel plus the TDC-value
histogram for any selected channel. No matplotlib; plots are drawn with
QPainter (numpy is the only scientific dep).

```bash
# Extract hits from an evio file (produces a 16-byte/hit binary)
./build/bin/tdc_dump /data/stage6/prad_023667/prad_023667.evio.00000 \
    -b /tmp/tagger_hits.bin           # add "-n 100000" to cap events

# Visualise
python scripts/tdc_viewer.py /tmp/tagger_hits.bin
```

`tdc_dump` also writes a CSV (`event_num,trigger_bits,roc_tag,slot,channel,edge,tdc`)
to stdout when no `-o` / `-b` is given, which is handy for `awk`/`head` checks.
