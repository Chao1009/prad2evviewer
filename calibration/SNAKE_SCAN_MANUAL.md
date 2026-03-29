# HyCal Snake Scan -- Operator Manual

**PRad-II, Jefferson Lab Hall B**

---

## Overview

Automates beam calibration of HyCal modules in a serpentine pattern.
All 1152 PbWO4 modules are scanned first, then surrounding PbGlass
sectors (0--2 layers, configurable). At each module the transporter
dwells for a configurable time (default 120 s), then advances.

The scan path is: **Center (PbWO4)** → **Bottom** → nearest **Side** →
**Top** → other **Side**, minimising slow y-axis travel.

A red dashed boundary on the map shows the transporter travel limits.
Moves outside this boundary are blocked automatically.

**Safety:** Only four EPICS PVs are written: `ptrans_{x,y}.VAL` and
`ptrans_{x,y}.SPMG`. All others are read-only.

---

## Launching

On **clonpc19**:

```bash
cd ~/prad2_daq/prad2evviewer && ./calibration/hycal_snake_scan.py --real
```

For simulation (no EPICS, no hardware):

```bash
python calibration/hycal_snake_scan.py
```

---

## Coordinate System

| Beam at HyCal (0,0) | ptrans_x = **-126.75** | ptrans_y = **10.11** |
|---|---|---|

```
ptrans_x = -126.75 + module_x       (x same direction)
ptrans_y =   10.11 - module_y       (y inverted)
```

Travel limits (symmetric about centre):
ptrans_x: **-582.65** to **329.15** mm,
ptrans_y: **-672.50** to **692.72** mm.

---

## Controls

| Control | Description |
|---------|-------------|
| **LG layers (0-2)** | PbGlass layers to include (0 = PbWO4 only). Locked during scan. |
| **Start** | Starting module (dropdown or click map) |
| **Count** | Number of modules to scan (0 = all from start to end) |
| **Dwell time** | Seconds per module (default 120) |
| **Pos. threshold** | Max position error in mm (default 0.5) |
| **Start Scan** | Begin scan |
| **Pause / Resume** | Pause/resume motors and dwell countdown |
| **Stop** | Abort scan, stop motors |
| **Skip Module** | Skip current dwell, advance to next module |
| **Ack Error** | Acknowledge position error, continue |
| **Move to Starting Point** | Move beam to selected module without scanning |
| **Reset to Beam Center** | Return to ptrans(-126.75, 10.11) |

---

## Running a Scan

1. Set **LG layers** (0 for PbWO4 only, 1--2 to include PbGlass).
2. Set **Start** module and **Count** (0 = scan all).
3. Set **Dwell time** and **Pos. threshold**.
4. Click **Start Scan**.

### Resume After Interruption

1. Find the last completed module (blue on map or in event log).
2. Select the next module as **Start**.
3. Click **Start Scan**.

### Position Error

Scan auto-pauses. Click **Ack Error** to skip and continue, or **Stop** to abort.

---

## Module Map

- Click a module to select it as start and see its info in the frame title.
- Click again to deselect.
- Path preview (blue line) shows the planned route from start to end.
- Red crosshair tracks actual motor position.
- Red dashed rectangle shows transporter travel limits.

| Color | Meaning |
|-------|---------|
| Dark grey (Todo) | Pending | Yellow | Moving to |
| Green | Dwelling | Blue | Done |
| Red | Position error | Orange | Selected start |
| Dim | Skipped (before start / after count) |

---

## Logging

All events are logged to `calibration/logs/snake_scan_YYYYMMDD_HHMMSS.log`
(one file per session, created on launch).

**Upload log files when the scan or shift ends.**

---

## Shift Checklist

During the scan, verify that:

1. The **red crosshair** position on this GUI matches the **occupancy
   plot** from the online event monitor — the current module should show
   the highest occupancy.
2. The **FADC scalers** are consistent with the beam being on the
   expected module.

If the occupancy does not match, **Stop** the scan, set **Start** to the
last good module, and try again. Log the observation and action on
**PRADLOG**. If the error persists, **contact the run coordinator**.

---

## Troubleshooting

If a problem occurs, **try again first**. If it persists, contact the
run coordinator.

| Symptom | Likely cause |
|---------|--------------|
| PVs not connecting | IOC down, network / firewall |
| Motors don't move | Interlocks, motor enable, SPMG not Go |
| Move blocked | Target outside travel limits (see log) |
| Frequent position errors | Motor speed, backlash, encoder |
| Move timeout (>300 s) | Motor stall, limit switch, IOC |
