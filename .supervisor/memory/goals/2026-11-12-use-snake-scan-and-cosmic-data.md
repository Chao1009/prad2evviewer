---
domain: science
summary: Use snake-scan and cosmic data to improve HyCal energy resolution; remind the user at session start, with the three prerequisites they must supply.
affects: [analysis/pyscripts/notes/calibration_runs_snake_cosmic_20261112.md]
verified: 2026-11-12
source: user 2026-11-12; analysis/pyscripts/notes/plan_energy_resolution_improvement_20261112.md
supersedes:
---
Use the located snake-scan runs (PbWO4 23699-23720, lead glass 23721-23723, 2108.85 MeV bremsstrahlung) and cosmic runs (23741, 23743-23745, 24154, 24164, 24185, 24387-24402, 24422, 24617) to improve HyCal energy resolution; the user plans to work on it on 2026-11-13 and asked to be reminded at the start of the next session.
Why: the runs are stored at JLab only, none local, and runs up to 23849 need `database/history/daq_map_20260424_til_run023849.json`, so work cannot start until the user supplies:
1. One EVIO file of a scan run and one of a cosmic run (e.g. `prad_023701.evio.00000`, `prad_023744.evio.00000`), or JLab farm access.
2. The scan logs `calibration/logs/snake_scan_2026041{8,9}.log`, or the logbook entries for runs 23699-23723.
3. Whether HyCal HV changed after April 19, and whether the tagger was read out in the `prad_v0_equa` trigger.
