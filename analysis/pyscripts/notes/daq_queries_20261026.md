# DAQ questions: run 24386, veto V1-V4, FADC clock vs RF

Data basis: 13 runs, 111 splits. H2 at 0.73 GeV (24590, 24603), 2.24 GeV (25228, 25257, 25260, 25299, 25320) and 3.49 GeV (24386, 24503, 24512); Ta X17 at 2.24 GeV (25800, 26138, 26230). Per-split dumps of PbWO4 matched-filter pulse times, V1-V4 pulses and RF TDC hits are in output/module_timing_20261025/ (dump/, rf/). RF is CEBAF RF/32 in TDC ROC 0x40 slot 16 ch 0 and 8.

## 1. Run 24386 (3.49 GeV H2, 11 splits)

The FADC-RF phase drifts between splits: 0.88 ns span. The span is 1-4 ps in 9 of 13 runs; the other three are 24503 (0.23 ns), 24512 (0.06 ns) and 26230 (0.04 ns). Within most splits the events sit 50/50 in two RF phases 1.33 ns apart (a diagnostic estimate), with no link to the TI timestamp. The veto coincidence peak (t_veto - t_event) is at +18.49 ns, against -33.70 ns in 24503 and 24512 (0.73 GeV: -33.42, 2.24 GeV: -33.32). The run is excluded from the timing studies. In the per-run RF calibration, 78 % of the modules measured in 24386 show a bimodal RF phase (output/rf_calib_20261026/rf_calib_runs.tsv, column frac_bimodal); no other run exceeds 5.51 %, and the run is written as unusable in that calibration.

1.1 Was there a DAQ restart, clock resync, TI or FADC reconfiguration, cabling change or veto change during 24386?

1.2 Does the run log, e-log, run database or shift summary record anything for 24386, or flag it as bad or special?

## 2. Veto counters V1-V4

FADC crate 3 slot 3 ch 0-3, 50 x 50 mm; (x, y) in mm: V1 (-630, 0), V2 (-680, +50), V3 (-730, 0), V4 (-680, -50) (database/hycal_map.json). Read out in the H2 runs; no veto pulses in the Ta X17 data. No offsets for V1-V4 in database/hycal_time_offsets/time_calib.json. Fitted V1 - V4: -2.36 ns (0.73 GeV), -2.33 ns (2.24 GeV), -2.13 ns (3.49 GeV), i.e. V4 is 2.1-2.4 ns later than V1 (fit residual rms 0.07-0.11 ns). At 0.73 GeV V1 has two timing populations, at about -37.0 and -33.8 ns against the event time; origin not established. No coincidence peak in 2.24 GeV runs 25228, 25257, 25260 (present in 25299, 25320). Their veto pulses still show the accidental beam-bunch comb (period within 2.1 sigma of the RF value), so the counters were recording pulses.

2.1 What are the cable lengths or delays for each of V1-V4?

2.2 Was V1 retimed or re-cabled, or its discriminator/threshold or HV changed, during the 0.73 GeV runs? Does V1 have a split signal path?

2.3 What changed for the veto in 25228, 25257 and 25260: counters moved, HV, timing, or readout window?

2.4 Do time offsets for these four channels exist anywhere outside time_calib.json?

## 3. FADC clock vs RF

The accidental-coincidence comb in module-veto time differences gives a beam-bunch period of 4.00755 +/- 0.00013 ns (2.24 GeV, 46 teeth). This is consistent with 128.247/32 = 4.0077 ns (-1.3 sigma) and 59.5 sigma from 4.000 ns, so the FADC 250 MHz sampling clock is not RF-locked. Between the X17 runs 25800 and 26230 the FADC–RF phase jumps by +1.99 ns for the whole detector: 80 of 83 FADC boards agree with this value within ±0.1 ns. Three boards in crate 3 (slots 15, 8 and 18) deviate from it by +1.58, +0.32 and +0.14 ns; they read out 23 of the 33 modules of the central HyCal block, and in run 26230 the 26 central modules are early by 0.85–2.10 ns relative to the database time offsets, while in runs 25800 and 26138 they agree with the offsets, apart from the four modules named next. Four modules (W525, W528, W597, W697) are early in all three X17 runs: by 0.9–1.35 ns at 25800 and 26138 and by 1.1–2.0 ns at 26230.

3.1 Is the FADC clock documented as free-running relative to the RF?

3.2 Is there a description of the clock topology (source, fan-out, distribution per crate and per board)?

3.3 Is a common phase change of the FADC clock relative to the RF at each run start expected? We measure +1.99 ns between runs 25800 and 26230, the same for 80 of 83 boards.

3.4 Did anything change for FADC crate 3 slots 15, 8 and 18 (clock input, firmware, cabling, re-sync) at or before run 26230? These boards alone do not follow the common phase change, and the central block they read out is shifted in that run only.

Details: notes/PbWO4_energy_resolution_20261025.md §9.2 (veto), §9.3 (bunch period), §9.5 (RF phase, run 24386), §9.7 (RF time calibration per run), §10 (open items).
