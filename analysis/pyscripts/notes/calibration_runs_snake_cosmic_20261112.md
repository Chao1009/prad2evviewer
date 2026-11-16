# PRad-II snake-scan and cosmic runs: inventory and calibration use (P13)

Date: 2026-11-12. Plan: `notes/plan_energy_resolution_improvement_20261112.md`, package P13 (U4). Read-only
research; no data were decoded.

Scope: find the PRad-II HyCal snake-scan and cosmic-ray runs (run numbers, dates, beam and trigger
configuration, event counts, data location, readability from this machine), state what each set supports
for U4, and what is missing to use it. Out of scope: running anything, code changes, other calibration
runs except where they bound the answer.

## 1. Answer

The runs exist and are identified. RCDB (public, no login) lists them with operator comments. None is on
this machine; the data are at JLab. P14 can start its research stage once one EVIO split of each set is
copied here or read on the JLab farm.

| Set | Runs | Date (2026) | Beam / target | Trigger config | Events | Notes |
|---|---|---|---|---|---|---|
| Gain equalisation (HV loop) | 23638-23661, 23689-23698 | Apr 15-18 (equaliser timestamps; RCDB rows from Apr 16) | 2108.85 MeV, "Harp Wire", 3 nA request | `prad_v0_equa` | 23638-61: 2.0e7, most runs 0; 23689-98: 3.1e8 | HV changed during the runs; the equaliser read online histograms. Not a calibration data set |
| **Snake scan, PbWO4** | **23699-23720** (22 runs; 23704, 23707, 23714 short or "bad spectrum") | Apr 18 15:49 - Apr 19 10:03 | 2108.85 MeV bremsstrahlung, "Harp Wire", 3 nA request | `prad_v0_equa` (run type `PRAD_VTP`) | 1.02e9 | comments "Scan of HyCal crystals", "HyCal scan W419-421", "Scan paused at W477" |
| **Snake scan, lead glass** | **23721-23723** | Apr 19 10:04-18:25 | same | `prad_v0_equa` | 1.38e8 | "Pb-Glass layer", "Lead glass ring" |
| Cosmic, HV iteration | 23555, 23561, 23570, 23575 | about Apr 7-9 | no beam | unknown | unknown | not in RCDB (RCDB jumps from 23065 to 23585); known from the repo only. Taken before equalisation, so at other HVs |
| **Cosmic, after the scan** | **23741, 23743, 23744, 23745** | Apr 20 13:45 - Apr 21 09:33 | empty target | `prad_v0_cosmic` | 2.06e7 (23744: 1.40e7 in 12.4 h) | closest in time to the scan, presumably same HV (unverified) |
| Cosmic, 3.49 GeV period | 24154, 24164, 24185 (`prad_v1_cosmic`); 24387, 24389, 24390, 24393, 24394, 24396, 24400-24402, 24422 (`prad_v1`) | Apr 30 - May 7 | empty / H2 cell, no beam | as listed | 24164: 1.67e7 (14.5 h); 24387-24402 set: about 7e7 | 24387 "only cosmics" with H2; 24390 "had also beam on tagger"; 24155/24156 failed (no sum trigger / only LMS) |
| Cosmic, 0.73 GeV period | 24617 | May 15 | empty | `prad_v1_cosmic` | 3.6e6 (4.4 h) | "COSMIC"; 24615/24616 DAQ wrong |
| Cosmic, X17 period | 25461 | Jun 27 | empty | `prad_v1_cosmic` | 9.3e3 | negligible |

Many further `prad_v1_cosmic` runs (24091-24172) are GEM tests; most are "trash, no rate" or have
unverified HyCal content.

## 2. Evidence (verified)

- RCDB run list, public: `https://clasweb.jlab.org/rcdb/runs/23585-24200` and `.../runs/24201-26500`
  (fetched 2026-11-12). Columns used: start time, duration, events, current, energy, target, trigger,
  comment. All rows in the table above come from these pages, except the four HV-iteration cosmics.
- RCDB run page `https://clasweb.jlab.org/rcdb/runs/info/23701`: 54,305,662 events, 81 EVIO files,
  165,323 MB, rate 16.96 kHz, `run_config /usr/clas12/release/2.0.0/parms/trigger/prad_v0_equa.cnf`,
  `run_type PRAD_VTP`, beam_current readback 0.046 (request 3), comment "Scan of HyCal crystals".
- RCDB run page `.../runs/info/23744`: 14,039,373 events, 76 files, 154,027 MB, 291.8 Hz,
  `prad_v0_cosmic.cnf`, target Empty, 12:21:18.
- RCDB `.../runs/info/23575`: "Run 23575 is not found ... Previous 23065, Next 23585".
- HV-iteration cosmics: `calibration/cosmic_gain/new_vset.C:8` (`cosmic_modules_run575.json`),
  `read_cosmic_json.C:5` (`run570`); git history of `read_cosmic_json.C` names `run555`, that of `new_vset.C` `run561`;
  `analysis/README.md:285-286` ("`-r 575` reads run 23575"). Seed HV `vset_iter0.json` is dated
  2026-04-07 14:55.
- Gain equaliser output: `database/gain_equalization_results.json`, 1144 PbWO4 modules, timestamps
  2026-04-15 05:00 to 2026-04-17 08:14, target edge ADC 1920 (1090 modules) or 1910 (54). The equaliser drives
  "the right edge of the Bremsstrahlung spectrum" to a target ADC (`calibration/gain_scanner.py:5-6`,
  `hycal_gain_equalizer.py:7-8`).
- Snake-scan tool: `calibration/hycal_snake_scan.py`; manual `calibration/docs/hycal_snake_scan_wiki.txt`
  (beam centred on each module, configurable dwell, log `calibration/logs/snake_scan_YYYYMMDD.log` on the
  counting-house account). Path profiles in `calibration/paths.json` (`snake-all-pwo`, `snake-all-pwo-r`,
  lead-glass layers).
- Derived cosmic products already in the repo: `database/gain_response_curve.json` (per-module
  gain(HV) fit from cosmics for 1134 of 1152 PbWO4 modules, plus LMS);
  `database/calibration/old/adc_to_mev_factors_cosmic.json` (added 2026-04-12, "a rough calibration factor
  file from cosmic data", base_energy 25 MeV, median factor 0.108 MeV/ADC).
- Tagger: E-counter `tage1` added to the readout in runs 23665-23667 (Apr 17, MOR trigger, 0.5 um Ni);
  tagger-HyCal correlation macros use runs 23667, 23671, 23686 (`analysis/scripts/tagger_hycal_correlation.C:31`,
  `scripts/tagger_viewer.py:29`); ROC 0x8E "tagger" in `database/daq_config.json:51`.
- Data location convention: `/mss/clas12/rg-o/data` (tape) and `/cache/clas12/rg-o/data/prad_<RUN>/`
  (`scripts/shell/replay_recon.md:73,91`); counting-house stage disk `/data/stage6/prad_0<RUN>/`
  (`analysis/scripts/tagger_hycal_correlation.C:31`).
- Local data: `/Users/L00338853/work/PRad/data/evio` holds splits of 18 production runs (24386-26230)
  only; no run below 24386. Free disk: 87 GiB.
- DAQ map: runs up to 23849 were taken with `database/history/daq_map_20260424_til_run023849.json`;
  the decoder and replay contain no run-based selection of the history maps (grep for `history/` and
  `til_run` in `prad2dec`, `prad2det`, `analysis`, `python`, `src` finds nothing). LMS1-3 sit in other
  slots in that map (`docs/research-lms-gain-tracking.md:129`).
- `database/runinfo/general.json`: the `from_run: 0` configuration (beam 2108.8 MeV) covers all April
  runs; production starts at 24185.
- Public run plan (Gasparian, indico.jlab.org/event/989, slides 6-7) planned gain equalising with a
  photon beam and a separate calibration scan with the MOR trigger, "three ranges of T- and E-counters
  (for linearity check)", ~1 min per module, HyCal X,Y stored from EPICS. The Hall B status talk of
  2026-06-29 (indico.jlab.org/event/1077, slide 12) calls the actual calibration "improvised".

## 3. Inferences (not verified)

- The scan used an untagged bremsstrahlung beam with endpoint 2108.85 MeV, as the equaliser did. Whether
  the tagger TDC (ROC 0x8E) was in the `prad_v0_equa` readout is unknown; RCDB does not say.
- HyCal X,Y per dwell come from EPICS (the scan GUI reads ptrans); whether EPICS events with ptrans are in
  the EVIO stream is unverified. Without them the module-to-time mapping needs the scan log.
- HV between the scan (Apr 18-19) and production (from May 2) is not documented in the repo. If HVs were
  reset, scan gains transfer to production only through LMS ratios, which need LMS events in both.

## 4. What each set supports (U4) and what is missing

| U4 use | Snake scan 23699-23723 | Cosmics (23741-45; 24164, 24387-402, 24617) |
|---|---|---|
| Absolute per-module gain independent of ep kinematics | Yes, from the bremsstrahlung edge at 2108.85 MeV per module, all PbWO4 plus lead glass. Edge position depends on the resolution and threshold model, so it is an edge fit, not a peak | No (MIP scale, ~25 MeV) |
| Outer PbWO4 modules poorly covered by ep and Møller | Intended to cover every module (~1e9 events); actual coverage needs the scan log | Relative only |
| Per-module non-linearity | Only if tagged energies exist (tagger in readout) or the spectrum shape is fitted at several energies; one endpoint alone gives no NL | Weak: MIP to GeV extrapolation |
| Position dependence across a module | Only with per-event beam position (EPICS ptrans plus beam-spot size, or GEMs; GEM presence during the scan unverified) | No |
| Fourth energy point for the three-term fit | Weak: 2.11 GeV is next to 2.24 GeV; a tagged sub-range would add points | No |
| Relative gain equalisation and stability | Yes, at one time | Yes: compare Apr 20-21, Apr 30-May 7 and May 15 per module (time drift across the 3.49 and 0.73 GeV periods) |

Missing to use the snake scan: (1) the EVIO data (one split, ~2 GB, to test; full PbWO4 scan ~3 TB, which
does not fit on this machine, so production processing belongs on the JLab farm); (2) the
`daq_map_..._til_run023849.json` map loaded for these runs (code has no switch); (3) the module-to-time
mapping (EPICS in the stream, or `snake_scan_20260418/19.log` and PRADLOG screenshots from the counting
house); (4) the HV history from Apr 19 to production; (5) whether the tagger was read out.

Missing to use the cosmics: the EVIO data (23744 is 154 GB); the history DAQ map for 23741-23745; the
existing `prad2ana_cosmic_test` (`analysis/tools/cosmic_test.cpp`) expects `prad_023<run>.000NN_raw.root`
replay files, so a replay step first.

## 5. Risks, open questions, coverage gaps

- RCDB comments are operator text; event counts of 0 in RCDB (runs 23639-23658) may mean the run end was
  not recorded, not that no data exist.
- RCDB has no entries for 23066-23584; the April 7-9 cosmics and any earlier scan tests are not listed.
- PRADLOG (logbooks.jlab.org) needs a login (HTTP 401); not read. The run-plan PDF and status talks are
  public; the "HyCal_calibration_proposal.pdf" was not fetched.
- The trigger `.cnf` files under `/usr/clas12/release/2.0.0/parms/trigger/` are not public; thresholds and
  prescales for the scan and cosmic triggers are unknown here.

## 6. What the user needs to provide

1. One EVIO split each of a scan run and a cosmic run, e.g. `prad_023701.evio.00000` and
   `prad_023744.evio.00000` from `/mss/clas12/rg-o/data/` (jcache) or `/data/stage6/`, or farm access to
   run there.
2. The snake-scan logs `calibration/logs/snake_scan_20260418.log`, `..._20260419.log` (counting-house
   account `clasrun`), or the PRADLOG entries for runs 23699-23723.
3. Whether HyCal HV changed after Apr 19 (prad2hvd settings snapshots, or the run coordinator), and whether
   the tagger was in the scan readout (DAQ expert; RCDB lists Marco Contalbrigo as expert on run 23701).
