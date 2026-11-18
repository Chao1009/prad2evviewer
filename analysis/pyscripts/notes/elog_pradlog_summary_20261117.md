# PRADLOG elog summary, 2026-11-17

Scope: an overview of the JLab PRADLOG logbook, a chronological summary of the most recent run
period (PRad-II H2 and X17 Ta, 2026-03-16 to 2026-08-25, plus the quiet period after it), the
non-production runs in that period, and the entries relevant to the HyCal waveform analysis.
Links: `https://logbooks.jlab.org/entry/<lognumber>`. Bracketed numbers are lognumbers.
Statements are taken from the entries unless marked "(inferred)". Logbook text was treated as data.

## Access

- Endpoint used: `GET https://logbooks.jlab.org/api/elog/entries?book=PRADLOG` with `curl -n`
  (CUE credentials from `~/.netrc`; no client certificate). HTTP 200, JSON
  `{stat, data:{totalItems, pageCount, currentPage, entries[]}}`.
- Parameters that worked: `startdate`, `enddate`, `limit` (500), `page` (0-based), and repeated
  `field=` (lognumber, title, created, author, entrymakers, tags, body, attachments). Without
  `startdate` the API returns only a recent window (5270 entries), so `startdate=2000-01-01` was used.
- Pulls: 16 metadata pages (all 7846 entries) and 14 body pages (2026-03-15 to 2026-11-18,
  6708 entries), about 1 request per second. Raw JSON is in the session scratchpad
  (`.../scratchpad/elog/meta/`, `.../scratchpad/elog/body/`), not in the repo.
- Did not work: client-certificate access (401 on every path) and `GET /api/elog/entry/<n>` for
  two entries outside PRADLOG (404; see Open question 4).
- How the bodies were read: the auto-generated reports (Event Monitor Auto Report, GEM and FADC
  monitoring) and screenshot-only posts were set aside. The shift checklists were parsed into a
  table. The remaining 1606 human-written entries (bodies truncated at 4000 characters) were read
  in five monthly slices. Run-type and configuration claims for the runs we hold locally were
  checked against the entry text.

## Overview

All 7846 entries reachable (lognumbers 3290264 to 4552399), dated 2014-07-28 to 2026-11-13.
PRad-I (2016) accounts for 834 entries and 2026 for 6803.

Entries per month, last 12 months:

| Month | Entries | Month | Entries |
|---|---:|---|---:|
| 2025-12 | 37 | 2026-06 | 647 |
| 2026-01 | 32 | 2026-07 | 2965 |
| 2026-02 | 58 | 2026-08 | 1228 |
| 2026-03 | 49 | 2026-09 | 2 |
| 2026-04 | 497 | 2026-10 | 0 |
| 2026-05 | 1320 | 2026-11 | 5 |

The API has two name fields: `author` (the account that submitted the entry; `clasrun` is the
shared counting-house account) and `entrymakers` (the people named on the entry; a multi-person
entry counts once per person).

| # | author | n | entrymakers | n |
|---:|---|---:|---|---:|
| 1 | clasrun | 4323 | wrightso | 132 |
| 2 | cpeng | 1211 | liyuan | 110 |
| 3 | wrightso | 313 | Zhaozhong Shi | 93 |
| 4 | nithya | 164 | vidura | 80 |
| 5 | liyuan | 145 | hauenst | 76 |
| 6 | xiongw | 107 | vimukthi | 73 |
| 7 | vidura | 91 | arunin | 59 |
| 8 | liye | 86 | xbai | 59 |
| 9 | xbai | 85 | tamangbm | 55 |
| 10 | yiyu | 64 | acschick, ddutta (tie) | 50 |

Tags (top 10; 3892 entries have no tag): Autolog 3726, DAQ 1081, ShiftSummary 66, EndOfRun 49,
Analysis 46, PRad 32, Detector 12, Hall B 10, Beamline 6, Readme 5.

## Run-period summary

The last 60 days (2026-09-18 to 2026-11-17) have 7 entries. Two (2026-09-21) post veto
scintillator noise and edge studies [4547923, 4547928]. Five (11-09 to 11-13) cover a HyCal
chiller test at 18 C for leaks: tank level stable, chiller off on 11-13, N2 left running
[4551947, 4551964, 4552251, 4552322, 4552399]. The summary below therefore covers the most recent
run period.

### Commissioning, 03-16 to 05-04 (runs ~23350-24310)

- 03-16 to 03-31: HyCal HV on and set to the PRad-I 2.2 GeV values. HV cable swaps were found and
  fixed [4477641, 4478453]. Crate sync aligned to ~2 ns [4480024]. The DAQ was re-cabled: the map
  changes at run 23353 and again after 03-31 [4481468, 4481937].
- 04-04: the FADC250 trigger pedestal register is 12 bits, so some pulse integrals sent to the VTP
  were 4096 too large. As a workaround `adchycal.cnf` was set to NSB=8, NSA=56 ns (16 samples)
  [4483623].
- 04-05 to 04-09: gain equalization with cosmics (target ~35 ADC peak). More than 90% of
  crystals ended within +-3 ADC. Gain-curve runs 23589 (-80 V) and 23590 (-40 V) [4485780, 4485808].
- 04-14 to 04-20: first beam, a 2.108 GeV photon beam at ~3 nA [checklists].
  - Crystal gain equalization with the bremsstrahlung edge, runs 23630-23648 [4487613, 4488060].
  - Snake scan with the raw-sum trigger, runs ~23692-23722 [4489057, 4489317].
  - Photon beam on a 4 um Ta foil, runs 23734-23740 [4489580].
- 04-20 to 04-24: GEMs installed, HyCal moved to its production position [4489905, 4489993].
  Crate timing improved from +-4 ns to +-2 ns [4490940]. The veto was moved to its own FADC board
  [4491513].
- 04-26 to 04-28: 2.108 GeV electrons on a 1 um carbon foil, runs 23892-24014 [4491621]. These
  include trigger and livetime tests, veto timing scans and a GEM latency scan. 04-29:
  calibration factors published [4492682].
- 05-04: NSA=128, NSB=8 ns set in `adchycal.cnf` to cover the whole peak [4494681].

### PRad-II H2 production, 05-05 to 06-24 (runs ~24365-25370)

- 05-05: first H2 (current scan run 24381). Raw-sum threshold set to 1200 MeV [4495147, 4495151].
- 05-07 to 05-10: 3.488 GeV H2 at 180-200 nA, with empty-target (Type B/C) cycles [4496279, 4497290].
- 05-11 to 05-18: 0.728 GeV with trigger `prad_v1_700MeV`. The raw-sum threshold went
  233 -> 400 -> 300 -> 250 MeV; 250 MeV from run 24603 on [4497644, 4498081, 4498915]. The run
  stopped on 05-15/16 for an unacceptable beam profile [4500894]. Plan completed 05-18 [4501395].
- 05-19 to 06-01: 3.488 GeV at 400-700 nA; trigger `prad_v2` from 05-21. Interruptions:
  - 05-20: 11 DAQ crashes in one shift [4501838].
  - 05-23 to 05-24: chiller leak trip, HyCal off ~22 h [4503309, 4503611].
  - 05-27 to 05-28: target blockage, ~19 h [4505976].
  - 3.5 GeV statistical goal met by 06-01 [4507193].
- 06-01 to 06-14: energy change. HV boards replaced (W187, W17) [4508097, 4509566]. adchycal5
  slot 14 FADC replaced, so pedestals for that crate need remeasuring [4511184].
- 06-15 to 06-23: 2.239 GeV at 320 nA (600 nA on 06-23); `prad_v3` [4511957, 4512110].
  - 06-16: FADC firmware `fadc250v3_PRAD_top_0F09` on adchycal1-7 [4512352].
  - 06-18: pedestals shifted, fixed by rebooting adchycal7 [4513127].
  - 06-20: FADC gain factors for 2.2 GeV [4514415].
  - PRad-II statistical goal reached on the evening of 06-23 [4515974].
- 06-23 to 06-24: "special tests and systematic studies" until the morning [4515974]. The detail is
  in two shift entries outside PRADLOG (4515733, 4515793) that could not be fetched.

### X17 transition and production, 06-24 to 08-14 (runs ~25549-26266)

- 06-24 to 06-30: kapton pipe, GEM shielding, VTP firmware aligning the cluster-sum and raw-sum
  trigger timing [4518386, 4518733].
- 07-02 to 07-03: beam at 2.2395 GeV. The FADC timing was fixed: SSP_PRAD_TRIGGER latency
  1540 -> 1668, FAV3_W_OFFSET 2700 -> 2828, "so the peaks would be in the same location as in
  PRad-II" [4519495]. A sequence of non-production runs followed (empty target, 50 um C
  calibration, trigger scans), then X17_v1 from run 25608 [4519304, 4519683].
- 07-04: HyCal gain fell ~15%, so the current was lowered to 3.5 nA [4519752]. First X17
  calibration: ep resolution 3.3%/sqrt(E) [4519744].
- 07-05 onward: production at 3.0 nA on 1 um Ta #1 (X17_v2, then X17_v3 from 07-09).
  - 3.5 nA from 07-20 to 07-23.
  - X17_v4 from run 26034 (07-30); X17_v4_RawSum950 from 08-02 [4520072, 4521520, 4528102,
    4534429, 4535999].
- FADC gain files followed the radiation damage, all relative to the start of X17:
  - `x17_25655` [4521328]; `x17_v3` on 07-14, inner modules 10-30% [4524659];
  - `x17_v4` on 08-08, inner modules 100-150% above the start, i.e. 50-70% damage [4539799].
- 08-10: HV of the innermost two crystal layers raised, giving 2-3x gain in LMS; gain file
  `x17_v4tov5` for those layers only [4540517, 4540647, 4540653]. Calibration run 26214 on 50 um C
  at 0.8 nA [4540550].
- 08-14: the inner-layer HV was lowered again "after noticing the energy resolution getting worse"
  [4543254]. FADC gains updated [4543286].
- 08-14: run officially ended [4543525]. Beam-off calibration and cosmic runs 26263 and 26266
  followed [4543262, 4543360]. Final cosmic run by 08-24, then HyCal HV and LMS off [4544511, 4544514].

Recurring problems, every phase: ROCs 100% busy (adchycal4 most often), runs marked "junk",
lost FADC/VTP readout from crate 5 slot 20 and crate 6 slot 15 [4532624], and accelerator
downtime.

## Entries relevant to the waveform analysis

Readout settings and timing:
- NSB/NSA: 8/56 ns from 04-04 [4483623], then 8/128 ns from 05-04 [4494681]. Since 04-05 the
  `adchycal.cnf` settings are stored in the EVIO files [4484048], so each run carries its own
  values.
- FADC firmware `fadc250v3_PRAD_top_0F09`, 06-16 [4512352]. Trigger latency and W_OFFSET change on
  07-03 [4519495]: X17 runs before about 25570 may have the pulse peak elsewhere in the window
  (inferred).
- Raw replay stores up to 200 raw samples plus firmware Mode 1/2/3 peaks. It flags peaks at the
  window boundary and NSB/NSA truncation [4493048].
- Timing:
  - Neighbour-crystal Δt width 0.4-0.5 ns at 3.5 GeV (run 24512), with secondary peaks near the
    absorber, possibly pile-up [4511316].
  - X17: low-amplitude (<50 ADC) signals have a >20 ns tail; some waveforms are oddly symmetric
    [4536882, 4538062].
  - Log-normal fit with 50% CFD gives ~0.14 ns per module, with side peaks at about +-1 ns
    [4532585, 4533549].
  - RF time is in the raw replay files [4495316].
- Trigger EnergySumIntWidth 32/24/16 ns, from random-trigger runs: 32 ns loses no good events,
  16 ns is too tight [4520653].

Pedestals and noise:
- FADC pedestal issues on 16 modules at commissioning [4480768]. W822 low pedestal [4495613].
- adchycal4 pedestals remeasured; W1039/W1074/W1107 read lower than seen in the event monitor
  [4505258, 4505274].
- adchycal7 pedestal shift on 06-18 [4513127]; pedestals loaded wrongly in run 25344 [4514877].

Gains, HV and damage:
- FADC firmware calibration is ~9% below the physics calibration [4493134].
- Gain files v2/v3/v4 (PRad-II) [4493355, 4494879, 4498651] and the X17 sequence above.
- LMS gain correction every 2000 LMS events in replay [4500901]. The LMS correction cannot fully
  follow the damage [4523754]. W355 shows a false LMS drop [4509268].
- Resolution worse by 0.2-0.3%/sqrt(E) with the first-layer gain loss [4523571, 4523754]; 3.5 nA
  slightly worse than 3.0 nA [4529939, 4531944].
- Run 26230 (08-12, which we use) falls between the inner-layer HV raise (08-10) and its partial
  reversal (08-14), the period in which the resolution was reported worse [4540517, 4543254].
- Dead or bad channels: W781 HV dead [4520940]; W1141 FADC channel dead [4503980]; W428 gain loss
  with HV suspected [4540826].

Rate and pile-up:
- 4 nA costs 2.5-3.0% livetime [4539980].
- Empty-cell run 24774 had 10-15 times the background of 24748 [4502801].
- Residual gas contamination 2% [4506389].

Resolution and other notes:
- 4 ns vs 16 ns clustering windows give similar resolution [4541597].
- Non-linearity from e-p and e-e peaks at three energies [4514001].
- Veto waveform pedestal method (first 30 samples, median) [4540355].

## Non-production runs

Types and run numbers as stated in the entries. "Junk" runs (DAQ crashes, busy ROCs) are not
listed one by one.

| Type | Runs (date) | Configuration stated | Lognumbers |
|---|---|---|---|
| Cosmics, beam off | 23356 (03-29); 23417 (04-01); 23555-23605 (04-06 to 04-10, HV iterations); 23856-23872 (04-25, `prad_v1_cosmic` 150 MeV); 24123 (04-28); 24173 (05-02); 24422 (05-07); 24617 (05-15); 25953-25960, 25983, 25992 (07-24 to 07-28, X17_ExternalsOnly); 26266 (08-15, `X17_v4_cosmic` 250 MeV); final cosmic run by 08-24 | target empty or home; GEM HV usually off | 4481504, 4482703, 4484501, 4485350, 4491395, 4491425, 4492490, 4493527, 4496207, 4499394, 4530901, 4532206, 4532864, 4543360, 4544511 |
| HV / gain-curve scan | 23589 (-80 V), 23590 (-40 V), 23591 (final) (04-09) | cosmics | 4485808 |
| Gain equalization with beam | 23532-23559 (04-05 to 04-08); 23630-23648 (04-15/16, row by row); 23660 | `prad_v0_equa`, bremsstrahlung edge, 3 nA photon | 4483977, 4487613, 4488060, 4488524 |
| Snake scan | ~23692-23718, redo 23719-23720, first glass layer 23721-23722 (04-17 to 04-19) | raw-sum trigger, ~3 nA photon beam | 4488783, 4489057, 4489317, 4489346 |
| Special target / energy calibration | 23734-23740 (04-19/20, photon on 4 um Ta); 23892-23915 (04-26, 1 um C, 1-3 nA, 2.108 GeV); 25202, 25286-25308 (06-15/18/19, 1 um C, ~5 nA, veto in); 25556, 25559-25562 (07-03, 50 um C, 0.5 nA, X17_v0_calib); 25654-25657 (07-04, 50 um C); 26214 (08-10, 50 um C, 0.8 nA) | as listed | 4489580, 4491621, 4511707, 4513238, 4513708, 4519304, 4519946, 4540550 |
| Random trigger | 24014 (04-27, 1 um C, 10 nA); 24575, 24578 (05-13/14, 0.7 GeV, `prad_random_700MeV_v0`, 2 h); 24601-24602 (05-14, `prad_random_thr_100_300`); end of 3.5 GeV (05-31/06-01, number not given); ~25635-25639 (07-04, 3.5 nA Ta, 32/24/16 ns sum windows) | as listed | 4492181, 4498531, 4498915, 4506827, 4506940, 4519727, 4520653 |
| Empty target / target out (Type B empty cell, C empty cell and chamber, D target out) | many, interleaved with production, e.g. 24426-24434, 24504-24531, 24606-24614, 24634-24657, 24714-24715, 24763-24774, 24874, 24941-24945, 24962, 25215, 25223, 25289-25291, 25309; X17: 25550-25555 (1 nA), 25667, 25672, 25864 (3 to 7 nA ramp) | H2 current of the period, or 1-7 nA in X17 | 4496335, 4497290, 4498915, 4500390, 4501870, 4502801, 4506389, 4506858, 4506940, 4512110, 4513932, 4519304, 4520072, 4526452 |
| Beam-off external-trigger HyCal runs (LMS/alpha, ~27 Hz) | 25650 (07-04, radiation-recovery study); ~25776, 25783 (07-10); 25833 (07-14/15); 25837; 25848 (07-17, "special run"); 26263 (08-14/15, "for HyCal calibrations") | PRAD_HYCAL or X17_ExternalsOnly, GEM HV off | 4519946, 4522039, 4522372, 4524816, 4525279, 4525851, 4543262 |
| LMS runs | 24491, 24493 (05-10, gain monitoring); 26103 (08-04, radiation recovery) | not stated | 4497033, 4537540 |
| Trigger / threshold scans | 23445, 23453 (04-03/04); 24203-24215 (05-03, veto timing, W_OFFSET 2500-3150); carbon total-sum 800/1200/1500/1800 MeV (05-05); 25575-25591 (07-03, E_max 2000-1100 MeV); 25592-25607 (07-03, raw sum 1400-900 MeV); 25658-25663 (07-05, prescale); 26171-26173 (08-08, raw sum 950/1050/1150 MeV) | as listed | 4483413, 4483454, 4493889, 4495161, 4519521, 4519543, 4519578, 4520082, 4539980 |
| Current scans / high current | 24381 (05-05, 5-100 nA H2); 06-23 at 600 nA (runs not given); 25862-25863 (2 nA); 25896 (3.5 nA); 25897 (4 nA) | as listed | 4495151, 4515486, 4526452, 4528130 |
| Veto / GEM tests | veto delay scans 23916-23934 and 25263-25285; GEM latency scans 23935-23962 and 25716-25717; GEM HV scan 25640-25649; GEM pedestals 24111, 24162, 24392, 24495, 24783-24787, 24872, 25026-25035 | GEM pedestal runs with random trigger read all crates, FADC included (25026-25035) | 4491780, 4491812, 4492610, 4495593, 4508906, 4512726, 4513238, 4519827, 4520526 |

What each type can do for the waveform analysis (inferred; none tested yet):
- **Random trigger** (24578, 24601-24602, 25635-25639): unbiased samples of the readout window give
  per-module pedestal and noise, out-of-time pulses for an unbiased pile-up rate per module, and
  baseline shifts under beam. 24578 is the cleanest such run, at 0.7 GeV.
- **Beam-off external-trigger runs and LMS runs** (25650, 25776-25848, 26263, 26103): pedestal and
  noise per module without beam, and fixed-amplitude LMS pulses for a pulse-shape template per
  module that is free of pile-up. Repeated over July and August, they track how damage changes
  amplitude and shape.
- **Cosmics** (23356-24617; 25953-25992 and 26266 in X17): low-rate, pile-up-free pulses of a few
  tens of MeV. Pile-up-free shapes at low amplitude can test the reported >20 ns tail of small
  signals [4536882]. 26266 has the X17 trigger timing.
- **HV gain-curve runs** (23589-23591) and the 08-10/08-14 inner-layer HV changes: dependence of
  amplitude and shape on HV. Runs 26214 to ~26262 sit at raised inner-layer HV; 26263 onward at
  the lowered HV.
- **Trigger and timing scans** (07-03 latency and W_OFFSET change; 24203-24215 W_OFFSET scan):
  locate the peak inside the readout window. Our in-time window of 140-156 ns should be checked
  against X17 runs before and after 07-03 (around runs 25550-25580).
- **Empty target / target out** (Type B/C/D cycles): beam-related background, occupancy and pile-up
  rate without the H2 target. They pair with our production runs at the same current.
- **Carbon calibration and current scans** (26214, 25559-25562, 25862-25863, 25897, 24381): extra
  resolution points at low rate (0.5-0.8 nA) versus 3-4 nA. They separate pile-up from intrinsic
  resolution.
- **Snake scan and photon Ta runs** (23692-23722, 23734-23740): position-dependent response and
  bremsstrahlung-edge amplitudes. They also give high-amplitude pulses for saturation checks.

Data availability and our usage:
- No entry states tape or cache status for any run. Only "junk" and "do not use" labels appear.
  Raw replay of 24578 on ifarm lustre is mentioned [4498717].
- Local evio files (`/Users/L00338853/work/PRad/data/evio/`) that are non-production by the logbook:
  - **24578** is a 2 h random-trigger run at 0.7 GeV [4498531], not H2 production as listed in our
    `data/evio/README.md`. It is good for pedestal and pile-up studies; it should not be pooled
    with the 0.7 GeV production runs.
  - **25308** is a 1 um carbon run with veto in at 5.3 nA, logged as "Production Carbon" in the
    shift table [4513868]. It was used for the cluster-trigger efficiency study [4514402]. Our
    README lists it as a carbon run, consistent with this.
  - **23701 and 23744** are not named in any entry. 23701 falls inside the snake-scan range and
    23744 just after the Ta photon runs (both inferred).
  - **24386** was taken on 05-05 at 21:03 during H2 commissioning, after the current scan in run
    24381 [4495151]. It was the first H2 run replayed and analysed (ep/Moller yields, background
    estimate) [4495336, 4495575]; 21:03 is the auto-report time, not the run start.
  - **24503** has no human entry.
  - **24590, 24937, 25800** are not named in any human entry (bodies truncated at 4000 characters;
    see Gaps).
- Production runs confirmed by the logbook:
  - **24603**: H2, 0.7 GeV, 25 nA, 250 MeV threshold [4498915].
  - **24512**: H2 at 3.5 GeV on 05-10, just before an empty-target run [4497290]; used for the
    peak-time study [4511316].
  - **24698**: H2 at 400 nA [4501711].
  - **24911**: H2, Configuration A [4506221, 4506439].
  - **25228, 25257, 25260, 25299, 25320**: H2 at 2.2 GeV, 320 nA [4512110, 4512559, 4512726,
    4513708, 4514290].
  - **26138** (08-06) and **26230** (08-12): X17 Ta. For 26230 see the HV note above.
- None of the beam-off external-trigger, LMS, cosmic, empty-target or carbon calibration runs is
  in our local set.

## Open questions for the user

1. Is a run-type list available from RCDB (clasweb.jlab.org/rcdb) that we can use instead of
   shift-log text? The logbook gives types for most but not all runs.
2. Which non-production runs should be fetched first? Suggested set:
   - 24578 (already local);
   - one beam-off external-trigger run per month in X17 (25650, 25833, 26263);
   - X17 cosmics 26266;
   - carbon calibration 26214;
   - empty-target run 25864 (3 to 7 nA ramp).
3. Does our energy-resolution analysis account for the inner-layer HV change of 08-10 to 08-14
   for run 26230?
4. Should the 06-23/24 special-test runs (detail in entries 4515733 and 4515793, outside
   PRADLOG) be pulled?

## Gaps

Not obtained on 2026-11-17:

- Bodies were read truncated at 4000 characters. Long shift summaries lose their run tables, so
  some run ranges are incomplete (for example runs 25995-26048).
- Attachments (spreadsheets, plots, HV and gain tables) were not opened.
- The two special-test shift reports outside PRADLOG (lognumbers 4515733, 4515793, 06-23/24) were
  not obtained.
- Data-file availability (cache vs tape) was not obtained: no entry states it for any run.
- Auto-generated monitoring reports (1075 runs, 24361-26275) were used only for run dates.
- Some run numbers are marked "inferred" or "not given" where the entries are unclear.

User decisions 2026-11-17: fetch RCDB run-type list (in progress); start with local run 24578, then one beam-off run and 26266 cosmics; re-examine run 26230 against the 08-10/08-14 HV change; pull the two special-test reports.
