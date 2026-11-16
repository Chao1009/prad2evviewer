# Plan P14: snake-scan and cosmic constants through LMS ratios (comparison against control v3) — 2026-11-14

Working dir `analysis/pyscripts` (branch `waveform_jz`); paths are relative to it unless they start with `/`. BM =
`waveform_analysis/energy_resolution_benchmark.py`; ER = `waveform_analysis/energy_resolution`; SS =
`waveform_analysis/snake_scan/` (new); V3 = `output/energy_resolution/control_phase2_v3_20261114/`; DB =
`/Users/L00338853/work/PRad/prad2evviewer/database`; EVIO = `/Users/L00338853/work/PRad/data/evio/`; CAL =
`/Users/L00338853/work/PRad/prad2evviewer/calibration`. Parent:
`notes/plan_energy_resolution_improvement_20261112.md` (package P14, U4; D3 comparison rule; P7 ring machinery).
Research: `notes/calibration_runs_snake_cosmic_20261112.md` (P13). Environment: `source setup_env.sh`
(`PRAD2_DATABASE_DIR`, `PYTHONPATH=build/python`), `MPLBACKEND=Agg`. Facts below marked [W1]-[W5] come from the
planner-worker probes of 2026-11-14 (first 3 k to 100 k events of `prad_023701.evio.00000-00080` and
`prad_023744.evio.00000`; scripts in the session scratchpad, nothing in the repo).

## Time estimate

| | estimate |
|---|---|
| Agent time | about 3 agent-days: P14.1 0.3, P14.2 0.6, P14.3 0.5, P14.4 0.3, P14.5 0.5, P14.6 0.5, plus review |
| CPU, local | 23701 dwell map 9 files x 330 s = 0.8 h (6 min wall at 9 lanes); edge pass 25 files x 330 s = 2.3 h (20 min wall at 9 lanes); LMS pass over the local production splits of the 5 run types, about 130 splits x 2-5 min = 5-10 h (1 h wall at 8 lanes); 23744 3 files x 20 min = 1 h; benchmark reports 2 x 45 min (2.24 GeV) + 2 x 20 min (3.49 GeV) |
| Wall | 3-4 working days with the 16 extra 23701 files cached locally (34 GB; the evio volume had 15 GiB free on 2026-11-14 evening, so about 20 GB must be freed or the files fetched in batches of about 6, processed and deleted); the full scan (P14.7) adds a farm pass of about 160 CPU-h |
| Critical path | P14.1 map + dwell table (0.5 d) -> P14.2 edge pass (0.5 d) and P14.3 LMS (parallel, 0.5 d) -> P14.4 transfer + pool validation (0.3 d; the go/no-go for P14.6) -> P14.6 hook and D3 comparison (0.5 d). P14.5 cosmic runs in parallel after P14.1 |

Cost basis: decode 253-414 s per EVIO file single-core on this machine (`output/energy_resolution/logs/scan_*.log`
"[done]" lines; 18 cores) [W1]; a peak-only Python pass is 1.8 ms/event, 330 s per 0.65 M-event file [W2]; a
WaveAnalyzer pass with a few channels per event ran at 545 events/s [W5]; a 2.24 GeV report takes about 45 min [W4].

## Open Questions

Material decisions for the user; each has a recommended answer. The plan proceeds on the recommendations if
approved as is.

- Q1 (full scan). Run 23701 covers about 38 modules on rows 17-18 of the left half (W545-W564, W580-W599; T, O, M, H
  and 7 pool modules) [W2]. O and T constants for the other three quadrants need the other scan runs (23699-23723,
  1.02e9 events, about 3.3 TB, JLab only). Recommended: approve P14.7 now as a farm job that stores only the
  per-dwell edge histograms and the LMS/alpha/pedestal series (a few MB per run), run after P14.4 shows the
  transfer spread is below the 2 % threshold (D5), and start with every 3rd file of each run (a module spans 2.5-3
  files, so every module keeps at least one file). Caching 22 runs locally is out (space).
- Q2 (extra 23701 files). Recommended: fetch the 16 files `prad_023701.evio.000{03,06,13,16,23,26,33,36,43,46,53,56,63,66,73,76}`
  (34 GB) so that every module of the run has at least one file in a clean dwell; files 01-02 and 79-80 only if the
  first and last modules matter [W2]. Without them 9 modules have usable dwells and 2 only ramps. Space: the evio
  volume (321 GB of EVIO splits) had 15 GiB free at planning time, so either free about 20 GB or fetch in batches
  of about 6 files, run P14.1-P14.3 on them, delete, repeat (the per-file npz outputs are small).
- Q3 (spot size). Unknown; the data show 97 % of events with the central module as the maximum and a neighbour
  above 30 % of it in 1.6-5 % [W2], consistent with a spot well inside 20 mm. Recommended: if known, record it in
  `P14_snake_scan.md`; otherwise no input is needed, D3 fits the 3x3 sum, which is insensitive at that level.

## Goal

Per-module gain constants for the PbWO4 modules the April 2026 snake scan covers, transferred to the production
HV through LMS ratios, validated on the elastic-pool modules the scan also covers (U4), and offered to the
benchmark as a comparison-only prior for the O and T rings (the regions the pool solve leaves on a scalar prior,
`prior_frac` 1.0, res_cb 1.88 % / 2.14 % at 2.24 GeV [W4]). Cosmic run 23744 gives a MIP point per module group
as a cross-check (P14.5).

Non-goals: changing any benchmark default; the absolute energy scale (stays with the elastic pool, user decision
2026-11-14); a per-module non-linearity (one endpoint energy only); position dependence across a module (no
spot-size knowledge); the production consumer of ring-scoped constants; processing the full scan locally.

Expected change size: one new script `SS/snake_scan.py` (about 500 lines) and one `SS/cosmic_mip.py` (about 200),
one derived map JSON, a benchmark hook of about 30 lines (BM:1632-1676, `OVERRIDABLE`), tables under `SS/out/`
and `ER/improvement/P14_*.md/.tsv`; no change to existing outputs.

## Approach

Read run 23701 with the existing `dec`/`det` bindings and a derived HyCal map for runs <= 23849 (D1). Reconstruct
the beam path from the EPICS `ptrans_x/y` snapshots and the per-block occupancy argmax (D2). For each dwell, fit
the bremsstrahlung endpoint edge of the 3x3 sum around the beam module, in MeV of the production calibration
file, with single-pulse (pile-up-free) events (D3); the edge ratio `s_i = E_edge/2108.85` carries the module's
gain error at the scan HV times a common containment factor. Measure the LMS amplitude per module in the scan
epoch (23701 LMS events plus 23744's, alpha-normalised references) and in each production run type from the local
splits (D4). Transfer: `g_T,i = N (L_scan,i / L_prod,i) / s_i`, with N fixed by the pool modules the scan covers,
so the containment and the global scale never need to be known (D5). The spread of `g_pool,i / g_T,i` over the
pool modules is the figure of merit and the gate for P14.6 (D5). P14.6 adds `GAIN_PRIOR_FILE` to the benchmark
(vector g0 and lambda under `TRAIN_REGIONS`, D6) and compares O/T region metrics against v3 by D3 (D7).

Rejected: the equaliser's 5 %-log-cumulative edge (`CAL/gain_scanner.py:52-140`), not the endpoint and
peak-height based [W4]; absolute ADC/MeV transfer without LMS (HV changed, user decision); the benchmark's own
scan path for the dwell data (matched filter per hit, 330 s/file, needs a template bank and `RUN_TYPE_OF` entries
in a file P1b edits concurrently) in favour of a standalone script on `WaveAnalyzer`, the same integral the
production calibration uses; `gain_response_curve.json` gain(HV) curves to bridge the HV change (no HV history,
no consumer in the repo [W3]).

## Decisions

- D1 (readability and map). Run 23701 decodes with the current build; `database/daq_map.json` does not exist:
  the DAQ address is the `daq` field of `DB/hycal_map.json`, loaded by `HyCalSystem::Init`
  (`prad2det/src/HyCalSystem.cpp:81,114-117,143`). The history file
  `DB/history/daq_map_20260424_til_run023849.json` is a flat name/crate/slot/channel list and is rejected by that
  loader, so the switch is a derived map: `hycal_map.json` with each module's `daq` replaced by the history entry
  of the same name, written to `SS/hycal_map_til_run023849.json` (not under `DB/`) and passed with
  `det.PipelineBuilder().set_hycal_map(path)` (`python/bind_det.cpp:1187`; BM `build()` at BM:808-812 does not
  call it). Only 6 entries differ: LMS1-3 and LMSP at crate 3 slot 3 ch 5-7/4 (production slot 13 ch 10-12), W155
  (3,3,14) and W188 (3,3,15) [W1]. Events: 1 per block, trigger bits 0x100 on 99.97 % of events (the benchmark's
  `PHYSICS_BITS`), sparse readout of 3-33 channels of 100 samples at 4 ns, lead glass included, run number 23701,
  `Ebeam` 2108.8 MeV from `DB/runinfo/general.json` (`from_run: 0`); EPICS events (tag 31, about 5 per 3 k events)
  carry `MBSY2C_energy` 2108.853, `ptrans_x`, `ptrans_y`, the HyCal transporter encoders and the IPM positions [W1].
  Rate 14-22 kHz from the TI timestamps; one file is 0.65 M events, 33-42 s [W2].
- D2 (beam path). Primary: the EPICS `ptrans_x/y` series, which is the HyCal position; the module under the beam
  is the one whose centre is nearest to the beam at the transporter setting (sign and offset calibrated once from
  the occupancy argmax of the clean dwells; executor records the fitted affine map and its residual, expected
  < 5 mm). Validation: per 5 k-event block, the module with the largest pedestal-subtracted peak sample
  (`max(x - median(x[:8]))`, no WaveAnalyzer) and its share; a dwell is a maximal run of blocks with the same
  argmax at share >= 90 %; ramp blocks (share < 90 %, or EPICS position moving) are excluded. Measured: dwells are
  97 % pure, transitions are slow ramps of about 18 s (350 k events), a module takes 86-108 s of beam (1.5-2 M
  events, 2.5-3 files), not the 60 k events of the brief; the run walks `snake-all-pwo-r` (`CAL/paths.json`,
  `CAL/scan_engine.py:31`, dwell default 120 s at :22) in descending index order, W599 (file 00) to W564
  (file 80), i.e. row 18 descending to W580, the turn through W579/W545/W546, row 17 ascending to W564 [W2]. The 9
  local files give 9 whole-file dwells (W599, W594, W590, W585, W580, W547, W552, W557, W564; W564 weak, 29 % of
  events with a peak above 300 ADC) and 2 ramp partners (W589, W586). The dwell table (`SS/out/dwell_23701.tsv`:
  file, first/last event, module, share, n_events, ptrans_x/y) is the input of every later package.
- D3 (per-module measurement). Spectrum: not a peak; a trigger turn-on at 4-6 k ADC, a plateau, a sharp
  endpoint fall (half-plateau crossing at 18.5 k ADC of the 3x3 integral sum for every probed module, i.e. 0.114
  MeV/ADC at 2108.85 MeV, 14 % below the calibration file's median W of 0.1325 with a different integral
  convention) and a pile-up tail to 45 k ADC [W2]. Fit quantity: the 3x3 sum around the dwell module of
  `WaveAnalyzer` integrals (`WaveAnalyzer.analyze(samples)` -> ped, rms, peaks; the production calibration's
  integral convention, `python/bind_dec.cpp:412`) converted to MeV with the run's calibration file W_i
  (`calibration_factor_3p5_June7_nonLinear.json` per `DB/runinfo/general.json`), events with exactly one peak on
  every member that fired (pile-up rejection; 17 kHz x 400 ns gives 0.7 % naive overlap but the per-module photon
  rate is higher, the tail shows about 10 % [W2]). Edge model: a linear plateau times an error-function step at
  E_edge with width sigma, fitted over [0.6, 1.3] x the half-plateau crossing; free: plateau level and slope,
  E_edge, sigma; the pile-up residue above the edge is a free constant. Output per dwell: `s_i = E_edge/2108.85`
  (containment C and the global scale inside), its fit error (statistical error negligible at 1.5 M events; the
  systematic from the window and the slope term is taken from the spread over window choices, recorded), and the
  central-module-only edge as the spot-size cross-check. Neighbour gains: the 3x3 members enter with W_j (production
  calibration); since the scan visits the neighbours, a second pass replaces W_j by W_j / s_j for members with a
  dwell, and the executor reports the change (expected < 1 %, 20 % neighbour weight). Why 3x3 and not the central
  module: a centred 2 GeV shower leaves about 80 % in the central module with a steep dependence on the spot
  offset and 93-95 % in the 3x3 [W4]; why not the full island cluster: clustering lives in the C++ pipeline and the
  report side, and C cancels in D5 anyway.
- D4 (LMS). LMS events: Physics events with bit 24 set (always `0x1000100`, `trigger_type` unreliable), nearly all
  channels firing (1567 channels, PbWO4 median peak about 370 ADC in both 23701 and 25320); production rate about
  9 Hz (61-72 per 8 s split), run 23701 only 1-8 per 35 s split, about 25 events in the 9 local files, about 250 in
  the whole run; alpha-source events bit 25 (`0x2000000`, 144 per split in 23701) fire the 3 reference PMTs
  [W3]. Per module and event: `WaveAnalyzer` integral with `npeaks == 1`; per event normalised by the mean of
  LMS2 and LMS3 integrals each divided by its alpha-source integral of the same epoch, the `GainCorrCompute`
  convention (`analysis/src/GainCorrCompute.cpp`; Python template `analysis/pyscripts/lms_gain_variation.py:
  extract_split`, which selects on `trigger_bits & (1<<24)` and normalises `W_median x mean_{j=2,3}(1/R_j)`) [W3].
  Reference channels are mapped by name through the derived map (D1), never by slot. Epoch values: L_scan,i =
  median over the 23701 LMS events, with the 23744 LMS events (about 200 in 3 files, 0.04 % of events [W5]) as a
  second estimate; the two are merged when they agree within 2 SE per module over the scanned modules, else
  23701 alone is used and the disagreement reported. L_prod,i,rt = median over the LMS events of the local splits
  of the run type's runs (`RUN_TYPE_OF`, BM:63-69; local splits are every 10th-100th file [W3]), one value per run
  type, because the prior is used in that run type's report. Error per module: 1.25 x 1.48 MAD/sqrt(n). Existing
  tables `DB/gain_factor/gain_correction/*.root` and `ref_gain/*_LMS.dat` are relative to a reference run and
  exist for only some pool runs [W3]; not used.
- D5 (transfer and validation; the units). Let g_i be the benchmark's dimensionless multiplier on W_i E
  (`gains.tsv` columns moller, nonlinear, run, module, g, sigma_g, n_seed, n_member, prior_frac; median 0.965,
  pool modules `prior_frac < 0.5` have sigma_g about 0.002 [W4]). At the scan the module needs the correction
  a_i = C/s_i; a PMT gain change multiplies the ADC of signal and LMS light alike, so the production correction is
  a_i x L_scan,i / L_prod,i. Transferred multiplier: `g_T,i = N_rt (L_scan,i / L_prod,i,rt) / s_i`, with
  `N_rt = median over the scanned pool modules of g_pool,i s_i L_prod,i,rt / L_scan,i` (absorbs C and the global
  scale; the absolute scale stays with the pool). Figure of merit: the robust spread (1.48 MAD) of
  `r_i = g_pool,i / g_T,i` over the scanned pool modules (7 in 23701: I1-I3 modules of rows 17-18, cols 12-14 and
  21; H modules, 12, are added as a second group since they carry solved gains in `ring/gains.tsv`), against the
  expected quadrature of sigma_g, the edge error and the LMS errors. Gate for P14.6: spread <= 2 % (the alternative
  for O/T today is a scalar g0 whose error is the module-to-module scatter, 0.90-1.00 in the free modules);
  sigma_g of the prior = max(spread, per-module propagated error). Above 2 %: P14.6 is not run; report the
  attribution (edge vs LMS by the two groups).
- D6 (benchmark hook; comparison-only). New constant `GAIN_PRIOR_FILE = ""` in `OVERRIDABLE` (BM:241-244): a TSV
  (module, g, sigma_g; relative path under `TABLE_DIR`) read once in `report()`. In `solve_gains` (BM:1611) build
  `g0v = full(n, g0)` and the per-module lambda next to the scalar prior at BM:1632-1633; for listed modules whose
  region is in `TRAIN_REGIONS`, and only when `TRAIN_REGIONS` is non-empty (the pool-only solve swaps it to `()` at
  BM:2217, so pool rows stay byte-identical to the control), set `g0v[i] = g` and `lamv[i] = 1/sigma_g^2` (the
  existing `GAIN_RING_PRIOR_X` multiplier at BM:1668 then still applies); use `g0v[uniq]` at BM:1673 and
  `gnew = g0v.copy()` at BM:1676; `solve_gains_seed` (BM:1694, copy at :1703-1704) gets the same vector. The
  variant selects O and T with `--set 'RING_CALIB_BY_RT={"2p24GeV_H2": ("H","M","O","T")}'` (BM:2767-2768 turns
  it into `TRAIN_REGIONS`, `GAIN_REGION_SCOPE`, `GAIN_COMBINE_RUNS`). Defaults unchanged: with `GAIN_PRIOR_FILE ""`
  every path is the v3 code path (checkpoint: a default 2.24 GeV report byte-identical to one made from the BM
  state immediately before the hook edit, which is V3 if P1b changed no default). The hook is written only after
  the P1b executor has landed its BM edits or with the user's go-ahead (Owns conflict, Scope Notes).
- D7 (comparison). Control `V3/v3_defaults/2p24GeV_H2`; variant `output/energy_resolution/variants/P14/prior/`
  from `report --out ... output/energy_resolution/scan/2p24GeV_H2 --set RING_CALIB_BY_RT=... --set GAIN_PRIOR_FILE=...`;
  `compare_reports.py <control> <variant> --paths g16` -> `delta.tsv`; the O and T rows are the target region
  (D3 of the parent: real |Δ| >= 2 SE in `res_cb` or `tail_hi3`/`leak_cut`, guards `eff_cut`, `tail_lo_f`, yield);
  replication at 3.49 GeV (`3p49GeV_H2`, its own L_prod). Pool rows must show Δ = 0 (byte-identical). Attribution:
  the P7 ring arm with `("H","M","O","T")` and no prior (`output/energy_resolution/variants/P7/`,
  `ER/improvement/P7_beyond_pool.md`, O/T null) is the second control that isolates the prior's effect from the
  O/T ring solve; cite its O/T rows, no new report.
- D8 (HV/LMS segments within the scan). Per file: the pedestal mean per module from `WaveAnalyzer`
  (`ped_nsamples` 30, `DB/daq_config.json:109-111`) on the physics events and the LMS/alpha series with the TI
  timestamp. A segment boundary is a step > 3 sigma of the file-to-file scatter in the LMS2/LMS3 alpha-normalised
  reference or in the median PbWO4 pedestal; the edge and LMS epochs are then split at it and reported. Expected:
  none inside 23701 (the equalisation ended with 23698 [P13]).
- D9 (cosmic, P14.5). Geometry `DB/hycal_module.txt:10,586`: PbWO4 20.77 x 20.75 x 180 mm, lead glass
  38.15 x 38.15 x 450 mm; cosmics cross a module transversely (2.05 cm, 16.97 g/cm2): mean MIP deposit 20.9 MeV,
  most probable about 19-21 MeV (Landau MPV below the mean; relativistic rise and slant +5-10 %); the repo's
  `adc_to_mev_factors_cosmic.json` uses 25 MeV with no documented origin [W5]. Quantities: the track selection
  uses the peak sample (`max(x - ped)`, MIP median 31.5 ADC, threshold 10 ADC); the MIP peak is fitted on the
  `WaveAnalyzer` integral, the edge's quantity (D3), expected MPV about 185 ADC at 0.114 MeV/ADC (the repo's
  cosmic target is 250 ADC for 25 MeV, `analysis/tools/cosmic_test.cpp`). The existing `cosmic_test.cpp`
  selection (5 < nch <= 70, every single-peak hit, no track requirement) is contaminated, so P14.5 selects vertical
  tracks: >= 3 contiguous modules in a column, <= 2 off-column hits, each hit > 10 ADC (0.7 % of events, median
  track length 7 modules, interior hits 1.05 per module per 30 k events [W5]). With 3 files (about 19 hits per
  module) P14.5 reports MIP peaks per region group (T, O, M, H, pool; >= 20 modules each) as `s_mip,group` and the
  group-level transfer ratio against D5's `r_i` groups; per-module MIP constants need the whole run (P14.7). Run
  23744 LMS events (0.04 %, bit 24 with nch > 900; 910 further nch > 900 events with bits 0x100 only are flagged and
  excluded) feed D4.
- D10 (standalone script, not the benchmark's scan path). `SS/snake_scan.py` with subcommands `map`, `path`,
  `edge`, `lms`, `transfer`, each a pure function over one EVIO file or over the tables, parallel over files with
  `xargs -P` as `run_energy_resolution.sh:27-30` does for runs; decoding through `_common.setup_pipeline` and the
  `dec` bindings as `lms_gain_variation.py` does. Reason: BM is under concurrent edit (P1b), needs a template bank
  per run type, and the dwell data need no matched filter.

## Outline

- `SS/snake_scan.py` (new, about 500 lines):
  - `derive_map(hycal_map: Path, history: Path, out: Path) -> dict` (D1): keys as `hycal_map.json` entries
    (`n`, `t`, `geo`, `daq`); 6 changed entries, LMSP added with `t: "LMS"`.
  - `dwell_table(evio: Path, map: Path, block: int = 5000) -> rows(file, ev_first, ev_last, module, share,
    n_events, ptrans_x, ptrans_y, ts_first, ts_last)` (D2); EPICS parsing as the worker probe did (tag 31 text,
    name after value).
  - `edge_pass(evio, dwell_rows, map, calib) -> per-module histogram npz (3x3 MeV, central MeV, n_single,
    n_pileup)` and `fit_edge(hist) -> (E_edge, sigma, plateau, slope, err, window_spread)` (D3); histogram
    binning 10 MeV from 0 to 4000 MeV.
  - `lms_pass(evio, map) -> npz (lms_I[n, module], lms_ref[n, 3], alpha_ref[n, 3], ts, ped[file, module])`
    (D4, D8), the `lms_gain_variation.py:extract_split` shape with the map and reference channels by name.
  - `transfer(edge.tsv, lms_scan.tsv, lms_prod_<rt>.tsv, gains.tsv, ring/gains.tsv) -> SS/out/transfer_<rt>.tsv
    (module, region, s, s_err, L_scan, L_prod, g_T, sigma_T, g_pool, r) and SS/out/prior_<rt>.tsv (module, g,
    sigma_g; O and T modules only)` (D5).
  - Invariants: modules addressed by name; every table carries a `# source:` line with the file list and the map
    path; no module outside the dwell table gets an edge.
- `SS/hycal_map_til_run023849.json` (derived, D1); `SS/out/` tables; `ER/improvement/P14_snake_scan.md/.tsv`
  (results note, same layout as `ER/improvement/P7_*`).
- `SS/cosmic_mip.py` (new, about 200 lines): `tracks(evio, map) -> hits(module, peak_adc, integral_adc, column_len)`,
  `mip_peaks(hits, groups) -> (group, mpv, err, n)`; Gaussian core fit in [0.7, 1.5] x the mode as
  `analysis/tools/cosmic_test.cpp` does, on the vertical-track sample (D9).
- BM hook (D6): `GAIN_PRIOR_FILE` constant next to `AMP_SCALE` (BM:240), `OVERRIDABLE` (BM:243), loader next to
  `_AMP` (BM:2372, set in `report()` before the workers fork), vector prior in `solve_gains` (BM:1632-1633, :1668,
  :1673, :1676) and `solve_gains_seed` (:1703-1704). Gotcha: `corrections()` sets `TRAIN_REGIONS = ()` for the
  pool-only solve (BM:2217); the prior must be inert there.
- Pattern to follow for the tables and note: `ER/improvement/P7_*`, `compare_reports.py` usage header.

## Slices

1. Map and dwell table (P14.1). Delivers `SS/hycal_map_til_run023849.json` and `SS/out/dwell_23701.tsv` for the
   9 local files, plus the EPICS-to-module affine map. Units: `derive_map`, `dwell_table`. Tests:
   `map_diff: hycal_map.json + history -> exactly 6 changed daq entries, names LMS1-3, LMSP, W155, W188`;
   `dwell_00: file 00000 -> one dwell W599, share >= 0.95, n >= 600 k`;
   `dwell_30: file 00030 -> ramp W586 then dwell W585, ramp blocks excluded`;
   `epics_module: ptrans of a clean dwell -> the dwell module within 5 mm`. Checkpoint:
   `python SS/snake_scan.py path EVIO/prad_023701.evio.00000` prints the dwell rows above.
2. Edge pass and fit (P14.2). Delivers `SS/out/edge_23701.tsv` (module, s, s_err, window_spread, n_single,
   n_pileup, s_central) for every dwell module. Units: `edge_pass`, `fit_edge`. Tests:
   `edge_w547: file 00050 -> E_edge within 3 % of 18.5 k ADC x W (about 2.1 GeV x s), sigma < 5 % of E_edge`;
   `edge_pileup: n_pileup / (n_single + n_pileup) between 2 and 20 %`;
   `edge_windows: three windows -> s spread < 0.5 %`;
   `edge_consistency: s of the 4 probed modules (W547, W552, W580, W599) agree within 1 % [W2 saw identical edges]`.
   Checkpoint: `python SS/snake_scan.py edge --dwell SS/out/dwell_23701.tsv EVIO/prad_023701.evio.00050`.
3. LMS and segments (P14.3). Delivers `SS/out/lms_scan.tsv` (module, L_scan, err, n, source 23701/23744/merged),
   `SS/out/lms_prod_<rt>.tsv` for the 5 run types, `SS/out/segments_23701.tsv` (D8). Units: `lms_pass`, a
   `segments` reducer. Tests:
   `lms_count_23701: file 00010 -> 8-9 LMS events, 1569 channels, PbWO4 median peak 300-450 ADC`;
   `lms_ref_by_name: LMS2 of 23701 resolves to (3,3,6), of 25320 to (3,13,11)`;
   `lms_prod_25320: one split -> 61-72 LMS events, L_prod for >= 1140 PbWO4 modules`;
   `segments_none: 9 files -> no step > 3 sigma (else listed)`. Checkpoint:
   `python SS/snake_scan.py lms EVIO/prad_023701.evio.00010` and the same on one 25320 split.
4. Transfer and pool validation (P14.4). Delivers `SS/out/transfer_2p24GeV_H2.tsv`, `prior_2p24GeV_H2.tsv`,
   the 3.49 GeV pair, and `ER/improvement/P14_snake_scan.md` with the spread, its attribution and the go/no-go.
   Unit: `transfer`. Tests: `transfer_units: a module with s = 1, L_scan = L_prod -> g_T = N`;
   `transfer_pool_median: median r over the pool modules = 1 by construction`;
   `transfer_spread: spread of r reported with its expected value; gate decided`. Checkpoint:
   `python SS/snake_scan.py transfer --rt 2p24GeV_H2` prints the spread and writes the prior.
5. Cosmic group check (P14.5, parallel to 2-4). Delivers `SS/out/mip_23744.tsv` (group, mpv, err, n_hits, n_modules)
   and the LMS events of 23744 into slice 3's input. Units: `tracks`, `mip_peaks`. Tests:
   `tracks_rate: 30 k events -> 150-300 vertical tracks`;
   `mip_mpv: interior hits -> peak-sample median 25-40 ADC (31.5 seen) and integral MPV 150-260 ADC, consistent across the 3 files within 2 SE`;
   `mip_groups: every group has >= 20 modules and >= 300 hits`. Checkpoint:
   `python SS/cosmic_mip.py EVIO/prad_023744.evio.00000 --map SS/hycal_map_til_run023849.json`.
6. Benchmark hook and comparison (P14.6; only if slice 4 passes the gate and P1b has landed). Delivers the
   `GAIN_PRIOR_FILE` hook, `output/energy_resolution/variants/P14/prior/` for 2.24 and 3.49 GeV, their
   `delta.tsv`, and the result section of `P14_snake_scan.md`. Unit: BM hook (one file). Tests:
   `hook_off: default report of 2.24 GeV -> byte-identical to the pre-hook default report (V3/v3_defaults if P1b changed no default; regions.tsv, gains.tsv)`;
   `hook_pool_inert: variant report -> pool rows of regions.tsv byte-identical to the control`;
   `hook_prior_applied: ring/gains.tsv O/T rows with a prior row -> g within sigma_g of the prior for modules
   with prior_frac > 0.5`;
   `delta_ot: delta.tsv O and T rows carry res_cb, tail_hi3, eff_cut deltas with SE`. Checkpoint: the two
   `compare_reports.py` calls of D7.

## Verification

After slice 6: the D3 verdict on O and T at 2.24 GeV with the 3.49 GeV replication, pool rows unchanged,
guards within bounds; `P14_snake_scan.md` lists the dwell coverage, s_i and L ratios per module, the pool spread,
the cosmic group check and the verdict; a default `report` of 2.24 GeV still reproduces the pre-hook default report (V3 when P1b changed no default) byte for byte.

## Risks

- Only about 25 scan-epoch LMS events locally: per-module L_scan error of order 0.5-1 % if the normalised
  event-to-event scatter is 2-5 %; 23744's 200 events help only if the HV was unchanged between Apr 19 and 20
  (D4 checks the agreement). The whole scan on the farm (P14.7) has about 250 per run.
- The HV change between the scan and production may change the PMT's LMS-to-shower response ratio
  (different light paths, PMT non-linearity at 370 ADC vs shower pulses up to 2 GeV); the pool validation (D5)
  measures the net effect and its spread, which is why the gate sits there.
- The edge model: a spectrum shape (bremsstrahlung times trigger acceptance) with a slope under the edge and
  pile-up above it; the window spread is recorded as the systematic; identical edges on 4 modules [W2] show the
  equalisation worked to a few per cent, so s_i is expected near one constant.
- Transporter sign/offset: if `ptrans` is not the HyCal position in beam coordinates, D2 falls back to the
  occupancy argmax alone (already sufficient for the 9 files).
- 7 pool modules in one run is a thin validation; the H group (12 modules) is the second group, and P14.7
  extends it.
- Concurrency: BM is edited by P1b; slice 6 waits (Scope Notes). `_common.setup_pipeline` and `lms_gain_variation.py`
  are read-only dependencies.
- Run 23744 file 00060 was still a partial download at planning time; slice 5 lists the files it used.

## Scope Notes

In: packages P14.1-P14.6 as written; P14.7 only after the user's answer to Q1. Out: any benchmark default
change; applying constants to the pool; absolute-scale transfer; NL or position dependence from the scan; the
production ring consumer; `pileup_rates.py` and `exotic_study/` (untouched); git (the user commits). User
decisions executor must not reopen: LMS-ratio transfer only, absolute scale with the pool, beam energy 2108.85
MeV from RCDB, no scan log (path from data), HV/LMS segments from the data (parent Log 2026-11-14 lines); the
comparison protocol D3 and U2 of the parent; memory goal entry
`.supervisor/memory/goals/2026-11-12-use-snake-scan-and-cosmic-data.md`. BM edits of slice 6 are made only after
the P1b executor has landed or the user says so; until then P14.6 Status is blocked.

## Assumptions

- The transporter EPICS `ptrans_x/y` are in mm in HyCal coordinates up to sign and offset; checked by the affine
  fit of slice 1 against the occupancy argmax.
- The same `WaveAnalyzer` settings (`DB/daq_config.json`) apply to the scan files as to production; checked by the
  pedestal rms and peak counts in slice 3.
- The LMS light reaches every PbWO4 module in the scan epoch as in production (1151 PbWO4 pulses per LMS event
  seen in both [W3]); checked by the count per event in slice 3.
- 23744 is at the scan HV; checked by the D4 agreement test.

## Work Packages

### P14.1: map and beam path
- Goal: derived map and dwell table for run 23701; done when slice 1 tests pass and the dwell table lists the
  11 modules of the 9 files with their event ranges
- Owns: `SS/snake_scan.py` (`derive_map`, `dwell_table`), `SS/hycal_map_til_run023849.json`, `SS/out/dwell_*.tsv`
- Deps: none
- Stages: execute
- Status: done
- Artifacts: `SS/out/dwell_23701.tsv`
- Cost: 0.3 agent-day; 9 files x 330 s (6 min wall at 9 lanes)

### P14.2: endpoint edge per module
- Goal: `s_i` with errors for every dwell module; done when slice 2 tests pass on the local files (and the 16
  extra files if Q2 is approved)
- Owns: `SS/snake_scan.py` (`edge_pass`, `fit_edge`), `SS/out/edge_*.tsv`, the per-module histogram npz
- Deps: P14.1
- Stages: execute
- Status: done
- Artifacts: `SS/out/edge_23701.tsv`
- Cost: 0.6 agent-day; 9-25 files x 330 s (6-20 min wall at 9 lanes)

### P14.3: LMS, alpha and pedestal series
- Goal: `L_scan,i`, `L_prod,i,rt` for the 5 run types, and the segment table; done when slice 3 tests pass
- Owns: `SS/snake_scan.py` (`lms_pass`, `segments`), `SS/out/lms_*.tsv`, `SS/out/segments_*.tsv`
- Deps: P14.1 (map); P14.5's 23744 LMS extraction feeds it but is not required
- Stages: execute
- Status: done
- Artifacts: `SS/out/lms_scan.tsv`, `SS/out/lms_prod_<rt>.tsv`
- Cost: 0.5 agent-day; about 130 production splits + 12 scan/cosmic files, 5-10 CPU-h (about 1 h wall at 8 lanes)

### P14.4: transfer, pool validation, prior tables
- Goal: `g_T,i`, the pool spread and the go/no-go for P14.6; done when `P14_snake_scan.md` states the spread and
  the gate decision and `prior_<rt>.tsv` exist for 2.24 and 3.49 GeV
- Owns: `SS/snake_scan.py` (`transfer`), `SS/out/transfer_*.tsv`, `SS/out/prior_*.tsv`, `ER/improvement/P14_snake_scan.md/.tsv`
- Deps: P14.2, P14.3
- Stages: execute, write
- Status: done
- Artifacts: `ER/improvement/P14_snake_scan.md`
- Cost: 0.3 agent-day; minutes of compute

### P14.5: cosmic MIP groups (run 23744)
- Goal: vertical-track MIP peaks per region group and the 23744 LMS events; done when slice 5 tests pass and the
  group ratios are in `P14_snake_scan.md`
- Owns: `SS/cosmic_mip.py`, `SS/out/mip_*.tsv`, `SS/out/lms_23744.npz`
- Deps: P14.1 (map)
- Stages: execute
- Status: done
- Artifacts: `SS/out/mip_23744.tsv`
- Cost: 0.5 agent-day; 3 files, about 1 h CPU

### P14.6: benchmark prior hook and D3 comparison
- Goal: `GAIN_PRIOR_FILE` hook and the O/T comparison against v3 at 2.24 and 3.49 GeV; done when slice 6 tests
  pass and the verdict is in `P14_snake_scan.md`
- Owns: BM (the hook only, after P1b lands), `output/energy_resolution/variants/P14/`
- Deps: P14.4 (gate passed)
- Stages: execute, write
- Status: blocked: execute: gate not met (4 pool modules, spread 2.29/3.03 %); user decision
- Artifacts: `output/energy_resolution/variants/P14/prior/delta.tsv`, `ER/improvement/P14_snake_scan.md`
- Cost: 0.5 agent-day; 2 x 45 min + 2 x 20 min reports

### P14.7: full scan on the farm (Q1)
- Goal: edge histograms and LMS/alpha/pedestal series for runs 23699-23723 and the whole of 23744, stored as
  small npz/tsv per run, giving O/T constants for all quadrants and per-module MIP peaks; done when
  `transfer_<rt>.tsv` covers every O and T module with a clean dwell
- Owns: a farm job script under `SS/farm/`, `SS/out/<run>/`
- Deps: P14.4 (gate passed), user answer to Q1
- Stages: plan, execute
- Status: blocked: user: submit farm pass on ifarm (package prepared)
- Artifacts: ->
- Cost: about 160 CPU-h on the farm for every 3rd file (about 570 files), 0.5 agent-day

## Log

- 2026-11-14 decision: plan approved (user); Q1 yes after the gate; Q2 first pass on the 9 files; Q3 unknown.
- 2026-11-14 P14.1: derived map (5 daq changes plus LMSP added), 9 dwells and 8 ramp rows on the 9 local 23701 files, `ptrans` is a set point that maps exactly onto module centres, 4/4 slice-1 tests pass; deferred: none.
- 2026-11-14 P14.2: `edge_pass.py` and the 9-module edge table done (s = 0.96-1.11, errors 0.3-1.0%); edge at 15.4k ADC not 18.5k, pile-up rejection 29-54%, all 4 slice-2 tests fail on plan criteria (convention, curvature, by-design); deferred: all-event edge cross-check.
- 2026-11-14 P14.3: LMS/alpha/pedestal pass done (25 scan events, L_scan 1151 modules at 0.82%; 131 production splits, L_prod 0.05-0.12%); no segment steps; 23744 sits 2.5% below 23701, so the primary L_scan is 23701-only and the merged variant is a sensitivity; 4/4 tests pass; deferred: none.
- 2026-11-14 P14.4: pool spread of r (4 modules W552/W557/W590/W599) is 2.29% at 2.24 GeV and 3.03% at 3.49 GeV, so the gate fails at both (merged variant 1.57% / 3.02%); the transfer does not beat the constant pool g (1.1-1.3%); priors exist for only W547 (O) and W580 (T); deferred: P14.6 and P14.7 (gate not met, user decision).
- 2026-11-14 P14.5: vertical-track MIP on 23744 done (4421 tracks, 31,080 interior hits; group MPV 192-203 ADC, 22.8 MeV at 0.114 and 26.6 MeV via W); `lms_23744.npz` holds 202 LMS events and 30,699 alpha events; `mip_mpv` fails at 2.08 SE between files; deferred: per-module MIP constants (need the whole run, P14.7).
- 2026-11-14 P14.7: farm package prepared in `SS/farm/` (single decode loop per scan file, merge calls the wave-1 functions unchanged, `test_farm.py` 10/10 against wave 1, `--lms-runs 23701` reproduces the primary transfer); about 500 scan jobs rather than 190; swif2/sbatch syntax, accounts and the ifarm build untested; deferred: the farm run itself and the transfer on the pooled scan (user submits).
- 2026-11-14 P14.2b: all-event edge check done (`P14_edge_pileup_check.md/.tsv`); in-time pile-up settings shift s by 0.05–0.46 % (pool about +0.2 to +0.3 %), `none_sum` unusable, pool spread 2.28–2.34 % for all in-time settings, so rejection explains none of the gate failure; keep wave 1's rule for the farm pass; deferred: pulses inside `peak_pileup_gap` are not bounded.
- 2026-11-14 amend: slice criteria corrected (see `## Amendments` below); 23744/23701 LMS ratio corrected 0.975 → 0.9935.

## Todo

- [x] P14.1 derive map; dwell table of the 9 files; EPICS affine map; slice 1 tests
- [x] Q2: superseded — first pass on 9 files (user)
- [x] P14.2 edge pass and fit; window-spread systematic; neighbour second pass; slice 2 tests
- [x] P14.3 LMS/alpha/pedestal pass on 23701, 23744 and the local production splits; segments; slice 3 tests
- [x] P14.5 vertical-track MIP groups on 23744; slice 5 tests
- [x] P14.4 transfer tables, pool spread, gate decision, `P14_snake_scan.md`; slice 4 tests
- [ ] P14.6 hook (after P1b lands), hook_off byte-identity, variant reports, `delta.tsv`, verdict; slice 6 tests
- [ ] Q1: P14.7 farm pass decision (user); then its plan

## Amendments (2026-11-14, from P14.7/P14.2b)
- Slice 2 `edge_consistency`: compare the ADC edge, not s. Edge is at ~15.4 k ADC (3×3 integral convention), not 18.5 k. ADC edges agree within 4.3 % across the 9 modules at scan HV. s = E_edge/2108.85 spans 0.96–1.11 by design (W is the production calibration) → drop the "within 1 %" test on s.
- Slice 2 pile-up window: replace "about 10 % rejected" with the measured 29–54 % for the any-member `len(peaks) >= 2` rule (the ~10 % figure corresponds to `loose2`, 4–9 %). Criterion: in-time settings agree within 0.5 % in s per module (measured 0.05–0.46 %).
- Slice 2 curvature test: fails by design on these edges; replace with the window-spread systematic already recorded by `fit_edge`.
- Slice 5 MIP split tolerance: 3 SE (was failing at 2.08 SE). 23744/23701 LMS ratio is 0.9935 (measured), not 0.975.
- Slice 4 gate: judge on the 4-module pool spread; note that the pool+H spread differs (5.4–5.6 %). P14.7 adds modules; the gate is only assessable after the farm pass.
