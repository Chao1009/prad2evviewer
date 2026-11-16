# Plan: exotic waveforms in the HyCal pile-up fit — a module-complete check, then the fit — 2026-11-12

Working dir `analysis/pyscripts` (branch `waveform_jz`); paths are relative to it. PR =
`waveform_analysis/pileup_analysis/pileup_rates.py`; S = `waveform_analysis/pileup_analysis/pileup_deconvolution_study.py`;
G = `waveform_analysis/pileup_analysis/pileup_generator_stage1.py`; BM = `waveform_analysis/energy_resolution_benchmark.py`;
ED = `waveform_analysis/pileup_analysis/exotic_events/` (event displays of the two events, done: `summary.txt`, one
PDF per event); XS = `waveform_analysis/pileup_analysis/exotic_study/` (new: this plan's scripts and small outputs);
PP = `notes/plan_energy_resolution_improvement_20261112.md` (program plan; its D1-D13, P0-P14 are cited by name).
Environment: `source setup_env.sh`, `MPLBACKEND=Agg`.

## Open Questions

None blocks the plan. Each has a default that holds unless changed at approval.

1. Rescan eight runs, two per run type, with the new per-channel variables into a new directory (slice 3; 8 CPUs,
   about 1.4 h in parallel)? Recommended: yes, runs 24590, 24603, 25320, 25228, 24698, 24911, 26138, 25800: two independent
   runs per type are what confirms a flagged module (X7). Alternative: four runs (one per type), with the two
   halves of each run as the independent samples.
2. Benchmark side: BM gets an acceptance or refit change only as one arm of P8's first test (`Amend parent:` in
   X12); the per-hit class word in BM is proposed at the P5 review, not added to P8. Recommended: yes.
3. Is the HV archive (`/data/prad2/hv_data`, read by `hv_channel_trace`, `../../python/bind_dec.cpp:1146`; not on
   this machine) reachable, and who are the contacts for base history and the fADC250 boards? Recommended: answer
   when slice 6 delivers the list of anomalous modules (X8).
4. Lead glass: its channels get table rows and maps where the scans hold kept channels, reported apart from PbWO4
   and without a width or tail measure (the template banks hold the 1152 PbWO4 modules only). Recommended: yes.

## Goal

Check every HyCal module the same way for waveforms the template cannot describe: how the algorithm identifies
them, how often each module shows each class in each run type, and what they are. W121 (run 24590, event 4084847)
and W493 (run 24590, event 2069998) are two validation cases: the procedure must flag both events with limits set
without them. Deliverables: a class word per fitted channel from measured quality variables; one table row per
module and run type with the class rates and shape measures, HyCal maps of each column, ranked lists of anomalous
modules from an outlier rule fixed in advance (X7), example pages for the top modules; a stability check across
runs and between H2 and Ta; hot-module, template-candidate and bad-channel lists; a ranked, tested choice of fit
change for the wide-pulse class. Non-goals: edits to BM, `mf_pulse.py`, the template banks, `../../database/`, the
C++ `WaveAnalyzer`, ED; the Gaussian noise budget (PP D9); crosstalk and saturation corrections (dropped by PP);
new P3 or P8 work. Expected change: about 70 lines added to PR (five variable functions, record fields, two
counters, sample reservoirs; `fit()` untouched), two new tools of about 400 lines each and one run script in XS,
one results file. No change to any fit result.

## Approach

Measure first, in two phases, every quantity for all modules. Phase A (slices 1-2, about one day, no scan): the 13
existing scan files hold the fit records (k, amplitudes, times, χ², baseline counters) of every kept channel. That
gives the all-module table, maps, outlier lists and run-to-run stability for the quantities those records hold
(extras by (dt, r) locus, the accidental proxy, k, χ², `far`, bad rows) and a first re-ranked hot list: the
detector-wide picture at the cheap stop point. Its numbers are reported before phase B starts. Phase B: slice 3
adds five quality variables and the stored `far` flag where samples and fit coexist (PR `flush()`), rescans two
runs per run type, and sets limits from the one-pulse channels; slices 4-5 build a truth bench and test the fit
options in order; slice 6 completes the table with classes C1-C6, widths and baseline episodes, and runs the
tests that say what each class is; slice 7 hands the lists to P3/P8, the bank work and the user. Nothing in the fit
or the templates changes inside this plan; a change is a recommendation.

## Decisions

**X1. Corrections to the request and the brief (code facts).**
- The 5-pulse fit is not `mf_pulse.py`: PR, S and BM each carry a copy of class `MF` (PR:210; `KMAX = 6`, PR:128).
- The fit variance already holds the pulse-shape scatter (diagonal `shape_sigma` term, PR:237). An extra pulse is
  kept when A ≥ max(5σ_A, 5 ADC), A ≥ 0.02 × every pulse within 96 ns (`R_MODEL`, PR:132) and Δχ² > 2 ln N_w
  (PR:140, :514). What the variance lacks is the sample-to-sample correlation of a shape deviation.
- χ² is taken inside the pulse windows only (PR:671) and grows with amplitude: χ²/dof > 10 in 0.1 % of channels at
  0.73 GeV, 0.5 % at 2.24 GeV, 2.3-3.2 % at 3.49 GeV. χ² alone cannot select a class.
- A 30-sample baseline more than 3 σ_ref from the run reference b_ref is replaced by b_ref (`far`, PR:614), on the
  assumption that a pulse sits in the pedestal window; `far` is counted per module, not stored per channel.
- Both events are in `pileup_rates_out/examples.pdf` (pages 6, 7), not in `hot_modules.pdf`. W121 is hot rank 1
  only at 3.49 GeV (19 of 268); W493 is rank 473 at 0.73 GeV and hot only in Ta (7.1 %).

**X2. The two validation cases and the detector-wide numbers behind them** (`output/pileup_rates/scan/run_24590.npz`,
`ED/summary.txt`).

| | W121 ev 4084847 | W493 ev 2069998 |
|---|---|---|
| fit | k = 4; A = 76, 89, 73, 55 ADC at 195, 241, 275, 311 ns; χ²/dof 128 | k = 5; A = 120, 41, 23, 16, 8 ADC at 154, 164, 179, 197, 230 ns; χ²/dof 1.15 |
| first 30 samples | mean 203.9 = b_ref + 24 ADC (25 σ_ref), rms 3.9 = 4.0 σ_ref, slope −2.8 ADC per 100 ns | mean 169.9, rms 0.76 (clean) |
| shape | hump from 180 to 380 ns, +64 above its own baseline, FWHM 111 ns (about 140 ns above b_ref; template 23.9 ns) | flat top over 4 samples, FWHM 36 ns (template 22.7 ns, ratio 1.6); +17 % of peak at 68 ns after it |
| consequence | `far`: the fit starts from b_ref, so the offset is fitted as pulses; E_main 83 MeV is not shower energy | extras in a decreasing series along the tail; raw peak 133 ADC, main pulse fitted 120 ADC (122 MeV) |
| event (ED) | 5 channels read out; no neighbour and no board mate among them; one-module cluster of 79 MeV beside the physics cluster (seed W564, 750 MeV, E/E_beam 1.03) | 8 channels; W493 is 26 % of a 4-module cluster (seed W494 269 MeV, E_cluster 532 MeV, E/E_beam 0.73); W494, W460, W459 have one pulse each (151.0, 150.2, 148.4 ns), nothing at the extras' times |

- W121 is not an outlier module by the existing counts: `far` occupancy 1.7e-6 to 1.1e-4 per event over 13 runs, up
  to 3.4 times the PbWO4 median (rank 220-670 of 1152; `far` is dominated by real early pulses); in run 24590 it is
  a median module in hits (930 kept channels, median 942) with a χ²/dof > 10 fraction 2-3 times the per-module
  median. In its EVIO split 4 of 232 read-outs have a noisy pre-window, and W86, W87, W154, W155 show 1-4 such
  events each, in other events: the baseline episode is a class of many modules, seen here in one instance.
- Tail-locus extras are detector-wide and not accidentals. Run 24590: 89 % of W493's 10 417 extras lie 8-48 ns
  after the main pulse with median r 0.056 (the "locus"; `waveform_analysis/opens.md` item 9), 4 % before it.
  Channels with a locus extra: all modules 1.5-2.5 % in H2 (W493 2.5-3.7 %), 2.7-4.2 % in Ta (W493 3.5-6.4 %);
  channels with an extra before the main pulse (dt < −8 ns, the proxy for accidentals): 0.13-0.5 % in H2, 10-15 %
  in Ta. Accidentals change about 80-fold, the locus fraction by at most 3: the fit splits the tail of 2-3 % of
  single pulses, and the W493 event is the multi-extra end of that population.
- W493 event: if the four extras belong to the main pulse, W493 carries up to 91 MeV more than the 122 MeV used.
  Its template is built from 1000 elastic-pool pulses at 306-584 ADC; a 120 ADC non-seed member is outside it.
- Run 24937 is anomalous for all modules (k ≥ 4 in 1.0 % of channels, 94 k events); the survey treats it apart.
- Correction (phase A, `XS/results.md`): the Ta/H2 ratio of the before-main fraction f_acc is 33 in rings 1-2 (beam-hole rings) but 1.7 in rings 11-16 (median over modules 3.4), not 20-115 detector-wide; the "80-fold" above holds for the W493 ring only. Slices 3-6 use the ratio per ring group (X7 label rule, slice 6 `rate` test), not these numbers; S3's limits do not depend on it.

**X3. Classes.** One bit per class C2-C6, several may be set; the primary class follows the order saturated, C5,
C4, C3, C6, C2, C1, C0. Each rule is one boolean in the X4 variables; a limit written L is an output of slice 3
(X5), slice 4 (L_r) or slice 6 (L_4, L_5). "Occupancy" = flagged channels per event among the rows PR keeps: not
`bad_waveform` (PR:652), raw maximum ≥ 8 ADC above the pedestal, fitted main ≥ 10 ADC. It measures a duty cycle
only for anomalies that pass these cuts and the DAQ readout threshold, which is not established (PP U6, owned by
P12; assumed to be the firmware 10 ADC). Two counters give the flagged rows these cuts remove (X4).

| class | origin | rule in one waveform | over time / in the cluster | differs from |
|---|---|---|---|---|
| C0 clean single | one shower | k = 1, no bit set | neighbours share one time | reference |
| C1 resolved pile-up | second particle, any time in the 400 ns window | k ≥ 2, no bit set: an extra before the main pulse, later than 48 ns, or with r > L_r | rate ∝ beam rate, flat in dt; members able to show it do (C_t), with shower-like sharing | C2: extras are not ordered along the tail |
| C2 wide / shape mismatch | open (X9): shape depends on amplitude or role, slow tail, module shoulder, pair closer than 8 ns | not C3, and w50 > L_w or tail > L_t or every extra on the locus (8 ≤ dt ≤ 48 ns, r ≤ L_r) | fraction per kept channel independent of beam rate (X2); same locus event after event | C1: locus; C3: baseline clean |
| C3 noisy baseline with slow hump | open (X8): PMT/base/HV episode | `far`, z_n > L_n, and no fitted pulse overlapping the 30-sample pedestal window | one module, episodic; alone it forms an isolated one-module cluster (79 MeV in the W121 event) | C5: board mates quiet; C6: no preceding pulse |
| C4 persistent bad channel | ringing or broken channel (W335 at 0.73 GeV, `waveform_analysis/opens.md` item 11) | C3 rule | C3 occupancy of the module above L_4 over a whole run or split | C3: duty cycle |
| C5 coherent pickup | board, crate or array-wide electronics | C3 rule | ≥ L_5 channels with the C3 bit on one fADC board in one event; ≥ 300 kept channels is already vetoed (PR:101) | C3: several channels at once |
| C6 afterpulse / baseline recovery | large earlier pulse in the same channel | not C3, and slope < −L_s or an extra at dt > 100 ns with r < L_r | follows a large pulse in the same channel; rate ∝ hit rate and amplitude | C3: slope and history |
| saturated; glitch or undershoot; crosstalk spike | amplitude above 4095; readout; neighbour channel | `SAT_NCLIP` (PR:107); `bad_waveform` (PR:652); crosstalk not flagged | — | existing handling, unchanged |

**X4. Quality variables.** Per kept channel, computed in PR `flush()` (the only place where samples and fit
coexist) by module-level functions that the tools import, stored as new `ch_*` arrays next to the existing k, a,
t, χ², nclip. `fit()` is not edited.
- `far` (flag): |b_30 − b_ref| > 3 σ_ref (returned by `fit()`, now stored). z_b = (b_30 − b_ref)/σ_ref. W121: +25.
- z_n = rms of the first 30 samples about their median, divided by σ_ref. W121: 4.0; W493: 0.85.
- slope = straight-line slope of the first 30 samples, ADC per 100 ns. W121: −2.8.
- w50 = full width at half maximum of the waveform around its largest sample, above the baseline the fit used
  (b_ref when `far`), divided by the same width of the template at the fitted phase. W121: about 5; W493: 1.6.
- tail = mean over [t_main + 40, t_main + 120 ns], clipped at the record end, of (waveform − baseline − main
  component of the existing fit, `comp`, PR:680)/A_main. W493: about +0.1.
- Per event: timestamp (`fe.info`, PR:876-888).
- Two counters per module, from z_b and z_n evaluated on all rows at the top of `flush()`: rows with `far` and
  z_n > 3 among the `bad_waveform` rows, and among the rows failing the 8 or 10 ADC cut. No arrays for rejected rows.
- Samples for the refits of slice 5: kept channels with `far`, k ≥ 4, χ²/dof > 10 or w50 > 1.3 (four loose fixed
  conditions, not class limits), 20 per module and condition, through the existing `reservoir()` (PR:698) with a
  generator of their own, so the present `ex` and `pm` reservoirs keep their content (PR:746, :795).

Derived by the survey tool, not stored (limits may move): extras (dt_k, r_k = A_k/A_main); the class word; board
occupancy N_b = kept channels of the same fADC board in the event (`hycal_map.json` daq fields) and how many carry
the C3 bit; cluster agreement C_t = members that show an extra within 4 ns of the same time, divided by members
large enough to show it (A_main·r ≥ 5 ADC), from the stored `ev` and `midx`; per module and class the occupancy and
its dependence on split and event block (Fano factor). Readout is zero-suppressed (5-8 channels per event), so a
waveform correlation between board mates cannot be formed in one event; coincidence counts over events replace it.

**X5. Limits come from a clean population, not from the two events.** L_n, L_w, L_t, L_s (on −slope): the value 0.1 % of
the k = 1, not-`far` channels of the rescan exceed, per amplitude band (20-50, 50-200, 200-1000, > 1000 ADC), per
module where it has ≥ 10 000 such channels in the band and pooled over its ring (ring index of
`waveform_analysis/template_extract.py:1219`) otherwise. Holdout: limits from even event numbers, applied to odd.
The template pools are not a reference: they start at 50 ADC and are selected against bumps and tails. L_r is set
on the truth bench (slice 4): true pairs flagged C2 in at most 2.5 % of the cell dt ≥ 12 ns, r ≥ 0.1, the cell of
the P8 guard. The two events are test cases, not tuning data.

**X6. What the chain does with each class.** While measuring: flag only; members failing a χ² cut carry real
energy (`notes/PbWO4_energy_resolution_20261030.md:103-107`) and PP dropped χ² cuts. Proposed actions, each to be
confirmed by its test: C1 as now. C2: locus extras are not counted as pile-up; amplitude estimator chosen on
truth (X10 O2). C3: channel flagged and not counted in pile-up rates; a one-module C3 cluster is flagged as fake;
inside a cluster the channel stays in the sum until the P8 arm shows the effect on the D1 figures. C4: per-run
bad-channel list as a proposed table (`../../database/` is not edited; the W335 decision stays with the user).
C5: event flag. C6: extras labelled, not counted as pile-up. No class is down-weighted.

**X7. Module-complete outputs, outlier rule, stability, and what they feed.**
- Table: `XS/modules_<run type>.tsv`, one row per module (1152 PbWO4; lead glass per Open Question 4), and
  `XS/modules_by_run.tsv` with one row per module and run. Tested columns: kept channels n; fractions with k ≥ 2,
  3, 4; locus-extra fraction f_loc; before-main extra fraction f_acc (accidental proxy); χ²/dof > 10 fraction per
  amplitude band; `far` per fitted row (counter `mod_n_base_ref`, counted before the keep cut, PR:771-775);
  `bad_waveform` rows per `digit_count`. Descriptive columns, not tested, from the bank of the run type:
  `fallback`, `single`, `pool_general`, `sigma_borrowed`, `n_pool`, `a_band`, `shape_sigma` at the peak,
  `fall_med_ns`. Phase B adds tested columns: fraction of each class C1-C6; medians of w50 and tail of k = 1
  channels per band; C3 episode rate and its Fano factor over splits; rejected-row counters.
- Maps: one HyCal map per column and run type in `XS/exotic_survey.pdf` (map code of PR's report as pattern).
- Outlier rule, fixed here, per tested column q and run. Reference group: modules of the same type (PbWO4, lead
  glass), ring index (`waveform_analysis/template_extract.py:1219`; rings with fewer than 30 modules merged with
  the next), run type and amplitude band; q_g and s_g = 1.4826·MAD are its median and spread. Fractions:
  z = (q_m − q_g)/√(s_g² + q_g(1 − q_g)/n_m), binomial error at the group median; where the expected count
  n_m·q_g is below 10, z is the Gaussian equivalent of the Poisson probability of at least the observed count.
  Medians (w50, tail): σ_stat = 1.253·rms/√n_m, tested two-sided. A module is flagged at |z| > 5 with n_m ≥ 100
  (`MIN_N_MOD`, PR:115). Printed per list: number of tests, the Gaussian expectation (tests × 2.9e-7 per side), and
  the null count: modules flagged in the even or the odd events of a run but not in both (two-sided columns: also
  the count at z < −5).
- Stability: a module is confirmed when flagged in at least two runs and in more than half of the runs of its run
  type (run 24937 not counted). At 2.24 GeV, the only energy with both targets, a confirmed module is labelled
  module-intrinsic when |q_Ta − q_H2| < 3σ of the difference in the same amplitude band, rate-dependent when
  q_Ta/q_H2 > 3 (f_acc differs 20-115-fold, X2), else undecided. Episodic classes also get the split history.
- Ranked lists: `XS/anomalous_<column or class>.tsv` (module, z per run, runs flagged, label), confirmed modules first.
- Examples: for the top 5 confirmed modules of a list, 2 events each (largest χ²/dof for the χ² list, most locus
  extras for f_loc and k lists, the flagged channels of largest z_b for class lists; ties by lowest event number),
  the survey writes `XS/examples/cmds.sh` calling `ED/event_display.py --run --event --module --out` (waveform
  with fit, neighbours, cluster map); ED is called, not edited. Event numbers are `ev_number[ch_ev]`. `far` and
  bad-row lists rest on counters and get examples in phase B only.
- Feeds. Hot modules: two lists replace the single ranking (PR:1314-1318, cut pile-up per kept channel, top 8):
  cut pile-up rate over C0 + C1 channels only (phase A: without locus extras), ranked by z, so a module with 5 of
  181 no longer outranks one with 2000 of 230 000; and the noisy-channel list (C3 + C4). Template bank:
  `XS/template_candidates.tsv`, confirmed modules whose w50 or tail median, or f_loc, is an outlier, with the
  dependence of w50 on amplitude (flat: rebuild; sloped: amplitude bands), handed to the bank work. Energy
  resolution: slice 7 (P3/P8, per-run bad-channel table). PR's report changes only on a slice 7 recommendation.

**X8. Baseline episodes (C3-C5): hypotheses, and what decides.** Slice 6 applies the tests to every module on
the C3 list; W121 is the worked case and its event display fills the last column.

| hypothesis | confirmed by | refuted by | observed (ED) |
|---|---|---|---|
| A. single-channel episode (PMT, base, HV contact) | flagged events have no flagged board mate; occupancy equal in all run types; episodes clustered in time (Fano ≫ 1) | board mates flagged in the same events | consistent: no board mate read out; noisy fraction 1 of 11 within ±5000 events against 2 of 220 outside |
| B. coherent pickup on the module's fADC board, large enough to be read out | board mates (for W121, crate 3 slot 4: W154, W88, W120, W86, W87, W53, W54, W20, W52, W18, W19, G166, G136, G106, G76) read out and flagged in the same events far above the product of their single rates | no coincidences | disfavoured if the readout threshold is the firmware 10 ADC (X3): a common 60 ADC hump would be read out on every channel, none is; W154's noisy events are other events |
| C. baseline recovery or afterpulse train | large pulse of the module in a recorded event shortly before (timestamps); occupancy ∝ the module's hit rate (Ta ≫ H2) | flat raised pre-window; occupancy independent of hit rate | open: pre-window falls by 2.8 ADC per 100 ns on a +24 ADC offset |
| D. class common to many modules | the module's C3 z below 5 (X7); other modules show the class | module confirmed as a C3 outlier (X7) | supported locally: W86, W87, W154, W155 show it, four other neighbours do not |

B below the readout threshold is the known common-mode noise (60 % of the baseline variance within a board,
`output/noise_diag_20261028/report.md`), not this class. If the coincidence counts of slice 6 are too small, the
complete boards of LMS events are the fallback sample (`noise_diag.py`). Hardware questions prepared for the user:
HV and base history of the confirmed C3/C4 modules (Open Question 3); shared HV group or cable bundle; whether a
slow positive excursion of 60-80 ADC over 150 ns on a raised, noisy baseline matches a known PMT or divider fault.

**X9. Tail-locus extras and wide pulses (C2): hypotheses, and what decides.** Slice 6 applies the tests per
module and ring; W493 is the worked case.

| hypothesis | confirmed by | refuted by | observed (ED) |
|---|---|---|---|
| 1. real later pulses or an unresolved pair shared by the shower | neighbours show the same time structure with shower-like sharing; C2 fraction rises with the accidental proxy | neighbours single; C2 fraction flat against the accidental proxy | refuted for this event: three neighbours, one pulse each, nothing at the extras' times; population: locus fraction flat while accidentals change 80-fold (X2) |
| 2. shape depends on amplitude or on the module's role (non-seed member, below 300 ADC) | w50 and tail of k = 1 and pool pulses rise as amplitude falls, or differ between seed and non-seed at equal amplitude | w50 flat in both | open; the event is in the untested corner (120 ADC, non-seed, 3-6 ns late) |
| 3. template too narrow for all pulses | median w50 of pool pulses > 1.02 at every amplitude; `waveform_analysis/opens.md` item 18 (raw rise 6.6 vs 5.35 ns) | median w50 = 1.00 ± 0.02 | disfavoured: χ² median 0.85; 97 % of W493 channels have no locus extra |
| 4. rare wide class of this module or ring (shoulder modules W493, W661, W663, W664, `opens.md` item 17; direct PMT hit; afterpulse) | w50 distribution has a separate population near 1.6; module-specific | population continuous and common to all modules | open |

**X10. Fit options for the wide-pulse case, in this order.** Data: real singles from the pools
(`waveform_analysis/pulse_templates/0p73GeV_H2/pulses/run_24590.npz`, `.../2p24GeV_H2/pulses/run_25320.npz`),
generator pairs and r = 0 controls, distorted singles (stretch 1.1-1.6, added tail 5-15 %), the C2 samples of the
rescan. Metrics: F = fake extras per true single; ε = pair efficiency in the worst cell (dt ≥ 12 ns, r ≥ 0.1,
pooled as in `waveform_analysis/pileup_analysis/efficiency.tsv`, now 97.5 %); bias and spread of the main
amplitude. The D1 figures (res_cb, sep_pi, tail_hi3, tail_lo3) are judged only by P8 on BM.

| option | first test | success | kill |
|---|---|---|---|
| O0 no change, measure | F on real pool singles per amplitude band; (dt, r) locus of fake extras; F on distorted singles | numbers exist | — |
| O1 acceptance with existing constants, one at a time: `R_MODEL` 0.02 → 0.05, 0.1, 0.3; `KMAX` 6 → 3, 2; `BIC_PENALTY` 1 → 3; baseline rows `MIN_SEP_NS` 8 (present) and 12 (P3 variant) | same samples per setting | F halves, ε ≥ 97.5 %, amplitude spread of clean singles not worse by > 2 SE | no setting meets both |
| O2 classify, then refit | C2 channels refitted by `MF` with `KMAX` = 1 (amplitude and time free, baseline as in PR); estimators main, one-pulse, main + locus extras compared with truth (distorted single: the amplitude of the undistorted pulse of equal integral) | amplitude bias of distorted singles < 0.5 % and spread below the K-pulse fit's | bias or spread worse than O1 |
| O3 one shape parameter with a prior (time stretch s of the template, or the first principal component of pool residuals) | fit s on pool pulses: its distribution against amplitude; share of residual variance in the first component | F halves on C2, ε ≥ 97.5 %, clean-single spread not worse by > 2 SE | first component holds < 30 % of residual variance, or ε falls |
| O4 amplitude-banded templates | only if hypothesis 2 holds: request to the bank work, with the measured w50(A) | — | hypothesis 2 refuted |

O4 is last: per-module f(A) made σ/E worse in production (`notes/PbWO4_energy_resolution_summary.md:113-117`) and
banks are rebuilt outside this plan. Of O3 only the first test runs here; an O3 fit is a recommendation to P8.

**X11. Template jitter: no.** Jitter is already in the model twice: the alignment residual is in the mean shape,
the pulse-to-pulse spread is in the variance (X1). Measured per module
(`waveform_analysis/pulse_templates/timing/timing_summary.tsv`, medians over modules per run type): `jitter_50`
0.040-0.046 ns, rise jitter 0.14-0.17 ns, fall jitter 2.1-2.4 ns (W493 2.05 ns, W121 2.90 ns). Turning 22.7 ns
FWHM into 36 ns needs a Gaussian σ of 12 ns, about 6 times the largest of these, so the W493 event is another
class, not the tail of the jitter, and a smeared template would be wrong for the 97 % of W493 channels without a
locus extra. Deciding test (slice 5): fit pool singles and the C2 samples with the template convolved by σ_j = 0,
0.5, 1, 2, 4, 8, 12 ns. Adopt only if F on C2 halves while the pool's χ² median and amplitude bias move by less
than 2 SE; expected outcome: pool χ² and bias degrade from 1-2 ns on. The better way is X10: classify, then O1-O3.

**Not to do, and why.** Smear or widen the template (X11). Cut on χ² (energy loss; χ² scales with amplitude, X1).
Inflate the diagonal variance globally (hides real pile-up; the ε guard would fail). Mask W121 or W493 on single
events (W121's `far` occupancy is ordinary, X2). Tune limits on the two events (X5). Remove in-time energy of a
flagged channel before the P8 arm. Edit BM, banks or `mf_pulse.py` here.

**X12. Relation to the program plan: separate plan, three hand-offs.** This plan works at waveform level in the
pile-up study chain (PR, S, G), which no PP package owns; PP works at cluster level in BM. No BM switch is needed:
the bench imports PR and sets its module constants, as `ED/event_display.py` imports PR. P0: no dependency. P3:
unchanged; slice 7 gives P3/P5 the artefact locus in (dt, r), so a width or tail excess in a `SECONDARY_DUMP` cell
can be read as shape artefact or real pulse. P8: receives the candidate rule as one arm of its first test. The
per-hit class word and χ² in BM (`waveform_analysis/opens.md` item 6) are a proposal for the P5 review, not P8.

`Amend parent:` P8 — `Deps: P3, P5` unchanged. `Artifacts` adds section "For P8" of `XS/results.md`, read when P8
starts if it exists. First test adds one arm: this plan's candidate acceptance or refit rule beside the P3 winner.
Gate D7 and the 97.5 % guard are unchanged.

**X13. Compute.** One PR scan is one process per run, 0.5-1.3 h for the eight runs (stored `wall_s`). Slice 3 runs
the eight in parallel when no PP scan is running (PP D11, 18 CPUs shared). Everything else reads npz files.

## Outline

- PR (edit, about 70 lines): module-level functions for z_b, z_n, slope, w50, tail (beside `pedestal`, PR:578);
  `REC` (PR:110) gains these five fields and `far`; `flush()` (PR:761) fills them, fills the two rejected-row
  counters, records the timestamp per event, and keeps the sample reservoirs of X4; `load()` and `report()`
  ignore the new arrays. `fit()` (PR:663) is not edited; fit results stay bit-identical.
- `XS/exotic_survey.py` (new): `records <scan dir> <out dir>` (slice 2: table, maps, outlier lists, stability and
  example commands from the existing files), `classes <scan_q dir> <out dir>` (limits, class word, full table,
  lists, slice 6 tests), `selftest`. One code path builds table, z and lists for both modes.
- `XS/exotic_bench.py` (new): imports PR (`Bank`, `MF`, quality functions) and G's waveform builder; modes `pool`
  (real singles), `pseudo` (pairs, controls, distorted singles; writes L_r), `options` (O1-O3, X11 test). G is
  edited only if its builder cannot be imported, by one keyword argument with the present behaviour as default.
- `XS/run_scan_q.sh` (new): the scan line of `waveform_analysis/pileup_analysis/run_pileup_scan.sh` for the eight
  runs, output `output/pileup_rates/scan_q/run_<run>.npz` (new directory; the existing `scan/` is not overwritten).
- Outputs in XS (all new): `results.md` (one section per slice, verdict tables, "For P8", "For the P5 review"),
  `modules_<run type>.tsv`, `modules_by_run.tsv`, `anomalous_<column or class>.tsv`, `template_candidates.tsv`,
  `thresholds.tsv`, `locus_limits.tsv`, `hot_pileup.tsv`, `noisy_channels.tsv`, `bench_options.tsv`,
  `examples/`, `exotic_survey.pdf`, `exotic_bench.pdf`.

## Slices

**S1. Record the two validation cases** (no code). Units: `XS/results.md` (new): the X2 table and the ED columns
of X8 and X9, with their source lines in `ED/summary.txt` and the two cluster maps cited by PDF.
- Checkpoint: every observed cell carries a number or "open, slice 6".

**S2. All-module survey from the existing scan files** (no scan). Units: `XS/exotic_survey.py records`.
Delivers the phase A columns of X7 for every module and run, the maps, outlier lists, stability labels, first hot
lists and the example commands; the example pages are made for the k, f_loc and χ² lists (X7).
- table_complete: 13 scans -> 1152 PbWO4 rows per run type, no empty cell where n ≥ 100; lead-glass rows counted.
- w121_counts: `run_24590.npz` -> W121 kept 930, χ²/dof > 10: 7, k ≥ 3: 6, far 5.
- w493_multi: `run_24590.npz` -> W493 k ≥ 2 3.55 %, locus extra 3.21 %, extra before main 0.15 %, k ≥ 4 0.018 %.
- event_lookup: reservoir `pm1`, W121 ev 4084847 -> k 4, χ²/dof 127.6; W493 ev 2069998 -> k 5, χ²/dof 1.15.
- outlier_control: every list -> tests, flags, Gaussian expectation, even/odd null count, confirmed modules printed.
- scaling: median over modules, Ta over 2.24 GeV H2 -> f_acc ratio 20-115, f_loc ratio below 3; modules off this printed.
- hot_rerank: 3.49 GeV -> W121 printed as 7.1 ± 1.6 % with its z and rank; list without locus extras beside it.
- Checkpoint: `python3 waveform_analysis/pileup_analysis/exotic_study/exotic_survey.py records
  output/pileup_rates/scan waveform_analysis/pileup_analysis/exotic_study` prints the checks as passed and writes
  the tables, lists and `exotic_survey.pdf`. Phase A ends here: its numbers are reported before S3 starts.

**S3. Quality variables, rescan, limits.** Units: PR, `XS/run_scan_q.sh`, `exotic_survey.py classes|selftest`.
- w121_vars: stored samples of W121 ev 4084847 -> far, z_b 20-30, z_n 3-5, w50 > 3, slope −2 to −4 ADC per 100 ns.
- w493_vars: stored samples of W493 ev 2069998 -> not far, |z_b| < 2, z_n 0.6-1.2, w50 1.4-1.8, tail > 0.05.
- template_self: noiseless W493 template, A = 300 ADC, eight phases -> w50 1.00 ± 0.02, tail 0 ± 0.005.
- fit_unchanged: `PILEUP_MAX_EVENTS=200000` scan of run 24590 before and after the edit -> `ch_k`, `ch_a`, `ch_t`,
  `ch_chi2` and the `ex`, `pm` reservoirs identical.
- rejected_rows: the same scan -> both counters per module are ≤ `bad_count` and ≤ the number of cut failures.
- display_still_runs: the ED command for one event with `--out` in a scratch directory -> PDF written.
- threshold_holdout: limits from even event numbers of the run 24590 rescan applied to odd -> flagged fraction
  0.1 % ± 0.05 % per band, modules pooled.
- Checkpoint: `exotic_survey.py selftest` passes; then `sh waveform_analysis/pileup_analysis/exotic_study/run_scan_q.sh`
  and `exotic_survey.py classes output/pileup_rates/scan_q waveform_analysis/pileup_analysis/exotic_study` write
  `thresholds.tsv` and the occupancies of all classes but C2 (L_r comes from slice 4) for the eight runs.

**S4. Truth bench.** Units: `XS/exotic_bench.py pool|pseudo`.
- control_reproduced: r = 0 controls at A1 = 300, 1000, 3000 ADC, pooled over five modules -> fake rate 0.05 %, 0,
  0 as in `efficiency.tsv`.
- guard_reproduced: pairs r = 0.1, dt = 12 ns pooled -> ε 97.5 % ± 0.3 %.
- pool_fake_rate: pool singles of runs 24590 and 25320 -> F per amplitude band and for W493, ring 2, the rest,
  with errors (a number that does not exist yet).
- distorted_flagged: singles stretched by 1.6 -> w50 1.6 ± 0.1; fraction with k ≥ 3 and fraction flagged C2 at the
  X5 working point (true pairs flagged ≤ 2.5 %) reported; L_r written to `locus_limits.tsv`.
- Checkpoint: `exotic_bench.py pool` and `pseudo` write `bench_options.tsv` rows for the present constants and
  `locus_limits.tsv`.

**S5. Fit options and the jitter test.** Units: `XS/exotic_bench.py options`.
- o1_settings: the present constants, `MIN_SEP_NS` 12 and the six one-at-a-time settings of X10 -> F, ε, bias,
  spread per setting; the settings that pass named.
- o2_estimators: distorted singles and C2 samples -> bias and spread of the three estimators; W493 ev 2069998
  as a printed case (120 ADC main, one-pulse amplitude, sum 209 ADC, and the cluster energy with each).
- o3_first: stretch fit and principal components on pool pulses of W493, ring 2, and the rest -> s(A) and the
  first component's variance share.
- jitter: σ_j = 0, 0.5, 1, 2, 4, 8, 12 ns -> pool χ² median, amplitude bias, F on C2; X11 confirmed or reversed.
- Checkpoint: `exotic_bench.py options` completes `bench_options.tsv`; `results.md` names the candidate or states
  that none passes.

**S6. Module-complete class survey and what the classes are.** Units: `exotic_survey.py classes`, `XS/results.md`.
Completes the X7 table, maps, lists and examples with the phase B columns for the rescanned runs.
- classes_table: `scan_q/` -> fraction of C1-C6 per module and run; confirmed modules per class listed.
- validation: W121 ev 4084847 -> C3; W493 ev 2069998 -> C2; the rank of both modules in each list printed.
- episodes: C3 occupancy per module against its reference group, kept rows and rejected-row counters side by
  side -> list with z; L_4 set from the distribution (X8 A or D per module).
- board_coincidence: events with a C3 channel, all boards -> number with a second C3 channel on the same board
  against the number expected from single rates; B confirmed or refuted per board and L_5 set, or "too few"
  stated with the counts (fallback: X8).
- episode_history: time to the previous recorded pulse above 500 ADC in the same channel, C3 against unflagged ->
  X8 C (valid only if the timestamp unit is established, see Assumptions).
- width: w50 and tail against amplitude, seed and non-seed members apart, per module and ring, pool and k = 1
  pulses -> X9 hypotheses 2-4 and `template_candidates.tsv`.
- rate: per module, f_loc and C2 fraction against f_acc over runs of one type and between H2 and Ta -> slope
  with error; flat confirms a shape artefact, rising supports X9 hypothesis 1.
- cluster_agreement: C_t for C1 and C2 channels in Ta and 0.73 GeV -> separation power stated as a number.
- Checkpoint: after `exotic_survey.py classes` on `scan_q/`, every row of X8 and X9 has a verdict with its number,
  L_4 and L_5 are in `thresholds.tsv`, and the hardware question list names the confirmed C3/C4 modules.

**S7. Hand-off.** Units: `XS/results.md` sections "For P8" (candidate rule, limits files, final hot lists), "For
the P5 review" (class word per hit in BM, proposed per-run bad-channel table), "For the bank work"
(`template_candidates.tsv`), comparison with the P3 table when it exists.
- handoff_complete: section "For P8" -> names setting, F, ε and estimator; the other two sections name their table.
- Checkpoint: the supervisor is told the `Amend parent:` line of X12 is ready to apply.

## Verification

Each slice's checkpoint above. End to end: every PbWO4 module has a row for every run type in the final tables; the two events come out of
`exotic_survey.py classes` as C3 (W121) and C2 (W493) with limits that were set without them; the fit records of the rescan equal the old scan's for the same
events; `python3 waveform_analysis/pileup_analysis/pileup_rates.py report output/pileup_rates/scan_q <scratch dir>`
still runs (a scratch output directory: `pileup_rates_out/` is not overwritten).

## Risks

- Readout and selection: the effective DAQ threshold is not established (PP U6, P12) and PR adds its 8 and 10 ADC
  cuts; an effect below them is invisible to the coincidence test and belongs to the noise study, not here.
- The k = 1 reference contains mild members of the classes; the 0.1 % quantile is a convention (holdout: stability
  only). w50 is coarse at 4 ns sampling and may carry no information below 50 ADC (band reported apart).
- The generator's noise is white and its shape smear is shared by both pulses (`opens.md` item 13); F is therefore
  measured on real pool singles first and the pseudo-data only supply pairs and distortions.
- Rare classes: W121 has 5 `far` channels and 10 `bad_waveform` rows in run 24590; verdicts A-D may rest on tens
  of events per module over eight runs. If the error does not separate them, the five other runs are rescanned (Open Question 1).
- ED decodes one EVIO split per call; S2 times the first call and scales the number of example pages to an hour.

## Scope Notes

- Out: BM, `mf_pulse.py`, banks, database, C++; a pass over LMS or pulser events, except as the X8 fallback; noise
  below the readout and selection cuts (PP D9); cosmic or no-beam runs (P13 locates them); lead-glass templates; rebuilding any template.
- The `far` rule itself (PR:614) is not changed here although it turns a raised baseline into pulses; the class
  word marks those channels, and a change of the rule is part of the "For P8" recommendation.

## Assumptions

- `fe.info.timestamp` is a monotonic trigger time usable for intervals; its unit is checked in S3 against the event
  rate. If not usable, test episode_history is dropped and hypothesis C rests on slope and rate scaling.
- The DAQ reads out a channel from 10 ADC above pedestal (`../../database/daq_config.json:90`) until P12 reports.
- G's waveform builder and bank loader can be imported without running its `main` (checked in S4).
- PR reads `R_MODEL`, `BIC_PENALTY`, `KMAX` and `MIN_SEP_NS` at call time, so setting them on the imported module
  takes effect (checked by control_reproduced and one changed setting); where it does not, the bench caps the
  pulse count after the fit.
- Lead-glass channels are 0.13 % of the kept channels (about 30 per module in run 24590); S2 counts them, and
  rows below n = 100 are listed but not tested.
- Board membership from `$PRAD2_DATABASE_DIR/hycal_map.json` (crate, slot) is valid for all 13 runs.

## Log

- 2026-11-13 lead: S1 and S2 (phase A) done: `exotic_study/results.md`, `exotic_survey.py records` (3.4 min), 4 module tables, 11 anomalous lists, hot re-rank, 25-page maps, 62 example pages; checks table_complete, w121_counts, w493_multi, event_lookup, hot_rerank pass; W493 not on the f_loc list (z 0.8, typical inner-ring module) and the Ta/H2 f_acc median ratio is 3.4, not 20-115 (33 in rings 1-2), reported not tuned; deferred: S3-S7 (phase B decision by the user)

## Todo

- [x] S1 (phase A) record event-display results as the two validation cases — ED done — 0.1 d — `XS/results.md`
- [x] S2 (phase A) all-module table, maps, outlier lists, stability, hot lists, examples from the 13 scans — none —
      1 d — `XS/exotic_survey.py`, `XS/modules_*.tsv`, `XS/anomalous_*.tsv`, `XS/exotic_survey.pdf`
- [ ] S3 quality variables in PR, rescan 8 runs, limits — after the S2 numbers are reported; free CPUs (PP D11) —
      1 d + 1.4 h scan — `output/pileup_rates/scan_q/`, `XS/thresholds.tsv`
- [ ] S4 truth bench: pool singles, pairs, distorted singles — after S3's functions; parallel to the rescan — 1 d —
      `XS/exotic_bench.py`, `XS/bench_options.tsv`, `XS/locus_limits.tsv`
- [ ] S5 options O1-O3 and jitter test — after S4 — 1 d — `XS/bench_options.tsv`, `XS/exotic_bench.pdf`
- [ ] S6 full class table, lists, class tests, template candidates, hardware questions — after S1, S2, S3, S4 — 1.5 d —
      `XS/results.md`, `XS/template_candidates.tsv`
- [ ] S7 hand-off to P3/P5/P8 and the bank work, `Amend parent:` applied by the supervisor — after S5, S6; P3 table if present — 0.5 d
