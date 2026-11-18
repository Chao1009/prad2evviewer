# crosstalk_study

Independent study of the crosstalk issue in the HyCal non-seed cluster members: dt = T_module - T_seed, amplitude and
matched-filter (MF) chi2/dof of every non-seed member of single-cluster e-p events. Run types: `0p73GeV_H2` (runs 24590, 24603; 10 + 6 splits),
`2p24GeV_H2` (run 25228, 11 splits) and `3p49GeV_H2` (runs 24386, 24503, 24512; 10 + 10 + 7 splits = the 27 splits of the reference chain).
There is no `2p24GeV_X17` run type: the chain's X17 type is `2p24GeV_Ta_X17` (Ta target), analysed per run (never pooled) as pseudo run types
`2p24GeV_Ta_X17_<run>` (runs 25800, 26138, 26230; 8, 4 and 13 splits), figures in `figures/2p24GeV_Ta_X17_<run>/`. All
commands are run from `analysis/pyscripts` after `source setup_env.sh`. Everything lives in this directory; the frozen benchmark,
`database/` and `runinfo/` are not modified. Nothing here is committed (cache and figures are regenerable).

## Dependencies
- Frozen `benchmark_hycal_timing_test_20261021.py` (sha1 3293ab96870cef4580e886eb2dda5ec37c8821d1), chain job lines
  `output/r5ri_20261104/jobs/final.txt` + `--cluster-times --pulse-dump`.
- prad2py with `dec.WaveAnalyzer.analyze_fit` / `Peak.time_algo` / `LogNormalFitResult`: this comes from an uncommitted change of
  `python/bind_dec.cpp` (see `README_prad2py_fit.md`; build `cd build && make prad2py`). Without it `--cfd-algo fit` stops with a message;
  `--cfd-algo all` and `--time-source mf` work without it.
- `database/hycal_time_offsets/time_calib.json` (sha1 447d610db1..., 1152 entries, 1150 resolved, reference W495 = 0).

## Conventions
- **Sign (changed 2026-11-06):** dt = T_module - T_seed, positive = module later than the seed (the convention of
  `output/tseed_dt_20261030` and the 10-30 note). Earlier versions of this study used T_seed - T_module (negative of dt);
  all figures were re-plotted with the new sign. The crosstalk bump of the 10-30 note sits at dt = +7...+15 ns (shaded on fig1a, fig1c as
  "crosstalk region (10-30 note)"), inside the -5...25 ns windows of fig1b, fig1c, fig2, fig2b.
- Cross-check of the raw shape against the 10-30 sample (raw legacy, no selection, MF time, sign-matched): old run 025320
  (N = 5.60 M): |dt| > 10 ns 0.788 %, dt in [+7,+15] 3.58 %, dt in [-15,-7] 0.87 %, median +0.501 ns; this study, run 25228
  (N = 15.12 M): 0.798 %, 3.60 %, 0.86 %, +0.508 ns. Same raw shape; calibration and gating change it.
- Hits: every non-seed member of the cluster (`members` column of `result_clusters.tsv`, first member = seed), legacy path (ungated
  cluster, the path of the calibration pool) unless `--path gated` (seed time window 16 ns, the chain's `--window 16`).
- Amplitude (ADC) and chi2/dof come from the MF row of the `--pulse-dump` matched to the member time (|dt| < 0.01 ns); they are the
  same for every time source.

## Selection (default `--selection ep`; `single` = previous behaviour)
Cut-flow, legacy path (counts of the cache; Moller-tag counts overlap):

| step | 2.24 GeV H2 (run 25228) | 3.49 GeV H2 (24386, 24503, 24512) |
|---|---|---|
| single-cluster events (exactly one cluster > 500 MeV; applied by the benchmark, its cluster TSV keeps only these) | 1,424,607 | 2,360,537 |
| tagged by `moller_pair` / `moller_kin` / `moller_rej` | 278,315 / 392,130 / 392,130 | 318,421 / 329,775 / 329,775 |
| not Moller-tagged (pair = kin = rej = 0) | 1,014,738 | 1,970,251 |
| in the elastic window, |r - 1| < 3 sigma_CB | 523,352 | 1,191,910 |
| sigma_CB / E (chain `result.tsv`, `sigma_over_E_cb` of the legacy row, weighted by `n_inner_events`) | 0.01805 | 0.01594 |
| r window | [0.9459, 1.0541] | [0.9522, 1.0478] |

- Moller tag (benchmark TSV columns, written because the chain runs `--moller-reject kin`): `moller_pair` = another cluster with
  |E1 + E2 - Ebeam| <= 0.10 Ebeam and |dphi - 180| <= 10 deg; `moller_kin` = |E - E_ee(theta)| < 3 sigma and |E - E'(theta)| > 2 sigma;
  `moller_rej` = the rejection actually used (= kin in this chain). An event is kept only if all three are 0 (any tag removes it).
  The cache stores the bits (1 pair, 2 kin, 4 rej) in `moller`.
- r = E_cluster / E'_ep(seed): E_cluster = `cluster_E_raw_MeV` of the TSV (already contains the chain gain file and pos-corr: with
  `--gain-refine off` `cluster_E_cal_MeV` is identical where it exists, and it exists only for the inner-ring seeds with a gain);
  E' = `expected_energy_elastic(theta_seed, Ebeam, m_p)` of the benchmark (imported, not re-implemented), theta from the seed-module
  centre (`seed_x_mm`, `seed_y_mm`) with the HyCal z of the benchmark's own `result_window.tsv`; Ebeam and target mass from the rerun's
  `result.tsv` header (2239.51 MeV, 3488.43 MeV; 938.272 MeV). The peak sits at r = 1.0007 (2.24 GeV), the window is centred on 1.
- Hits of kept events only; N_ev counts kept events (events whose members all lack an MF or WaveAnalyzer match still count).

## Time source and algorithm (`--time-source cfd|mf`, `--cfd-algo fit|all`)
- `mf` (the previous basis): MF pulse time t = 50 % leading edge of the template fit.
- `cfd` (default): WaveAnalyzer time `pk.time` (ns) = method-0 row of the pulse dump (the CFD pick-off), the peak nearest in time to
  the MF pulse of the same (event, module) within 10 ns; hits without such a peak (12.7 % at 2.24 GeV, 11.1 % at 3.49 GeV, small pulses
  below the WaveAnalyzer threshold; median MF amplitude 9.9 ADC) are dropped from CFD figures. Median (t_WA - t_MF) = -0.60 ns (2.24),
  -0.54 ns (3.49); |t_WA - t_MF| > 1 ns for about 29 % of matched hits.
- `time_algo` and the log-normal fit chi2/dof come from a per-channel `analyze_fit` pass over the EVIO (about 75 s per split, run
  automatically when the cache lacks them, `--jobs` parallel) matched to the CFD time (|dt| < 0.01 ns). Algorithm breakdown of the
  matched hits (e-p sample): seeds 100.00 % log-normal fit (2); non-seed hits 87.3 % fit (2), 0.05 % linear CFD (1),
  0 % peak (0), 12.7 % unmatched at 2.24 GeV; 88.8 / 0.05 / 0 / 11.1 % at 3.49 GeV. So the fit-only restriction costs nothing beyond
  the 0.05 % linear-CFD hits; the any-algorithm set (`--cfd-algo all`, directory `<slug>_cfd_all`) differs by that amount.
- fig2b (dt vs log-normal fit chi2/dof of the module) exists for the fit-only CFD sets.

## Time calibration (`--tcal off|auto|<json>`, default auto)
- `auto` = `database/hycal_time_offsets/time_calib.json`, applied as t = t - offset to the seed and the module (benchmark
  :1616-1619), for MF and CFD alike. Basis: the table is derived by `analysis/tools/hycal_module_time_calib.cpp` from the stored
  `ev->peak_time` (`hycal_module_time_calib.cpp:213-250`), which Replay fills from the WaveAnalyzer peak time
  (`Replay.cpp:523-532`, `WaveAnalyzer.cpp:572`), raw and without offset (`plan_20261022.md:257`): it is native to the CFD time. For the MF
  time (t0 + template t50) it is what the benchmark has always applied but it is cross-basis.
- Every figure title carries "tcal: ..."; with offsets the raw figures are written to the sibling `<out>_raw`. fig1a / fig1c print the core
  width (robust sigma = half of the 16-84 % range; sigma of a Gaussian fit to median +- 3 ns, 0.1 ns bins), raw and calibrated.
- Self-referenced offsets (`--selfref`, cross-check only, `tcal/<slug>_{mf,cfd}_selfref_offsets.tsv`: module, n, offset_ns, se): per
  (module, seed) pair the median dt over the e-p sample, solved as a weighted least-squares pair network (gauge: median 0), windows
  re-centred 20/5/3/3 ns. They reproduce the core width of the official table (Gauss 0.676 vs 0.679 ns at 2.24 GeV CFD on four
  splits) but only the local differences agree with it (RMS over seed-neighbour pairs 0.27 ns CFD, 0.41 ns MF at 2.24 GeV;
  0.50 / 0.46 ns at 3.49 GeV); the module-by-module RMS is 1.7 (MF) / 1.8 (CFD) ns at 2.24 GeV and 1.9 / 2.4 ns at 3.49 GeV
  (smooth large-scale modes the local pair network does not fix, median removed).

## Cache (`cache/<run-type>/<run>_<split>_<path>.npz`, per hit unless noted)
`ev, seed_t, seed_amp, seed_E, t, amp, chi2, E` (legacy columns, kept), `mod_i, seed_mod_i` (index into `module_names`),
`cl_E, eprime, r, moller` (cluster energy, E', r, Moller bits), `seed_t_wa, t_wa` (WaveAnalyzer times, nan = no match),
`algo, lnchi2, seed_algo, seed_lnchi2` (added by the WaveAnalyzer pass; algo -1 = unmatched); per cluster `cl_ev, cl_r, cl_moller`;
scalars `sigma_cb, n_inner, ebeam, eprime_hz, n_clusters, n_nomatch, n_events_dump`. Inventory: `cache/0p73GeV_H2` (16 legacy splits), `cache/2p24GeV_H2` (11 legacy + 11 gated), `cache/3p49GeV_H2` (27 legacy); about 0.3 to 0.8 GB per run type.

## Commands
Benchmark outputs are kept with `--work` so that later passes do not rerun the benchmark. For 3.49 GeV the 27 benchmark jobs were
run beforehand in `<scratch>/ct_3p49_bench/<run>_<split>/` with exactly the job lines this script builds (checked line by line:
27 of 27 identical); the script ingests a directory when its `rc` is 0 and `pulse.npz` exists.

    source setup_env.sh
    W=<scratch>/ct_2p24_bench      # benchmark outputs (about 2 CPU-min per split, sequential about 22 min for 11 splits; --jobs 4 shortens it)
    # first run / regenerate: benchmark + cache + WaveAnalyzer pass + figures (main set: e-p, CFD fit-only, calibrated; raw sibling _raw)
    python3 crosstalk_study/crosstalk_study.py --run-type 2p24GeV_H2 --out crosstalk_study/figures/2p24GeV_H2 --jobs 4 --work $W --selfref
    # variants (cache only, minutes); same for --run-type 3p49GeV_H2 (default splits = all 27)
    python3 crosstalk_study/crosstalk_study.py --run-type 2p24GeV_H2 --reuse --out crosstalk_study/figures/2p24GeV_H2_cfd_all --cfd-algo all
    python3 crosstalk_study/crosstalk_study.py --run-type 2p24GeV_H2 --reuse --out crosstalk_study/figures/2p24GeV_H2_mf --time-source mf
    python3 crosstalk_study/crosstalk_study.py --run-type 2p24GeV_H2 --reuse --out crosstalk_study/figures/2p24GeV_H2_single --selection single
    python3 crosstalk_study/crosstalk_study.py --run-type 2p24GeV_H2 --reuse --out crosstalk_study/figures/2p24GeV_H2_gated --path gated
    # 3.49 GeV, benchmark outputs from the pre-run:
    python3 crosstalk_study/crosstalk_study.py --run-type 3p49GeV_H2 --out crosstalk_study/figures/3p49GeV_H2 --jobs 6 --work <scratch>/ct_3p49_bench --selfref
    # single-cluster selection (no e-p cuts) for any run type, e.g. Ta per run:
    python3 crosstalk_study/crosstalk_study.py --run-type 2p24GeV_Ta_X17_25800 --reuse --no-raw --selection single --out crosstalk_study/figures/2p24GeV_Ta_X17_25800_single
    # regenerate only the main sets from the cache (no raw sibling), optional crosstalk window for fig4:
    python3 crosstalk_study/crosstalk_study.py --run-type 2p24GeV_H2 --reuse --no-raw --out crosstalk_study/figures/2p24GeV_H2 --ct-window 7 15 --core-window 2
    # 0.73 GeV H2 (16 splits; benchmark pre-run into <scratch>/ct_0p73_bench, then the WaveAnalyzer pass and figures):
    python3 crosstalk_study/crosstalk_study.py --run-type 0p73GeV_H2 --out crosstalk_study/figures/0p73GeV_H2 --no-raw --jobs 8 --work <scratch>/ct_0p73_bench --selfref

Default `--splits` = all splits of the run type (`TYPE_SPLITS`); without a known run type, `25228/50`.

`--reuse` and the cache: the cache is used when it has all of `r, moller, t_wa, seed_t_wa, mod_i, seed_mod_i, module_names`
(`CACHE_KEYS`); otherwise (e.g. a cache of an earlier version) the split is re-reduced, which reuses the benchmark output in `--work`
if present and reruns the benchmark otherwise (a message is printed). The WaveAnalyzer pass (`algo` columns) runs only for
`--time-source cfd --cfd-algo fit` and is skipped when the columns exist. All flag combinations of the commands above work with `--reuse`
once the cache is complete; `--selection single` needs no extra data, `--tcal` needs only `mod_i`.

## Figure sets (`figures/`; each has a `_raw` sibling without time calibration; titles state selection, N_ev, N_hits, time source, tcal)
| directory | selection | time | path |
|---|---|---|---|
| `<slug>` | e-p | WaveAnalyzer CFD, log-normal fit only | legacy |
| `<slug>_cfd_all` | e-p | WaveAnalyzer CFD, any time_algo | legacy |
| `<slug>_mf` | e-p | matched filter | legacy |
| `<slug>_single` | single-cluster only | WaveAnalyzer CFD, fit only | legacy |
| `2p24GeV_H2_gated` | e-p | WaveAnalyzer CFD, fit only | gated (W = 16 ns) |

Figures (pure 2D histograms; the median-per-bin profile overlays were removed 2026-11-06): fig1a_dt_hist (+-60 ns, log y), fig1b_dt_vs_amp
(-5...25 ns, MF amplitude), fig1c_dt_hist_zoom (linear and log y panels, one short title line with run type and N_ev, axis labels only, the
`--ct-window` band shaded, no legend or text boxes), fig2_dt_vs_chi2 (MF chi2/dof, linear y clipped at the 99 % quantile of the window),
fig2b_dt_vs_lnchi2 (log-normal fit chi2/dof), fig3_amp_vs_amp (seed vs module MF amplitude), fig4_amp2d_crosstalk (same axes and binning as
fig3, only hits with `--ct-window LO HI` <= dt <= HI, default +7...+15 ns; text box "N = ..."), fig4b_amp2d_core (|dt| < 5 ns, identical
axes). fig3, fig4, fig4b have no ratio lines and no identity line (dropped: the log axes are identical, the identity is not needed to read them)
and fixed axes 5...5000 ADC (hits outside are clipped into the edge bins), identical for all run types. The shading on fig1a also stays.

The main sets `figures/0p73GeV_H2`, `figures/2p24GeV_H2` and `figures/3p49GeV_H2` (e-p, calibrated, CFD fit-only) are regenerated in this version
(`--reuse --no-raw`). Fig4 and fig4b are present in these, in the per-run Ta sets `figures/2p24GeV_Ta_X17_<run>`, and in the `_single` (100 MeV),
`_single_e500`, `_nmod20` and `_single_nmod20` sets (checked 2026-11-07 for every run type where the set exists). The sibling sets `_raw`, `_mf`,
`_mf_raw`, `_cfd_all`, `_cfd_all_raw`, `_single_raw`, `_gated` and `_gated_raw` have no fig4 or fig4b: they were left as they were and predate the
median removal (they still show the median overlay on fig1b / fig2 / fig2b). Regenerating one needs only its command from "Commands".

### fig4 / fig4b amplitude statistics (e-p, calibrated CFD fit-only)
Window = +7...+15 ns (fig4); core = |dt| < 5 ns (fig4b), not the |dt| < 2 ns core of the ratio tables below; MF amplitude; ratio = median non-seed / seed amplitude.

| sample | selection | window N | window median amp (ADC) | window median amp ratio | window frac < 10 ADC | core median amp (ADC) | core median ratio |
|---|---|---|---|---|---|---|---|
| 0.73 GeV H2 | e-p | 6,690 | 9.53 | 1.64 % | 61.8 % | n/a | n/a |
| 2.24 GeV H2 | e-p | 446,909 | 9.46 | 0.60 % | 56.7 % | 32.5 | 2.14 % |
| 3.49 GeV H2 | e-p | 1,803,823 | 11.49 | 0.50 % | 35.4 % | 35.0 | 1.58 % |
| Ta X17 (2p24GeV_Ta_X17_25800, _26138, _26230; range over the runs) | e-p | n/a | 9.2...10.3 | 0.57...0.71 % | 47...60 % | n/a | n/a |

`n/a` = not reported (for Ta the note gives only the window ranges). Core N (no column): 2.24 GeV 4,891,526; 3.49 GeV 14,208,472;
core frac < 10 ADC 0.20 % for both H2 rows.

Source: script `[fig4]`/`[fig4b]` stdout as reported in the 2026-11-06/07 run logs.

### Crosstalk and core ratios (`tables/crosstalk_ratio.md`, written by the script; main sample, legacy path)
(a) events with >= 1 non-seed hit in the region / total e-p events (N_ev, events without fit-only hits included); (b) hits in the region /
total non-seed hits; h/ev = region hits / N_ev. Crosstalk window `--ct-window` = +7...+15 ns; core `--core-window` = |dt| < 2 ns.
| run type | N_ev | non-seed hits | (a) window | (b) window | h/ev window | (a') core | (b') core | h/ev core |
|---|---|---|---|---|---|---|---|---|
| 0.73 GeV H2 | 1,080,125 | 4,720,970 | 0.61 % | 0.14 % | 0.006 | 99.98 % | 94.38 % | 4.125 |
| 2.24 GeV H2 | 523,352 | 5,512,007 | 66.58 % | 8.11 % | 0.854 | 99.95 % | 81.23 % | 8.555 |
| 3.49 GeV H2 | 1,191,910 | 16,571,558 | 87.07 % | 10.89 % | 1.513 | 100.00 % | 78.54 % | 10.920 |
Ta X17 (per run, same recipe; sigma_CB/E from each run's chain `result.tsv`):
| run | single-cluster events | not tagged | sigma_CB/E, r window | N_ev | non-seed hits | (a) window | (b) window | h/ev | (a') core | (b') core | h/ev core |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 25800 | 88,098 | 81,684 | 0.02049, [0.9385, 1.0615] | 61,332 | 601,787 | 51.79 % | 6.56 % | 0.643 | 100.00 % | 79.66 % | 7.816 |
| 26138 | 19,304 | 17,669 | 0.02121, [0.9364, 1.0636] | 14,222 | 130,083 | 42.10 % | 5.45 % | 0.498 | 100.00 % | 80.92 % | 7.401 |
| 26230 | 50,721 | 42,110 | 0.02125, [0.9362, 1.0638] | 32,150 | 367,824 | 68.19 % | 9.25 % | 1.058 | 99.99 % | 69.11 % | 7.907 |
Job lines: 25 of 25 identical to `final.txt` + the two flags, rc 0.

0.73 GeV H2 check of the threshold hypothesis: 14.67 % of the non-seed hits of the e-p sample have no WaveAnalyzer peak (median MF amplitude 10.1 ADC,
matched hits 27.1 ADC; median seed amplitude 496 ADC against about 1700 ADC at 2.24 GeV). The MF-time sample (`--time-source mf`, row "MF" in
`tables/crosstalk_ratio.md`) gives the same picture (window ratios 0.67 % / 0.13 %, 7,410 hits): its window hits have a median amplitude of 6.7 ADC and
90 % are below 10 ADC, i.e. they sit at the MF floor (min amplitude 5 ADC); the expected crosstalk of about 0.6 % of the seed (3 ADC) is below it. So
the data are consistent with a threshold effect but this sample cannot separate it from an absence.

### Single-cluster selection (`--selection single`), 500 x kE MeV definition (`figures/<slug>_single_e500/`), ratios beside the e-p rows
Selection side by side: **e-p** = single-cluster, not Moller-tagged (pair|kin|rej = 0), |r - 1| < 3 sigma_CB; **single** = every event with exactly
one cluster above 500 MeV (the benchmark's own sample), no tag veto, no r cut. Both: legacy path, calibrated WaveAnalyzer CFD (time_algo 2 at seed
and module), window +7...+15 ns, core |dt| < 2 ns. The full table with a `selection` column is `tables/crosstalk_ratio.md` (one json per sample in `tables/`).
| sample | N_ev single | non-seed hits | (a) window | (b) window | h/ev | (a') core | (b') core | h/ev core |
|---|---|---|---|---|---|---|---|---|
| 0.73 GeV H2 | 2,893,353 | 11,273,453 | 0.31 % | 0.09 % | 0.003 | 99.10 % | 95.65 % | 3.727 |
| 0.73 GeV H2, MF time | 2,893,353 | 13,100,155 | 0.37 % | 0.09 % | 0.004 | 99.37 % | 95.28 % | 4.314 |
| 2.24 GeV H2 | 1,424,607 | 13,204,084 | 44.88 % | 5.89 % | 0.546 | 99.64 % | 84.61 % | 7.842 |
| 3.49 GeV H2 | 2,360,537 | 30,273,407 | 75.50 % | 9.54 % | 1.223 | 99.99 % | 80.33 % | 10.302 |
| Ta X17 25800 | 88,098 | 822,643 | 42.37 % | 5.72 % | 0.534 | 99.84 % | 79.26 % | 7.401 |
| Ta X17 26138 | 19,304 | 170,366 | 35.83 % | 4.93 % | 0.435 | 99.89 % | 80.03 % | 7.063 |
| Ta X17 26230 | 50,721 | 565,621 | 66.42 % | 9.36 % | 1.044 | 99.94 % | 67.70 % | 7.550 |
The single-cluster window ratios (a) are lower than the e-p ones (2.24 GeV: 44.9 % vs 66.6 % of events), the hit fractions (b) differ less (5.9 % vs 8.1 %); the cause was not studied. The wedge in fig3 (seed amplitude below about
500 ADC with module amplitude comparable to the seed) is present in all `_single` sets.

### Summary figure: dt (CFD fit) vs MF chi2/dof, three energies (`summary_dtcfd_vs_mfchi2.py`, from the caches only)
`python3 crosstalk_study/summary_dtcfd_vs_mfchi2.py` writes `figures/summary/dtcfd_vs_mfchi2.{png,pdf}` and `tables/dtcfd_vs_mfchi2.md`: 2 x 3 panels, rows e-p
(`cache/<slug>`, 500 x kE threshold) and single-cluster > 100 MeV (`cache/<slug>_e100`), columns 0.73, 2.24, 3.49 GeV H2; x = calibrated CFD log-normal-fit
dt = T_module - T_seed (-5...25 ns, +7...+15 ns band), y = MF chi2/dof (linear, common range 0..81 = 99 % quantile over all six panels in -5...25 ns, clipped
into the top row), log colour scale. It is the same quantity as the per-sample `fig2_dt_vs_chi2` of `figures/<slug>/` (e-p) and `figures/<slug>_single/` (100 MeV)
(there the y range is per sample). Thresholds 2 and 3 are the good/bad cuts of waveform_good_bad.py.
Window +7..+15 ns, core |dt| < 2 ns; hits = non-seed hits with time_algo 2 at seed and module; N_ev = events of the selection.

Window +7..+15 ns, core |dt| < 2 ns; hits = non-seed hits with time_algo 2 at seed and module; N_ev = events of the selection.

Core fractions = false-rejection rate if chi2/dof were used as a crosstalk cut.

| sample | selection | N_ev | N hits | window hits | core hits | median window | median core | window > 2 | > 3 | > 5 | > 10 | core > 2 | > 3 | > 5 | > 10 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.73 GeV H2 | e-p | 1,080,125 | 4,720,970 | 6,690 | 4,455,810 | 1.948 | 1.097 | 47.4 % | 29.8 % | 25.8 % | 14.0 % | 30.8 % | 22.0 % | 15.7 % | 8.3 % |
| 2.24 GeV H2 | e-p | 523,352 | 5,512,007 | 446,909 | 4,477,305 | 8.114 | 1.250 | 96.4 % | 88.3 % | 75.2 % | 36.2 % | 35.7 % | 25.7 % | 16.8 % | 8.4 % |
| 3.49 GeV H2 | e-p | 1,191,910 | 16,571,558 | 1,803,823 | 13,015,430 | 13.381 | 1.315 | 98.7 % | 94.0 % | 79.7 % | 59.6 % | 36.6 % | 26.5 % | 17.8 % | 9.1 % |
| 0.73 GeV H2 | single 100 MeV | 2,713,040 | 10,638,306 | 9,207 | 10,166,096 | 1.946 | 1.066 | 47.6 % | 31.7 % | 27.6 % | 16.2 % | 30.2 % | 21.6 % | 15.6 % | 8.2 % |
| 2.24 GeV H2 | single 100 MeV | 968,883 | 9,155,444 | 574,705 | 7,705,556 | 7.697 | 1.235 | 96.1 % | 87.5 % | 72.1 % | 34.9 % | 35.4 % | 25.4 % | 16.9 % | 8.5 % |
| 3.49 GeV H2 | single 100 MeV | 1,723,192 | 22,950,608 | 2,349,590 | 18,220,828 | 12.667 | 1.316 | 98.5 % | 93.3 % | 78.7 % | 57.4 % | 36.7 % | 26.6 % | 17.8 % | 9.1 % |

Note on the earlier headline table: its "chi2/dof in / out" column was the median for |dt| < 5 ns / |dt| > 20 ns (calibrated; N = 4,891,526 / 5,029 at 2.24 GeV e-p), not |dt| < 10 / > 10 ns. The +7...+15 ns window (median 8.1) lies in neither bin, so the two statistics do not conflict.

### Profile: MF chi2/dof median and 16-84 % band vs CFD dt (`profile_mfchi2_vs_dtcfd.py`)
`figures/summary/mfchi2_profile_vs_dtcfd.{png,pdf}`: two panels (left linear y 0-20, right log y), three curves (0.73, 2.24, 3.49 GeV H2; single-cluster > 100 MeV,
`cache/<slug>_e100`, calibrated CFD fit-only dt), 0.5 ns bins from -5 to 25 ns, bins with < 100 hits skipped, +7...+15 ns band shaded. The curve is the per-bin median
MF chi2/dof; the filled band is the 16th-84th percentile interval (asymmetric, not a Gaussian sigma). Binned values: `tables/mfchi2_profile_vs_dtcfd.tsv`
(energy, dt centre, N, median, p16, p84). Pooled values (median [p16, p84]):
| sample | core |dt| < 2 ns | window +7...+15 ns |
|---|---|---|
| 0.73 GeV H2 | 1.066 [0.569, 4.739] | 1.946 [1.322, 10.702] (N = 9,207) |
| 2.24 GeV H2 | 1.235 [0.620, 5.410] | 7.697 [3.427, 33.047] (N = 574,705) |
| 3.49 GeV H2 | 1.316 [0.651, 5.860] | 12.667 [4.263, 52.129] (N = 2,349,590) |

### Examples: MF chi2/dof > 5 with |dt(CFD)| < 2 ns (`examples_mfchi2_gt5.py`)
`figures/examples/mfchi2_gt5_core_2p24GeV_H2.pdf` (caption page + 9 pages, 2 x 3 panels): 50 random non-seed hits of 2.24 GeV H2 (single-cluster > 100 MeV, fit-only calibrated
|dt| < 2 ns, MF chi2/dof > 5; population 1,302,480 hits; numpy default_rng seed 20261107; 4-5 per split, at most one per event and module). Each panel: raw samples, MF fit
(benchmark bank, nsigma 5, min_amp 5, kmax 4, r_model 0.3, mad noise; fit window 32 samples = 128 ns from 3 samples before ceil(t0/4); A and t0 free, dof = valid samples - 2), residual, and the markers
t_MF, t_CFD (module) and t_CFD (seed). The recomputed chi2/dof equals the cached one exactly (max |diff| = 0 over the 50; all have K = 1). By-eye tally of the 50 (rough): 23 low-amplitude pulses
(A below about 30 ADC; residuals are sample-to-sample scatter and ADC-step noise against a small sigma_ped), 24 clean single pulses of larger amplitude with structured residuals at the rise/peak/tail (template-shape mismatch, residuals
+-5 to +-20 ADC; chi2/dof 5-27), 3 marginal ones (A 30-60 ADC, chi2/dof 5-9, noise-like). None showed a second pulse, clipping or a baseline step.

### Fired-module requirement n_mod < 20 (cache-only; `--nmod-max 20`, figures `<slug>_nmod20/` (e-p) and `<slug>_single_nmod20/`)
Source: the benchmark has no per-event fired-module count, but its `--pulse-dump` contains every fired PbWO4 channel of every kept event (not only cluster members), so the count is taken from the
dump in the benchmark output dirs (`--work`), no rerun. `n_mod` (cache columns `nmod_wa20`, `nmod_mf5`, cluster-level `cl_nmod_*`): **`wa20` (default, the N_FIRED = 20 criterion of waveform_good_bad.py / waveform_examples.py with H_MIN = 20 ADC)** = number of distinct PbWO4
modules with a WaveAnalyzer peak (method-0 row) of height >= 20 ADC; the quality/overflow flags are not in the dump (they are rare); `mf5` = distinct modules with an MF pulse (method-1 row, amplitude >= min_amp 5 ADC, `--nmod-def mf5`).
Selection: single-cluster > 100 MeV AND n_mod < 20 (then the e-p cut). Ratio rows carry an `n_mod` column (`all`, `<20 (wa20)`, `<20 (mf5)`); the figures use `wa20`. 
Events before -> after (single-cluster > 100 MeV; n_mod median / p90 in brackets, wa20 | mf5): 0.73 GeV 2,713,040 -> 2,712,927 (4 / 6 | 6 / 9) | 2,710,615; 2.24 GeV 968,883 -> 965,913 (8 / 10 | 14 / 17) | 927,201;
3.49 GeV 1,723,192 -> 1,715,862 (11 / 13 | 17 / 21) | 1,407,100; Ta 25800 69,837 -> 65,234 (11 / 18 | 19 / 30) | 36,319; Ta 26138 15,217 -> 14,351 (10 / 17 | 19 / 30) | 8,142; Ta 26230 39,218 -> 30,419 (15 / 23 | 24 / 36) | 10,496.
e-p (100 MeV, n_mod all -> wa20): 0.73 GeV 1,078,406 -> 1,078,387; 2.24 GeV 519,622 -> 517,957; 3.49 GeV 1,184,917 -> 1,179,923; Ta 51,721 -> 48,234, 11,828 -> 11,154, 25,921 -> 19,876.
With the wa20 definition the cut removes under 1 % of events at H2 (it matters for Ta: 7-22 %), so the ratios barely move; with mf5 it removes 0.1 % (0.73 GeV), 4 % (2.24 GeV) and 18 % (3.49 GeV) of the H2 events and 46-73 % of the Ta events; the ratios change by up to about 4 points (Ta).
Full rows (all samples, both definitions, with and without the cut, e-p and single, 100 MeV): `tables/crosstalk_ratio.md`. Summary figures with the cut:
`figures/summary/dtcfd_vs_mfchi2_nmod20.{png,pdf}` (`tables/dtcfd_vs_mfchi2_nmod20.md`) and `figures/summary/mfchi2_profile_vs_dtcfd_nmod20.{png,pdf}` (`tables/mfchi2_profile_vs_dtcfd_nmod20.tsv`), made with `--nmod 20`; both rows of the grid come from `cache/<slug>_e100`.
Commands: `python3 crosstalk_study/crosstalk_study.py --run-type <slug> --cluster-thr 100 --reuse --no-raw --work <scratch>/ct_<slug>_e100 --selection ep --nmod-max 20 --out crosstalk_study/figures/<slug>_nmod20`
(`--selection single` for `<slug>_single_nmod20`; `--nmod-def mf5` for the alternative definition).

### Seed-module MF chi2/dof (`seed_chi2.py`; e-p, single-cluster > 100 MeV, calibrated CFD fit-only; `tables/seed_chi2.md`)
Seed chi2/dof = chi2/dof of the seed's own MF pulse (method-1 row of the pulse dump that gives seed_t / seed_amp), cache column `seed_chi2` (hit-level) in `cache/<slug>_e100`;
missing for 1 of 12.4 M hits (0.73 GeV), 316 of 10.5 M (2.24 GeV), 13 of 25.9 M (3.49 GeV), none for Ta (no MF row within 0.01 ns of seed_t). Figures in `figures/<slug>/`:
`fig5_seedchi2_vs_amp` (x = non-seed MF amplitude, log 5-5000 ADC, y = seed chi2/dof linear clipped at the 99 % quantile), `fig5b_seedchi2_vs_amp_logy`, `fig5c_seedchi2_vs_dt` (x = CFD dt, band +7..+15 ns).
One entry per (event, non-seed member): an event with 10 members enters 10 times at the same y.
# Seed-module MF chi2/dof (e-p, single-cluster > 100 MeV, calibrated CFD fit-only hits)

Seed chi2/dof = chi2/dof of the seed's own MF pulse (benchmark pulse dump). Event-level values: one per event with at least one selected non-seed hit; hit-level: one per (event, non-seed member). 'seeds without MF chi2' counts hits of the whole cache (before the e-p cut).

| sample | events | hits | seeds without MF chi2 | seed chi2 events: median [p16, p84] | window hits +7..+15 ns | core hits abs(dt)<2 ns | events > 2 | > 3 | > 5 | Spearman(seed chi2, amp) | Spearman(seed chi2, dt) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.73 GeV H2 | 1,078,325 | 4,712,544 | 1 of 12,370,747 | 1.173 [0.637, 2.829] | 0.996 [0.580, 2.200] (N=6,575) | 1.180 [0.639, 2.845] (N=4,448,837) | 26.2 % | 14.6 % | 6.0 % | +0.075 | -0.011 |
| 2.24 GeV H2 | 519,539 | 5,467,302 | 316 of 10,480,934 | 1.006 [0.458, 2.803] | 0.999 [0.458, 2.718] (N=443,366) | 1.001 [0.458, 2.774] (N=4,445,876) | 24.5 % | 14.6 % | 7.0 % | +0.003 | +0.026 |
| 3.49 GeV H2 | 1,184,915 | 16,463,532 | 13 of 25,873,795 | 2.077 [0.559, 10.373] | 2.066 [0.557, 9.720] (N=1,793,451) | 2.003 [0.552, 10.116] (N=12,936,494) | 50.9 % | 41.4 % | 30.6 % | +0.005 | +0.019 |
| Ta X17 25800 | 51,721 | 497,265 | 0 of 754,584 | 1.081 [0.476, 3.112] | 1.066 [0.469, 2.937] (N=32,339) | 1.057 [0.469, 3.035] (N=404,279) | 27.7 % | 16.8 % | 8.3 % | +0.008 | +0.016 |
| Ta X17 26138 | 11,828 | 105,770 | 0 of 154,767 | 0.952 [0.461, 2.704] | 1.010 [0.472, 2.754] (N=5,625) | 0.935 [0.457, 2.638] (N=87,457) | 23.5 % | 13.9 % | 6.3 % | -0.000 | +0.042 |
| Ta X17 26230 | 25,919 | 286,110 | 0 of 491,724 | 0.903 [0.451, 2.461] | 0.925 [0.453, 2.594] (N=26,594) | 0.892 [0.447, 2.398] (N=205,682) | 20.4 % | 12.6 % | 7.0 % | +0.003 | +0.027 |
Reading: in every sample the seed chi2/dof is uncorrelated with the non-seed amplitude and with dt (|Spearman| <= 0.075; window hits have the same seed chi2 as core hits), so the seeds of crosstalk-window hits are not fitted worse.
Seed chi2 is about 1.0 (24-28 % of events above 2) at 2.24 GeV and Ta, 1.17 at 0.73 GeV, and larger at 3.49 GeV (median 2.1, 31 % of events above 5).

### Single-cluster threshold: 100 MeV (current `figures/<slug>_single/`) and 500 MeV (`figures/<slug>_single_e500/`, kept for reference)
The single-cluster cut is the benchmark option `--signal-min` (default 500 x kE MeV: 162.7 MeV at 0.73 GeV, 500 MeV at 2.24 GeV and Ta, 779 MeV at
3.49 GeV). The benchmark writes only events with exactly one cluster above it, so a different threshold needs a benchmark rerun: `--cluster-thr 100` appends
`--signal-min 100` to the chain job line (79 of 79 job lines = chain line + `--cluster-times --pulse-dump --signal-min 100`, rc 0; about 73 min wall at about
20 concurrent jobs including the WaveAnalyzer pass). Caches: `cache/<slug>_e100/`; benchmark outputs in the scratchpad (`ct_<slug>_e100`). More events have a
second cluster above 100 MeV than events gain a sole 100..500 MeV cluster, so N_ev drops. `tables/crosstalk_ratio.md` has a `cluster_thr_MeV` column
("500 x kE" = old rows).
| sample | N_ev 500 x kE | N_ev 100 MeV | non-seed hits (100) | (a) window | (b) window | h/ev | (a') core | (b') core | h/ev core |
|---|---|---|---|---|---|---|---|---|---|
| 0.73 GeV H2 | 2,893,353 | 2,713,040 | 10,638,306 | 0.32 % | 0.09 % | 0.003 | 99.15 % | 95.56 % | 3.747 |
| 0.73 GeV H2, MF time | 2,893,353 | 2,713,040 | 12,370,747 | 0.39 % | 0.09 % | 0.004 | 99.39 % | 95.19 % | 4.341 |
| 2.24 GeV H2 | 1,424,607 | 968,883 | 9,155,444 | 47.73 % | 6.28 % | 0.593 | 99.51 % | 84.16 % | 7.953 |
| 3.49 GeV H2 | 2,360,537 | 1,723,192 | 22,950,609 | 80.99 % | 10.24 % | 1.364 | 99.98 % | 79.39 % | 10.574 |
| Ta X17 25800 | 88,098 | 69,837 | 639,551 | 43.15 % | 5.78 % | 0.529 | 99.90 % | 81.25 % | 7.441 |
| Ta X17 26138 | 19,304 | 15,217 | 131,398 | 35.81 % | 4.88 % | 0.422 | 99.92 % | 82.10 % | 7.090 |
| Ta X17 26230 | 50,721 | 39,218 | 420,538 | 67.31 % | 9.60 % | 1.029 | 99.94 % | 70.54 % | 7.564 |
Command: `python3 crosstalk_study/crosstalk_study.py --run-type <slug> --cluster-thr 100 --selection single --no-raw --jobs 6 --work <scratch>/ct_<slug>_e100 --out crosstalk_study/figures/<slug>_single`.

Window hits per event that has a window hit: 1.02 (0.73), 1.28 (2.24), 1.74 (3.49). The 2.24 GeV X17 type does not exist (see top).

### 0.73 GeV H2 selection (same recipe, nothing adapted)
Cut-flow: 2,893,353 single-cluster events; tagged pair / kin / rej 9,246 / 595,388 / 595,388; not tagged 2,296,454; sigma_CB/E = 0.03416,
window |r - 1| < 0.1025 (r in [0.8975, 1.1025]): 1,080,125 events. Flags: (1) the chain's 0.73 GeV job lines use
`calibration_factor_3p5_June7_nonLinear.json` (3.5 GeV calibration; the benchmark prints `calib_match=mismatch`); they were used unchanged
(16 of 16 job lines identical to `final.txt` + the two flags), not the runinfo/general.json calibration. (2) The e-p window holds only 47 % of
the untagged events (mode of r at 0.99), the Moller kin tag removes 20.6 %. (3) The late (+7...+15 ns) population is absent at 0.73 GeV
(0.14 % of hits), unlike 2.24 and 3.49 GeV; a different regime (small clusters, 4.4 hits per event in the core), not a selection artefact as far as
checked.

## Numbers (legacy path unless noted; raw -> calibrated with time_calib.json; N_ev and hits as in the figure titles)
fraction |dt| > 10 ns [%]; late / early = number of hits with dt > +10 ns (module later) / dt < -10 ns; core = robust sigma / Gauss sigma [ns];
median chi2/dof (MF) for |dt| < 5 ns / |dt| > 20 ns (calibrated); amp = median non-seed/seed MF amplitude. The fraction above the fig2
chi2/dof limit is 0.98-1.0 % for every set by construction (limit = 99 % quantile of the window; 52-120 depending on set).

| sample | raw: >10 ns, late/early, core | cal: >10 ns, late/early, core | chi2/dof, |dt|<5 / |dt|>20 ns | amp |
|---|---|---|---|---|
| 2.24 e-p, CFD fit-only (5.51 M hits) | 4.75, 257,496/4,148, 3.31/2.37 | 4.67, 253,802/3,769, 1.58/0.68 | 1.33 / 1.47 | 0.0185 |
| 2.24 e-p, MF (6.32 M) | 1.00, 57,935/5,084, 2.91/2.44 | 0.74, 42,097/4,622, 1.27/0.56 | 1.36 / 1.39 | 0.0157 |
| 2.24 single, CFD fit-only (13.2 M) | 3.43, 438,317/14,205, 3.03/2.33 | 3.40, 435,835/13,504, 1.17/0.67 | 1.29 / 1.55 | 0.0258 |
| 2.24 single, MF (15.12 M) | 0.80, 103,178/17,523, 2.79/2.43 | 0.59, 72,072/16,808, 1.05/0.56 | 1.30 / 1.44 | 0.0218 |
| 2.24 e-p, CFD, gated | 4.65, 255,043/859, 3.30/2.37 | 4.57, 251,366/481, 1.57/0.68 | 1.33 / 4.4 (N=26) | 0.0185 |
| 3.49 e-p, CFD fit-only (16.57 M) | 6.52, 1,074,735/6,084, 3.78/2.50 | 6.50, 1,071,431/4,824, 2.23/0.67 | 1.41 / 1.54 | 0.0130 |
| 3.49 e-p, MF (18.66 M) | 1.26, 226,786/7,559, 3.14/2.61 | 0.89, 160,741/5,786, 1.48/0.58 | 1.44 / 1.42 | 0.0111 |
| 3.49 single, CFD fit-only (30.27 M) | 5.85, 1,756,276/15,531, 3.58/2.39 | 5.67, 1,703,480/13,696, 1.78/0.66 | 1.39 / 2.14 | 0.0161 |
| 3.49 single, MF (34.16 M) | 1.16, 376,210/21,171, 3.07/2.53 | 0.78, 248,317/18,826, 1.33/0.57 | 1.42 / 2.30 | 0.0137 |

Reading: with MF times the late population is small (0.9-1.3 % beyond +-10 ns raw); with WaveAnalyzer CFD times it is 4-7 %, a shoulder in
the crosstalk region dt = +7...+15 ns that the time calibration does not remove (it only shifts per-module offsets); the CFD core
narrows from 2.4 to 0.68 ns (Gauss) with the official table, i.e. the table is native to that basis. The long flat plateau at
dt = -58...+52 ns (about 10-100 hits / 0.5 ns, 8 ns modulation) is present in every set, also after calibration. The median module amplitude
falls from about 20 ADC at dt = +3 ns to 7-10 ADC at dt = +10...+18 ns (fig1b, e-p CFD), i.e. the late hits are the smallest pulses. The
fig3 wedge (seed amplitude below about 500 ADC with a module amplitude comparable to or above the seed, and a tail with seed above
3000 ADC) is present in the single-cluster samples and absent in the e-p samples, which keep seed amplitudes of about 600-3000 ADC.
