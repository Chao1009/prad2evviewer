# crosstalk_study

Matched-filter chi-squared per degree of freedom of the seed and neighbour modules of HyCal clusters, in
slices of cluster energy, for the three H2 run types (0.73, 2.24, 3.49 GeV), plus the constant-fraction-discriminator
time of every fitted channel. One standalone script, `crosstalk_study.py`, read straight from the EVIO files; the
scripts, figures, tables and caches of the previous study (2026-11-06/07) are kept under `previous_20261107/`,
`figures/`, `tables/`, `cache/` and `tcal/` and are not used.

## Chain

`scan <out.npz> <evio>...` (one run, all its splits):

1. Events: physics trigger gate of the run type (`bits == 0x100`; X17: `bits & 0xf00 != 0`) and light-monitoring veto
   (`bits & 0x1000000`). No other event selection.
2. Channels: every PbWO4 module with 100 samples, except the run's dead modules (`database/runinfo/general.json`,
   `hycal_dead_modules`). Lead glass is not read.
3. Bad waveforms (any sample equal to 0 or more than 20 ADC below the baseline) are dropped and counted per module
   (`mod_n_bad`). Pre-cut (`PRECUT_ADC`): raw maximum minus the baseline (median of the first 30 samples) >= 8 ADC. Noise reference (per file, from its first 20 000 events): per module, sigma_ref (the reference noise sigma: 1.4826 times the median absolute deviation (MAD) of the residuals of the waveforms about their own baselines) and b_ref (the median of those baselines); a module with fewer than 20 reference waveforms takes the median sigma_ref of the other modules in the file and no b_ref. b_ref replaces the start of the baseline fit when the 30-sample baseline is more than 3 sigma_ref away.
4. Matched filter, stages A (candidates), B (joint fit) and C (shoulders), with the run type's template bank
   `waveform_analysis/pulse_templates/<run type>/templates_<run type>.npz` and its sigma(t) shape scatter. The baseline offset is a fit parameter, clipped to ±1 sigma_ref (as in `pileup_rates.py`). A channel is
   kept if its largest fitted amplitude (`KEEP_ADC`) is >= 10 ADC (about 9.5 MeV at the median PbWO4 gain of 0.95 MeV per ADC,
   16-84 % of modules 8.9-10.3 MeV; 8.7 MeV for X17).
   A kept-size channel whose main pulse has its 32-sample fit window cut by the end of the waveform is dropped and
   counted per module (`mod_n_trunc`).
5. Stored per kept channel: module, event, number of pulses k, amplitudes and onsets of all pulses (time order; the main
   pulse is the largest), chi-squared per degree of freedom, pedestal, pedestal sigma, constant-fraction-discriminator
   time of the main pulse, number of clipped samples. Per module: the template's 50 % rise time after onset. Per
   cluster-energy slice: up to 12 random seeds with chi-squared per degree of freedom > 5, with raw waveform and fit.

`report <scan dir> <out dir>` (all `run_*.npz`, pooled per run type):

1. Every cluster of the event. In-window pulses have onset in 130-200 ns; a seed must be a local maximum (no larger
   in-window channel in its 3 x 3 neighbourhood); seed = the free (not yet clustered) local maximum with the largest
   in-window pulse; members = seed and the free channels within 1.5 x the mean module size (the 3 x 3 neighbourhood)
   with an in-window pulse; repeated until no free local maximum is left. The first cluster follows the same
   seed-and-neighbours rule as the leading cluster of `waveform_analysis/pileup_analysis/pileup_rates.py`. All
   clusters are used, not only the leading one: under the physics trigger the leading cluster is almost never below 100 MeV (30 000-event test: 0, 2
   and 0 leading clusters in 10-50 MeV at 0.73, 2.24, 3.49 GeV), so the low slices would be empty.
2. Cluster energy = sum over members of gain x A, A = the member's largest in-window pulse, gain (MeV per ADC) =
   calibration factor of the run (runinfo calibration file) x `iint_wa` of the module's template. Members without a
   calibration factor are left out; clusters whose seed has none are dropped (counted in `summary.txt`).
3. Slices of cluster energy (MeV): [10, 50), [50, 100), [100, 200), [200, 500), [500, 1000), [1000, 2000),
   [2000, 3000), [3000, 4000).

Definitions:

- Chi-squared per degree of freedom (as `pileup_rates.py`): sum of squared residuals over the union of the pulses'
  32-sample fit windows (clipped samples excluded), divided by the sigma(t) variance = per-file reference noise sigma_ref squared (not the channel pedestal sigma) + 1/12 +
  (sum of A_k x template shape scatter x sqrt(a_band / max(A_k, 20 ADC)))^2; degrees of freedom = samples - 2k - 1 (the baseline offset is a fit parameter).
- Constant-fraction-discriminator time: the largest pedestal-subtracted sample within 12 ns of the main pulse's template
  peak, 50 % of it, linear interpolation between the last sample below 50 % and the next one on the rising edge. Not
  computed for clipped channels.
- Matched-filter 50 % time t_mf = onset of the main pulse + the template's 50 % rise time, so a pulse that matches its
  template gives t_cfd - t_mf near 0.

## Constants

Constants are at the top of `crosstalk_study.py`, each with a one-line comment; `CROSSTALK_MAX_EVENTS=<n>`
(environment) stops the scan after n decoded events, for tests.

The matched filter (class `MF` and its constants) is copied verbatim from
`waveform_analysis/pileup_analysis/pileup_deconvolution_study.py`. Class `Bank` is copied from `pileup_rates.py`, the
data version of the study's class (reads the run-type bank and uses the bank's average template for modules without
their own); the study's own `Bank` needs every module to have its own template and stops on the 26-50 PbWO4 modules per
bank that do not. `run_config`, `module_gain`, `neighbours` and the decoder loop are those of `pileup_rates.py`.

## Outputs (`out/`)

| file | content |
|---|---|
| `seed_chi2_vs_cluster_energy.png` | the requested figure: one row per H2 run type, one panel per slice, histogram of the seed's chi-squared per degree of freedom (log x, 0.1-100), entries, median (solid) and 90th percentile (dashed) |
| `chi2_summary_vs_cluster_energy.png` | median and 90th percentile vs slice, three run types; seed and neighbour (non-seed member) panels |
| `chi2_slices.tsv` | entries, median and 90th percentile per run type, seed or neighbour, slice (every scanned run type, X17 included) |
| `chi2_vs_dt_cfd_mf.png`, `.tsv` | median and 16-84 % band of chi-squared per degree of freedom vs t_cfd - t_mf (0.5 ns bins), seeds and neighbours of 10-4000 MeV clusters |
| `seed_vs_neighbour_chi2.png` | 2D histogram, seed chi-squared per degree of freedom vs that of each neighbour |
| `examples_seed_chi2_gt5_2p24GeV_H2.pdf` | one page per slice, 6 random seeds with chi-squared per degree of freedom > 5: raw samples, pedestal + fit, matched-filter onsets, constant-fraction-discriminator time |
| `summary.txt` | per run type: runs, splits, scan time, events, clusters, dropped clusters, constant-fraction-discriminator coverage |
| `late_component.pdf`, `late_component.tsv` | late constant-fraction component study: partner-channel cross-talk (produced by one-off analyses, not by the script) |
| `seed_chi2_two_peaks_10_50MeV.png`, `seed_chi2_two_peaks_10_50MeV.tsv` | seed chi-squared per degree of freedom of 10-50 MeV clusters, noise-sigma quantisation check (produced by one-off analyses, not by the script) |

Percentiles come from histograms with 200 log bins per decade (0.01-10^4), interpolated inside a bin.

## Commands

From `analysis/pyscripts` after `source setup_env.sh`:

```
python3 crosstalk_study/crosstalk_study.py scan output/crosstalk_study/scan/run_25228.npz /path/prad_025228.evio.*
python3 crosstalk_study/crosstalk_study.py report output/crosstalk_study/scan crosstalk_study/out
bash crosstalk_study/run_crosstalk_scan.sh          # all H2 runs, NPROC = 6; X17 runs listed, add them to RUNS to include
```

Test on part of a run: `CROSSTALK_MAX_EVENTS=30000 python3 crosstalk_study/crosstalk_study.py scan ...`.

## Previous study (2026-11-06/07; `previous_20261107/`)

Inputs: the frozen benchmark `benchmark_hycal_timing_test_20261021.py` rerun with the r5ri chain's job lines plus
`--cluster-times --pulse-dump`, a WaveAnalyzer pass with an uncommitted log-normal fit binding, and
`database/hycal_time_offsets/time_calib.json`. Run types 0.73 GeV H2 (24590, 24603), 2.24 GeV H2 (25228), 3.49 GeV H2
(24386, 24503, 24512) and 2.24 GeV Ta X17 per run (25800, 26138, 26230). Findings: `notes/crosstalk_findings_20261107.md`.

What it looked at:

- Selections: e-p (single cluster above 500 x beam energy in GeV MeV or 100 MeV, no Moller tag, |E / E_ep - 1| < 3
  sigma) and single-cluster; a fired-module cut n_mod < 20 as a variant.
- dt = T_module - T_seed of every non-seed cluster member, with constant-fraction-discriminator (WaveAnalyzer
  log-normal fit) or matched-filter times, with and without per-module time offsets; crosstalk window +7...+15 ns and
  core |dt| < 2 ns; event and hit ratios in both (`tables/crosstalk_ratio.md`).
- Amplitude of window and core hits, seed vs module amplitude (fig3, fig4).
- dt vs matched-filter chi-squared per degree of freedom (fig2, `summary_dtcfd_vs_mfchi2.py`), its median and 16-84 %
  profile (`profile_mfchi2_vs_dtcfd.py`).
- 50 examples of core hits with chi-squared per degree of freedom > 5 (`examples_mfchi2_gt5.py`).
- Seed chi-squared per degree of freedom vs non-seed amplitude and dt (`seed_chi2.py`).
- Self-referenced time offsets compared with the official table (`tcal/`).

Kept, re-implemented on the new chain:

- Chi-squared per degree of freedom vs a timing difference, now per hit t_cfd - t_mf (the discriminant the findings
  note proposed) instead of dt to the seed, seeds and neighbours separately.
- Seed vs neighbour chi-squared per degree of freedom.
- Examples with chi-squared per degree of freedom > 5, now seeds per cluster-energy slice, with the fit drawn.
- New: seed and neighbour chi-squared per degree of freedom per cluster-energy slice.

Dropped, and why:

- dt to the seed, the crosstalk-window and core ratios, and the amplitude figures: they need the benchmark cluster
  reconstruction and the per-module time offsets of the old chain; the new chain has neither (onsets are template
  onsets, no time calibration), and the request is about chi-squared vs cluster energy.
- e-p and Moller selections: they need the benchmark's elastic kinematics and Moller tags; the new report uses every
  cluster with no kinematic selection, so the low-energy slices are populated.
- The WaveAnalyzer log-normal fit time: it needs the uncommitted `python/bind_dec.cpp` change; replaced by a 50 %
  constant-fraction-discriminator time computed from the samples.
- The n_mod < 20 variant and the self-referenced time offsets: tied to the old cluster and time-offset inputs.
