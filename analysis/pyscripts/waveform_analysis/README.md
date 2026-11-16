# waveform_analysis

Run every command below from `analysis/pyscripts` after `source setup_env.sh`.

## Scripts

- `template_extract.py` builds per-module sampled (spline) PbWO4 pulse templates from the pulses that pass the quality cuts, pooled per module: `select` turns the EVIO splits of one run into one candidate npz that stores every such pulse, and `build` pools the pulses with e_rough >= 0.5 E_beam per module (`POOL_EFRAC`; at most `POOL_MAX` = 1000 at random); the pools are tried in order and the first whose fit succeeds is used: the elastic pool if it has at least `POOL_MIN_EP` = 20 pulses, the pulses with amplitude `A_LO` = 200 to `A_HI` = 1000 ADC if at least `POOL_MIN_GEN` = 20 (`pool_general`), the elastic pool of 3 to 19 pulses, the single largest pulse, the run-type average; `build` fits the templates for one run type from its runs' npz files. Run types: `0p73GeV_H2`, `2p24GeV_H2`, `2p24GeV_C`, `2p24GeV_Ta_X17`, `3p49GeV_H2`. `run_template_select.sh` runs `select` for every run of the four run types (8 runs in parallel; `pulse_templates/<RT>/pulses/run_<run>.npz`, logs in `pulse_templates/<RT>/logs/`), then `build` for each run type.

```bash
python3 waveform_analysis/template_extract.py select <out.npz> <run type> <evio file> [<evio file> ...]
python3 waveform_analysis/template_extract.py build <run type>
nohup waveform_analysis/run_template_select.sh > output/run_template_select.log 2>&1 &
```

- `cmp_templates.py` computes the 10-90% rise and 90-10% fall times per module from the template banks (`waveform_analysis/pulse_templates/<RT>/templates_<RT>.npz`) and writes tables and distributions to `waveform_analysis/pulse_templates/timing/` (`timing_summary`, `timing_modules`, `timing_maxdiff` including the rise/fall standard deviation over run types (at least 3 run types), `timing_sigma`, `sigma_modules`, `jitter_modules`; figures `timing_dist` and `timing_sigma`, each as png and pdf) and HyCal maps to `waveform_analysis/figures/` (`timing_map_<RT>`, `timing_map_all`, `timing_maxdiff_map`, `timing_std_map`, `sigma_map_<RT>`/`_all`, `jitter_map_<RT>`/`_all`, `rise_jitter_map_<RT>`/`_all`, `fall_jitter_map_<RT>`/`_all`).

```bash
python3 waveform_analysis/cmp_templates.py
```

- `pulse_shape_study.py` classifies the tails of good pulses as smooth, plateau or bumpy (a dip followed by a secondary rise >= max(3 ADC, 3 % of A); plateau is not counted as a bump), and writes the PbWO4-only bumpy-rate heat map, `bumpy_rate_hist.png`, `examples.pdf`, `modules.tsv`, `hot_modules.tsv` (bumpy rate > 10 % AND >= 20 pulses), `hot_examples.pdf` (2 smooth + 2 bumpy pulses per hot module) and `summary.txt`; the report keeps A >= 200 ADC.

`run_pulse_shape_scan.sh` scans every run of the four run types (one npz per run, `output/pulse_shape_study/scan/run_<run>.npz`, 8 runs in parallel) and then writes the report to `output/pulse_shape_study/`.

```bash
python3 waveform_analysis/pulse_shape_study.py scan <out.npz> <evio file> [<evio file> ...]
python3 waveform_analysis/pulse_shape_study.py report <scan dir> <out dir>
nohup waveform_analysis/run_pulse_shape_scan.sh > output/run_pulse_shape_scan.log 2>&1 &
```

Each driver sets `DATA_DIR`, `OUT_DIR`, `NPROC` and the run lists at its top. One driver job is one Python call per run, so on the JLab ifarm `DATA_DIR` changes (and the paths in `setup_env.sh`), and a batch system can submit one job per run.

- `energy_resolution_benchmark.py` compares the PbWO4 energy resolution of legacy and time-gated island clustering (standalone: NumPy, SciPy, Matplotlib and `prad2py`). `scan` stores the matched-filter hits of one run without clustering; `report` clusters them for the legacy path and each time gate given on the command line and writes tables and a PDF to `energy_resolution/<RT>/` (recipe and outputs in `energy_resolution/README.md`). `run_energy_resolution.sh` scans every run of the five run types (the four above plus `3p49GeV_H2_old`, the previous 3.49 GeV runs; 8 runs in parallel), then runs `report` per run type with its `GATES`; a new gate list needs only `report`.

```bash
python3 waveform_analysis/energy_resolution_benchmark.py scan <out.npz> <evio file> [<evio file> ...]
python3 waveform_analysis/energy_resolution_benchmark.py report <scan dir> <gate_ns> [<gate_ns> ...]
python3 waveform_analysis/energy_resolution_benchmark.py summary waveform_analysis/energy_resolution   # summary.pdf, <RT>/runs.tsv
nohup waveform_analysis/run_energy_resolution.sh > output/run_energy_resolution.log 2>&1 &
```

- `compare_old_vs_new_3p49.py` is a one-off comparison of the old 3p49GeV_H2 runs (24386, 24503, 24512) with the new ones; it reads `pulse_templates/3p49GeV_H2` and `pulse_templates/3p49GeV_H2_old` and writes `pulse_templates/3p49GeV_H2_old/compare_old_vs_new.pdf` and `compare.txt`.

```bash
python3 waveform_analysis/compare_old_vs_new_3p49.py
```

### pileup_analysis/

- `pileup_analysis/pileup_generator_stage1.py` writes fully synthetic 2.24 GeV H2 two-pulse waveforms (100 samples, uint16) for five second-ring PbWO4 modules from their own spline templates, cycling over main amplitude (300, 1000, 3000 ADC), secondary/main amplitude ratio (0.05-2) and signed time gap dt = t2 - t1 (symmetric grid, -128 to -8 and 8 to 128 ns; for dt < 0 the secondary arrives first, the main pulse stays the A1 one at t1), with exact truth for matched-filter deconvolution tests; each event's template is smeared by the bank's `shape_sigma` scaled by sqrt(`a_band` / A1) (the measured 1 / sqrt(amplitude) law; `SMEAR_K = 1`; one smooth random spline perturbation per event shared by both pulses, rms-calibrated per module to `shape_sigma`, stored as `smear_k` and per-event `smear_rms`); no waveform reaches 4095 ADC (grid combinations whose noiseless peak exceeds 4095 - 4 sigma are dropped and listed at startup, any event that still reaches 4095 is redrawn); output defaults to `output/pileup_stage1/pileup_stage1.npz`.

```bash
python3 waveform_analysis/pileup_analysis/pileup_generator_stage1.py --n-events 19510000 output/pileup_stage1/pileup_stage1.npz
```

- `pileup_analysis/pileup_deconvolution_study.py` is standalone (NumPy, SciPy, Matplotlib only): it evaluates the generator's own spline templates exactly as `pileup_generator_stage1.py` does and runs its own three-stage matched filter on a stage-1 file: a one-pass peak search on the smoothed matched-filter amplitude (a TSpectrum-like search without ROOT), a joint fit of all candidates, then residual shoulders added one at a time and kept only when the BIC improves (at most 6 pulses, `MIN_SEP_NS` = 8 ns; the module docstring's Algorithm has the details). It counts the secondary as found when a fitted pulse other than the one nearest t1 lies within 4 ns of t2 with amplitude within 30 % of A2. It writes `efficiency.tsv` (per module, A1, ratio and gap, with binomial errors, the r = 0 false-positive rate, resolutions and the median chi2/dof against the generator's noise sigma, the main-pulse amplitude pull and `ped_bias`, the median (fitted - true pedestal) / sigma, which shows the effect of an early secondary in the pedestal window), `efficiency_tables.pdf` (one ratio x gap table of efficiency ± binomial error per module pooled over A1, then split by A1) and `pileup_examples.pdf` (example fits over a ratio x gap grid, negative and positive dt, for one module and A1, set by the `EX_*` constants) to its own directory, `waveform_analysis/pileup_analysis/` (8 worker processes).

```bash
python3 waveform_analysis/pileup_analysis/pileup_deconvolution_study.py output/pileup_stage1/pileup_stage1.npz
```

- `pileup_analysis/pileup_rates.py` is standalone (NumPy, SciPy, Matplotlib and the `prad2py` decoder; database files from `$PRAD2_DATABASE_DIR`, which `setup_env.sh` exports) and measures the pile-up rate in the data for the four run types (`0p73GeV_H2`, `2p24GeV_H2`, `3p49GeV_H2` and `2p24GeV_Ta_X17`; other runs are skipped), one run per scan. It runs the same three-stage matched filter as `pileup_deconvolution_study.py` (an identical `MF` class) with the run type's own spline bank on the HyCal channels that pass the bad-waveform flag (any sample equal to 0 or more than 20 ADC below the baseline; removed and counted per module, `bad_count`) and have raw max - median(first 30 samples) >= 8 ADC (`PRECUT_ADC`; only these are fitted) of the events that pass the trigger gate and the LMS veto. A channel counts (is kept) if its largest fitted amplitude is >= 10 ADC (`KEEP_ADC`, about 9.5 MeV at the median gain), and it counts as pile-up if the fit finds two or more pulses. The fit includes a baseline offset (`dped`, clipped at ±1 sigma (one standard deviation of the pedestal noise); stored, its distribution is in `summary.txt`); a kept channel whose main pulse is truncated by the end of the waveform (`trunc_main`) stays in the raw k counts and k maps but is excluded from every E_main-based quantity. `report` quotes the primary rates with noise-suppressing cuts: a channel is pile-up (cut) when E_main = A_main x gain >= 0.1 x E_beam of its run (E_beam from `database/runinfo/general.json`, per run: 3485.4 MeV for the new and 3488.4 MeV for the old 3p49GeV_H2 runs; gain = calibration factor of the run from `database/runinfo/general.json` -> `database/calibration/*.json` x I_int (`iint_wa`) of the module's spline template, as computed in `pileup_rates.py`; modules without a factor, including lead glass, are left out and counted in `summary.txt`) and another fitted pulse has A_other / A_main > 0.1 and |dt| >= 8 ns, the channel rate is over kept channels with E_main >= 0.1 E_beam, the cluster rate is over leading clusters with seed E_used >= 0.1 E_beam (a member at E_used >= 0.1 E_beam with such an extra pulse; all extras and in-window extras), and the per-module rates and maps (rate_cut, modules with >= 100 channels at E_main >= 0.1 E_beam) use the same cut; a channel that is the only kept channel in its event is vetoed from all cut rates, maps, medians and clusters (in a few modules these lone channels are self-triggered bursts, not pile-up) and counted in `n_single_veto` and `n_cut_single_veto`; channels of events with 300 or more kept channels (array-wide pulses) are vetoed the same way (`n_array_veto`, `n_cut_array_veto`); saturated channels (3 or more clipped samples, i.e. raw samples >= 4095 ADC) and bad fits (chi-squared per degree of freedom > 1000) are excluded like bad waveforms (from the cut rates, clusters, slices, medians, maps and examples) and counted per module only; the uncut k >= 2 rates stay as comparison columns. It writes `pileup_rates.tsv` and `summary.txt` (channel-level and cluster-level rates with binomial errors, plus the clipped rates: kept channels and clusters with any sample at 4095 ADC, and the pile-up fraction among clipped channels; the cluster is one leading cluster per event, a seed plus the modules within 1.5 x the mean module size, using only pulses with onset in the (130, 200) ns pre-window), `pileup_grid_<RT>.tsv` (two-pulse channels per time-gap x amplitude-ratio cell), `pileup_modules_<RT>.tsv` (per-module raw rates for modules with > 100 kept channels, `rate_k2_raw` = N(k >= 2) / N(kept before the vetoes), `rate_k2` ... `rate_k6` = N(k >= n) / N(kept after the vetoes), cut rates for modules with >= 100 channels at E_main >= 0.1 E_beam; column `n_emain`), `pileup_rate_map_<RT>.png`, `pileup_rate_slices_<RT>.png` and `pileup_slices.tsv` (the slice rate, the cut rate in 6 slices of E_main / E_beam, 0.1-0.3 ... 1.1-1.3: N(cut channels in the slice) / N(all kept channels of the module after vetoes and exclusions) per module, drawn for modules with >= 10 channels in the slice, and the pooled rate over the run type with N(all kept channels of the run type) as denominator; one shared logarithmic colour scale) and `pileup_rates.pdf` (tables, PbWO4 maps, 20 pile-up examples, one per time-gap x ratio cell as far as populated cells allow, and 10 no-pile-up example fits; `examples.pdf` (the X17 beam-hole ring page, 16 spread examples per run type on 2 pages, the excluded-channels map page (bad waveforms, saturated and bad fits per module) and the N(k >= n) / N(kept channels after the vetoes) map pages for n = 2 ... 6, k >= 2 first, 5 pages) and `hot_modules.pdf` (the 8 modules per run type with the highest cut rate: a summary table, then one page per module with 6 pile-up and 2 one-pulse panels)) to `waveform_analysis/pileup_analysis/pileup_rates_out/`. `pileup_analysis/run_pileup_scan.sh` scans every run of the four run types (one npz per run, `output/pileup_rates/scan/run_<run>.npz`, 8 runs in parallel) and then runs the report.

```bash
python3 waveform_analysis/pileup_analysis/pileup_rates.py scan <out.npz> <evio file> [<evio file> ...]
python3 waveform_analysis/pileup_analysis/pileup_rates.py report <scan dir> <out dir>
nohup waveform_analysis/pileup_analysis/run_pileup_scan.sh > output/run_pileup_scan.log 2>&1 &
```

## Notes

Baseline alternative (not implemented; considered as an alternative to the free baseline offset of the joint fit): take, per module and per file, a reference baseline (the median of the per-event 30-sample baselines) and allow the fitted baseline to differ from it by at most 3 sigma of the pedestal noise; where the 30-sample window is itself contaminated, fall back to the lowest-mean 15-sample window of the waveform.
