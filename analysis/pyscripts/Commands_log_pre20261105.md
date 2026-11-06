# Commands Log

Running log of commands for the PRad-II waveform analysis. Each entry has the
command and a concise explanation.

**Environment (run once per shell):**
```bash
export PYTHONPATH=$HOME/work/PRad/prad2evviewer/build/python
export PRAD2_DATABASE_DIR=$HOME/work/PRad/prad2evviewer/database
export DAQ_CONFIG=$PRAD2_DATABASE_DIR/daq_config.json
export HC_MAP=$PRAD2_DATABASE_DIR/hycal_map.json
```
`$PRAD2_DATABASE_DIR` is required for scripts that resolve default calibration paths (e.g. `module_dt_diagnostic.py`'s `--time-calib`, `--hycal-geometry` flags).

**Run → beam energy reference:**
| Run | Beam energy | Target | Notes |
|---|---|---|---|
| 024590 | 728.918 MeV (0.7 GeV) | — | Energy scan point |
| 025320 | 2239.51 MeV (2.2 GeV) | H2 | ep elastic, 60 nA |
| 024512 | 3488.43 MeV (3.5 GeV) | — | Energy scan point |
| 025308 | 2239.51 MeV (2.2 GeV) | Carbon | Carbon target |
| 026138 | 2239.51 MeV (2.2 GeV) | X17 | X17 physics, high pile-up |

---

## Template extraction

**Extract a per-material template for one run (physics population):**
```bash
PYTHONPATH=$HOME/work/PRad/prad2evviewer/build/python \
python3 fit_pulse_template.py \
    ~/work/PRad/data/evio/prad_<RUN>.evio.* \
    -o output/trigger_scan/pulse_templates_<RUN>_SSP_RawSum_h500_t0min22.0.json \
    --max-events 0 --max-pulses-per-channel 0 \
    --height-min 500 --model-err-floor 0.03 --t0-min 22.0 \
    --trigger-event-type SSP_RawSum \
    --daq-config $DAQ_CONFIG --hc-map-file $HC_MAP
```
Extracts the clean physics-population PbWO4 template. `--t0-min 22` selects
Population A (physics) over Population B (background). `--model-err-floor 0.03`
corrects the chi2 amplitude-scaling bias.

**Trigger + cut scan (parallel batch):**
```bash
bash run_trigger_scan.sh \
    --runs 024590,025320,024512 \
    --triggers SSP_RawSum \
    --height-mins 100,300,500,700 \
    --t0-mins 20.0,22.0,24.0
```
Runs the full grid of runs x triggers x height-cuts x t0-cuts in parallel.
Output: output/trigger_scan/pulse_templates_<RUN>_<TRIG>_h<H>_t0min<T0>p<D>.json

**Extended cut scan (low + high amplitude, wider t0 range):**
```bash
./run_scan.sh \
  --runs 024590,025320,024512 \
  --triggers SSP_RawSum \
  --height-mins 20,50,900,1200,1500 \
  --t0-mins 16.0,18.0,26.0,28.0
```
Extends the cut grid to probe extremes: very low amplitude (20, 50 ADC —
background-dominated), very high amplitude (900-1500 ADC — tests PMT
saturation shape distortion), and t0 cuts below (16, 18) and above (26, 28)
the physics window. Runs 024590/025320/024512 = 0.7/2.2/3.5 GeV beam energies.

---

## Diagnostics

**Peaks-per-channel distribution (three denominators + event-level pile-up):**
```bash
python3 npeaks_distribution.py \
    ~/work/PRad/data/evio/prad_025308.evio.00050 \
    ~/work/PRad/data/evio/prad_025320.evio.00099 \
    ~/work/PRad/data/evio/prad_026138.evio.00099 \
    --max-events 10000 --daq-config $DAQ_CONFIG \
    --out output/npeaks_distribution.json
```
Measures pile-up three ways: channel-weighted, pulse-weighted, event-level.

**Compare templates across the full cut grid:**
```bash
python3 compare_templates.py \
    --dir output/trigger_scan \
    --runs 024590,025320,024512 \
    --height-mins 20,50,100,300,500,700,900,1200,1500 \
    --t0-mins 16.0,18.0,20.0,22.0,24.0,26.0,28.0 \
    --material PbWO4 --plot \
    --out output/compare_templates_full.json
```
Compares PbWO4 tau_r/tau_f across all 9 height cuts x 7 t0 cuts x 3 energies
(0.7/2.2/3.5 GeV). Missing cells show "—". Produces vs-energy, vs-hmin, and
vs-t0min plots. Key questions: (1) does tau_r drop at high height-min
(900-1500 ADC) indicating PMT saturation? (2) is tau_r energy-independent
at fixed amplitude cut?

---

## Pile-up benchmarking

**Generate synthetic pile-up:**
```bash
python3 pileup_generator.py \
    ~/work/PRad/data/evio/prad_025308.evio.* \
    --template output/trigger_scan/pulse_templates_025308_SSP_RawSum_h500_t0min22.0.json \
    --out-dir output/pileup_synth_025308 \
    --n-per-config 50 --dt-grid 8,12,16,20,24,32,48,64,96,128 \
    --ratio-grid 0.1,0.2,0.3,0.5,1.0,2.0 --max-events 0 \
    --residual-veto-sigma 999 \
    --daq-config $DAQ_CONFIG --hc-map-file $HC_MAP
```
2-pulse synthetic pile-up on a (dT, ratio) grid. `--residual-veto-sigma 999`
effectively disables the residual veto (too strict for busy runs).

**Benchmark the deconvolver:**
```bash
python3 benchmark_deconv.py output/pileup_synth_025308 \
    --template output/trigger_scan/pulse_templates_025308_SSP_RawSum_h500_t0min22.0.json \
    --daq-config $DAQ_CONFIG --plot
```
Runs WaveAnalyzer.deconvolve on synthetic pile-up, compares recovered
amplitudes to truth.

**Benchmark the matched-filter peak finder:**
```bash
python3 benchmark_matched_filter.py output/pileup_synth_025308 \
    --template output/trigger_scan/pulse_templates_025308_SSP_RawSum_h500_t0min22.0.json \
    --daq-config $DAQ_CONFIG --plot-examples 10
```
Compares matched-filter detection vs. WaveAnalyzer.findPeaks on synthetic
pile-up.

---

## Timing gate (item 4)

**Legacy vs. timing-gated clustering:**
```bash
python3 benchmark_hycal_timing.py \
    ~/work/PRad/data/evio/prad_025308.evio.00050 \
    output/timing_gate_025308 \
    --window 8.0 --pre-window 100 200 --signal-min 500 \
    --max-events 50000 --e-range 500 3000 --calibrate --target-peak 2239 \
    --hc-calib-file $HOME/work/PRad/prad2evviewer/database/calibration/calibration_factor_2p2_fadc.json \
    --daq-config $DAQ_CONFIG --hc-map-file $HC_MAP
```
Compares HyCal cluster energy resolution with and without the seed-time gate.
`--window` is the timing gate half-width in ns.

---

## HyCal timing coincidence diagnostic (item 5)

**Raw (seed, neighbour, ΔT) extraction:**
```bash
cd $HOME/work/PRad/prad2evviewer/analysis/pyscripts

python3 study_hycal_timing.py \
    "$HOME/work/PRad/data/evio/prad_025320.evio.*" \
    output/hcal_dt_025320.tsv \
    --max-qdist 5.0 --pre-window 100 200 --require-trigger 0 \
    --hc-calib-file $PRAD2_DATABASE_DIR/calibration/calibration_factor_2p2_fadc.json
```
Dumps one TSV row per (seed, neighbour) pulse pair; 14 columns including `dt`, `dx_q`, `dy_q`, energies.
`--require-trigger 0` disables the default SSP0 strict-equality gate for an unbiased ΔT landscape.
`--hc-calib-file` override needed on run 25320 because runinfo chains it to the wrong (3.5 GeV) calibration; pass the 2.2 GeV one explicitly.
Quote the glob so Python expands it, not the shell. Output overwrites. Expect multi-GB TSV for a full run; use `--max-events N` for a smoke test.

**Per-module ΔT diagnostic (streaming consumer):**
```bash
python3 module_dt_diagnostic.py \
    output/hcal_dt_025320.tsv \
    output/dt_diagnostic_025320/ \
    --time-calib $PRAD2_DATABASE_DIR/hycal_time_offsets/time_calib.json
```
Streams the TSV in `--chunk-rows` batches; applies `time_calib.json` offsets Python-side (TSV `dt` is raw).
Produces: per-module fits PDF, 2D (module, ΔT) heatmap PNG+PDF, summary CSV, markdown flag report, and `hycal_secondary_peak_map.{png,pdf}` (three-panel: binary flag / μ₂ / significance).
Primary-fit flags: `low-stats`, `no-peak`, `offset-peak`, `wide-peak`, `sideband-overlap`. Informational: `secondary-peak`.
`--no-spatial-map` skips the HyCal spatial figure only; secondary detection always runs. Run `--self-test` first on a new install to confirm the fit and offset code paths.

**Self-test (sanity check):**
```bash
python3 module_dt_diagnostic.py --self-test
```
Three cases: synthetic Gaussian+flat (checks density-convention fit normalization), time-calib JSON parse (checks sign of `dt_corr = dt_raw - (off_nbr - off_seed)`), and synthetic double-peak (checks secondary-peak detection and flagging).
Expected output: `[self-test] Results: 3 passed, 0 failed.` with exit code 0.
