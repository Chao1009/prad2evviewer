# P12: readout threshold and below-threshold loss (U6; KTeV)

2026-11-13. Measurement only; no production threshold changed. Tools: `tools/p12_threshold.py` (raw decode, edge table), `tools/p12_table.py`
(loss table), `output/energy_resolution/variants/P12/run_p12.sh` (the 12 reports). Outputs: `output/energy_resolution/variants/P12/`
(`raw/`, `m1|m8|m12|m16/<run type>/`, `p12_loss_vs_cut.tsv`, `p12_loss_fit.tsv`, `logs/`). Reports: gates 8 and 12 ns, `EVAL_REGIONS=True`,
`CLUSTER_SET={'min_module_energy': c}`, on the existing H2 scans; the table uses path `g12`, default arm (the `none` arm peak from `legacy`).

## 1. Threshold from the data and the configuration

- The readout is zero-suppressed. In 40000 physics events per energy (run 24590 split 10, 25228 split 50, 24698 split 50) only 7.1 / 16.0 / 20.9
  PbWO4 channels per event are present (0.73 / 2.24 / 3.49 GeV) and every present channel carries the full 100 samples. No event has zero
  channels. Channels without a pulse above threshold are not read out.
- Configuration: `database/daq_config.json:90` firmware TET = 10 ADC above the per-channel pedestal, NSAT = 4; `database/daq_cnf/*.cnf` hold only per-channel FADC gains (`FAV3_ALLCH_GAIN`), no threshold.
- Data: pulse height (max minus median of the first 8 samples) of the present channels: 1st percentile 9.5 ADC at all three energies, median 24-26;
  1.0-1.3 % are below 10 (pedestal-estimate noise; 0.1-0.2 % below 8). Lower edge = the 10 ADC TET.
- Fit acceptance (BM:90, A >= max(5 sigma_A, 5 ADC)) is not binding: the matched-filter amplitude of stored hits has its 0.5 / 2 / 5 percentile at
  7.3-7.7 / 8.4-8.7 / 9.3-9.5 ADC (per-module medians over 1148 modules), below the TET; the energy edge is dominated by the hardware.
- MeV per module: k = median(h_e/h_a) = 0.95-0.97 MeV/ADC (16-84 % over modules 0.90-1.03). The TET is therefore **9.5-9.7 MeV median over modules**
  (0.73 / 2.24 / 3.49 GeV: 9.62 / 9.49 / 9.67), 16-84 % 8.9-10.4 MeV, full range 5.3-13.3 MeV. Stored hit energy, 2 % percentile per module:
  8.0-8.3 MeV median. A fixed ADC threshold is a 5-13 MeV threshold across modules.

## 2. Loss against cluster cut (cut c applied on top of the readout threshold, 1 MeV = control, no effect below the stored 6 MeV)

res_cb = Crystal Ball sigma/E (%), g12, default arm; edge/centre = `pool_edge` / `pool_centre`; mu_none = legacy CB peak / E' - 1 (%) before gains
(mean over splits, SE of that mean 0.006-0.28).

| run type | c (MeV) | res_cb | edge | centre | tail_lo3 (%) | edge | centre | mu_none (%) |
|---|---|---|---|---|---|---|---|---|
| 0.73 | 1 | 3.370 | 3.440 | 3.255 | 44.8 | 46.0 | 43.6 | -1.73 |
| 0.73 | 8 | 3.378 | 3.438 | 3.275 | 44.7 | 45.9 | 43.5 | -1.90 |
| 0.73 | 12 | 3.525 | 3.535 | 3.473 | 44.1 | 44.7 | 43.4 | -3.24 |
| 0.73 | 16 | 3.800 | 3.768 | 3.795 | 42.8 | 42.9 | 43.0 | -4.94 |
| 2.24 | 1 | 1.864 | 1.915 | 1.832 | 27.9 | 29.9 | 26.7 | 2.13 |
| 2.24 | 8 | 1.869 | 1.920 | 1.835 | 27.9 | 30.0 | 26.7 | 1.95 |
| 2.24 | 12 | 1.936 | 1.987 | 1.901 | 28.6 | 31.0 | 27.3 | 0.86 |
| 2.24 | 16 | 2.011 | 2.067 | 1.978 | 29.5 | 33.2 | 27.1 | -0.26 |
| 3.49 | 1 | 1.571 | 1.645 | 1.527 | 11.1 | 13.3 | 9.8 | 4.18 |
| 3.49 | 8 | 1.575 | 1.648 | 1.530 | 11.1 | 13.3 | 9.8 | 4.07 |
| 3.49 | 12 | 1.613 | 1.692 | 1.565 | 11.0 | 13.2 | 9.7 | 3.23 |
| 3.49 | 16 | 1.669 | 1.752 | 1.626 | 10.8 | 13.0 | 9.5 | 2.26 |

c = 1 and 8 agree within 0.004-0.008 pp (the readout edge sits above 8 for most modules); the response falls and the width rises from c = 12.
Edge-minus-centre gap of res_cb stays 0.08-0.12 pp (2.24, 3.49 GeV) at every c; at 0.73 GeV it closes (0.19 to -0.03 pp): the loss acts at the
centre there. tail_lo3 at the edge rises with c at 2.24 GeV (+3.3 pp from c = 1 to 16), not at the centre (+0.4 pp).

## 3. Size of the loss at the readout edge (model-dependent)

The cut lifts the effective threshold to t = max(c, T_HW). Fit on t = T_HW (c = 1 and 8 pooled), 12, 16, extrapolated to 0 (`p12_loss_fit.tsv`):
res^2(t) = r0^2 + b t^2, mu_none(t) = mu0 - s t.

| run type | T_HW (MeV) | res_cb ctrl | r0 (zero threshold) | loss term at T_HW (pp, quadrature) | loss / res_cb | same if variance ~ t (t^3) | mean loss at T_HW |
|---|---|---|---|---|---|---|---|
| 0.73 | 9.62 | 3.37 | 3.12 | 1.31 | 0.39 | 0.64 (0.27) | 4.6 % = 34 MeV |
| 2.24 | 9.49 | 1.86 | 1.80 | 0.54 | 0.29 | 0.48 (0.20) | 3.3 % = 74 MeV |
| 3.49 | 9.67 | 1.57 | 1.52 | 0.42 | 0.27 | 0.44 (0.18) | 2.8 % = 98 MeV |

**Verdict: visible under the extrapolation; robust only at 0.73 GeV on the stricter reading.** Reading used: loss term / res_cb (amplitude ratio) >= 0.10,
which holds at all three energies for every variance power tried (0.18-0.64). On the variance share (that ratio squared) it is 0.41 / 0.15 / 0.07 at 0.73 GeV
and 0.23 / 0.08 / 0.04 at 2.24 GeV and 0.19 / 0.07 / 0.03 at 3.49 GeV for variance ~ t, t^2, t^3, so at 2.24 and 3.49 GeV only the t and (just) t^2 models pass.
Largest at 0.73 GeV as expected (U6).

Caveats. (a) The extrapolation from t = 9.5-16 MeV to 0 is unconstrained by any data point. (b) The cut-scan rise in res^2 also contains changes to the
log-weighted position and to the refitted gains (a 12 MeV module still carries weight at 0.73 GeV), so it is an upper bound on the loss term. (c) No tail
power fits the onset: the c = 12 point lies above every fit, so t = max(c, median T_HW) is only approximate (module edges span 5-13 MeV). (d) The mean loss of
3-5 % is a linear extrapolation; a flat hit density near threshold would make it grow like t^2, so it is an upper bound too. (e) NSAT = 4 is not visible in
the data: the peak-height edge (9.5 ADC) is that of a single-sample 10 ADC crossing. Not measurable here: loss below the hardware edge itself.

## 4. Consequence

A KTeV-style below-threshold correction meets the plan's criterion and becomes an `Amend parent:` candidate: a per-module (not per-ADC) threshold from
the 1148-module edge table (`raw/edge_*.npz`), an expected-loss term from the neighbours' energies, tested in P1-style paired reports. Before that,
the extrapolation should be tightened: a scan with a lower effective hardware threshold does not exist, so the only cross-check is Geant4-free
module-wise (select modules with edge < 8 vs > 11 MeV and compare res_cb and mu_none; not run, outside this block).
