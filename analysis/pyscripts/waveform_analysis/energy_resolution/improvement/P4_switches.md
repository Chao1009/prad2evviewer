# P4: pool attribution switches (rows a, b)

Control: row a against C0 (`baseline_20261112/c0/2p24GeV_H2`, run 25320, 8 splits); row b against `baseline_20261112/pre/<RT>`.
Variants: `output/energy_resolution/variants/P4/{a_phase_off,b_logw_r7ri}/<RT>/` (`delta.tsv` per D3; `table.py` prints the rows below).
Δ = variant − control, mean ± SE over paired splits, ep pool, gate 12 ns (legacy in the TSV). `bank_sha1`: a 7f125b75 (rescan and C0 equal); b the existing scans' banks (0.73: 5d588cce, 2.24: 15dfd867, 3.49 and 3.49_old: 294184b0, same on control and variant); override strings in the first line of each `delta.tsv`.
Units: res_cb, tail, mu in pp or %, leak_cut and eff_cut in %.

## a: `PHASE_CORR=False` (rescan of run 25320, `scan_variants/phase_off/2p24GeV_H2/`)

| path | Δres_cb | Δres_cb oos | Δtail_hi3 | Δleak_cut | Δeff_cut | Δtail_lo3 | Δmu | yield loss % | verdict |
|---|---|---|---|---|---|---|---|---|---|
| legacy | +0.058 ± 0.007 | +0.055 ± 0.007 | −0.006 ± 0.004 | +0.433 ± 0.036 | −0.061 ± 0.036 | −0.180 ± 0.023 | −0.012 ± 0.007 | −0.19 | res_cb real, leak real |
| g12 | +0.063 ± 0.007 | +0.062 ± 0.009 | −0.003 ± 0.003 | +0.440 ± 0.039 | −0.060 ± 0.028 | −0.198 ± 0.032 | −0.014 ± 0.007 | −0.18 | res_cb real, leak real |

Reading: removing the phase correction makes the core 0.06 pp wider (g12 1.863 → 1.926 %) and leak_cut worse by 0.44 %; the phase correction is a gain, not the −0.04 pp "production only" effect expected. The tail_lo3 fall (−0.2) is a guard that moves in the favourable direction. Assumption check (mean energy change < 0.05 %): Δmu = −0.012 ± 0.007 %, holds. P10 gate (|Δ| ≥ 2 SE and ≥ 0.02 pp): met; single run, so by D5 it is extended to the other 2.24 GeV runs and 3.49 GeV before any verdict (the gate decides P10 starts; the verdict on the switch is "keep PHASE_CORR=True"). Not a candidate for adoption: it is the production default.

## b: `CLUSTER_SET` `log_weight_thres` = R7ri value (0.73 / 2.24 / 3.49 / 3.49_old: 3.4101 / 3.6576 / 3.7073 / 3.7073) vs 5.2

| RT | path | n | Δres_cb | Δres_cb oos | Δtail_hi3 | Δleak_cut | Δeff_cut | Δtail_lo3 | Δmu | yield loss % |
|---|---|---|---|---|---|---|---|---|---|---|
| 0p73GeV_H2 | g12 | 16 | +0.006 ± 0.003 | −0.003 ± 0.003 | −0.001 ± 0.001 | +0.000 ± 0.000 | −0.004 ± 0.010 | −0.049 ± 0.026 | +0.002 ± 0.004 | 0.36 |
| 2p24GeV_H2 | g12 | 42 | −0.011 ± 0.002 | −0.010 ± 0.002 | −0.001 ± 0.002 | −0.098 ± 0.013 | +0.059 ± 0.009 | +0.014 ± 0.009 | −0.004 ± 0.002 | 0.09 |
| 3p49GeV_H2 | g12 | 16 | −0.026 ± 0.003 | −0.034 ± 0.004 | −0.013 ± 0.003 | −0.371 ± 0.032 | +0.128 ± 0.022 | +0.069 ± 0.017 | +0.007 ± 0.002 | 0.16 |
| 3p49GeV_H2_old | g12 | 28 | −0.033 ± 0.003 | −0.033 ± 0.004 | −0.012 ± 0.004 | −0.501 ± 0.055 | +0.120 ± 0.014 | +0.126 ± 0.023 | +0.012 ± 0.003 | 0.17 |

Legacy path agrees in sign and size (3.49: res_cb −0.028 ± 0.003, leak_cut −0.38 ± 0.05; 3.49_old −0.036, −0.56). The 0.73 GeV change is null (res_cb +0.006 ± 0.003, below 0.01 pp; g8 shows +0.008 ± 0.003, one isolated "real").

D3 verdict: real at 2.24, 3.49 and 3.49_old in res_cb, leak_cut (and tail_hi3 at 3.49), same sign on both 3.49 GeV types, out of sample unchanged (no U2 gap problem; gap within 1 SE). Adoption test not passed as written: `tail_lo3` rises at 3.49 GeV (+0.069 ± 0.017, +0.126 ± 0.023; > 2 SE) and the peak scale shifts (Δmu +0.007 / +0.012 %); the other guards hold (eff_cut up, yield loss ≤ 0.4 %). The res_cb gain at 3.49 GeV (−0.026 to −0.033 pp) matches the expected bound (≤ 0.03 pp), and the leak_cut gain is the larger effect. Decision for the user and the P5 gate: whether the low-side tail rise is acceptable against the leak_cut gain (the tail_lo3 guard is the single failed adopt condition); no default change is made here.
