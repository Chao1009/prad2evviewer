# P16: below-threshold (KTeV-style) correction

2026-11-13. Control: phase-2 control v2 (`control_phase2_v2.md`), rerun as `output/energy_resolution/variants/P16/ctrl/` with the cluster dump
on (`--set THR_DUMP=True`). Its `regions_splits.tsv` pool rows equal control v2 at the control gate (max |diff| 0 on all four H2 run types).
Code sha1 4764b432 (`variants/P16/code.sha1`). Scripts: `variants/P16/run_p16.sh` (control and variant), `run_p16_cross.sh` (cross-fitted
table). Tools: `tools/p16_common.py`, `p16_slice1.py`, `p16_table.py`, `p16_d3.py`. Tables: `P16_threshold_corr.tsv` (D12),
`variants/P16/p16_slice1.tsv`, `p16_d3_kt.tsv`, `p16_d3_kt_cross.tsv`. Units: `res_cb` is the Crystal Ball σ/E of the ep pool (%); tails,
`leak_cut` and `eff_cut` are shares in %; Δ is variant minus control, paired per (run, split), mean ± SE (D3).

## Slice 1: is the threshold loss visible in the data?

Module threshold: T_m = TET × k_m (k_m = MeV per amplitude ADC from the P12 edge table; for 3.49 GeV old runs from `variants/P16/raw/`),
multiplied by the run's gain of the module (calibrated MeV). The brief's split (modules below 8 MeV against modules above 11 MeV) cannot be
applied per cluster. Only 9-13 of 1148 modules sit below 8 MeV and 52-68 above 11 MeV (raw thresholds). A cluster reads 8-24 modules. Each cluster therefore gets one weighted threshold T_w: the mean T_m over its 5 × 5 window, weighted by
P(1 − P), where P is the share of modules at the same distance d that are read out. The weight picks the modules that sit near threshold:
d ≈ 1.0-1.5 module sizes at 0.73 GeV, 1.5-2.5 at 2.24 GeV and 1.75-2.75 at 3.49 GeV. The classes are lo (T_w below its 25 % quantile) and
hi (above its 75 % quantile), with the quantiles taken in 25 mm rings of seed radius so that both classes have the same radii. The lo and
hi class means are 9.5 / 10.3, 9.2 / 9.7, 9.1 / 9.7 and 9.2 / 9.8 MeV (0.73 / 2.24 / 3.49 / 3.49 old).

| run type | Δ res_cb oos (hi − lo) | Δ mu_none | Δ mu, gains only | within-seed: absent modules per MeV | within-seed: read-out ring energy (MeV per MeV) | within-seed: E (%/MeV), gains only |
|---|---|---|---|---|---|---|
| 0.73 | +0.018 ± 0.046 | −0.33 ± 0.12 | −0.13 ± 0.03 | −0.11 ± 0.02 | +5.8 ± 0.5 | −1.62 ± 0.21 |
| 2.24 | +0.015 ± 0.012 | +0.18 ± 0.02 | +0.12 ± 0.01 | **+0.52 ± 0.03** | **−4.0 ± 0.9** | +0.98 ± 0.05 |
| 3.49 | +0.034 ± 0.027 | −0.09 ± 0.14 | +0.15 ± 0.02 | **+0.45 ± 0.05** | **−6.3 ± 1.9** | +1.01 ± 0.07 |
| 3.49 old | +0.000 ± 0.018 | −1.40 ± 0.35 | +0.08 ± 0.01 | **+0.39 ± 0.05** | **−4.7 ± 2.1** | +0.90 ± 0.06 |

The within-seed columns are slopes on T_w per split, with the mean of each seed subtracted and u, v, u², v², uv as nuisance terms. The
read-out ring is the modules at d ≥ 1.25. Edge and centre rows, tails and `eff_cut` are in `p16_slice1.tsv`.

**Verdict: the loss is real at the module level, but it does not show in the class split.**
- Module level. At 2.24 GeV, 3.49 GeV and the 3.49 GeV old runs, a higher T_w means more absent near-threshold modules (8-17 SE) and less
  read-out energy in the ring (2.2-4.6 SE). The ring loses 4-6 MeV per MeV of threshold, 40-50 % of the P12 extrapolation (about 10 and
  13 MeV per MeV). So the mechanism is real at ≥ 2 SE at 2 energies and replicates in the old runs.
- 0.73 GeV. T_w does not predict absence (the sign is reversed), so the 0.73 GeV module-level test is inconclusive.
- Class split. The cluster-level response does not follow the threshold. `res_cb` (hi − lo) is null at every energy. The total energy moves
  with T_w in the wrong direction (+0.9 to +1.0 %/MeV within seed at 2.24 and 3.49 GeV). A module-property confound therefore dominates the class split. It survives the gains and an
  energy-weighted near-neighbour regressor (`*_thr2` and `*_near2` rows). Reading the classes as a threshold effect would give the wrong
  sign.
- Power. The class contrast is 0.5-0.8 MeV. With the P12 term, the expected Δ `res_cb` is 0.01-0.03 pp, below the split's 1 SE. A class
  split is the wrong instrument for the width term. Slice 2 measures the same term directly, from which modules are absent.
- By the brief's literal criterion (class split: `res_cb` or mu at ≥ 2 SE with a consistent sign at ≥ 2 energies) the loss is not
  confirmed. Slice 2 was run on the module-level criterion instead; this substitution is flagged for the supervisor.

## Slice 2: the correction

**Method.** For each cluster, every usable PbWO4 module of the 5 × 5 window around the seed that has no stored hit in the event (any time,
any path) gets a virtual member with energy c(d, T class) × T_m. Usable means the threshold is known and the module is not dead. d is the
distance of the module centre from the cluster position, in module sizes, in bins of 0.125. The T classes split at 9.0 and 9.8 MeV.

The table (`p16_table.py`, one per run type, 3 × 24 bins) comes from the control's peak clusters (|E/E′ − 1| < 2 × the core width). It
works in units of x = e/T_m, which do not depend on the gains. In each bin, the true deposits below x* = 1.2 (full efficiency above x*)
are modelled as a truncated exponential. Its mass is the absent share plus the read-out share below x*. Its edge density is the measured
density at x*. The expected energy of an absent module is c = (mass × truncated mean − read-out energy below x*) / absent share. Bins
with an absent share below 0.005 or above 0.95 get c = 0.

Resulting c: 0.6-0.7 T (6-7 MeV) where few modules are absent, falling to 0.25-0.35 T where most are absent. The added energy per peak
cluster is 24 / 33 / 31 MeV (rms 5-7 MeV), and mu_none rises by +3.1 / +1.5 / +0.9 %. Before the run, the variance a presence-based
correction can remove was estimated as Σ P(1 − P) c² = 17 / 42 / 49 MeV². This bound predicts −0.04 / −0.026 / −0.011 pp. The observed oos gains exceed it at 0.73 GeV (5×) and 3.49 GeV
(1.7×), so part of the gain is not the presence term (see edge and centre below).

The virtual members enter before every correction. Gains, the position correction and the P9 shape factor are refitted in the report:
`kt/<rt>/poscorr_extra.tsv` coefficients differ from the control's at 2.24 and 3.49 GeV. `shape_feats` counts real members only, so its
count feature keeps its meaning. At 0.73 GeV the control has no position correction and no shape factor, so the correction acts alone
after the gains.

**Overfitting (U2).** The table has 72 numbers per run type, built from 1.8-7.3 × 10⁵ peak clusters (10⁶-10⁷ module entries).
- Even- and odd-half tables agree to a weighted rms of 1.4-3.5 % of c.
- In the cross-fitted run, each split is scored with the table built on the other half of the splits (gains and position correction are
  out of sample as usual). It reproduces the full-table Δ to ≤ 0.004 pp.

**D3, cross-fitted (`kt_cross`), pool, control gate:**

| run type | Δ res_cb in | Δ res_cb oos | Δ tail_hi_f | Δ tail_lo_f | Δ tail_hi3 | Δ leak_cut | Δ eff_cut | peak shift (mu) | yield loss |
|---|---|---|---|---|---|---|---|---|---|
| 0.73 (control 3.365 / 3.397) | −0.196 ± 0.006 | **−0.215 ± 0.012** | −0.001 ± 0.001 (null) | **−1.49 ± 0.07** | +0.008 ± 0.003 | −0.000 ± 0.001 (null) | +0.03 ± 0.03 (null) | −0.002 ± 0.012 (null) | −2.46 % (gain) |
| 2.24 (1.790 / 1.802) | −0.024 ± 0.002 | **−0.026 ± 0.002** | −0.001 ± 0.001 (null) | **+0.053 ± 0.006** | +0.002 ± 0.002 (null) | **−0.18 ± 0.01** | +0.05 ± 0.01 | +0.016 ± 0.001 | −0.02 % |
| 3.49 (1.501 / 1.526) | −0.019 ± 0.002 | **−0.019 ± 0.003** | −0.012 ± 0.004 | −0.004 ± 0.010 (null) | +0.001 ± 0.002 (null) | **−0.26 ± 0.06** | +0.06 ± 0.02 | +0.014 ± 0.002 | −0.05 % |
| 3.49 old (1.479 / 1.500) | −0.018 ± 0.002 | **−0.018 ± 0.002** | −0.007 ± 0.004 (null) | −0.04 ± 0.02 | +0.002 ± 0.003 (null) | **−0.17 ± 0.05** | +0.08 ± 0.02 | +0.012 ± 0.002 | +0.05 % |

Tails, `leak_cut` and `eff_cut` are out of sample.

**Edge and centre (oos, cross-fitted):**

| run type | Δ res_cb edge | Δ res_cb centre | Δ eff_cut edge | Δ eff_cut centre | Δ tail_lo_f edge | Δ tail_lo_f centre | Δ peak edge | Δ peak centre |
|---|---|---|---|---|---|---|---|---|
| 0.73 | −0.183 ± 0.015 | −0.216 ± 0.015 | −0.05 ± 0.04 | +0.07 ± 0.07 | −1.36 ± 0.08 | −1.56 ± 0.07 | +0.143 ± 0.013 | −0.162 ± 0.018 |
| 2.24 | −0.021 ± 0.004 | −0.027 ± 0.003 | +0.07 ± 0.01 | +0.03 ± 0.01 | +0.057 ± 0.009 | +0.051 ± 0.008 | +0.015 ± 0.003 | +0.017 ± 0.002 |
| 3.49 | −0.019 ± 0.004 | −0.018 ± 0.002 | +0.05 ± 0.03 | +0.06 ± 0.02 | −0.006 ± 0.017 | −0.001 ± 0.014 | +0.007 ± 0.003 | +0.016 ± 0.002 |
| 3.49 old | −0.020 ± 0.003 | −0.013 ± 0.005 | +0.15 ± 0.02 | +0.08 ± 0.03 | −0.05 ± 0.02 | −0.04 ± 0.02 | +0.004 ± 0.003 | +0.011 ± 0.003 |

The edge and centre rows are refitted from the dumps (the same `fit_path`; the dump fit of the pool reproduces the official pool `res_cb`
to 6e-6). At 0.73 GeV the correction moves the edge peak up relative to the centre peak. The edge − centre peak offset falls from
−0.415 % to −0.100 %, and edge and centre each narrow by about 0.18 pp (in sample). The correction therefore partly acts as the physical position
correction that 0.73 GeV lacks, beyond the presence term.

The same result holds for the other estimators (splits.tsv, control gate, in sample, full table): `res_g2` Δ −0.201 / −0.027 / −0.017 / −0.014.

## Recommendation (U7, D3)

- **Adopt at 0.73 GeV, 3.49 GeV and 3.49 GeV old.** Both figures of merit gain or hold, and every guard passes:
  - `res_cb` oos improves by 7-18 SE at all three H2 energies, the 3.49 GeV gain replicates in the old runs, and 0.73 GeV gains most
    (−0.215 pp, `sep_pi` +6 %).
  - `leak_cut` −0.17 to −0.26 at 2.24 / 3.49 GeV.
  - `tail_hi_f` null or lower. At 0.73 GeV, `tail_hi3` +0.008 ± 0.003 is the 3σ line moving in (σ-narrowing; the amended guard uses
    `tail_hi_f`, which is null).
  - `eff_cut` up or null, `tail_lo_f` null or lower at these three, and the yield changes by at most 0.05 % or rises (0.73 GeV: +2.5 %).
- **2.24 GeV: one guard fails formally.** `tail_lo_f` rises by +0.053 ± 0.006 (9 SE). That is 0.2 % of its 26.9 % control value, while
  `eff_cut` rises by +0.05 and `res_cb` and `leak_cut` gain at 13 and 17 SE. Adopting at 2.24 GeV needs a ruling on this guard.
- Mechanism of the 2.24 GeV `tail_lo_f` rise: not measured. Plausibly the table is built at the elastic energy and gives the low-energy
  tail clusters (Møller and radiative, inside the fit range at 2.24 GeV) less relative energy than the peak.
- The switch stays default off, as briefed. Turning it on means setting `BELOW_THR` to a directory holding the tables of the adopted run
  types: a run type without a table is unchanged.
- **Ta (2.24 GeV X17).** Not measured. Ta borrows the 2.24 GeV position correction from the report root. If the 2.24 GeV report is built
  with the correction, Ta needs its own table first, or Ta's borrowed map no longer matches its energies.

## Benchmark edits (additive, default off; flagged P16 in the file)

- Two new switches, both in `OVERRIDABLE`:
  - `THR_DUMP` (False) writes `thr_clusters.npz`: the control-gate pool clusters with e (in sample), e_oos, er, position, raw members, the
    used flags, and the 5 × 5 window modules with a stored-hit flag.
  - `BELOW_THR` ("") names the directory of tables `<run type>.npz`.
- `thr_window()` and `below_thr()` sit before `cluster_events`. `split_job` calls `below_thr` right after `cluster_events`, when
  `BELOW_THR` is set, and builds the dump.
- `report` writes the dump. `shape_feats` subtracts `nv`, the number of virtual members. Without `nv` (the default) the result is
  `c / 9.0`, as before.
- Default regression: the control rerun with `THR_DUMP` on reproduces control v2 exactly; `BELOW_THR` "" leaves the code path untouched.

Open: the 2.24 GeV `tail_lo_f` guard ruling; a Ta table; whether the 0.73 GeV gain partly replaces a position correction (a 0.73 GeV
position-correction arm on top of the correction is untested); the slice-1 module-property confound (rising energy with T_w) is not
identified.
