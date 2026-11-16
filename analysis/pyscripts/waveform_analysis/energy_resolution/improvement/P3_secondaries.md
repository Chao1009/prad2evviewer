# P3: secondary pulses in cluster members (2026-11-13)

Paths relative to `analysis/pyscripts/`. Tool: `waveform_analysis/energy_resolution/improvement/tools/p3_cells.py`.
Control for (A): `output/energy_resolution/baseline_20261112/pre`; for (B): `.../baseline_20261112/c0/2p24GeV_H2` (same bank sha1 7f125b75 on all arms).

## A. (dt, A2/A1) cells (`output/energy_resolution/variants/P3/<RT>/`, `SECONDARY_DUMP=True`, gates 8 12; full table `output/energy_resolution/P3/cells.tsv`)
Default-report outputs are identical to `pre` (compare_reports: all d = 0). Quantities are scan-level E/E' (no gains, no position
correction) of legacy pool clusters in 0.6-1.4; reference = clusters with no secondary on a member module; mu, sigma of the reference peak
(Gaussian window); tails at mu -/+ 3 sigma_ref. Cells are not independent of cluster geometry (more members, inner modules), so a width
excess is an association, not a cause.

| RT | flagged share (in range) | ref sigma | with-secondary sigma | tail_lo3 ref / sec | tail_hi3 ref / sec | A2/A1 >= 0.25: share, tail_hi3 |
|---|---|---|---|---|---|---|
| 2p24GeV_H2 | 25.5 % | 1.918 % | 2.139 % | 28.6 / 28.8 % | 0.105 / 0.255 % | 1.4 %, 2.91 % |
| 3p49GeV_H2 | 34.8 % | 1.620 % | 1.788 % | 10.4 / 14.4 % | 0.085 / 0.148 % | 1.9 %, 1.08 % |
| 2p24GeV_Ta_X17 | 43.4 % | 2.086 % | (see cells.tsv) | 52.4 | 0.39 | pre-pulse (dt < 0) cells +10-18 pp hi |

- Width: almost every cell with n >= 500 has a peak-sigma excess of +0.1 to +0.4 pp (>= 2 SE); the largest are dt 16-30 ns with A2/A1 0.04-0.12.
- Median shift of each cell is about -(median e_extra / beam): -0.5 to -1.1 % for e_extra about 1 % of the beam. The legacy path keeps one
  pulse per module, so the extra pulse's energy is dropped from the cluster; this is the signature of real small pulses, not of a split of one pulse.
- Low side: dt 16-30 ns, A2/A1 0.04-0.12 carry +4 to +17 pp tail_lo3 (3.49 GeV: +10 to +17); dt < 16 and A2/A1 >= 0.25 carry a deficit (-10 to -14 pp, selection effect of large-amplitude members).
- High side: excess only where A2/A1 >= 0.25 or the pulse precedes the main (dt < 0): 2.24 GeV +8.6 pp (dt < 0), +0.5-0.7 pp (dt 0-16);
  3.49 GeV +5.6 pp (dt < 0). These cells hold 331 of 1557 (2.24) and 77 of 379 (3.49) clusters above mu + 3 sigma_ref. Clusters with any secondary hold 45 % / 48 % of that tail.
- Ta: dt < 0 cells are much larger (share 10-15 % of clusters) and carry +10 to +18 pp hi tail; not interpreted (inelastic physics, U-guard only).

## B. Paired rescans, run 25320 (8 splits), `output/energy_resolution/scan_variants/{minsep12,rmodel0p3}/{2p24GeV_H2,report}`; delta.tsv per arm

| Variant | flag share (legacy pu) | extra hits | d res_cb | d res_cb_oos | d tail_hi3 | d leak_cut | d tail_lo3 | d eff_cut | yield loss |
|---|---|---|---|---|---|---|---|---|---|
| C0 | 23.94 % | 400694 | - | - | - | - | - | - | - |
| `MIN_SEP_NS=12` | 23.51 % | 394990 | -0.0004 +- 0.0018 | -0.0022 +- 0.0023 | +0.0010 +- 0.0014 | +0.004 +- 0.004 | -0.008 +- 0.009 | +0.017 +- 0.006 | 0.006 % |
| `R_MODEL=0.3` | 0.34 % | 10586 | -0.036 +- 0.008 (real) | -0.037 +- 0.009 (real) | +0.014 +- 0.007 | -0.160 +- 0.034 (real) | -0.32 +- 0.04 | +0.63 +- 0.04 | -0.4 % (gain) |

Gated paths (d res_cb, same C0): g8 -0.030 +- 0.008, g12 -0.020 +- 0.006; leak_cut g8 -0.13 +- 0.04, g12 -0.11 +- 0.03; U2 gap change null (-0.001 to -0.003 +- 0.003-0.006).
`MIN_SEP_NS=12` is null in every figure of merit (the open to-do of `notes/PbWO4_energy_resolution_20261108.md:16` is answered: 8 ns costs nothing and gains nothing at 25320).
`R_MODEL=0.3` removes 97 % of the extra pulses and narrows the peak by 0.036 pp (4.6 SE; 0.02 pp at gate 12); only run 25320 and 8 splits, so by D5 not yet a verdict.

## Verdict (D3 / D7 gate)
Not "real small pulses, close": (A) cells with a width and tail excess exist (about 2 SE or more, many cells) and (B) one rescan variant moves a figure of merit by more than 2 SE.
P8 gate condition ("a cell carries a width or tail excess >= 2 SE and a rescan variant confirms it") is met on run 25320.
Caveats: (i) the narrowing could be removal of real small pulses' effect in the Crystal Ball fit rather than of spurious pulses (the median shift of the cells says the pulses carry about 1 % energy that legacy drops); (ii) `R_MODEL=0.3` also changes the main-pulse fit; (iii) D5 requires extension to the other 2.24 GeV runs and the 3.49 GeV runs (about 9.5 CPU-h for both arms) and, for a passing variant, the Ta guard (run 26138, 4 splits, gated); neither was run (outside this block).
