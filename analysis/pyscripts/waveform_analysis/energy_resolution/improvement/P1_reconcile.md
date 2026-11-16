# P1: reconcile the 2.24 GeV benchmark; Møller rejection and NL comparisons (2026-11-13)

Control: `output/energy_resolution/baseline_20261112/pre/<RT>`. Variants: `output/energy_resolution/variants/P1/{a_moller_reject,b_moller_reject_1cb,g_nonlinear}/<RT>`
(`delta.tsv` per run type, logs in `variants/P1/logs`). Same scans, same splits, one switch per variant; reports run with
gates `8 12`. Scan banks (sha1): 0.73 GeV 5d588cce, 2.24 GeV H2 15dfd867 (bank no longer on disk, D5), 3.49 GeV (both) 294184b0,
Ta 8faf5d6f. Full table: `P1_reconcile.tsv`. Values are legacy path (g12 / g8 in the TSV), in-sample `res_cb` in %, Δ = variant − control ± paired SE.
Both switches stay default OFF (U3).

## Row a: `MOLLER_REJECT=True` (all five run types)

| Run type | control | variant | Δ res_cb | Δ oos | Δ tail_hi3 (pp) | Δ leak_cut (pp) | Δ eff_cut (pp) | n_used ctrl → var | D3 |
|---|---|---|---|---|---|---|---|---|---|
| 0p73GeV_H2 | 3.389 | 3.556 | +0.167 ± 0.013 | +0.152 ± 0.010 | +0.011 ± 0.004 | +0.009 | +34.6 | 443k → 316k | real, worse |
| 2p24GeV_H2 | 1.876 | 1.885 | +0.009 ± 0.001 | +0.010 ± 0.001 | +0.037 ± 0.001 | +1.22 | +17.3 | 1224k → 1001k | real, slightly worse |
| 3p49GeV_H2 | 1.585 | 1.585 | 0 | 0 | +0.001 | +0.24 | +0.9 | 382k → 375k | null on res_cb |
| 3p49GeV_H2_old | 1.559 | 1.559 | 0 | 0 | +0.001 | +0.24 | +0.9 | 579k → 569k | null on res_cb |
| 2p24GeV_Ta_X17 (g8) | 2.035 | 2.017 | −0.018 ± 0.009 | −0.018 ± 0.008 | +0.092 ± 0.012 | +0.67 | +7.7 | 178k → 169k | res_cb about 2 SE better; tail_hi3 worse |

- It removes 18 % (2.24 GeV) and 29 % (0.73 GeV) of the clusters, so tail, `eff_cut` and `leak_cut` change population as well as shape (the
  low tail loses 18-31 pp: Møller events). Those columns are not like for like; only `res_cb` is.
- Rejection does not help the H2 core at any energy (+0.17 pp at 0.73, +0.009 at 2.24, no effect at 3.49 where there is no Møller pool).
  Adoption fails D3 at every H2 energy. Gap (in-sample vs oos) unchanged (|Δgap| < 0.02 pp, inside 2 SE).

## Row b (= `ref`): row a + `TWO_PEAK_WITH_REJECT=False`, `2p24GeV_H2`, 42 splits

`ref` in-sample `res_cb` = 1.8822 ± 0.0045 % (oos 1.8952; g12 1.8697). Against control 1.8759: Δ +0.006 ± 0.001. The single CB is 0.003 pp
better than the two-peak fit with rejection (1.8853), neither is better than no rejection. `ref` ≤ 1.81 % is **not reached**.

## Attribution of 1.876 (BM) vs 1.794 % (old note), 2.24 GeV H2 legacy in-sample

| Step | BM (this table) | Old chain (`opens.md:38`) |
|---|---|---|
| no gains, no position correction | 2.0045 | 2.015 (same selection; pulse chain equivalent) |
| gains (c, read from `poscorr_comparison.tsv`) | 1.9241 (−0.080) | −0.18 |
| position correction (c) | 1.8759 (−0.048) | −0.05 |
| Møller rejection (a / b) | +0.009 / +0.006 | −0.01 |
| fixed gain windows (d, read) | not in BM; old +0.004 | +0.004 |
| `ref` | 1.8822 | 1.794 |

Residual 1.8822 − 1.794 = 0.088 pp (about 14 combined SE; the 1.81 criterion misses by 0.072 pp). Per-step: the gain step is smaller here
by about 0.10 pp (−0.080 vs −0.18), the Møller step has the opposite sign (+0.006 vs −0.010, 0.016 pp), the position step matches. Measured
by P1: that the gap sits in the gain solve and the rejection step. Not measurable without the fallback (D6: rebuild the old chain for run
25320, about 3 CPU-h plus a best-50 bank): why the old gain solve reached −0.18 pp (its best-50 template bank, per-crystal solver and
selection are not in BM; the old chain's 1.794 note-only value cannot be rechecked, its outputs are gone). Row f (`R_MODEL` 0.3) is not run
here; it is P3's rescan and cannot touch the gain solve. The step (c) in-sample number also includes the out-of-sample gap, 0.011 pp, too small to matter.

## Row g: `NONLINEAR=True`

`2p24GeV_H2` and `2p24GeV_Ta_X17`: the 2.24 GeV calibration has nl = 0, the switch is an identity, no run.

| Run type | control | variant | Δ res_cb | Δ oos | Δ tail_hi3 (pp) | Δ leak_cut (pp) | Δ eff_cut (pp) | D3 |
|---|---|---|---|---|---|---|---|---|
| 0p73GeV_H2 | 3.389 | 3.355 | −0.033 ± 0.011 | −0.038 ± 0.016 | +0.001 ± 0.003 | +0.002 | −0.9 | res_cb real in sample (3 SE), oos 2.3 SE; g12 null (−0.014 ± 0.016); n_used 443k → 481k |
| 3p49GeV_H2 | 1.585 | 1.681 | +0.095 ± 0.003 | +0.095 ± 0.004 | −0.017 ± 0.006 | +1.01 | −0.9 | real, worse |
| 3p49GeV_H2_old | 1.559 | 1.662 | +0.103 ± 0.002 | +0.106 ± 0.003 | +0.019 ± 0.006 | +1.22 | −0.8 | real, worse (replicates) |

NL as applied by the clusterer's factor widens the 3.49 GeV core by about 0.1 pp and costs `eff_cut` 0.9 pp, with `mu` unchanged
(legacy |Δmu| < 0.004 pp); a gain at 0.73 GeV is of 2-3 SE and does not carry over to g12. The 0.73 GeV `res_g2` moves by +1.35 pp, a symptom of the
unstable Gaussian fit there (`opens.md:5`), not used. No adoption on D3 at any energy; kill for the H2 energies where it matters.

## Verdicts

- Møller rejection: no adoption (worse on H2 core at 0.73 and 2.24 GeV, null at 3.49 GeV; Ta marginal, with a tail_hi3 worsening).
- NL: no adoption (worse at 3.49 GeV, replicated; a small 0.73 GeV gain that fails g12 and replication).
- Both default OFF stands; the user decides on the table (U3).
- Open for the user: close P1 on the residual (0.088 pp, mostly the old gain solve) or authorise the D6 fallback (about 3 CPU-h).
