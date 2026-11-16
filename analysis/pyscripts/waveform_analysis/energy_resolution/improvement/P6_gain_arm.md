# P6: gain solve as a comparison arm (2026-11-13)

Control: `output/energy_resolution/control_phase2_20261113/<RT>` (TIME_OFFSETS on, gates 20/16/16/16/8, log weight, 24937 dropped). Variants:
`output/energy_resolution/variants/P6/<variant>/<RT>/` (`delta.tsv` per run type; logs and `run1.sh`, `batch*.sh` in `variants/P6/`). Gates 12 and the control gate; same
scans, same splits, one switch per variant; scan bank sha1 0.73 GeV 5d588cce, 2.24 GeV 15dfd867, 3.49 GeV 294184b0, Ta 8faf5d6f. Benchmark sha1 at the last run 364f967e.
Full table: `P6_gain_arm.tsv` (built by `tools/p6_table.py`; res_cb in %, Δ = variant − control in pp, paired SE over splits; gains-only variant `g1` is against `g0`, not the control).
Refinement stays default OFF (U2); the control's member ridge solve is unchanged.

## What the old chain did that the benchmark does not

The old chain's 1.794 % applied two gain layers: the per-crystal member solve (`percrystal_calib.py`, same ridge solve as `solve_gains`, fixed windows) and then the
**seed-indexed refinement** `fit_elastic_window_gains` (OLD:823): one gain per seed module = E' / median of that module's in-window clusters, refitted by the benchmark
itself **per job, per path, on the events it then evaluates** (`PbWO4_energy_resolution_20261025.md:20`, `:863`). The benchmark has the member solve only. The member solve is
the same (iteration count, fixed windows, prior, seed fraction and low-side window are all null, below); the difference is the second layer.

## Benchmark edits (additive, `--set`-switchable, defaults unchanged, `energy_resolution_benchmark.py`)

- New constants `GAIN_MODE` ("member"; "seed" replaces the member solve by the seed-indexed one), `GAIN_SEED_STEP` ("" off; "run": seed-indexed step per run stacked after member gains
  and position correction, in sample = all splits of the run, out of sample = the other half; "split": refit on the pool at hand per path and split, the old chain; "splithalf": that, cross-fitted over even / odd clusters),
  `GAIN_FIXED_WINDOWS`, `GAIN_MIN_SEED_FRAC`, `GAIN_WIN_LO_NSIGMA`; `OVERRIDABLE` gains those plus `GAIN_SEED_FRAC`, `GAIN_FREE_FRAC`, `GAIN_MAX_ITER`, `GAIN_N_REFIT`, `PC_SETTING`.
- New `solve_gains_seed`; branches in `solve_gains` (training sample switches), `gains_job`, `apply_arm`, `corrections` (run-level step: seed gains ride behind the member gains in the gains array, so the data flow is unchanged).
- Default check: `g_iter12`/`f_iter1` and the null rows agree with the control to ≤ 0.001 pp, and the control rows are untouched (no default path reads a new constant).

## Gain-step attribution at 2.24 GeV H2, legacy path, in sample (res_cb %)

| step | value | note |
|---|---|---|
| no gains, no position correction | 2.0044 | `res_cb_none` of the control report |
| member gains only (`g0`) | 1.9241 (−0.080) | reproduces P1's number |
| seed-indexed gains only instead of member (`g1`) | 1.9482 (−0.056) | worse than member by +0.023 ± 0.002 out of sample |
| + position correction (control) | 1.8648 | |
| + seed-indexed step per run (`s1`) | 1.8623 | in-sample −0.003 ± 0.001, out of sample +0.001 ± 0.001 |
| + seed-indexed step per split, own events (`s2`, old chain) | **1.8358** | −0.029 ± 0.002 in sample |
| `s2` cross-fitted within the split (`s3`) | 1.8884 | +0.024 ± 0.002 in sample, +0.017 ± 0.003 "oos" |
| `s2` with Møller rejection, single CB (`s4`, old chain's `ref` settings) | 1.8420 | rejection adds +0.006 as in P1 |
| old chain, note-only | 1.794 | |

The 0.088 pp residual of P1 (`ref` 1.882 vs 1.794): the old chain's in-sample seed refit on its own events accounts for 0.029 pp (s2). The rest, about 0.06 pp (0.042 left after s2 on the control
basis, plus the Møller step, which here has the opposite sign, +0.006 vs −0.010), is not attributable without the D6 rebuild. The gain step itself (−0.080 vs −0.18) is not a solver
difference: the member solve is insensitive to everything the old chain set (table below).

## Arm table, 2.24 GeV H2 (control gate 16 ns, 42 splits); other run types in the tsv

| variant (`--set`) | Δ in | Δ oos | Δ tail_hi3 oos | Δ leak_cut oos | Δ eff_cut | Δ tail_lo_f | Δ gap | yield loss % | verdict (D3 / U2) |
|---|---|---|---|---|---|---|---|---|---|
| a_seed (`GAIN_MODE="seed"`) | +0.003 ± 0.002 | +0.004 ± 0.002 | +0.001 ± 0.001 | +0.04 ± 0.01 | -0.02 ± 0.01 | +0.01 ± 0.00 | +0.001 ± 0.002 | -0.03 | null |
| b_fixedwin (`GAIN_FIXED_WINDOWS=True`) | -0.000 ± 0.001 | +0.001 ± 0.001 | -0.001 ± 0.001 | +0.01 ± 0.01 | +0.00 ± 0.00 | -0.00 ± 0.00 | +0.002 ± 0.001 | -0.04 | null |
| c_prior_x4w (`GAIN_SEED_FRAC=0.32`) | +0.001 ± 0.001 | +0.002 ± 0.001 | -0.001 ± 0.001 | -0.00 ± 0.01 | +0.00 ± 0.00 | +0.00 ± 0.00 | +0.002 ± 0.001 | -0.00 | null |
| d_prior_x16w (`GAIN_SEED_FRAC=0.16`) | +0.001 ± 0.001 | +0.001 ± 0.001 | -0.001 ± 0.001 | -0.01 ± 0.01 | +0.01 ± 0.01 | -0.00 ± 0.00 | +0.001 ± 0.001 | -0.01 | null |
| e_prior_x4s (`GAIN_SEED_FRAC=1.28`) | +0.000 ± 0.001 | +0.002 ± 0.002 | +0.001 ± 0.001 | +0.02 ± 0.01 | -0.00 ± 0.00 | +0.00 ± 0.00 | +0.002 ± 0.002 | -0.02 | null |
| f_iter1 (`GAIN_MAX_ITER=1`) | -0.000 ± 0.001 | +0.001 ± 0.001 | -0.001 ± 0.001 | +0.01 ± 0.01 | +0.00 ± 0.00 | -0.00 ± 0.00 | +0.002 ± 0.001 | -0.04 | null |
| g0_ctrl_gonly (`PC_SETTING={"0p73GeV_H2":("none",0),"2p24GeV_H2":("none",0),"3p49GeV_H2":("none",0),"3p49GeV_H2_old":("none",0),"2p24GeV_Ta_X17":("none",0)}`) | +0.054 ± 0.003 | +0.054 ± 0.003 | +0.002 ± 0.002 | +0.47 ± 0.02 | -0.19 ± 0.02 | +0.03 ± 0.01 | -0.000 ± 0.002 | -0.24 | worse; guard: eff_cut down |
| g1_seed_gonly (`PC_SETTING={"0p73GeV_H2":("none",0),"2p24GeV_H2":("none",0),"3p49GeV_H2":("none",0),"3p49GeV_H2_old":("none",0),"2p24GeV_Ta_X17":("none",0)} GAIN_MODE="seed"`) | +0.025 ± 0.002 | +0.022 ± 0.002 | +0.008 ± 0.002 | +0.25 ± 0.01 | -0.03 ± 0.01 | -0.00 ± 0.00 | -0.003 ± 0.002 | -0.03 | worse |
| g_iter12 (`GAIN_MAX_ITER=12 GAIN_N_REFIT=6`) | +0.000 ± 0.001 | +0.001 ± 0.001 | -0.001 ± 0.000 | +0.00 ± 0.00 | +0.00 ± 0.00 | +0.00 ± 0.00 | +0.001 ± 0.001 | +0.00 | null |
| h_seedfrac50 (`GAIN_MIN_SEED_FRAC=0.5`) | +0.001 ± 0.001 | +0.001 ± 0.001 | -0.001 ± 0.001 | -0.00 ± 0.01 | -0.00 ± 0.01 | -0.01 ± 0.00 | +0.001 ± 0.001 | +0.03 | null |
| i_winlo1 (`GAIN_WIN_LO_NSIGMA=1.0`) | -0.006 ± 0.003 | -0.003 ± 0.003 | +0.003 ± 0.001 | -0.07 ± 0.01 | +0.07 ± 0.01 | -0.07 ± 0.00 | +0.002 ± 0.002 | +0.72 | null |
| j_winlo15 (`GAIN_WIN_LO_NSIGMA=1.5`) | +0.000 ± 0.003 | -0.000 ± 0.002 | -0.000 ± 0.001 | -0.03 ± 0.01 | +0.04 ± 0.01 | -0.03 ± 0.00 | -0.001 ± 0.002 | +0.23 | null |
| s1_seedstep_run (`GAIN_SEED_STEP="run"`) | -0.004 ± 0.001 | +0.003 ± 0.001 | -0.000 ± 0.001 | +0.02 ± 0.01 | -0.00 ± 0.01 | +0.01 ± 0.00 | +0.006 ± 0.002 | +0.02 | null |
| s2_seedstep_split (`GAIN_SEED_STEP="split"`) | -0.028 ± 0.002 | -0.033 ± 0.002 | +0.016 ± 0.002 | -0.06 ± 0.01 | -0.00 ± 0.01 | -0.01 ± 0.00 | -0.004 ± 0.002 | +0.38 | better |
| s3_seedstep_splithalf (`GAIN_SEED_STEP="splithalf"`) | +0.022 ± 0.003 | +0.017 ± 0.003 | +0.002 ± 0.002 | +0.12 ± 0.01 | -0.01 ± 0.01 | -0.03 ± 0.00 | -0.005 ± 0.002 | +0.02 | worse |
| s4_split_ref (`GAIN_SEED_STEP="split" MOLLER_REJECT=True TWO_PEAK_WITH_REJECT=False`) | -0.023 ± 0.003 | -0.028 ± 0.002 | +0.050 ± 0.002 | +1.09 ± 0.02 | +17.40 ± 0.04 | -18.39 ± 0.04 | -0.005 ± 0.002 | +0.40 | better |

`g0`, `g1` are gains-only (`PC_SETTING` none) rows, Δ against `g0` for `g1` and against the control for `g0`. In `s2` and `s4` the "oos" column is not out of sample (the step refits on the evaluated pool);
only `s1` (run level, other-half fit) and `s3` (cross fit) have a held-out column. Fitted factors: member solve 533 per run, 284 free (prior fraction < 0.5), 4109 / 498 events per free module (median / min) at 2.24 GeV;
seed-indexed step 264 factors per run on a median 218 (min 47) in-window clusters per module per run, so about 23 per module in a per-split fit (the old chain's setting), below `MIN_PER_MODULE` 30.

## Other run types (control gate; seed-indexed step per run `s1` and per split `s2`)

| run type | s1 Δ in | s1 Δ oos | s2 Δ in | s2 Δ tail_hi3 | s2 Δ tail_lo_f | s2 yield loss % |
|---|---|---|---|---|---|---|
| 0.73 GeV (g20) | +0.000 ± 0.006 | −0.012 ± 0.004 | −0.040 ± 0.010 | +0.019 ± 0.005 | −0.24 ± 0.06 | +1.30 |
| 3.49 GeV (g16) | −0.006 ± 0.002 | −0.002 ± 0.003 | −0.027 ± 0.004 | +0.008 ± 0.010 | −0.06 ± 0.01 | +0.56 |
| Ta X17 (g8) | −0.014 ± 0.008 | −0.019 ± 0.010 | −0.084 ± 0.022 | +0.065 ± 0.017 | +0.35 ± 0.08 | +0.24 |

3.49 GeV old runs not run (3.49 GeV already null; same code path). Ta `s2` in sample reproduces the old chain's −0.07 to −0.09 pp; its tail_hi3 and tail_lo_f both worsen.

## Member-solve settings (2.24 GeV, control gate): all null

| variant | Δ oos | verdict |
|---|---|---|
| prior strength λ/4 (`GAIN_SEED_FRAC=0.32`), λ/16 (0.16), 4λ (1.28) | +0.002 ± 0.001, +0.001 ± 0.001, +0.002 ± 0.002 | null |
| one pass at the g = 1 windows (`GAIN_MAX_ITER=1`, = fixed windows) / 12 passes | +0.001 ± 0.001 / +0.001 ± 0.001 | null |
| training clusters with seed fraction ≥ 0.5 | +0.001 ± 0.001 | null |
| lower training-window edge 1.0 σ / 1.5 σ (`GAIN_WIN_LO_NSIGMA`) | −0.003 ± 0.003 / −0.000 ± 0.002 | null (1.0 σ: yield loss +0.72 %, `eff_cut` +0.07, `tail_lo_f` −0.07: the cut moves the population) |

## Verdicts (D3, U2)

- Seed-indexed refinement, in sample, per split (the old chain): −0.03 pp at 0.73, 2.24, 3.49 GeV and −0.08 pp on Ta; it **worsens `tail_hi3`** (+0.017 ± 0.004 at 0.73, +0.016 ± 0.002 at 2.24, +0.008 ± 0.010 at 3.49 GeV, Ta +0.065 ± 0.017) and cross-fitted inside the split it is **+0.017 pp worse**. A training artefact of fitting about 23 clusters per module on the events then scored.
- Seed-indexed refinement at run level, held-out half (the judge): Δ oos +0.003 ± 0.001 (2.24), −0.002 ± 0.003 (3.49), −0.019 ± 0.010 (Ta, null), −0.012 ± 0.004 (0.73 GeV, 3 SE but 0.012 pp, in-sample 0.000). Gap +0.006 ± 0.002 pp at 2.24 GeV (within max(2 SE, 0.02)). D7 kill rule met at 2.24 and 3.49 GeV: no out-of-sample Δ ≤ −2 SE at the H2 energies where it counts.
- Member solve: prior strength, passes, windows, seed-fraction sample: all null; no setting beats the control.
- D6 rebuild not run: it could only restate a note-only in-sample number; the decision (no out-of-sample gain) does not depend on it.

## Recommendation

Keep gain refinement a comparison arm, default OFF. Do not promote it: its apparent −0.03 to −0.18 pp is in-sample (per-split self-fit); held out it is 0 at 2.24 and 3.49 GeV and ≤ 0.02 pp elsewhere. Only the 0.73 GeV run-level −0.012 pp deserves a replication (second 0.73 GeV run set) before anyone cites it, and it is far below the 0.04 pp decision scale. U2 stays the user's call; no default changed.
