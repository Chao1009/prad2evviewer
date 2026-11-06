# Commands
Earlier command log (through 2026-11-03) preserved in Commands_log_pre20261105.md.

Command reference for the waveform chains: template extraction, matched filter, energy
resolution, timing/pile-up, non-linearity. Every command is taken from the scripts or the
chain drivers; anything reconstructed from logs is marked `(reconstructed)`. Paths are
relative to `analysis/pyscripts/` unless absolute. Run everything from that directory.

Frozen: `benchmark_hycal_timing_test_20261021.py` (shasum `3293ab96870cef4580e886eb2dda5ec37c8821d1`,
enforced by `output/r5ri_20261104/tools/run1_ri.sh`) and the chain tools under `output/*_20261104/tools/`.
Live (cleaned) versions: `spline_template_extract.py`, `mf_pulse.py`.

## 0. Environment

```bash
cd /Users/L00338853/work/PRad/prad2evviewer/analysis/pyscripts && source setup_env.sh && export MPLBACKEND=Agg
# setup_env.sh exports PRAD2_DATABASE_DIR=$HOME/work/PRad/prad2evviewer/database, DAQ_CONFIG, HC_MAP and
# prepends $HOME/work/PRad/prad2evviewer/build/python to PYTHONPATH; it must be sourced, not executed.
```

Data: one EVIO split per job, `/Users/L00338853/work/PRad/data/evio/prad_<run6>.evio.<split5>`.
The job lists hard-code `/Users/L00338853/Software/opt/miniforge3/bin/python3`; the drivers call `python3`.

Calibration comes from `$PRAD2_DATABASE_DIR/runinfo/general.json` through `calib_resolve.py`:
`resolve_calib(run)` returns the absolute path of the `calibration.file` in force for `run`, or exits
(`[ERROR] runinfo has no calibration file for run N`, `[ERROR] calibration file for run N not found: PATH`);
`check_calib_json(path)` exits when fewer than 100 modules have factor > 0, `check_pipeline(p)` when fewer
than half of the PbWO4 modules are calibrated. A script resolves it as (`benchmark_hycal_timing_scan.py:405-407`):

```python
calib = resolve_calib(run)
print(f"[calib] run {run}: {calib} ({check_calib_json(calib)} modules)", flush=True)
cmd.extend(["--hc-calib-file", calib])
```

```bash
# which file a run resolves to (derived one-liner; calib_resolve.py has no CLI)
python3 -c "from calib_resolve import resolve_calib, check_calib_json; p = resolve_calib(24386); print(p, check_calib_json(p))"
```

`--hc-calib-file` (from `_common.py`, default `""` = via runinfo) overrides it; the ri job lists pass the
resolved path explicitly (`ri_jobs.py` calls `resolve_calib(run)` when it writes the lines).

## 1. Template extraction chain (`spline_template_extract.py`)

Reference run: `output/templates_spline_best50_p1r3_20261104/run_chain.sh` with jobs `jobs/j000.sh..j237.sh`
(identical to `select_cmds.txt`). Arms `A` = `baseline_subbank`, `p1shape_r3`; run types
`T="0p73GeV_H2 2p24GeV_H2 2p24GeV_C 2p24GeV_Ta_X17 3p49GeV_H2"` (16/42/8/25/28 splits); `N` = chain dir.
Each arm holds `<arm>/<type>/{pulses,logs,check}/`, `templates_<type>.{npz,pdf,root}`, `summary.tsv`,
`build.log`, `check.log`, and `<arm>/ampcorr.json`. Stage order in `run_chain.sh`: all selects
(`ls $N/jobs/j*.sh | xargs -P 14 -n 1 bash`), 10 builds in parallel, ampcorr per arm, 10 checks in
parallel, then `CHAIN_DONE` in `chain.out`.

```bash
# select: one EVIO split -> candidate pulses npz (jobs/j000.sh; baseline arm). Ta_X17 uses --trigger physics-any.
/usr/bin/time -p python3 spline_template_extract.py select /Users/L00338853/work/PRad/data/evio/prad_024590.evio.00010 --trigger physics --mf-bank output/templates_spline_best50_20261026/0p73GeV_H2/templates_0p73GeV_H2.npz --per-module-best 50 --emax-frac 1.2 --mf-strict --max-modules-per-event 20 --exclude-trigger-bits 0x1000000 --out output/templates_spline_best50_p1r3_20261104/baseline_subbank/0p73GeV_H2/pulses/024590_00010.npz > output/templates_spline_best50_p1r3_20261104/baseline_subbank/0p73GeV_H2/logs/select_024590_00010.log 2>&1
# p1shape_r3 arm: same line with --out .../p1shape_r3/... plus the shape cut
#   --p1-shape-cut --p1-min-rising 3      (default --p1-min-rising is 4; output/templates_spline_best50_p1shape_20261104 uses the default)
# No --calib-file: select resolves the calibration from runinfo (resolve_calib) and runs check_calib_json.

# build: candidate npz files -> templates_<type>.npz/.root/.pdf + summary.tsv (run_chain.sh:8,12), per arm A and type R
BX="--pool-order energy --n-max 50 --fit-points direct --shape-sigma 0.012 --lambda-select heldout --lambda-chi2-sigma full --lambda-min 1e2 --lambda-fallback 1e3 --converge-tol 0.01 --converge-df 5e-3 --max-iterations 10 --osc-escalate"
python3 spline_template_extract.py build --pulses "$PWD/$N/$A/$R/pulses/*.npz" --run-type $R $BX --out $PWD/$N/$A/esc_stage/$R/templates_$R > $N/$A/esc_stage/$R/build.log 2>&1
# --osc-escalate: modules whose template has spurious wiggles get lambda x10 and a rebuilt pool until none; still wiggly at 1e4 -> run-type average.

# ampcorr: compare the run types' f(A); write one global f(A) and rewrite the banks in place (run_chain.sh:18)
python3 spline_template_extract.py ampcorr --out $PWD/$N/$A/esc_stage/ampcorr.json --templates $(for R in $T; do echo $PWD/$N/$A/esc_stage/$R/templates_$R.npz; done)
# run_chain.sh then moves esc_stage/$R/* to $N/$A/$R/ and ampcorr.json to $N/$A/ampcorr.json.

# check: held-out gates G1-4..G1-9 and the PASS 1 cross-check -> $N/$A/$R/check/{gates.json,gates.tsv,g14-g19 png/tsv} (run_chain.sh:24)
python3 spline_template_extract.py check --pulses "$PWD/$N/$A/$R/pulses/*.npz" --templates $PWD/$N/$A/$R/templates_$R.npz --out $PWD/$N/$A/$R/check/
# plot: approximate re-render of the module pages from the bank npz, no rebuild (no pooled-point scatter, grid-interpolated bins; the build PDF is the reference; covers kept from the existing --out; writes onset_abs_ns.json)
python3 spline_template_extract.py plot --templates $PWD/$N/$A/$R/templates_$R.npz --out $PWD/$N/$A/$R/templates_$R.pdf
```

One-pulse variant (3p49GeV_H2 only, f(A) = 1; inputs are the 28 `p1shape_r3/3p49GeV_H2` splits of
`select_cmds.txt`; output `output/templates_wassim1_20261104/3p49GeV_H2/templates_3p49GeV_H2.{npz,pdf,root}`):

```bash
# driver as run (output/templates_wassim1_20261104/3p49GeV_H2/README.md); --workers is its only option (default 6)
MPLBACKEND=Agg nice -n 10 python3 wassim_one_pulse_template.py --workers 6
# what it issues per split, then once in-process (wassim_one_pulse_template.py:48-49, 69-70)
python3 spline_template_extract.py select <evio> --trigger physics --wassim-ref --out output/templates_wassim1_20261104/3p49GeV_H2/pulses/<run>_<split>.npz
python3 spline_template_extract.py build --single-pulse --pulses "output/templates_wassim1_20261104/3p49GeV_H2/pulses/*.npz" --run-type 3p49GeV_H2 --out output/templates_wassim1_20261104/3p49GeV_H2/templates_3p49GeV_H2
```

`--wassim-ref` cannot be combined with `--per-module-best`, `--mf-strict`, `--max-modules-per-event`,
`--emax-frac`, `--exclude-trigger-bits`, `--mf-bank` or `--p1-shape-cut`; `--single-pulse` cannot be combined with
`--pool-order`, `--lambda-select`, `--osc-escalate`, `--bands`, `--iterations`, `--converge-tol`,
`--max-iterations` or `--fit-points`.

## 2. Matched filter (`mf_pulse.py`)

`mf_pulse.py` has no CLI (`python3 mf_pulse.py --self-test` runs the self-test; no argument prints the
docstring). API, as called in `module_timing_dump.py:111-119`:

```python
bank, info = mf_pulse.build_bank([m.name for m in hy], pwo4, tmpl,
    dict(smooth_order=int(wcfg.smooth_order), int_tail_ratio=float(wcfg.int_tail_ratio), tail_break_n=int(wcfg.tail_break_n)),
    column=mf_pulse.resolve_column(tmpl, None, prefer="25320"), int_tail_ratio=None)   # tmpl: .npz (spline bank) or two-tau table
rec = mf_pulse.MFReconstructor(bank, nsigma=5.0, min_amp=5.0, kmax=4, r_model=0.3, ped_nsamples=int(wcfg.ped_nsamples))
res = rec.run(Y, midx)          # Y[n, ns] ADC samples, midx[n] module indices -> MFResult
```

The constructor also takes `noise_model="mad"` (choices `"mad"`, `"nps"`, `"nps_ped"`); no script or flag sets it.

On data the MF is driven by:

```bash
# module timing dump: one npz per split -> output/module_timing_ri_20261104/dump/<run6>_<split5>.npz (+ dump/logs/ with --parallel)
# (reconstructed: the line was not recorded; dump_run.log shows 17 splits of 26138/26230, templates_spline_best50_20261026 bank, parallel driver)
python3 module_timing_dump.py --all --runs 26138 26230 --parallel 2 --templates-root output/templates_spline_best50_20261026 --out output/module_timing_ri_20261104/dump
# splits that already have their npz are skipped, so rerunning the same line resumes.
```

In the benchmark (`benchmark_hycal_timing_test_20261021.py`) the switch is `--pulse-method mf` with
`--mf-templates PATH` (a `spline_templates_v1` .npz, or a two-tau table), `--mf-scale` (1.010), `--mf-r-model` (0.3), `--mf-nsigma` (5.0), `--mf-min-amp` (5.0),
`--mf-kmax` (4), `--mf-template-column`, `--mf-int-tail-ratio`, `--mf-ped-floor`, `--mf-shape`,
`--mf-nb-amp`, `--single-pulse {off,seed,cluster}`. The ri chains use
`--pulse-method mf --mf-scale 1.01 --mf-r-model 0.3 --mf-templates <bank>.npz --mf-scale 1.0` on every job
line (section 3); banks: `output/templates_spline_best50_20261026/<type>/templates_<type>.npz` for H2 and C,
`output/templates_spline_ta_fix_20261028/n200/templates_2p24GeV_Ta_X17.npz` for Ta.

## 3. Energy-resolution chain (reference R5nr-ri, `output/r5ri_20261104/tools/`)

Stages `dump -> gains -> map -> maps -> final -> repro`. Gates G0-G6 are defined in
`plans/plan_ri_20261104.md`; G2-G4 and G6 are checked inline by `chain_ri.sh` (`GATE FAILED (<chain>): ...`
in `logs/chain.out`, exit 1), G5 is not. Job lists: dump 119, map 118, final 118, repro 1.

```bash
# fresh chain dir (jobs/ and logs/P are not created by the tools; same for r7ri and fAri below), then the
# job lists -> output/<chain>_20261104/jobs/{dump,map,final,repro}.txt; --chain r5ri|r7ri|fAri
mkdir -p output/r5ri_20261104/{jobs,logs} && echo 4 > output/r5ri_20261104/logs/P
python3 output/r5ri_20261104/tools/ri_jobs.py --chain r5ri
# whole chain, one stage at a time; P (parallel jobs) is re-read from output/<chain>_20261104/logs/P before each stage
nohup bash output/r5ri_20261104/tools/drive_ri.sh r5ri > output/r5ri_20261104/logs/chain.out 2>&1 &
# single stage(s): STAGES is a subset of "dump gains map maps final repro"
P=4 STAGES=final bash output/r5ri_20261104/tools/chain_ri.sh r5ri
# resume a stopped chain: rerun the same drive_ri.sh/chain_ri.sh line. Job stages re-queue only dirs whose
# benchmark.log lacks rc=0; gains skips runs with gains/<run>/gains_<run>_all.tsv; maps always reruns.
# one job by hand (DIR CMD...): run1_ri.sh checks the benchmark sha (rc=97 on mismatch), then run1.sh writes DIR/job.cmd and DIR/benchmark.log (rc= line)
bash output/r5ri_20261104/tools/run1_ri.sh $(sed -n 1p output/r5ri_20261104/jobs/final.txt)
# tables after the three chains: tables/ri.{md,tsv}, ri_splits.tsv (pairs scans per (run, split)), then tables/ri_table.md
bash output/r5ri_20261104/tools/post.sh
```

`post.sh` runs `python3 compare_scan_dirs.py --scans R5nr-ri=$O/scan/splits R7ri=output/r7ri_20261104/scan/splits fAri=output/fAri_20261104/scan/splits --gauss-window 2sigma --out $O/tables/ri`
then `python3 $O/tools/ri_table.py $O/tables/ri_splits.tsv $O/tables/ri_table.md` (third argument defaults to `jobs/final.txt`).
Inside `chain_ri.sh` the gains stage runs `percrystal_calib.py fit --dumps "$D" --splits all --fixed-windows --out $O/gains_pass0/$R/`
then the same with `--moller-reject kin --moller-kin-n 3 --moller-kin-m 2 --moller-gains $O/gains_pass0/$R/gains_${R}_all.tsv --out $O/gains/$R/`;
the maps stage runs `position_corr.py map --scan-dir $S --out $M --splits all --workers 4` and `position_corr.py fit --map $M/map_${s}_all.tsv --form bilinear` (plus `fit --per-module ... --fallback` for 3p49GeV_H2, 2p24GeV_H2, 2p24GeV_C).

Full final-stage job line (`jobs/final.txt:1`; the first token is the job DIR, the rest is CMD):

```bash
output/r5ri_20261104/scan/splits/scan_024386/split_00010 /Users/L00338853/Software/opt/miniforge3/bin/python3 /Users/L00338853/work/PRad/prad2evviewer/analysis/pyscripts/output/merge_check_20261028/tools/merge_cfg_wrap.py --hycal-cfg output/r5_20261028/configs/reconstruction_config_r5.json /Users/L00338853/work/PRad/data/evio/prad_024386.evio.00010 output/r5ri_20261104/scan/splits/scan_024386/split_00010/result --calibrate --exclude-rowcol 15 20 15 20 --pre-window 100.0 200.0 --inner-radius 200 --min-per-module 30 --lms-corr off --hc-calib-file /Users/L00338853/work/PRad/prad2evviewer/database/calibration/calibration_factor_3p5_June7_nonLinear.json --window 16.0 --time-calib off --pulse-method mf --mf-scale 1.01 --mf-r-model 0.3 --mf-templates output/templates_spline_best50_20261026/3p49GeV_H2/templates_3p49GeV_H2.npz --mf-scale 1.0 --gain-file output/r5ri_20261104/gains/024386/gains_024386_all.tsv --pos-corr output/r5ri_20261104/poscorr/maps/3p49GeV_H2/poscorr_r5ri_3p49GeV_H2_all.json --write-cluster-tsv --moller-reject kin --moller-kin-n 3 --moller-kin-m 2 --moller-reject-gains on --fit-select chi2 --cb-two-peak --cluster-profile geant4 --gain-refine off --nl-members on
```

`merge_cfg_wrap.py` strips `--hycal-cfg`/`--runinfo` and imports the frozen benchmark. Flags that matter:
`--nl-members on`, `--hc-calib-file` (from runinfo via `ri_jobs.py`), `--gain-file` (gains stage output),
`--pos-corr` (maps stage output), `--mf-templates` (bank). Split pairing is not a flag: `compare_scan_dirs.py`
pairs scans per (run, split) by directory layout. Dump lines end with
`--cal-method percrystal --cluster-dump auto --cluster-profile geant4 --gain-refine off --nl-members on`
instead of the `--gain-file ... --nl-members on` tail; map lines have `--gain-file` but no `--pos-corr`.

Variants (same `drive_ri.sh`/`chain_ri.sh` unless stated; the final line differs from the one above only
in the tokens listed):

```bash
# R7ri: energy-dependent W0 via per-type config (log_weight_thres), not a benchmark flag;
#   configs from `python3 output/r7_20261104/tools/r7_jobs.py` -> output/r7_20261104/configs/hycal_<slug>.json, tables/w0.md
python3 output/r5ri_20261104/tools/ri_jobs.py --chain r7ri && nohup bash output/r5ri_20261104/tools/drive_ri.sh r7ri > output/r7ri_20261104/logs/chain.out 2>&1 &
#   final.txt:1 tokens: --hycal-cfg output/r7_20261104/configs/hycal_3p49GeV_H2.json ; DIR/result output/r7ri_20261104/scan/splits/scan_024386/split_00010 ;
#   --gain-file output/r7ri_20261104/gains/024386/gains_024386_all.tsv ; --pos-corr output/r7ri_20261104/poscorr/maps/3p49GeV_H2/poscorr_r7ri_3p49GeV_H2_all.json

# fAri: no f(A) = banks with the ampcorr_* keys removed (f(A) = 1); no benchmark flag, only the bank path changes
python3 output/fA_ablation_20261104/tools/make_banks.py && python3 output/fA_ablation_20261104/tools/verify_banks.py   # -> output/fA_ablation_20261104/banks/<slug>/templates_<slug>.npz, exit 0 on PASS
python3 output/r5ri_20261104/tools/ri_jobs.py --chain fAri && nohup bash output/r5ri_20261104/tools/drive_ri.sh fAri > output/fAri_20261104/logs/chain.out 2>&1 &
#   final.txt:1 tokens: --mf-templates output/fA_ablation_20261104/banks/3p49GeV_H2/templates_3p49GeV_H2.npz ; fAri paths for DIR/result, --gain-file, --pos-corr (poscorr_fAri_3p49GeV_H2_all.json)

# w1ri: one-pulse bank (section 1), 3p49GeV_H2 only (28/27/27/1 jobs); own copies of the tools
mkdir -p output/w1ri_20261104/{jobs,logs} && python3 output/w1ri_20261104/tools/ri_jobs_w1.py   # no arguments; BANK = output/templates_wassim1_20261104/3p49GeV_H2/templates_3p49GeV_H2.npz
nohup env P=3 STAGES=dump bash output/w1ri_20261104/tools/chain_w1.sh w1ri > output/w1ri_20261104/logs/chain.out 2>&1 &   # dump stage (reconstructed launch)
nohup bash output/w1ri_20261104/tools/after_dump.sh > /dev/null 2>&1 &   # waits for CHAIN w1ri DONE in logs/chain.out (exits 1 on GATE FAILED), then runs gains map maps final repro at P=4
#   final.txt:1 tokens: --mf-templates output/templates_wassim1_20261104/3p49GeV_H2/templates_3p49GeV_H2.npz ; w1ri paths (poscorr_w1ri_3p49GeV_H2_all.json)
sed -n 1p output/{r7ri,fAri,w1ri}_20261104/jobs/final.txt             # the full lines; they differ from R5nr-ri only in the tokens above
```

`chain_r7.sh`/`chain_fa.sh` (sha guard `run1_sha.sh`) are not used by the ri variants.

## 4. Timing / pile-up tools (as run by `output/tseed_dt_ri_20261104/tools/chain.sh`)

`chain.sh` takes no arguments (`bash output/tseed_dt_ri_20261104/tools/chain.sh`): it waits for
`[driver] finished` in `output/module_timing_ri_20261104/dump_run.log` (the dump of section 2), then runs,
with `T=output/tseed_dt_ri_20261104`, `O=output/module_timing_ri_20261104`, `P=output/pileup_rates_ri_20261104`:

```bash
# tseed_dt benchmark jobs (17 raw + 17 cal splits of 26138/26230) -> $T/splits/{raw,cal}/scan_<run>/split_<n>/result*; cal.txt = raw.txt with --time-calib <db>/hycal_time_offsets/time_calib.json
cat $T/jobs/raw.txt $T/jobs/cal.txt | xargs -P 2 -L1 bash output/r5_20261028/tools/run1.sh > $T/logs/run.out 2>&1
# pile-up rates: EVIO splits -> $P/dump/<run6>_<split5>.npz, then tables/maps/report.md in $P
python3 pileup_rates.py dump /Users/L00338853/work/PRad/data/evio/prad_026138.evio.00099 /Users/L00338853/work/PRad/data/evio/prad_026138.evio.00299 /Users/L00338853/work/PRad/data/evio/prad_026138.evio.00399 /Users/L00338853/work/PRad/data/evio/prad_026230.evio.00100 /Users/L00338853/work/PRad/data/evio/prad_026230.evio.00600 /Users/L00338853/work/PRad/data/evio/prad_026230.evio.01300 --parallel 2 --max-events 60000 --out $P/dump > $P/dump.log 2>&1
python3 pileup_rates.py analyse --dump-dir $P/dump --out $P --max-events-note 60000 > $P/analyse.log 2>&1
#   (chain.sh writes `--max-events 60000` here; analyse only has --max-events-note and accepts the abbreviation)
# module timing plots from the section-2 dump -> module_timing{,_tcal,_nogate,_typealign}.pdf, module_timing_veto{,_tcal,_typealign}.tsv, module_timing_nogate.tsv, veto_summary*.tsv, module_timing_vs_integral{,_V4}.{pdf,tsv}, *_veto_offsets{,_V4}.json
A="--dump-dir $O/dump --out-dir $O --jobs 2"
python3 plot_module_timing.py $A > $O/run_veto.log 2>&1
python3 plot_module_timing.py $A --time-calib > $O/run_veto_tcal.log 2>&1
python3 plot_module_timing.py $A --veto-gate none > $O/run_nogate.log 2>&1
python3 plot_module_timing.py $A --veto-align type > $O/run_veto_typealign.log 2>&1
python3 plot_module_timing_vs_integral.py $A > $O/run_vs_integral.log 2>&1
python3 plot_module_timing_vs_integral.py $A --veto-counter V4 > $O/run_vs_integral_V4.log 2>&1
# seed-time dt analysis (no arguments, paths hard-coded) -> $T/figs/tseed_dt.pdf, $T/tables/tseed_dt.md, $T/cache/*.npz
python3 $T/tools/tseed_dt.py > $T/logs/tseed_dt.log 2>&1
```

```bash
# benchmark_hycal_timing_scan.py: not used by the 20261104 chains; derived from its argparse, not from a driver (it calls resolve_calib(run) per run)
python3 benchmark_hycal_timing_scan.py --output-dir <DIR> --parallel 2 --pulse-method mf --mf-templates-map <PATTERN>
```

## 5. Non-linearity tools (`output/nonlin_20261104/tools`, no arguments except step3_load)

```bash
python3 output/nonlin_20261104/tools/step1_june7_factors.py > output/nonlin_20261104/step1/june7_factors.txt   # June7 non-linear factors to stdout (redirect reconstructed)
nice -n 10 python3 output/nonlin_20261104/tools/step2_ratios.py   # reads output/r5_20261028/gains -> step2/{step2_numbers.txt,ratio_vs_radius.png,ratio_map_xy.png,ratio_spread_pull.png,seed_vs_gain_and_lms.png}
python3 output/nonlin_20261104/tools/step3_load.py 024386 024503 024512 025228 025257 025260 025299 025308 025320 024590 024603   # (reconstructed run list from step3/reduced_*.npz) -> step3/reduced_<run>.npz; then step3_beamoffset.py, step3_tables.py, step3_fits.py (-> beam_offsets.tsv, binned_by_r.tsv, *_slopes.tsv, fit_summary.txt, fig_*.png)
```
