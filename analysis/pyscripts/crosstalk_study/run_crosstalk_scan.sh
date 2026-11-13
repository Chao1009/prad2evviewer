#!/bin/bash
# crosstalk_study.py: one scan per run (all its EVIO splits, NPROC runs in parallel), then the report over all runs.
# ifarm: set DATA_DIR (and the ls pattern if the file names differ); with swarm/slurm submit
# one job per run, each job being the one python call of the xargs line for that run; run the
# report step once all jobs have finished.
HERE=$(cd "$(dirname "$0")" && pwd)
PYDIR=$(cd "$HERE/.." && pwd)                    # analysis/pyscripts
DATA_DIR=/Users/L00338853/work/PRad/data/evio
OUT_DIR=$PYDIR/output/crosstalk_study
NPROC=6
RUNS_0p73GeV_H2=(24590 24603)
RUNS_2p24GeV_H2=(25228 25257 25260 25299 25320)
RUNS_3p49GeV_H2=(24698 24911 24937)
RUNS_2p24GeV_Ta_X17=(25800 26138 26230)
RUNS=("${RUNS_0p73GeV_H2[@]}" "${RUNS_2p24GeV_H2[@]}" "${RUNS_3p49GeV_H2[@]}")   # add "${RUNS_2p24GeV_Ta_X17[@]}" to include X17 (not in the figures)
source "$PYDIR/setup_env.sh"
export MPLBACKEND=Agg PYDIR DATA_DIR OUT_DIR
mkdir -p "$OUT_DIR/scan" "$OUT_DIR/logs"
printf '%s\n' "${RUNS[@]}" | xargs -P "$NPROC" -n 1 sh -c 'files=$(ls "$DATA_DIR"/prad_0"$1".evio.*); python3 "$PYDIR/crosstalk_study/crosstalk_study.py" scan "$OUT_DIR/scan/run_$1.npz" $files > "$OUT_DIR/logs/run_$1.log" 2>&1' _
python3 "$PYDIR/crosstalk_study/crosstalk_study.py" report "$OUT_DIR/scan" "$HERE/out" > "$OUT_DIR/report.log" 2>&1
echo DONE $(date)
