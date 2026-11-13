#!/bin/bash
# Template candidates for the four run types: one `select` job per run, then one `build` per run type.
# ifarm: set DATA_DIR (and the ls pattern if the file names differ); with swarm/slurm submit
# one job per run, each job being the one python call of the xargs line for that run; run the
# build step once all jobs have finished.
HERE=$(cd "$(dirname "$0")" && pwd)
PYDIR=$(cd "$HERE/.." && pwd)                    # analysis/pyscripts
DATA_DIR=/Users/L00338853/work/PRad/data/evio
OUT_DIR=$HERE/pulse_templates                    # must equal OUT_ROOT of template_extract.py
NPROC=8
RUN_TYPES=(0p73GeV_H2 2p24GeV_H2 3p49GeV_H2 2p24GeV_Ta_X17)
RUNS_0p73GeV_H2=(24590 24603)
RUNS_2p24GeV_H2=(25228 25257 25260 25299 25320)
# RUNS_3p49GeV_H2=(24386 24503 24512)    # previous 3p49GeV_H2 runs (replaced 2026-11-10)
RUNS_3p49GeV_H2=(24698 24911 24937)
RUNS_2p24GeV_Ta_X17=(25800 26138 26230)
source "$PYDIR/setup_env.sh"
export MPLBACKEND=Agg PYDIR DATA_DIR OUT_DIR
for rt in "${RUN_TYPES[@]}"; do mkdir -p "$OUT_DIR/$rt/pulses" "$OUT_DIR/$rt/logs"; done
for rt in "${RUN_TYPES[@]}"; do
    runs="RUNS_${rt}[@]"
    for run in "${!runs}"; do printf '%s %s\n' "$rt" "$run"; done
done | xargs -P "$NPROC" -n 2 sh -c 'files=$(ls "$DATA_DIR"/prad_0"$2".evio.*); python3 "$PYDIR/waveform_analysis/template_extract.py" select "$OUT_DIR/$1/pulses/run_$2.npz" "$1" $files > "$OUT_DIR/$1/logs/run_$2.log" 2>&1' _
for rt in "${RUN_TYPES[@]}"; do
    python3 "$PYDIR/waveform_analysis/template_extract.py" build "$rt" > "$OUT_DIR/$rt/build.log" 2>&1
done
echo DONE $(date)
