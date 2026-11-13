#!/bin/bash
# energy_resolution_benchmark.py: one `scan` per run (all its EVIO splits, NPROC runs in parallel), then one
# `report` per run type over the time gates GATES (the report reuses the scans; rerun it alone for other gates).
# Report order: 2p24GeV_H2 before 2p24GeV_Ta_X17 (the Ta position correction borrows the 2.24 GeV H2 polynomials and map).
# ifarm: set DATA_DIR (and the ls pattern if the file names differ); with swarm/slurm submit
# one job per run, each job being the one python call of the xargs line for that run; run the
# report step once all jobs have finished.
HERE=$(cd "$(dirname "$0")" && pwd)
PYDIR=$(cd "$HERE/.." && pwd)                    # analysis/pyscripts
DATA_DIR=/Users/L00338853/work/PRad/data/evio
OUT_DIR=$PYDIR/output/energy_resolution         # scan/<RT>/run_<run>.npz, logs/; tables and PDF: $HERE/energy_resolution/<RT>/
NPROC=8
GATES="2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10 10.5 11 11.5 12 12.5 13 13.5 14 14.5 15 15.5 16 16.5 17 17.5 18 18.5 19 19.5 20 20.5 21 21.5 22 22.5 23 23.5 24"   # time-gate scan: 2–24 ns in 0.5 ns steps
RUN_TYPES=(0p73GeV_H2 2p24GeV_H2 3p49GeV_H2 3p49GeV_H2_old 2p24GeV_Ta_X17)
RUNS_0p73GeV_H2=(24590 24603)
RUNS_2p24GeV_H2=(25228 25257 25260 25299 25320)
RUNS_3p49GeV_H2=(24698 24911 24937)
RUNS_3p49GeV_H2_old=(24386 24503 24512)         # previous 3.49 GeV runs (E_beam 3488.43 MeV); scan uses the 3p49GeV_H2 bank
RUNS_2p24GeV_Ta_X17=(25800 26138 26230)
source "$PYDIR/setup_env.sh"
export MPLBACKEND=Agg PYDIR DATA_DIR OUT_DIR
mkdir -p "$OUT_DIR/logs"
for rt in "${RUN_TYPES[@]}"; do mkdir -p "$OUT_DIR/scan/$rt"; done
for rt in "${RUN_TYPES[@]}"; do
    runs="RUNS_${rt}[@]"
    for run in "${!runs}"; do printf '%s %s\n' "$rt" "$run"; done
done | xargs -P "$NPROC" -n 2 sh -c 'files=$(ls "$DATA_DIR"/prad_0"$2".evio.*); python3 "$PYDIR/waveform_analysis/energy_resolution_benchmark.py" scan "$OUT_DIR/scan/$1/run_$2.npz" $files > "$OUT_DIR/logs/scan_$2.log" 2>&1' _
for rt in "${RUN_TYPES[@]}"; do
    python3 "$PYDIR/waveform_analysis/energy_resolution_benchmark.py" report "$OUT_DIR/scan/$rt" $GATES > "$OUT_DIR/logs/report_$rt.log" 2>&1
done
python3 "$PYDIR/waveform_analysis/energy_resolution_benchmark.py" summary "$HERE/energy_resolution" > "$OUT_DIR/logs/summary.log" 2>&1
echo DONE $(date)
