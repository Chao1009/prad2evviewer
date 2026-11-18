#!/bin/bash
# Full raw-sum scan: as run_pileup_scan.sh (13 runs, all EVIO splits of each, one process per run, 8 runs in parallel; no report step) with PILEUP_RAWSUM=1.
# Output: $OUT_DIR/rawsum_full/run_<run>.npz, logs in $OUT_DIR/rawsum_full/logs.  Progress: bash waveform_analysis/pileup_analysis/rawsum_full_status.sh
# Then: python3 waveform_analysis/pileup_analysis/integral_vs_rawsum_full.py <that dir> --stamp all_<date>
HERE=$(cd "$(dirname "$0")" && pwd)
PYDIR=$(cd "$HERE/../.." && pwd)                    # analysis/pyscripts
DATA_DIR=${DATA_DIR:-/Users/L00338853/work/PRad/data/evio}
OUT_DIR=${OUT_DIR:-$PYDIR/output/pileup_rates}/rawsum_full
NPROC=${NPROC:-8}
RUNS=(25228 25800 24698 24911 24590 25257 25260 25299 25320 24603 24937 26138 26230)   # the run_pileup_scan.sh set
export MPLBACKEND=Agg PYDIR DATA_DIR OUT_DIR PILEUP_RAWSUM=1
source "$PYDIR/setup_env.sh"
mkdir -p "$OUT_DIR/logs"
printf '%s\n' "${RUNS[@]}" | xargs -P "$NPROC" -n 1 sh -c 'files=$(ls "$DATA_DIR"/prad_0"$1".evio.*); python3 "$PYDIR/waveform_analysis/pileup_analysis/pileup_rates.py" scan "$OUT_DIR/run_$1.npz" $files > "$OUT_DIR/logs/run_$1.log" 2>&1' _
echo DONE $(date)
