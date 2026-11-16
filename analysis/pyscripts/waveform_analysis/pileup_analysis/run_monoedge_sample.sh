#!/bin/bash
# Bounded old-vs-new sample for the monotonic-edge rule and the noisy flag (pileup_rates.py options PILEUP_MONO_EDGE, PILEUP_NOISE_FLAG).
# 2 EVIO splits per run of run_pileup_scan.sh (choice: pileup_monoedge.py splits, seed fixed), the same splits for both variants:
#   old = rule off, noisy flag on;  new = rule and noisy flag on.   Output: $OUT_DIR/monoedge_sample/{old,new}/run_<run>.npz
# Then: python3 waveform_analysis/pileup_analysis/pileup_monoedge.py report
HERE=$(cd "$(dirname "$0")" && pwd)
PYDIR=$(cd "$HERE/../.." && pwd)                    # analysis/pyscripts
DATA_DIR=${DATA_DIR:-/Users/L00338853/work/PRad/data/evio}
OUT_DIR=${OUT_DIR:-$PYDIR/output/pileup_rates}/monoedge_sample
NPROC=${NPROC:-16}
source "$PYDIR/setup_env.sh"
export MPLBACKEND=Agg PYDIR OUT_DIR
mkdir -p "$OUT_DIR/old" "$OUT_DIR/new" "$OUT_DIR/logs"
python3 "$HERE/pileup_monoedge.py" splits --evio-dir "$DATA_DIR" > "$OUT_DIR/splits.txt"
# job lines "<variant> <run> <file> <file>", the heaviest runs first
awk '{for (v = 0; v < 2; v++) print (v ? "new" : "old"), $0}' "$OUT_DIR/splits.txt" | sort -k2,2r | \
  xargs -P "$NPROC" -L 1 sh -c 'v=$0; run=$1; shift; if [ "$v" = new ]; then M=1; else M=0; fi
    PILEUP_MONO_EDGE=$M PILEUP_NOISE_FLAG=1 python3 "$PYDIR/waveform_analysis/pileup_analysis/pileup_rates.py" scan "$OUT_DIR/$v/run_$run.npz" "$@" > "$OUT_DIR/logs/${v}_$run.log" 2>&1'
echo DONE $(date)
