#!/bin/bash
# Noisy-baseline sample (pileup_rates.py option PILEUP_NOISY_BASE=1; MONO_EDGE and NOISE_FLAG off): the same 2 EVIO splits per run as run_monoedge_sample.sh
# (monoedge_sample/splits.txt).  Output: $OUT_DIR/noisybase_sample/run_<run>.npz.  Then: python3 pileup_noisybase.py report
HERE=$(cd "$(dirname "$0")" && pwd)
PYDIR=$(cd "$HERE/../.." && pwd)
OUT_DIR=${OUT_DIR:-$PYDIR/output/pileup_rates}
NPROC=${NPROC:-16}
source "$PYDIR/setup_env.sh"
export MPLBACKEND=Agg
mkdir -p "$OUT_DIR/noisybase_sample/logs"
sort -k1,1r "$OUT_DIR/monoedge_sample/splits.txt" | \
  xargs -P "$NPROC" -L 1 sh -c 'run=$0; PILEUP_NOISY_BASE=1 python3 "'"$HERE"'/pileup_rates.py" scan "'"$OUT_DIR"'/noisybase_sample/run_$run.npz" "$@" > "'"$OUT_DIR"'/noisybase_sample/logs/$run.log" 2>&1'
echo DONE $(date)
