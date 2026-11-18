#!/bin/bash
# P3 of the exotic-events plan: the 2-file sample (same splits as run_rawsum_sample.sh, same splits.txt from pileup_monoedge.py), default fit, PILEUP_EXOTIC=1 only
# (PILEUP_BASE_SLOPE, PILEUP_CLIP_SINGLE, PILEUP_RAWSUM off).  Output: $OUT_DIR/exotic_sample/run_<run>.npz.  Then: python3 exotic_p3_sample.py
HERE=$(cd "$(dirname "$0")" && pwd)
PYDIR=$(cd "$HERE/../../.." && pwd)                 # analysis/pyscripts
PADIR=$(cd "$HERE/.." && pwd)                       # pileup_analysis
DATA_DIR=${DATA_DIR:-/Users/L00338853/work/PRad/data/evio}
OUT_DIR=${OUT_DIR:-$PYDIR/output/pileup_rates}/exotic_sample
NPROC=${NPROC:-8}
source "$PYDIR/setup_env.sh"
export MPLBACKEND=Agg PYDIR OUT_DIR PADIR PILEUP_EXOTIC=1
unset PILEUP_BASE_SLOPE PILEUP_CLIP_SINGLE PILEUP_RAWSUM PILEUP_MONO_EDGE PILEUP_NOISE_FLAG
mkdir -p "$OUT_DIR/logs"
python3 "$PADIR/pileup_monoedge.py" splits --evio-dir "$DATA_DIR" > "$OUT_DIR/splits.txt"
sort -k1,1r "$OUT_DIR/splits.txt" | \
  xargs -P "$NPROC" -L 1 sh -c 'run=$0; python3 "$PADIR/pileup_rates.py" scan "$OUT_DIR/run_$run.npz" "$@" > "$OUT_DIR/logs/$run.log" 2>&1'
echo DONE $(date)
