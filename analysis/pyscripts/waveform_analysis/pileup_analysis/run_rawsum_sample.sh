#!/bin/bash
# Raw-sum sample: the 2 EVIO splits per run of run_monoedge_sample.sh (same splits.txt, seed fixed), default fit (switches off), PILEUP_RAWSUM=1.
# Output: $OUT_DIR/rawsum_sample/run_<run>.npz.  Then: python3 waveform_analysis/pileup_analysis/integral_vs_rawsum_full.py <that dir>
HERE=$(cd "$(dirname "$0")" && pwd)
PYDIR=$(cd "$HERE/../.." && pwd)                    # analysis/pyscripts
DATA_DIR=${DATA_DIR:-/Users/L00338853/work/PRad/data/evio}
OUT_DIR=${OUT_DIR:-$PYDIR/output/pileup_rates}/rawsum_sample
NPROC=${NPROC:-16}
source "$PYDIR/setup_env.sh"
export MPLBACKEND=Agg PYDIR OUT_DIR PILEUP_RAWSUM=1
mkdir -p "$OUT_DIR/logs"
python3 "$HERE/pileup_monoedge.py" splits --evio-dir "$DATA_DIR" > "$OUT_DIR/splits.txt"
sort -k1,1r "$OUT_DIR/splits.txt" | \
  xargs -P "$NPROC" -L 1 sh -c 'run=$0; python3 "$PYDIR/waveform_analysis/pileup_analysis/pileup_rates.py" scan "$OUT_DIR/run_$run.npz" "$@" > "$OUT_DIR/logs/$run.log" 2>&1' 
echo DONE $(date)
