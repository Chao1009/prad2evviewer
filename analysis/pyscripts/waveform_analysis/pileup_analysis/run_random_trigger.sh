#!/bin/bash
# random_trigger_study.py: scan of the local EVIO splits of run 24578 (random trigger, all events of each split) and of two splits of run 24590
# (production, first 60000 events each; same code, same selection), 4 jobs in parallel; existing outputs are skipped; then the report (PDF + TSV in this directory).
HERE=$(cd "$(dirname "$0")" && pwd)
PYDIR=$(cd "$HERE/../.." && pwd)                    # analysis/pyscripts
DATA_DIR=/Users/L00338853/work/PRad/data/evio
OUT_DIR=$PYDIR/output/random_trigger_24578
NPROC=4
STAMP=20261117
SPLITS_24590=(00010 00060)
source "$PYDIR/setup_env.sh"
export MPLBACKEND=Agg PYDIR DATA_DIR OUT_DIR HERE
mkdir -p "$OUT_DIR/scan" "$OUT_DIR/ref24590" "$OUT_DIR/logs"
{ for f in "$DATA_DIR"/prad_024578.evio.*; do echo "24578 $f $OUT_DIR/scan/split_${f##*.}.npz 0 --ref-from $PYDIR/output/pileup_rates/scan/run_24590.npz"; done
  for s in "${SPLITS_24590[@]}"; do echo "24590 $DATA_DIR/prad_024590.evio.$s $OUT_DIR/ref24590/split_$s.npz 60000 --max-events 60000"; done; } |
  xargs -P "$NPROC" -L 1 sh -c '[ -s "$2" ] || python3 "$HERE/random_trigger_study.py" scan "$2" "$0" "$1" --max-events "$3" $4 $5 > "$OUT_DIR/logs/$(basename "$2" .npz)_$0.log" 2>&1' 
python3 "$HERE/random_trigger_study.py" report "$OUT_DIR" "$HERE/random_trigger_24578_$STAMP" > "$OUT_DIR/logs/report.log" 2>&1
echo DONE $(date)
