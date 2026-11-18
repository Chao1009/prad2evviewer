#!/bin/bash
# Continuation of run_exotic_full.sh with LEAN=1 (disk space): launches the files of $OUT_DIR/exotic_full/files.txt that have no log yet, keeping at most NPROC scan
# processes running in total; after each scan the npz is reduced to run_<run>_s<split>.p4.pkl (exotic_p4_module_rates.py --reduce) and, once that pickle loads, the npz is removed.
HERE=$(cd "$(dirname "$0")" && pwd)
PYDIR=$(cd "$HERE/../../.." && pwd); PADIR=$(cd "$HERE/.." && pwd)
OUT_DIR=${OUT_DIR:-$PYDIR/output/pileup_rates}/exotic_full
NPROC=${NPROC:-16}
source "$PYDIR/setup_env.sh" > /dev/null
export MPLBACKEND=Agg PYDIR OUT_DIR PADIR HERE PILEUP_EXOTIC=1
unset PILEUP_BASE_SLOPE PILEUP_CLIP_SINGLE PILEUP_RAWSUM PILEUP_MONO_EDGE PILEUP_NOISE_FLAG PILEUP_SHAPE_SCALE
task() {
  f=$1; b=$(basename "$f"); run=$(echo "$b" | sed "s/prad_0*\([0-9]*\)\.evio\.0*\([0-9]*\)/\1/"); sp=$(echo "$b" | sed "s/prad_0*\([0-9]*\)\.evio\.0*\([0-9]*\)/\2/"); sp=${sp:-0}
  o="$OUT_DIR/run_${run}_s${sp}.npz"
  python3 "$PADIR/pileup_rates.py" scan "$o.part.npz" "$f" > "$OUT_DIR/logs/run_${run}_s${sp}.log" 2>&1 && mv "$o.part.npz" "$o" \
    && python3 "$HERE/exotic_p4_module_rates.py" --reduce "$o" >> "$OUT_DIR/logs/run_${run}_s${sp}.log" 2>&1 \
    && python3 -c "import pickle,sys; pickle.load(open(sys.argv[1],'rb'))" "${o%.npz}.p4.pkl" && rm -f "$o"
}
while read -r f; do
  b=$(basename "$f"); run=$(echo "$b" | sed "s/prad_0*\([0-9]*\)\.evio\.0*\([0-9]*\)/\1/"); sp=$(echo "$b" | sed "s/prad_0*\([0-9]*\)\.evio\.0*\([0-9]*\)/\2/"); sp=${sp:-0}
  [ -e "$OUT_DIR/logs/run_${run}_s${sp}.log" ] && continue
  while [ "$(pgrep -f 'pileup_rates.py scan .*exotic_full' | wc -l)" -ge "$NPROC" ]; do sleep 20; done
  task "$f" &
  sleep 5
done < "$OUT_DIR/files.txt"
wait
while [ "$(pgrep -f 'pileup_rates.py scan .*exotic_full' | wc -l)" -gt 0 ]; do sleep 30; done   # tasks of the original xargs driver
echo DONE $(date)
