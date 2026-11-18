#!/bin/bash
# P4 of the exotic-events plan: full scan of the 13 runs of run_rawsum_full.sh (all EVIO splits), default fit, PILEUP_EXOTIC=1 only
# (PILEUP_BASE_SLOPE, PILEUP_CLIP_SINGLE, PILEUP_MONO_EDGE, PILEUP_NOISE_FLAG, PILEUP_RAWSUM, SHAPE_SCALE unset).
# One task per EVIO file (the noise reference is per file, so rates equal the per-run scan; balances the load), Ta X17 files first (longest).
# Output: $OUT_DIR/exotic_full/run_<run>_s<split>.npz, logs in $OUT_DIR/exotic_full/logs.  Already finished files are skipped (restartable).
# Disk note: the run was continued by run_exotic_full_rest.sh (LEAN: each finished npz reduced to run_<run>_s<split>.p4.pkl by exotic_p4_module_rates.py --reduce, npz removed; only 25800 s100-s600
# and 26230 s100-s1300 kept as npz) because the volume had ~10 GB free.  Progress: bash exotic_full_status.sh.  Then: python3 exotic_p4_module_rates.py
HERE=$(cd "$(dirname "$0")" && pwd)
PYDIR=$(cd "$HERE/../../.." && pwd)                 # analysis/pyscripts
PADIR=$(cd "$HERE/.." && pwd)                       # pileup_analysis
DATA_DIR=${DATA_DIR:-/Users/L00338853/work/PRad/data/evio}
OUT_DIR=${OUT_DIR:-$PYDIR/output/pileup_rates}/exotic_full
NPROC=${NPROC:-16}
RUNS=(26230 25800 26138 24698 24911 24937 25228 25257 25260 25299 25320 24590 24603)   # Ta X17 first, then H2 (the run_rawsum_full.sh set)
source "$PYDIR/setup_env.sh"
export MPLBACKEND=Agg PYDIR OUT_DIR PADIR PILEUP_EXOTIC=1
unset PILEUP_BASE_SLOPE PILEUP_CLIP_SINGLE PILEUP_RAWSUM PILEUP_MONO_EDGE PILEUP_NOISE_FLAG PILEUP_SHAPE_SCALE
mkdir -p "$OUT_DIR/logs"
for r in "${RUNS[@]}"; do ls "$DATA_DIR"/prad_0"$r".evio.*; done > "$OUT_DIR/files.txt"
echo "$(wc -l < "$OUT_DIR/files.txt") files, NPROC=$NPROC, start $(date)"
xargs -P "$NPROC" -n 1 sh -c 'f=$0; b=$(basename "$f"); run=$(echo "$b" | sed "s/prad_0*\([0-9]*\)\.evio\.0*\([0-9]*\)/\1/"); sp=$(echo "$b" | sed "s/prad_0*\([0-9]*\)\.evio\.0*\([0-9]*\)/\2/"); o="$OUT_DIR/run_${run}_s${sp:-0}.npz"; [ -s "$o" ] && exit 0; python3 "$PADIR/pileup_rates.py" scan "$o.part.npz" "$f" > "$OUT_DIR/logs/run_${run}_s${sp:-0}.log" 2>&1 && mv "$o.part.npz" "$o"' < "$OUT_DIR/files.txt"
echo DONE $(date)
