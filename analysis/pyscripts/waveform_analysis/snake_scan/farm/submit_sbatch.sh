#!/bin/bash
# submit_sbatch.sh CACHE_LIST OUTDIR [--submit]   -- README step 3 (Slurm alternative to swif2). Run on ifarm.
# Splits the list into scan files and 23744 files (different time/memory), writes OUTDIR/slurm_scan_<tag>.list and slurm_cosmic_<tag>.list and
# the sbatch command lines in OUTDIR/sbatch_cmds_<tag>.sh (<tag> = timestamp, so a resubmission of a missing.txt does not overwrite the lists of
# running array jobs, which read their list at task start); --submit runs them. Array task k processes line k of its list;
# products that already exist are skipped (resubmit the same list to finish after failures).
set -eu
LIST="${1:?CACHE_LIST}"; OUTDIR="${2:?OUTDIR}"; DO="${3:-}"
TAG="${TAG:-$(date +%Y%m%d_%H%M%S)}"
FARM_DIR="${FARM_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
ACCOUNT="${ACCOUNT:-hallb}"; PARTITION="${PARTITION:-production}"   # defaults of scripts/shell/replay_common.sh; confirm with `sacctmgr show assoc user=jyzhou` / `sinfo`
CONC="${CONC:-100}"                                                  # max simultaneous array tasks
SCAN_TIME="${SCAN_TIME:-02:00:00}"; SCAN_MEM="${SCAN_MEM:-2000}"
COS_TIME="${COS_TIME:-01:00:00}"; COS_MEM="${COS_MEM:-1500}"
mkdir -p "$OUTDIR/logs"
grep -v 'prad_023744\.' "$LIST" > "$OUTDIR/slurm_scan_$TAG.list" || true
grep 'prad_023744\.' "$LIST" > "$OUTDIR/slurm_cosmic_$TAG.list" || true
CMDS="$OUTDIR/sbatch_cmds_$TAG.sh"; : > "$CMDS"
for kind in scan cosmic; do
    L="$OUTDIR/slurm_${kind}_$TAG.list"; n=$(wc -l < "$L" | tr -d " "); [ "$n" -gt 0 ] || continue
    if [ $kind = scan ]; then T=$SCAN_TIME; M=$SCAN_MEM; else T=$COS_TIME; M=$COS_MEM; fi
    echo "sbatch --account=$ACCOUNT --partition=$PARTITION --job-name=p14_${kind}_$TAG --array=1-$n%$CONC --cpus-per-task=1 --mem=${M}M --time=$T" \
         "--output=$OUTDIR/logs/p14_${kind}_${TAG}_%A_%a.out --export=ALL,FARM_DIR=$FARM_DIR $FARM_DIR/farm_array.sbatch $L $OUTDIR" >> "$CMDS"
done
cat "$CMDS"
[ "$DO" = "--submit" ] && bash "$CMDS" && echo "submitted; monitor: squeue -u jyzhou"
exit 0
