#!/bin/bash
# submit_swif2.sh CACHE_LIST OUTDIR [--submit]   -- README step 3 (swif2). Run on ifarm.
#   CACHE_LIST : one /cache evio path per line (jcache_list_cache.txt from make_filelist.sh; any subset, e.g. farm_missing.sh output)
#   OUTDIR     : directory for the products (<OUTDIR>/<run>/<file>.npz) and logs (<OUTDIR>/logs); must be on /work or /volatile
# Without --submit the commands are only written to $OUTDIR/swif2_jobs.sh (review, then `bash` it or rerun with --submit).
# One job per file, 1 core. The 23744 files (cosmic scan, about 5x slower per file) get a longer time/ram request.
# Resources (measured on the Mac: scan file 320-590 s wall, max RSS 676 MB, product <= 19.7 MB; 23744 file 95 s, 540 MB; requests are
# 3-4x the measurement, see README.md): override via env. Options -account/-partition as the Slurm settings of scripts/shell/replay_common.sh
# (hallb / production); the repo has no other swif2 usage, so confirm the option names and the -time/-ram/-disk syntax with
# `swif2 add-job -help` on ifarm before --submit (the default run only writes the file).
set -eu
LIST="${1:?CACHE_LIST}"; OUTDIR="${2:?OUTDIR}"; DO="${3:-}"
FARM_DIR="${FARM_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
WORKFLOW="${WORKFLOW:-p14_snake_scan}"
ACCOUNT="${ACCOUNT:-hallb}"           # swif2 -account, same default as SLURM_ACCOUNT in scripts/shell/replay_common.sh
PARTITION="${PARTITION:-production}"  # swif2 -partition, same default as SLURM_PARTITION there
SCAN_TIME="${SCAN_TIME:-2h}"; SCAN_RAM="${SCAN_RAM:-2000MB}"; SCAN_DISK="${SCAN_DISK:-2GB}"
COS_TIME="${COS_TIME:-1h}"; COS_RAM="${COS_RAM:-1500MB}"; COS_DISK="${COS_DISK:-1GB}"
mkdir -p "$OUTDIR/logs"
JOBS="$OUTDIR/swif2_jobs.sh"
{
  echo "#!/bin/bash"
  echo "# generated $(date) from $LIST"
  echo "swif2 create -workflow $WORKFLOW 2>/dev/null || true      # fine if it already exists"
  while IFS= read -r f; do
      [ -n "$f" ] || continue
      b=$(basename "$f"); r=$((10#$(echo "$b" | sed -E 's/^prad_0*([0-9]+)\..*/\1/'))); n="p14_${r}_${b##*.}"
      if [ "$r" = 23744 ]; then T=$COS_TIME; M=$COS_RAM; D=$COS_DISK; else T=$SCAN_TIME; M=$SCAN_RAM; D=$SCAN_DISK; fi
      echo "swif2 add-job -workflow $WORKFLOW -account $ACCOUNT -partition $PARTITION -cores 1 -ram $M -disk $D -time $T -name $n" \
           "-stdout $OUTDIR/logs/$n.out -stderr $OUTDIR/logs/$n.err" \
           "$FARM_DIR/farm_job.sh $f $OUTDIR"
  done < "$LIST"
  echo "swif2 run -workflow $WORKFLOW"
} > "$JOBS"
echo "wrote $(grep -c '^swif2 add-job' "$JOBS") add-job lines to $JOBS"
if [ "$DO" = "--submit" ]; then bash "$JOBS"; echo "submitted; monitor: swif2 status -workflow $WORKFLOW"; fi
