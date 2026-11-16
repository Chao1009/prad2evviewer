#!/bin/bash
# farm_array.sbatch LIST OUTDIR   -- Slurm array task: process line $SLURM_ARRAY_TASK_ID of LIST (a /cache evio path). Started by submit_sbatch.sh.
# FARM_DIR is exported by submit_sbatch.sh (a spooled copy of this script cannot find its own directory).
set -eu
LIST="$1"; OUTDIR="$2"
EVIO=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
[ -n "$EVIO" ] || { echo "no line ${SLURM_ARRAY_TASK_ID} in $LIST" >&2; exit 2; }
exec "${FARM_DIR:?FARM_DIR not exported}/farm_job.sh" "$EVIO" "$OUTDIR"
