#!/bin/bash
# farm_job.sh EVIO OUTDIR   -- one EVIO split -> OUTDIR/<run>/<file>.npz (P14.7; see README.md)
#   scan runs 23699-23723: dwell blocks + EPICS + edge events + LMS/alpha series   (farm_pass.py scan mode)
#   run 23744            : cosmic MIP hits + LMS/alpha                            (farm_pass.py cosmic mode)
# Only the small npz product is written; nothing else is kept (scratch goes to $TMPDIR and is removed).
# Extra arguments are passed to farm_pass.py (e.g. --force, --max-events N, --calib FILE).
#
# ---- USER-FILLED SETTINGS (environment variables override; edit the defaults once) -------------------------------
# Repository checkout on ifarm: this script lives in <PRAD2_ROOT>/analysis/pyscripts/waveform_analysis/snake_scan/farm
# (the whole snake_scan/ directory, _common.py and database/ must be there: see README step 0).
FARM_DIR="${FARM_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
# swif2 jobs inherit no environment: put the settings below in $FARM_DIR/farm_env.sh (see farm_env.sh.example); it is sourced first.
[ -f "$FARM_DIR/farm_env.sh" ] && source "$FARM_DIR/farm_env.sh"
PRAD2_ROOT="${PRAD2_ROOT:-$(cd "$FARM_DIR/../../../../.." && pwd)}"
PRAD2_DATABASE_DIR="${PRAD2_DATABASE_DIR:-$PRAD2_ROOT/database}"
# directory holding prad2py*.so built for the Python below (this Mac: $PRAD2_ROOT/build/python, as setup_env.sh)
PRAD2_BINDINGS="${PRAD2_BINDINGS:-$PRAD2_ROOT/build/python}"
PYTHON="${PYTHON:-python3}"                 # needs numpy and scipy (scipy only in farm_merge.py)
PY_SETUP="${PY_SETUP:-}"                    # optional line evaluated first, e.g. 'source /work/.../venv/bin/activate' or 'module load python'
# -----------------------------------------------------------------------------------------------------------------
set -euo pipefail
if [ $# -lt 2 ]; then echo "usage: $0 EVIO OUTDIR [farm_pass.py options]" >&2; exit 2; fi
EVIO="$1"; OUTDIR="$2"; shift 2
[ -n "$PY_SETUP" ] && eval "$PY_SETUP"
export PRAD2_DATABASE_DIR PRAD2_ROOT
export PYTHONPATH="$PRAD2_BINDINGS${PYTHONPATH:+:$PYTHONPATH}"
export MPLBACKEND=Agg OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1   # one core per job
export TMPDIR="${TMPDIR:-/tmp}"
SS="$FARM_DIR/.."
for f in "$FARM_DIR/farm_pass.py" "$SS/snake_scan.py" "$SS/edge_pass.py" "$SS/lms_pass.py" "$SS/cosmic_mip.py" \
         "$SS/../../_common.py" "$PRAD2_DATABASE_DIR/runinfo/general.json" "$FARM_DIR/hycal_map_til_run023849.json"; do
    [ -e "$f" ] || { echo "[farm_job] ERROR: missing $f (PRAD2_ROOT=$PRAD2_ROOT)" >&2; exit 2; }
done
[ -e "$EVIO" ] || { echo "[farm_job] ERROR: $EVIO not found (not on /cache yet? run jcache.sh --check)" >&2; exit 2; }
echo "[farm_job] $(date -Is 2>/dev/null || date) host $(hostname) evio $EVIO out $OUTDIR python $($PYTHON --version 2>&1)"
exec "$PYTHON" "$FARM_DIR/farm_pass.py" "$EVIO" "$OUTDIR" --map "$FARM_DIR/hycal_map_til_run023849.json" "$@"
