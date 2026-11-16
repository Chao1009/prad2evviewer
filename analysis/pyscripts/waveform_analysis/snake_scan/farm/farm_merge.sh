#!/bin/bash
# farm_merge.sh [farm_merge.py options]   -- README step 5/6: farm_merge.py in the environment of farm_env.sh
# (PRAD2_ROOT, PRAD2_BINDINGS, PYTHON, PY_SETUP ... as farm_job.sh and check_build.sh).
FARM_DIR="${FARM_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
[ -f "$FARM_DIR/farm_env.sh" ] && source "$FARM_DIR/farm_env.sh"
PRAD2_ROOT="${PRAD2_ROOT:-$(cd "$FARM_DIR/../../../../.." && pwd)}"
PRAD2_DATABASE_DIR="${PRAD2_DATABASE_DIR:-$PRAD2_ROOT/database}"
PRAD2_BINDINGS="${PRAD2_BINDINGS:-$PRAD2_ROOT/build/python}"
PYTHON="${PYTHON:-python3}"
PY_SETUP="${PY_SETUP:-}"
set -eu
[ -n "$PY_SETUP" ] && eval "$PY_SETUP"
export PRAD2_DATABASE_DIR PRAD2_ROOT PYTHONPATH="$PRAD2_BINDINGS${PYTHONPATH:+:$PYTHONPATH}" MPLBACKEND=Agg
exec "$PYTHON" "$FARM_DIR/farm_merge.py" "$@"
