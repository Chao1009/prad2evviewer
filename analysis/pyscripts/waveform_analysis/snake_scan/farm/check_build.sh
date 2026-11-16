#!/bin/bash
# check_build.sh [EVIO] -- README step 1 (1a; with EVIO also 1b, the golden-number check on that file): are the bindings, packages, database and the derived map usable in this environment?
# Same settings as farm_job.sh (PRAD2_ROOT, PRAD2_BINDINGS, PYTHON, PY_SETUP ...). Prints OK or the first failure.
# EVIO (optional): /cache path of prad_023701.evio.00000; runs check_golden.py in this same environment.
EVIO_GOLD="${1:-}"
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
echo "PRAD2_ROOT=$PRAD2_ROOT"; echo "PRAD2_BINDINGS=$PRAD2_BINDINGS"; echo "PRAD2_DATABASE_DIR=$PRAD2_DATABASE_DIR"
"$PYTHON" - <<PY
import sys, hashlib, pathlib
print("python", sys.version.split()[0])
import numpy, scipy
print("numpy", numpy.__version__, "scipy", scipy.__version__)
from prad2py import dec, det
print("prad2py bindings: OK", dec.__name__, det.__name__)
for n in ("EvChannel", "WaveAnalyzer"):
    assert hasattr(dec, n), n
sys.path.insert(0, "$FARM_DIR/..")
sys.path.insert(0, "$FARM_DIR/../../..")
import _common, snake_scan, edge_pass, lms_pass, cosmic_mip, transfer
print("study modules import: OK")
a = pathlib.Path("$FARM_DIR/hycal_map_til_run023849.json"); b = pathlib.Path("$FARM_DIR/../hycal_map_til_run023849.json")
if b.exists():
    ha, hb = hashlib.md5(a.read_bytes()).hexdigest(), hashlib.md5(b.read_bytes()).hexdigest()
    print("derived map copy equals snake_scan/ map:", ha == hb)
    if ha != hb:
        sys.exit("FAIL: farm/hycal_map_til_run023849.json differs from snake_scan/hycal_map_til_run023849.json; copy the current one")
print("calibration file of run 23701:", edge_pass.resolve_calib(23701))
print("OK")
PY
if [ -n "$EVIO_GOLD" ]; then
    echo "--- check_golden.py $EVIO_GOLD"
    "$PYTHON" "$FARM_DIR/check_golden.py" "$EVIO_GOLD"
fi
