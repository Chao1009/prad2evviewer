# setup_env.sh — Environment setup for PRad-II Python analysis scripts
#
# USAGE
#   source setup_env.sh        (bash / zsh)
#   . setup_env.sh             (POSIX sh)
#
# This file MUST be sourced, NOT executed.  Executing it as a standalone
# script (./setup_env.sh or bash setup_env.sh) runs it in a child shell:
# the variables are exported there, the child exits, and all exports are
# immediately discarded in the calling shell — i.e., nothing happens.
#
# Re-sourcing in the same shell is safe.  PYTHONPATH gets a duplicate
# prepend entry, but duplicate entries in PYTHONPATH are harmless.

# ---------------------------------------------------------------------------
# Guard: detect execution rather than sourcing and abort with a clear message.
# BASH_SOURCE is a bash-only array; skip the check on plain sh (sourcing
# works correctly there regardless).
# ---------------------------------------------------------------------------
if [ -n "${BASH_SOURCE+x}" ] && [ "${BASH_SOURCE[0]}" = "${0}" ]; then
    echo "ERROR: setup_env.sh must be sourced, not executed." >&2
    echo "Usage: source setup_env.sh" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Exports
# ---------------------------------------------------------------------------
export PRAD2_DATABASE_DIR=$HOME/work/PRad/prad2evviewer/database
export DAQ_CONFIG=$PRAD2_DATABASE_DIR/daq_config.json
export HC_MAP=$PRAD2_DATABASE_DIR/hycal_map.json

# Prepend to PYTHONPATH; if already set, insert a ':' separator so existing
# entries are preserved.
export PYTHONPATH=$HOME/work/PRad/prad2evviewer/build/python${PYTHONPATH:+:$PYTHONPATH}

# ---------------------------------------------------------------------------
# Confirmation
# ---------------------------------------------------------------------------
echo "[setup_env] PRAD2_DATABASE_DIR=$PRAD2_DATABASE_DIR"
echo "[setup_env] PYTHONPATH=$PYTHONPATH"

# ---------------------------------------------------------------------------
# Sanity checks — warn if expected JSON config files are absent.
# Sourcing still succeeds; the user may have a non-standard layout.
# ---------------------------------------------------------------------------
[ -f "$DAQ_CONFIG" ] || echo "[setup_env] WARNING: DAQ_CONFIG not found: $DAQ_CONFIG" >&2
[ -f "$HC_MAP" ]     || echo "[setup_env] WARNING: HC_MAP not found: $HC_MAP" >&2
