#!/usr/bin/env bash
# run_trigger_scan.sh — batch fit_pulse_template.py across runs × trigger event types
#
# PURPOSE:
#   Produce one template-extraction JSON per (run, trigger-event-type) combination,
#   all using the same physics-selection cuts, for downstream comparison of how
#   template shape depends on trigger event type across runs.
#
#   Filtering uses --trigger-event-type (fadc_evt.info.trigger_type, a single
#   per-event uint8 TI event-type code) rather than the trigger bitmask.
#   Known event-type names: SSP_RawSum, Pulser, LMS, Alpha, MasterOR, SSP_Cluster.
#
# ESTIMATED RUNTIME WARNING:
#   Each (run, trigger) combination may take 30 min – several hours depending
#   on EVIO file size and the number of pulses collected.  The full 12-
#   combination grid (3 runs × 4 triggers) can easily exceed 24 hours on a
#   single core.  With parallel execution the wall-clock time is reduced
#   proportionally to the number of jobs run simultaneously.  Test with a
#   subset first, e.g.:
#       ./run_trigger_scan.sh --dry-run --runs 025308 --triggers SSP_RawSum
#
# BASH VERSION NOTE:
#   Written to be compatible with bash 3.2 (default on macOS).  Associative
#   arrays (bash 4+) are intentionally avoided; state is tracked with plain
#   delimited strings and POSIX-compatible constructs only.
#   wait -n (bash 4.3+) is NOT used; a portable PID-polling loop is used
#   instead so the script works on macOS with bash 3.2.
#
# Usage: run_trigger_scan.sh [OPTIONS]
#
#   --runs RUN1,RUN2,...           Comma-separated run numbers
#                                  (default: 025308,025320,026138)
#   --triggers T1,T2,...           Comma-separated trigger event-type names
#                                  (passed to --trigger-event-type in fit_pulse_template.py)
#                                  (default: SSP_RawSum,Pulser,LMS,SSP_Cluster)
#   --out-dir DIR                  Output directory
#                                  (default: $HOME/work/PRad/prad2evviewer/analysis/pyscripts/output/trigger_scan/)
#   --evio-dir DIR                 EVIO input directory
#                                  (default: $HOME/work/PRad/data/evio/)
#   --dry-run                      Print commands without executing
#   --force                        Overwrite existing output (DEFAULT)
#   --no-force                     Skip combinations where output JSON exists
#   --max-events N                 Cap events per run/trigger (default: 0 = all)
#   --jobs N                       Max parallel jobs (default: total_cpus - 4, min 1)
#   -h, --help                     Show this help
#
# EXAMPLES:
#   # Dry-run a single combination
#   ./run_trigger_scan.sh --dry-run --runs 025308 --triggers SSP_RawSum
#
#   # Run only two runs and two trigger event types, unlimited events
#   ./run_trigger_scan.sh --runs 025308,025320 --triggers SSP_RawSum,Pulser
#
#   # Skip existing outputs and run only missing combinations, using 8 parallel jobs
#   ./run_trigger_scan.sh --no-force --jobs 8

set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DEFAULT_OUT_DIR="$HOME/work/PRad/prad2evviewer/analysis/pyscripts/output/trigger_scan"
DEFAULT_EVIO_DIR="$HOME/work/PRad/data/evio"
DEFAULT_RUNS="025308,025320,026138"
DEFAULT_TRIGGERS="SSP_RawSum,Pulser,LMS,SSP_Cluster"

PY_SCRIPT="$HOME/work/PRad/prad2evviewer/analysis/pyscripts/fit_pulse_template.py"
DAQ_CONFIG="$HOME/work/PRad/prad2evviewer/database/daq_config.json"
HC_MAP_FILE="$HOME/work/PRad/prad2evviewer/database/hycal_map.json"

HEIGHT_MIN=500
MODEL_ERR_FLOOR=0.03
T0_MIN=22.0

# ---------------------------------------------------------------------------
# Detect CPU count — macOS and Linux compatible
# ---------------------------------------------------------------------------
if command -v nproc &>/dev/null; then
    TOTAL_CPUS=$(nproc)
elif command -v sysctl &>/dev/null; then
    TOTAL_CPUS=$(sysctl -n hw.ncpu 2>/dev/null || echo 4)
else
    TOTAL_CPUS=4
fi
MAX_JOBS=$((TOTAL_CPUS - 4))
if [ "$MAX_JOBS" -lt 1 ]; then
    MAX_JOBS=1
fi

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
RUNS="$DEFAULT_RUNS"
TRIGGERS="$DEFAULT_TRIGGERS"
OUT_DIR="$DEFAULT_OUT_DIR"
EVIO_DIR="$DEFAULT_EVIO_DIR"
DRY_RUN=0
FORCE=1
MAX_EVENTS=0
USER_JOBS=""   # empty means "use auto-detected MAX_JOBS"

usage() {
    # Print the block of comment lines starting at "# Usage:" until the first
    # non-comment line.  Uses only basic sed/awk constructs for portability
    # with macOS BSD sed and bash 3.2.
    awk '
        /^# Usage:/   { in_usage=1 }
        in_usage && /^[^#]/ { exit }
        in_usage { sub(/^# ?/, ""); print }
    ' "$0"
    exit 0
}

while [ $# -gt 0 ]; do
    case "$1" in
        --runs)
            RUNS="$2"; shift 2 ;;
        --triggers)
            TRIGGERS="$2"; shift 2 ;;
        --out-dir)
            OUT_DIR="$2"; shift 2 ;;
        --evio-dir)
            EVIO_DIR="$2"; shift 2 ;;
        --dry-run)
            DRY_RUN=1; shift ;;
        --force)
            FORCE=1; shift ;;
        --no-force)
            FORCE=0; shift ;;
        --max-events)
            MAX_EVENTS="$2"; shift 2 ;;
        --jobs)
            USER_JOBS="$2"; shift 2 ;;
        -h|--help)
            usage ;;
        *)
            echo "[ERROR] Unknown option: $1" >&2
            usage ;;
    esac
done

# Apply user override for --jobs
if [ -n "$USER_JOBS" ]; then
    MAX_JOBS="$USER_JOBS"
    if [ "$MAX_JOBS" -lt 1 ]; then
        MAX_JOBS=1
    fi
fi

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Format seconds as HHh MMm SSs
format_elapsed() {
    local total=$1
    local h=$(( total / 3600 ))
    local m=$(( (total % 3600) / 60 ))
    local s=$(( total % 60 ))
    printf "%02dh %02dm %02ds" "$h" "$m" "$s"
}

# IFS-safe split: split_csv <string> <separator>  — prints each token on its own line
split_csv() {
    local str="$1"
    local sep="$2"
    # Replace sep with newline then echo
    printf '%s' "$str" | tr "$sep" '\n'
}

# ---------------------------------------------------------------------------
# Job-slot throttling — bash 3.2 compatible PID polling
# ---------------------------------------------------------------------------
# RUNNING_PIDS: space-separated list of active background PIDs
# RUNNING_NAMES: space-separated list of "run:trigger" tokens, parallel to PIDs
# N_RUNNING: count of currently active background jobs
RUNNING_PIDS=""
RUNNING_NAMES=""
N_RUNNING=0

# wait_for_slot — blocks until the number of running jobs drops below MAX_JOBS.
# Sets global REAPED_EXIT to the exit code of the reaped job (0 if none reaped).
wait_for_slot() {
    REAPED_EXIT=0
    while [ "$N_RUNNING" -ge "$MAX_JOBS" ]; do
        local new_pids=""
        local new_names=""
        local count=0
        local found=0
        # Walk through tracked PIDs; reap any that have finished
        local i=1
        for pid in $RUNNING_PIDS; do
            name=$(echo "$RUNNING_NAMES" | cut -d' ' -f$i)
            if kill -0 "$pid" 2>/dev/null; then
                # Still running — keep it
                new_pids="$new_pids $pid"
                new_names="$new_names $name"
                count=$((count + 1))
            else
                # Finished — reap it
                wait "$pid" 2>/dev/null
                REAPED_EXIT=$?
                found=1
            fi
            i=$((i + 1))
        done
        # Trim leading spaces
        RUNNING_PIDS="${new_pids# }"
        RUNNING_NAMES="${new_names# }"
        N_RUNNING=$count
        if [ "$found" -eq 1 ]; then
            return 0
        fi
        # No job finished yet — wait a moment before polling again
        sleep 1
    done
    return 0
}

# wait_all — wait for all remaining background jobs to finish; accumulates results.
# Appends to global N_SUCCESS, N_FAILED, FAILED_LIST based on STATUS files.
wait_all() {
    for pid in $RUNNING_PIDS; do
        wait "$pid" 2>/dev/null || true
    done
    RUNNING_PIDS=""
    RUNNING_NAMES=""
    N_RUNNING=0
}

# ---------------------------------------------------------------------------
# Prepare output directories
# ---------------------------------------------------------------------------
mkdir -p "$OUT_DIR"
mkdir -p "$OUT_DIR/logs"

LOG_DIR="$OUT_DIR/logs"

# ---------------------------------------------------------------------------
# Count total combinations
# ---------------------------------------------------------------------------
# Count runs and triggers from comma-separated lists
count_tokens() {
    local str="$1"
    # Count commas + 1
    local n
    n=$(printf '%s' "$str" | tr -cd ',' | wc -c)
    echo $(( n + 1 ))
}

N_RUNS=$(count_tokens "$RUNS")
N_TRIGGERS=$(count_tokens "$TRIGGERS")
TOTAL=$(( N_RUNS * N_TRIGGERS ))

# ---------------------------------------------------------------------------
# State tracking (bash 3.2 compatible — plain strings)
# ---------------------------------------------------------------------------
COUNT=0          # combinations launched so far (including skipped)
N_SUCCESS=0
N_SKIPPED=0
N_FAILED=0
FAILED_LIST=""   # newline-separated "run trigger exitcode logpath"

WALL_START=$(date +%s)

# ---------------------------------------------------------------------------
# Banner / setup block
# ---------------------------------------------------------------------------
echo "========================================"
echo " PRad trigger scan — fit_pulse_template"
echo "========================================"
echo "  Runs:        $RUNS"
echo "  Triggers:    $TRIGGERS"
echo "  Out dir:     $OUT_DIR"
echo "  EVIO dir:    $EVIO_DIR"
echo "  Max events:  $MAX_EVENTS"
echo "  Dry run:     $DRY_RUN"
echo "  Force:       $FORCE"
echo "  Total combos: $TOTAL"
echo "[setup] CPUs=$TOTAL_CPUS  max_parallel_jobs=$MAX_JOBS"
if [ "$DRY_RUN" -eq 1 ]; then
    echo "  [DRY RUN — commands will be printed but not executed]"
fi
echo "========================================"
echo ""

# ---------------------------------------------------------------------------
# Main loop — iterate over all (run, trigger) combinations
# ---------------------------------------------------------------------------
# Use tr to convert commas to newlines, then use while read to iterate.
# This avoids any bash 4+ array syntax.
#
# For each combination we either:
#   (dry-run)  print the command and continue
#   (live)     launch a background subshell, throttled to MAX_JOBS

while IFS= read -r RUN; do
    [ -z "$RUN" ] && continue

    while IFS= read -r TRIGGER; do
        [ -z "$TRIGGER" ] && continue

        COUNT=$(( COUNT + 1 ))
        OUT_JSON="$OUT_DIR/pulse_templates_${RUN}_${TRIGGER}_h${HEIGHT_MIN}_t0min${T0_MIN}.json"
        LOG_FILE="$LOG_DIR/${RUN}_${TRIGGER}_h${HEIGHT_MIN}_t0min${T0_MIN}.log"
        STATUS_FILE="$LOG_DIR/${RUN}_${TRIGGER}_h${HEIGHT_MIN}_t0min${T0_MIN}.status"
        PLOT_DIR="$OUT_DIR/template_plots_${RUN}_${TRIGGER}_h${HEIGHT_MIN}_t0min${T0_MIN}"

        # ---- Build command string (always, for dry-run and CMD display) ----
        CMD="PYTHONPATH=\$HOME/work/PRad/prad2evviewer/build/python \\"
        CMD="$CMD
python3 $PY_SCRIPT \\"
        CMD="$CMD
    $EVIO_DIR/prad_${RUN}.evio.* \\"
        CMD="$CMD
    -o $OUT_JSON \\"
        CMD="$CMD
    --plot-dir $PLOT_DIR \\"
        CMD="$CMD
    --max-events $MAX_EVENTS \\"
        CMD="$CMD
    --max-pulses-per-channel 0 \\"
        CMD="$CMD
    --height-min $HEIGHT_MIN \\"
        CMD="$CMD
    --model-err-floor $MODEL_ERR_FLOOR \\"
        CMD="$CMD
    --t0-min $T0_MIN \\"
        CMD="$CMD
    --trigger-event-type $TRIGGER \\"
        CMD="$CMD
    --daq-config $DAQ_CONFIG \\"
        CMD="$CMD
    --hc-map-file $HC_MAP_FILE"

        # ---- Dry-run: print and continue; never background anything ----
        if [ "$DRY_RUN" -eq 1 ]; then
            echo "------------------------------------------------------------"
            echo "[RUN $COUNT/$TOTAL] run=$RUN  trigger=$TRIGGER"
            echo "  output:  $OUT_JSON"
            echo "  plot-dir: $PLOT_DIR"
            echo "  log:     $LOG_FILE"
            echo ""
            echo "[CMD]"
            echo "$CMD"
            echo ""
            echo "[DRY RUN] (not executed)"
            continue
        fi

        # ---- Skip check (in parent, before grabbing a slot) ----
        if [ -f "$OUT_JSON" ] && [ "$FORCE" -eq 0 ]; then
            echo "[SKIP] run=$RUN trigger=$TRIGGER — $OUT_JSON (exists)"
            N_SKIPPED=$(( N_SKIPPED + 1 ))
            continue
        fi

        # ---- EVIO existence check (in parent, fast) ----
        EVIO_PATTERN="$EVIO_DIR/prad_${RUN}.evio.*"
        set +e
        EVIO_FILES=$(ls $EVIO_PATTERN 2>/dev/null | head -1)
        set -e
        if [ -z "$EVIO_FILES" ]; then
            echo "[ERROR] no EVIO files for run $RUN matching: $EVIO_PATTERN"
            N_FAILED=$(( N_FAILED + 1 ))
            FAILED_LIST="${FAILED_LIST}  $RUN $TRIGGER (EVIO_MISSING, log at $LOG_FILE)\n"
            continue
        fi

        # ---- Throttle: block until a job slot is free ----
        wait_for_slot

        # ---- Launch subshell in background ----
        # Capture loop vars so the subshell sees the right values
        _RUN="$RUN"
        _TRIGGER="$TRIGGER"
        _OUT_JSON="$OUT_JSON"
        _PLOT_DIR="$PLOT_DIR"
        _LOG_FILE="$LOG_FILE"
        _STATUS_FILE="$STATUS_FILE"
        _COUNT="$COUNT"

        (
            echo "[START] run=$_RUN trigger=$_TRIGGER  combo=$_COUNT/$TOTAL  (pid=$$)"

            START_TIME=$(date +%s)

            set +e
            PYTHONPATH="$HOME/work/PRad/prad2evviewer/build/python" \
            python3 "$PY_SCRIPT" \
                $EVIO_DIR/prad_${_RUN}.evio.* \
                -o "$_OUT_JSON" \
                --plot-dir "$_PLOT_DIR" \
                --max-events "$MAX_EVENTS" \
                --max-pulses-per-channel 0 \
                --height-min "$HEIGHT_MIN" \
                --model-err-floor "$MODEL_ERR_FLOOR" \
                --t0-min "$T0_MIN" \
                --trigger-event-type "$_TRIGGER" \
                --daq-config "$DAQ_CONFIG" \
                --hc-map-file "$HC_MAP_FILE" \
                > "$_LOG_FILE" 2>&1
            EXIT_CODE=$?
            set -e

            # If fit succeeded, also generate the 2D HyCal maps
            if [ "$EXIT_CODE" -eq 0 ]; then
                python3 "$SCRIPT_DIR/plot_template_2d_map.py" \
                    "$_OUT_JSON" \
                    --out-dir "$_PLOT_DIR" \
                    --materials PbWO4 \
                    >> "$_LOG_FILE" 2>&1 || true
                # || true: don't fail the whole combo if 2D plot has issues
            fi

            END_TIME=$(date +%s)
            ELAPSED=$(( END_TIME - START_TIME ))
            ELAPSED_FMT=$(printf "%02dh %02dm %02ds" \
                $((ELAPSED/3600)) $(((ELAPSED%3600)/60)) $((ELAPSED%60)))

            if [ "$EXIT_CODE" -eq 0 ]; then
                echo "OK" > "$_STATUS_FILE"
                echo "[OK]   run=$_RUN trigger=$_TRIGGER  elapsed=${ELAPSED_FMT}"
            else
                echo "FAIL:$EXIT_CODE" > "$_STATUS_FILE"
                echo "[FAIL] run=$_RUN trigger=$_TRIGGER  exit=$EXIT_CODE  elapsed=${ELAPSED_FMT}  log=$_LOG_FILE"
            fi

            exit "$EXIT_CODE"
        ) &

        _PID=$!
        RUNNING_PIDS="$RUNNING_PIDS $_PID"
        RUNNING_NAMES="$RUNNING_NAMES ${RUN}:${TRIGGER}"
        N_RUNNING=$(( N_RUNNING + 1 ))

    done <<EOF
$(split_csv "$TRIGGERS" ',')
EOF

done <<EOF
$(split_csv "$RUNS" ',')
EOF

# ---------------------------------------------------------------------------
# Wait for all remaining background jobs
# ---------------------------------------------------------------------------
if [ "$DRY_RUN" -eq 0 ]; then
    echo ""
    echo "[wait] All combinations launched — waiting for remaining $N_RUNNING job(s)..."
    wait_all
fi

# ---------------------------------------------------------------------------
# Tally results from status files
# ---------------------------------------------------------------------------
if [ "$DRY_RUN" -eq 0 ]; then
    while IFS= read -r RUN; do
        [ -z "$RUN" ] && continue
        while IFS= read -r TRIGGER; do
            [ -z "$TRIGGER" ] && continue

            OUT_JSON="$OUT_DIR/pulse_templates_${RUN}_${TRIGGER}_h${HEIGHT_MIN}_t0min${T0_MIN}.json"
            LOG_FILE="$LOG_DIR/${RUN}_${TRIGGER}_h${HEIGHT_MIN}_t0min${T0_MIN}.log"
            STATUS_FILE="$LOG_DIR/${RUN}_${TRIGGER}_h${HEIGHT_MIN}_t0min${T0_MIN}.status"

            # Combinations that were skipped (no status file written, output exists)
            if [ ! -f "$STATUS_FILE" ]; then
                # Could be skipped or EVIO-missing (already counted above); skip here
                continue
            fi

            STATUS_CONTENT=$(cat "$STATUS_FILE")
            if [ "$STATUS_CONTENT" = "OK" ]; then
                N_SUCCESS=$(( N_SUCCESS + 1 ))
            else
                EXIT_CODE_FROM_FILE=$(echo "$STATUS_CONTENT" | cut -d: -f2)
                N_FAILED=$(( N_FAILED + 1 ))
                FAILED_LIST="${FAILED_LIST}  $RUN $TRIGGER (exit code ${EXIT_CODE_FROM_FILE:-?}, log at $LOG_FILE)\n"
            fi

        done <<EOF2
$(split_csv "$TRIGGERS" ',')
EOF2
    done <<EOF2
$(split_csv "$RUNS" ',')
EOF2
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
WALL_END=$(date +%s)
WALL_ELAPSED=$(( WALL_END - WALL_START ))
WALL_FMT=$(format_elapsed "$WALL_ELAPSED")

echo ""
echo "==== summary ===="

if [ "$DRY_RUN" -eq 1 ]; then
    echo "Dry run: $COUNT command(s) printed (nothing executed)"
else
    echo "Successful:  $N_SUCCESS/$TOTAL combinations"
    echo "Skipped:     $N_SKIPPED combinations (existing output)"
    echo "Failed:      $N_FAILED combinations"
fi

echo "Total time:  $WALL_FMT"
echo "Max parallel jobs: $MAX_JOBS"

if [ "$DRY_RUN" -eq 0 ] && [ "$N_FAILED" -gt 0 ]; then
    echo ""
    echo "Failed combinations:"
    printf '%b' "$FAILED_LIST"
fi

echo ""
