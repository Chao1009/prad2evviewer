#!/bin/bash
# submit_replay_recon_m.sh
# Submit multiple single-run PRad replay pipeline Slurm jobs.
#
# Prompts once for a run range or list and the shared settings, then submits
# one job per cached run and requests jcache staging for the others.  The CPU
# count is entered once and reused by replay_recon -j, replay_filter -t,
# quick_check -j, and Slurm --cpus-per-task in each generated job.  Each run
# writes all products directly under <output_base>/prad_<RUN>/.

set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/replay_common.sh" ||
    { echo "ERROR: replay_common.sh not found next to this script"; exit 1; }

echo "Submit multiple PRad replay/recon Slurm jobs"
echo ""

RUNS=()
read -rp "Enter start run number, or press Enter to input a run list: " START_RUN
if [[ -z "${START_RUN}" ]]; then
    read -rp "Enter run list separated by spaces or commas: " RUN_LIST
    if [[ -z "${RUN_LIST}" ]]; then
        echo "ERROR: run list cannot be empty."
        exit 1
    fi
    RUN_LIST="${RUN_LIST//,/ }"
    for run in ${RUN_LIST}; do
        if [[ ! "${run}" =~ ^[0-9]+$ ]]; then
            echo "ERROR: invalid run number: ${run}"
            exit 1
        fi
        RUNS+=("${run}")
    done
else
    if [[ ! "${START_RUN}" =~ ^[0-9]+$ ]]; then
        echo "ERROR: invalid start run number: ${START_RUN}"
        exit 1
    fi
    read -rp "Enter end run number [${START_RUN}]: " END_RUN
    END_RUN="${END_RUN:-${START_RUN}}"
    if [[ ! "${END_RUN}" =~ ^[0-9]+$ ]]; then
        echo "ERROR: invalid end run number: ${END_RUN}"
        exit 1
    fi
    if (( 10#${END_RUN} < 10#${START_RUN} )); then
        echo "ERROR: end run number must be greater than or equal to start run number."
        exit 1
    fi

    WIDTH="${#START_RUN}"
    if (( ${#END_RUN} > WIDTH )); then
        WIDTH="${#END_RUN}"
    fi
    for ((run=10#${START_RUN}; run<=10#${END_RUN}; run++)); do
        RUNS+=("$(printf "%0${WIDTH}d" "${run}")")
    done
fi

prompt_replay_settings
prompt_job_settings
submit_replay_runs "${RUNS[@]}"
