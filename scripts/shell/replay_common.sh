# replay_common.sh — defaults, prompts and the replay pipeline used by the
# PRad replay scripts in this directory.  Sourced, not run; it locates the
# checkout from its own path, so it must stay next to those scripts.
#
# Everything here is safe under `set -euo pipefail` and without it.

# --- Defaults (every one can be overridden from the environment) ---
# PRAD2_SOFT defaults to the checkout holding this file.  PRAD2_BIN and
# DEFAULT_CUTS follow the prad2evviewer directory entered at the prompt unless
# they were set explicitly.
PRAD2_BIN_WAS_SET=0
DEFAULT_CUTS_WAS_SET=0
[[ -n "${PRAD2_BIN:-}" ]] && PRAD2_BIN_WAS_SET=1
[[ -n "${DEFAULT_CUTS:-}" ]] && DEFAULT_CUTS_WAS_SET=1

PRAD2_SOFT="${PRAD2_SOFT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
PRAD2_BIN="${PRAD2_BIN:-${PRAD2_SOFT}/build/bin}"
CACHE_BASE="${CACHE_BASE:-/cache/clas12/rg-o/data}"
MSS_BASE="${MSS_BASE:-/mss/clas12/rg-o/data}"
OUTPUT_BASE="${OUTPUT_BASE:-./}"
REPLAY_CORES="${REPLAY_CORES:-15}"
REPLAY_ZERO_SUPPRESS="${REPLAY_ZERO_SUPPRESS:-5}"
REPLAY_MAX_FILES="${REPLAY_MAX_FILES:-10000}"
REPLAY_MERGE_FILES="${REPLAY_MERGE_FILES:-62}"
DEFAULT_CUTS="${DEFAULT_CUTS:-${PRAD2_SOFT}/analysis/cuts/prad2_default.json}"
SLURM_ACCOUNT="${SLURM_ACCOUNT:-hallb}"
SLURM_PARTITION="${SLURM_PARTITION:-production}"
SLURM_TIME="${SLURM_TIME:-12:00:00}"
SLURM_MEM_PER_CPU="${SLURM_MEM_PER_CPU:-1500}"
ROOT_SETUP="${ROOT_SETUP:-}"

# --- Small helpers ---
prompt_default() {
    local prompt="$1"
    local default="$2"
    local value
    read -rp "${prompt} [${default}]: " value
    if [[ -n "${value}" ]]; then
        printf '%s\n' "${value}"
    else
        printf '%s\n' "${default}"
    fi
}

to_abs_path() {
    local path="$1"
    local base="$2"
    if [[ "${path}" == /* ]]; then
        printf '%s\n' "${path}"
    else
        printf '%s\n' "${base}/${path#./}"
    fi
}

shell_quote() {
    printf "'%s'" "$(printf '%s' "$1" | sed "s/'/'\\\\''/g")"
}

# True when the run's cache directory holds at least one regular file.
run_is_cached() {
    local run_dir="${CACHE_BASE}/prad_$1"
    [[ -d "${run_dir}" && -n "$(find "${run_dir}" -maxdepth 1 -type f -print -quit)" ]]
}

# --- Interactive prompts ---
prompt_run_number() {
    read -rp "Enter run number (e.g. 024650): " RUN_NUMBER
    if [[ -z "${RUN_NUMBER}" ]]; then
        echo "ERROR: run number cannot be empty."
        exit 1
    fi
}

# Replay mode, GEM hit branches, software/data paths, replay options and the
# replay_filter cut JSON.  Sets REPLAY_MODE_FLAG, REPLAY_MODE_NAME,
# REPLAY_GEM_HIT_FLAG and CUT_JSON and updates the defaults above.
prompt_replay_settings() {
    local answer
    while true; do
        read -rp "Enter replay mode (prad2, x17, x17_full, prad1, or random) [prad2]: " answer
        case "${answer,,}" in
            ""|prad2|-prad2)    REPLAY_MODE_FLAG="";                REPLAY_MODE_NAME="PRad2" ;;
            x17|-x17)           REPLAY_MODE_FLAG="-x17";            REPLAY_MODE_NAME="X17" ;;
            x17_full|-x17_full) REPLAY_MODE_FLAG="-x17 --x17_full"; REPLAY_MODE_NAME="X17 Full" ;;
            prad1|-prad1)       REPLAY_MODE_FLAG="-prad1";          REPLAY_MODE_NAME="PRad1" ;;
            random|-random)     REPLAY_MODE_FLAG="-random";         REPLAY_MODE_NAME="Random" ;;
            *) echo "ERROR: enter prad2, x17, x17_full, prad1, or random."; continue ;;
        esac
        break
    done
    echo "Replay mode: ${REPLAY_MODE_NAME}"

    while true; do
        read -rp "Include GEM hit-level branches in output? [no]: " answer
        case "${answer,,}" in
            ""|n|no|false)               REPLAY_GEM_HIT_FLAG="" ;;
            y|yes|true|gem_hit|-gem_hit) REPLAY_GEM_HIT_FLAG="-gem_hit" ;;
            *) echo "ERROR: enter yes/no or gem_hit."; continue ;;
        esac
        break
    done
    if [[ -n "${REPLAY_GEM_HIT_FLAG}" ]]; then
        echo "GEM hit-level branches will be included in the output"
    fi

    PRAD2_SOFT="$(prompt_default "Enter prad2evviewer directory" "${PRAD2_SOFT}")"
    if [[ "${PRAD2_BIN_WAS_SET}" -eq 0 ]]; then
        PRAD2_BIN="${PRAD2_SOFT}/build/bin"
    fi
    if [[ "${DEFAULT_CUTS_WAS_SET}" -eq 0 ]]; then
        DEFAULT_CUTS="${PRAD2_SOFT}/analysis/cuts/prad2_default.json"
    fi
    PRAD2_BIN="$(prompt_default "Enter executable directory" "${PRAD2_BIN}")"
    CACHE_BASE="$(prompt_default "Enter EVIO cache base directory" "${CACHE_BASE}")"
    OUTPUT_BASE="$(prompt_default "Enter output base directory" "${OUTPUT_BASE}")"
    REPLAY_CORES="$(prompt_default "Enter number of parallel jobs (-j)" "${REPLAY_CORES}")"
    REPLAY_ZERO_SUPPRESS="$(prompt_default "Enter GEM zero suppression (-z)" "${REPLAY_ZERO_SUPPRESS}")"
    REPLAY_MAX_FILES="$(prompt_default "Enter max number of files to process (-f)" "${REPLAY_MAX_FILES}")"
    REPLAY_MERGE_FILES="$(prompt_default "Enter replay merge group size (-m, 0 disables)" "${REPLAY_MERGE_FILES}")"

    echo "Enter cut JSON file for replay_filter (path to cuts.json, or 'default' to use ${DEFAULT_CUTS}):"
    read -rp "Cut JSON [default]: " answer
    if [[ -z "${answer}" || "${answer}" == "default" ]]; then
        CUT_JSON="${DEFAULT_CUTS}"
    else
        CUT_JSON="${answer}"
    fi
}

# ROOT setup and Slurm settings, then make every path absolute: the batch job
# runs from its output directory, not from where it was submitted.
prompt_job_settings() {
    ROOT_SETUP="$(prompt_default "Enter ROOT setup script, or 'none' to rely on submitted environment" "${ROOT_SETUP:-none}")"
    if [[ "${ROOT_SETUP}" == "none" ]]; then
        ROOT_SETUP=""
    fi

    SLURM_ACCOUNT="$(prompt_default "Enter Slurm account" "${SLURM_ACCOUNT}")"
    SLURM_PARTITION="$(prompt_default "Enter Slurm partition (production or priority)" "${SLURM_PARTITION}")"
    SLURM_TIME="$(prompt_default "Enter Slurm time limit" "${SLURM_TIME}")"
    SLURM_MEM_PER_CPU="$(prompt_default "Enter Slurm mem-per-cpu MB" "${SLURM_MEM_PER_CPU}")"

    PRAD2_SOFT="$(to_abs_path "${PRAD2_SOFT}" "${PWD}")"
    PRAD2_BIN="$(to_abs_path "${PRAD2_BIN}" "${PWD}")"
    CACHE_BASE="$(to_abs_path "${CACHE_BASE}" "${PWD}")"
    CUT_JSON="$(to_abs_path "${CUT_JSON}" "${PWD}")"
    if [[ -n "${ROOT_SETUP}" ]]; then
        ROOT_SETUP="$(to_abs_path "${ROOT_SETUP}" "${PWD}")"
    fi
    OUTPUT_BASE="$(to_abs_path "${OUTPUT_BASE}" "${PWD}")"
}

# --- jcache tape staging ---
# Ask for a notification email and request staging of each run's EVIO files.
# Returns 1 when no email is given or any request fails.
request_jcache() {
    local email run
    local failed=()
    echo "No EVIO files under ${CACHE_BASE} for run(s) $*; the data may still be on tape (MSS)."
    read -rp "Enter your email address for jcache notification: " email
    if [[ -z "${email}" ]]; then
        echo "ERROR: email cannot be empty."
        return 1
    fi

    echo ""
    for run in "$@"; do
        echo "jcache get ${MSS_BASE}/prad_${run}/* -e ${email}"
        if ! jcache get "${MSS_BASE}/prad_${run}"/* -e "${email}"; then
            echo "WARNING: jcache request failed for run ${run}."
            failed+=("${run}")
        fi
    done

    echo ""
    echo "Wait for staging, then re-run this script for these run(s):"
    echo "$*"
    if [[ "${#failed[@]}" -gt 0 ]]; then
        echo "jcache failed for these run(s); check manually:"
        echo "${failed[*]}"
        return 1
    fi
}

# --- Replay pipeline: replay_recon -> replay_filter -> live_charge -> quick_check ---
# Print and run one pipeline step; exit with its status if it fails.
run_step() {
    local label="$1"
    local rc
    shift
    echo "Command: $*"
    echo ""
    "$@" || { rc=$?; echo ""; echo "ERROR: ${label} exited with code ${rc}."; exit "${rc}"; }
    echo ""
    echo "${label^} finished."
}

# replay_filter -o for the given recon inputs.  Several inputs take the output
# directory; a single input takes a file name, spelled the way
# filtered_output_path in replay_filter.cpp names it in directory mode.
filter_output_path() {
    local out_dir="$1"
    shift
    if [[ "$#" -gt 1 ]]; then
        printf '%s\n' "${out_dir}"
        return
    fi
    local base
    base="$(basename "$1")"
    if [[ "${base}" =~ ^(.+)_recon_([^/]+)\.root$ ]] ||
       [[ "${base}" =~ ^(.+)\.evio\.([0-9]+)_recon\.root$ ]] ||
       [[ "${base}" =~ ^(.+)\.([0-9]+)_recon\.root$ ]]; then
        printf '%s\n' "${out_dir}/${BASH_REMATCH[1]}_filter_${BASH_REMATCH[2]}.root"
    else
        printf '%s\n' "${out_dir}/${base%.root}_filter.root"
    fi
}

# Run the whole pipeline for one run in OUT_DIR, which must exist.  Reads
# RUN_NUMBER RUN_DIR OUT_DIR PRAD2_BIN CUT_JSON REPLAY_CORES REPLAY_ZERO_SUPPRESS
# REPLAY_MAX_FILES REPLAY_MERGE_FILES REPLAY_MODE_FLAG REPLAY_GEM_HIT_FLAG.
# Exits on the first failing step.  A missing filter/live_charge/quick_check
# executable or cut JSON only skips that step; live_charge reads the filter
# outputs, so it is skipped with the filter.
run_replay_pipeline() {
    local replay_cmd="${PRAD2_BIN}/prad2ana_replay_recon"
    local filter_cmd="${PRAD2_BIN}/prad2ana_replay_filter"
    local live_charge_cmd="${PRAD2_BIN}/prad2ana_live_charge"
    local quick_check_cmd="${PRAD2_BIN}/prad2ana_quick_check"
    local report="${OUT_DIR}/prad_${RUN_NUMBER}_filter_report.json"
    local lc_json="${OUT_DIR}/prad_${RUN_NUMBER}_live_charge.json"
    local qc_output="${OUT_DIR}/prad_${RUN_NUMBER}_quick_check.root"
    local mode_args=() gem_hit_args=() recon_inputs=() lc_inputs=()
    local filter_out

    if [[ ! -x "${replay_cmd}" ]]; then
        echo "ERROR: executable not found: ${replay_cmd}"
        exit 1
    fi
    if [[ -n "${REPLAY_MODE_FLAG:-}" ]]; then
        read -ra mode_args <<< "${REPLAY_MODE_FLAG}"
    fi
    if [[ -n "${REPLAY_GEM_HIT_FLAG:-}" ]]; then
        read -ra gem_hit_args <<< "${REPLAY_GEM_HIT_FLAG}"
    fi

    echo ""
    echo "Starting replay..."
    run_step "replay" "${replay_cmd}" "${RUN_DIR}" -o "${OUT_DIR}" -j "${REPLAY_CORES}" \
        -z "${REPLAY_ZERO_SUPPRESS}" -f "${REPLAY_MAX_FILES}" -m "${REPLAY_MERGE_FILES}" \
        "${mode_args[@]}" "${gem_hit_args[@]}"

    # Merged prad_<RUN>_recon_<NNN>.root, or the per-split recon files with -m 0.
    mapfile -t recon_inputs < <(find "${OUT_DIR}" -maxdepth 1 -type f -name "prad_${RUN_NUMBER}_recon_*.root" | sort)
    if [[ "${#recon_inputs[@]}" -eq 0 ]]; then
        mapfile -t recon_inputs < <(find "${OUT_DIR}" -maxdepth 1 -type f -name "prad_${RUN_NUMBER}.*_recon.root" | sort)
    fi
    if [[ "${#recon_inputs[@]}" -eq 0 ]]; then
        echo "ERROR: no reconstructed ROOT files found in ${OUT_DIR}."
        exit 1
    fi
    echo "Downstream input ROOT file(s): ${#recon_inputs[@]}"

    if [[ ! -x "${filter_cmd}" ]]; then
        echo "WARNING: executable not found: ${filter_cmd}, skipping replay filter."
    elif [[ ! -f "${CUT_JSON:-}" ]]; then
        echo "WARNING: cut JSON not found: ${CUT_JSON:-}, skipping replay filter."
    else
        filter_out="$(filter_output_path "${OUT_DIR}" "${recon_inputs[@]}")"
        echo ""
        echo "Running replay filter..."
        run_step "replay filter" "${filter_cmd}" "${recon_inputs[@]}" -o "${filter_out}" \
            -c "${CUT_JSON}" -j "${report}" -t "${REPLAY_CORES}"
        echo "Filtered out  : ${filter_out}"
        # Only multi-input mode writes the run-level slow-control ROOT.
        if [[ "${#recon_inputs[@]}" -gt 1 ]]; then
            echo "Slow ROOT     : ${OUT_DIR}/prad_${RUN_NUMBER}_epics.root"
        fi
        echo "Filter report : ${report}"

        mapfile -t lc_inputs < <(find "${OUT_DIR}" -maxdepth 1 -type f -name "prad_${RUN_NUMBER}_filter*.root" | sort)
        if [[ "${#lc_inputs[@]}" -eq 0 ]]; then
            echo "ERROR: no filtered ROOT files found in ${OUT_DIR}; live_charge requires prad_${RUN_NUMBER}_filter*.root inputs."
            exit 1
        fi
    fi

    if [[ ! -x "${live_charge_cmd}" ]]; then
        echo "WARNING: executable not found: ${live_charge_cmd}, skipping live charge."
    elif [[ "${#lc_inputs[@]}" -eq 0 ]]; then
        echo "WARNING: replay filter did not run, skipping live charge."
    else
        echo ""
        echo "Running live charge calculation..."
        run_step "live charge" "${live_charge_cmd}" "${lc_inputs[@]}" -j "${lc_json}"
        echo "Live charge JSON: ${lc_json}"
    fi

    if [[ ! -x "${quick_check_cmd}" ]]; then
        echo "WARNING: executable not found: ${quick_check_cmd}, skipping quick check."
    else
        echo ""
        echo "Running quick check..."
        run_step "quick check" "${quick_check_cmd}" "${recon_inputs[@]}" -o "${qc_output}" -j "${REPLAY_CORES}"
        echo "Quick check ROOT: ${qc_output}"
    fi

    echo ""
    echo "Replay pipeline finished."
    echo "Recon inputs: ${#recon_inputs[@]}"
    echo "Output dir: ${OUT_DIR}"
}

# --- Slurm submission ---
# Write <OUTPUT_BASE>/prad_<RUN>/replay_recon_<RUN>.sbatch and submit it.  The
# batch script embeds the pipeline functions, so it stays self-contained and
# can be resubmitted by hand.  Every step is checked explicitly because callers
# run this where set -e is suspended.
submit_replay_job() {
    local run="$1"
    local out_dir="${OUTPUT_BASE}/prad_${run}"
    local sbatch_script="${out_dir}/replay_recon_${run}.sbatch"

    if [[ ! -x "${PRAD2_BIN}/prad2ana_replay_recon" ]]; then
        echo "ERROR: executable not found: ${PRAD2_BIN}/prad2ana_replay_recon"
        return 1
    fi
    mkdir -p "${out_dir}" || return 1

    cat > "${sbatch_script}" <<EOF || return 1
#!/bin/bash
#SBATCH --job-name=recon${run}
#SBATCH --account=${SLURM_ACCOUNT}
#SBATCH --partition=${SLURM_PARTITION}
#SBATCH --output=${out_dir}/replay_recon_${run}-%j.out
#SBATCH --error=${out_dir}/replay_recon_${run}-%j.err
#SBATCH --mail-user=$(whoami)@jlab.org
#SBATCH --time=${SLURM_TIME}
#SBATCH --mem-per-cpu=${SLURM_MEM_PER_CPU}
#SBATCH --cpus-per-task=${REPLAY_CORES}

set -euo pipefail

RUN_NUMBER=$(shell_quote "${run}")
RUN_DIR=$(shell_quote "${CACHE_BASE}/prad_${run}")
OUT_DIR=$(shell_quote "${out_dir}")
PRAD2_BIN=$(shell_quote "${PRAD2_BIN}")
CUT_JSON=$(shell_quote "${CUT_JSON}")
REPLAY_CORES=$(shell_quote "${REPLAY_CORES}")
REPLAY_ZERO_SUPPRESS=$(shell_quote "${REPLAY_ZERO_SUPPRESS}")
REPLAY_MAX_FILES=$(shell_quote "${REPLAY_MAX_FILES}")
REPLAY_MERGE_FILES=$(shell_quote "${REPLAY_MERGE_FILES}")
REPLAY_MODE_FLAG=$(shell_quote "${REPLAY_MODE_FLAG:-}")
REPLAY_GEM_HIT_FLAG=$(shell_quote "${REPLAY_GEM_HIT_FLAG:-}")
REPLAY_MODE_NAME=$(shell_quote "${REPLAY_MODE_NAME:-PRad2}")
ROOT_SETUP=$(shell_quote "${ROOT_SETUP}")

$(declare -f run_step filter_output_path run_replay_pipeline)

echo "Host: \$(hostname)"
echo "Work dir: \$(pwd)"
echo "Run: \${RUN_NUMBER}"
echo "Replay mode: \${REPLAY_MODE_NAME}"
echo "Input: \${RUN_DIR}"
echo "Output: \${OUT_DIR}"
date

if [[ -n "\${ROOT_SETUP}" ]]; then
    echo "Sourcing ROOT setup: \${ROOT_SETUP}"
    source "\${ROOT_SETUP}"
fi

if ! command -v root >/dev/null 2>&1; then
    echo "ERROR: root is not found in PATH. Set ROOT_SETUP or submit from a ROOT-ready environment."
    exit 1
fi
echo "ROOT is available: \$(root-config --version 2>/dev/null || root --version 2>&1 | head -1)"

if [[ ! -d "\${RUN_DIR}" ]]; then
    echo "ERROR: input directory not found: \${RUN_DIR}"
    exit 1
fi
FILE_COUNT=\$(find "\${RUN_DIR}" -maxdepth 1 -type f | wc -l)
if [[ "\${FILE_COUNT}" -eq 0 ]]; then
    echo "ERROR: no input files found in \${RUN_DIR}. Stage them with jcache first."
    exit 1
fi

mkdir -p "\${OUT_DIR}"
cd "\${OUT_DIR}"
echo "Job work dir: \$(pwd)"

run_replay_pipeline
date
EOF

    chmod +x "${sbatch_script}" || return 1
    echo ""
    echo "Prepared Slurm script: ${sbatch_script}"
    echo "Submitting job..."
    sbatch "${sbatch_script}"
}

# Submit one job per cached run, then request jcache staging for the rest.
# Returns 1 when any submission or staging request failed.
submit_replay_runs() {
    local run
    local status=0
    local cached=() missing=() failed=()

    echo ""
    echo "Runs requested: $*"
    echo "Checking cache under: ${CACHE_BASE}"
    for run in "$@"; do
        if run_is_cached "${run}"; then
            cached+=("${run}")
        else
            missing+=("${run}")
        fi
    done
    echo "Ready to submit: ${#cached[@]} run(s)"
    if [[ "${#cached[@]}" -gt 0 ]]; then
        echo "  ${cached[*]}"
    fi
    echo "Need jcache staging: ${#missing[@]} run(s)"
    if [[ "${#missing[@]}" -gt 0 ]]; then
        echo "  ${missing[*]}"
    fi

    echo ""
    for run in "${cached[@]}"; do
        echo "Submitting replay/recon job for run ${run}..."
        if ! submit_replay_job "${run}"; then
            echo "WARNING: submit failed for run ${run}."
            failed+=("${run}")
        fi
        echo ""
    done

    if [[ "${#missing[@]}" -gt 0 ]]; then
        request_jcache "${missing[@]}" || status=1
    elif [[ "${#failed[@]}" -eq 0 ]]; then
        echo "All requested runs were submitted; no jcache staging is needed."
    fi

    if [[ "${#failed[@]}" -gt 0 ]]; then
        echo ""
        echo "Submit failed for these cached run(s); check the messages above:"
        echo "${failed[*]}"
        status=1
    fi
    return "${status}"
}
