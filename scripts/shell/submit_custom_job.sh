#!/bin/bash
# submit_custom_job.sh
# Prompt for a full shell command, then submit it as a Slurm batch job.
#
# Usage:
#   ./submit_custom_job.sh
#
# The script asks for:
#   - Slurm account/partition/time/memory/cpus
#   - working directory
#   - optional ROOT setup script (or blank)
#   - a full command line, for example:
#         source /path/to/setup.sh && my_exe input.root -o output.root
#
# The command is saved to a helper shell script and then run inside the job.

set -euo pipefail

read -rp "Enter Slurm account [hallb]: " SLURM_ACCOUNT
SLURM_ACCOUNT="${SLURM_ACCOUNT:-hallb}"

read -rp "Enter Slurm partition [production]: " SLURM_PARTITION
SLURM_PARTITION="${SLURM_PARTITION:-production}"

read -rp "Enter Slurm time [12:00:00]: " SLURM_TIME
SLURM_TIME="${SLURM_TIME:-12:00:00}"

read -rp "Enter memory per CPU in MB [1500]: " SLURM_MEM_PER_CPU
SLURM_MEM_PER_CPU="${SLURM_MEM_PER_CPU:-1500}"

read -rp "Enter cpus-per-task [1]: " SLURM_CPUS_PER_TASK
SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-1}"

read -rp "Enter job name [custom_job]: " JOB_NAME
JOB_NAME="${JOB_NAME:-custom_job}"

read -rp "Enter working directory [$(pwd)]: " WORK_DIR
WORK_DIR="${WORK_DIR:-$(pwd)}"
if [[ ! -d "${WORK_DIR}" ]]; then
    echo "ERROR: working directory does not exist: ${WORK_DIR}"
    exit 1
fi
WORK_DIR="$(cd "${WORK_DIR}" && pwd)"

read -rp "Enter output directory [${WORK_DIR}/slurm_out]: " OUTPUT_DIR
OUTPUT_DIR="${OUTPUT_DIR:-${WORK_DIR}/slurm_out}"
mkdir -p "${OUTPUT_DIR}"

read -rp "Enter ROOT setup script, or 'none' to rely on the current environment [none]: " ROOT_SETUP
ROOT_SETUP="${ROOT_SETUP:-none}"
if [[ "${ROOT_SETUP}" == "none" ]]; then
    ROOT_SETUP=""
elif [[ ! -f "${ROOT_SETUP}" ]]; then
    echo "ERROR: ROOT setup file not found: ${ROOT_SETUP}"
    exit 1
fi

read -rp "Enter the full command to run: " COMMAND
if [[ -z "${COMMAND}" ]]; then
    echo "ERROR: command cannot be empty."
    exit 1
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
JOB_SCRIPT_PATH="${OUTPUT_DIR}/${JOB_NAME}_${STAMP}.sbatch"
CMD_SCRIPT_PATH="${OUTPUT_DIR}/${JOB_NAME}_${STAMP}.command.sh"

printf '%s\n' "${COMMAND}" > "${CMD_SCRIPT_PATH}"
chmod +x "${CMD_SCRIPT_PATH}"

cat > "${JOB_SCRIPT_PATH}" <<EOF
#!/bin/bash
#SBATCH --job-name=${JOB_NAME}
#SBATCH --account=${SLURM_ACCOUNT}
#SBATCH --partition=${SLURM_PARTITION}
#SBATCH --output=${OUTPUT_DIR}/${JOB_NAME}-%j.out
#SBATCH --error=${OUTPUT_DIR}/${JOB_NAME}-%j.err
#SBATCH --time=${SLURM_TIME}
#SBATCH --mem-per-cpu=${SLURM_MEM_PER_CPU}
#SBATCH --cpus-per-task=${SLURM_CPUS_PER_TASK}

set -euo pipefail

echo "Host: \\$(hostname)"
echo "Working directory: ${WORK_DIR}"
echo "Command script: ${CMD_SCRIPT_PATH}"

cd "${WORK_DIR}"

if [[ -n "${ROOT_SETUP}" ]]; then
    echo "Sourcing ROOT setup: ${ROOT_SETUP}"
    source "${ROOT_SETUP}"
fi

bash "${CMD_SCRIPT_PATH}"
EOF

chmod +x "${JOB_SCRIPT_PATH}"

echo ""
echo "Prepared Slurm job script: ${JOB_SCRIPT_PATH}"
echo "Saved command script: ${CMD_SCRIPT_PATH}"
echo "Submitting job..."
sbatch "${JOB_SCRIPT_PATH}"
