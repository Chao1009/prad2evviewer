#!/bin/bash
# submit_replay_recon.sh
# Submit one PRad replay pipeline job to JLab ifarm Slurm.
#
# The generated job runs:
#   1. prad2ana_replay_recon, with optional grouped hadd merging via -m
#   2. prad2ana_replay_filter over all recon ROOTs
#   3. prad2ana_live_charge over all filtered ROOTs together
#   4. prad2ana_quick_check over all recon ROOTs
#
# All outputs are written directly under <output_base>/prad_<RUN>/.
# The requested CPU count is shared by replay_recon -j, replay_filter -t,
# quick_check -j, and Slurm --cpus-per-task.

set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/replay_common.sh" ||
    { echo "ERROR: replay_common.sh not found next to this script"; exit 1; }

echo "Submit one PRad replay/recon Slurm job"
echo ""

prompt_run_number
prompt_replay_settings
prompt_job_settings
submit_replay_runs "${RUN_NUMBER}"
