#!/bin/bash
# replay_recon.sh — run the single-run PRad replay pipeline on JLab ifarm
#
# Usage: bash /path/to/prad2evviewer/scripts/shell/replay_recon.sh
#   Run it from the checkout; it sources replay_common.sh next to it.
#   Prompts for run number and settings, then:
#     1. Checks the EVIO files in <cache_base>/prad_<RUN>/ (jcache if missing)
#     2. Runs prad2ana_replay_recon; by default, recon ROOTs are merged in
#        groups of 62 into prad_<RUN>_recon_<NNN>.root
#     3. Runs prad2ana_replay_filter on all recon ROOTs and writes the
#        corresponding *_filter.root files directly under the run output dir
#     4. Runs prad2ana_live_charge over all filtered ROOTs together
#     5. Runs prad2ana_quick_check over all recon ROOTs
#
# The one "parallel jobs" prompt is reused for replay_recon -j,
# replay_filter -t, and quick_check -j. All products are written under
# <output_base>/prad_<RUN>/; no filter subdirectory is created.

source "$(dirname "${BASH_SOURCE[0]}")/replay_common.sh" ||
    { echo "ERROR: replay_common.sh not found next to this script"; exit 1; }

# Check ROOT environment
if ! command -v root &>/dev/null; then
    echo "ERROR: 'root' is not found in PATH. Please set up the ROOT environment first."
    echo "       e.g.: source /path/to/root/bin/thisroot.sh"
    exit 1
fi
if ! root -l -q &>/dev/null; then
    echo "ERROR: 'root -l' failed. The ROOT installation may be incomplete or misconfigured."
    exit 1
fi
echo "ROOT is available: $(root-config --version 2>/dev/null || root --version 2>&1 | head -1)"

prompt_run_number
prompt_replay_settings

RUN_DIR="${CACHE_BASE}/prad_${RUN_NUMBER}"
OUT_DIR="${OUTPUT_BASE}/prad_${RUN_NUMBER}"

echo ""
echo "Checking input directory: ${RUN_DIR}"
if ! run_is_cached "${RUN_NUMBER}"; then
    request_jcache "${RUN_NUMBER}"
    exit
fi

mkdir -p "${OUT_DIR}"
echo "Output directory: ${OUT_DIR}"
run_replay_pipeline
