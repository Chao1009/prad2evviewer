# Source this file to set up the PRad2 environment (bash/zsh).
#   source <prefix>/bin/prad2_setup.sh

PRAD2_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

export PATH="${PRAD2_DIR}/bin${PATH:+:$PATH}"

# Prepend both lib64 and lib — evio + prad2py land in lib64 on RHEL-family
# systems (GNUInstallDirs default), our static libs land in lib.  Listing
# both keeps the install layout tolerant of either convention.
export LD_LIBRARY_PATH="${PRAD2_DIR}/lib64:${PRAD2_DIR}/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PYTHONPATH="${PRAD2_DIR}/lib64/prad2py:${PRAD2_DIR}/lib/prad2py${PYTHONPATH:+:$PYTHONPATH}"

export PRAD2_DATABASE_DIR="${PRAD2_DIR}/share/prad2evviewer/database"
export PRAD2_RESOURCE_DIR="${PRAD2_DIR}/share/prad2evviewer/resources"
