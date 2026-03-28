# Source this file to set up the PRad2 environment (csh/tcsh).
#   source <prefix>/bin/setup.csh

set PRAD2_DIR = `dirname $0`
set PRAD2_DIR = `cd "$PRAD2_DIR/.." && pwd`

setenv PATH "${PRAD2_DIR}/bin${PATH:+:$PATH}"
if ($?LD_LIBRARY_PATH) then
    setenv LD_LIBRARY_PATH "${PRAD2_DIR}/lib:${LD_LIBRARY_PATH}"
else
    setenv LD_LIBRARY_PATH "${PRAD2_DIR}/lib"
endif
setenv PRAD2_DATABASE_DIR "${PRAD2_DIR}/share/prad2evviewer/database"
setenv PRAD2_RESOURCE_DIR "${PRAD2_DIR}/share/prad2evviewer/resources"
