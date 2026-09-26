# python/CMakeLists.txt — prad2py (pybind11 bindings for prad2dec + prad2det)
#
# Built when BUILD_PYTHON=ON (root CMakeLists).  Fetches pybind11 via
# FetchContent so no system install is required.  The resulting
# extension module lands in ${CMAKE_BINARY_DIR}/python/ — add that
# directory to PYTHONPATH to import it, e.g.:
#
#   export PYTHONPATH=$PWD/build/python:$PYTHONPATH
#   python -c "import prad2py; print(prad2py.__version__)"

include(FetchContent)

FetchContent_Declare(pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG        v2.13.6
    GIT_SHALLOW    TRUE
)
# Tell pybind11 to skip its own Python search — use the interpreter resolved
# by the top-level probe.
set(PYBIND11_FINDPYTHON ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(pybind11)

pybind11_add_module(prad2py MODULE
    prad2py.cpp
    bind_dec.cpp
    bind_det.cpp
)

target_link_libraries(prad2py PRIVATE
    prad2dec
    prad2det
    nlohmann_json::nlohmann_json
)

set_target_properties(prad2py PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python"
)

# prad2_setup.sh and the generated bin/ wrappers put
# <prefix>/<libdir>/prad2py on PYTHONPATH.
install(TARGETS prad2py
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}/prad2py"
)

message(STATUS "Python bindings: output -> ${CMAKE_BINARY_DIR}/python/")
