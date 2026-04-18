# gem/CMakeLists.txt — GEM diagnostic tools and helper scripts.
#
# Contents:
#   gem_dump.cpp           — C++ command-line inspector (installed)
#   gem_event_viewer.py    — PyQt6 event-by-event viewer (installed via wrapper)
#   gem_cluster_view.py    — static cluster plotter for evdump JSON (wrapper)
#   gem_layout.py          — strip-layout visualizer (wrapper)
#   gem_strip_map.py       — thin wrapper over prad2py.det.map_strip
#   gem_view.py            — rendering library used by the above
#   check_strip_map.py     — cross-check against reference implementations
#   README.md              — overview + GEM detector facts

# --- C++ binary --------------------------------------------------------------
# add_prad_tool() links prad2dec, adds DATABASE_DIR compile-def, sends the
# binary to <build>/bin/.  Defaults SRC to <name>.cpp relative to this
# CMakeLists — here that's gem/gem_dump.cpp.
add_prad_tool(gem_dump LINKS prad2det nlohmann_json::nlohmann_json)

install(TARGETS gem_dump RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# --- Python scripts + README -------------------------------------------------
# Installed under <prefix>/share/prad2evviewer/gem/.  Sibling imports work
# because Python puts the script's own directory on sys.path when it runs
# ``python <path>/gem_event_viewer.py``.  The bin wrapper for each tool
# points at this directory.
install(DIRECTORY .
    DESTINATION ${CMAKE_INSTALL_DATADIR}/prad2evviewer/gem
    FILES_MATCHING
        PATTERN "*.py"
        PATTERN "README.md"
        PATTERN "__pycache__" EXCLUDE
        PATTERN "CMakeLists.txt" EXCLUDE
)
