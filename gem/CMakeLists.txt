# gem/CMakeLists.txt — GEM diagnostic tools and helper scripts (see README.md).

# --- C++ binary --------------------------------------------------------------
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
