# test/CMakeLists.txt — small command-line dump / diagnostic tools.

# --- Installed tools (user-facing; install rule in the top-level CMakeLists.txt)
foreach(_t evio_dump ped_calc)
    add_prad_tool(${_t} LINKS nlohmann_json::nlohmann_json)
endforeach()

# gem_dump lives in gem/ (see gem/CMakeLists.txt).

# --- Dev-only tools (not installed) — sources under test/dev/
add_prad_tool(ts_dump   SRC dev/ts_dump.cpp   LINKS nlohmann_json::nlohmann_json)
add_prad_tool(livetime  SRC dev/livetime.cpp  LINKS nlohmann_json::nlohmann_json)
add_prad_tool(dsc_scan  SRC dev/dsc_scan.cpp  LINKS nlohmann_json::nlohmann_json)
add_prad_tool(hycal_energy_bias_test SRC dev/hycal_energy_bias_test.cpp
    LINKS prad2det nlohmann_json::nlohmann_json)
add_test(NAME prad2det.hycal_energy_bias COMMAND hycal_energy_bias_test)
add_prad_tool(gem_quality_test SRC dev/gem_quality_test.cpp
    LINKS prad2det nlohmann_json::nlohmann_json)
add_test(NAME prad2det.gem_quality COMMAND gem_quality_test)

# ET-only dev tools (compiled only when WITH_ET is on in prad2dec).
if(WITH_ET)
    add_prad_tool(evc_scan   SRC dev/evc_scan.cpp)
    add_prad_tool(et_feeder  SRC dev/et_feeder.cpp)
    add_prad_tool(evet_diff  SRC dev/evet_diff.cpp)
endif()
