# test/CMakeLists.txt — small command-line dump / diagnostic tools.
#
# Each tool links prad2dec (via add_prad_tool).  LINKS adds extra private
# dependencies; SRC overrides the default source path (test/<name>.cpp).
#
# Layout:
#   test/        — installed, user-facing diagnostic tools.
#   test/dev/    — dev-only test scaffolding (smoketests, ET plumbing,
#                  livetime calc, timestamp dump).  Built but never installed.

# ---------------------------------------------------------------------------
# Installed tools (user-facing)
# ---------------------------------------------------------------------------

# Tools that need nlohmann/json (for --daq-config parsing).
foreach(_t evio_dump ped_calc)
    add_prad_tool(${_t} LINKS nlohmann_json::nlohmann_json)
endforeach()

# gem_dump lives in gem/ (see gem/CMakeLists.txt).

# ---------------------------------------------------------------------------
# Dev-only tools (not installed) — sources under test/dev/
# ---------------------------------------------------------------------------
add_prad_tool(ts_dump   SRC dev/ts_dump.cpp   LINKS nlohmann_json::nlohmann_json)
add_prad_tool(livetime  SRC dev/livetime.cpp  LINKS nlohmann_json::nlohmann_json)

# ET-only dev tools (compiled only when WITH_ET is on in prad2dec).
#   evc_scan  — smoketest reader (evio file or ET ring).
#   et_feeder — feeds an evio file into an ET ring.
#   evet_diff — reads both an evio file and an ET ring in parallel and
#               diffs their raw buffers; pairs with et_feeder.
if(WITH_ET)
    add_prad_tool(evc_scan   SRC dev/evc_scan.cpp)
    add_prad_tool(et_feeder  SRC dev/et_feeder.cpp)
    add_prad_tool(evet_diff  SRC dev/evet_diff.cpp)
endif()
