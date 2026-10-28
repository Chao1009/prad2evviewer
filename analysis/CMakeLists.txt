# analysis/ — replay and physics analysis tools (requires ROOT, found by the
# top-level CMakeLists before this directory is added)

# Shared analysis library, linked by every tool and by the ACLiC scripts.
add_library(prad2ana STATIC
    src/Replay.cpp
    src/GainCorrCompute.cpp
    src/PhysicsTools.cpp
    src/MatchingTools.cpp
    src/SlowControl.cpp
    src/ToolUtils.cpp)

target_include_directories(prad2ana
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include/prad2ana>)

target_link_libraries(prad2ana PUBLIC
    prad2dec prad2det nlohmann_json::nlohmann_json
    ROOT::Core ROOT::Tree ROOT::TreePlayer ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Gpad ROOT::Spectrum ROOT::Physics)

# Tools: tools/<name>.cpp -> bin/prad2ana_<name>.  The CMake target name
# stays clean for internal referencing; the built/installed binary gets a
# `prad2ana_` prefix so multiple projects can coexist in the same bindir.
set(ANA_TOOLS
    replay_rawdata replay_gainCorr online_gain_monitor_base replay_recon
    replay_filter replay_add_hv live_charge
    cosmic_test quick_check GEM_matching det_align sim2replay sim2replay_hc
    gain_monitor gain_fitter refGain_produce
    hycal_rf_offset_calib nonlinearity calib_5by5 hycal_module_time_calib
    hycal_shower_profile time_recon_check physics_calib leakage_correction_check
    hycal_cluster_density hycal_Erecon_check hycal_energy_bowl_shape)
foreach(_t IN LISTS ANA_TOOLS)
    add_prad_tool(${_t} SRC tools/${_t}.cpp LINKS prad2ana)
    set_target_properties(${_t} PROPERTIES OUTPUT_NAME prad2ana_${_t})
endforeach()

install(TARGETS ${ANA_TOOLS} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# Python tools — staged into <build>/bin/ alongside the C++ executables and
# installed under the same prad2ana_ name so users do not need to know which
# language a tool happens to be implemented in.
set(_REPLAY_REPORT_VIEWER_SRC
    ${CMAKE_SOURCE_DIR}/scripts/replay_report_viewer.py)
add_custom_target(replay_report_viewer ALL
    COMMAND ${CMAKE_COMMAND} -E copy
            ${_REPLAY_REPORT_VIEWER_SRC}
            ${CMAKE_BINARY_DIR}/bin/prad2ana_replay_report_viewer
    DEPENDS ${_REPLAY_REPORT_VIEWER_SRC}
    COMMENT "Staging prad2ana_replay_report_viewer")
install(PROGRAMS ${_REPLAY_REPORT_VIEWER_SRC}
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    RENAME prad2ana_replay_report_viewer)

# Install the static analysis library so installed-mode ACLiC scripts can
# link against analysis::* symbols (no shared lib needed since tools also
# link statically).
install(TARGETS prad2ana ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Install analysis/include headers so installed-tree ACLiC scripts can
# include "PhysicsTools.h" / "MatchingTools.h" / "Replay.h" / "ConfigSetup.h"
# the same way they do from the build tree.
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/prad2ana
        FILES_MATCHING PATTERN "*.h")

# Install the ACLiC macro scripts and the script_helpers.h they include for
# the install-tree workflow described in rootlogon.C:
# `root -l <prefix>/share/prad2evviewer/analysis/scripts/rootlogon.C`.
install(DIRECTORY scripts/ DESTINATION ${CMAKE_INSTALL_DATADIR}/prad2evviewer/analysis/scripts
        FILES_MATCHING PATTERN "*.C" PATTERN "*.h")

# Install pyscripts/ — both the standalone analysis tools and the shared
# _common.py helper imported by scripts/gem_hycal_match_viewer.py.
install(DIRECTORY pyscripts/ DESTINATION ${CMAKE_INSTALL_DATADIR}/prad2evviewer/analysis/pyscripts
        FILES_MATCHING PATTERN "*.py")
