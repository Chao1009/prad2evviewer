# cmake/PradHelpers.cmake — shared helpers for the prad2evviewer build
#
# Provides three tiny utilities used throughout the build:
#
#   prad_option(NAME DEFAULT DESCRIPTION)
#       Drop-in replacement for option() that also records whether the user
#       explicitly passed -D${NAME}=<val> on the command line.  The result
#       is available afterwards as ``${NAME}_EXPLICIT``.
#
#   prad_optional_probe(NAME FOUND <flag-var> [HINT <msg>])
#       Inspect the result of a prior find_package(… QUIET) call.  If the
#       dependency was found, do nothing.  If not:
#         - ${NAME}_EXPLICIT = TRUE  → FATAL_ERROR with HINT + escape route.
#         - otherwise                → WARNING + set ${NAME} OFF in cache.
#
#   add_prad_tool(NAME [SRC <path>] [LINKS <libs…>])
#       Build a command-line tool executable (test/, gem/, analysis/).
#       Links prad2dec by default, adds DATABASE_DIR compile definition,
#       and drops the binary in ${CMAKE_BINARY_DIR}/bin.

include_guard(GLOBAL)


# prad_option — option() + remember explicit cmdline usage
macro(prad_option _name _default)
    # Remaining args are the human-readable description.
    set(_desc "${ARGN}")

    # Was the variable set BEFORE we create the option default?  The only
    # way that happens is via -D on the command line (or an earlier
    # CMakeLists, which we don't have).
    if(DEFINED ${_name})
        set(${_name}_EXPLICIT TRUE)
    else()
        set(${_name}_EXPLICIT FALSE)
    endif()
    option(${_name} "${_desc}" ${_default})
endmacro()


# prad_optional_probe — fail-or-disable based on find_package result
macro(prad_optional_probe _name)
    # HINT is multiValue so callers can supply several "…\n" fragments
    # without worrying about semicolon handling.
    cmake_parse_arguments(_PROBE "" "FOUND" "HINT" ${ARGN})

    if(NOT _PROBE_FOUND)
        message(FATAL_ERROR "prad_optional_probe: missing FOUND <var>")
    endif()

    # Collapse HINT fragments into one string (joined by empty — each
    # fragment should carry its own trailing \n if the author wants one).
    set(_PROBE_HINT_JOINED "")
    foreach(_frag IN LISTS _PROBE_HINT)
        string(APPEND _PROBE_HINT_JOINED "${_frag}")
    endforeach()

    if(NOT ${${_PROBE_FOUND}})
        if(${_name}_EXPLICIT)
            message(FATAL_ERROR
                "${_name}=ON but the required dependency was not found.\n"
                "${_PROBE_HINT_JOINED}\n"
                "Re-configure with -D${_name}=OFF to skip.")
        else()
            message(WARNING
                "Required dependency not found — disabling ${_name}. "
                "Core targets still build.\n${_PROBE_HINT_JOINED}\n"
                "Pass -D${_name}=OFF explicitly to silence this message.")
            set(${_name} OFF CACHE BOOL "" FORCE)
        endif()
    endif()
endmacro()


# add_prad_tool — one-liner executable for command-line tools
function(add_prad_tool _name)
    set(_oneValue SRC)
    set(_multiValue LINKS)
    cmake_parse_arguments(_T "" "${_oneValue}" "${_multiValue}" ${ARGN})

    # Default source: <name>.cpp relative to the caller's CMakeLists.  From
    # test/CMakeLists.txt this resolves to test/<name>.cpp automatically.
    if(NOT _T_SRC)
        set(_T_SRC ${_name}.cpp)
    endif()

    add_executable(${_name} ${_T_SRC})
    target_link_libraries(${_name} PRIVATE prad2dec ${_T_LINKS})
    target_compile_definitions(${_name} PRIVATE DATABASE_DIR="${DATABASE_DIR}")
    set_target_properties(${_name} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
endfunction()
