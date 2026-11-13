# cmake/WebDeps.cmake — FetchContent for the HTTP/WebSocket stack
#
# Creates three imported/interface targets:
#
#   nlohmann_json::nlohmann_json   — header-only JSON library
#   asio_standalone                — non-Boost Asio, header-only
#   websocketpp_lib                — header-only WebSocket server, links asio
#
# These are all header-only, so the cost at configure time is a shallow
# git clone and at build time it's whatever the consumers compile.

include_guard(GLOBAL)

include(FetchContent)
set(FETCHCONTENT_QUIET OFF CACHE BOOL "" FORCE)

FetchContent_Declare(json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG        v3.11.3
    GIT_SHALLOW    TRUE
)
FetchContent_Declare(websocketpp
    GIT_REPOSITORY https://github.com/zaphoyd/websocketpp.git
    GIT_TAG        0.8.2
    GIT_SHALLOW    TRUE
)
FetchContent_Declare(asio
    GIT_REPOSITORY https://github.com/chriskohlhoff/asio.git
    GIT_TAG        asio-1-30-2
    GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(json websocketpp asio)

# nlohmann's own install() is off as a sub-project, but the installed
# prad2 headers include <nlohmann/json.hpp> (ACLiC macros built against
# the install tree, see analysis/scripts/rootlogon.C), so ship it.
install(DIRECTORY ${json_SOURCE_DIR}/include/nlohmann
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Asio has no CMake target of its own; wrap the include dir as an
# interface library and set the standalone define here so consumers
# don't have to know.
if(NOT TARGET asio_standalone)
    add_library(asio_standalone INTERFACE)
    target_include_directories(asio_standalone INTERFACE
        ${asio_SOURCE_DIR}/asio/include)
    target_compile_definitions(asio_standalone INTERFACE ASIO_STANDALONE)
endif()

# WebSocketPP is header-only and transitively needs asio_standalone.
if(NOT TARGET websocketpp_lib)
    add_library(websocketpp_lib INTERFACE)
    target_include_directories(websocketpp_lib INTERFACE
        ${websocketpp_SOURCE_DIR})
    target_link_libraries(websocketpp_lib INTERFACE asio_standalone)
endif()
