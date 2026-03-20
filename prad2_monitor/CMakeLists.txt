cmake_minimum_required(VERSION 3.14)

find_package(Qt5 REQUIRED COMPONENTS Widgets WebEngineWidgets)

add_executable(prad2_monitor main.cpp)
target_link_libraries(prad2_monitor PRIVATE
    Qt5::Widgets
    Qt5::WebEngineWidgets
)
set_target_properties(prad2_monitor PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    AUTOMOC ON
)
