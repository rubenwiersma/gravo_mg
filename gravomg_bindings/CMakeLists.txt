cmake_minimum_required(VERSION 3.16)
project(gravomg_bindings)

add_subdirectory(deps/pybind11)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

add_subdirectory(../gravomg gravomg-build)

# Python Bindings
pybind11_add_module(gravomg_bindings
    src/cpp/core.cpp
)

find_package(OpenMP REQUIRED)
include_directories(gravomg_bindings ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp)
target_link_libraries(gravomg_bindings PRIVATE gravomg)