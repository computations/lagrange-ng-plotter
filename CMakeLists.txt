cmake_minimum_required(VERSION 3.13)
project(lagrange-ng-plotter)
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

add_subdirectory(libs)
add_subdirectory(src)
