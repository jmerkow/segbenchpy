cmake_minimum_required(VERSION 2.8)
project(seg-bench)
set(PACKAGE_VERSION "0.1")

include_directories(csa util Benchmark)
add_subdirectory(csa)
add_subdirectory(util)
add_subdirectory(benchmark)

add_subdirectory(python)
