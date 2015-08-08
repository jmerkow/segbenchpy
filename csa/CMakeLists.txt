cmake_minimum_required(VERSION 2.8)

project(csa++)

set(srcs csa.cc csaAssign.cc)
set(hdrs csa.hh csa_types.h csa_defs.h)

add_library(csa SHARED ${srcs})
add_executable(test test.cc)
target_link_libraries(test csa)