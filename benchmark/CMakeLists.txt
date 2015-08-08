add_library(benchmark SHARED correspondPixels.cc match.cc)
target_link_libraries(benchmark util csa)