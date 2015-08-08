set(srcs Exception.cc String.cc Random.cc Timer.cc Matrix.cc kofn.cc)

add_library(util SHARED ${srcs})
target_link_libraries(util -lblas -lm)