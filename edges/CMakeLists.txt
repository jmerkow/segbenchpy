

add_library(edges SHARED edgesNms2d.cxx edgesNms3d.cxx)

if(GLOG_FOUND)
    target_link_libraries(edges ${GLOG_LIBRARY})
endif()

if(GLOG_FOUND)
    add_definitions("-DHAVE_GLOG")
    include_directories(${GLOG_INCLUDE_DIR})
endif()


add_executable(testNMS3D edgesNms3d.cxx)
target_compile_definitions(testNMS3D PUBLIC "TEXT_EXE")

if(GLOG_FOUND)
    target_link_libraries(testNMS3D ${GLOG_LIBRARY})
endif()