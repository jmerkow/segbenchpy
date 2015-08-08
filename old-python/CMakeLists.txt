find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

find_package(PythonLibs)
include_directories(${PYTHON_INCLUDE_PATH})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})


foreach (_wrap csa)


include_directories(../${_wrap})


SET_SOURCE_FILES_PROPERTIES(${_wrap}.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(${_wrap}.i PROPERTIES SWIG_FLAGS "-I/usr/include -I${CMAKE_CURRENT_SOURCE_DIR} -includeall")
# SET_SOURCE_FILES_PROPERTIES(csa.i PROPERTIES SWIG_FLAGS "-includeall")
SWIG_ADD_MODULE(${_wrap} python ${_wrap}.i)
SWIG_LINK_LIBRARIES(${_wrap} ${PYTHON_LIBRARIES} ${_wrap})

endforeach()