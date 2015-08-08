find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

find_package(PythonLibs)
include_directories(${PYTHON_INCLUDE_PATH})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})



file(REMOVE ${CMAKE_CURRENT_BINARY_DIR}/segbench/__init__.py)
foreach (_wrap csa benchmark)


include_directories(../${_wrap})


SET_SOURCE_FILES_PROPERTIES(${_wrap}.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(${_wrap}.i PROPERTIES SWIG_FLAGS "-I/usr/include -I${CMAKE_CURRENT_SOURCE_DIR} -includeall")
# SET_SOURCE_FILES_PROPERTIES(csa.i PROPERTIES SWIG_FLAGS "-includeall")
SWIG_ADD_MODULE(${_wrap} python ${_wrap}.i)
SWIG_LINK_LIBRARIES(${_wrap} ${PYTHON_LIBRARIES} util ${_wrap} util)
file(APPEND ${CMAKE_CURRENT_SOURCE_DIR}/segbench/__init__.py "from ${_wrap} import *\n")

# if(UNIX OR APPLE)
    set(__linkname "${CMAKE_CURRENT_SOURCE_DIR}/segbench/_${_wrap}.so")
    set(__pynamegen "${CMAKE_CURRENT_BINARY_DIR}/${_wrap}.py")
    set(__pyname "${CMAKE_CURRENT_SOURCE_DIR}/segbench/${_wrap}.py")
    add_custom_command(TARGET _${_wrap} POST_BUILD
                       COMMAND ln -sf $<TARGET_LINKER_FILE:_${_wrap}> "${__linkname}"
                       COMMAND echo ln -sf $<TARGET_LINKER_FILE:_${_wrap}> "${__linkname}"
                       COMMENT "Creating symlink ${__linkname} -> ${CMAKE_CURRENT_BINARY_DIR}/segbench/_${_wrap}${CMAKE_SHARED_LIBRARY_SUFFIX}")
    add_custom_command(TARGET _${_wrap} POST_BUILD
                       COMMAND ln -sf ${__pynamegen} "${__pyname}"
                       COMMENT "Creating symlink ${__pyname} -> ${__pynamegen}")
# endif()

endforeach()