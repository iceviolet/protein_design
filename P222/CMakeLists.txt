# The Board of Trustees of the University of Pennsylvania. Copyright (C) 2013-2014. All Rights Reserved.

##
## Example CMake file for your own application
##

cmake_minimum_required (VERSION 2.8)

set (PROJ_NAME P222)

##
## Dependencies
##

set (MAINFOLDER ${PROJECT_SOURCE_DIR})
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/tools/modules/")
find_package (VERGIL)

##
## Build project
##

if (VERGIL_FOUND)

  file (GLOB_RECURSE PROJ_SRC_FILES src/*.cpp src/*.h)
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set (PROJ_COMPILE_FLAG "-Wall -W -Wno-uninitialized -fomit-frame-pointer -fno-exceptions -g -fstrict-aliasing -fopenmp")
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set (PROJ_COMPILE_FLAG "-Wall -W -Wno-uninitialized -fomit-frame-pointer -fno-exceptions -g -fstrict-aliasing")
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set (PROJ_COMPILE_FLAG "-Wall -W -Wno-uninitialized -fno-exceptions -g -fstrict-aliasing -openmp")
  endif ()  

  include_directories (${VERGIL_INCLUDE_DIRS})
  foreach(SOURCE_FILE ${PROJ_SRC_FILES}) 
  	get_filename_component(EXEC_NAME ${SOURCE_FILE} NAME_WE)
  	add_executable (${EXEC_NAME} ${SOURCE_FILE})
  	target_link_libraries (${EXEC_NAME} ${VERGIL_LIBRARIES})
  	target_link_libraries (${EXEC_NAME} ${MPI_CXX_LIBRARIES})
  	target_link_libraries (${EXEC_NAME} ${IPOPT_LIBRARIES})
  	set_target_properties(${EXEC_NAME} PROPERTIES COMPILE_FLAGS ${PROJ_COMPILE_FLAG} LINK_FLAGS ${VERGIL_LINK_FLAGS})

  	add_custom_target ( "run_${EXEC_NAME}" DEPENDS ${EXEC_NAME} COMMENT "Executing target ${EXEC_NAME}..." VERBATIM SOURCES ${SOURCE_FILE})
 	add_custom_command (TARGET "run_${EXEC_NAME}" POST_BUILD COMMAND "${CMAKE_CURRENT_BINARY_DIR}/${EXEC_NAME}")
  endforeach( SOURCE_FILE ${PROJ_SRC_FILES} )	

endif ()

##
## Clean all command
## http://stackoverflow.com/questions/9680420/looking-for-a-cmake-clean-command-to-clear-up-cmake-output
##

add_custom_target(clean-all COMMAND ${CMAKE_BUILD_TOOL} clean COMMAND ${CMAKE_COMMAND} -P ${MAINFOLDER}/tools/clean-all.cmake 
  COMMENT "Removing all compiled binaries, libraries, and associated CMake files in build directory")
