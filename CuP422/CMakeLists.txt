# The Board of Trustees of the University of Pennsylvania. Copyright (C) 2013-2014. All Rights Reserved.

##
## Example CMake file for your own application
##

cmake_minimum_required (VERSION 2.8)

set (PROJ_NAME example)

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
  set (PROJ_COMPILE_FLAG "-Wall -W -Wno-uninitialized -fomit-frame-pointer -fno-exceptions -g -fstrict-aliasing")

  include_directories (${VERGIL_INCLUDE_DIRS})
  add_executable (${PROJ_NAME} ${PROJ_SRC_FILES})
  target_link_libraries (${PROJ_NAME} ${VERGIL_LIBRARIES})
  target_link_libraries (${PROJ_NAME} ${MPI_CXX_LIBRARIES})
  target_link_libraries (${PROJ_NAME} ${IPOPT_LIBRARIES})
  set_target_properties(${PROJ_NAME} PROPERTIES COMPILE_FLAGS ${PROJ_COMPILE_FLAG} LINK_FLAGS ${VERGIL_LINK_FLAGS})

  add_custom_target ( "run_${PROJ_NAME}" DEPENDS ${PROJ_NAME} COMMENT "Executing target ${PROJ_NAME}..." VERBATIM SOURCES ${PROJ_SRC_FILES})
  add_custom_command (TARGET "run_${PROJ_NAME}" POST_BUILD COMMAND "${CMAKE_CURRENT_BINARY_DIR}/${PROJ_NAME}")

endif ()

##
## Clean all command
## http://stackoverflow.com/questions/9680420/looking-for-a-cmake-clean-command-to-clear-up-cmake-output
##

add_custom_target(clean-all COMMAND ${CMAKE_BUILD_TOOL} clean COMMAND ${CMAKE_COMMAND} -P ${MAINFOLDER}/tools/clean-all.cmake 
  COMMENT "Removing all compiled binaries, libraries, and associated CMake files in build directory")
