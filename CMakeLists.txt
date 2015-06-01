# The Board of Trustees of the University of Pennsylvania. Copyright (C) 2013-2014. All Rights Reserved.

##
## Includes and flags
##

include_directories (${COMMON_INCLUDES})
include_directories (${MPI_CXX_INCLUDE_PATH})

set (WORKING_BUILD_FLAGS "${OPENMP_FLAG} ${CXXDEBUG_FLAG} ${OPTIMIZE_FLAG} ${WARNING_FLAG} ${MPI_CXX_COMPILE_FLAGS}")
set (WORKING_LINK_FLAGS "${OPENMP_FLAG} ${LINK_FLAGS} ${MPI_CXX_LINK_FLAGS}")

##
## Create Working Executable
##

function(create_working_target WORKING_EXEC_NAME SRC_FILE)  

  add_executable (${WORKING_EXEC_NAME} ${SRC_FILE})
  target_link_libraries (${WORKING_EXEC_NAME} ${PROJECT_NAME})
  target_link_libraries (${WORKING_EXEC_NAME} ${MPI_CXX_LIBRARIES})
  set_target_properties(${WORKING_EXEC_NAME} PROPERTIES COMPILE_FLAGS ${WORKING_BUILD_FLAGS} LINK_FLAGS ${WORKING_LINK_FLAGS})

  add_custom_target ( "run_${WORKING_EXEC_NAME}" DEPENDS ${WORKING_EXEC_NAME} COMMENT "Executing target ${WORKING_EXEC_NAME}..." VERBATIM SOURCES ${SRC_FILE})
  add_custom_command (TARGET "run_${WORKING_EXEC_NAME}" POST_BUILD COMMAND "${CMAKE_CURRENT_BINARY_DIR}/${WORKING_EXEC_NAME}")
  
endfunction()

##
## Add working targets
##
  
create_working_target(example src/example.cpp)
create_working_target(NSC_calculation src/NSC_example.cpp)
create_working_target(terbium src/terbium/Terbium.cpp)
create_working_target(helix src/helix.cpp)
create_working_target(smo src/SMO.cpp)
create_working_target(training_set src/training_set_2.cpp)
create_working_target(CuP422 src/CuP422.cpp)
create_working_target(Crosslinker src/crosslinker.cpp)