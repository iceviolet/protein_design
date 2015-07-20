# The Board of Trustees of the University of Pennsylvania. Copyright (C) 2013-2014. All Rights Reserved.
# Author: Chris Von Bargen

# Created:
# VERGIL_INCLUDE_DIRS - Directories to include to use Vergil
# VERGIL_LIBRARIES    - Default library to link against to use Vergil
# VERGIL_LINK_FLAGS   - Flags to be added to linker's options
# VERGIL_FOUND        - If false, don't try to use Vergil

# in linux if the env var VERGIL_DIR is not set, use /usr/local
set (VERGIL_DIR_TEST $ENV{VERGIL_DIR})
if (VERGIL_DIR_TEST)
  set (VERGIL_DIR $ENV{VERGIL_DIR} CACHE PATH "Path to VERGIL build directory")
else ()
  set (VERGIL_DIR /usr/local CACHE PATH "Path to VERGIL build directory")
endif()

# Find MPI, otherwise use the STUBs library
find_package (MPI)
if (MPI_CXX_FOUND)
  include_directories (${MPI_CXX_INCLUDE_PATH})
else (MPI_CXX_FOUND)
	find_library (MPI_CXX_LIBRARIES libmpi_stubs.a ${VERGIL_DIR}/lib)
	message ("-- Unable to find MPI, using Vergil MPI STUBS library...")
endif (MPI_CXX_FOUND)

# Find IPOPT
find_package (IPOPT)

# Set include directory and shared library object
set (VERGIL_INCLUDE_DIRS ${VERGIL_DIR}/include)
if(APPLE)
	find_library (VERGIL_LIBRARIES libvergil.dylib ${VERGIL_DIR}/lib)
ELSEIF(UNIX)
	find_library (VERGIL_LIBRARIES libvergil.so ${VERGIL_DIR}/lib)
ENDIF()

# Set linking flags
set (LINK_FLAGS "-O -g")
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set (OPENMP_FLAG "-fopenmp")
  set (LINK_FLAGS "${LINK_FLAGS} -mpc64")
endif () 
set (VERGIL_LINK_FLAGS "${OPENMP_FLAG} ${LINK_FLAGS} ${MPI_CXX_LINK_FLAGS}")

# Establish success/failure
if (VERGIL_LIBRARIES)
   set (VERGIL_FOUND TRUE)
	 message ("-- Found Vergil: ${VERGIL_DIR}")
	 message ("-- Vergil Build Flags: ${VERGIL_LINK_FLAGS}")
else ()
   set (VERGIL_FOUND FALSE)
   set (VERGIL_INCLUDE_DIRS "")
   set (VERGIL_LIBRARIES "")
   set (VERGIL_LINK_FLAGS "")
	 message("Vergil not found. Please install Vergil to use this target.")
endif()
