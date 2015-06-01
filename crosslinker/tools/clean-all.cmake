##
## Clean all routine
## http://stackoverflow.com/questions/9680420/looking-for-a-cmake-clean-command-to-clear-up-cmake-output
##

set(cmake_generated ${CMAKE_BINARY_DIR}/CMakeCache.txt
                    ${CMAKE_BINARY_DIR}/cmake_install.cmake  
                    ${CMAKE_BINARY_DIR}/Makefile
                    ${CMAKE_BINARY_DIR}/CMakeFiles
                    ${CMAKE_BINARY_DIR}/benchmark
                    ${CMAKE_BINARY_DIR}/docs
                    ${CMAKE_BINARY_DIR}/src
                    ${CMAKE_BINARY_DIR}/test
                    ${CMAKE_BINARY_DIR}/tools
                    ${CMAKE_BINARY_DIR}/working
                    ${CMAKE_BINARY_DIR}/vergil_config.h
)

foreach(file ${cmake_generated})
  if (EXISTS ${file})
     file(REMOVE_RECURSE ${file})
  endif()
endforeach(file)
