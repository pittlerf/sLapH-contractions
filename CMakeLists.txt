# Copyright © 2017-2018 Martin Ueding <dev@martin-ueding.de>

# Enforce a version that is actually available on QBiG.
cmake_minimum_required(VERSION 2.8)

project(sLapH-contractions CXX C)

include_directories("include" "include/GlobalData/")

###############################################################################
#                          Git Revision Information                           #
###############################################################################

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(GetGitRevisionDescription)

get_git_head_revision(git_refspec git_sha1)
git_local_changes(git_changes)

message(STATUS "git sha1: ${git_sha1}")
message(STATUS "git refspec: ${git_refspec}")
message(STATUS "git state: ${git_changes}")

execute_process(COMMAND
  hostname -f
  OUTPUT_VARIABLE hostname
  ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

set(git_user $ENV{USER})
set(git_host ${hostname})

message(STATUS "user: ${git_user}")
message(STATUS "host: ${git_host}")

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/include/git.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/include/git.h"
  @ONLY)
list(APPEND SOURCES "${CMAKE_CURRENT_BINARY_DIR}/include/git.h")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/include")

###############################################################################
#                             External Libraries                              #
###############################################################################

# HDF5
if(NOT DEFINED HDF5_INCLUDE_DIRS OR NOT DEFINED HDF5_LIBRARIES)
  find_package(HDF5 REQUIRED COMPONENTS C CXX)
endif()
include_directories(${HDF5_INCLUDE_DIRS})

# Boost
if(NOT DEFINED Boost_INCLUDE_DIRS OR NOT DEFINED Boost_LIBRARIES)
  find_package(Boost REQUIRED COMPONENTS filesystem system program_options)
endif()
include_directories(${Boost_INCLUDE_DIRS})

# Eigen
if(NOT DEFINED EIGEN3_INCLUDE_DIR)
  find_package(Eigen3 REQUIRED)
endif()
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

# Lime
if(NOT DEFINED LIME_INCLUDE_DIRS)
  message(FATAL_ERROR "Must specify -DLIME_INCLUDE_DIRS on the command line to CMake.")
endif()
include_directories(SYSTEM ${LIME_INCLUDE_DIRS})

# Lime
if(NOT DEFINED LIME_LIBRARIES)
  message(FATAL_ERROR "Must specify -DLIME_LIBRARIES on the command line to CMake.")
endif()

# OpenMP
find_package(OpenMP)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

###############################################################################
#                                 Build Flags                                 #
###############################################################################

# We want warnings, but not all of them.
add_definitions(-Wall -pedantic --std=c++11)

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-fdiagnostics-color has_diagnostics_color)
if(has_diagnostics_color)
  if (CMAKE_GENERATOR STREQUAL "Ninja")
    add_definitions("-fdiagnostics-color=always")
  else()
    add_definitions("-fdiagnostics-color=auto")
  endif()
endif()

# Default to “Release” build type.
message(STATUS "Build Type: '${CMAKE_BUILD_TYPE}'")
if(CMAKE_BUILD_TYPE STREQUAL "")
  message(STATUS "No CMAKE_BUILD_TYPE, assuming release and enabling all optimizations.")
  set(CMAKE_BUILD_TYPE Release)
endif()

# Tune for the local architecture because this code will be used where it is
# compiled, in most cases.
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")

###############################################################################
#                                 Executables                                 #
###############################################################################

### contract

add_executable(contract

  main/contract.cpp
  modules/Correlators.cpp
  modules/Diagram.cpp
  modules/dilution-iterator.cpp
  modules/EigenVector.cpp
  modules/Gamma.cpp
  modules/GaugeField.cpp
  modules/GlobalData/global_data.cpp
  modules/GlobalData/global_data_build_IO_names.cpp
  modules/GlobalData/global_data_input_handling.cpp
  modules/GlobalData/global_data_input_handling_utils.cpp
  modules/GlobalData/init_lookup_tables.cpp
  modules/h5-wrapper.cpp
  modules/io_utils.cpp
  modules/OperatorsForMesons.cpp
  modules/Perambulator.cpp
  modules/QuarkLineBlock2.cpp
  modules/RandomVector.cpp
  modules/ranlxs.cpp
  modules/StopWatch.cpp

  )

target_link_libraries(contract 
  ${Boost_LIBRARIES}
  ${HDF5_LIBRARIES}
  ${LIME_LIBRARIES}
  )

install(TARGETS contract DESTINATION bin)

### test-iterator

add_executable(test-dilution-iterator
  main/test-iterator.cpp
  modules/dilution-iterator.cpp
  )

### test-derivative
add_executable(test-derivative
  main/test-derivative.cpp
  modules/io_utils.cpp
  modules/GaugeField.cpp
  modules/EigenVector.cpp
  )

target_link_libraries(test-derivative
  ${LIME_LIBRARIES}
  )

###############################################################################
#                               Subdirectories                                #
###############################################################################

add_subdirectory(tests/gtest)

###############################################################################
#                                    Tests                                    #
###############################################################################

enable_testing()

add_test(NAME unit-dilution-iterator
  COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/tests/dilution-iterator/run-test" "${CMAKE_CURRENT_SOURCE_DIR}")

add_test(NAME integration-L4
  COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/tests/integration-L4/run-integration-test" "${CMAKE_CURRENT_BINARY_DIR}"
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/tests/integration-L4")
