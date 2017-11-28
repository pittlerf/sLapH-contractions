# Copyright © 2017 Martin Ueding <dev@martin-ueding.de>

# Enforce a version that is actually available on QBiG.
cmake_minimum_required(VERSION 2.8)

project(sLapH-contractions CXX C)

include_directories("include" "include/GlobalData/")

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

# OpenMP
find_package(OpenMP)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

###############################################################################
#                                 Build Flags                                 #
###############################################################################

# We want warnings, but not all of them.
add_definitions(-Wall -pedantic --std=c++11)
add_definitions(-Wno-unused-variable -Wno-unused-local-typedefs)

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-fdiagnostics-color has_diagnostics_color)
if(has_diagnostics_color)
    add_definitions("-fdiagnostics-color=always")
endif()


# Default to “Release” build type.
message(STATUS "Build Type: '${CMAKE_BUILD_TYPE}'")
if(CMAKE_BUILD_TYPE STREQUAL "")
    message(STATUS "No CMAKE_BUILD_TYPE, assuming release and enabling all optimizations.")
    set(CMAKE_BUILD_TYPE Release)
endif()

# Tune for the local architecture because this code will be used where it is
# compiled, in most cases.
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")

###############################################################################
#                                 Executables                                 #
###############################################################################

### contract

add_executable(contract

    main/contract.cpp
    modules/Correlators.cpp
    modules/dilution-iterator.cpp
    modules/EigenVector.cpp
    modules/GlobalData/global_data.cpp
    modules/GlobalData/global_data_build_IO_names.cpp
    modules/GlobalData/global_data_input_handling.cpp
    modules/GlobalData/global_data_input_handling_utils.cpp
    modules/GlobalData/init_lookup_tables.cpp
    modules/OperatorsForMesons.cpp
    modules/Perambulator.cpp
    modules/QuarkLineBlock.cpp
    modules/Quarklines.cpp
    modules/Quarklines_one_t.cpp
    modules/RandomVector.cpp
    modules/ranlxs.cpp
    modules/StopWatch.cpp

    )

target_link_libraries(contract 
    ${Boost_LIBRARIES}
    ${HDF5_LIBRARIES}
    )

if(cmake_build_type_lower STREQUAL "debug")
    target_link_libraries(contract)
endif()

install(TARGETS contract DESTINATION bin)

### test-iterator

add_executable(test-dilution-iterator
    main/test-iterator.cpp
    modules/dilution-iterator.cpp
    )

###############################################################################
#                                    Tests                                    #
###############################################################################

enable_testing()

add_test(NAME unit-dilution-iterator
    COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/tests/dilution-iterator/run-test" "${CMAKE_CURRENT_SOURCE_DIR}")

add_test(NAME integration-L4
    COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/tests/integration-L4/run-integration-test" "${CMAKE_CURRENT_BINARY_DIR}"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/tests/integration-L4")
