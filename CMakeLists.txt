# Copyright © 2017 Martin Ueding <dev@martin-ueding.de>
# Licensed under the MIT/Expat license.

# Enforce a version that is actually available on QBiG.
cmake_minimum_required(VERSION 2.8)

project(sLapH-contractions CXX C)

include_directories("include" "include/GlobalData/")

if(NOT DEFINED HDF5_INCLUDE_DIRS OR NOT DEFINED HDF5_LIBRARIES)
    find_package(HDF5 REQUIRED COMPONENTS C CXX)
endif()
include_directories(${HDF5_INCLUDE_DIRS})

if(NOT DEFINED Boost_INCLUDE_DIRS OR NOT DEFINED Boost_LIBRARIES)
    find_package(Boost REQUIRED COMPONENTS filesystem system program_options)
endif()
include_directories(${Boost_INCLUDE_DIRS})

if(NOT DEFINED EIGEN3_INCLUDE_DIRS)
    find_package(Eigen3 REQUIRED)
endif()
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIRS})

find_package(OpenMP)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

add_definitions(-Wall -pedantic --std=c++11)
add_definitions(-Wno-unused-variable -Wno-unused-local-typedefs -Wno-sign-compare)

add_executable(contract
    modules/RandomVector.cpp
    modules/Correlators/Correlators.cpp
    modules/EigenVector.cpp
    modules/Quarklines_one_t.cpp
    modules/ranlxs.cpp
    modules/Quarklines.cpp
    modules/Perambulator.cpp
    modules/GlobalData/init_lookup_tables.cpp
    modules/GlobalData/global_data_input_handling_utils.cpp
    modules/GlobalData/global_data_input_handling.cpp
    modules/GlobalData/global_data_build_IO_names.cpp
    modules/GlobalData/global_data.cpp
    modules/OperatorsForMesons.cpp
    main/contract.cpp
    )

target_link_libraries(contract 
    ${Boost_LIBRARIES}
    ${HDF5_LIBRARIES}
    )

install(TARGETS contract DESTINATION bin)
