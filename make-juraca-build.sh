#!/bin/bash
# Copyright © 2017 Martin Ueding <dev@martin-ueding.de>
# Modified by Christian Jost
# Licensed under The MIT/Expat License

# Sets up a debug and a release build outside of the source code with the paths
# that are needed on QBiG. The debug build will automatically set the `-g` flag
# for debugging, the release version will be performance and has `-O2 -DNDEBUG`
# set.

set -e
set -u
set -x

module load CMake
# Load standard compiler
module load defaults
module load HDF5
module load Boost
module load Eigen

sourcedir=$(pwd)

for buildtype in release; do
    builddir=../sLapH-contractions-$buildtype
    rm -rf "$builddir"
    mkdir "$builddir"
    pushd "$builddir"

    cmake \
        "$sourcedir" \
         -DEIGEN3_INCLUDE_DIR='${EBROOTEIGEN}/include/' \
        -DCMAKE_C_COMPILER=mpicc \
        -DCMAKE_CXX_COMPILER=mpicxx \
        -DCMAKE_BUILD_TYPE=$buildtype

    make -j $(nproc)
    popd
done
