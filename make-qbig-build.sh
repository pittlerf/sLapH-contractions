#!/bin/bash
# Copyright © 2017 Martin Ueding <dev@martin-ueding.de>
# Licensed under The MIT/Expat License

# Sets up a debug and a release build outside of the source code with the paths
# that are needed on QBiG. The debug build will automatically set the `-g` flag
# for debugging, the release version will be performance and has `-O2 -DNDEBUG`
# set.

set -e
set -u
set -x

for buildtype in release debug; do
    builddir=../sLapH-contractions-$buildtype
    rm -rf "$builddir"
    mkdir "$builddir"
    pushd "$builddir"

    cmake \
        ../cntr.v0.1 \
        -DEIGEN3_INCLUDE_DIRS='/hiskp2/werner/libraries/eigen-3.2.7' \
        -DHDF5_INCLUDE_DIRS=/hiskp2/knippsch/hdf5-1.8.17/include \
        -DHDF5_LIBRARIES='-L/hiskp2/knippsch/hdf5-1.8.17/lib -lhdf5_cpp -lhdf5 -lsz -lz' \
        -DCMAKE_CXX_COMPILER=g++-4.7 \
        -DCMAKE_BUILD_TYPE=$buildtype

    make -j $(nproc)
    popd
done
