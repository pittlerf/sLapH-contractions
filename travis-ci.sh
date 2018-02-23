#!/bin/bash
# Copyright © 2017-2018 Martin Ueding <dev@martin-ueding.de>
# Licensed under the MIT/Expat license.

set -e
set -u
set -x

sourcedir="$(pwd)"
builddir=build-contraction

###############################################################################
#                              Install Packages                               #
###############################################################################

ubuntu_packages=(
    gcc-7 g++-7
    ccache
    cmake
    libhdf5-dev 
    hdf5-tools
    libeigen3-dev
    libboost-filesystem-dev libboost-system-dev libboost-program-options-dev
    libgtest-dev
)
sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
sudo apt-get update
sudo apt-get install -y "${ubuntu_packages[@]}"

###############################################################################
#                              Fix Eigen Library                              #
###############################################################################

sudo updatedb
locate FindEigen3.cmake

cd ..

mkdir cmake-module
cp $(locate FindEigen3.cmake) cmake-module

###############################################################################
#                              Build Google Test                              #
###############################################################################

# https://www.eriksmistad.no/getting-started-with-google-test-on-ubuntu/

mkdir "gtest"
pushd "gtest"
cmake /usr/src/gtest
make -j $(nproc)
sudo cp *.a /usr/lib/
popd

###############################################################################
#                          Build sLapH-contractions                           #
###############################################################################

rm -rf "$builddir"
mkdir -p "$builddir"
pushd "$builddir"

CXX=$(which g++-7)

# Compile gtest
# Modified from https://www.eriksmistad.no/getting-started-with-google-test-on-ubuntu/
pushd /usr/src/gtest
sudo cmake CMakeLists.txt -DCMAKE_CXX_COMPILER="$CXX"
sudo make -j $(nproc)
sudo cp *.a /usr/lib
popd

cmake "$sourcedir" -DCMAKE_MODULE_PATH=../cmake-module -DCMAKE_CXX_COMPILER="$CXX"
make -j $(nproc)

ctest --output-on-failure
