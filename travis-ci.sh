#!/bin/bash
# Copyright © 2017 Martin Ueding <dev@martin-ueding.de>
# Licensed under the MIT/Expat license.

set -e
set -u
set -x

ubuntu_packages=(
    gcc-7 g++-7
    ccache
    cmake
    libhdf5-dev 
    hdf5-tools
    libeigen3-dev
    libboost-filesystem-dev libboost-system-dev libboost-program-options-dev
)
sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
sudo apt-get update
sudo apt-get install -y "${ubuntu_packages[@]}"

sourcedir="$(pwd)"
builddir=build-contraction

sudo updatedb
locate FindEigen3.cmake

cd ..


mkdir cmake-module
cp $(locate FindEigen3.cmake) cmake-module


rm -rf "$builddir"
mkdir -p "$builddir"
cd "$builddir"

CXX=$(which g++-7)

cmake "$sourcedir" -DCMAKE_MODULE_PATH=../cmake-module -DCMAKE_CXX_COMPILER="$CXX"
make -j $(nproc)

ctest --output-on-failure
