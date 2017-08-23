#!/bin/bash
# Copyright © 2017 Martin Ueding <dev@martin-ueding.de>
# Licensed under the MIT/Expat license.

set -e
set -u
set -x

ubuntu_packages=(
    gcc-6 g++-6
    ccache
    cmake
    libhdf5-dev 
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

# Download the `local` directory if desired.
wget http://qphix.martin-ueding.de/local.tar.gz
tar -xzf local.tar.gz


mkdir cmake-module
cp $(locate FindEigen3.cmake) cmake-module


rm -rf "$builddir"
mkdir -p "$builddir"
cd "$builddir"

# Get lime flags
cmake "$sourcedir" -DCMAKE_MODULE_PATH=../cmake-module \
    -DLIME_INCLUDE_DIRS='../local/include/' \
    -DLIME_LIBRARIES='-L../local/lib/ -llime'
make -j $(nproc) VERBOSE=1
