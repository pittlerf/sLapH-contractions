#!/bin/bash
# Copyright © 2017 Martin Ueding <dev@martin-ueding.de>

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

cd ..
mkdir build
cd build

cmake "$sourcedir"
