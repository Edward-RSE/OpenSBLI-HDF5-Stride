#!/bin/bash

export CC=gcc
export CXX=g++
export OMPI_CC=$CC
export OMPI_CXX=$CXX

OPS_INSTALL_DIR="$HOME/codes/OpenSBLI_JAXA/OPS-INSTALL"
HDF5_INSTALL_DIR="$HOME/codes/OpenSBLI_JAXA/HDF5"

rm -f *.h5
rm -rf build && mkdir build && cd build || { echo "Failed to initialise build directory"; exit; }

cmake .. -DOPS_INSTALL_DIR=$OPS_INSTALL_DIR -DCMAKE_BUILD_TYPE=Release -DHDF5_ROOT=$HDF5_INSTALL_DIR -DGPU_ARCH=80
cmake --build . -j
