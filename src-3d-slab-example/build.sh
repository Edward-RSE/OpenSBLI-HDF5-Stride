#!/bin/bash

source ~/modules/opensbli.lmod

export CC=gcc-14
export CXX=g++-14
export OMPI_CC=$CC
export OMPI_CXX=$CXX

OPS_INSTALL_DIR="$HOME/srsg-projects/opensbli/OPS-INSTALL-DEBUG"
HDF5_INSTALL_DIR="$HOME/srsg-projects/opensbli/HDF5"

cp ../include/io_strided3d.h io_strided.h && \
    rm -rf build && mkdir build && cd build || { echo "Failed to initialise build directory"; exit; }

cmake .. -DOPS_INSTALL_DIR=$OPS_INSTALL_DIR -DCMAKE_BUILD_TYPE=Debug -DHDF5_ROOT=$HDF5_INSTALL_DIR
cmake --build . -j
