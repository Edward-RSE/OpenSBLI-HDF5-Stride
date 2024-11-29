#!/bin/bash

export CC=/opt/homebrew/bin/gcc-14
export CXX=/opt/homebrew/bin/g++-14
export OMPI_CC=$CC
export OMPI_CXX=$CXX

OPENSBLI_INSTALL_DIR="$HOME/srsg-projects/opensbli"
OPS_INSTALL_DIR="$OPENSBLI_INSTALL_DIR/OPS-INSTALL"
HDF5_INSTALL_DIR="$OPENSBLI_INSTALL_DIR/HDF5"

rm -rf build && mkdir build && cd build || { echo "Failed to initialise build directory"; exit; }
cmake .. -DOPS_INSTALL_DIR=$OPS_INSTALL_DIR -DCMAKE_BUILD_TYPE=Debug -DHDF5_ROOT=$HDF5_INSTALL_DIR
cmake --build . -j
