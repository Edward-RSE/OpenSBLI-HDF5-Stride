#!/bin/bash

OPS_INSTALL_DIR="$HOME/srsg-projects/opensbli/OPS-INSTALL"
HDF5_INSTALL_DIR="$HOME/srsg-projects/opensbli/HDF5"

rm -rf build && mkdir build && cd build || { echo "Failed to initialise build directory"; exit; }

cmake .. -DOPS_INSTALL_DIR=$OPS_INSTALL_DIR -DCMAKE_BUILD_TYPE=Debug -DHDF5_ROOT=$HDF5_INSTALL_DIR $OPTIMISATION
cmake --build . -j
