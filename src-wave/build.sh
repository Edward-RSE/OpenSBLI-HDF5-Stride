#!/bin/bash

rm -r build && mkdir build && cd build || { echo "Failed to initialise build directory"; exit; }

OPS_INSTALL_DIR="$HOME/srsg-projects/opensbli/OPS-INSTALL"
HDF5_INSTALL_DIR="$HOME/srsg-projects/opensbli/HDF5"
OPTIMISATION="-DCFLAG=-ftree-vectorize -DCXXFLAG=-ftree-vectorize"

set -x
cmake .. -DOPS_INSTALL_DIR=$OPS_INSTALL_DIR -DCMAKE_BUILD_TYPE=Release -DHDF5_ROOT=$HDF5_INSTALL_DIR $OPTIMISATION
cmake --build . -j
