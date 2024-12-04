#!/bin/bash

rm -f *.h5 \
    && ./build.sh \
    && mpirun -n 4 ./build/slab-example_mpi \
