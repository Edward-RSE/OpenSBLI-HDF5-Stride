#!/bin/bash

#   && ./build/kernel-example_seq \
#   && mpirun -n 4 ./build/kernel-example_mpi \

rm -f *.h5 \
    && ./build.sh \
    && mpirun -n 4 ./build/kernel-example_mpi \
    && python ../scripts/compare-2d-kernel-example.py opensbli_output.h5 opensbli_output-strided.h5 2 2 --halo_size 5