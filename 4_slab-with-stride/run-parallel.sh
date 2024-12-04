#!/bin/bash

rm -f *.h5 \
    && ./build.sh \
    && mpirun -n 4 ./build/slab-example_mpi \
    && python ../scripts/compare-slab-example.py opensbli_output_000004.h5 opensbli_output-slab_strided.h5 1 1 1 \
