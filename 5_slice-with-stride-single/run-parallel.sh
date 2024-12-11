#!/bin/bash

rm -f *.h5 \
    && ./build.sh \
    && mpirun -n 4 ./build/slice-example_mpi \
    && python ../scripts/compare-slice-example.py opensbli_output_000002.h5 K75.h5 2 2 2 \
