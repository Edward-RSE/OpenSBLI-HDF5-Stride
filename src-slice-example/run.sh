#!/bin/bash

NP=4

rm -f *.h5 \
    && ./build.sh \
    && ./build/slice-example_seq \
    && python ../scripts/compare-slice-example.py opensbli_output_000001.h5 K75.h5 2 2 1 \

# &&  mpirun -np $NP ./build/slice-example_mpi \