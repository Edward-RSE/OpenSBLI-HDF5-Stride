#!/bin/bash

NP=4

rm -f *.h5 \
    && ./build.sh \
    &&  mpirun -np $NP ./build/slice-example_mpi \
    && python ../scripts/compare-strided_output.py opensbli_output_000001.h5 K75.h5 2 2 1 \
