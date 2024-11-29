#!/bin/bash

rm -f *.h5 \
    && ./build.sh \
    && mpirun -n 4 ./build/slice-example_cuda \
    && python ../scripts/compare-slice-example.py opensbli_output_000001.h5 K75.h5 2 2 2 \
