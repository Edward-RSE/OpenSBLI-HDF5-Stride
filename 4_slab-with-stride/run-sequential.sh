#!/bin/bash

rm -f *.h5 \
    && ./build.sh \
    && ./build/slab-example_seq \
    && python ../scripts/compare-slab-example.py opensbli_output.h5 rho_B0-slab.h5 1 1 1