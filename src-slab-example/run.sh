#!/bin/bash

NP=2

rm -f *.h5 \
    && ./build.sh \
    && ./build/slab-example_seq
    # && python ../scripts/compare-slice-example.py opensbli_output_000001.h5 K75.h5 2 2 1 \
