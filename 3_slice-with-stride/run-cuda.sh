#!/bin/bash

rm -f *.h5 \
    && ./build.sh \
    && ./build/slice-example_cuda \
    && python ../scripts/compare-slice-example.py opensbli_output_000002.h5 K75.h5 2 2 2 \
