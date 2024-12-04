#!/bin/bash

rm -f *.h5 \
    && ./build.sh \
    && ./build/slab-example_seq \
    && python ../scripts/compare-slab-example.py opensbli_output_000004.h5 opensbli_output-slab_strided.h5 1 1 1 \