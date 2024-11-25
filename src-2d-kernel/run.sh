#!/bin/bash

rm -f *.h5 \
    && ./build.sh \
    && ./build/kernel-example_seq \
    && python ../scripts/compare-2d-kernel-example.py opensbli_output.h5 opensbli_output-strided.h5 2 2
