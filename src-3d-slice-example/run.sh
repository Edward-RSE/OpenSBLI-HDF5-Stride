#!/bin/bash

rm -f *.h5 \
    && ./build.sh \
    && ./build/slice-example_seq \
    && python ../scripts/compare-strided_output.py opensbli_output_000005.h5 K16.h5 2 2 2 \
    && rm io_strided.h
