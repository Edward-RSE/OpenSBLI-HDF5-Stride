#!/bin/bash

rm -f *.h5 \
    && ./build.sh \
    && ./build/slice-example_cuda \
    && python ../scripts/compare-strided_output.py opensbli_output_000500.h5 K150.h5 2 2 2 \
    && rm io_strided.h
