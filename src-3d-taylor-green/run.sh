#!/bin/bash

rm -f *.h5 \
    && ./build.sh \
    && ./build/taylorgreen_seq \
    && python ../scripts/compare-taylor-green.py opensbli_output.h5 opensbli_output-strided.h5 2 2 2 --halo_size 5 \
    && python ../scripts/compare-taylor-green.py opensbli_output.h5 opensbli_output-strided-single-precision.h5 2 2 2 --halo_size 5 \
    && rm io_strided.h
