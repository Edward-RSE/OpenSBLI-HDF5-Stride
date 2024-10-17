#!/bin/bash

rm -f *.h5 \
    && ./build.sh \
    && ./build/taylorgreen_seq \
    && python ../scripts/compare-taylor-green.py opensbli_output_000001.h5 opensbli_output-strided.h5
    # && python ../scripts/compare-taylor-green.py ../data/taylor-green-reference.h5 opensbli_output-strided.h5
