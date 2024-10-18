#!/bin/bash

rm -f opensbli_output-strided.h5 \
    && ./build.sh \
    && ./build/wave_seq \
    && python ../scripts/compare-wave-halos.py ../data/wave-reference.h5 opensbli_output-strided.h5 \
    && rm io_strided.h