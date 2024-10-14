#!/bin/bash

rm -f opensbli_output-strided.h5 && \
    ./build.sh && \
    ./build/wave_seq && \
    python ../scripts/compare-wave.py ../data/wave-reference.h5 opensbli_output-strided.h5
