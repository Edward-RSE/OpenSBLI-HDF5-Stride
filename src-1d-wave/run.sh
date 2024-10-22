#!/bin/bash

rm -f *.h5 \
    && ./build.sh \
    && ./build/wave_seq \
    && python ../scripts/compare-wave-halos.py opensbli_output.h5 opensbli_output-strided.h5 \
    && python ../scripts/compare-wave-halos.py opensbli_output.h5 opensbli_output-strided-single-precision.h5 \
    && rm io_strided.h