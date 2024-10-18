#!/bin/bash

rm -f *.h5 \
    && ./build.sh \
    && ./build/taylorgreen_seq \
    && python ../scripts/compare-taylor-green-halos.py opensbli_output.h5 opensbli_output-strided.h5 \
    && rm io_strided.h
