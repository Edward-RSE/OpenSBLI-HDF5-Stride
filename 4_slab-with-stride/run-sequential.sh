#!/bin/bash

rm -f *.h5 \
    && ./build.sh \
    && time ./build/slab-example_seq \
    && rm  -f ../scripts/fig/*.png \
    && python ../scripts/compare-slab-example.py opensbli_output_000001.h5 opensbli_output-slab_strided.h5 1 1 1 \