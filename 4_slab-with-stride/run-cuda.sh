#!/bin/bash

rm -f *.h5 \
    && ./build.sh \
    && ./build/slab-example_cuda \
