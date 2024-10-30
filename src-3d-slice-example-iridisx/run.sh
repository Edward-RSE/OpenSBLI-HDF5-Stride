#!/bin/bash

rm -f *.h5 \
    && ./build.sh \
    && ./build/slice-example_cuda \
    && rm io_strided.h
