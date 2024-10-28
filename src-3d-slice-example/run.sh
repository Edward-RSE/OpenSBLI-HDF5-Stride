#!/bin/bash

rm -f *.h5 \
    && ./build.sh \
    && ./build/slice-example_seq\
    && rm io_strided.h
