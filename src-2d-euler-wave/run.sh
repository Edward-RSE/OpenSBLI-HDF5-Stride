#!/bin/bash

rm -f *.h5 \
    && ./build.sh \
    && ./build/eulerwave_seq \
    && rm io_strided.h
