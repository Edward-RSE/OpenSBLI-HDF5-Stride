#!/bin/bash

./build.sh && ./build/wave_seq && python ../scripts/compare-wave.py ../data/wave-reference.h5 opensbli_output.h5
