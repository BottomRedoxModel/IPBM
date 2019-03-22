#!/bin/bash

mkdir -p build
rm -r build/*
CWD=`pwd`
cp ../data_spbm_laptev/baseline-0.001/* build &&
cd build && cmake $CWD/src -DFABM_BASE=$FABMDIR \
                           -DFABM_NIVA_BASE=$BROMDIR \
                           -DCMAKE_BUILD_TYPE=Release \
                           #-DERSEM_USE_IRON=ON \
