#!/bin/bash

mkdir -p build
rm -r build/*
CWD=`pwd`
cp ../SPBM_data/laptev-2/* build &&
cd build && cmake $CWD/src -DFABM_BASE=$FABMDIR \
                           -DFABM_ERSEM_BASE=$ERSEMDIR \
                           -DFABM_NIVA_BASE=$BROMDIR \
                           -DCMAKE_BUILD_TYPE=Debug \
                           #-DERSEM_USE_IRON=ON \
