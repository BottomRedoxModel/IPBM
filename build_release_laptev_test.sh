#!/bin/bash

#mkdir -p build
#rm -r build/*
 
make_build () {
path=$1
srcpath=$HOME/SPBM/src
cd $path
#cp ../data_spbm_laptev/baseline-B/* build &&
#cd build &&
cmake $srcpath -DFABM_BASE=$FABMDIR \
               -DFABM_ERSEM_BASE=$ERSEMDIR \
               -DFABM_NIVA_BASE=$BROMDIR \
               -DCMAKE_BUILD_TYPE=Release \
               #-DERSEM_USE_IRON=ON \
make
		   }

scenarios="baseline-B baseline-O"
for scen in $scenarios
    do 	
        make_build $HOME'/data_spbm_laptev/'$scen
    done
