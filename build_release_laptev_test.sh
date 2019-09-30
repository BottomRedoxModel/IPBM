#!/bin/bash

make_build () {
path=$1
srcpath=$HOME/SPBM/src
cd $path
cmake $srcpath -DFABM_BASE=$FABMDIR \
               -DFABM_ERSEM_BASE=$ERSEMDIR \
               -DFABM_NIVA_BASE=$BROMDIR \
               -DCMAKE_BUILD_TYPE=Release \
               #-DERSEM_USE_IRON=ON \
make
		   }

#change scenarios according to you case
scenarios="baseline-B baseline-O"

# change directory according to your case 
dir='/data_spbm_laptev/'

for scen in $scenarios
    do 	
        make_build $HOME$dir$scen
    done
