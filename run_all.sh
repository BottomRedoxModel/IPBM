#!/bin/sh
for i in $(find $HOME"/data_spbm_laptev/" -maxdepth 2 -name "SPBM") 
    do
	    dir=$(dirname $i)
	    cd $dir 
	    echo "pwd $(pwd)"
	$i
    done
