#!/usr/bin/env bash

args=(-DFABM_BASE=$FABMDIR -DFABM_ERSEM_BASE=$ERSEMDIR -DFABM_NIVA_BASE=$BROMDIR -DCMAKE_BUILD_TYPE=Release)

if [[ $# -ne 1 ]] ; then
    echo 'Please specify arguments north or north-2'
	exit 1
fi

if [[ ( "$1" = 'north' || "$1" = 'north-2' ) ]] ; then
	mkdir -p build
	rm -r build/* 2> /dev/null
	CWD=$(pwd)
	cp ../SPBM_data/$1/* build && cd build && cmake $CWD/src "${args[@]}"
else
	echo 'Wrong arguments'
	exit 1
fi
