#!/bin/bash

mkdir -p build
rm -r build/*
CWD=`pwd`
cp data/* build &&
cd build && cmake $CWD/code -DFABM_BASE=$FABMDIR -DCMAKE_BUILD_TYPE=Debug
