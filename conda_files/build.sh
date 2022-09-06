#!/bin/bash

export C_INCLUDE_PATH=${PREFIX}/include
export LIBRARY_PATH=${PREFIX}/lib
echo $PREFIX

pip install .
bash build.sh -l $LIBRARY_PATH
cp -r build/hapog ${PREFIX}/bin/hapog_bin
