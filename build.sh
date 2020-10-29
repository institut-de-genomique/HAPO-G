#!/bin/bash
set -e
set -o pipefail

mkdir build
cd build

cmake ../src/
make
cd ..

mkdir bin
ln -s ../build/hapog bin/hapog

echo "HAPoG was successfully built!"
