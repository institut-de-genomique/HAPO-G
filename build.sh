#!/bin/bash
set -e
set -o pipefail


help()
{
   echo "Hapo-G building script"
   echo
   echo "Syntax: bash build.sh -l PATH_TO_HTSLIB"
   echo "options:"
   echo "l     Path to htslib. As an example, if samtools is installed in /home/user/samtools, htslib is probably in /home/user/samtools/htslib"
   echo "h     Print this help."
   echo
}


while getopts "hl:" flag
do
    case "${flag}" in
        l) 
            htslib_root=${OPTARG}
            ;;
        h) 
            help
            exit
            ;;
        \?) 
            echo "Invalid option: -$OPTARG" >&2
            ;;
    esac
done

mkdir hapog_build
cd hapog_build

export HTSLIB_ROOT=${htslib_root}
cmake ../src/
make
cd ..

mkdir bin
ln -s ../hapog_build/hapog bin/hapog

echo "HAPoG was successfully built!"
