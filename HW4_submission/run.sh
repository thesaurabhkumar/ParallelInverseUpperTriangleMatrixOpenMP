#!/bin/bash

module load intel/2017A
export OMP_NUM_THREADS=8
ulimit -S -s 1048576

icc -qopenmp -o par.exe inverse.cpp

echo "Running in Parallel ..."
bsub < grid.job




