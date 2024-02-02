#!/usr/bin/bash

#set -v

module purge
module load intel/2023.1.0.x
module load impi/2021.9
module load mkl/2023.1


FileName=gemm_big_matrices
icc -DUSE_MKL -fopenmp -qmkl=parallel -o $FileName $FileName.c
