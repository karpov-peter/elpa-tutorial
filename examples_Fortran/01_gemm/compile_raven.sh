#!/usr/bin/bash

#set -v

module purge
module load intel/2023.1.0.x
module load impi/2021.9
module load mkl/2023.1


FileName=gemm_big_matrices
ifort -fopenmp -qmkl=parallel -I${MKLROOT}/include/intel64/lp64 -o $FileName $FileName.f90 -L${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a
