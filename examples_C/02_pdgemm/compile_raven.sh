#!/usr/bin/bash

#set -v

module purge
module load intel/2023.1.0.x
module load impi/2021.9
module load mkl/2023.1


FileName=pdgemm
mpiicc -qmkl=sequential -o $FileName $FileName.c -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64
