#!/usr/bin/bash

set -v

module purge
module load intel/2023.1.0.x
module load impi/2021.9
module load mkl/2023.1
module load elpa/mpi/standard/2023.11.001

echo $ELPA_HOME
echo $ELPA_LIBS

ELPA_INCLUDE_DIR=$ELPA_HOME/include/elpa-2023.11.001
#ELPA_INCLUDE_DIR=$ELPA_HOME/include/elpa_openmp-2023.11.001

FileName="eigenproblem_elpa_autotuning"

mpiicc -std=c11 -I$ELPA_INCLUDE_DIR $FileName.c -o $FileName $ELPA_LIBS
