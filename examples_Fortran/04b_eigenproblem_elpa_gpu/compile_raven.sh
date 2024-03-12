#!/usr/bin/bash

set -v

module purge
module load gcc/11 cuda/11.4 impi/2021.11
module load elpa/mpi/standard/gpu/2023.11.001

echo $ELPA_HOME
echo $ELPA_LIBS

ELPA_INCLUDE_DIR=$ELPA_HOME/include/elpa-2023.11.001
#ELPA_INCLUDE_DIR=$ELPA_HOME/include/elpa_openmp-2023.11.001

CUDA_LIBS="-L$CUDA_HOME/lib64 -lcudart -lcuda -lm -Wl,-rpath=$CUDA_HOME/lib64"

FileName="eigenproblem_elpa_gpu"

mpigfortran -I$ELPA_INCLUDE_DIR/modules $FileName.f90 -o $FileName $ELPA_LIBS $CUDA_LIBS
