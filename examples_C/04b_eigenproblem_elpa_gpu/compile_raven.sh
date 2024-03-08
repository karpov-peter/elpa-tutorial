#!/usr/bin/bash

set -v

module purge
module load gcc/11 cuda/11.4 impi/2021.11 mkl/2023.1
module load elpa/mpi/standard/gpu/2023.11.001

echo $ELPA_HOME
echo $ELPA_LIBS

ELPA_INCLUDE_DIR=$ELPA_HOME/include/elpa-2023.11.001

CUDA_LIBS="-L$CUDA_HOME/lib64 -lcudart -lcuda -lm -Wl,-rpath=$CUDA_HOME/lib64"

FileName="eigenproblem_elpa_gpu"

mpigcc -I$ELPA_INCLUDE_DIR -I$MKL_HOME/include -I$CUDA_HOME/include $FileName.c -o $FileName $ELPA_LIBS $CUDA_LIBS
