#!/usr/bin/bash

set -v

module purge
module load gcc/14
module load openmpi/5.0
module load mkl/2025.1
module load elpa/mpi/standard/2025.01.001

echo $ELPA_HOME
echo $ELPA_LIBS

#ELPA_INCLUDE_DIR=$ELPA_HOME/include/elpa-2023.11.001
#ELPA_INCLUDE_DIR=$ELPA_HOME/include/elpa_openmp-2023.11.001
ELPA_INCLUDE_DIR=$ELPA_HOME/include/elpa-2025.01.001

FileName="eigenproblem_elpa"

#mpiicc -std=c11 -I$ELPA_INCLUDE_DIR $FileName.c -o $FileName $ELPA_LIBS
mpicc -std=c11 -I$ELPA_INCLUDE_DIR $FileName.c -o $FileName $ELPA_LIBS -lm
