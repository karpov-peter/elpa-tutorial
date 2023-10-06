#!/usr/bin/bash

set -v

module purge

# module load LUMI/21.12  partition/C
# module load ELPA/2021.11.001-cpeGNU-21.12
# echo ELPA_INCLUDE_DIR=$ELPA_INCLUDE_DIR

module load PrgEnv-gnu
module load cray-libsci
ELPA_HOME=/scratch/project_465000539/ELPA/ELPA_2023.11_CPU
ELPA_LIBS="-L$ELPA_HOME/lib -lelpa -Wl,-rpath=$ELPA_HOME/lib"
ELPA_INCLUDE_DIR=$ELPA_HOME/include/elpa-2023.05.001
#ELPA_INCLUDE_DIR=$ELPA_HOME/include/elpa_openmp-2023.05.001

FileName="eigenproblem_elpa"

cc -std=c11 -I$ELPA_INCLUDE_DIR $FileName.c -o $FileName $ELPA_LIBS
