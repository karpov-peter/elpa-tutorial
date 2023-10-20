#!/usr/bin/bash

set +v

module purge
module load intel/21.2.0  
module load impi/2021.2
module load mkl/2021.2
module load cuda/11.2

set -v

module list -l

echo $MKLROOT

ELPA_HOME=$HOME/soft/elpa_single_precision_cuda
ELPA_LIBS="-L$ELPA_HOME/lib -lelpa -Wl,-rpath=$ELPA_HOME/lib"


MKL_LIBS="-L$MKL_HOME/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm -Wl,-rpath,$MKL_HOME/lib/intel64"

echo $ELPA_HOME
echo $ELPA_LIBS

MKLINCLUDE=$MKLROOT/include
MKLLIBINTEL64=$MKLROOT/lib/intel64

MPLAPACK_HOME=$HOME/soft/MPLAPACK/install
MPLAPACK_LIBS="-L$MPLAPACK_HOME/lib -lmpblas__Float128 -lmpblas_dd -lqd -lmpblas__Float64x -Wl,--rpath,$MPLAPACK_HOME/lib"

CUDA_LIBS="-L$CUDA_HOME/lib64 -lcudart -lcuda -Wl,-rpath=$CUDA_HOME/lib64"


filename="main"


mpiicc -DUSE_CUDA -I$CUDA_HOME/include -I$MPLAPACK_HOME/include -I$MPLAPACK_HOME/include/mplapack -I$ELPA_HOME/include/elpa-2022.11.001.rc2 -I$MKLROOT/include  $filename.c -o $filename $MPLAPACK_LIBS $ELPA_LIBS $MKL_LIBS $CUDA_LIBS
#mpiicc -std=c11 -I${ChASEROOT}/include -I$ELPA_HOME/include/elpa-2022.05.001 -I$MKLROOT/include  $filename.c -o $filename $ELPA_LIBS $MKL_LIBS


# debug
#mpiicpc -g -std=c++11 -I${ChASEROOT}/include -I$ELPA_HOME/include/elpa-2022.05.001 -I$MKLROOT/include  $filename.cpp -o $filename $ELPA_LIBS $MKL_LIBS

#expand preprocessor commands
#mpiicpc -std=c++11 -I$ELPA_HOME/include/elpa-2022.05.001 -I$MKLROOT/include  -E $filename.cpp > $filename-expanded.cpp 
