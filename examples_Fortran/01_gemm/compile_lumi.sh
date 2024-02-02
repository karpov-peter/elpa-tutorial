#set -v

module load PrgEnv-gnu
module load cray-libsci

FileName=gemm_big_matrices
ftn -fopenmp -o $FileName $FileName.f90
