#set -v

module load cray-libsci
module load PrgEnv-gnu

FileName=gemm_big_matrices
ftn -fopenmp -o $FileName $FileName.f90
