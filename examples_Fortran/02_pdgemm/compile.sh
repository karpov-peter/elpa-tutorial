#set -v

module load PrgEnv-gnu
module load cray-libsci

FileName=pdgemm
ftn  $FileName.f90 -o $FileName
