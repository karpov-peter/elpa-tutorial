#set -v

module load PrgEnv-gnu
module load cray-libsci

FileName=pdgemm
cc  $FileName.c -o $FileName
