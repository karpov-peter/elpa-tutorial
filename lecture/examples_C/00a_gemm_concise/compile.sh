set -v

module load cray-libsci
echo CRAY_LIBSCI_DIR=$CRAY_LIBSCI_DIR

#module load gcc/...
module load PrgEnv-gnu
#gcc -o  00_gemm_concise 00_gemm_concise.c -llapack -lblas -lcblas

LIBS_CRAY="-lsci_cray_mpi"
LIBS_CRAY="-L$CRAY_LIBSCI_DIR/CRAY/9.0/x86_64/lib -lsci_cray_mpi -Wl,-rpath,$CRAY_LIBSCI_DIR/CRAY/9.0/x86_64/lib"

FileName=gemm_concise
cc $FileName.c -o $FileName $LIBS_CRAY
