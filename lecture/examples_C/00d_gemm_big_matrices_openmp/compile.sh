#set -v

module load cray-libsci
echo CRAY_LIBSCI_DIR=$CRAY_LIBSCI_DIR

#module load gcc/...
module load gcc/11.2.0
#gcc -o  00_gemm_concise 00_gemm_concise.c -llapack -lblas -lcblas

LIBS_CRAY="-L$CRAY_LIBSCI_DIR/CRAY/9.0/x86_64/lib -lsci_cray_mp -Wl,-rpath,$CRAY_LIBSCI_DIR/CRAY/9.0/x86_64/lib"

FileName=gemm_big_matrices
cc -fopenmp -o $FileName $FileName.c $LIBS_CRAY
