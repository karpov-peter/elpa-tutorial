#include <stdio.h>
#include <stdlib.h>
#include <hip/hip_runtime.h>
#include <hipblas/hipblas.h>
#include "error_checks.h"

#include <omp.h> // for omp_get_wtime() 

int main(int argc, char *argv[]) 
{
    if (argc != 2) 
    {
        printf("Usage: %s <matrix_size>\n", argv[0]);
        return 1;
    }

    long long int i, j;
    long long int N = atol(argv[1]);

    hipblasHandle_t handle;
    CHECK_HIP_ERROR( hipblasCreate(&handle) );

    // Dynamically allocate matrices A, B, and C
    double *A = (double *)malloc(N * N * sizeof(double));
    double *B = (double *)malloc(N * N * sizeof(double));
    double *C = (double *)malloc(N * N * sizeof(double));

    double *dA;
    double *dB;
    double *dC;

    CHECK_HIP_ERROR( hipMalloc((void **) &dA, N * N * sizeof(double)) );
    CHECK_HIP_ERROR( hipMalloc((void **) &dB, N * N * sizeof(double)) );
    CHECK_HIP_ERROR( hipMalloc((void **) &dC, N * N * sizeof(double)) );

    for (j = 0; j < N; j++)
      {
      for (i = 0; i < N; i++) 
        {
        A[i + N*j] = i+1;
        B[i + N*j] = j+1;
        }
      }

    CHECK_HIP_ERROR( hipblasSetPointerMode(handle, HIPBLAS_POINTER_MODE_HOST) );

    // Parameters for dgemm
    int m = N;
    int n = N;
    int k = N;
    double alpha = 1.0/N;
    double beta = 0.0;
    int lda = N;
    int ldb = N;
    int ldc = N;

    const hipblasOperation_t transA = HIPBLAS_OP_N;
    const hipblasOperation_t transB = HIPBLAS_OP_N;

    CHECK_HIP_ERROR( hipMemcpy(dA, A, N * N * sizeof(double), hipMemcpyHostToDevice) );
    CHECK_HIP_ERROR( hipMemcpy(dB, B, N * N * sizeof(double), hipMemcpyHostToDevice) );

    // First execution is typically slow, so let's perform warm-up run
    // Perform matrix multiplication: C = alpha * A * B + beta * C
    CHECK_HIP_ERROR( hipblasDgemm(handle, transA, transB, m, n, k, &alpha, dA, lda, dB, ldb, &beta, dC, ldc) );
    hipDeviceSynchronize(); 
    double t1 = omp_get_wtime();

    // Perform matrix multiplication: C = alpha * A * B + beta * C
    CHECK_HIP_ERROR( hipblasDgemm(handle, transA, transB, m, n, k, &alpha, dA, lda, dB, ldb, &beta, dC, ldc) );
    hipDeviceSynchronize();

    double t2 = omp_get_wtime();
    CHECK_HIP_ERROR( hipMemcpy(C, dC, N * N * sizeof(double), hipMemcpyDeviceToHost) );

    // Print the resulting C matrix
    int print_size = 10;
    printf("Upper-left %dx%d corner of the resulting C matrix:\n", print_size, print_size);
    for (i = 0; i < print_size; i++) 
      {
      for (j = 0; j < print_size; j++) 
        printf("%f\t", C[i * n + j]);
      
      printf("\n");
      }


    printf("\nLower-right %dx%d corner of the resulting C matrix:\n", print_size, print_size);
    for (i = N-print_size; i < N; i++) 
      {
      for (j = N-print_size; j < N; j++) 
        printf("%e\t", C[i * (long long int)n + j]);
      
      printf("\n");
      }

    printf("\nGEMM time: %f\n", t2-t1);
    
    // Free dynamically allocated memory
    free(A);
    free(B);
    free(C);

    hipFree(dA);
    hipFree(dB);
    hipFree(dC);

    hipblasDestroy(handle);

    return 0;
}

