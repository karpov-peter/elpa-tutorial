#include <stdio.h>
#include <stdlib.h>
#include "cblas.h"

int main(int argc, char *argv[]) 
{
    if (argc != 2) 
    {
        printf("Usage: %s <matrix_size>\n", argv[0]);
        return 1;
    }

    int i, j;
    int N = atoi(argv[1]);

    // Dynamically allocate matrices A, B, and C
    double *A = (double *)malloc(N * N * sizeof(double));
    double *B = (double *)malloc(N * N * sizeof(double));
    double *C = (double *)malloc(N * N * sizeof(double));

    for (j = 0; j < N; j++)
      {
      for (i = 0; i < N; i++) 
        {
        A[i + N*j] = i+1;
        B[i + N*j] = j+1;
        }
      }

    // Parameters for dgemm
    int m = N;
    int n = N;
    int k = N;
    double alpha = 1.0/N;
    double beta = 0.0;
    int lda = N;
    int ldb = N;
    int ldc = N;

    // Perform matrix multiplication: C = alpha * A * B + beta * C
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);

    // Print the resulting C matrix
    int print_size = 10;
    printf("Upper-left %dx%d corner of the resulting C matrix:\n", print_size, print_size);
    for (i = 0; i < print_size; i++) 
      {
      for (j = 0; j < print_size; j++) 
        printf("%f\t", C[i * n + j]);
      
      printf("\n");
      }

    printf("Lower-right %dx%d corner of the resulting C matrix:\n", print_size, print_size);
    for (i = N-print_size; i < N; i++) 
      {
      for (j = N-print_size; j < N; j++) 
        printf("%f\t", C[i * n + j]);
      
      printf("\n");
      }

    // Free dynamically allocated memory
    free(A);
    free(B);
    free(C);

    return 0;
}

