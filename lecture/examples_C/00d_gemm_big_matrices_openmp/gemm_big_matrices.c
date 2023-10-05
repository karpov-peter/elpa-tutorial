#include <stdio.h>
#include <stdlib.h>
#include "cblas.h"

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

    double t1 = omp_get_wtime();

    // Perform matrix multiplication: C = alpha * A * B + beta * C
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);

    double t2 = omp_get_wtime();

    // Print the resulting C matrix
    int print_size = 10;
    printf("Upper-left %dx%d corner of the resulting C matrix:\n", print_size, print_size);
    for (i = 0; i < print_size; i++) 
      {
      for (j = 0; j < print_size; j++) 
        printf("%f\t", C[i * n + j]);
      
      printf("\n");
      }

    
    printf("\nLower-right %dx%d corner of the A matrix:\n", print_size, print_size);
    for (i = N-print_size; i < N; i++) 
      {
      for (j = N-print_size; j < N; j++) 
        printf("%e\t", A[i * (long long int)n + j]);
      
      printf("\n");
      }

    printf("\nLower-right %dx%d corner of the B matrix:\n", print_size, print_size);
    for (i = N-print_size; i < N; i++) 
      {
      for (j = N-print_size; j < N; j++) 
        printf("%e\t", B[i * (long long int)n + j]);
      
      printf("\n");
      }

    printf("\nLower-right %dx%d corner of the resulting C matrix:\n", print_size, print_size);
    for (i = N-print_size; i < N; i++) 
      {
      for (j = N-print_size; j < N; j++) 
        printf("%e\t", C[i * (long long int)n + j]);
      
      printf("\n");
      }


    for (long long int k=0; k<N*N; k++)
      {
      if (C[k] < 0.01)
        {
        printf("C[%lld] = %e\n", k, C[k]);
        break;
        }
      }

    long long int bigValue = (long long int) N * (long long int) N - 1;
    printf("%lld\n", bigValue);
    printf("A last=%e\n", A[(long long int) N * (long long int) N - 1]);
    printf("B last=%e\n", B[(long long int) N * (long long int) N - 1]);
    printf("C last=%e\n", C[(long long int) N * (long long int) N - 1]);
    printf("\nGEMM time: %f\n", t2-t1);
    
    // Free dynamically allocated memory
    free(A);
    free(B);
    free(C);

    return 0;
}

