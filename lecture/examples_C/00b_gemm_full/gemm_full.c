#include <stdio.h>
#include "cblas.h"

int main() 
{
int i, j;

// Define matrices A, B, and C
double A[6] = {1.0, 1.0, 
               2.0, 2.0, 
               3.0, 3.0}; // 3x2 matrix

double B[6] = {1.0, 2.0, 3.0, 
               1.0, 2.0, 3.0}; // 2x3 matrix
               
double C[9] = {0.0, 0.0, 0.0, 
               0.0, 0.0, 0.0, 
               0.0, 0.0, 0.0}; // 3x3 matrix

// Parameters for dgemm
int m = 3; // Number of rows in matrices A and C
int n = 3; // Number of columns in matrices B and C
int k = 2; // Number of columns in matrix A; number of rows in matrix B
double alpha = 0.5;
double beta = 0.0;
int lda = 2; // Leading dimension of matrix A
int ldb = 3; // Leading dimension of matrix B
int ldc = 3; // Leading dimension of matrix C

// Perform matrix multiplication: C = alpha * A * B + beta * C
cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);

// Print the resulting C matrix
printf("Resulting C matrix:\n");
for (i = 0; i < m; i++) 
    {
    for (j = 0; j < n; j++) 
        printf("%f ", C[i * n + j]);
  
    printf("\n");
    }


return 0;
}
