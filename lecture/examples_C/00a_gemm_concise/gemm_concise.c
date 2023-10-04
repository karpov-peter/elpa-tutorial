#include "cblas.h"

int main() 
{
double a[6] = {1, 2, 1, 4, 5, 6}; 
double b[6] = {1, 2, 3, 4, 5, 6}; 
double c[9] = {0}; 

cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 2, 1.0, a, 2, b, 3, 0.0, c, 3);

return 0;
}
