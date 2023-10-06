#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
//#include "scalapack.h"


static int max( int a, int b ){if (a>b) return(a); else return(b);}
static int min( int a, int b ){if (a<b) return(a); else return(b);}

//__________________________________________________________________________________
// default values	

int N = 1000; // matrix size
int NB = 32; // block size

int debug_mode=0;

//__________________________________________________________________________________
// blacs functions

void   Cblacs_pinfo( int* mypnum, int* nprocs);
void   Cblacs_get( int context, int request, int* value);
int    Cblacs_gridinit( int* context, char * order, int np_row, int np_col);
void   Cblacs_gridinfo( int context, int*  np_row, int* np_col, int*  my_row, int*  my_col);
void   Cblacs_gridexit( int context);
void   Cblacs_exit( int error_code);

// scalapack functions
// find size of the local matrix (NUMber of Rows Or Columns), m_loc or n_loc
int  numroc_( int *N, int *NB, int *iproc, int *isrcproc, int *nprocs);

// initialize the matrix descriptor
void descinit_( int *desc, int *M, int *N, int *mb, int *NB, int *irsrc, int *icsrc, int *ictxt, int *lld, int *info);

// redistribute the matrix from a global one to the local ones
void pdgemr2d_ (int *m , int *n , double *a , int *ia , int *ja , int *desca , 
                                  double *b , int *ib , int *jb , int *descb , int *ictxt );

// parallel matrix-matrix multiplication
void pdgemm_ (const char *transa, const char *transb, const int *m , const int *n, const int *k, const double *alpha, 
              const double *a, const int *ia, const int *ja, const int *desca, 
              const double *b, const int *ib, const int *jb, const int *descb, const double *beta , 
              double *c , const int *ic , const int *jc , const int *descc );

//__________________________________________________________________________________

int main(int argc, char **argv) 
{
//____________________________________________ 
// Initialize MPI
int world_rank, world_size; // MPI

MPI_Init( &argc, &argv);
MPI_Comm_size(MPI_COMM_WORLD, &world_size);
MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

//____________________________________________ 

int iam, nprocs; // blacs
int ictxt, nprow, npcol, myrow, mycol;
int np, nq;
int m_loc, n_loc;
int info, itemp, seed;
int descA[9], descB[9], descC[9], descGlobal[9]; // matrix descriptors
long long int I_gl, J_gl;

int izero=0, ione=1;

//____________________________________________ 
// Read command line arguments (assuming arguments are passed as N, NB)

if (argc==1) // one argument was provided: filename (default)
   {
   N = 4;
   NB = 2;
   }
   
if (argc==2) // two arguments were provided: filename (default), N
   {
   N   = atoi(argv[1]);
   NB = 2;
   }

if (argc==3) // three arguments were provided: filename (default), N, NB
   {
   N   = atoi(argv[1]);
   NB = atoi(argv[2]);
   }

//____________________________________________ 

if (N>100) debug_mode=0; // turn off debug mode for large matrices
	
//____________________________________________ 
// ! Compute grid size

// try to set square grid
nprow = sqrt(world_size);
npcol = sqrt(world_size);

// if the grid is not square, try to find the "most square" rectangular grid
for(; nprow>0; nprow--)
  {
  if (world_size%nprow==0)
    {
    npcol = world_size/nprow;
    break;
    }
  }


if (NB>N/max(nprow,npcol)) NB = N/max(nprow,npcol); // case of small matrices

if (nprow*npcol!=world_size)
	{
  if (world_rank==0)
  printf("ERROR:  wrong grid size \n");
  MPI_Finalize(); 
	return(1);
  }

if (debug_mode) printf("world_rank= %i \n", world_rank);
if (world_rank==0) printf("world_size= %i \n", world_size);
if (world_rank==0) printf("world_size=%i, nprow=%i, npcol=%i, N=%i, NB=%i \n", world_size, nprow, npcol, N, NB);

//____________________________________________ 

// Setup blacs grid
//Cblacs_pinfo( &iam, &nprocs ) ;
Cblacs_get( -1, 0, &ictxt );
Cblacs_gridinit( &ictxt, "Row", nprow, npcol );
Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );

// Compute the size of the local matrices (thanks to numroc)
m_loc = numroc_( &N, &NB, &myrow, &izero, &nprow );
n_loc = numroc_( &N, &NB, &mycol, &izero, &npcol );
if (debug_mode) printf("myrow=%i, mycol=%i, m_loc=%i, n_loc=%i, \n", myrow, mycol, m_loc, n_loc);

// Initialize the array descriptor for the distributed matrices A/A_loc, B/B_loc, C/C_loc
itemp = max( 1, m_loc );
descinit_( descA,  &N, &N, &NB, &NB, &izero, &izero, &ictxt, &itemp, &info );
descinit_( descB,  &N, &N, &NB, &NB, &izero, &izero, &ictxt, &itemp, &info );
descinit_( descC,  &N, &N, &NB, &NB, &izero, &izero, &ictxt, &itemp, &info );

// Allocate memory for the local matrices
double *A_loc = (double *)malloc(m_loc * n_loc * sizeof(double));
double *B_loc = (double *)malloc(m_loc * n_loc * sizeof(double));
double *C_loc = (double *)malloc(m_loc * n_loc * sizeof(double));

//____________________________________________ 

// Fill global matrices A and B. This is not a very efficient way to do it, 
// it's better to initialize the local matrices in parallel -- this is covered in the next example
// For saving memory, we will fill the global matrices only on 0-th MPI rank
double *A, *B, *C;
if (world_rank==0)
  {
  A = (double *)malloc(N * N * sizeof(double));
  B = (double *)malloc(N * N * sizeof(double));

  for (int J_gl = 0; J_gl < N; J_gl++) 
    {
    for (int I_gl = 0; I_gl < N; I_gl++) 
      {
      A[I_gl + N*J_gl] = I_gl+1;
      B[I_gl + N*J_gl] = J_gl+1;
      }
    }
  }

// create a desccriptor for the global matrix
descinit_( descGlobal, &N, &N, &N, &N, &izero, &izero, &ictxt, &N, &info );

// redestribute the global matrices to the local ones
pdgemr2d_(&N, &N, A,     &ione, &ione, descGlobal,
                  A_loc, &ione, &ione, descA,  &ictxt);
pdgemr2d_(&N, &N, B,     &ione, &ione, descGlobal,
                  B_loc, &ione, &ione, descB,  &ictxt);

// we don't need the global matrices A and B anymore and can deallocate them to save memory
free(A);
free(B);

//____________________________________________ 
// perform parallel matrix-matrix multiplication

double t_start = MPI_Wtime();

char no = 'N'; // C = A * B, both matrices A and B are not transposed
double alpha = 1.0/N;
double beta = 0.0;
pdgemm_(&no, &no, &N, &N, &N, &alpha, A_loc, &ione, &ione, descA, B_loc, &ione, &ione, descB, &beta, C_loc, &ione, &ione, descC);

double t_stop = MPI_Wtime();

//____________________________________________ 
// redistribute the local matrix C_loc to the global matrix C
C = (double *)malloc(N * N * sizeof(double));

pdgemr2d_(&N, &N, C_loc,&ione, &ione, descC,
                  C,    &ione, &ione, descGlobal,  &ictxt);

// Output C matrix for the process with row 0 and col 0
if (world_rank==0)
  {
  int print_size = min(10, N);
  printf("\nUpper left corner of the global matrix C:\n");
  for (int I_gl = 0; I_gl < print_size; I_gl++) 
      {
      for (int J_gl = 0; J_gl < print_size; J_gl++) 
        {
        printf("%f\t", C[I_gl + J_gl*N]);
        }
      printf("\n");
      }

  printf("\nLower-right corner of the global matrix C:\n");
  for (int I_gl = N-print_size; I_gl < N; I_gl++) 
      {
      for (int J_gl = N-print_size; J_gl < N; J_gl++) 
        {
        printf("%f\t", C[I_gl + J_gl*N]);
        }
      printf("\n");
      }
  
  printf("\nPDGEMM time (sec): \n%g\n", t_stop-t_start);
  }

free(A_loc);
free(B_loc);
free(C_loc);
free(C);

Cblacs_gridexit(ictxt);
MPI_Finalize();

return 0;
}
