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
int  numroc_( int *N, int *NB, int *iproc, int *isrcproc, int *nprocs);
void descinit_( int *desc, int *M, int *N, int *mb, int *NB, int *irsrc, int *icsrc, int *ictxt, int *lld, int *info);

void pdgemm_ (const char *transa, const char *transb, const int *m , const int *n, const int *k, const double *alpha, 
              const double *a, const int *ia, const int *ja, const int *desca, 
              const double *b, const int *ib, const int *jb, const int *descb, const double *beta , 
              double *c , const int *ic , const int *jc , const int *descc );

//__________________________________________________________________________________

int main(int argc, char **argv) 
{
//____________________________________________ 
// setup MPI
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
int descA[9], descB[9], descC[9]; // matrix descriptors
  
int izero=0, ione=1;

//____________________________________________ 

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
// determine the grid size

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

// setup Blacs
//Cblacs_pinfo( &iam, &nprocs ) ;
Cblacs_get( -1, 0, &ictxt );
Cblacs_gridinit( &ictxt, "Row", nprow, npcol );
Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );

// Compute the size of the local matrices (thanks to numroc)
m_loc = numroc_( &N, &NB, &myrow, &izero, &nprow );
n_loc = numroc_( &N, &NB, &mycol, &izero, &npcol );
if (debug_mode) printf("myrow=%i, mycol=%i, m_loc=%i, n_loc=%i, \n", myrow, mycol, m_loc, n_loc);

// Initialize the array descriptor for the distributed matrices A, Z
itemp = max( 1, m_loc );
descinit_( descA,  &N, &N, &NB, &NB, &izero, &izero, &ictxt, &itemp, &info );
descinit_( descB,  &N, &N, &NB, &NB, &izero, &izero, &ictxt, &itemp, &info );
descinit_( descC,  &N, &N, &NB, &NB, &izero, &izero, &ictxt, &itemp, &info );

double *A = (double *)malloc(m_loc * n_loc * sizeof(double));
double *B = (double *)malloc(m_loc * n_loc * sizeof(double));
double *C = (double *)malloc(m_loc * n_loc * sizeof(double));

// Fill A and B with sample data (for brevity, just set all elements to 1.0)
for (int i = 0; i < m_loc * n_loc; i++) 
  {
  A[i] = 1.0;
  B[i] = 1.0;
  C[i] = 0.0;
  }

char no = 'N';
double alpha = 1.0;
double beta = 0.0;
pdgemm_(&no, &no, &N, &N, &N, &alpha, A, &ione, &ione, descA, B, &ione, &ione, descB, &beta, C, &ione, &ione, descC);

// Output C matrix for the process with row 0 and col 0
int print_size = min(10, N);
printf("Upper left corner of C:\n");
if (myrow == 0 && mycol == 0) {
    for (int i = 0; i < print_size; i++) {
        for (int j = 0; j < print_size; j++) {
            printf("%f ", C[i + j * m_loc]);
        }
        printf("\n");
    }
}

free(A);
free(B);
free(C);
Cblacs_gridexit(ictxt);
MPI_Finalize();

return 0;
}
