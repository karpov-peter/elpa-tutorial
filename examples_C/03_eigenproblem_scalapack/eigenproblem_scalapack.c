#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "mpi.h"

#include <assert.h>

static int max( int a, int b ){if (a>b) return(a); else return(b);}
static int min( int a, int b ){if (a<b) return(a); else return(b);}

//__________________________________________________________________________________
	
int N = 1000;
int nev=64;
int NB = 32; 

int debug_mode=0;

enum matrix_types {Symmetrized_Clement, Random_0_1};
enum matrix_types matrix_type = Symmetrized_Clement;

enum diagonalization_methods {scalapack_pdsyev, scalapack_pdsyevd, scalapack_pdsyevr, scalapack_pdsyevx};
enum diagonalization_methods diagonalization_method = scalapack_pdsyev;

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

// Driver routines: 
// they all do first tridiagonalization with Housholder transformations and then they differ only in a way how they solve a tridiagonal matrix
// ~4/3 N^3 flops if only eigenvalues are computed, ~8/3 N^3 flops if both eigenvalues and eigenvectors are computed

// Сomputes all eigenvalues and, optionally, all or no eigenvectors of a real symmetric matrix 
// uses the QR Algorithm (p.108 of ScaLAPACK User Guide)
void pdsyev_  (char *jobz , char *uplo , int *n , double *a , int *ia , int *ja , int *desca , double *w , double *z , int *iz , int *jz , int *descz , double *work , int *lwork , int *info ); 

// Computes all eigenvalues and eigenvectors of a real symmetric matrix
// Uses divide and conquer algorithm.
void pdsyevd_ (char *jobz , char *uplo , int *n , double *a , int *ia , int *ja , int *desca , double *w , double *z , int *iz , int *jz , int *descz , double *work , int *lwork , int *iwork , int *liwork , int *info ); 

// Computes selected eigenvalues and, optionally, eigenvectors of a real symmetric matrix
// Uses Relatively Robust Representation.
void pdsyevr_(char* jobz, char* range, char* uplo, int* n, double* a, int* ia, int* ja, int* desca, double* vl, double* vu, int* il, int* iu, int* m, int* nz, double* w, double* z, int* iz, int* jz, int* descz, double* work, int* lwork, int* iwork, int* liwork, int* info);

// Computes selected eigenvalues and, optionally, eigenvectors of a symmetric matrix.
// Uses bisection method (for eigenvalues) and inverse iteration (for eigenvectors)
void pdsyevx_ (char *jobz , char *range , char *uplo , int *n , double *a , int *ia , int *ja , int *desca , double *vl , double *vu , int *il , int *iu , double *abstol , int *m , int *nz , double *w , double *orfac , double *z , int *iz , int *jz , int *descz , double *work , int *lwork , int *iwork , int *liwork , int *ifail , int *iclustr , double *gap , int *info );
//__________________________________________________________________________________

double get_global_matrix_element (int I_gl, int J_gl)
   {
   double matrix_element=0;
   
   if (matrix_type==Symmetrized_Clement)
      {
      if (I_gl==J_gl+1 || J_gl==I_gl+1) 
         {
         int K = min(I_gl, J_gl);
         matrix_element = (double) sqrtl((K+1) * (N - K-1));
         }
      }
   else if (matrix_type==Random_0_1)
      {
      matrix_element = (double) (rand())/RAND_MAX;
      }
   else
      {
      printf("matrix_type=%d is not supported", matrix_type);
      exit(1);
      }
   
   return matrix_element;
   }

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
int descA[9], descZ[9]; // matrix descriptors
  
int izero=0, ione=1;

//____________________________________________ 

if (argc==1) // one argument was provided: filename (default)
   {
   N = 10;
   nev = N;
   }
   
if (argc==2) // two arguments were provided: filename (default), N
   {
   N   = atoi(argv[1]);
   nev = N;
   }

if (argc==3) // three arguments were provided: filename (default), N, nev
   {
   N   = atoi(argv[1]);
   nev = atoi(argv[2]);
   }

if (argc==4) // four arguments were provided: filename (default), N, nev, NB
   {
   N   = atoi(argv[1]);
   nev = atoi(argv[2]);
   NB  = atoi(argv[3]);
   }

//____________________________________________ 

if (N>100) debug_mode=0;

//srand(time(NULL));	
srand(1); // seeding
	
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

// try some other grids if you want
//if (world_size==72) {nprow=12; npcol=6;}

// specialy treat a case of small matrices
if (NB>N/max(nprow,npcol)) NB = N/max(nprow,npcol);

if (nprow*npcol!=world_size)
	{
  if (world_rank==0)
  printf("ERROR:  wrong grid size \n");
  MPI_Finalize(); 
	return(1);
  }

if (debug_mode) printf("world_rank= %i \n", world_rank);
if (world_rank==0) printf("world_size= %i \n", world_size);
if (world_rank==0) printf("world_size=%i, nprow=%i, npcol=%i, N=%i, nev=%i, NB=%i \n", world_size, nprow, npcol, N, nev, NB);

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

// Initialize the array descriptor for the distributed matrices A, Z
itemp = max( 1, m_loc );
descinit_( descA,  &N, &N, &NB, &NB, &izero, &izero, &ictxt, &itemp, &info );
descinit_( descZ,  &N, &N, &NB, &NB, &izero, &izero, &ictxt, &itemp, &info );


double *A_loc, *Z_loc, *Eigenvalues;

// Allocate the matrices A_loc, Z_loc (matrix of eigenvectors), Eigenvalues (vector of eigenvalues)
A_loc = (double *)calloc(m_loc*n_loc, sizeof(double)); // matrix to be diagonalized, local
Z_loc = (double *)calloc(m_loc*n_loc,sizeof(double));  // matrix of eigenvectors, local
Eigenvalues = (double *)calloc(N,sizeof(double)); // vector of eigenvalues, global
  

// fill the local part of matrix A 
int i_loc, j_loc;
for (i_loc = 0; i_loc < m_loc; i_loc+=1)
  {
  int l_1 = i_loc/NB; // local coord of the (NBxNB) block among other blocks
  int x_1 = i_loc%NB; // local coord within the block
  int I_gl= (l_1*nprow + myrow)*NB + x_1; // gl-global; nprow = "P_r"; myrow="p_r"  (quoted-ScaLAPACK userguide notation, p.88-90)
  for (j_loc = 0; j_loc < n_loc; j_loc+=1) 
    {
    int l_2 = j_loc/NB;
    int x_2 = j_loc%NB;
    int J_gl= (l_2*npcol + mycol)*NB + x_2;

    A_loc[i_loc + j_loc*m_loc] = get_global_matrix_element(I_gl, J_gl);
    }
  }

// print the local part of matrix A
if (debug_mode==1 && world_rank==0)
      {
      //cout << "A_loc:" << endl;
      printf("A_loc:\n");
      int i_loc, j_loc;
      for(i_loc=0; i_loc<m_loc; i_loc++) 
          {
          for(j_loc=0; j_loc<n_loc; j_loc++) 
              {
              printf("%g\t", A_loc[i_loc+j_loc*m_loc]);
              }
          //cout << endl;
          printf("\n");
          }
      //cout << endl;
      printf("\n");
      }


// auxillary array for pdsyev, pdsyevd, pdsyevr, pdsyevx
double* work = (double *)calloc(2,sizeof(double));
int lwork=-1;

// auxillary array for pdsyevd, pdsyevr, pdsyevx
int liwork=-1;
int* iwork = (int *)calloc(2,sizeof(int));

double t0, t1;

if (diagonalization_method == scalapack_pdsyev)
  {
  // "empty" diagonalization run for finding lwork, which is a required size of the array "work"
  pdsyev_( "V", "U", &N, A_loc, &ione, &ione, descA, Eigenvalues, Z_loc, &ione, &ione, descZ, work, &lwork, &info );
  lwork= (int) work[0];
    
  if (debug_mode) printf("lwork=%i", lwork);
  free(work);
  work = (double *)calloc(lwork, sizeof(double));
    
  // ED-main run
  t0 = MPI_Wtime();
  pdsyev_( "V", "U", &N, A_loc, &ione, &ione, descA, Eigenvalues, Z_loc, &ione, &ione, descZ, work, &lwork, &info );
  // arguments:
  // char jobz; //  "N" - only eigenvalues, "V" - eigenvalues and eigenvectors; A-??
  // char uplo; //  "U" - upper,  "L" - lower triangular matrix
  }
else if (diagonalization_method == scalapack_pdsyevd)
  {
  // "empty" diagonalization run for finding lwork and liwork, which are required sizes of the array "work" and "iwork"
  pdsyevd_("V", "U", &N, A_loc, &ione, &ione, descA, Eigenvalues, Z_loc, &ione, &ione, descZ, work, &lwork, iwork, &liwork, &info);
  lwork = (int)work[0];
  liwork = iwork[0];
  if (debug_mode) printf("world_rank=%d, lwork=%d, liwork=%d\n", world_rank, lwork, liwork);

  // Allocate work and iwork arrays with appropriate sizes
  work  = (double *)calloc( lwork, sizeof(double));
  iwork = (int *)   calloc(liwork, sizeof(int));

  // ED-main run
  t1 = MPI_Wtime();
  pdsyevd_("V", "U", &N, A_loc, &ione, &ione, descA, Eigenvalues, Z_loc, &ione, &ione, descZ, work, &lwork, iwork, &liwork, &info);
  
  }
else
  {
  printf("diagonalization_method=%d is not supported", diagonalization_method);
  MPI_Finalize();
  return(1);
  }

double t2 = MPI_Wtime();


if (world_rank==0 && info!=0) printf("info pdsyev(d)=%i, error occured!!!\n", info); 

if (world_rank==0)
  {
  printf("Diagonalization is done\n");
  
  int N_eigenvalues = min(10, nev);
  printf("First %d eigenvalues: \n", N_eigenvalues);
  int i; // for old Intel compilers
  for (i=0; i<N_eigenvalues; i++) printf("%g ",Eigenvalues[i]);
  printf("\n");
  
  if (matrix_type==Symmetrized_Clement)
    {
    printf("Absolute diff of eigenvalues wrt to analytic ones: \n");
    for (i=0; i<N_eigenvalues; i++) printf("%.16g ", Eigenvalues[i] - (-N+1+2*i) );
    printf("\n");
    }
  
  printf("Diagonalization time  (sec): \n%g\n", t2-t1);
  }


// cleanup and finalize

free(A_loc);
free(Z_loc);
free(Eigenvalues);

free(work);
free(iwork);


Cblacs_gridexit(0);

MPI_Finalize();

return(0);
}
