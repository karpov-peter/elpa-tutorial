// Description: This code solves the eigenproblem using ELPA library with GPU acceleration.
// In this example, the device pointers are passed to the ELPA solver function
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "mpi.h"

#include <cuda.h>
#include <cuda_runtime.h>

// Step 1: include the ELPA header file
#include <elpa/elpa.h>

#include <assert.h>
#define assert_elpa_ok(x) assert(x == ELPA_OK)

static int max( int a, int b ){if (a>b) return(a); else return(b);}
static int min( int a, int b ){if (a<b) return(a); else return(b);}

//__________________________________________________________________________________

int N   = 1000;
int nev = 64;
int NB  = 32; 

int debug_mode=0;

enum matrix_types {Symmetrized_Clement, Random_0_1};
enum matrix_types matrix_type = Symmetrized_Clement;

//__________________________________________________________________________________
// BLACS routines
void   Cblacs_pinfo( int* mypnum, int* world_size);
void   Cblacs_get( int context, int request, int* value);
int    Cblacs_gridinit( int* context, char * order, int np_row, int np_col);
void   Cblacs_gridinfo( int context, int*  np_row, int* np_col, int*  my_row, int*  my_col);
void   Cblacs_gridexit( int context);
void   Cblacs_exit( int error_code);

// ScaLAPACK routines (Fortran-style, hence the trailing underscore is needed)
int    numroc_( int *n, int *NB, int *iproc, int *isrcproc, int *world_size);
void   descinit_( int *desc, int *m, int *n, int *mb, int *NB, int *irsrc, int *icsrc, int *ictxt, int *lld, int *info);

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
// Setup MPI
int world_rank, world_size; // MPI

// For hybrid MPI+OpenMP
//int thread_level;
//MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &thread_level); 

// For pure MPI
MPI_Init( &argc, &argv);

MPI_Comm_size(MPI_COMM_WORLD, &world_size);
MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

//____________________________________________ 

int ictxt, nprow, npcol, myrow, mycol;
int np, nq;
int m_loc, n_loc;
int info, itemp, seed, lwork;
int descA[9], descZ[9];
int izero=0;

//____________________________________________ 
// Parse command line arguments

if (argc==1) // one argument was provided: filename (default)
   {
   N = 10;
   nev = N;
   }
   
if (argc==2) // two arguments were provided: filename (default), N (matrix size)
   {
   N   = atoi(argv[1]);
   nev = N;
   }

if (argc==3) // three arguments were provided: filename (default), N (matrix size), nev
   {
   N   = atoi(argv[1]);
   nev = atoi(argv[2]);
   }

if (argc==4) // four arguments were provided: filename (default), N (matrix size), nev, NB ()
   {
   N   = atoi(argv[1]);
   nev = atoi(argv[2]);
   NB  = atoi(argv[3]);
   }

//____________________________________________ 
// Determine the grid size

nprow = sqrt(world_size); 
npcol = sqrt(world_size); 
if (world_size==2)  {nprow=2; npcol=1;}
if (world_size==6)  {nprow=3; npcol=2;}

if (world_size==10) {nprow=2; npcol=5;}
if (world_size==20) {nprow=4; npcol=5;}
if (world_size==40) {nprow=5; npcol=8;}

if (world_size==72) {nprow=8; npcol=9;}

if (world_size==72*4) {nprow=8*2; npcol=9*2;}

if (NB>N/max(nprow,npcol)) NB = N/max(nprow,npcol);

if (nprow*npcol!=world_size)
	{
  if (world_rank==0)
  printf("ERROR:  wrong grid \n");
  MPI_Finalize(); 
	return(1);
  }

if (debug_mode) printf("world_rank= %i \n", world_rank);
if (world_rank==0) printf("world_size=%i, nprow=%i, npcol=%i, N=%i, nev=%i, NB=%i \n", world_size, nprow, npcol, N, nev, NB);

//____________________________________________ 
	
// Setup BLACS
Cblacs_pinfo( &world_rank, &world_size ) ;
Cblacs_get( -1, 0, &ictxt );
Cblacs_gridinit( &ictxt, "Col", nprow, npcol ); // "Row" or "Col" is the ordering of the processes in the grid. ELPA works with either of them
Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );

// Compute the size of the local matrices A_loc, Z_loc (thanks to numroc)
m_loc = numroc_( &N, &NB, &myrow, &izero, &nprow );
n_loc = numroc_( &N, &NB, &mycol, &izero, &npcol );
if (debug_mode) printf("myrow=%i, mycol=%i, m_loc=%i, n_loc=%i, \n", myrow, mycol, m_loc, n_loc);

// Initialize the array descriptors descA and descZ for the distributed matrices A and Z
itemp = max( 1, m_loc );
descinit_( descA,  &N, &N, &NB, &NB, &izero, &izero, &ictxt, &itemp, &info );
descinit_( descZ,  &N, &N, &NB, &NB, &izero, &izero, &ictxt, &itemp, &info );

//____________________________________________ 
// Setup ELPA 

// Step 2: define a handle for an ELPA object
elpa_t handle;
int status;

// Step 3: initialize the ELPA library
status = elpa_init(20170403);   
if (status != ELPA_OK) 
	{
	fprintf(stderr, "Error: ELPA API version not supported");
  MPI_Finalize();
	exit(1);
	}

// Step 4: allocate the ELPA object
handle = elpa_allocate(&status);
if (status != ELPA_OK) 
	{
	fprintf(stderr, "Error: cannot allocate the ELPA object");
  MPI_Finalize();
	exit(1);
	}

//____________________________________________ 
// Step 5: set mandatory parameters describing the matrix and its MPI distribution
// This can be done only once per ELPA object (handle)
elpa_set(handle, "na", N, &status); // matrix dimension
if (world_rank==0 || debug_mode) printf("na=%i is set, status=%i\n", N, status);

elpa_set(handle, "nev", nev, &status); // number of eigenvectors to be calculated (all eigenvalues are calculated regardless this setting)
if (world_rank==0 || debug_mode) printf("nev=%i is set, status=%i\n", nev, status);

elpa_set(handle, "local_nrows", m_loc, &status); // m_loc
if (world_rank==0 || debug_mode) printf("local_nrows=%i is set, status=%i\n", m_loc, status);

elpa_set(handle, "local_ncols", n_loc, &status); // n_loc
if (world_rank==0 || debug_mode) printf("local_ncols=%i is set, status=%i\n", n_loc, status);

elpa_set(handle, "nblk", NB, &status); // NB
if (world_rank==0 || debug_mode) printf("nblk=%i is set, status=%i\n", NB, status);

elpa_set(handle, "mpi_comm_parent", MPI_Comm_c2f(MPI_COMM_WORLD), &status);
if (world_rank==0 || debug_mode) printf("Fortran MPI communicator is set, status=%i\n", status);

elpa_set(handle, "process_row", myrow, &status); // myrow
if (world_rank==0 || debug_mode) printf("process_row=%i is set, status=%i\n", myrow, status);

elpa_set(handle, "process_col", mycol, &status); // mycol
if (world_rank==0 || debug_mode) printf("process_col=%i is set, status=%i\n", mycol, status);

// Step 6: Finalize the setup of mandatory parameters of the ELPA object (handle)
status = elpa_setup(handle);
printf("ELPA setup done, status=%i\n", status);
if (status!=ELPA_OK)
  {
  printf("elpa_setup() failed, status=%i\n", status);
  MPI_Finalize();
  return(1);
  }


// Step 7: set ELPA runtime options. They can be changed between different ELPA runs, e.g. elpa_eigenvectors() calls
elpa_set(handle, "solver", ELPA_SOLVER_1STAGE, &status);
if (world_rank==0 || debug_mode) printf("solver=ELPA_SOLVER_1STAGE is set, status=%i\n", status);

elpa_set(handle, "nvidia-gpu", 1, &status); // 1=on, 0=off
if (status!=ELPA_OK)
  {
  printf("elpa_setup(nvidia-gpu) failed, status=%i\n", status);
  MPI_Finalize();
  return(1);
  }

// Finalize the GPU setup, needed only when using GPUs
status = elpa_setup_gpu(handle);
if (status!=ELPA_OK)
  {
  printf("elpa_setup_gpu() failed, status=%i\n", status);
  MPI_Finalize();
  return(1);
  }
// End of Setup ELPA 
//____________________________________________ 


// Allocate the matrices A_loc, Z_loc(vector of eigenvalues), Eigenvalues
double *A_loc, *Z_loc, *work, *Eigenvalues;
A_loc = (double *)calloc(m_loc*n_loc,sizeof(double)) ;
Z_loc = (double *)calloc(m_loc*n_loc,sizeof(double)) ;
Eigenvalues = (double *)calloc(N,sizeof(double)) ;
  

// Fill the local part of matrix A 
int k = 0;
int i_loc, j_loc;
for(int i_loc=0; i_loc<m_loc; i_loc++) // iteration over rows of A_loc
    {
    int l_1 = i_loc/NB; // local coord of the (NBxNB) block among other blocks
    int x_1 = i_loc%NB; // local coord within the block
    int I_gl= (l_1*nprow + myrow)*NB + x_1; // gl-global; nprow = "P_r"; myrow="p_r"  (quoted-ScaLAPACK userguide notation, p.88-90)

    for(int j_loc=0; j_loc<n_loc; j_loc++) //  iteration over columns of A_loc
        {
        int l_2 = j_loc/NB; // local coord of the (NBxNB) block among other blocks
        int x_2 = j_loc%NB; // local coord within the block
        int J_gl= (l_2*npcol + mycol)*NB + x_2;

        // ELPA assumes column-major matrix layout
        A_loc[i_loc+j_loc*m_loc] = get_global_matrix_element(I_gl, J_gl);
        }
    }

//____________ CUDA-part ________________

// Allocate device memory and perform host-to-device copying
double *A_loc_dev, *Z_loc_dev, *Eigenvalues_dev;

cudaError_t error_cuda = cudaMalloc((void **) &A_loc_dev , m_loc*n_loc*sizeof(double));
if (error_cuda != cudaSuccess){    
   fprintf(stderr, "Error in cudaMalloc(A_loc_dev)\n");
   exit(1);
   }

error_cuda = cudaMalloc((void **) &Z_loc_dev , m_loc*n_loc*sizeof(double));
if (error_cuda != cudaSuccess){    
   fprintf(stderr, "Error in cudaMalloc(Z_loc_dev)\n");
   exit(1);
   }

error_cuda = cudaMalloc((void **) &Eigenvalues_dev , N*sizeof(double));
if (error_cuda != cudaSuccess){
   fprintf(stderr, "Error in cudaMalloc(Eigenvalues_dev)\n");
   exit(1);
   }

// copy
error_cuda = cudaMemcpy(A_loc_dev, A_loc,  m_loc*n_loc*sizeof(double), cudaMemcpyHostToDevice);
if (error_cuda != cudaSuccess){   
   fprintf(stderr, "Error in cudaMemcpy(A_loc_dev, A_loc)\n");
   exit(1);
   }  

//____________________________________________ 
// Print ELPA settings
//elpa_print_settings(handle, &status);

//____________________________________________ 
// Step 8: Perform the diagonalization using ELPA

MPI_Barrier(MPI_COMM_WORLD); // for timing
double t_start = MPI_Wtime();

// A_loc=local part of the matrix to be diagonalized, input
// Eigenvalues=global array of eigenvalues, output
// Z_loc=matrix of eigenvectors (contained as columns), output
elpa_eigenvectors_double(handle, A_loc_dev, Eigenvalues_dev, Z_loc_dev, &status);

if (status!=ELPA_OK)
  {
  printf("elpa_eigenvectors() failed, status=%i\n", status);
  MPI_Finalize();
  return(1);
  }

MPI_Barrier(MPI_COMM_WORLD);
double t_stop = MPI_Wtime();

//____________________________________________ 
// copying back from device to host

error_cuda = cudaMemcpy(Z_loc, Z_loc_dev, m_loc*n_loc*sizeof(double), cudaMemcpyDeviceToHost);
if (error_cuda != cudaSuccess){    
   fprintf(stderr, "Error in cudaMemcpy(Z_loc, Z_loc_dev)\n");
   exit(1);
   }

error_cuda = cudaMemcpy(Eigenvalues, Eigenvalues_dev, N*sizeof(double), cudaMemcpyDeviceToHost);
if (error_cuda != cudaSuccess){    
   fprintf(stderr, "Error in cudaMemcpy(Eigvalues, Eigenvalues_dev)\n");
   exit(1);
   }

//____________________________________________ 
// Print the results	

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
  
  printf("Diagonalization time (sec): \n%g\n", t_stop-t_start);
  }

//____________________________________________ 
// Step 9: Clean up ELPA

elpa_deallocate(handle, &status);
elpa_uninit(&status);

//____________________________________________ 
// deallocating device pointers
error_cuda = cudaFree(A_loc_dev);
if (error_cuda != cudaSuccess){    
   fprintf(stderr, "Error in cudaFree(A_loc_dev)\n");
   exit(1);
   } 

error_cuda = cudaFree(Z_loc_dev);
if (error_cuda != cudaSuccess){    
   fprintf(stderr, "Error in cudaFree(Z_loc_dev)\n");
   exit(1);
   } 

error_cuda = cudaFree(Eigenvalues_dev);
if (error_cuda != cudaSuccess){    
   fprintf(stderr, "Error in cudaFree(Eigenvalues_dev)\n");
   exit(1);
   } 

//____________________________________________ 
// Clean up the rest and finalize

free(A_loc);
free(Z_loc);
free(Eigenvalues);

Cblacs_gridexit(0);
MPI_Finalize();

return(0);
}
