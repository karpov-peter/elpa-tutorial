#include <mpi.h>
#include "mkl.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include <cuda.h>
#include <cuda_runtime.h>

#include <elpa/elpa.h>

int ELPA_SOLVER = ELPA_SOLVER_1STAGE;

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))


//______________________________________________

int debug_mode=0;

int N; //problem size
int nev; //number of eigenpairs to be computed

enum matrix_types {Symmetrized_Clement, Random_0_1, Random_m1_1, Sequence};
enum matrix_types matrix_type = Symmetrized_Clement;

/************   MPI  ***************************/

int world_rank, world_size;

/************  BLACS ***************************/

int NB=64; // Global column block size for partitioning the global matrix
int ictxt; // ordinal number i of the "context",
int nprow; // Number of total process rows    in the process grid (equivalent to Pr)
int npcol; // Number of total process columns in the process grid (equivalent to Pc)

int myrow; // The calling process's row     coordinate in the process grid
int mycol; // The calling process's column coordinate in the process grid.
int info=0; // Output integer argument of driver and computational routines indicating the success(0) or failure(1) of the routine

int m_loc; // size of local matrices is m_loc x n_loc for NxN global matrix; size of local vectors is m_loc x 1
int n_loc;

int m_loc_reduced; // size of local matrices is m_loc_reduced x n_loc_reduced for nev x nev global matrix
int n_loc_reduced;

int desc_N_N[9], desc_N_nev[9], desc_nev_nev[9];
#ifdef __cplusplus
extern "C" 
	{
#endif
	// Cblacs declarations
	void Cblacs_pinfo(int*, int*);
	void Cblacs_get(int, int, int*);
	void Cblacs_gridinit(int*, const char*, int, int);
	void Cblacs_pcoord(int, int, int*, int*);
	void Cblacs_gridexit(int);
	void Cblacs_barrier(int, const char*);
	int numroc_(int*, int*, int*, int*, int*);
	void Cblacs_gridinfo(int, int*, int*, int*, int*);
	void descinit_ (MKL_INT *desc, const MKL_INT *m, const MKL_INT *n, const MKL_INT *mb, const MKL_INT *nb, const MKL_INT *irsrc, const MKL_INT *icsrc, const MKL_INT *ictxt, const MKL_INT *lld, MKL_INT *info);
	
   // P BLAS declarations
   void pdgemm_ ( const char *transa , const char *transb , const MKL_INT *m , const MKL_INT *n , const MKL_INT *k , const double *alpha , const double *a , const MKL_INT *ia , const MKL_INT *ja , const MKL_INT *desca , const double *b , const MKL_INT *ib , const MKL_INT *jb , const MKL_INT *descb , const double *beta , double *c , const MKL_INT *ic , const MKL_INT *jc , const MKL_INT *descc );
   void pdscal_ ( const MKL_INT *n , const double *a , double *x , const MKL_INT *ix , const MKL_INT *jx , const MKL_INT *descx , const MKL_INT *incx );
   void pdaxpy_ ( const MKL_INT *n , const double *a , const double *x , const MKL_INT *ix , const MKL_INT *jx , const MKL_INT *descx , const MKL_INT *incx , double *y , const MKL_INT *iy , const MKL_INT *jy , const MKL_INT *descy , const MKL_INT *incy );
   void pdnrm2_ ( const MKL_INT *n , double *norm2 , const double *x , const MKL_INT *ix , const MKL_INT *jx , const MKL_INT *descx , const MKL_INT *incx );
	
   //void cblas_dgemm ( const CBLAS_LAYOUT Layout , const CBLAS_TRANSPOSE transa , const CBLAS_TRANSPOSE transb , const MKL_INT m , const MKL_INT n , const MKL_INT k , const double alpha , const double *a , const MKL_INT lda , const double *b , const MKL_INT ldb , const double beta , double *c , const MKL_INT ldc );
#ifdef __cplusplus
	}
#endif

//____________________________________________________________________________________________

int LoadCblacs(int N) // return 0 (for successfull) and 1 (for unsuccessful) initialization
	{
	if (world_size==72)
      {
      nprow = 9;
      npcol = 8;
      }
   else
      {      
      nprow=sqrt(world_size); // Number of process rows     in the process grid (equivalent to Pr)
      if (nprow>1 && NB*nprow>N) nprow = N/NB;
      npcol=nprow; // Number of process columns in the process grid (equivalent to Pc)
      }
      
	if (nprow==1 && npcol==1) NB=N;
	if (world_rank==0) printf("N= %d, nprow=%d, npcol=%d, NB=%d \n", N, nprow, npcol, NB);

	Cblacs_pinfo( &world_rank, &world_size ) ; // Routine is used when some initial system information is required before the BLACS are set up
	Cblacs_get( -1, 0, &ictxt );
	Cblacs_gridinit( &ictxt, "Row", nprow, npcol );
	Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );

	if (!((myrow>-1)&(mycol>-1)&(myrow<nprow)&(mycol<npcol))) return 1;

	// NUMROC computes the NUMber of Rows Or Columns of a distributed (small local) matrix owned by the process indicated by myrow,mycol
	int iZERO=0;
	m_loc = numroc_( &N, &NB, &myrow, &iZERO, &nprow ); // size of local matrix is m_loc x n_loc
	n_loc = numroc_( &N, &NB, &mycol, &iZERO, &npcol );
   
  m_loc_reduced = numroc_( &nev, &NB, &myrow, &iZERO, &nprow );
	n_loc_reduced = numroc_( &nev, &NB, &mycol, &iZERO, &npcol );
   
	//Cblacs_barrier(ictxt, "All");

	// DESCINIT initializes the descriptor vector with the 8 input arguments: M, N, MB, NB, IRSRC, ICSRC, ICTXT, LLD.
	descinit_(desc_N_N,   &N, &N  ,    &NB, &NB,   &iZERO, &iZERO, &ictxt, &m_loc, &info); // m_loc -- is local leading dimension
	descinit_(desc_N_nev, &N, &nev,    &NB, &NB,   &iZERO, &iZERO, &ictxt, &m_loc, &info); // m_loc -- is local leading dimension
	descinit_(desc_nev_nev, &nev, &nev,    &NB, &NB,   &iZERO, &iZERO, &ictxt, &m_loc_reduced, &info); // m_loc_reduced -- is local leading dimension
	//descinit_(descVector,  &N, &iONE, &NB, &iONE, &iZERO, &iZERO, &ictxt, &m_loc, &info);

	// print parameters of the processes
	for(int iter_rank = 0; iter_rank < world_size; iter_rank++)
		{
		if (iter_rank==world_rank) printf("load_cblacs.h: Proc-%d, m_loc= %d, n_loc=%d, m_loc_reduced= %d, n_loc_reduced=%d  \n", world_rank, m_loc, n_loc, m_loc_reduced, n_loc_reduced);
		//Cblacs_barrier(ictxt, "All");
		}

  return 0;
	}

int initialize_elpa(elpa_t *handle_ptr, int N, int nev_elpa){
    // Setup ELPA 
    //*
    int error;
        
    //elpa_set(handle, "na", n, &error); // matrix dimension
    elpa_set(*handle_ptr, "na", N, &error); // matrix dimension
    printf("na=%d is set, error=%d \n", N, error);

    elpa_set(*handle_ptr, "nev", nev_elpa, &error); // ?? number of eigenvalues/eigenvectors?
    printf("nev=%d is set, error=%d \n", nev_elpa, error);

    elpa_set(*handle_ptr, "local_nrows", m_loc, &error); // m_loc
    printf("local_nrows=%d is set, error=%d \n", m_loc, error);

    elpa_set(*handle_ptr, "local_ncols", n_loc, &error); // n_loc
    printf("local_ncols=%d is set, error=%d \n", n_loc, error);

    elpa_set(*handle_ptr, "nblk", NB, &error); // NB
    printf("nblk=%d is set, error=%d \n", NB, error);

    elpa_set(*handle_ptr, "mpi_comm_parent", MPI_Comm_c2f(MPI_COMM_WORLD), &error);
    printf("Fortran MPI communicator is set, error=%d \n", error);

    elpa_set(*handle_ptr, "process_row", myrow, &error); // myrow
    printf("process_row=%d is set, error=%d \n", myrow, error);
	
    elpa_set(*handle_ptr, "process_col", mycol, &error); // mycol
    printf("process_col=%d is set, error=%d \n", mycol, error);
    
    elpa_set(*handle_ptr, "timings", 1, &error);
    printf("timings=1 is set, error=%d \n", error);
   
    error = elpa_setup(*handle_ptr);
    printf("ELPA setup done, error=%d", error);

  //if(iam == 0) printf("ELPA setup done\n");
  // Set ELPA tunables 
    elpa_set(*handle_ptr, "solver", ELPA_SOLVER, &error);
    printf("elpa_set solver done, error=%d \n", error);
    if (ELPA_SOLVER == ELPA_SOLVER_1STAGE) printf("ELPA_SOLVER_1STAGE \n");
    if (ELPA_SOLVER == ELPA_SOLVER_2STAGE) printf("ELPA_SOLVER_2STAGE \n");
  
      
   elpa_set(*handle_ptr, "nvidia-gpu", 1, &error);
   printf("use nvidia-gpu is set, error=%d \n", error);
   
   if (ELPA_SOLVER == ELPA_SOLVER_2STAGE)
      {
      elpa_set(*handle_ptr, "real_kernel", ELPA_2STAGE_REAL_NVIDIA_GPU, &error);
      printf("setting real_kernel=ELPA_2STAGE_REAL_NVIDIA_GPU, error=%d \n", error);
      }

    return 0;
}   


//___________________________________________________________________________

double function_set_matrix (int I_gl, int J_gl)
   {
   double matrix_element=0;
   
   if (matrix_type==Symmetrized_Clement)
      {
      if (I_gl==J_gl+1 || J_gl==I_gl+1) 
         {
         int K = MIN(I_gl, J_gl);
         matrix_element = (double) sqrtl((K+1) * (N - K-1));
         }
      }
   else if (matrix_type==Sequence)
      {
      if (J_gl>I_gl) return function_set_matrix (J_gl, I_gl);
      return I_gl + N*J_gl + 1;
      }
   else if (matrix_type==Random_0_1)
      {
      matrix_element = (double) (rand())/RAND_MAX;
      //if (abs(I_gl-J_gl)!=1) matrix_element = double(rand())/RAND_MAX;
      //if (!(I_gl==0&&J_gl==1 || I_gl==1&&J_gl==0)) matrix_element = double(rand())/RAND_MAX;
      }
   else if (matrix_type==Random_m1_1)
      {
      matrix_element = 2*((double)(rand())/RAND_MAX - 0.5);
      }
   else
      {
      printf("matrix_type=%d is not supported", matrix_type);
      exit(1);
      }
   
   return matrix_element;
   }
//___________________________________________________________________________


int main(int argc, char** argv)
{

/************  MPI ***************************/
MPI_Init( &argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
MPI_Comm_size(MPI_COMM_WORLD, &world_size);
if (world_rank==0) printf("world_size=%d \n", world_size);

/************  --- ***************************/
// CUDA
int nDevices;
cudaGetDeviceCount(&nDevices);
printf("Number of Devices found: %d\n\n", nDevices);
//____________________________________________ 

if (argc==1) // one argument was provided: filename (default)
   {
   N = 10;
   nev = N;
   }
   
if (argc==2) // two arguments were provided: filename (default), N
   {
   N = (int) strtol(argv[1], NULL, 10);
   nev = N;
   }

if (N>100) debug_mode=0;

//srand(time(NULL));	
srand(1); // seeding
	
//----------------------------------------------------

// 0. Setup ELPA
LoadCblacs(N);
if (!((myrow>-1)&(mycol>-1)&(myrow<nprow)&(mycol<npcol))) // checks whether the process is in the grid.
	{
	printf("bye-bye from process-%d \n", world_rank); 
	MPI_Finalize();
	return 0;
	}

if (elpa_init(20170403) != ELPA_OK) 
  {
  printf("Error: ELPA API version not supported \n");
  return 1;
  }

int error;
elpa_t handle = elpa_allocate(&error);
  
initialize_elpa(&handle, N, nev);
 
// 1. Fill the local matrix
//vector<double> H_loc = vector<double>(m_loc*n_loc, 0);
double *H_loc;
H_loc  = (double *) calloc(m_loc*n_loc, sizeof(double));
	
for(int i_loc=0; i_loc<m_loc; i_loc++) // iteration over the first state
	{
	int l_1 = i_loc/NB; // local coord of the (NBxNB) block among other blocks
	int x_1 = i_loc%NB; // local coord within the block
	int I_gl= (l_1*nprow + myrow)*NB + x_1; // gl-global; nprow = "P_r"; myrow="p_r"  (quoted-ScaLAPACK userguide notation, p.88-90)
		
	for(int j_loc=0; j_loc<n_loc; j_loc++) // iteration over the second state
		{
      int l_2 = j_loc/NB; // local coord of the (NBxNB) block among other blocks
      int x_2 = j_loc%NB; // local coord within the block
      int J_gl= (l_2*npcol + mycol)*NB + x_2;

		H_loc[i_loc+j_loc*m_loc] = function_set_matrix(I_gl, J_gl);
		}
	} 
   
// 2. Allocating host memory for eigenvectors and eigenvalues

double *Eigvec_loc;
double *Eigval;
Eigvec_loc  = (double *) calloc(m_loc*n_loc, sizeof(double)); // must be always dimensioned to full size (N,N) even if only part of eigenvalues is needed [ELPA 2014 paper]
Eigval  = (double *) calloc(N, sizeof(double));

// analytic eigenvalues for testing
double *Eigval_analyt;
Eigval_analyt  = (double *) calloc(N, sizeof(double));
for (int i=0; i<N; i++) Eigval_analyt[i] = -N+1+2*i;

//____________ CUDA-part ________________

// 3. Allocating device memory and host-device copying
double *H_loc_dev, *Eigvec_loc_dev, *Eigval_dev;

cudaError_t error_cuda = cudaMalloc((void **) &H_loc_dev , m_loc*n_loc*sizeof(double));
if (error_cuda != cudaSuccess){    
   fprintf(stderr, "Error in cudaMalloc(H_loc_dev)\n");
   exit(1);
   }

error_cuda = cudaMalloc((void **) &Eigvec_loc_dev , m_loc*n_loc*sizeof(double));
if (error_cuda != cudaSuccess){    
   fprintf(stderr, "Error in cudaMalloc(Eigvec_loc_dev)\n");
   exit(1);
   }

error_cuda = cudaMalloc((void **) &Eigval_dev , N*sizeof(double));
if (error_cuda != cudaSuccess){    
   fprintf(stderr, "Error in cudaMalloc(Eigval_dev)\n");
   exit(1);
   }

// copy
error_cuda = cudaMemcpy(H_loc_dev, H_loc,  m_loc*n_loc*sizeof(double), cudaMemcpyHostToDevice);
if (error_cuda != cudaSuccess){   
   fprintf(stderr, "Error in cudaMemcpy(H_loc_dev, H_loc)\n");
   exit(1);
   }  
   
// 4. Diagonalization with ELPA

cudaDeviceSynchronize(); // for correct measurement of time
clock_t t0_elpa = clock();

elpa_eigenvectors_double(handle, H_loc_dev, Eigval_dev, Eigvec_loc_dev, &error);

cudaDeviceSynchronize();
clock_t t1_elpa = clock();

// 5. copying back from host to device

error_cuda = cudaMemcpy(Eigvec_loc, Eigvec_loc_dev, m_loc*n_loc*sizeof(double), cudaMemcpyDeviceToHost);
if (error_cuda != cudaSuccess){    
   fprintf(stderr, "Error in cudaMemcpy(Eigvec_loc, Eigvec_loc_dev)\n");
   exit(1);
   }

error_cuda = cudaMemcpy(Eigval, Eigval_dev, N*sizeof(double), cudaMemcpyDeviceToHost);
if (error_cuda != cudaSuccess){    
   fprintf(stderr, "Error in cudaMemcpy(H_loc, H_loc_dev)\n");
   exit(1);
   }

// 6. deallocating device pointers
error_cuda = cudaFree(H_loc_dev);
if (error_cuda != cudaSuccess){    
   fprintf(stderr, "Error in cudaFree(H_loc_dev)\n");
   exit(1);
   } 

error_cuda = cudaFree(Eigvec_loc_dev);
if (error_cuda != cudaSuccess){    
   fprintf(stderr, "Error in cudaFree(Eigvec_loc_dev)\n");
   exit(1);
   } 

error_cuda = cudaFree(Eigval_dev);
if (error_cuda != cudaSuccess){    
   fprintf(stderr, "Error in cudaFree(Eigval_dev)\n");
   exit(1);
   } 
   

if (world_rank==0) printf("ELPA diagonalization time: %f sec.\n", (double)(t1_elpa - t0_elpa)/CLOCKS_PER_SEC);

	
if (world_rank==0)
	{
	printf("ELPA diagonalization is done, Emin=%.16g\n", Eigval[0]);

	int N_eigval=MIN(10, N);
	
	printf("First %d eigenvalues: \n", N_eigval);
	for (int i=0; i<N_eigval; i++) printf("%.16g ", Eigval[i]);
	
	printf("\nAbsolute diff of eigenvalues wrt to analytic ones: \n");
	for (int i=0; i<N_eigval; i++) printf("%.16g ", Eigval[i] - Eigval_analyt[i]);
	}


free(H_loc);
free(Eigvec_loc);
free(Eigval);

if ((myrow>-1)&(mycol>-1)&(myrow<nprow)&(mycol<npcol)) Cblacs_gridexit(0); // leave grid if we were in one
  
MPI_Finalize();
}


