#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "mpi.h"

#include <elpa/elpa.h>
#include <assert.h>

#define assert_elpa_ok(x) assert(x == ELPA_OK)

static int max( int a, int b ){if (a>b) return(a); else return(b);}
static int min( int a, int b ){if (a<b) return(a); else return(b);}

//__________________________________________________________________________________

int N = 1000;
int nev=64;
int NB = 32; 

int debug_mode=0;

enum matrix_types {Symmetrized_Clement, Random_0_1};
enum matrix_types matrix_type = Symmetrized_Clement;

int ELPA_SOLVER = ELPA_SOLVER_2STAGE;

//__________________________________________________________________________________
// blacs
void   Cblacs_pinfo( int* mypnum, int* world_size);
void   Cblacs_get( int context, int request, int* value);
int    Cblacs_gridinit( int* context, char * order, int np_row, int np_col);
void   Cblacs_gridinfo( int context, int*  np_row, int* np_col, int*  my_row, int*  my_col);
void   Cblacs_gridexit( int context);
void   Cblacs_exit( int error_code);

// scalapack
int    numroc_( int *n, int *NB, int *iproc, int *isrcproc, int *world_size);
void   descinit_( int *desc, int *m, int *n, int *mb, int *NB, int *irsrc, int *icsrc, int *ictxt, int *lld, int *info);

//__________________________________________________________________________________


int Setup_ELPA()
	{
	return 0;	
	}

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

// for hybrid MPI+OpenMP
//int thread_level;
//MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &thread_level); 

MPI_Init( &argc, &argv);
MPI_Comm_size(MPI_COMM_WORLD, &world_size);
MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

//____________________________________________ 

//int world_rank, world_size; // blacs
int ictxt, nprow, npcol, myrow, mycol;
int np, nq;
int m_loc, n_loc;
int info, itemp, seed, lwork;
int descA[9], descZ[9];
int izero=0;
//double *A_loc; // this is the matrix we diagonalize here
//double *Z_loc, *work, *Eigenvalues;
  


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

nprow = sqrt(world_size); 
npcol = sqrt(world_size); 
if (world_size==2)  {nprow=2; npcol=1;}
if (world_size==6)  {nprow=3; npcol=2;}

if (world_size==10) {nprow=2; npcol=5;}
if (world_size==20) {nprow=4; npcol=5;}
if (world_size==40) {nprow=5; npcol=8;}

if (world_size==72) {nprow=8; npcol=9;}
//if (world_size==72) {nprow=9; npcol=8;}
//if (world_size==72) {nprow=6; npcol=12;}
//if (world_size==72) {nprow=12; npcol=6;}

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
	
// setup Blacs
Cblacs_pinfo( &world_rank, &world_size ) ;
Cblacs_get( -1, 0, &ictxt );
Cblacs_gridinit( &ictxt, "Row", nprow, npcol );
Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );

// Compute the size of the local matrices (thanks to numroc)
m_loc = numroc_( &N, &NB, &myrow, &izero, &nprow );
n_loc = numroc_( &N, &NB, &mycol, &izero, &npcol );
if (debug_mode) printf("myrow=%i, mycol=%i, m_loc=%i, n_loc=%i, \n", myrow, mycol, m_loc, n_loc);

//____________________________________________ 

/* Setup ELPA */
///*
int error, status;
elpa_t handle;
   
if (elpa_init(20170403) != ELPA_OK) 
	{
	fprintf(stderr, "Error: ELPA API version not supported");
	exit(1);
	}

handle = elpa_allocate(&error);
assert_elpa_ok(error);

//elpa_set(handle, "na", n, &error); // matrix dimension
elpa_set_integer(handle, "na", N, &error); // matrix dimension
if (world_rank==0 || debug_mode) printf("na=%i is set, error=%i\n", N, error);

elpa_set_integer(handle, "nev", nev, &error); // ?? number of eigenvalues/eigenvectors?
if (world_rank==0 || debug_mode) printf("nev=%i is set, error=%i\n", nev, error);

elpa_set_integer(handle, "local_nrows", m_loc, &error); // m_loc
if (world_rank==0 || debug_mode) printf("local_nrows=%i is set, error=%i\n", m_loc, error);

elpa_set_integer(handle, "local_ncols", n_loc, &error); // n_loc
if (world_rank==0 || debug_mode) printf("local_ncols=%i is set, error=%i\n", n_loc, error);

elpa_set_integer(handle, "nblk", NB, &error); // NB
if (world_rank==0 || debug_mode) printf("nblk=%i is set, error=%i\n", NB, error);

elpa_set_integer(handle, "mpi_comm_parent", MPI_Comm_c2f(MPI_COMM_WORLD), &error);
if (world_rank==0 || debug_mode) printf("Fortran MPI communicator is set, error=%i\n", error);

elpa_set_integer(handle, "process_row", myrow, &error); // myrow
if (world_rank==0 || debug_mode) printf("process_row=%i is set, error=%i\n", myrow, error);

elpa_set_integer(handle, "process_col", mycol, &error); // mycol
if (world_rank==0 || debug_mode) printf("process_col=%i is set, error=%i\n", mycol, error);

int error1 = elpa_setup(handle);
printf("ELPA setup done, error=%i\n", error1);
if (error1!=ELPA_OK)
  {
  printf("ELPA setup failed, error=%i\n", error1);
  MPI_Finalize();
  return(1);
  }
//if(world_rank == 0) printf("ELPA setup done\n");
// Set ELPA tunables 
elpa_set_integer(handle, "solver", ELPA_SOLVER, &error);
if (world_rank==0 || debug_mode) 
  {
  printf("elpa_set solver done, error=%d \n", error);
  if (ELPA_SOLVER == ELPA_SOLVER_1STAGE) printf("ELPA_SOLVER_1STAGE \n");
  if (ELPA_SOLVER == ELPA_SOLVER_2STAGE) printf("ELPA_SOLVER_2STAGE \n");
  }
  
MPI_Barrier(MPI_COMM_WORLD);   
// End Setup ELPA 


// parameters of the matrix


 
// Work only the process in the process grid

double *A_loc, *Z_loc, *work, *Eigenvalues;
//vector<double> A_loc; // this is the matrix we diagonalize here
//vector<double> Z_loc, work, Eigenvalues;


// Allocate the matrices A_loc, Z_loc(vector of eigenvalues), Eigenvalues
A_loc = (double *)calloc(m_loc*n_loc,sizeof(double)) ;
Z_loc = (double *)calloc(m_loc*n_loc,sizeof(double)) ;
Eigenvalues = (double *)calloc(N,sizeof(double)) ;
  

// fill the local part of matrix A 
int k = 0;
int i_loc, j_loc;
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

        A_loc[i_loc+j_loc*m_loc] = get_global_matrix_element(I_gl, J_gl);
        }
    }

      
// Initialize the array descriptor for the distributed matrices xA, U and VT
itemp = max( 1, m_loc );
descinit_( descA,  &N, &N, &NB, &NB, &izero, &izero, &ictxt, &itemp, &info );
descinit_( descZ,  &N, &N, &NB, &NB, &izero, &izero, &ictxt, &itemp, &info );



double t_start = MPI_Wtime();

elpa_eigenvectors(handle, A_loc, Eigenvalues, Z_loc, &error); // Eigenvalues=eigenvector matrix, Z_loc =eigenvalue vector

double t_stop = MPI_Wtime();
	

	

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

// clean up

free(A_loc);
free(Z_loc);
free(Eigenvalues);

Cblacs_gridexit(0);

elpa_deallocate(handle, &error);
elpa_uninit(&error);
		
MPI_Finalize();

return(0);
}
