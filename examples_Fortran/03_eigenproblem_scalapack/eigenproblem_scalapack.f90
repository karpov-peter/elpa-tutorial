! Include MPI and BLACS modules
program eigenproblem_scalapack
  use mpi
  !use blacs_f77
  implicit none
  
  ! Declare global variables
  integer, parameter :: N_DEFAULT = 1000, NEV_DEFAULT = 64, NB_DEFAULT = 32
  integer :: N = N_DEFAULT, nev = NEV_DEFAULT, NB = NB_DEFAULT
  integer :: debug_mode = 0
  integer, parameter :: scalapack_pdsyev = 1, scalapack_pdsyevd = 2, scalapack_pdsyevr = 3, scalapack_pdsyevx = 4
  integer :: diagonalization_method = scalapack_pdsyev
  
  integer :: world_rank, world_size
  integer :: iam, nprocs, ictxt, nprow, npcol, myrow, mycol
  integer :: np, nq, m_loc, n_loc, info, itemp, seed
  integer :: descA(9), descZ(9)
  integer :: izero = 0, ione = 1
  integer :: i_arg, argc
  character(len=32) :: argv
  real(8) :: t0, t1, t2
  real(8), allocatable :: A_loc(:,:), Z_loc(:,:), Eigenvalues(:)
  integer :: I_loc, J_loc, I_global, J_global
  real(8), allocatable :: WORK(:)
  integer :: LWORK
  integer :: ierr
  integer :: l_1, x_1, I_gl
  integer :: l_2, x_2, J_gl
  ! Assuming the size of A_loc is m_loc by n_loc


  ! Initialize MPI
  call MPI_Init(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, world_size, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, world_rank, ierr)
  
  ! Process command-line arguments
  argc = command_argument_count()
  
  if (argc == 1) then
    N = 10
    nev = N
  end if
  
  if (argc == 2) then
    call get_command_argument(1, argv)
    read(argv, *) N
    nev = N
  end if
  
  if (argc == 3) then
    call get_command_argument(1, argv)
    read(argv, *) N
    call get_command_argument(2, argv)
    read(argv, *) nev
  end if
  
  if (argc == 4) then
    call get_command_argument(1, argv)
    read(argv, *) N
    call get_command_argument(2, argv)
    read(argv, *) nev
    call get_command_argument(3, argv)
    read(argv, *) NB
  end if
  
  if (N > 100) debug_mode = 0
  
  ! determine the grid size
  nprow = nint(sqrt(real(world_size)))
  nprow = nint(sqrt(real(world_size)))
  
  ! ... More of the grid size calculation ... 
! BLACS initialization
  call blacs_pinfo(iam, nprocs)
  call blacs_get(0, 0, ictxt)
  call blacs_gridinit(ictxt, 'row', nprow, npcol)
  call blacs_gridinfo(ictxt, nprow, npcol, myrow, mycol)
  
  ! Compute local size
  CALL NUMROC(np, N, myrow, 0, nprow, NB)
  CALL NUMROC(nq, N, mycol, 0, npcol, NB)
  
  m_loc = np
  n_loc = nq
  
  ! Initialize array descriptors for distributed matrix A_loc and Z_loc
  CALL DESCINIT(descA, N, N, NB, NB, 0, 0, ictxt, m_loc, info)
  CALL DESCINIT(descZ, N, N, NB, NB, 0, 0, ictxt, m_loc, info)
  
  ! Allocate local storage for distributed matrix A_loc and Z_loc
  allocate(A_loc(m_loc, n_loc))
  allocate(Z_loc(m_loc, n_loc))
  
  ! Fill in the local pieces of A_loc
  do i_loc = 1, m_loc
    l_1 = (i_loc-1)/NB
    x_1 = mod(i_loc-1, NB)
    I_gl = (l_1*nprow + myrow)*NB + x_1 + 1
  
    do j_loc = 1, n_loc
      l_2 = (j_loc-1)/NB
      x_2 = mod(j_loc-1, NB)
      J_gl = (l_2*npcol + mycol)*NB + x_2 + 1
  
      A_loc(i_loc, j_loc) = get_global_matrix_element(I_gl, J_gl, N)
    end do
  end do
  
  ! Now, compute the eigenvalues using ScaLAPACK, depending on the chosen method
  allocate(Eigenvalues(N))
  
  ! Other necessary workspace and lwork computation will be needed here. 
  ! This is a general placeholder and might need adjustments based on the chosen method.
  LWORK = m_loc * n_loc
  allocate(WORK(LWORK))
  
  if (diagonalization_method == scalapack_pdsyev) then
    CALL PDSYEV('V', 'U', N, A_loc, ione, ione, descA, Eigenvalues, Z_loc, ione, ione, descZ, WORK, LWORK, info)
  elseif (diagonalization_method == scalapack_pdsyevd) then
    ! ... PDSYEVD call and necessary adjustments ...
  elseif (diagonalization_method == scalapack_pdsyevr) then
    ! ... PDSYEVR call and necessary adjustments ...
  elseif (diagonalization_method == scalapack_pdsyevx) then
    ! ... PDSYEVX call and necessary adjustments ...
  else
    print *, "Unknown diagonalization method!"
    CALL EXIT(1)
  end if
  
  ! You may want to print the eigenvalues or any other results here
  
  ! Clean up and finalize
  deallocate(A_loc, Z_loc, Eigenvalues, WORK)
  CALL BLACS_GRIDEXIT(ictxt)
  CALL MPI_FINALIZE(ierr)
  

    
contains

  ! Function to compute the global matrix element of "Symmetrized_Clement" matrix
function get_global_matrix_element(i_gl, j_gl, n) result(matrix_element)
  implicit none
  integer, intent(in) :: i_gl, j_gl, n
  real(8) :: matrix_element
  integer :: k

  matrix_element = 0.0d0
  if (i_gl == j_gl + 1 .or. j_gl == i_gl + 1) then
      k = min(i_gl, j_gl)
      matrix_element = sqrt(real((k + 1) * (n - k - 1)))
  end if

end function get_global_matrix_element

!   ! Function to compute the global matrix element
! real(8) function get_global_matrix_element(I_gl, J_gl)
! integer, intent(in) :: I_gl, J_gl
! integer :: K

! get_global_matrix_element = 0.0d0

! if (matrix_type == Symmetrized_Clement) then
!   if (I_gl == J_gl+1 .or. J_gl == I_gl+1) then
!     K = min(I_gl, J_gl)
!     get_global_matrix_element = sqrt(dble(K+1) * dble(N - K-1))
!   end if
! else if (matrix_type == Random_0_1) then
!   get_global_matrix_element = dble(rand())/dble(RAND_MAX)
! else
!   print *, "matrix_type=", matrix_type, " is not supported"
!   call exit(1)
! end if

! end function get_global_matrix_element

end program