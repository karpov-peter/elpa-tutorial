#define assert_elpa_ok(error_code) call x_ao(error_code, stringify(error_code), __FILE__, __LINE__)

program eigenproblem_elpa
  use mpi
  ! Step 1: use the ELPA module
  use elpa
  implicit none
  
  ! Explicit interface for numroc
  interface 
    integer function numroc(N, NB, iproc, isrcproc, nprocs)
      integer, intent(in) :: N, NB, iproc, isrcproc, nprocs
    end function numroc
  end interface
  
  ! Step 2: define a handle for an ELPA object
  class(elpa_t), pointer :: elpaInstance
  class(elpa_autotune_t), pointer :: elpaTuneState
  integer :: status

  ! Declare global variables
  integer, parameter :: N_DEFAULT = 1000, NEV_DEFAULT = 500, NB_DEFAULT = 32
  integer :: N=N_DEFAULT, nev=NEV_DEFAULT, NB=NB_DEFAULT
  integer :: debug_mode = 0
  
  integer :: world_rank, world_size
  integer :: iam, nprocs, ictxt, nprow, npcol, myrow, mycol
  integer :: np, nq, m_loc, n_loc, info, itemp, seed
  integer :: descA(9), descZ(9)
  integer :: izero = 0, ione = 1
  integer :: i_arg, argc
  character(len=32) :: argv
  real(8) :: t0, t1, t2
  real(8), allocatable :: A_loc(:,:), A_loc_copy(:,:), Z_loc(:,:), Eigenvalues(:)
  integer :: I_loc, J_loc, I_global, J_global
  integer :: ierr
  integer :: l_1, x_1, I_gl
  integer :: l_2, x_2, J_gl
  integer :: iter, iter_max, N_eigenvalues_print
  character(len=5) :: iter_string
  logical :: unfinished

  ! ____________________________________________
  ! Set up MPI
  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, world_size, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, world_rank, ierr)
  
  ! ____________________________________________ 
  ! Process command-line arguments
  argc = command_argument_count()
  
  if (argc == 1) then
    call get_command_argument(1, argv)
    read(argv, *) N
    nev = N
    NB = min(N, NB_DEFAULT)
  end if
  
  if (argc == 2) then
    call get_command_argument(1, argv)
    read(argv, *) N
    call get_command_argument(2, argv)
    read(argv, *) nev
    NB = min(NB_DEFAULT, N)
  end if
  
  if (argc >= 3) then
    call get_command_argument(1, argv)
    read(argv, *) N
    call get_command_argument(2, argv)
    read(argv, *) nev
    call get_command_argument(3, argv)
    read(argv, *) NB
  end if

  if (N > 100) debug_mode = 0

  ! ____________________________________________ 
  ! determine the blacs grid size
  nprow = nint(sqrt(real(world_size)))
  npcol = world_size / nprow
  
  ! if the grid is not square, try to find the "most square" rectangular grid
  do while (nprow > 0)
    if (mod(world_size, nprow) == 0) then
        npcol = world_size / nprow
        exit
    end if
    nprow = nprow - 1
  end do
  
  ! try some other grids if you want
  ! if (world_size==72) {nprow=12; npcol=6;}

  ! specially treat a case of small matrices
  if (NB > N/max(nprow, npcol)) NB = N/max(nprow, npcol)

  if (nprow*npcol /= world_size) then
    if (world_rank == 0) then
        print *, 'Error:  wrong grid size'
    end if
    call MPI_Finalize(ierr) 
    stop 1
  end if

  print *, "N=", N, ", nev=", nev, ", NB=", NB

  ! ____________________________________________ 
  ! Set up blacs grid
  call blacs_pinfo(iam, nprocs)
  call blacs_get(0, 0, ictxt)
  call blacs_gridinit(ictxt, 'row', nprow, npcol)
  call blacs_gridinfo(ictxt, nprow, npcol, myrow, mycol)
  
  ! Compute the size of the local matrices using numroc
  m_loc = numroc(N, NB, myrow, 0, nprow)
  n_loc = numroc(N, NB, mycol, 0, npcol)
  if (debug_mode == 0) write(*,*) 'myrow=', myrow, ', mycol=', mycol, ', m_loc=', m_loc, ', n_loc=', n_loc

  ! Initialize array descriptors for distributed matrices A and Z
  itemp = max( 1, m_loc );
  print *, "m_loc=", m_loc, ", n_loc=", n_loc, ", itemp=", itemp
  call descinit(descA, N, N, NB, NB, 0, 0, ictxt, itemp, info)
  call descinit(descZ, N, N, NB, NB, 0, 0, ictxt, itemp, info)

  ! ____________________________________________ 
  ! Setup ELPA 
  
  ! Step 3: initialize the ELPA library
  status = elpa_init(20231705)
  if (status /= ELPA_OK) then
      print *, "ELPA API version not supported"
    stop 1
  endif
  
  ! Step 4: allocate the ELPA object
  elpaInstance => elpa_allocate(status)
  ! Check status code, e.g. with
  if (status /= ELPA_OK) then
    print *, "Could not allocate ELPA instance"
    stop 1
  endif

  ! Step 5: set mandatory parameters describing the matrix and its MPI distribution
  ! This can be done only once per ELPA object (handle)
  call elpaInstance%set("na", N, status)
  if (status /= ELPA_OK) then
    print *, "Could not set parameter na"
    ! Handle this error in your application
  endif
  
  call elpaInstance%set("nev", nev, status)
  ! Check status code ...
  
  call elpaInstance%set("local_nrows", m_loc, status)
  ! Check status code ...
  
  call elpaInstance%set("local_ncols", n_loc, status)
  ! Check status code ...
  
  call elpaInstance%set("nblk", NB, status)
  ! Check status code ...
  
  call elpaInstance%set("mpi_comm_parent", MPI_COMM_WORLD, status)
  ! Check status code ...
  
  call elpaInstance%set("process_row", myrow, status)
  ! Check status code ...
  
  call elpaInstance%set("process_col", mycol, status)
  ! Check status code ...

  ! Step 6: set up the elpa object, finalize setting of mandatory parameters
  status = elpaInstance%setup()
  if (status /= ELPA_OK) then
    print *, "Could not setup the ELPA object"
    ! Handle this error in your application
  endif

  ! Step 7: set runtime options, e.g. GPU settings
  call elpaInstance%set("solver", ELPA_SOLVER_2STAGE, status)
  ! Check status code ...

  ! End of Setup ELPA 
  ! ____________________________________________ 
  ! Allocate the matrices A_loc, Z_loc (matrix of eigenvectors), and vector of Eigenvalues
  allocate(A_loc(m_loc, n_loc)) ! matrix to be diagonalized, local
  allocate(A_loc_copy(m_loc, n_loc))
  allocate(Z_loc(m_loc, n_loc)) ! matrix of eigenvectors,    local
  allocate(Eigenvalues(N))     ! vector of eigenvalues,     global

  ! Fill the local matrix A_loc 
  do i_loc = 1, m_loc
    l_1 = (i_loc-1)/NB     ! local coord of the (NBxNB) block among other blocks
    x_1 = mod(i_loc-1, NB) ! local coord within the block
    I_gl = (l_1*nprow + myrow)*NB + x_1 + 1 ! gl-global; nprow = "P_r"; myrow="p_r"  (quoted-ScaLAPACK userguide notation, p.88-90)
  
    do j_loc = 1, n_loc
      l_2 = (j_loc-1)/NB
      x_2 = mod(j_loc-1, NB)
      J_gl = (l_2*npcol + mycol)*NB + x_2 + 1
  
      A_loc(i_loc, j_loc) = get_global_matrix_element(I_gl, J_gl, N)
      A_loc_copy(i_loc, j_loc) = get_global_matrix_element(I_gl, J_gl, N)
    end do
  end do
  
  ! ____________________________________________ 
  ! print the local part of matrix A
  if (debug_mode == 1 .and. world_rank == 0) then
    write(*,*) "A_loc:"
    do i_loc = 1, m_loc
       do j_loc = 1, n_loc
          write(*,'(F6.2,X)', advance='no'), A_loc(i_loc, j_loc)
       end do
       write(*,*)
    end do
    write(*,*)
  end if

  ! ____________________________________________ 
  ! Setup autotune

  ! if you want to use the new autotuning implentation
  !call elpaInstance%autotune_set_api_version(20211125, status) ! THIS CAUSES ERROR!
  call elpaInstance%autotune_set_api_version(20170403, status) ! THIS CAUSES ERROR!
  
  call assert(status == ELPA_OK)
  ! if you want to use the old one, either do not set autotune_set_api_version
  ! or set autotune_set_api_version to a supported api version < 20211125
  elpaTuneState => elpaInstance%autotune_setup(ELPA_AUTOTUNE_MEDIUM, ELPA_AUTOTUNE_DOMAIN_REAL, status)
  call assert(status == ELPA_OK)

  !____________________________________________ 
  ! Print ELPA settings
  call elpaInstance%print_settings(status)

  ! ____________________________________________ 
  ! Step 8: Perform the diagonalization using ELPA
  iter_max=100

  do iter = 1, iter_max
    unfinished = elpaInstance%autotune_step(elpaTuneState, status)
    call assert(status == ELPA_OK)

    if (unfinished == .false.) then
      if (world_rank == 0) then
          print *, 'ELPA autotuning finished in ', iter-1, ' steps'
      end if
      exit   
    end if
  
    ! Solve EV problem 
    call elpaInstance%eigenvectors(A_loc, Eigenvalues, Z_loc, status)
    call assert(status == ELPA_OK)
    
    ! restore the matrix
    A_loc(:,:) = A_loc_copy(:,:)

    if (world_rank == 0) then
      print *, ""
      print *, "iter ", iter
      call elpaInstance%autotune_print_state(elpaTuneState)
    endif

  end do

  ! set and print the autotuned-settings
  call elpaInstance%autotune_set_best(elpaTuneState, status)
  call assert(status == ELPA_OK)

  if (world_rank .eq. 0) then
    print *, ""
    print *, "The best combination found by the autotuning:"
    call elpaInstance%autotune_print_best(elpaTuneState, status)
    call assert(status == ELPA_OK)
  endif

  ! ____________________________________________ 
  ! Print the results
  if (world_rank == 0) then
    write(*,*) "Autotuning is done"
  end if

  ! ____________________________________________ 
  ! Step 9: Clean up ELPA

  call elpa_deallocate(elpaInstance, status)
  ! Check status code ...
  call elpa_uninit()

  ! ____________________________________________ 
  ! Clean up the rest and finalize
  deallocate(A_loc) 
  deallocate(A_loc_copy) 
  deallocate(Z_loc)
  deallocate(Eigenvalues)
  
  call blacs_gridexit(ictxt)
  call MPI_Finalize(ierr)
  
contains

! Function to compute the global matrix element of "Symmetrized_Clement" matrix
function get_global_matrix_element(I_gl, J_gl, n) result(matrix_element)
  implicit none
  integer, intent(in) :: I_gl, J_gl, n
  real(8) :: matrix_element
  integer :: k

  matrix_element = 0.0d0
  if (I_gl == J_gl+1 .or. J_gl == I_gl+1) then
    K = min(I_gl, J_gl)
    matrix_element = dsqrt(dble(k * (n - k)))
  end if
end function get_global_matrix_element

subroutine assert(condition)
  logical, intent(in) :: condition
  if (.not. condition) then
    stop 1
  endif
end subroutine assert

end program