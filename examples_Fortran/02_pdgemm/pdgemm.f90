program pdgemm_example
  use mpi

  implicit none

  ! Explicit interface for numroc
  interface 
    integer function numroc(N, NB, iproc, isrcproc, nprocs)
      integer, intent(in) :: N, NB, iproc, isrcproc, nprocs
    end function numroc
  end interface

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: N_default = 1000, NB_default = 32
  integer :: N = n_default, NB = nb_default
  integer :: world_rank, world_size, iam, nprocs
  integer :: ictxt, nprow, npcol, myrow, mycol
  integer :: m_loc, n_loc, info, itemp, izero = 0, ione = 1
  integer :: descA(9), descB(9), descC(9), descGlobal(9)
  integer :: i, j, print_size, ierror
  real(dp) :: alpha, beta = 0.0_dp
  real(dp), allocatable :: A(:), B(:), C(:)
  real(dp), allocatable :: A_loc(:), B_loc(:), C_loc(:)
  character(len=20) :: arg_str   ! temporary string to hold the argument
  integer :: I_gl, J_gl

  ! Initialize MPI
  call mpi_init(info)
  call mpi_comm_size(mpi_comm_world, world_size, info)
  call mpi_comm_rank(mpi_comm_world, world_rank, info)

  !  Read command line arguments (assuming arguments are passed as N, NB)
  if (command_argument_count() >= 1) then
    call get_command_argument(1, arg_str)
    read(arg_str, *) N
    if (command_argument_count() == 2) then
      call get_command_argument(2, arg_str)
      read(arg_str, *) NB
    end if
  end if

  ! Compute grid size
  nprow = int(sqrt(real(world_size)))
  npcol = world_size / nprow
  if (NB > N / max(nprow, npcol)) NB = N / max(nprow, npcol)

  ! Setup blacs grid
  call blacs_get(-1, 0, ictxt)
  call blacs_gridinit(ictxt, 'Row', nprow, npcol)
  call blacs_gridinfo(ictxt, nprow, npcol, myrow, mycol)

  ! Compute the size of the local matrices (thanks to numroc)
  m_loc = numroc(N, NB, myrow, izero, nprow)
  n_loc = numroc(N, NB, mycol, izero, npcol)

  ! Initialize the array descriptor for the distributed matrices A/A_loc, B/B_loc, C/C_loc
  itemp = max(1, m_loc)
  call descinit(descA, N, N, NB, NB, izero, izero, ictxt, itemp, info)
  call descinit(descB, N, N, NB, NB, izero, izero, ictxt, itemp, info)
  call descinit(descC, N, N, NB, NB, izero, izero, ictxt, itemp, info)

  ! Allocate memory for the local matrices
  allocate(A_loc(m_loc * n_loc))
  allocate(B_loc(m_loc * n_loc))
  allocate(C_loc(m_loc * n_loc))
  

  !____________________________________________ 

  ! Fill global matrices A and B. This is not a very efficient way to do it, 
  ! it's better to initialize the local matrices in parallel -- this is covered in the next example
  ! For saving memory, we will fill the global matrices only on 0-th MPI rank

  if (world_rank == 0) then
      allocate(A(N*N))
      allocate(B(N*N))

      do J_gl = 1, N
          do I_gl = 1, N
              A(I_gl+ (J_gl-1)*N) = dble(I_gl)
              B(I_gl+ (J_gl-1)*N) = dble(J_gl)
          end do
      end do
  end if

  ! create a desccriptor for the global matrix
  call descinit( descGlobal, N, N, N, N, izero, izero, ictxt, N,  info );
  
  ! redestribute the global matrices to the local ones
  call pdgemr2d(N, N, A,     ione, ione, descGlobal, &
                      A_loc, ione, ione, descA,  ictxt);
  call pdgemr2d(N, N, B,     ione, ione, descGlobal, &
                      B_loc, ione, ione, descB,  ictxt);
  
  ! we don't need the global matrices A and B anymore and can deallocate them to save memory
  deallocate(A);
  deallocate(B);
  !____________________________________________ 

  ! Call pdgemm
  alpha = 1.0_dp/dble(N)
  beta = 0.0_dp
  call pdgemm('N', 'N', N, N, N, alpha, A_loc, ione, ione, desca, B_loc, ione, ione, descb, beta, C_loc, ione, ione, descc)


  !____________________________________________ 
  ! redistribute the local matrix C_loc to the global matrix C
  allocate(C(N*N))

  call pdgemr2d(N, N, C_loc, ione, ione, descC, &
                  C,     ione, ione, descGlobal, ictxt);

  ! Print matrix C
  print_size = min(10,N)
  if (myrow == 0 .and. mycol == 0) then
    print *, "Upper-left ", print_size, "x", print_size, " corner of the resulting C matrix:"
    do i = 1, print_size
        do j = 1, print_size
            write(*,'(e12.4)', advance="no") C(i+(j-1)*N)
        end do
        print *   ! Move to the next line after each row
    end do
    print *  ! Empty line
  end if

  ! Cleanup
  deallocate(A_loc)
  deallocate(B_loc)
  deallocate(C_loc)
  deallocate(C)
  call blacs_gridexit(ictxt)
  call mpi_finalize(ierror)

end program
