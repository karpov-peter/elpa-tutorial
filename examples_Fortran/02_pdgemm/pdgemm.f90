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
  integer :: desca(9), descb(9), descc(9)
  integer :: i, j, print_size, ierror
  real(dp) :: alpha = 1.0_dp, beta = 0.0_dp
  real(dp), allocatable :: A(:), B(:), C(:)
  character(len=20) :: arg_str   ! temporary string to hold the argument
  

  ! Initialize MPI
  call mpi_init(info)
  call mpi_comm_size(mpi_comm_world, world_size, info)
  call mpi_comm_rank(mpi_comm_world, world_rank, info)

  ! Argument handling (assuming arguments are passed as N, NB)
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

  ! Setup blacs
  call blacs_get(-1, 0, ictxt)
  call blacs_gridinit(ictxt, 'Row', nprow, npcol)
  call blacs_gridinfo(ictxt, nprow, npcol, myrow, mycol)

  ! Compute local matrix size
  m_loc = numroc(N, NB, myrow, izero, nprow)
  n_loc = numroc(N, NB, mycol, izero, npcol)

  ! Initialize descriptors
  itemp = max(1, m_loc)
  call descinit(desca, N, N, NB, NB, izero, izero, ictxt, itemp, info)
  call descinit(descb, N, N, NB, NB, izero, izero, ictxt, itemp, info)
  call descinit(descc, N, N, NB, NB, izero, izero, ictxt, itemp, info)

  ! Allocate and initialize matrices
  allocate(a(m_loc * n_loc), b(m_loc * n_loc), c(m_loc * n_loc))
  A = 1.0_dp
  B = 1.0_dp
  C = 0.0_dp

  ! Call pdgemm
  call pdgemm('N', 'N', N, N, N, alpha, a, ione, ione, desca, b, ione, ione, descb, beta, c, ione, ione, descc)


  ! Print matrix C
  print_size = 10
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
  deallocate(a, b, c)
  call blacs_gridexit(ictxt)
  call mpi_finalize(ierror)

end program
