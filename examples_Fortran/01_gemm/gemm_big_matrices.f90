program gemm_big_marices
  use omp_lib, only: omp_get_wtime
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  integer :: i, j, print_size, argc
  integer(kind=8) :: N, k
  real(dp), allocatable :: A(:), B(:), C(:)
  real(dp) :: t1, t2, alpha, beta
  character(len=32) :: argv
  integer :: ierr

  ! Get command line arguments
  argc = COMMAND_ARGUMENT_COUNT()
  if (argc /= 1) then
      print *, 'Usage: <executable_name> <matrix_size>'
      stop
  end if
  
  call GET_COMMAND_ARGUMENT(1, argv, status=ierr)
  
  ! Convert command line argument to integer
  if (ierr == 0) then
      read(argv, *) N
  else
      print *, "Error reading command line argument."
      stop
  end if

  ! Allocate matrices
  allocate(A(N*N), B(N*N), C(N*N))

  ! Initialize matrices
  do j = 1, N
      do i = 1, N
          A(i + (j-1)*N) = dble(i)
          B(i + (j-1)*N) = dble(j)
      end do
  end do

  alpha = 1.0_dp / dble(N)
  beta = 0.0_dp

  t1 = omp_get_wtime()

  ! Matrix multiplication using dgemm
  call dgemm('N', 'N', N, N, N, alpha, A, N, B, N, beta, C, N)

  t2 = omp_get_wtime()

  ! Print matrices
  print_size = 10
  print *, "Upper-left ", print_size, "x", print_size, " corner of the resulting C matrix:"
  do i = 1, print_size
      do j = 1, print_size
          write(*,'(e12.4)', advance="no") C(i+(j-1)*N)
      end do
      print *   ! Move to the next line after each row
  end do
  print *  ! Empty line

  ! print *, 'GEMM time (sec): ', t2-t1
  write(*,'("GEMM time (sec): ",f8.3)') t2-t1
  

  ! Deallocate matrices
  deallocate(A, B, C)

end program matrix_mult
