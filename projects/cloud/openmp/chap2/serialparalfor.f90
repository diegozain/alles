program main
use omp_lib
implicit none
! ------------------------------------------------------------------------------
! parallelize 2 for loops. first one sequentially, second in parallel.
!
! gfortran -fopenmp file.f90
! ------------------------------------------------------------------------------
! ⬛⬜ variables
integer, parameter :: max_threads = 4
integer ::  i, id, numthreads, nthread
integer, parameter :: num_steps = 100000000

! ⌚
real*8 :: start_time, run_time

! for-🍭 variables
integer, parameter :: nfreq=5
! with nab=16000000 there is an
! order of magnitude difference between serial & parallel.
integer, parameter :: nab=16
integer :: ifreq, iab
! ------------------------------------------------------------------------------
! this call can be avoided by declaring a bash variable like so:
! export OMP_NUM_THREADS=4
! this means the 'max_threads' variable can be avoided and
! re-compiling when changing the num of threads too.
call omp_set_num_threads(max_threads)
! ------------------------------------------------------------------------------
! ⌚
start_time = omp_get_wtime()
! ------------------------------------------------------------------------------
!                            ⬛⬜⬛⬜
! ------------------------------------------------------------------------------
do ifreq=1,nfreq
  !$omp parallel
  !$omp do
  do iab=1,nab
    nthread = omp_get_thread_num()
    ! ✏️
    write(*,101) ifreq, iab, nthread
    101 format('i𝙵 ', i2, ' . i𝒜ℬ ', i2, ' •', i2)
  enddo
  !$omp end do
  !$omp end parallel
  print*,''
enddo
! ------------------------------------------------------------------------------
! ⌚
run_time = omp_get_wtime() - start_time
! ------------------------------------------------------------------------------
! ✏️
write(*,102) run_time
102 format('run time ', f15.8)
! ------------------------------------------------------------------------------
end program main
