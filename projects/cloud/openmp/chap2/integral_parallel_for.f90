program main
use omp_lib
implicit none
! ------------------------------------------------------------------------------
! parallelize the numerical integral of 4/(1+x^2) between 0 and 1.
! this integral is equal to pi.
! 
! this approach is the cyclic parallelization of the integral for loop.
!
! it implements 'critical barriers' that solve the problem of 'false sharing'. 
! ------------------------------------------------------------------------------

! parallel variables
integer, parameter :: max_threads=4
integer :: i
integer, parameter :: num_steps=100000000

! integral variables
real*8 :: pi, step, full_sum, x
real*8 :: start_time, run_time
real*8 :: partial_sum

full_sum = 0.0
step = 1.0 / num_steps

call omp_set_num_threads(max_threads)
start_time = omp_get_wtime()
! ------------------------------------------------------------------------------
!$omp parallel private(i,partial_sum,x)

partial_sum = 0.0

! the construct '$omp do ... $omp end do' takes care of all the
! ID/num of threads/cyclic or block/ stuff, and has the compiler (via ompenMP)
! do it for you.

!$omp do
do i = 1, num_steps
   x = (i+0.5)*step
   partial_sum = partial_sum + 4.0/(1.0+x*x)
enddo
!$omp end do
! ------------------------------------------------------------------------------
! the construct '$omp do ... $omp end do' also takes care of the
! critical barrier.

full_sum = full_sum + partial_sum
! ------------------------------------------------------------------------------
!$omp end parallel
! ------------------------------------------------------------------------------
! this just scales the full sum
pi = step * full_sum

run_time = omp_get_wtime() - start_time

write(*,100) pi, run_time
100     format('pi is ',f15.8,' in ',f8.3,'secs')

end program main
