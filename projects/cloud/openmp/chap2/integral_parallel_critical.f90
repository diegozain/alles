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
integer :: i, j, id, numthreads, nthreads
integer, parameter :: num_steps=100000000
integer, parameter :: max_threads=4

! integral variables
real*8 :: pi, real_sum, step, full_sum, x
real*8 :: start_time, run_time
real*8 :: sum(0:max_threads-1)
real*8 :: partial_sum

full_sum = 0.0
step = 1.0 / num_steps

call omp_set_num_threads(max_threads)
full_sum = 0.0
start_time = omp_get_wtime()
! ------------------------------------------------------------------------------
!$omp parallel private(i,id,numthreads,partial_sum,x)
id = omp_get_thread_num()
numthreads = omp_get_num_threads()
partial_sum = 0.0

if (id == 0)  nthreads = numthreads

do i = id, num_steps-1, numthreads
   x = (i+0.5)*step
   partial_sum = partial_sum + 4.0/(1.0+x*x)
enddo
! ------------------------------------------------------------------------------
!$omp critical
full_sum = full_sum + partial_sum
!$omp end critical
! ------------------------------------------------------------------------------
!$omp end parallel
! ------------------------------------------------------------------------------
! this just scales the full sum
pi = step * full_sum

run_time = omp_get_wtime() - start_time

write(*,100) pi, run_time, nthreads
100     format('pi is ',f15.8,' in ',f8.3,'secs and ',i3,' threads')

end program main
