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
integer :: i
integer, parameter :: num_steps=100000000
integer, parameter :: nthreads=4

! integral variables
real*8 :: pi, step, sum_, x
real*8 :: start_time, run_time

sum_ = 0.0
step = 1.0 / num_steps
! ------------------------------------------------------------------------------
! starts chronometer
start_time = omp_get_wtime()

! sets number of threads
call omp_set_num_threads(nthreads)

!$omp parallel private(i,x)
! ------------------------------------------------------------------------------
! reduction(op:vars)     -- this will enable 'sum_' to be spread privately to 
!                            each tread, and then get them all.
!                            it does the same as breaking up the for loop into 
!                            partial sums, but is much cleaner code.
! schedule(static,ichunk) -- this will break up the for loop into chunks of 
!                            size ichunks to each thread.
!                            it is OPTIONAL, but allows for more control.
!
!                            schedule(static,1) is cyclic
!                            schedule(static)   is block

!$omp do reduction(+:sum_) schedule(static)
do i = 1, num_steps
   x = (i+0.5)*step
   sum_ = sum_ + 4.0/(1.0+x*x)
enddo
!$omp end do
! ------------------------------------------------------------------------------
!$omp end parallel

! stops chronometer
run_time = omp_get_wtime() - start_time
! ------------------------------------------------------------------------------
! this just scales the full sum_
pi = step * sum_

write(*,100) pi, run_time
100     format('pi is ',f15.8,' in ',f8.3,'secs')

end program main
