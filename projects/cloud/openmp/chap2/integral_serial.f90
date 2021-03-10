program main
use omp_lib
implicit none
! ------------------------------------------------------------------------------
! compute the numerical integral of 4/(1+x^2) between 0 and 1.
! this integral is equal to pi.
! 
! this approach is the serial integral for loop.
! ------------------------------------------------------------------------------

integer :: i
integer, parameter :: num_steps = 100000000
real*8 :: x, pi, sum, step
real*8 :: start_time, run_time

sum = 0.0

step = 1.0 / num_steps
start_time = omp_get_wtime()

do i = 1, num_steps
   x = (i - 0.5) * step
   sum = sum + 4.0 / (1.0 + x * x)
enddo
pi = step * sum

run_time = omp_get_wtime() - start_time

write(*,100) pi, num_steps,  run_time
100 format('pi = ', f15.8, ',', i14, ' steps,',f8.3,' secs')

end program main 
