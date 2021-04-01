program main
use omp_lib

! ------------------------------------------------------------------------------
! parallelize the numerical integral of 4/(1+x^2) between 0 and 1.
! this integral is equal to pi.
! 
! this approach is the block parallelization of the integral for loop.
! ------------------------------------------------------------------------------



! parallel variables
integer, parameter :: max_threads = 4
integer ::  i, id, numthreads, nthreads
integer, parameter :: num_steps = 100000000

! integral variables
real*8 ::  pi, real_sum, step, x
real*8 :: start_time, run_time
real*8 :: sum(0:max_threads-1)

! ------------------------------------------------------------------------------
! this part is just to visualize where the parallel part is placing values
! for each thread.
integer, parameter :: id_ = 0
integer, parameter :: numthreads_= 4
integer, parameter :: num_steps_ = 10
do i = id_, num_steps_ - 1, numthreads_
   print '(i14)', i
enddo
! ------------------------------------------------------------------------------

step = 1.0 / num_steps

! this call can be avoided by declaring a bash variable like so:
! export OMP_NUM_THREADS=4
! this means the 'max_threads' variable can be avoided and
! re-compiling when changing the num of threads too.
call omp_set_num_threads(max_threads)
start_time = omp_get_wtime()

! ------------------------------------------------------------------------------
!$omp parallel private(id,x,numthreads)
id = omp_get_thread_num()
numthreads = omp_get_num_threads()
sum(id) = 0.0

! this is done to avoid "data racing".
! if not done, then 'nthreads' would update for every thread,
! and Fortran doesn't like that.
if (id == 0)  then
   nthreads = numthreads
endif

istart = id * num_steps / numthreads + 1
iend = (id+1) * num_steps / numthreads
if (id == (numthreads - 1)) then 
   iend = num_steps
endif

! this is the breakup of the integral into various threads.
! the idea is that each thread sums in blocks every '1' position.
do i = istart, iend
   x = (i - 0.5) * step
   sum(id) = sum(id) + 4.0 / (1.0 + x * x)
enddo
!$omp end parallel
! ------------------------------------------------------------------------------
! this collects all the threads-sums into one sum
pi = 0.0
do i = 0, nthreads-1
   pi = pi + sum(i)
enddo
pi = step * pi

run_time = omp_get_wtime() - start_time

write(*,100) pi, num_steps, run_time, nthreads
100 format('pi = ', f15.8, ',', i14, ' steps,',f8.3,' secs',i14,' nthreads')

end program main
