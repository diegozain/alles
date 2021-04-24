# openMP
diego domenzain

March 2021 @ Colorado School of Mines

## The OpenMP Common Core - Chapter II

Parallelize the integral:

π = ∫ 4 / (1+x^2) dx

Compile with,
```bash
gfortran -fopenmp file.f90
./a.out
```

### Simple openMP construct

```fortran90
use omp_lib

integer :: it, sum_
integer, parameter :: num_steps=100
real*8 :: start_time, run_time
integer, parameter :: nthreads=4

! starts chronometer 
start_time = omp_get_wtime()

! sets number of threads
call omp_set_num_threads(nthreads)

! tell compiler what vars are no-go
!$omp parallel private(it)

! start parallel block
!$omp do reduction (+:sum_) schedule(static)
do it = 1, num_steps
   sum_ = sum_ + 1
enddo
!$omp end do
!$omp end parallel

! stop chronometer
run_time = omp_get_wtime() - start_time

print *, ' run-time was ',run_time,'secs'
```

### Walking through the chapter

1. do it in **serial**: ```integral_serial.f90```
1. do it in **parallel cyclic**: ```integral_parallel_cyclic.f90```
1. do it in **parallel block**: ```integral_parallel_block.f90```
1. do it in **parallel cyclic** with **critical** synchronization: ```integral_parallel_critical.f90```
    * cyclic and block suffer from **false sharing**
    * this is fixed by placing a **critical construct**:
      * this forces one thread at a time to execute this block of code,
      * so if another thread wants to do this, it must wait for the previous thread to finish.
1. all these parallel constructs do similar stuff:
    * define ```omp parallel``` block,
    * get **ID** and **num. of threads**,
    * use **ID** and **num. of threads** to break up for loop in **cyclic** or **block**,
    * place a **critical construct** to avoid **false sharing**.
1. the **ID/num of threads/cyclic or block/critical** declaration can be handled by openMP on its own: ```integral_parallel_for.f90```
1. what about the **partial sums** that had to be defined to make the parallel scheme work? : ```integral_parallel_redu.f90```
    * well, you can **get rid** of them by using ```reduction(op:vars)```
1. can we schedule how many chunks go to each thread? : ```integral_parallel_redu.f90```
    * yes, with ```schedule(static or dynamic, ichunks)```
    * ```schedule(static,1)``` is **cyclic**
    * ```schedule(static)```   is **block**