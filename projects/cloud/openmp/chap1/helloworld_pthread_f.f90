! This is a Fortran wrapper to call C Pthreads function
! To compile:
! % gfortran -c helloworld_pthread_f.f90 
! % gcc -c pthread_c.c 
! % gfortran helloworld_pthread_f.o pthread_c.o
! 
! this function is cool.
! it runs the C function "pthreads_c_" from the file pthread_c.c
! (more precisely from the compiled object pthread_c.o)
! 
! NOTE the underscore of "pthreads_c_" in pthread_c.c, but not in this file.

program main
   implicit none
   external :: pthread_c
   call pthreads_c()
end program main
