program umfpack_cool
 ! -----------------------------------------------------------------------------
 ! diego domenzain @ CSM 2021
 ! 
 ! solve a sparse system of equations using umfpack: Ax=b
 ! 
 ! -----------------------------------------------------------------------------
 ! for compiling:
 ! 
 ! gfortran -c umfpack_cool.f90
 ! gcc -c ../../../../SuiteSparse/UMFPACK/Demo/umf4_f77wrapper.c
 ! gfortran umfpack_cool.o umf4_f77wrapper.o -L../../../../SuiteSparse/lib -lumfpack 
 ! -----------------------------------------------------------------------------
 implicit none
 ! -----------------------------------------------------------------------------
 ! sparse matrix A
 integer, parameter :: nz=12, n=5
 integer Ai(nz), Ap(n+1)
 double precision Ax(nz)
 
 ! vector b
 double precision b(n)
 ! vector x
 double precision x(n)
 
 ! umfpack stuff (from umfpack_userguide.pdf chap 7)
 integer symbolic, numeric, filenum, status
 double precision control(20), info(90)
 
 ! printing
 integer i
 ! -----------------------------------------------------------------------------
 ! fill A & b
 Ap = [ 0,2,5,9,10,12 ]
 Ai = [ 0,1,0,2,4,1,2,3,4,2,1,4 ]
 Ax = [ 2,3,3,-1,4,4,-3,1,2,2,6,1 ]
 
 b = [ 8,45,-3,3,19 ]
 ! -----------------------------------------------------------------------------
 ! print funny stuff so user is amused
 write(*,*) ' Ap = ' , (Ap(i),i=1,n+1)
 write(*,*) ' Ai = ' , (Ai(i),i=1,nz)
 write(*,*) ' Ax = ' , (Ax(i),i=1,nz)
 
 write(*,*) ' b = ' , (b(i),i=1,n)
 ! -----------------------------------------------------------------------------
 ! solve using Umfpack.f
 call umf4sym(n, n, Ap, Ai, Ax, symbolic, control, info)
 call umf4num(Ap, Ai, Ax, symbolic, numeric, control, info)
 call umf4sol(0, x, b, numeric, control, info)
 ! -----------------------------------------------------------------------------
 ! print solution
 write(*,*) ' x = ' , (x(i),i=1,n)
 ! -----------------------------------------------------------------------------
end program umfpack_cool
