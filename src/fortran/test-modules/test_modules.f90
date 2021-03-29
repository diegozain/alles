program test_modules
 ! -----------------------------------------------------------------------------
 ! diego domenzain @ CSM 2021
 ! 
 ! test modules.
 !
 ! compile and execute:
 !
 !         $> gfortran -c ../calculus.f90 test_modules.f90
 !         $> gfortran calculus.o test_modules.o
 !         $> ./a.out
 !
 !         $> python vis_test.py
 !
 !         $> rm *.o *.mod *.out *.dat fort.*
 ! -----------------------------------------------------------------------------
 use calculus
 implicit none
 ! -----------------------------------------------------------------------------
 double precision, parameter :: pi=3.1415926535898
 
 integer, parameter :: nt=10000
 integer :: it
 double precision, parameter :: tmin=0, tmax=2*pi
 double precision :: t(nt),y(nt),y_diff(nt),y_inte(nt),dt
 ! -----------------------------------------------------------------------------
 ! print some funny shit so user is amussed 
 print *,' '
 print *,'   alright sir, Im gonna test if your modules work.'
 print *,' '
 print *,'      If you so desire, run:'
 print *,'         python vis_test.py'
 print *,'      to see the output of this.'
 print *,' '
 ! -----------------------------------------------------------------------------
 t = linspace(tmin,tmax,nt)
 dt = t(2)-t(1)
 ! -----------------------------------------------------------------------------
 do it=1,nt
  y(it) = dsin(t(it))
  y_diff(it) = y(it)
  y_inte(it) = y(it)
 enddo
 
 call differentiate_o6(y_diff,dt,nt) ! sixth order derivative
 ! call differentiate(y_diff,dt,nt) ! second order derivative
 call integrate(y_inte,dt,nt) ! second order antiderivative
 ! -----------------------------------------------------------------------------
 ! write data
 ! -- data transform
 open(1, file = 't.dat', status='unknown')
 open(2, file = 'y.dat', status='unknown')
 open(3, file = 'y_diff.dat', status='unknown')
 open(4, file = 'y_inte.dat', status='unknown')
 do it = 1,nt
  write(1,"(E15.7)") t(it)
  write(2,"(E15.7)") y(it)
  write(3,"(E15.7)") y_diff(it)
  write(4,"(E15.7)") y_inte(it)
 end do  
 close(1)
 close(2)
 close(3)
 close(4)
 ! -----------------------------------------------------------------------------
end program test_modules