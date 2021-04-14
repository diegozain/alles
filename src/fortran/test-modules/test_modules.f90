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
 open(1,file = 't.dat', status='unknown',form='unformatted',access='stream')
 open(2,file = 'y.dat', status='unknown',form='unformatted',access='stream')
 open(3,file = 'y_diff.dat',status='unknown',form='unformatted',access='stream')
 open(4,file = 'y_inte.dat',status='unknown',form='unformatted',access='stream')

 write(1) (t(it),it=1,nt)
 write(2) (y(it),it=1,nt)
 write(3) (y_diff(it),it=1,nt)
 write(4) (y_inte(it),it=1,nt)
 
 close(1)
 close(2)
 close(3)
 close(4)
 ! -----------------------------------------------------------------------------
 ! in Matlab
 ! fid=fopen('t.dat','r');
 ! t=fread(fid,'double');
 ! fclose(fid);
 ! -----------------------------------------------------------------------------
end program test_modules