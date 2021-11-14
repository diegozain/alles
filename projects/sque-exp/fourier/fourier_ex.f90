program fourier_ex
 ! -----------------------------------------------------------------------------
 ! diego domenzain @ CSM 2021
 !
 ! use the ** crazy ** fft function named dfork from Thomas Frobriger's lamb.f,
 ! and plot it in Python.
 !
 ! this is supposed to be a simple example on how to plot data produced by
 ! Fortran.
 !
 ! the idea is to use modules too:
 !
 !         $> gfortran -c ../../../src/fortran/fourier.f90 fourier_ex.f90
 !         $> gfortran fourier.o fourier_ex.o
 !         $> ./a.out
 !
 !         $> python visualize_f.py
 !
 !         $> rm *.o *.mod *.out *.dat fort.*
 ! -----------------------------------------------------------------------------
 use fourier
 implicit none
 ! -----------------------------------------------------------------------------
 double precision, parameter :: pi=3.1415926535898

 integer, parameter :: nt_pow2 = 4096
 integer, parameter :: nt = 4096
 complex*16 g(nt_pow2)
 complex*16 g_(nt_pow2)

 double precision HIN, RUECK
 parameter(HIN=-1.D0, RUECK=1.D0)

 double precision t, to_, source, dt, df
 source(t,to_)=dsin(pi*t/to_)**3

 integer i
 ! -----------------------------------------------------------------------------
 ! print some funny shit so user is amussed
 print *,' '
 print *,'   alright sir, Im gonna do the fft for you.'
 print *,'      If you so desire, run visualize.py to_ see the output of this.'
 print *,' '
 ! -----------------------------------------------------------------------------
 ! initialize g and g_ (of size a power of two!!)
 do i=1,nt_pow2
  g(i) = dcmplx(1.D0,0.D0)
  g_(i)= dcmplx(1.D0,0.D0)
 enddo

 to_= 0.5
 dt = 1e-3
 df = 1/(dt*dble(nt_pow2))

 ! initialize g and g_
 ! (on their own domain - not necessarily of size a power of 2)
 do i=1,nt
  t = dble(i-1)*dt
  g(i) = dcmplx(source(t,to_),0.D0)
  g_(i)= g(i)
 enddo
 ! -----------------------------------------------------------------------------
 ! run the fft,
 ! size of signal has to be a power of 2!!!!
 ! HIN  = -1 ---> frequency domain
 ! RUECK=  1 ---> time domain
 call dfork(nt,g_,HIN)
 ! call dfork(nt,g_,RUECK)
 ! -----------------------------------------------------------------------------
 ! write data
 ! -- data transform
 open(1, file = 'data_t.dat', status='unknown')
 do i = 1,nt_pow2
  write(1,"(E15.7,E15.7)") g(i)
 end do
 close(1)
 ! -- fourier transform
 open(2, file = 'data_f.dat', status='unknown')
 do i = 1,nt_pow2
  ! write(2,*) g_(i)
  write(2,"(E15.7,E15.7)") g_(i)
 end do
 close(2)
 ! -- time vector
 open(3, file = 'time.dat', status='unknown')
 do i = 1,nt_pow2
  write(3,*) dble(i-1)*dt
 end do
 close(3)
 ! -- frequency vector
 open(4, file = 'frequency.dat', status='unknown')
 do i = -nt_pow2/2,(nt_pow2/2 - 1)
  write(4,*) dble(i)*df
 end do
 close(4)
 ! -----------------------------------------------------------------------------
end program fourier_ex
