program lamb_solu
 ! -----------------------------------------------------------------------------
 ! diego domenzain @ CSM 2021
 ! 
 ! inspired by Thomas Frobriger's lamb.f, but easy to read and modify. 
 ! Honestly, his code is a mess.
 ! https://git.scc.kit.edu/Seitosh/Seitosh/-/blob/master/src/synt/misc/lamb.f
 !
 ! Uses the module 'fourier'. Results are plotted in Python.
 !
 !
 !         $> gfortran -c ../../../src/fortran/fourier.f90 lamb_solu.f90
 !         $> gfortran fourier.o lamb_solu.o
 !         $> ./a.out
 !
 !         $> python visualize_f.py
 !
 !         $> rm *.o *.mod *.out *.dat fort.*
 ! -----------------------------------------------------------------------------
 ! output is:
 ! vx_rx , vz_rx : shot-gathers on receivers at the free-surface on a 
 !                 horizontal line. Receivers are inline with the source.
 ! sz_time, time, frequency, rx : source signature on vz in time, ... ,
 !                                receiver coordinates on x.
 ! -----------------------------------------------------------------------------
 use fourier
 implicit none
 ! -----------------------------------------------------------------------------
 double precision, parameter :: pi=3.1415926535898
 complex*16, parameter :: IMAGINE=(0.D0,1.D0)
 
 integer, parameter :: nt = 512 ! 4096 ! a power of 2 so fft plays nice
 integer, parameter :: nrx= 10
 double precision, parameter :: drx = 10
 
 complex*16 rx(nrx)
 complex*16 sz(nt)
 complex*16 sz_time(nt)
 complex*16 vx_rx(nt , nrx)
 complex*16 vz_rx(nt , nrx)
 
 ! model parameters
 double precision, parameter :: vp= 500 ! 6000.D0 ! m/s
 double precision, parameter :: vs= 300 ! 4000.D0 ! m/s
 double precision, parameter :: rho= 1800 ! 2700.D0 ! kg/m^3
 double precision, parameter :: dt=0.001D0 ! s 0.5
 double precision mu, f(nt), time(nt), df, ftap
 complex*16 vp_c, vs_c, a, b, vx_num, vz_num, wvlt
 
 ! source parameters
 double precision, parameter :: to=0.032D0 ! s
 double precision, parameter :: Fo=1e4 ! kg*m/s^2
 
 ! slowness integration (s/m)
 double precision, parameter :: umin=0.D0    ! minimal slowness being calculated
 double precision, parameter :: uwil=0.1D-10 ! left cosine-taper slowness
 double precision, parameter :: uwir=1D-2  ! right cosine-taper slowness
 double precision, parameter :: umax=1.5D-2  ! maximal slowness being calculated
 integer, parameter :: nu=1000 ! 10000
 double precision u, du, ddu, utap, omegau, bess_arg
 complex*16 bess_o, bess_1
 complex*16 Ivx(nt,nrx), Ivz(nt,nrx)
 
 ! freq integration (Hz)
 double precision :: fmin=0.D0    ! minimal freq being calculated
 double precision, parameter :: fwil=0.0D0 ! left cosine-taper freq
 double precision :: fwir  ! right cosine-taper freq
 double precision :: fmax  ! maximal freq being calculated
 integer iminf, imaxf, nf
 
 ! rayleigh pole nonsense
 double precision, parameter :: raylim=1.D-100
 double precision reray, imray, absray
 complex*16 rayleigh
 
 ! quality nonsense
 double precision :: qvp=0.D0 ! 200.D0 ! qvp=0 -> elastic p-waves
 double precision :: qvs=0.D0 ! 200.D0 ! qvs=0 -> elastic s-waves
 
 ! fft stuff
 double precision HIN, RUECK
 parameter(HIN=-1.D0, RUECK=1.D0)
 
 ! source signature stuff
 double precision t, to_, source_impulse
 ! legit hammer
 source_impulse(t,to_)=dsin(pi*t/to_)**3
 ! ! derivative of hammer
 ! source_impulse(t,to_)=(3*pi/to_)*(dsin(pi*t/to_)**2)*dcos(pi*t/to_)
 ! ! antiderivative of hammer (by wolframalpha)
 ! source_impulse(t,to_)=(to_*dcos((pi*t)/to_)*(dcos((2*pi*t)/to_) - 5))/(6*pi)
 
 ! for loops
 integer it,irx,iu
 ! -----------------------------------------------------------------------------
 ! print funny shit so user is amused
 print *,' '
 print *,'   solution to Lamb problem: '
 print *,'            elastic wave propagation on a free-surface'
 print *,'            of a homogeneous half-space.'
 print *,' '
 print *,'   if you want to see model parameters, look inside me (lamb_solu.f90).'
 print *,' '
 ! -----------------------------------------------------------------------------
 ! init constant values
 mu = rho*vs*vs
 
 df = 1/(dt*dble(nt))
 nf = nt/2          ! this is only for one side of freq
 fmax = dble(nf)*df ! this is the actual max frequency given by dt and nt
 fwir = fmax*0.95   ! max taper value to be 95% of fmax
 
 ! recalculate frequency bandwidth from index-borders
 iminf=int(fmin/df)+1
 imaxf=int(fmax/df)+2
 imaxf=max(imaxf,nf)
 fmin=df*dble(iminf-1)
 fmax=df*dble(imaxf-1)
 
 du=(umax-umin)/dble(nu)
 
 if (qvp.eq.0.D0) then
   vp_c = dcmplx(vp)
 else
   vp_c = dcmplx(vp)*dcmplx(1.D0,0.5D0/qvp)
 endif
 if (qvs.eq.0.D0) then
   vs_c = dcmplx(vs)
 else
   vs_c = dcmplx(vs)*dcmplx(1.D0,0.5D0/qvs)
 endif
 
 ! receivers away from source
 do irx=1,nrx
   rx(irx) = dble(irx) * drx
 end do
 ! -----------------------------------------------------------------------------
 print *, '   largest frequency to sample =',fmax,'Hz'
 print *,' '
 ! -----------------------------------------------------------------------------
 ! initialize arrays
 do it=1,nt
  ! use impulse-response first
  sz(it)=dcmplx(1.D0,0.D0)
  ! set frequency
  f(it)=df*dble(it-1)
  ! set time
  time(it)=dt*dble(it-1)
  do irx=1,nrx
    ! arrays for sum-results
    vx_rx(it,irx)=dcmplx(0.D0,0.D0)
    vz_rx(it,irx)=dcmplx(0.D0,0.D0)
  end do
 end do
 
 ! true source init 
 if (to.gt.(3.D0*dt)) then
   do it=1,nt
     t=dble(it-1)*dt
     sz(it)=dcmplx(0.D0,0.D0)
     if (t.le.to) sz(it)=dcmplx(source_impulse(t,to),0.D0)
     sz_time(it) = dcmplx(Fo,0.D0)*sz(it)
   end do
   call dfork(nt,sz,HIN)
 endif
 ! -----------------------------------------------------------------------------
 ! slowness integral
 print *,' '
 print *,'   * slowness integral '
 ! -----------------------------------------------------------------------------
 do iu=0,nu
   u = dble(iu)*du+umin
   
   a = sqrt( (1/vp_c**2) - u**2 )
   b = sqrt( (1/vs_c**2) - u**2 )
   
   ddu=du
   if ((iu.eq.0).or.(iu.eq.nu)) ddu=du*0.5D0
   
   ! slowness taper
   utap = costap(umin,uwil,uwir,umax,u)
   
   ! rayleigh function
   rayleigh = (2*u**2 - 1/vs_c**2)**2 + (4*u**2)*a*b
   reray = dble(rayleigh)
   imray = dimag(rayleigh)
   absray= abs(rayleigh)
   ! rayleigh dipole nonsense
   if (absray.eq.0.D0) then
     rayleigh = raylim
   elseif (absray.lt.raylim) then
     rayleigh = dcmplx((reray/abs(reray)),(imray/abs(imray)))*raylim
   endif
   
   ! numerator for green fncs
   vx_num = dcmplx(u**2)*((1/vs_c**2)-dcmplx(2.D0*u**2)-dcmplx(2.D0,0.D0)*a*b)
   vz_num = dcmplx(u)*(1/vs_c**2)*a*IMAGINE
   
   ! now we loop on frequencies
   ! only positive and between fmin and fmax!!
   do it=iminf,imaxf
     omegau=2.D0*pi*f(it)*u
     
     ! now we loop over receivers
     do irx=1,nrx
       bess_arg=omegau*rx(irx)
       
       ! bessel functions of the first kind of order zero and one
       bess_o = dcmplx(bessel_jn(0,bess_arg))
       bess_1 = dcmplx(bessel_jn(1,bess_arg))
       
       Ivx(it,irx) = bess_1 * (vx_num/rayleigh*utap*ddu)
       Ivz(it,irx) = bess_o * (vz_num/rayleigh*utap*ddu)
       
       vx_rx = vx_rx + Ivx
       vz_rx = vz_rx + Ivz
     end do ! recs
   end do ! freqs
 end do ! slowness
 ! -----------------------------------------------------------------------------
 ! now comes the convolution with the source wavelet
 print *,'   * convolution with the source wavelet '
 ! -----------------------------------------------------------------------------
 ! loop on frequencies
 ! only positive and between fmin and fmax!!
 do it=iminf,imaxf
   ftap=costap(fmin,fwil,fwir,fmax,f(it))
   wvlt=dcmplx(ftap*(-Fo)*f(it)/mu)*sz(it)
   do irx=1,nrx
     vx_rx(it,irx)=vx_rx(it,irx)*wvlt
     vz_rx(it,irx)=vz_rx(it,irx)*wvlt
   end do
 end do
 ! -----------------------------------------------------------------------------
 ! now back to the time domain
 print *,'   * back to the time domain '
 ! -----------------------------------------------------------------------------
 do irx=1,nrx
   
   ! conjugate so time series is real
   do it=2,nf-1
     vx_rx(nt-it+2,irx)=dconjg(vx_rx(it,irx))
     vz_rx(nt-it+2,irx)=dconjg(vz_rx(it,irx))
   end do
   
   ! vx wavefield
   do it=1,nt
     sz(it)=vx_rx(it,irx)
   end do
   call dfork(nt, sz, RUECK)
   do it=1,nt
     vx_rx(it,irx)=sz(it)
   end do
   
   ! vz wavefield
   do it=1,nt
     sz(it)=vz_rx(it,irx)
   end do
   call dfork(nt, sz, RUECK)
   do it=1,nt
     vz_rx(it,irx)=sz(it)
   end do
   
 end do ! receivers
 ! -----------------------------------------------------------------------------
 ! write data
 print *,'   * write data '
 print *,' '
 ! -----------------------------------------------------------------------------
 ! -- time, frequency, and sz on time
 open(1, file = 'time.dat', status='unknown')
 open(2, file = 'frequency.dat', status='unknown')
 open(3, file = 'sz_time.dat', status='unknown')
 do it = 1,nt
  write(1,"(E15.7,E15.7)") time(it)
  write(2,"(E15.7,E15.7)") f(it)
  write(3,"(E15.7,E15.7)") sz_time(it)
 end do  
 close(1)
 close(2)
 close(3)
 ! -- receivers
 open(1, file = 'rx.dat', status='unknown')
 do irx=1,nrx
   write(1,"(E15.7,E15.7)") rx(irx)
 end do 
 close(1)
 ! -- vx and vz on receivers
 open(1, file = 'vx_rx.dat', status='unknown')
 open(2, file = 'vz_rx.dat', status='unknown')
 do irx=1,nrx ! irx=1,nrx  it=1,nt
   do it=1,nt ! it=1,nt    irx=1,nrx
     write(1,"(E15.7,E15.7)") vx_rx(it,irx)
     write(2,"(E15.7,E15.7)") vz_rx(it,irx)
   end do
 end do  
 close(1)
 close(2)
 ! -----------------------------------------------------------------------------
end program lamb_solu