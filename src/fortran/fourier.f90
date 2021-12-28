module fourier
! ------------------------------------------------------------------------------
! diego domenzain @ CSM 2021
!
! this module is about fourier.
!
! * 'dfork' is a discrete fast fourier transform by Thomas Forbriger.
! * 'costap' is a cosine taper by Thomas Forbriger.
!
! Thomas Forbriger's code: https://git.scc.kit.edu/Seitosh/Seitosh
! ------------------------------------------------------------------------------
use mesher3d
use calculus
contains
! ------------------------------------------------------------------------------
! the following fast fourier transform is just copied from
! the program refseis from j. ungerer
!
! this code was originally published by gerhard müller
! in his lecture notes on digital signal processing:
! gerhard mueller, 1999. digitale signalverarbeitung. skriptum zur
! gleichnamigen vorlesung. institut für meteorologie und geophysik,
! universität frankfurt.
!
! the original algorithm appears to be due to claerbout, j.f.,
! "fundamentals of geophysical data processing with applications
! to petroleum prospecting", mcgraw-hill, 1976.
!
!
! 23456789012345678901234567890123456789012345678901234567890123456789012
! *** subroutine fuer fouriertransformation ************************** !
! die von gerherd mueller verwendetet schnelle fouriertransformation   !
! fork wurde umgeschrieben fuer double complex                         !
! es muessen implementiert sein: double complex,dcmplx,cdexp           !
!                                                                      !
! zum verfahren der schnellen fouriertransformation(fft) und zur ar-   !
! beitsweise von fork siehe g.mueller: digitale signalverarbeitung i,  !
! vorlesungsmanuskript.                                                !
!                                                                      !
! variablen:                                                           !
!    lx       seismogrammlaenge bzw. anzahl der stuetzstellen,abtast-  !
!             werte des seismogramms/spektrums.muss eine zeier-potenz  !
!             sein.                                                    !
!    cx(lx)   feld auf dessen realteil die funktionswerte der zeit-    !
!             funktion stehen und nach transformation ihre fourierko-  !
!             effizienten.                                             !
!    signi    signi=-1.d0 bedeutet berechnung der fourierkoeffizienten !
!             signi=+1.d0 bedeutet ruecktransformation                 !
! ******************************************************************** !
      subroutine dfork(lx,cx,signi)
      integer     i,istep,j,l,lx,m
      real*8      sc,pi,signi
      complex*16  cx(lx),carg,cw,ctemp

      pi=3.14159265358979d0
      j=1
      sc=1.d0/dble(lx)
      sc=dsqrt(sc)
      do 5  i=1,lx
      if(i.gt.j) goto 2
      ctemp=cx(j)*sc
      cx(j)=cx(i)*sc
      cx(i)=ctemp
2     m=lx/2
3     if(j.le.m) goto 5
      j=j-m
      m=m/2
      if(m.ge.1) goto 3
5     j=j+m
      l=1
6     istep=2*l
      do 8  m=1,l
      carg=dcmplx(0.,1.)*(pi*signi*dble(m-1))/dble(l)
      cw=cdexp(carg)
      do 8 i=m,lx,istep
      ctemp=cw*cx(i+l)
      cx(i+l)=cx(i)-ctemp
8     cx(i)=cx(i)+ctemp
      l=istep
      if(l.lt.lx) goto 6
      return
      end subroutine dfork
! ------------------------------------------------------------------------------
double precision function costap(min,wil,wir,max,val)
  double precision min,wil,wir,max,val
  double precision pi
  parameter(pi=3.1415926535898D0)
  if (val.lt.min) then
   costap=0.D0
  elseif(val.le.wil) then
   if (wil.eq.min) then
      costap=0.D0
   else
      costap=0.5D0-0.5D0*dcos((val-min)*pi/(wil-min))
   endif
  elseif(val.gt.max) then
   costap=0.D0
  elseif(val.ge.wir) then
   if (max.eq.wir) then
      costap=0.D0
   else
      costap=0.5D0+0.5D0+dcos((val-wir)*pi/(max-wir))
   endif
  else
   costap=1.D0
  endif
  return
end function costap
! ------------------------------------------------------------------------------
subroutine fft_3d(sigm,nx,ny,nz,signi)
  ! 𝓕(σ)       : signi=-1.D0
  ! 𝓕^-1(𝓕(σ)) : signi= 1.D0
  !
  !⚠️ although input σ can be thought as being in ℝ (for 𝓕),
  !⚠️ it must be inputed as type ℂ.
  real*8, intent(in)  :: signi
  integer, intent(in) :: nx, ny, nz
  complex*16, intent(in out) :: sigm(ny,nx,nz)

  do ix=1,nx
    do iy=1,ny
      call dfork(nz,sigm(iy,ix,:),signi)
    enddo
  enddo

  do iz=1,nz
    do iy=1,ny
      call dfork(nx,sigm(iy,:,iz),signi)
    enddo
  enddo

  do iz=1,nz
    do ix=1,nx
      call dfork(ny,sigm(:,ix,iz),signi)
    enddo
  enddo
end subroutine fft_3d
! ------------------------------------------------------------------------------
subroutine fftshift_real(x,nx)
  ! ⚠️ nx HAS to be even.
  integer, intent(in) :: nx
  double precision, intent(in out) :: x(nx)
  double precision :: x_
  integer :: ix, ix_, ixshift

  ixshift= floor(nx*0.5);

  do ix=1,ixshift
    ix_   = ix+ixshift;
    x_    = x(ix_);
    x(ix_)= x(ix);
    x(ix) = x_;
  enddo
end subroutine fftshift_real
! ------------------------------------------------------------------------------
subroutine fftshift_cmplx(x,nx)
  ! ⚠️ nx HAS to be even.
  integer, intent(in) :: nx
  complex*16, intent(in out) :: x(nx)
  complex*16 :: x_
  integer :: ix, ix_, ixshift

  ixshift= floor(nx*0.5);

  do ix=1,ixshift
    ix_   = ix+ixshift;
    x_    = x(ix_);
    x(ix_)= x(ix);
    x(ix) = x_;
  enddo
end subroutine fftshift_cmplx
! ------------------------------------------------------------------------------
subroutine smooth_cube(sigm,nx,ny,nz,ax,ay,az)
  !  f      f.shifted     𝓕(σ)    𝓕^-1    Re
  ! 🎲   →     🎲    ⊙   🎲  →   🎲   →  🎲
  ! ℝ          ℝ          ℂ       ℂ        ℝ
  !
  !⚠️ although input σ is thought as being in ℝ,
  !⚠️ it must be inputed as type ℂ.
  integer, intent(in)          :: nx, ny, nz
  double precision, intent(in) :: ax, ay, az
  complex*16, intent(in out)   :: sigm(ny,nx,nz)
  double precision             :: filter(ny,nx,nz), x(nx), y(ny), z(nz)
  double precision             :: ax_, ay_, az_
  integer                      :: ix, iy, iz, iyxz
  ! ----------------------------------------------------------------------------
  ! get 𝓕(σ)
  call fft_3d(sigm,nx,ny,nz,-1.D0)
  ! ----------------------------------------------------------------------------
  ! build stuff
  x = linspace(0.D0,dble(nx-1),nx)
  y = linspace(0.D0,dble(ny-1),ny)
  z = linspace(0.D0,dble(nz-1),nz)

  ax_ = dble(nx)*(ax*0.5)*dsqrt(2.D0)
  ay_ = dble(ny)*(ay*0.5)*dsqrt(2.D0)
  az_ = dble(nz)*(az*0.5)*dsqrt(2.D0)
  ! ----------------------------------------------------------------------------
  ! build filter
  do iyxz=1,nx*ny*nz
    call get_ixyz(ix,iy,iz,iyxz,nx,ny,nz)
    ! 🔅 low-pass gaussian filter
    filter(iy,ix,iz) = exp(-( ((x(ix)-x(nx/2 +1))/ax_)**2 + ((y(iy)-y(ny/2 +1))/ay_)**2 + ((z(iz)-z(nz/2 +1))/az_)**2 ))
    ! ! 🔆 high-pass gaussian filter
    ! filter(iy,ix,iz) = 1 - exp(-( ((x(ix)-x(nx/2))/ax_)**2 + ((y(iy)-y(ny/2))/ay_)**2 + ((z(iz)-z(nz/2))/az_)**2 ))
  enddo
  ! ----------------------------------------------------------------------------
  ! fftshift the filter
  do ix=1,nx
    do iy=1,ny
      call fftshift_real(filter(iy,ix,:),nz)
    enddo
  enddo
  do iz=1,nz
    do iy=1,ny
      call fftshift_real(filter(iy,:,iz),nx)
    enddo
  enddo
  do iz=1,nz
    do ix=1,nx
      call fftshift_real(filter(:,ix,iz),ny)
    enddo
  enddo
  ! ----------------------------------------------------------------------------
  ! filter ⊙ σ
  do iyxz=1,nx*ny*nz
    call get_ixyz(ix,iy,iz,iyxz,nx,ny,nz)
    sigm(iy,ix,iz) = dcmplx(filter(iy,ix,iz),0.D0) * sigm(iy,ix,iz)
  enddo
  ! ----------------------------------------------------------------------------
  ! 𝓕^-1
  call fft_3d(sigm,nx,ny,nz,1.D0)
end subroutine smooth_cube
! ------------------------------------------------------------------------------
subroutine nextpow2(nt_,nt)
  integer, intent(in) :: nt
  integer, intent(in out) :: nt_
  double precision :: ntdbl_

  ntdbl_ = log(dble(2*nt)) / log(2.0)
  nt_ = int(floor(ntdbl_))
end subroutine nextpow2
! ------------------------------------------------------------------------------
subroutine vecpad2(vec_,vec,nt_,nt)
  integer, intent(in) :: nt, nt_
  complex*16, intent(in) :: vec(nt)
  complex*16, intent(in out) :: vec_(nt_)

  do it=1,nt
    vec_(it) = vec(it)
  enddo
  do it=1,nt_-nt
    vec_(nt+it) = 0 ! vec(nt) its actually better with =0 😮
  enddo
end subroutine vecpad2
! ------------------------------------------------------------------------------
subroutine vectrim2(vec,vec_,nt,nt_)
  integer, intent(in) :: nt, nt_
  complex*16, intent(in) :: vec_(nt_)
  complex*16, intent(in out) :: vec(nt)

  do it=1,nt
    vec(it) = vec_(it)
  enddo
end subroutine vectrim2
! ------------------------------------------------------------------------------
subroutine matpad2(mat_,mat,nrows_,ncols_,nrows,ncols)
  integer, intent(in) :: nrows_,ncols_,nrows,ncols
  complex*16, intent(in) :: mat(nrows,ncols)
  complex*16, intent(in out) :: mat_(nrows_,ncols_)

  ! fill main matrix
  do icol=1,ncols
    do irow=1,nrows
      mat_(irow,icol) = mat(irow,icol)
    enddo
  enddo
  ! fill right side
  do icol=1,ncols_-ncols
    do irow=1,nrows
      mat_(irow,ncols+icol) = 0 ! mat(irow,ncols) its actually better with =0 😮
    enddo
  enddo
  ! fill bottom side
  do icol=1,ncols_
    do irow=1,nrows_-nrows
      mat_(nrows+irow,icol) = 0 ! mat(nrows,icol) its actually better with =0 😮
    enddo
  enddo
end subroutine matpad2
! ------------------------------------------------------------------------------
subroutine matrim2(mat,mat_,nrows,ncols,nrows_,ncols_)
  integer, intent(in) :: nrows_,ncols_,nrows,ncols
  complex*16, intent(in) :: mat_(nrows_,ncols_)
  complex*16, intent(in out) :: mat(nrows,ncols)

  ! fill main matrix
  do icol=1,ncols
    do irow=1,nrows
      mat(irow,icol) = mat_(irow,icol)
    enddo
  enddo
end subroutine matrim2
! ------------------------------------------------------------------------------
subroutine cubpad2(cub_,cub,nx_,ny_,nz_,nx,ny,nz)
  integer, intent(in) :: nx_,ny_,nz_,nx,ny,nz
  complex*16, intent(in) :: cub(ny,nx,nz)
  complex*16, intent(in out) :: cub_(ny_,nx_,nz_)

  ! fill main matrix
  do ix=1,nx
    do iy=1,ny
      do iz=1,nz
        cub_(iy,ix,iz) = cub(iy,ix,iz)
      enddo
    enddo
  enddo
  ! fill x-side
  do ix=1,nx_-nx
    do iy=1,ny
      do iz=1,nz
        cub_(iy,nx+ix,iz) = 0 ! cub(iy,nx,iz) its actually better with =0 😮
      enddo
    enddo
  enddo
  ! fill z-side
  do ix=1,nx_
    do iy=1,ny
      do iz=1,nz_-nz
        cub_(iy,ix,nz+iz) = 0 ! cub(iy,ix,nz) its actually better with =0 😮
      enddo
    enddo
  enddo
  ! fill y-side
  do ix=1,nx_
    do iy=1,ny_-ny
      do iz=1,nz_
        cub_(ny+iy,ix,iz) = 0 ! cub(ny,ix,iz) its actually better with =0 😮
      enddo
    enddo
  enddo
end subroutine cubpad2
! ------------------------------------------------------------------------------
subroutine cubtrim2(cub,cub_,nx,ny,nz,nx_,ny_,nz_)
  integer, intent(in) :: nx_,ny_,nz_,nx,ny,nz
  complex*16, intent(in) :: cub_(ny_,nx_,nz_)
  complex*16, intent(in out) :: cub(ny,nx,nz)

  ! fill main matrix
  do ix=1,nx
    do iy=1,ny
      do iz=1,nz
        cub(iy,ix,iz) = cub_(iy,ix,iz)
      enddo
    enddo
  enddo
end subroutine cubtrim2
! ------------------------------------------------------------------------------
end module fourier
