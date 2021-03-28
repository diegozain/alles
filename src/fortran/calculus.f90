module calculus
! ------------------------------------------------------------------------------
! diego domenzain @ CSM 2021
!
! this module is about simple calculus stuff.
! It is based on my own Matlab code.
!
! * 'linspace' like Matlab.
! * 'differentiate' is a 1D discrete derivative.
! * 'integrate' is a 1D discrete antiderivative with constant of integration 
!               zero.
!
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
function linspace(tmin,tmax,nt) result(t)
   double precision :: tmin, tmax
   integer :: nt,it
   double precision :: t(nt), dt
   
   dt = abs(tmax-tmin)/dble(nt)
   
   do it=1,nt
      t(it) = tmin + dble(it-1)*dt
   enddo

end function linspace
! ------------------------------------------------------------------------------
subroutine differentiate(y,dt,nt)
   
   integer, intent(in) :: nt
   real*8, intent(in) :: dt
   real*8, intent(in out) :: y(nt)
   
   real*8 :: y_, y__, a,b,c,d
   
   a=-1.5/dt
   b=2/dt
   c=-0.5/dt
   d=2/dt
   
   y_ = y(1);
   y(1) = a*y(1) + b*y(2) + c*y(3);
   
   do it=2,(nt-2),2
    y__    = y(it)
    y(it)  = c*y_ - c*y(it+1)
    y_     = y(it+1)
    y(it+1)= c*y__- c*y(it+2)
   enddo
   
   ! NOTE: this part with ifs is ugly. 
   ! I should fix it one day.
   if (mod(nt,2)==1) then
    y__ = y(nt-1)
    y(nt-1) = c*y_ - c*y(nt)
    y(nt) = -c*y_ - d*y__ - a*y(nt)
   else
    y(nt) = -c*y__ - d*y_ - a*y(nt)
   endif
   
end subroutine differentiate
! ------------------------------------------------------------------------------
subroutine integrate(y,dt,nt)
   
   integer, intent(in) :: nt
   real*8, intent(in) :: dt
   real*8, intent(in out) :: y(nt)
   
   real*8 :: y_, y__, a,b
   
   a = 0.5*dt;
   b = 2*dt;
   
   y_ = y(3);
   y(3) = b * y(2);
   y(2) = a * (y(1) + y(2));
   do it=4,nt
    y__  = y(it);
    y(it)= y(it-2) + b*y_;
    y_   = y__;
   enddo
   y(1) = 0;

end subroutine integrate
! ------------------------------------------------------------------------------
end module calculus