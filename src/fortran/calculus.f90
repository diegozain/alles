module calculus
! ------------------------------------------------------------------------------
! diego domenzain @ CSM 2021
!
! this module is about simple calculus stuff.
! It is based on my own Matlab code.
!
! * 'linspace' like Matlab.
! * 'differentiate' is a 1D discrete derivative (second order accurate).
! * 'differentiate_o6' is a 1D discrete derivative (sixth order accurate).
! * 'integrate' is a 1D discrete antiderivative with constant of integration 
!               zero (second order accurate).
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
   
   real*8 :: y_, y__, a,b,c
   
   a=-1.5/dt
   b=2/dt
   c=-0.5/dt
   
   y_ = y(1)
   y(1) = a*y(1) + b*y(2) + c*y(3)
   
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
    y(nt) = -c*y_ - b*y__ - a*y(nt)
   else
    y(nt) = -c*y__ - b*y_ - a*y(nt)
   endif
   
end subroutine differentiate
! ------------------------------------------------------------------------------
subroutine differentiate_o6(y,dt,nt)
   
   integer, intent(in) :: nt
   real*8, intent(in) :: dt
   real*8, intent(in out) :: y(nt)
   
   real*8 :: y_(9), y__(9), fwd(7), ctd(7), bwd(7)
   integer :: it
   
   fwd(1)= -49./20. * (1/dt)
   fwd(2)= 6. * (1/dt)
   fwd(3)= -15./2. * (1/dt)
   fwd(4)= 20./3. * (1/dt)
   fwd(5)= -15./4. * (1/dt)
   fwd(6)= 6./5. * (1/dt)
   fwd(7)= -1./6. * (1/dt)

   ctd(1)= -1./60. * (1/dt)
   ctd(2)= 3./20. * (1/dt)
   ctd(3)= -3./4. * (1/dt)
   ctd(4)= 0.
   ctd(5)= 3./4. * (1/dt)
   ctd(6)= -3./20. * (1/dt)
   ctd(7)= 1./60. * (1/dt)

   bwd(1)= 1./6. * (1/dt)
   bwd(2)= -6./5. * (1/dt)
   bwd(3)= 15./4. * (1/dt)
   bwd(4)= -20./3. * (1/dt)
   bwd(5)= 15./2. * (1/dt)
   bwd(6)= -6. * (1/dt)
   bwd(7)= 49./20. * (1/dt)
   
   y_(1:9) = y(1:9)
   y(1) = dot_product(fwd,y_(1:7))
   y(2) = dot_product(fwd,y_(2:8))
   y(3) = dot_product(fwd,y_(3:9))
   
   do it=4,(nt-8),3
      y__(1:9)= y(it:(it+8))

      y(it)  = dot_product(ctd,y_(1:7))
      y(it+1)= dot_product(ctd,y_(2:8))
      y(it+2)= dot_product(ctd,y_(3:9))

      y_=y__
   enddo
   
   ! NOTE: this part with ifs is ugly. 
   ! I should fix it one day.
   if (mod(nt,3)==2) then
      it=nt-10

      y__(1:9)= y((it+2):(it+10))

      y(it+3)= dot_product(ctd,y_(1:7))
      y(it+4)= dot_product(ctd,y_(2:8))
      y(it+5)= dot_product(ctd,y_(3:9))

      y(it+6)= dot_product(bwd,y_(1:7))
      y(it+7)= dot_product(bwd,y_(2:8))
      y(it+8)= dot_product(bwd,y_(3:9))

      y_=y__

      y(nt-1)= dot_product(bwd,y_(2:8))
      y(nt)  = dot_product(bwd,y_(3:9))
   elseif (mod(nt,3)==1) then
      it=nt-9

      y__(1:9)= y((it+1):(it+9))

      y(it+3)= dot_product(ctd,y_(1:7))
      y(it+4)= dot_product(ctd,y_(2:8))
      y(it+5)= dot_product(ctd,y_(3:9))

      y(it+6)= dot_product(bwd,y_(1:7))
      y(it+7)= dot_product(bwd,y_(2:8))
      y(it+8)= dot_product(bwd,y_(3:9))

      y_=y__

      y(nt)  = dot_product(bwd,y_(3:9))
   elseif (mod(nt,3)==0) then
      it=nt-8

      y(it+3)= dot_product(ctd,y_(1:7))
      y(it+4)= dot_product(ctd,y_(2:8))
      y(it+5)= dot_product(ctd,y_(3:9))

      y(it+6)= dot_product(bwd,y_(1:7))
      y(it+7)= dot_product(bwd,y_(2:8))
      y(it+8)= dot_product(bwd,y_(3:9))
   endif
   
end subroutine differentiate_o6
! ------------------------------------------------------------------------------
subroutine integrate(y,dt,nt)
   
   integer, intent(in) :: nt
   real*8, intent(in) :: dt
   real*8, intent(in out) :: y(nt)
   
   real*8 :: y_, y__, a,b
   
   a = 0.5*dt
   b = 2*dt
   
   y_ = y(3)
   y(3) = b * y(2)
   y(2) = a * (y(1) + y(2))
   do it=4,nt
    y__  = y(it)
    y(it)= y(it-2) + b*y_
    y_   = y__
   enddo
   y(1) = 0

end subroutine integrate
! ------------------------------------------------------------------------------
end module calculus