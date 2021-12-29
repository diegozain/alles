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
function logspace(tmin,tmax,nt) result(t)
   double precision :: tmin, tmax
   integer:: nt,it
   double precision :: t_lin(nt), t(nt)

   t_lin = linspace(tmin,tmax,nt)

   do it=1,nt
      t(it) = 10**t_lin(it)
   enddo

end function logspace
! ------------------------------------------------------------------------------
subroutine flip(x,nx)
  integer, intent(in) :: nx
  double precision, intent(in out) :: x(nx)

  integer :: ix, nx_
  double precision :: x_

  if (mod(nx,2)==0) then
    nx_ = nx/2
  elseif (mod(nx,2)==1) then
    nx_ = (nx-1) / 2
  endif

  do ix=1,nx_
    x_ = x(ix)
    x(ix) = x(nx-ix+1)
    x(nx-ix+1) = x_
  enddo
end subroutine flip
! ------------------------------------------------------------------------------
subroutine normali(array,narray)
  integer, intent(in) :: narray
  double precision, intent(in out) :: array(narray)

  integer :: iarray(narray)
  double precision :: maxarray, array_(narray)

  do iarr=1,narray
    array_(iarr) = abs(array(iarr))
  enddo

  call quicksort(array_,iarray,narray)
  maxarray = array_(narray)

  do iarr=1,narray
    array(iarr) = array(iarr) / maxarray
  enddo
end subroutine normali
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
subroutine quicksort(array,iarray,narray)
  ! 'array'  is the array to sort,
  ! 'iarray' are the indexes in the order they were sorted,
  ! 'narray' is the length of the array.
  !
  ! This version maintains its own stack, to avoid needing to call
  ! itself recursively. By always pushing the larger "half" to the
  ! stack, and moving directly to calculate the smaller "half",
  ! it can guarantee that the stack needs no more than log_2(N)
  ! entries.
  integer, intent(in) :: narray
  integer, intent(in out) :: iarray(narray)
  double precision, intent(in out) :: array(narray)
  double precision :: temp, pivot
  integer :: i, j, left, right, low, high, itemp
  ! If your compiler lacks storage_size(), replace
  ! storage_size(i) by 64
  integer :: stack(2,storage_size(i)), stack_ptr

  low=1
  stack_ptr=1
  high=narray

  do i=1,narray
    iarray(i)=i
  enddo

  do
     if (high-low.lt.50) then ! use insertion sort on small arrays
        do i=low+1,high
           temp=array(i)
           itemp=iarray(i)
           do j=i-1,low,-1
              if (array(j).le.temp) exit
              array(j+1)=array(j)
              iarray(j+1)=iarray(j)
           enddo
           array(j+1)=temp
           iarray(j+1)=itemp
        enddo
        ! now pop from stack
        if (stack_ptr.eq.1) return
        stack_ptr=stack_ptr-1
        low=stack(1,stack_ptr)
        high=stack(2,stack_ptr)
        cycle
     endif

     ! find median of three pivot
     ! and place sentinels at first and last elements
     temp=array((low+high)/2)
     array((low+high)/2)=array(low+1)

     iarray((low+high)/2)=iarray(low+1)
     itemp=iarray((low+high)/2)
     if (temp.gt.array(high)) then
        array(low+1)=array(high)
        array(high)=temp

        iarray(low+1)=iarray(high)
        iarray(high) =itemp
     else
       array(low+1)=temp
       iarray(low+1)=itemp
     endif
     if (array(low).gt.array(high)) then
        temp=array(low)
        array(low)=array(high)
        array(high)=temp

        itemp=iarray(low)
        iarray(low)=iarray(high)
        iarray(high)=itemp
     endif
     if (array(low).gt.array(low+1)) then
        temp=array(low)
        array(low)=array(low+1)
        array(low+1)=temp

        itemp=iarray(low)
        iarray(low)=iarray(low+1)
        iarray(low+1)=itemp
     endif
     pivot=array(low+1)

     left=low+2
     right=high-1
     do
        do while(array(left).lt.pivot)
           left=left+1
        enddo
        do while(array(right).gt.pivot)
           right=right-1
        enddo
        if (left.ge.right) exit
        temp=array(left)
        array(left)=array(right)
        array(right)=temp

        itemp=iarray(left)
        iarray(left)=iarray(right)
        iarray(right)=itemp

        left=left+1
        right=right-1
     enddo
     if (left.lt.(low+high)/2) then
        stack(1,stack_ptr)=left
        stack(2,stack_ptr)=high
        stack_ptr=stack_ptr+1
        high=left-1
     else
        stack(1,stack_ptr)=low
        stack(2,stack_ptr)=left-1
        stack_ptr=stack_ptr+1
        low=left
     endif
  enddo
end subroutine quicksort
! ------------------------------------------------------------------------------
subroutine mean(x_,x,nx)
  integer, intent(in) :: nx
  double precision, intent(in) :: x(nx)
  double precision, intent(in out) :: x_

  integer :: ix
  ! ----------------------------------------------------------------------------
  x_=0
  do ix=1,nx
    x_ = x_ + x(ix)
  enddo
  x_ = x_ / nx
end subroutine mean
! ------------------------------------------------------------------------------
subroutine std(x_,x,nx)
  integer, intent(in) :: nx
  double precision, intent(in) :: x(nx)
  double precision, intent(in out) :: x_

  double precision :: x__
  integer :: ix
  ! ----------------------------------------------------------------------------
  call mean(x__,x,nx)
  x_ = 0
  do ix=1,nx
    x_ = x_ + (x(ix)-x__)**2
  enddo
  x_ = dsqrt(x_ / dble(nx))
end subroutine std
! ------------------------------------------------------------------------------
subroutine window_mean(v,u,nt,nw)
  integer, intent(in) :: nt, nw
  double precision, intent(in) :: u(nt)
  double precision, intent(in out) :: v(nt)

  integer :: nb_, nb__, it, it_
  double precision :: v_
  ! ----------------------------------------------------------------------------
  if (mod(nw,2)==0) then
    nb_ = nw/2
    ! last block starts at (-1)
    nb__ = nt - nb_
  else
    nb_ = ceiling(dble(nw/2))
    ! last block starts at (-1)
    nb__ = nt - nb_ + 1
  endif
  ! ----------------------------------------------------------------------------
  ! starting block
  do it=1,(nb_-1)
    ! v(it) = (ones(1,it+nb_-1)/(it+nb_-1)) * u(1:(it+nb_-1));
    v_ = 0
    do it_=1,it+nb_-1
      v_ = v_ + u(it_)
    enddo
    v(it) = v_ / dble(it+nb_-1)
  enddo
  ! middle block
  do it=nb_,nb__
    ! v(it) = (ones(1,nw)/nw) * u((it-nb_+1):(it-nb_+nw));
    v_ = 0
    do it_=(it-nb_+1),(it-nb_+nw)
      v_ = v_ + u(it_)
    enddo
    v(it) = v_ / dble(nw)
  enddo
  ! ending block
  do it=(nb__+1),nt
    ! v(it) = (ones(1,(nt-it+nb_))/(nt-it+nb_)) * u((it-nb_+1):nt);
    v_ = 0
    do it_=(it-nb_+1),nt
      v_ = v_ + u(it_)
    enddo
    v(it) = v_ / dble(nt-it+nb_)
  enddo
end subroutine window_mean
! ------------------------------------------------------------------------------
end module calculus
