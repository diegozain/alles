program qrfits_ie
 ! ----------------------------------------------------------------------
 ! 🎛️🖥️
 !
 !                           🟪🔴 = 🔷🔺
 !
 !
 ! Ax = b
 ! AP = QR
 ! -----------------------------------------------------------------------
 use qrfits
 implicit none
 ! -----------------------------------------------------------------------
 integer, parameter :: nd=5
 double precision :: x(nd), y(nd)
 double precision :: bl(2), bp(3)

 double precision :: A(3,3), b(3)
 double precision :: C(4,3), d(4), Cd(3)

 integer :: ii
 ! --------------------------------------------------------------------------
 x(1) = 0
 x(2) = 1
 x(3) = 3
 x(4) = 4
 x(5) = 5

 y(1) = 5
 y(2) = 2
 y(3) = 10
 y(4) = 20
 y(5) = 50
 ! --------------------------------------------------------------------------
 A(1,1) = 8
 A(2,1) = 3
 A(3,1) = 4

 A(1,2) = 1
 A(2,2) = 5
 A(3,2) = 9

 A(1,3) = 6
 A(2,3) = 7
 A(3,3) = 2

 b(1) = 28
 b(2) = 34
 b(3) = 28
 ! --------------------------------------------------------------------------
 C(1,1) = 8
 C(2,1) = 3
 C(3,1) = 4
 C(4,1) = 1

 C(1,2) = 1
 C(2,2) = 5
 C(3,2) = 9
 C(4,2) = 2

 C(1,3) = 6
 C(2,3) = 7
 C(3,3) = 2
 C(4,3) = 9

 d(1) = 28
 d(2) = 34
 d(3) = 28
 d(4) = 31
 ! --------------------------------------------------------------------------
 !
 !
 !                    🎸 lets rock and roll 💃
 !
 !
 ! ---------------------------------------------------------------------------
 call linefit(nd,x,y,bl)
 call parafit(nd,x,y,bp)
 ! ⭐🌟🌠
 ! remember that: xmin = - bp(2) / (2*bp(1))
 call linreg(3,A,b)
 call linreg_(4,3,C,d,Cd)
 ! ---------------------------------------------------------------------------
 print *, ''
 print *, ' the line ax + b is given by'
 do ii=1,2
   write(*,*) ' ', bl(ii)
 enddo
 print *, ''
 print *, ' the parabola ax² + bx + c is given by'
 do ii=1,3
   write(*,*) ' ', bp(ii)
 enddo
 print *, ''
 print *, ' x = A\b'
 do ii=1,3
   write(*,*) ' ', b(ii)
 enddo
 print *, ''
 print *, ' x = C\d'
 do ii=1,3
   write(*,*) ' ', Cd(ii)
 enddo
 print *, ''
 ! ---------------------------------------------------------------------------
end program qrfits_ie