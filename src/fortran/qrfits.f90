module qrfits
 ! -----------------------------------------------------------------------------
 ! diego domenzain
 !
 ! 🟪🔴 = 🔷🔺 wrapper of qr routines in lapack-oneapi-mkl
 ! -----------------------------------------------------------------------------
 include 'mkl.fi'
 contains
 ! -----------------------------------------------------------------------------
 subroutine linefit(nd,x,y,b)
  integer, intent(in) :: nd
  double precision, intent(in) :: x(nd), y(nd)
  double precision, intent(in out) :: b(2)

  double precision :: P(nd,2)
  double precision :: A(2,2)
  integer :: ii
  ! dtrsm & dgemm
  double precision, parameter :: alph=1, beta=0

  ! 👷 P
  do ii=1,nd
   P(ii,1) = x(ii)
   P(ii,2) = 1
  enddo
  ! 👷 A ⟵ P'P     pg. 119
  call dgemm('T', 'N', 2, 2, nd, alph, P, nd, P, nd, beta, A, 2)
  ! 👷 b ⟵ P'y     pg. 73
  call dgemv('T', nd, 2, alph, P, nd, y, 1, beta, b, 1)
  ! 🟪🔴 = 🔷🔺
  call linreg(2,A,b)
 end subroutine linefit
 ! -----------------------------------------------------------------------------
 subroutine parafit(nd,x,y,b)
  ! --------------------------------------------------------------------------
  ! ⭐🌟🌠
  ! remember that
  !               xmin = - b(2) / (2*b(1))
  ! --------------------------------------------------------------------------
  integer, intent(in) :: nd
  double precision, intent(in) :: x(nd), y(nd)
  double precision, intent(in out) :: b(3)

  double precision :: P(nd,3)
  double precision :: A(3,3)
  integer :: ii
  ! dtrsm & dgemm
  double precision, parameter :: alph=1, beta=0

  ! 👷 P
  do ii=1,nd
   P(ii,1) = x(ii)**2
   P(ii,2) = x(ii)
   P(ii,3) = 1
  enddo
  ! 👷 A ⟵ P'P     pg. 119
  call dgemm('T', 'N', 3, 3, nd, alph, P, nd, P, nd, beta, A, 3)
  ! 👷 b ⟵ P'y     pg. 73
  call dgemv('T', nd, 3, alph, P, nd, y, 1, beta, b, 1)
  ! 🟪🔴 = 🔷🔺
  call linreg(3,A,b)
 end subroutine parafit
 ! -----------------------------------------------------------------------------
 subroutine linreg(n,A,b)
   integer, intent(in) :: n
   double precision, intent(in) :: A(n,n)
   double precision, intent(in out) :: b(n)

   ! 🟪🔴 = 🔷🔺
   ! dgeqp3
   integer :: jpvt(n)
   double precision :: b_(n)
   double precision, dimension(:), allocatable :: work
   double precision :: tau(n)
   integer :: lwork, info
   ! dtrsm
   double precision, parameter :: alph=1

   ! 🟪🔴 = 🔷🔺
   do ii=1,n
     jpvt(ii) = 0
   enddo
   lwork = -1
   allocate(work(3*n+1))
   ! AP = QR
   call dgeqp3(n, n, A, n, jpvt, tau, work, lwork, info)
   lwork=work(1)
   deallocate(work)
   allocate(work(lwork))
   call dgeqp3(n, n, A, n, jpvt, tau, work, lwork, info)
   ! b ⟵ Q'b
   lwork=-1
   call dormqr('L','T', n, 1, n, A, n, tau, b, n, work, lwork, info)
   lwork=work(1)
   deallocate(work)
   allocate(work(lwork))
   call dormqr('L','T', n, 1, n, A, n, tau, b, n, work, lwork, info)
   ! Rx = Q'b
   ! x  = R \ (Q'b)
   call dtrsm('L', 'U', 'N', 'N', n, 1, alph, A, n, b, n)
   ! x ⟵ Px
   do ii=1,n
    b_(ii) = b(ii)
   enddo
   do ii=1,n
    b(jpvt(ii)) = b_(ii)
   enddo
   ! 🧼
   deallocate(work)
 end subroutine linreg
 ! -----------------------------------------------------------------------------
 subroutine linreg_(m,n,A,b,Ab)
  integer, intent(in) :: m,n
  double precision, intent(in) :: A(m,n), b(m)
  double precision, intent(in out) :: Ab(n)

  ! dgemm & dgemv
  double precision, parameter :: alph=1, beta=0
  double precision :: AA(n,n)

  ! 👷 AA ⟵ A'A     pg. 119
  call dgemm('T', 'N', n, n, m, alph, A, m, A, m, beta, AA, n)
  ! 👷 Ab ⟵ A'b     pg. 73
  call dgemv('T', m, n, alph, A, m, b, 1, beta, Ab, 1)
  ! 🟪🔴 = 🔷🔺
  call linreg(n,AA,Ab)
 end subroutine linreg_
 ! -----------------------------------------------------------------------------
end module qrfits
