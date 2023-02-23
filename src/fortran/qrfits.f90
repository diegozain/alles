module qrfits
 ! ------------------------------------------------------------------------------
 ! diego domenzain
 !
 ! 🟪🔴 = 🔷🔺 wrapper of qr routines in lapack-oneapi-mkl
 ! ------------------------------------------------------------------------------
 include 'mkl.fi'
 contains
 ! ------------------------------------------------------------------------------
 subroutine linefit(nd,x,y,b)
  integer, intent(in) :: nd
  double precision, intent(in) :: x(nd), y(nd)
  double precision, intent(in out) :: b(2)

  double precision :: P(nd,2)
  double precision :: A(2,2)
  integer :: ii
  ! dgeqp3
  integer :: jpvt(2)
  double precision :: b_(2)
  double precision, dimension(:), allocatable :: work
  double precision :: tau(2)
  integer :: lwork, info
  ! dtrsm & dgemm
  double precision, parameter :: alph=1, beta=0
  double precision :: C(1)

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
  do ii=1,2
   jpvt(ii) = 0
  enddo
  lwork = -1
  allocate(work(2*2+1))
  ! AP = QR
  call dgeqp3(2, 2, A, 2, jpvt, tau, work, lwork, info)
  lwork=work(1)
  deallocate(work)
  allocate(work(lwork))
  call dgeqp3(2, 2, A, 2, jpvt, tau, work, lwork, info)
  lwork=-1
  ! b ⟵ Q'b
  call dormqr('L','T', 2, 1, 2, A, 2, tau, b, 2, work, lwork, info)
  lwork=work(1)
  deallocate(work)
  allocate(work(lwork))
  call dormqr('L','T', 2, 1, 2, A, 2, tau, b, 2, work, lwork, info)
  ! x = R \ Q'b
  call dtrsm('L', 'U', 'N', 'N', 2, 1, alph, A, 2, b, 2)
  ! x ⟵ Px
  do ii=1,2
   b_(ii) = b(ii)
  enddo
  do ii=1,2
   b(jpvt(ii)) = b_(ii)
  enddo
  ! 🧼
  deallocate(work)
 end subroutine linefit
 ! ------------------------------------------------------------------------------
 subroutine parafit(nd,x,y,b)
  ! ---------------------------------------------------------------------------
  ! ⭐🌟🌠
  ! remember that
  !               xmin = - b(2) / (2*b(1))
  ! ---------------------------------------------------------------------------
  integer, intent(in) :: nd
  double precision, intent(in) :: x(nd), y(nd)
  double precision, intent(in out) :: b(3)

  double precision :: P(nd,3)
  double precision :: A(3,3)
  integer :: ii
  ! dgeqp3
  integer :: jpvt(3)
  double precision :: b_(3)
  double precision, dimension(:), allocatable :: work
  double precision :: tau(3)
  integer :: lwork, info
  ! dtrsm & dgemm
  double precision, parameter :: alph=1, beta=0
  double precision :: C(1)

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
  do ii=1,3
   jpvt(ii) = 0
  enddo
  lwork = -1
  allocate(work(3*3+1))
  ! AP = QR
  call dgeqp3(3, 3, A, 3, jpvt, tau, work, lwork, info)
  lwork=work(1)
  deallocate(work)
  allocate(work(lwork))
  call dgeqp3(3, 3, A, 3, jpvt, tau, work, lwork, info)
  lwork=-1
  ! b ⟵ Q'b
  call dormqr('L','T', 3, 1, 3, A, 3, tau, b, 3, work, lwork, info)
  lwork=work(1)
  deallocate(work)
  allocate(work(lwork))
  call dormqr('L','T', 3, 1, 3, A, 3, tau, b, 3, work, lwork, info)
  ! x = R \ Q'b
  call dtrsm('L', 'U', 'N', 'N', 3, 1, alph, A, 3, b, 3)
  ! x ⟵ Px
  do ii=1,3
   b_(ii) = b(ii)
  enddo
  do ii=1,3
   b(jpvt(ii)) = b_(ii)
  enddo
  ! 🧼
  deallocate(work)
 end subroutine parafit
 ! ------------------------------------------------------------------------------
 ! subroutine linreg()

 ! end subroutine linreg
end module qrfits