program linreg
 ! ----------------------------------------------------------------------
 ! 🎛️🖥️
 !
 !                           🟪🔴 = 🔷🔺
 !
 !
 ! Ax = b
 ! AP = QR
 ! -----------------------------------------------------------------------
 !
 ! ••• compiling in 💩
 !
 ! cmd.exe "/K" '"C:\Program Files (x86)\Intel\oneAPI\setvars.bat" && powershell'
 ! ifort /Qmkl /c "C:\Program Files (x86)\Intel\oneAPI\mkl\latest\include\lapack.f90" linreg.f90
 ! ifort /Qmkl linreg.obj lapack.obj
 ! .\linreg.exe
 !
 !
 ! ••• compiling in 🚀
 !
 ! source /opt/intel/oneapi/setvars.sh intel64
 ! rm -f *.o *.mod
 ! ifort -qmkl -c /opt/intel/oneapi/mkl/2023.0.0/include/lapack.f90 linreg.f90
 ! ifort -qmkl linreg.o lapack.o
 ! rm -f *.o *.mod
 ! mv a.out linreg.out
 ! ./linreg.out
 ! -----------------------------------------------------------------------
 ! AX = B
 !
 ! ?geqp3 ⟶ AP = QR
 ! ?ormqr ⟶ C  = Q.'B (for ℝ)
 ! unmqr  ⟶ C  = Q'B  (for ℂ)
 ! trsm   ⟶ X  = R\C
 !
 ! page 922, 927, 140 of oneapi-mkl.pdf
 ! -----------------------------------------------------------------------
 !
 !     _________     _____________     _____________
 !    |         |   |             |   |             |
 ! m  |   A     |   |     X       |   |      B      |
 !    |         | * |             | = |             | m
 !    |         |   |_____________|   |             |
 !    |_________|         p           |_____________|
 !         n                                 p
 !
 !
 ! ⟶ call dgeqp3(m, n, A, lda, jpvt, tau, work, lwork, info)
 !
 ! 📥
 ! A is of size m by n
 ! lda is a ℕ the leading dimension of A. at least max(1, m)
 ! jpvt is an array ℤ size at least max(1, n)
 ! work is an array (real, double, or cmplx) of size max(1, lwork)
 ! lwork must be at least:
 !                max(1, 3*n+1) ℝ
 !                max(1, n+1)   ℂ
 !
 ! 📤
 ! A
 ! tau is an array (real, double, or cmplx) of size at least max(1, min(m,n))
 ! jpvt the columns of AP are the columns of A like so:
 !                jpvt(1), ..., jpvt(n)
 ! info is an ℤ. 0 is good, -i the ith param is invalid.
 !
 ! -----------------------------------------------------------------------
 ! 👨🏻‍🏫
 ! lwork = -1 finds the best value for lwork and writes it in work(1).
 ! -----------------------------------------------------------------------
 ! matlab check
 !
 ! A = magic(3);
 ! [Q,R,P] = qr(A);
 ! b = [28 ; 34 ; 28];
 ! x = A\b;
 !
 ! x =
 ! 1
 ! 2
 ! 3
 !
 ! A =
 !     8     1     6
 !     3     5     7
 !     4     9     2
 !
 ! P =
 ! 0     1     0
 ! 1     0     0
 ! 0     0     1
 !
 ! Q =
 !    -0.0967    0.9912    0.0901
 !    -0.4834    0.0323   -0.8748
 !    -0.8701   -0.1281    0.4760
 !
 ! R =
 !   -10.3441   -5.7037   -5.7037
 !          0    7.5145    5.9176
 !          0         0   -4.6314
 ! -----------------------------------------------------------------------
 implicit none
 include 'mkl.fi'
 ! -----------------------------------------------------------------------
 integer, parameter :: n = 3, m = 3, p = 1

 double precision :: A(m,n)
 double precision :: b(m)
 integer :: ii

 ! dgeqp3
 integer :: lda
 integer :: jpvt(n)
 double precision, dimension(:), allocatable :: work
 double precision, dimension(:), allocatable :: tau ! size is min(m,n)
 integer :: lwork, info

 ! dormqr
 integer :: ldc
 ! k = "# of elementary reflectors whose product defines the matrix Q." wtf
 integer :: k = m
 double precision :: b_(m)


 ! dtrsm
 double precision, parameter :: alph=1
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
 !
 !
 !                    🎸 lets rock and roll 💃
 !
 !
 ! ---------------------------------------------------------------------------
 lda = m ! max(m,n)
 if (n>m) then
   lda=n
 endif
 ldc = m ! max(m,p)
 if (p>m) then
   ldc=p
 endif
 do ii=1,n
   jpvt(ii) = 0 ! 1 or 0
 enddo
 lwork = -1
 allocate(work(3*n+1))
 allocate(tau(n))

 ! lets factorize A as AP = QR
 ! find the best value for lwork
 call dgeqp3(m, n, A, lda, jpvt, tau, work, lwork, info)
 lwork=work(1)
 deallocate(work)
 allocate(work(lwork))
 ! now run again with optimal size
 call dgeqp3(m, n, A, lda, jpvt, tau, work, lwork, info)

 ! now A is factorized (in a weird way)
 ! time to multiply Q'b
 ! but first more optimal shit
 lwork=-1
 call dormqr('L','T', m, p, k, A, lda, tau, b, ldc, work, lwork, info)
 lwork=work(1)
 deallocate(work)
 allocate(work(lwork))
 call dormqr('L','T', m, p, k, A, lda, tau, b, ldc, work, lwork, info)

 ! ok, now the ending:
 ! Rx = Q'b
 ! x  = R \ (Q'b)
 call dtrsm('L', 'U', 'N', 'N', m, p, alph, A, lda, b, ldc)

 ! now unscramble Px
 do ii=1,m
  b_(ii) = b(ii)
 enddo
 do ii=1,m
  b(jpvt(ii)) = b_(ii)
 enddo
 ! ---------------------------------------------------------------------------
 ! 🧼
 ! ---------------------------------------------------------------------------
 deallocate(work)
 deallocate(tau)
 ! ---------------------------------------------------------------------------
 print *, ''
 print *, ' x'
 do ii=1,n
   write(*,*) ' ', b(ii)
 enddo
 print *, ''
 print *, ' jpvt'
 print *, jpvt
 print *, ''
 ! ---------------------------------------------------------------------------
end program linreg
