program linreg_
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
 ! ifort /Qmkl /c "C:\Program Files (x86)\Intel\oneAPI\mkl\latest\include\lapack.f90" linreg_.f90
 ! ifort /Qmkl linreg_.obj lapack.obj
 ! .\linreg_.exe
 !
 !
 ! ••• compiling in 🚀
 !
 ! source /opt/intel/oneapi/setvars.sh intel64
 ! rm -f *.o *.mod
 ! ifort -qmkl -c /opt/intel/oneapi/mkl/2023.0.0/include/lapack.f90 linreg_.f90
 ! ifort -qmkl linreg_.o lapack.o
 ! rm -f *.o *.mod
 ! mv a.out linreg_.out
 ! ./linreg_.out
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
 !     _________     ____    ____
 !    |         |   |   |   |   |
 ! n  |   A     |   | x |   | b |
 !    |         | * |   | = |   | n
 !    |_________|   |___|   |___|
 !         n
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
 ! jvpt the columns of AP are the columns of A like so:
 !                jvpt(1), ..., jvpt(n)
 ! info is an ℤ. 0 is good, -i the ith param is invalid.
 !
 ! -----------------------------------------------------------------------
 ! 👨🏻‍🏫
 ! lwork = -1 finds the best value for lwork and writes it in work(1).
 ! -----------------------------------------------------------------------
 ! matlab check
 !
 !
 ! -----------------------------------------------------------------------
 use omp_lib, only : omp_get_wtime
 implicit none
 include 'mkl.fi'
 ! -----------------------------------------------------------------------
 integer, parameter :: n = 1019 ! <1019 ⇒ no bueno.
 double precision :: A(n,n), b(n)
 integer :: ii, jj
 real :: r

 ! 🟪🔴 = 🔷🔺
 ! dgeqp3
 integer :: jpvt(n)
 double precision, dimension(:), allocatable :: work
 double precision :: tau(n)
 integer :: lwork, info
 ! dtrsm
 double precision, parameter :: alph=1
 double precision :: b_(n)

 ! ⌚
 real :: start_time, end_time
 ! --------------------------------------------------------------------------
 do ii=1,n
   do jj=1,n
     call random_number(r)
     A(ii,jj) = dble(r)
   enddo
   b(ii) = 0
 enddo
 b(1) = 1
 ! --------------------------------------------------------------------------
 ! print *, ''
 ! print *, 'A'
 ! print *, A
 ! print *, ''
 ! print *, 'b'
 ! print *, b
 ! --------------------------------------------------------------------------
 !
 !
 !                    🎸 lets rock and roll 💃
 !
 !
 ! ---------------------------------------------------------------------------
 print *,''
 ! ⌚
 start_time = omp_get_wtime()
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
 print *,'   😃  just did AP = QR'
 ! b ⟵ Q'b
 lwork=-1
 call dormqr('L','T', n, 1, n, A, n, tau, b, n, work, lwork, info)
 lwork=work(1)
 deallocate(work)
 allocate(work(lwork))
 call dormqr('L','T', n, 1, n, A, n, tau, b, n, work, lwork, info)
 print *,'   😃  just did b  ⟵   Qᵀb'
 ! Rx = Q'b
 ! x  = R \ (Q'b)
 call dtrsm('L', 'U', 'N', 'N', n, 1, alph, A, n, b, n)
 print *,'   😃  just did Rx = Qᵀb'
 ! x ⟵ Px
 do ii=1,n
  b_(ii) = b(ii)
 enddo
 do ii=1,n
  b(jpvt(ii)) = b_(ii)
 enddo
 ! ---------------------------------------------------------------------------
 ! 🧼
 ! ---------------------------------------------------------------------------
 deallocate(work)
 ! ---------------------------------------------------------------------------
 ! ⌚
 end_time = omp_get_wtime()
 print *, ''
 print *, 'elapsed time: ', (end_time - start_time), 'seconds'
 print *, ''
 ! ---------------------------------------------------------------------------
 ! print *, ''
 ! print *, 'x'
 ! print *, b
 ! print *, ''
 ! ---------------------------------------------------------------------------
end program linreg_
