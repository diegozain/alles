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
 ! compiling in 💩
 !
 ! $> cmd.exe "/K" '"C:\Program Files (x86)\Intel\oneAPI\setvars.bat" && powershell'
 ! $> ifort /Qmkl [name of file].f
 ! $> .\[name of file].exe
 ! -----------------------------------------------------------------------
 ! compiling in 🚀
 !
 ! source /opt/intel/oneapi/setvars.sh intel64
 ! rm -f *.o *.mod
 ! ifort -qmkl -c /opt/intel/oneapi/mkl/2023.0.0/include/lapack.f90 linreg.f90
 ! ifort -qmkl lapack.o linreg.o
 ! rm -f *.o *.mod
 ! mv a.out linreg.out
 ! ./linreg.out
 ! -----------------------------------------------------------------------
 ! AX = B 
 !
 ! ?geqp3 ⟶ AP = QR
 ! ormqr  ⟶ C  = Q.'B (for ℝ)
 ! unmqr  ⟶ C  = Q'B  (for ℂ)
 ! trsm   ⟶ X  = R\C
 !
 ! page 922 of oneapi-mkl.pdf
 ! -----------------------------------------------------------------------
 !
 !     _________     _____________     _____________
 !    |         |   |             |   |             |
 ! m  |   A     |   |     X       |   |      B      | 
 !    |         | * |             | = |             | m
 !    |         |   |_____________|   |             |
 !    |_________|         m           |_____________|
 !         n                                 m
 !
 !
 ! ⟶ call dgeqp3(m, n, A, lda, jpvt, tau, work, lwork, info)
 !
 ! 📥
 ! lda is a ℕ the leading dimension of A. at least max(1, m)
 ! jpvt is an array (integer) size at least max(1, n)
 ! work is an array (real, double, or cmplx) of size max(1, lwork)
 ! lwork must be at least:
 !                max(1, 3*n+1) ℝ
 !                max(1, n+1)   ℂ
 ! 
 ! 📤
 ! A
 ! tau is an array of size at least max(1, min(m,n))
 ! jvpt the columns of AP are the columns of A like so:
 !                jvpt(1), ..., jvpt(n)
 ! info is an ℤ. 0 is good, -i the ith param is invalid.
 ! 
 ! -----------------------------------------------------------------------
 ! lwork = -1 finds the best value for lwork and writes it in work(1).
 ! -----------------------------------------------------------------------
 ! matlab check
 !
 ! A = magic(3);
 ! [Q,R,P] = qr(A);
 ! b = [1 ; 2 ; 3];
 ! x = A\b;
 !
 ! x = 
 ! 0.05
 ! 0.3
 ! 0.05
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
 
 ! --------------------------------------------------------------------------
 !
 !
 !                    🎸 lets rock and roll 💃
 !
 !
 ! ---------------------------------------------------------------------------
 ! ! find the best value for lwork
 ! call dgeqp3(m, n, A, lda, jpvt, tau, work, lwork, info)
 ! ! now run again with optimal size
 ! call dgeqp3(m, n, A, lda, jpvt, tau, work, lwork, info)
 ! ---------------------------------------------------------------------------

 ! ---------------------------------------------------------------------------
 ! ---------------------------------------------------------------------------
end program linreg