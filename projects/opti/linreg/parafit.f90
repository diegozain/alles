program parafit
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
 ! ifort /Qmkl /c "C:\Program Files (x86)\Intel\oneAPI\mkl\latest\include\lapack.f90" parafit.f90
 ! ifort /Qmkl parafit.obj lapack.obj
 ! .\parafit.exe
 !
 !
 ! ••• compiling in 🚀
 !
 ! source /opt/intel/oneapi/setvars.sh intel64
 ! rm -f *.o *.mod
 ! ifort -qmkl -c /opt/intel/oneapi/mkl/2023.0.0/include/lapack.f90 parafit.f90
 ! ifort -qmkl parafit.o lapack.o
 ! rm -f *.o *.mod
 ! mv a.out parafit.out
 ! ./parafit.out
 ! -----------------------------------------------------------------------
 !
 !    fit a parabola to some points.
 !
 ! -----------------------------------------------------------------------
 ! matlab check
 !
 ! x = [0;1;3;4;5];
 ! y = [5;2;10;20;50];
 !
 ! nd= numel(x);
 ! P = [x.^2 x ones(nd,1)];
 !
 ! P'*P
 ! P.'*y
 ! b = (P'*P)\(P.'*y)
 !
 ! xx=linspace(min(x)*1.1,max(x)*1.1,100);
 ! yy=b(1)*xx.^2 + b(2)*xx + b(3);
 !
 ! figure;
 ! hold on;
 ! plot(xx,yy,'linewidth',3)
 ! plot(x,y,'.','markersize',20)
 ! hold off;
 ! axis tight;
 ! axis square;
 ! -----------------------------------------------------------------------
 implicit none
 include 'mkl.fi'
 ! -----------------------------------------------------------------------
 integer, parameter :: nd=5

 double precision :: P(nd,3), x(nd), y(nd)
 double precision :: A(3,3)
 double precision :: b(3)
 integer :: ii
 double precision :: xmin

 ! dgeqp3
 integer :: jpvt(3)
 double precision, dimension(:), allocatable :: work
 double precision :: tau(3)
 integer :: lwork, info
 ! dtrsm & dgemm
 double precision, parameter :: alph=1, beta=0
 double precision :: C(1)
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
 !
 !
 !                    🎸 lets rock and roll 💃
 !
 !
 ! ---------------------------------------------------------------------------
 ! 👷👷
 ! ---------------------------------------------------------------------------
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
 ! ---------------------------------------------------------------------------
 ! 🟪🔴 = 🔷🔺
 ! --------------------------------------------------------------------------
 do ii=1,3
   jpvt(ii) = 1
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
 ! ---------------------------------------------------------------------------
 ! 🧼
 ! ---------------------------------------------------------------------------
 deallocate(work)
 ! ---------------------------------------------------------------------------
 ! ⭐🌟🌠
 ! ---------------------------------------------------------------------------
 xmin = - b(2) / (2*b(1))
 ! ---------------------------------------------------------------------------
 print *, ''
 print *, ' the parabola ax² + bx + c is given by'
 do ii=1,3
   write(*,*) ' ', b(ii)
 enddo
 print *, ''
 print *, ' minimum of the parabola is achieved at this x'
 print *, xmin
 print *, ''
 ! ---------------------------------------------------------------------------
end program parafit
