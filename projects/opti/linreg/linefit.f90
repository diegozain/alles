program linefit
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
 ! ifort /Qmkl /c "C:\Program Files (x86)\Intel\oneAPI\mkl\latest\include\lapack.f90" linefit.f90
 ! ifort /Qmkl linefit.obj lapack.obj
 ! .\linefit.exe
 !
 !
 ! ••• compiling in 🚀
 !
 ! source /opt/intel/oneapi/setvars.sh intel64
 ! rm -f *.o *.mod
 ! ifort -qmkl -c /opt/intel/oneapi/mkl/2023.0.0/include/lapack.f90 linefit.f90
 ! ifort -qmkl linefit.o lapack.o
 ! rm -f *.o *.mod
 ! mv a.out linefit.out
 ! ./linefit.out
 ! -----------------------------------------------------------------------
 !
 !    fit a line to some points.
 !
 ! -----------------------------------------------------------------------
 ! matlab check
 !
 ! x = [0;1;3;4;5];
 ! y = [5;2;10;20;50];
 
 ! nd= numel(x);
 ! P = [x ones(nd,1)];
 
 ! P'*P
 ! P.'*y
 ! b = (P'*P)\(P.'*y)
 
 ! xx=linspace(min(x)*1.1,max(x)*1.1,100);
 ! yy=b(1)*xx + b(2);
 
 ! figure;
 ! hold on;
 ! plot(xx,yy,'linewidth',2)
 ! plot(x,y,'.','markersize',20)
 ! hold off;
 ! axis tight;
 ! axis square;
 ! -----------------------------------------------------------------------
 implicit none
 include 'mkl.fi'
 ! -----------------------------------------------------------------------
 integer, parameter :: nd=5

 double precision :: P(nd,2), x(nd), y(nd)
 double precision :: A(2,2)
 double precision :: b(2)
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
   P(ii,1) = x(ii)
   P(ii,2) = 1
 enddo
 ! 👷 A ⟵ P'P     pg. 119
 call dgemm('T', 'N', 2, 2, nd, alph, P, nd, P, nd, beta, A, 2)
 ! 👷 b ⟵ P'y     pg. 73
 call dgemv('T', nd, 2, alph, P, nd, y, 1, beta, b, 1)
 ! ---------------------------------------------------------------------------
 ! 🟪🔴 = 🔷🔺
 ! --------------------------------------------------------------------------
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
 ! ---------------------------------------------------------------------------
 ! 🧼
 ! ---------------------------------------------------------------------------
 deallocate(work)
 ! ---------------------------------------------------------------------------
 print *, ''
 print *, ' the line ax + b is given by'
 do ii=1,2
   write(*,*) ' ', b(ii)
 enddo
 print *, ''
 ! ---------------------------------------------------------------------------
end program linefit
