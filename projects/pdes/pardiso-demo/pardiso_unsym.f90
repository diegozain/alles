program pardiso_unsym
! -----------------------------------------------------------------------
! diego domenzain @ AU 2021
!
! solve a sparse system of equations using pardiso: a_=b
!
! -----------------------------------------------------------------------
! for compiling in 💩
!
! $> cmd.exe "/K" '"C:\Program Files (x86)\Intel\oneAPI\setvars.bat" && powershell'
! $> ifort /Qmkl [name of file].f
! $> .\[name of file].exe
! -----------------------------------------------------------------------
! for compiling in 🚀
!
! source /opt/intel/oneapi/setvars.sh intel64
! ifort -qmkl pardiso_unsym.f90
! mv a.out pardiso_unsym.out
! ./pardiso_unsym.out
! -----------------------------------------------------------------------
implicit none
include 'mkl_pardiso.fi'
! -----------------------------------------------------------------------
! -- sparse matrix A
integer, parameter :: nnz=12, n=5
integer ja(nnz), ia_(n+1)
double precision a_(nnz)

! -- vector b
double precision b(n)

! -- vector x
! clever case when solution is written in b
! iparm(6) = 1 (dummy)
double precision x(1)
! ! normal case when solution is written in a new x
! ! iparm(6) = 0 (not dummy)
! double precision x(n)
! -----------------------------------------------------------------------
!
!                         🌴 pardiso stuff 🌴
!
! from onemkl-developerreference-fortran.pdf,
!                            chapter 'Sparse Solver Routines', page 2598,
! and also from pardiso_unsym.f.
!
! -----------------------------------------------------------------------
! solves X for the problem:
!
!     _________     _____________     _____________
!    |         |   |             |   |             |
! n  |   A     |   |     X       |   |      B      | n
!    |         | * |             | = |             |
!    |_________|   |_____________|   |_____________|
!         n              m                  m
!
! -----------------------------------------------------------------------
!
! pardiso(pt, maxfct, mnum, mtype, phase, n, a_, ia_, ja, perm, m, iparm, prt_msg, b, x, error)
!
! pt     :: internal solver memory pointer for 64 & 32 bit architectures.
! maxfct :: fancy memory storage option. keep it =1.
! mnum   :: idk, just keep it =1.
! mtype  :: type of matrix: real, complex, symmetric, etc.
!           =11 real    and nonsymmetric ⬜
!           =13 complex and nonsymmetric 🔲
! phase  :: controls the execution of the solver. keep it as shown here.
! n      :: integer, size of matrix A (n × n)
! a_     :: array. value entries of A (nnz)
! ia_    :: array. index row entries of the CSR vector.
! ja     :: array. index column entries.
! perm   :: very complicated. here it is just a dummy array of size 1.
! m      :: size of x and b (n × m).
! iparam :: 😵 VERY complicated. see manual, page 2662. 🚑
! prt_msg:: =1 print messages from pardiso. =0 no print.
! b      :: values of b of size (n × m).
! x      :: initialize the solution x of size (n × m).
!           if iparm(6) = 0 then solution is written in b! 🔥
! error  :: error flag.
!
! -----------------------------------------------------------------------
type(MKL_PARDISO_HANDLE) pt(64)
integer :: m=1, maxfct=1, mnum=1, prt_msg=0, mtype=11, error=0
integer iparm(64), idummy(1)
double precision  ddum(1)
integer i, phase
! -----------------------------------------------------------------------
! initialize iparam.
!
! 😵 VERY complicated. see manual. 🚑
!
! page 2622 of onemkl-developerreference-fortran.pdf
! -----------------------------------------------------------------------
do i = 1, 64
    iparm(i) = 0
end do

iparm(1) = 1 ! no solver default
iparm(2) = 2 ! fill-in reordering from METIS
iparm(3) = 1 ! numbers of processors
iparm(4) = 0 ! no iterative-direct algorithm
iparm(5) = 0 ! no user fill-in reducing permutation
! iparm(6) = 0 ! solution on the first n components of x
iparm(6) = 1 ! solution x is written in b... less memory 🔥
iparm(8) = 9 ! numbers of iterative refinement steps
iparm(10)= 13 ! perturb the pivot elements with 1E-13
iparm(11)= 1 ! use nonsymmetric permutation and scaling MPS
iparm(13)= 1 ! maximum weighted matching algorithm is ON

! iparm(14)  Output: number of perturbed pivots
! iparm(15)  Output: Peak memory on symbolic factorization
! iparm(17)  Output: Peak memory on numerical factorization and solution.
! iparm(18)  Output: number of nonzeros in the factor LU
! iparm(19)  Output: Mflops for LU factorization
! iparm(20)  Output: Numbers of CG Iterations
! -----------------------------------------------------------------------
!
!
!                    🎸 lets rock and roll 💃
!
!
! -----------------------------------------------------------------------
! fill A & b
ia_= [ 1,3,6,9,10,13 ]
ja = [ 1,2,1,3,5,2,3,4,3,2,3,5 ]
a_ = [ 2,3,3,4,6,-1,-3,2,1,4,2,1 ]

b = [ 8,45,-3,3,19 ]
! -----------------------------------------------------------------------
! print funny stuff so user is amused
write(*,*) 'ia_= '
do i = 1, n+1
    write(*,*) ' ', ia_(i)
end do

write(*,*) 'ja= '
do i = 1, nnz
    write(*,*) ' ', ja(i)
end do

write(*,*) 'a_= '
do i = 1, nnz
    write(*,*) ' ', a_(i)
end do

write(*,*) 'b = '
do i = 1, n
    write(*,*) ' ', b(i)
end do
! -----------------------------------------------------------------------
! Initialize the internal solver memory pointer. This is only
! necessary for the FIRST call of the PARDISO solver.
do i = 1, 64
    pt(i)%DUMMY = 0
end do
! -----------------------------------------------------------------------
! Reordering and Symbolic Factorization, This step also allocates
!   all memory that is necessary for the factorization
phase = 11 ! only reordering and symbolic factorization
call pardiso(pt, maxfct, mnum, mtype, phase, n, a_, ia_, ja, idummy, m, iparm, prt_msg, ddum, ddum, error)
! -----------------------------------------------------------------------
! Factorization.
phase = 22 ! only factorization
call pardiso(pt, maxfct, mnum, mtype, phase, n, a_, ia_, ja,idummy, m, iparm, prt_msg, ddum, ddum, error)
! -----------------------------------------------------------------------
! Back substitution and iterative refinement
iparm(8) = 2 ! max numbers of iterative refinement steps
phase    = 33 ! only factorization
call pardiso(pt, maxfct, mnum, mtype, phase, n, a_, ia_, ja,idummy, m, iparm, prt_msg, b, x, error)
! -----------------------------------------------------------------------
! print solution to terminal
write(*,*) 'x = '
do i = 1, n
    write(*,*) ' ', b(i)
end do
! -----------------------------------------------------------------------
! Termination and release of memory
phase = -1 ! release internal memory
call pardiso(pt, maxfct, mnum, mtype, phase, n, ddum, idummy,idummy, idummy, m, iparm, prt_msg, ddum, ddum, error)
! -----------------------------------------------------------------------
end program pardiso_unsym
