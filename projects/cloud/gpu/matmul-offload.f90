! ==============================================================
! Copyright ©Intel Corporation
!
! SPDX-License-Identifier: MIT
! =============================================================
!
! Sample Matrix Multiply using OpenMP offload using 'target' and 'map'
! compile with ifx -qopenmp -fopenmp-targets=spir64 -xhost -o matmul-offload matmul-offload.f90
!
! =============================================================
! see here for more insight:
!
! https://www.intel.com/content/www/us/en/develop/documentation/get-started-with-cpp-fortran-compiler-openmp/top.html
!
! https://enccs.github.io/openmp-gpu/introduction/
!
! =============================================================
program matrix_multiply
use omp_lib
implicit none
integer, parameter :: N=1000
integer :: i, j, k, my_thread_id
real, allocatable, dimension(:,:) :: a, b, c, c_validate

    allocate( a(N,N), b(N,N), c(N,N), c_validate(N,N))

    !... Initialize a,b matrices to some values
    do j=1,N
        do i=1,N
            a(i,j) = i - j + 6.0
            b(i,j) = i - j + 7.0
        enddo
    enddo
    !... initialize the output matrix c
    !... initialize c_validate to hold expected values for validation of c
    c = 0.0
    c_validate = 0.0

    !... offload data and compute the matrix multiply on the GPU
    !... note we send 'a' and 'b' but do not move them back (no change)
    !... but 'c' needs to go to the GPU and brought back from GPU (changed)
    !$omp target map(to: a, b ) map(tofrom: c )
    !$omp parallel do private(j,i,k)
    do j=1,N
        do i=1,N
            do k=1,N
                c(i,j) = c(i,j) + a(i,k) * b(k,j)
            enddo
        enddo
    enddo
    !$omp end parallel do
    !$omp end target

    !... on host compute matrix multiplication with results in array 'c_validate'
    do j=1,N
        do i=1,N
            do k=1,N
                c_validate(i,j) = c_validate(i,j) + a(i,k) * b(k,j)
            enddo
        enddo
    enddo

    !... verify device computed array 'c' matches host computed array 'c_validate'
    do j=1,N
        do i=1,N
            if (c_validate(i,j) .ne. c(i,j)) then
            write(*,*) " VALIDATION FAILED, i, j, c_validate(i,j), c(i,j) ",&
                           i, j, c_validate(i,j), c(i,j)
            stop
            endif
        enddo
    enddo

    write(*,*) " VALIDATION PASSED"

end program matrix_multiply
