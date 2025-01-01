module scale_array_mod
    implicit none
contains
    subroutine scale_array(arr, n, scalar)
        real(8), intent(inout) :: arr(:)   ! Accepts a 1D array
        integer, intent(in) :: n          ! Size of the array
        real(8), intent(in) :: scalar     ! Scalar value
        integer :: i

        do i = 1, n
            arr(i) = arr(i) * scalar
        end do
    end subroutine scale_array
end module scale_array_mod

