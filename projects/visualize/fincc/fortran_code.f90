! fortran_code.f90
module fortran_module
  use, intrinsic :: iso_c_binding ! ⟶ c_int, c_double, c_double_complex
  implicit none
contains
  subroutine add_numbers(a, b, result) bind(c, name="add_numbers")
    real(c_double), intent(in)  :: a, b
    real(c_double), intent(out) :: result
    result = a + b
  end subroutine add_numbers
end module fortran_module

