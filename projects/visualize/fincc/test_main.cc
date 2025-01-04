#include <iostream>

// Declare the function provided by the C++ shared library
extern "C" void call_fortran();

int main() {
  std::cout << "Calling Fortran function from C++ shared library..." << std::endl;
  call_fortran();
  return 0;
}

