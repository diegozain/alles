#include <stdio.h>

// Declare the function provided by the C shared library
void call_fortran(double a, double b);

int main() {
  printf("Calling Fortran function from C shared library...\n");
  call_fortran(3.0, 4.0);
  return 0;
}

