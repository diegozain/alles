#include <stdio.h>
#include <stdlib.h>
#include <dlfcn.h>

typedef void (*add_numbers_t)(double*, double*, double*);

void call_fortran_function(double a, double b) {
  // Load the Fortran shared library
  void* handle = dlopen("./libfortran.so", RTLD_LAZY);
  if (!handle) {
    fprintf(stderr, "Cannot load library: %s\n", dlerror());
    exit(EXIT_FAILURE);
  }

  // Resolve the symbol for the Fortran function
  add_numbers_t add_numbers = (add_numbers_t)dlsym(handle, "add_numbers");
  const char* dlsym_error = dlerror();
  if (dlsym_error) {
    fprintf(stderr, "Cannot load symbol 'add_numbers': %s\n", dlsym_error);
    dlclose(handle);
    exit(EXIT_FAILURE);
  }

  // Call the Fortran function
  double result = 0.0;
  add_numbers(&a, &b, &result);

  printf("Result: %.2f\n", result);

  // Close the shared library
  dlclose(handle);
}

// Exported function for external use
void call_fortran(double a, double b) {
  call_fortran_function(a, b);
}

