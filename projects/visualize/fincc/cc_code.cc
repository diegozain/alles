#include <iostream>
#include <cstdlib>
#include <dlfcn.h>

typedef void (*add_numbers_t)(double*, double*, double*);

void call_fortran_function() {
  // Load the Fortran shared library
  void* handle = dlopen("./libfortran.so", RTLD_LAZY);
  if (!handle) {
    std::cerr << "Cannot load library: " << dlerror() << '\n';
    exit(EXIT_FAILURE);
  }

  // Resolve the symbol
  dlerror(); // Clear any existing errors
  add_numbers_t add_numbers = (add_numbers_t)dlsym(handle, "add_numbers");
  const char* dlsym_error = dlerror();
  if (dlsym_error) {
    std::cerr << "Cannot load symbol 'add_numbers': " << dlsym_error << '\n';
    dlclose(handle);
    exit(EXIT_FAILURE);
  }

  // Use the Fortran function
  double a = 3.0, b = 4.0, result = 0.0;
  add_numbers(&a, &b, &result);

  std::cout << "Result: " << result << std::endl;

  // Close the shared library
  dlclose(handle);
}

// Export the function from the shared library
extern "C" void call_fortran() {
  call_fortran_function();
}

