#include <stddef.h>

// gcc -shared -o matrix_ops.so -fPIC matrix_ops.c

void scale_array(double* arr, size_t size, double scalar) {
  for (size_t i = 0; i < size; i++) {
    arr[i] *= scalar;
  }
}
