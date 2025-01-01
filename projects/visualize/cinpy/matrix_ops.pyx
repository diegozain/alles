# Import required modules
import numpy as np
cimport numpy as cnp

# Declare the C function
cdef extern void scale_array(double* arr, size_t size, double scalar)

# Define a Python-accessible function
def scale_numpy_array(cnp.ndarray[cnp.float64_t, ndim=1] arr, double scalar):
  if not arr.flags['C_CONTIGUOUS']:
    raise ValueError("Array must be C-contiguous")

  # Get the size of the array
  cdef size_t size = arr.shape[0]

  # Call the C function
  scale_array(<double*> arr.data, size, scalar)
