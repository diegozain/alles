import numpy as np
import scale_array

# Create a NumPy array
arr = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64)

print("Before scaling:", arr)

# Call the Fortran function
scale_array.scale_array(arr, arr.size, 2.0)

print("After scaling:", arr)

