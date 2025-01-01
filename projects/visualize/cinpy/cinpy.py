import numpy as np
from matrix_ops import scale_numpy_array

arr = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64)
print("Before scaling:", arr)

scale_numpy_array(arr, 2.0)

print("After scaling:", arr)