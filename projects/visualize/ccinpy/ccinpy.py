import numpy as np
import matrix_ops

arr = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64)
print("Before scaling:", arr)

matrix_ops.scale_array(arr, 2.0)

print("After scaling:", arr)

