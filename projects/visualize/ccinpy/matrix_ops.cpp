#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <iostream>

namespace py = pybind11;

// Scale a NumPy array by a scalar
void scale_array(py::array_t<double> input_array, double scalar) {
    // Request buffer information
    py::buffer_info buf = input_array.request();

    // Ensure the input is a 1D array
    if (buf.ndim != 1)
        throw std::runtime_error("Input array must be 1-dimensional");

    // Get a pointer to the data
    double* ptr = static_cast<double*>(buf.ptr);

    // Apply scalar multiplication
    for (size_t i = 0; i < buf.shape[0]; i++) {
        ptr[i] *= scalar;
    }
}

// Expose the function to Python
PYBIND11_MODULE(matrix_ops, m) {
    m.def("scale_array", &scale_array, "Scale a NumPy array by a scalar");
}

