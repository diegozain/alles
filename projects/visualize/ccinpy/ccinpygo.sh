#!/bin/bash
pip install pybind11
#
c++ -O3 -Wall -shared -std=c++11 -fPIC \
    $(python3 -m pybind11 --includes) \
    -o matrix_ops$(python3-config --extension-suffix) matrix_ops.cpp
