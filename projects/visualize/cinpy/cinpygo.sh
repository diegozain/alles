#!/bin/bash
gcc -shared -o matrix_ops.so -fPIC matrix_ops.c
python setup.py build_ext --inplace
