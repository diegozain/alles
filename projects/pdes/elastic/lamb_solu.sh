#!/bin/bash

gfortran -c ../../../src/fortran/fourier.f90 lamb_solu.f90
gfortran fourier.o lamb_solu.o
./a.out
python visualize_f.py

rm *.o *.mod *.out *.dat fort.*
