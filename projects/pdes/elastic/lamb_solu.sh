#!/bin/bash

# # no openmp
# gfortran -c ../../../src/fortran/fourier.f90 lamb_solu.f90
# gfortran fourier.o lamb_solu.o

# with openmp
gfortran -c ../../../src/fortran/fourier.f90 -fopenmp lamb_solu.f90
gfortran -fopenmp fourier.o lamb_solu.o

# execute fortran code
./a.out

# plot
python visualize_f.py

# clean shit
rm *.o *.mod *.out *.dat fort.*



