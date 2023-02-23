#!/bin/bash
# ----------------------------------------------------------------------------
#
#                ifort in 🚀
#
#   diego 2022
# ----------------------------------------------------------------------------
# 🚀
# source /opt/intel/oneapi/setvars.sh intel64
# ----------------------------------------------------------------------------
printf "\n       •••••••••••••••  compiling  ••••••••••••••••\n"
# ----------------------------------------------------------------------------
rm -f *.o *.mod
# ----------------------------------------------------------------------------
ifort -qmkl -qopenmp -c -traceback -heap-arrays /opt/intel/oneapi/mkl/2023.0.0/include/lapack.f90 linreg.f90 linreg_.f90 parafit.f90 linefit.f90
#
ifort -qmkl linreg.o lapack.o
mv a.out linreg.out
ifort -qmkl -qopenmp linreg_.o lapack.o
mv a.out linreg_.out
ifort -qmkl parafit.o lapack.o
mv a.out parafit.out
ifort -qmkl linefit.o lapack.o
mv a.out linefit.out
# ----------------------------------------------------------------------------
ifort -qmkl -qopenmp -c -traceback -heap-arrays ../../../src/fortran/qrfits.f90 qrfits_ie.f90
#
ifort -qmkl -qopenmp qrfits_ie.o lapack.o qrfits.o
mv a.out qrfits_ie.out
# ----------------------------------------------------------------------------
rm -f *.o *.mod
# ----------------------------------------------------------------------------
printf "       •••••••• author: diego domenzain  ••••••••••\n\n"
# ----------------------------------------------------------------------------
