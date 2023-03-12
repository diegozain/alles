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
#
ifort -qmkl -qopenmp -c -traceback -heap-arrays ../../../src/fortran/calculus.f90 ../../../src/fortran/readfiles.f90 ../../../src/fortran/harmodenoiser.f90 harmodenoibin.f90
#
ifort -qmkl -qopenmp harmodenoibin.o calculus.o readfiles.o harmodenoiser.o
mv a.out harmodenoibin.out
#
rm -f *.o *.mod
# ----------------------------------------------------------------------------
printf "       •••••••• author: diego domenzain  ••••••••••\n\n"
# ----------------------------------------------------------------------------
