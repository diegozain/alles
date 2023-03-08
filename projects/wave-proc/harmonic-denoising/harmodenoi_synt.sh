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
ifort -qmkl -qopenmp -c -traceback -heap-arrays ../../../src/fortran/calculus.f90 ../../../src/fortran/readfiles.f90 ../../../src/fortran/harmodenoiser.f90 harmodenoi_synt.f90
#
ifort -qmkl -qopenmp harmodenoi_synt.o calculus.o readfiles.o harmodenoiser.o
mv a.out harmodenoi_synt.out
#
ifort -qmkl -qopenmp -c -traceback -heap-arrays harmodenoi_synt_.f90
#
ifort -qmkl -qopenmp harmodenoi_synt_.o calculus.o readfiles.o harmodenoiser.o
mv a.out harmodenoi_synt_.out
#
rm -f *.o *.mod
# ----------------------------------------------------------------------------
printf "       •••••••• author: diego domenzain  ••••••••••\n\n"
# ----------------------------------------------------------------------------
