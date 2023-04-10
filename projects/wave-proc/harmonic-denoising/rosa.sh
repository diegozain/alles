#!/bin/bash
# source /opt/intel/oneapi/setvars.sh intel64
printf "\n       •••••••••••••••  compiling  ••••••••••••••••\n"
rm -f *.o *.mod
#
ifort rosa.f90
#
mv a.out rosa.out
#
rm -f *.o *.mod

