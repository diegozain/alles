#!/bin/bash
gfortran -fPIC -shared -o libfortran.so fortran_code.f90
rm *.mod
g++ -fPIC -shared -o libcc.so cc_code.cc -ldl
g++ -o test_main test_main.cc -L. -lcc -Wl,-rpath,.


