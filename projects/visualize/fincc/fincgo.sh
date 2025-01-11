#!/bin/bash
gfortran -fPIC -shared -o libfortran.so fortran_code.f90
rm *.mod
gcc -fPIC -shared -o libc.so c_code.c -ldl
gcc -o test_main test_main.c -L. -lc -Wl,-rpath,. 


