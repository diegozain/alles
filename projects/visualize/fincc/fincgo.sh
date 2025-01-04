#!/bin/bash
gfortran -fPIC -shared -o libfortran.so fortran_code.f90
gcc -fPIC -shared -o libc.so c_code.c -ldl
gcc -o test_main test_main.c -L. -lc -Wl,-rpath,. -lc
# /usr/bin/ld: /usr/lib/gcc/x86_64-linux-gnu/11/../../../x86_64-linux-gnu/Scrt1.o: undefined reference to symbol '__libc_start_main@@GLIBC_2.34'
# /usr/bin/ld: /lib/x86_64-linux-gnu/libc.so.6: error adding symbols: DSO missing from command line
# collect2: error: ld returned 1 exit status


