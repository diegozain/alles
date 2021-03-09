# OpenMP
diego domenzain
March 2021 @ Colorado School of Mines

## Implementation of OpenMP 

Communicating with multiple processors (MP) using the openMP framework.

__These functions are examples on how to use openMP.__

They are taken from the book [The openMP common core](https://mitpress.mit.edu/books/openmp-common-core).

The book itself comes with code in *C* and *Fortran*. The code can be found [here](http://ompcore.com/) and [here](https://github.com/tgmattso/OmpCommonCore/tree/master/Book/).

This repo is just a way for me to learn this stuff. The code is to be compiled and run here in a simple way.

For *C*,
```shell
gcc -fopenmp file.c
./a.out
```

For *Fortran*,
```shell
gfortran -fopenmp file.f90
./a.out
```

My old Mac doesn't have openMP enabled for *C*. *Homebrew* won't even download ```gcc```. 

*Fortran* seems to be doing ok though.

---

