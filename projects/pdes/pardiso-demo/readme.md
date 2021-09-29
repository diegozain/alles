# *Pardiso* via *oneAPI* 🏝

diego domenzain

summer 2021 @ Aarhus University

```
$> cmd.exe "/K" '"C:\Program Files (x86)\Intel\oneAPI\setvars.bat" && powershell'
```
---

```
pardiso-demo
      |--- *you should be here*

$> ifort /Qmkl [name of file].f
$> .\[name of file].exe
```
For the magic function ```pardiso( args )```, the arguments are explained in page **2598** & **2617** of ```onemkl-developerreference-fortran.pdf```.

For the crazy 😵 ```iparam``` array see page **2622**.

## Simple solver ⬜🔷 = 🔶

```fortran90
! solves X for the problem:
!
!     _________     _____________     _____________
!    |         |   |             |   |             |
! n  |   A     |   |     X       |   |      B      | n
!    |         | * |             | = |             |
!    |_________|   |_____________|   |_____________|
!        n              m                  m
```

Specifically,
```fortran90
A = 
     2     3     0     0     0
     3     0     4     0     6
     0    -1    -3     2     0
     0     0     1     0     0
     0     4     2     0     1
     
b =

     8
    45
    -3
     3
    19
! ----------------------------------------------------------
nnz = 12
n   = 5
m   = 1
! ----------------------------------------------------------
a_ = a_(nnz)
ja = ja(nnz)
ia = ia(nnz)
! ia_ is the CRS array
ia_= ia_(n+1)

b = b(n,m)
! ----------------------------------------------------------
ia = [ 1,1,2,2,2,3,3,3,4,5,5,5 ]
ij = [ 1,2,1,3,5,2,3,4,3,2,3,5 ]
a_ = [ 2,3,3,4,6,-1,-3,2,1,4,2,1 ]
! ia_ is the CRS array
ia_= [ 1,3,6,9,10,13 ]
! ----------------------------------------------------------
```

After running, the solution should be

```fortran90
x =

    1.0000
    2.0000
    3.0000
    4.0000
    5.0000
```
---
