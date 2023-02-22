# 🟪🔴 = 🔷🔺



Solves ```Ax = b```  by QR magic using the *lapack* in **oneapi-mkl**.

```
Ax  = b
A   = QR
QRx = b
Rx  = Q'b
```
In general, it solves this scenario:

```
     _________     _____________     _____________
    |         |   |             |   |             |
 m  |   A     |   |     X       |   |      B      |
    |         | * |             | = |             | m
    |         |   |_____________|   |             |
    |_________|         p           |_____________|
         n                                 p
```

but it is written for ```p=1```.

## ⚠ 

When ```m≠n```  it is not working.  

Not really a problem, just do ```A⟵A'A``` & ```b⟵A'b``` using ```dtrsm``` & ```dgemm```.

## 🗃

* ```linreg.f90``` minimal example and well explained.
* ```linreg_.f90``` huge dense matrix 🤪.
* ```parafit.f90``` fit a parabola to a set of points.
