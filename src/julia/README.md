# My code in Julia

This folder is a _Julia_ clone of the _Matlab_ code.

It is a work in progress. The idea is to fill it in _sloooowly_.

The next functions on the list are:

 * ```DDx_DDz.m```
 * ```geom_median.m```
 * ```shape_curve.m```

## How to use these functions

_Julia_ needs to import libraries and stuff.

```julia
pwd()
cd("path/to/files/")

using SparseArrays, LinearAlgebra, FFTW, MAT
using PyPlot
using Plots

include("file.jl")

# like linspace:
t=LinRange(0, 2*pi, 1000);
y=sin.(t);
dt=t[2]-t[1];

# add dimension so it is (nt x 1)
y=reshape(y, (size(y)...,1));

# to plot:
using Plots
plot(x,y)
plot!(x,y_)

# full matrix:
Matrix(a)

# to imagesc:
using PyPlot
# pcolormesh(M)
imshow(M)
colorbar()
```
