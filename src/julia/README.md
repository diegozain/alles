# My code in Julia

This folder is a _Julia_ clone of the _Matlab_ code.

It is a work in progress. The idea is to fill it in _sloooowly_.

The next functions on the list are:

 * ```DDx_DDz.m```
 * ```geom_median.m```
 * ```shape_curve.m```

## Generic stuff

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

## A specific example

Let's run the functions ```differentiate_cube.jl``` and ```integrate_cube.jl```.

They both act on a cube matrix, and they both do so by reference.

*Matlab* (at least my version) __cannot pass by reference__. I checked with these same functions and looked at the task monitor. At least an extra copy of the cube was stored.

*Julia* __can pass by reference__, but is about __4 times slower__ than *Matlab*.

```julia
# change path to wherever alles is
cd("alles/src/julia/")

# this includes the functions we want to use.
# every time you change the function you need to include it again.
include("differentiate_cube.jl")
include("integrate_cube.jl")

# numbers without the .0 are ints automatically
nx=100;
nz=100;

# make nt so that the cube is of size 1Gb
nt=round(1/(nx*nz*8*1e-9));
nt=convert(Int,nt);

# some fake data
dt=0.1;
v=rand(nz,nx,nt);

# some idiot decided it was a good idea to
# remove the tic toc functions in julia.
t1=time_ns();

# the exclamation point ! makes the function to pass by reference, 
# so no extra memory is used.
differentiate_cube!(v,dt);

t2=time_ns();
elapsed_time=(t2-t1)/1.0e9

# now the integration
t1=time_ns();
integrate_cube!(v,dt);
t2=time_ns();
elapsed_time=(t2-t1)/1.0e9
```