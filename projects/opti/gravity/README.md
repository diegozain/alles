# Gravity inversion
diego domenzain
September 2020 @ Colorado School of Mines

## What is gravity inversion?

Gravity is dependent on density. By measuring the gravity field, it is possible to solve for (invert) density.

__This script is an example of density inversion at depth, jointly using gravity gradient data on x and z.__

```matlab
the gravity field u = (ux,uy,uz) at point ro is given by,

u = G int_V( rho(r) * (r-ro) / ||r-ro||^3 )

G is the gravitational constant
where V is the integration volume
r is the integration variable

Assume the data u has been normalized by G.

At receiver locations the data is given by,

d_x = M*Lx*rho
d_z = M*Lz*rho

M is the measuring operator
Lx and Lz are matrices depending only on geometry
rho is the density in vector form (rho=rho(:))

The data sensitivity with respect to the density (jacobian) is,

Jx = M*Lx

and the gradient is

g_x = Jxt*e_x

where e_x is the residual ( e_x = d_x - d_x_observed ). t is for transposed.
```

## References
1. [__3-D inversion of gravity data__](https://library.seg.org/doi/abs/10.1190/1.1444302). Yaoguo Li, and Douglas W. Oldenburg. *Geophysics* 1998.

The inversion itself is different than that of Yaoguo's,
* We do not use "depth weighting" because it is not necessary.
* We explicitly write the gradients.

Depth resolution *is* an issue nonetheless. 

---

Below is an example of the true and recovered density. The colormap is the same for both pictures.

[![](../pics/gravity_inversion.png)](./)
