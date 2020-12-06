# Distanced but connected
diego domenzain
September 2020 @ Colorado School of Mines

In the times of COVID19 we are told to stay 2m apart from other people.

If a large group of people were to meet, 

__how could people be as close together while still distancing themselves?__

A similar problem called _circle packing_ has been solved in the literature.
This problem aims to enclose _n_ circles of the same radii in the smallest possible circle.

[Circle packing in Wikipedia](https://en.wikipedia.org/wiki/Circle_packing_in_a_circle)
[Awesome and long list of known circle packings](http://hydra.nat.uni-magdeburg.de/packing/cci/)
Graham, Ronald L., et al. _Dense packings of congruent circles in a circle._ Discrete Mathematics 181.1-3 (1998): 139-154.

The approach for solving the _circle packing_ problem is to 

1. perform a non-linear optimization on a potential-minimizer type objective function,
1. fine tune this solution by tightening the contacts between the _n_ circles. 
   This is in turn another non-linear optimization.

The code presented here is a bit different:

1. we do not assume a bounded domain (e.g. a circle) for our solution to be in,
1. we do not tighten the contacts after step 1 above.

Furthermore, the approach presented here uses three objective functions rather than just one.
This approach enables us to exit local minima and achieve a qualitative much better solution compared to just one objective function.

The role of these objective functions is:

1. Contract all nodes (e.g. people) to their geometric median
1. Expand all nodes as if they were repelling electric charges
1. Minimize their distance to 2m

__The contraction-expansion-minimization strategy helps the nodes to find their most compact form.__

---

Below is the output of 19 people gathering in an (approximate) optimal sense while keeping 2m distance.

[![](../pics/covid_19-people.png)](./)
