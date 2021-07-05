# 3d DC & IP data visualization
diego domenzain
June 2021 @ Aarhus University 

## 👓 Sensitivity
For an **abmn** quadrupole, the sensitivity is:

👓 = ∇φ ⋅ ∇ v 

where,

-∇ ⋅ ∇ φ = source a - sink b

-∇ ⋅ ∇ v = source m - sink n

See Line's [.pdf](https://github.com/LineMeldgaardMadsen/ERT-XBH-sensitivity/blob/main/Catalogue_2D%20ERT_sensitivity.pdf) for more info.

## 👀 Visualize cross-borehole data 

### Localizing the sensitivities 👓 ⤍ ⚫

The location is the [geometric median](https://en.wikipedia.org/wiki/Geometric_median) of the **abmn** quadrupole positions in 3d space.

[![](../pics/pseudo-sensitivity.png)](./)

### Data structures 📚

```electrodes```
is a matrix of size ```nelectrodes```× 3. With entries of *xyz* coordinates, and indexes pointing to all electrodes.

```abmn```
is a matrix of size ```nabmn```× 4. With entries of indexes of ```electrodes```, and indexes pointing to all ABMN quadrupoles.

```pseud```
is a matrix of size ```nabmn```× 3. With entries of *xyz* coordinates, and indexes pointing to ```abmn```.

```data```
is a matrix of size ```nabmn```× 1. With entries of data (V, Ωm, ...), and indexes pointing to ```abmn```.

```klusters_```
is a cell of size ```nklu```. With entries of indexes of ```abmn``` (and so of ```pseud``` as well), and indexes pointing to clusters in ```pseud```.

### Data cleanup 📚⤍📘📗📕

Someone is like *dude, electrode #3 was shit. You have to remove it from the dataset*.

And then you're like *no problem,  it's cool*.

Let ```ielectrode_[a,b,m,n]``` be the electrodes indexes where ```ielectrode``` appears in ```abmn``` as a, b, m, or n.

```
ielectrode_a = abmn(find(abmn(:,1)==ielectrode_a),:);
ielectrode_b = abmn(find(abmn(:,2)==ielectrode_b),:);
ielectrode_m = abmn(find(abmn(:,3)==ielectrode_m),:);
ielectrode_n = abmn(find(abmn(:,4)==ielectrode_n),:);
```
Now you do this,

* delete ```ielectrode``` from ```electrodes```, 
* delete the collection ```ielectrode_[a,b,m,n]``` from: 
  * ```abmn```
  * ```pseud```
  * ```data```
* recompute ```klusters_```.

### AB.MN quadrupoles 🍀

[![](../pics/example-sensitivities.png)](./)

[![](../pics/pseudo-14electrodes.png)](./)

### Code 📝

* ```dc_pseudo_vis.m``` simple example script. Assumes all electrodes in Tx borehole can be **ab** pairs.
* ```dc_pseudo_vis_.m``` fancier example script. It can handle case when only *some* electrodes of Tx borehole can be **ab** pairs.
* ```dc_kaergaard.m``` simple example to visualize with field data.

## 🎯 Analytical forward model vs AarhusInv

The method for computing the analytical solutions is the method of images.

### Code 📝

* ```dc_analy_data.m``` homogeneous and two layered solution in 3d (only a 2d slice is computed). The entire 2d xz-slice is computed. *This script is a proof of concept*.
* ```dc_analy_data.m``` two layered solution in 3d (only a 2d slice is computed). Just the data points are computed. *This script is a proof of concept*.
* ```dc_analy_rip.m``` two layered solution in 3d (only a 2d slice is computed). 
  * Takes as input a modified version of the ```.rip``` (or ```.rap```) file, 
  * then it compares the ```.rip``` (or ```.rap```) output to the analytical solution for **one** quadrupole.
* ```dc_analy_rip_.m``` two layered solution in 3d (only a 2d slice is computed). 
  * Takes as input a modified version of the ```.rip``` (or ```.rap```) file, 
  * then it compares the ```.rip``` (or ```.rap```) output to the analytical solution for **all** quadrupoles.
* ```abmn2aainv.m```
* ```dc_aainv_vis.m```

## Visualize IP data
