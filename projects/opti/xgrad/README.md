# Cross-gradients
diego domenzain
September 2020 @ Colorado School of Mines

## What is the cross-gradient constraint?

Two (or more) physical parameters can share structure. Here, the parameters are two-dimensional matrices.

The question is then, how can we transfer the structure of one parameter onto the other?

__This script is an example of cross-gradient inversion.__

## The parameters

Here, they are called _a_ and _b_. They do not represent any physical parameters. However, they do emulate the mathematical description of a physical parameter.

## The inversions

There are three possibilities,

* _a_ and _b_ get to be like each other,
* _a_ gets to be like _b_ (and _b_ is fixed),
* _b_ gets to be like _a_ (and _a_ is fixed).

These inversions are done with gradient descent. 

## The scripts

* ```xgrad_ex_circs.m``` perform these inversions on two discs for _a_ and _b_.
* ```xgrad_ex_circs.m``` performs these inversions on a disc and a square for _a_ and _b_.