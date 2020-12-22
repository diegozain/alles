# Bayesian Optimization
diego domenzain
September 2020 @ Colorado School of Mines

Suppose you have an objective function with many hyper-parameters. 

__How do you find the best hyper-parameters?__

You _sample_ the objective function many times with different hyper-parameters.
You then use Gaussian Kernels to _grow_ in value an approximate of the objective function at these hyper-parameter locations.

---

The initial samples are in yellow, true minimum in blue the , the solution path in green, the recovered minimum in red. 

_True_ and _Approximate_ are the true and recovered objective functions respectively.

[![](../pics/bayes-opti-ex.png)](./)
