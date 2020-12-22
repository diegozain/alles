# Bayesian Optimization
diego domenzain
September 2020 @ Colorado School of Mines

Suppose you have an objective function with many hyper-parameters. 

__how do you find the best hyper-parameters?__

You _sample_ the objective function many times with different hyper-parameters.
You then use Gaussian Kernels to _grow_ in value an approximate of the objective function at these hyper-parameter locations.

In yellow are the initial samples, in blue the true minimum, in green the solution path, in red the recovered minimum. 

_True_ and _Approximate_ are the true and recovered objective functions respectively.

[![](../pics/bayes-opti-ex.png)](./)
