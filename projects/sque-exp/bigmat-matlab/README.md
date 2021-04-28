# Big Mat
diego domenzain
April 2021 @ Colorado School of Mines

## Assign by reference in *Matlab*

Let's say you are working on a big matrix stored in memory.

There is no memory left, so you write the next action on the matrix assigning *by reference*.

You run the program and it turns out *Matlab* did not really assign by reference. Let's skip the details on why *Matlab* didn't do that.

This trick will let you implement *by reference* no matter what.

1. Write a class handle like ```packer.m```.
1. Write your *by reference* code in a function designed to act on an instance of the ```packer.m``` class, and with no output. For example, see ```integrate_line_.m```.
1. When you encounter the big matrix, pack it into a ```packer.m``` instance.
1. Clear the big matrix.
1. Pass *by reference* on the ```packer.m``` instance.

## ```bigmat.m``` 

Very simple example. Two variables are declared, ```t``` and ```y(t)```. Then, the antiderivative of ```y``` is computed numerically.

Without this trick, *Matlab* would make a copy of ```y``` inside the ```integrate_line.m``` function.

With this trick, the script only uses twice the size of ```t``` in memory.

## ```bigmat_.m``` 

What if the big matrix already lives inside a structure?

This example shows that case. Turns out it is as efficient as before.

At any given point, the script only uses twice the size of ```t``` in memory.

---