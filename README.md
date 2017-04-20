# Finite Differences

I'm trying to implement the finite difference method for solving differential equations.

* You can read the [documentation](http://laurent.claessens-donadello.eu/finitediff/html/index.html) by Doxygen (work in progress).
* The mathematics is developed [here](http://laurent.claessens-donadello.eu/pdf/lefrido.pdf) as I'm learning it.

## Purpose

The official aim of this program is to implement the [finite difference method ](https://en.wikipedia.org/wiki/Finite_difference_method)  for solving differential equations.

Non-official purposes include
* A first experiment of implementing numerical methods.
* I (believe to) know that the finite difference method requires solving system of equations with sparse matrices. However I'm implementing the PLU decomposition for solving the system while I (believe to) know that there are more efficient ways to solve for sparse matrices.
* Learning the [cppunit framework](https://en.wikipedia.org/wiki/Cppunit) for my unit tests.
* Programming in C++ with no raw pointers.

## Latest news

* The PLU decomposition is done.
* Next step : solving linear systems. See the branch `alpha`.

## deploy.sh

The script `deploy.sh` does the following things :
* show `git status`
* make clean
* make all 
* launch the unit tests

The test are covering more or less everything, but the most interesting ones are in `plu_unit_tests.cpp`. There are many matrices for which the PLU decomposition are tested against the one provided by [sage](https://sagemath.org).
