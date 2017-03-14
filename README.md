# Finite Differences

I'm trying to implement the finite difference method for solving differential equations.

The mathematics is developed [here](http://laurent.claessens-donadello.eu/pdf/lefrido.pdf) as I'm learning it.

## deploy.sh

The script `deploy.sh` is the only working stuff up to now. It does the following things :
* show `git status`
* make clean
* make all 
* launch the unit tests

## Purpose

The official aim of this program is to implement the [finite difference method ](https://en.wikipedia.org/wiki/Finite_difference_method)  for solving differential equations.

Non-official purposes include
* A first experiment of implementing numerical methods.
* I (believe to) know that the finite difference method requires solving system of equations with sparse matrices. However I'm implementing the PLU decomposition for solving the system while I (believe to) know that there are more efficient ways to solve for sparse matrices.
* Learning the [cppunit framework](https://en.wikipedia.org/wiki/Cppunit) for my unit tests.
* Programming in C++ with no raw pointers.

## The classes

There are many different matrices types because there are many special cases : permutation, gaussian transformation, upper/lower triangular.

On the long run, I'm going to create a class for each and give them all the matrix semantics. In particular, define the multiplication operator for each (ordered) pair of matrix types.

The code for these matrices are in [./src/SNmatrices](src/SNmatrices/README.md).

