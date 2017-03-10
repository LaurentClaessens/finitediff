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

### `SNmatrix`

This is a template class that represent a fixed size square numerical matrix (SN=Square Numerical).

```C++
template <class T,unsigned int tp_size>
class SNmatrix
```

* `T` is the type of values contained in the matrix. This is expected to be a numeric type in the sense that it has to accept comparison, absolute value (from cmath) and other operations like that.
* `tp_size` is the size of the matrix.


### `SNelement`

This is a template class that represent an element of a `SNmatrix`.

```C++
template <class T,unsigned int tp_size>
class SNelement
```

The template arguments `T` and `tp_size` have the same meaning. 

A `SNelement` contains a reference to its matrix. 

