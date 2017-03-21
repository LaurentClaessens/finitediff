# Finite Differences

I'm trying to implement the finite difference method for solving differential equations.

The mathematics is developed [here](http://laurent.claessens-donadello.eu/pdf/lefrido.pdf) as I'm learning it.

## Purpose

The official aim of this program is to implement the [finite difference method ](https://en.wikipedia.org/wiki/Finite_difference_method)  for solving differential equations.

Non-official purposes include
* A first experiment of implementing numerical methods.
* I (believe to) know that the finite difference method requires solving system of equations with sparse matrices. However I'm implementing the PLU decomposition for solving the system while I (believe to) know that there are more efficient ways to solve for sparse matrices.
* Learning the [cppunit framework](https://en.wikipedia.org/wiki/Cppunit) for my unit tests.
* Programming in C++ with no raw pointers.

## Latest news

* I'm currently implementing the PLU decomposition.
* U is done. P is almost done (I have the list of line swap, but not yet a matrix). L is still to be done.

## deploy.sh

The script `deploy.sh` does the following things :
* show `git status`
* make clean
* make all 
* launch the unit tests

## The matrices

There are many different matrices types because there are many special cases : permutation, gaussian transformation, upper/lower triangular, etc.

### Hierarchy

The base type is `SNgeneric` (pure virtual) and the most general matrix type is `SNmatrix`.

The code for these matrices are in `src/SNmatrices/README.md`.

### Template

The matrices have two template arguments :

* `class T` is supposed to be a numeric class like `double`. It is the type of the matrix elements. 
* `unsigned int tp_size` is the size of the matrix (well, the size is the square of that you know).

### Operators

* Operators like equality and multiplication are defined for each (ordered) pair of these special cases.
* Many of these operators are returning the answer by value and are relying on the RVO.


## Exceptions

### `SNoutOfRangeException`

Raised when one ask by `get` or `at` an element that is out of the range of the matrix. 

```C++
SNmatrix<int,4> A;
A.get(5,1);     // raises
```

### `SNchangeNotAllowed`

Raised when one uses the method `at` for requesting an element that cannot be changed (use `get` instead).

```C++
SNlowerTriangularMatrix<int,4> A;
A.at(1,2);     // raises
```
A lower diagonal matrix will always have element `(1,2)` set to zero.

### `IncompatibleMatrixSizeException`

Raised when trying to perform operation with matrices with incompatible sizes.

```C++
SNlowerTriangularMatrix<int,4> A;
SNmatrix<double,2> B;
bool a=(A==B);        // raise
```

```C++
SNmatrix<double,2> A;
SNmatrix<double,3> B;
auto C=A*B;     //raise
```
