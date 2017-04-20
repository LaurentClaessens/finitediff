# Documentation's matin page


## The matrices

There are many different matrices types because there are many special cases : permutation, gaussian transformation, upper/lower triangular, etc.

### Hierarchy

The base type is `SNgeneric` (pure virtual) and the most general matrix type is `SNmatrix`.

The code for these matrices are in [./src/SNmatrices](src/SNmatrices/README.md).

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
SNlowerTriangular<int,4> A;
A.at(1,2);     // raises
```
A lower diagonal matrix will always have element `(1,2)` set to zero.

### `IncompatibleMatrixSizeException`

Raised when trying to perform operation with matrices with incompatible sizes.

```C++
SNlowerTriangular<int,4> A;
SNmatrix<double,2> B;
bool a=(A==B);        // raise
```

```C++
SNmatrix<double,2> A;
SNmatrix<double,3> B;
auto C=A*B;     //raise
```
