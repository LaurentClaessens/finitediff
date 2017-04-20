# Documentation's matin page


## The matrices

There are many different matrices types because there are many special cases : permutation, gaussian transformation, upper/lower triangular, etc.

### Hierarchy

The base type is `SNgeneric` (pure virtual) and the most general matrix type is `SNmatrix`.


### Template

The matrices have two template arguments :

* `class T` is supposed to be a numeric class like `double`. It is the type of the matrix elements. 
* `unsigned int tp_size` is the size of the matrix (well, the size is the square of that you know).

### Operators

* Operators like equality and multiplication are defined for each (ordered) pair of these special cases.
* Many of these operators are returning the answer by value and are relying on the RVO.
* They are defined in `SNoperators.h`

### More ...

Read more [here](src/SNmatrices/README.md).

## Exceptions

We define some exceptions adapted to our matrices manipulations. Read more  [here](src/exceptions/README.md).
