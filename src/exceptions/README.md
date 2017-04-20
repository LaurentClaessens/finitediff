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
