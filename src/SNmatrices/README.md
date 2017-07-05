#### Creating a matrix

Internally, `SNmatrix' uses a `std::array` to store the matrix entries, but you should not try to take advantage of that.

The correct way to populate your matrix is :


#### const ?

Many methods are not 'const'. Here is the rationale.

- When you extract an element from the matrix with the method 'getElement', 
  the element contains a reference (in fact, a pointer) to the original matrix.

- The matrix itself is supposed to change with time. Many operations are "in-place" replacement. Examples :
   * during the 'PLU' decomposition, the matrix itself becomes the 'U' matrix and th original is lost.
   * swaping lines is in place


 Note : the first point is maybe useless and some 'const' could be added in the future.

 Note : you cannot instantiate  
```c++
SNmatrix<const double,7> sn;
```

If you want to know the value of a matrix element, you can use two methods :
```c++
T& at(unsigned int,unsigned int);
T get(unsigned int,unsigned int) const;
```
The first one returns a reference ans should only be used to populate the matrix
(or more generally to modify the matrix) and the second one (get) should be
used whenever you just want to know the value.


#### Why 'subtract' instead of 'operator-'  ?

The operation
```c++
A.subtract(B)
```
changes A into A-B.

The operation (not yet implemented)
```c++
C=A-B
```
creates a new matrix that has the value A-B.


#### "PLU" decomposition

P : permutation matrix,
L : lower triangular matrix with `1` on the diagonal,
U : upper triangular matrix.

The mathematics can be found (in French) here :
http://laurent.claessens-donadello.eu/pdf/mazhe.pdf
