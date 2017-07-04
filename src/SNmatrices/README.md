
### `SNmatrix`

#### Gauss's PLU decomposition

An invertible matrix `A` can be written as `A=PLU`  where
* `P` is a permutation matrix,
* `L` is a lower triangular matrix with unital diagonal,
* `U` is an upper triangular matrix.

When `A` is a `SNmatrix`, we can use the method
```C++
template <class T,unsigned int tp_size>
void SNmatrix<T,tp_size>::makeUpperTriangular()
```
to turn the matrix `A` into its upper triangular matrix. 

I'm currently working on the other two. The interface will be the following (I'm writing here to clarify my own ideas).
* Let `A` be a `SNmatrix`
* The call to `A.PLU()` returns an object with `P`, `L` and `U` fields.
* The `U` field is a reference to the original matrix.
* The original matrix is no more available because it became the `U` matrix.

The reference trick is to avoid the creation a creation and a copy of matrix. 

### `SNelement`

This is a template class that represent an element of a `SNmatrix`.

```C++
template <class T,unsigned int tp_size>
class SNelement
```

The template arguments `T` and `tp_size` have the same meaning. 


### Permutations

* `Mpermutation` has one template parameter : `tp_size` which is the size of the set on which it acts.
* `SNpermutation` is the associated matrix. It is obviously composed of many zero's and some ones. But these are `0` and `1` of the template type `T`.
 
### `SNlowerTriangular`

There is a subtlety with the `at` method.

* Usually, the `get` method returns *by value* the value of the requested entry. This is the usual behaviour and it is the behaviour here.
* The `at` method returns *by reference* the value of the requested entry. This is useful for populate or modify the matrix. However asking for the element `1,3` on a `SNlowerTriangular` will raise a `SNnonAllowedChange` exception.

## Operators

We define some operators between different types of matrices : mainly `==` and `*`. For optimization reasons, the implementation depend on the precise types.

Remark : many operators depend on the "return value optimization" because they are building a local object and then returning it by value. A (simplified) example with the matrix product :


```C++
template <class U,class V,unsigned int s,unsigned int t>
SNmatrix<U,s> operator*(const SNmatrix<U,s>& A, const SNmatrix<V,t>& B)
{
    SNmatrix<U,s> ans;
    // populate ans
    return ans      // do not copy by RVO 
}
```

### Product operator `*`

The declaration of the product `SNmatrix*SNmatrix` is as follow : 

```C++
template <class U,class V,unsigned int s,unsigned int t>
SNmatrix<U,s> operator*(const SNmatrix<U,s>& A, const SNmatrix<V,t>& B)
```

The matching of the size if checked but the return type will be the one of the left matrix.

Let the example :
```C++
SNmatrix<int,2> A; 
SNmatrix<double,2> B; 
// populate the matrices A and B
auto C=A*B   // inacurate !
```
Here the type of `C` will be `SNmatrix<int,2>` and the product will be inaccurate  because of the products `int*double` are casted to `int`.

By the way, the matrix product is not commutative, so do not try to workaround with
```C++
auto C=B*A    // instead of A*B
```


## Implementation details

### SNmatrix<class T,unsigned int tp_size>

#### Template arguments

- `class T` : the type of the matrix entries (typically 'double')
- `unsigned int tp_size` : the size of the matrix (the number of entries is the square of that).


#### Creating a matrix

Internally, `SNmatrix' uses a `std::array` to store the matrix entries, but you should not try to take advantage of that.

The correct way to populate your matrix is :

```C++
SNmatrix<int,2> sn;
sn.at(0,0)=1;
sn.at(0,1)=2;
sn.at(1,0)=3;
sn.at(1,1)=4;
```

Notice that the elements are numbered from `0` to `tp_size-1`. Not from `1`.


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
