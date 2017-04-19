This file is mainly for personal usage. It describes the purpose of each created branch.


# remove_snmatrix_reference

The class `SNline` has a reference to the matrix it is issued from. We remove that reference.

# operator_member

On this branch we work on [The Decision between Member and Non-member](http://stackoverflow.com/questions/4421706/operator-overloading/4421729#4421729).

We will follow the three rules of thumb proposed there.

# multigaussian

We develop the class `SNmultiGaussian` and its tests. This class will be used for the PLU decomposition.

# no_reference_get


In `SNgeneric` the methods `at` and `get` are prototyped as

```c++
virtual T& at(const m_num,const m_num) final;        
virtual T get(const m_num&,const m_num&) const final;
```

For consistency we pass the arguments of `get` by value.

# m_num_increment

Rationalize how `i++` and `++i` work when `i` is of type `m_num`.

# getColumn

For `SNgaussian`, we set the attribute `column` as private and add the 'get' and 'set' methods.

# utilities

* Regroup the small functions here and there in `Utilities.h` and `Utilities.cpp`
* Develop and test the function `value_length`

# prod_perm

Work on the product of permutations.

* define Mpermutation * MelementaryPermutation
* Make non-member the other one.
