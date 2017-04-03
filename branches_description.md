This file is mainly for personal usage. It describes the purpose of each created branch.


# remove_snmatrix_reference

The class `SNline` has a reference to the matrix it is issued from. We remove that reference.


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
