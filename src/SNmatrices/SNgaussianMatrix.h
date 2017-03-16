/*
Copyright 2017 Laurent Claessens
contact : laurent@claessens-donadello.eu

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __SNgaussianMatrix_H__174236
#define __SNgaussianMatrix_H__174236

#include <array>

#include "SNgeneric.h"
#include "../SNexceptions.cpp"


#include "../DebugPrint.h"
DebugPrint debug_printGS;

/*
   This represents a gaussian matrix.

   A gaussian matrix is the multiplicative matrix which serves to eliminate
   a column under the diagonal.

   Let A be a matrix and 'c' be an index. Let 'm=A_{c,c}'  The gaussian matrix for 'A' on the index 'c' is full of zeros apart of
   - 1 on the diagonal
   - A_{i,c}/m on position (i,c) when i>c. 

   That is there is something non trivial only on the column 'c', under the diagonal.
 *
 *
 * TODO : giving 'c' as a template argument, one can in average spare half of the memory by allowing 'c' as size in 'data' instead of 'tp_size'.
*/

// forward definition

template <class T,unsigned int tp_size>
class SNgeneric;
template <class T>
class SpecialValue;


// THE CLASS HEADER -----------------------------------------

template <class T,unsigned int tp_size>
class SNgaussianMatrix : public SNgeneric<T,tp_size>
{

    private:
        std::array<T,tp_size> data;     // see implementation of "_at"
        SpecialValue<T> checkForSpecialElements(const unsigned int& i,const unsigned int& j) const;

        // populate the matrix from the elements of the given matrix
        template <class U,unsigned int s>
        void populate_from(const SNgeneric<U,s>&);
        T _get(unsigned int,unsigned int) const;
        T& _at(unsigned int,unsigned int);
    public :
        const unsigned int column;

        template <class U,unsigned int s>
        SNgaussianMatrix(const SNgeneric<U,s>& , const unsigned int&);

        
};

// CONSTRUCTOR  ---------------------------------------

template <class T,unsigned int tp_size> 
template <class U, unsigned int s>
void SNgaussianMatrix<T,tp_size>::populate_from(const SNgeneric<U,s>& A)
{
    if (s!=tp_size)
    {
        throw IncompatibleMatrixSizeException(tp_size,s);
    }

    T m = A.get(column,column);
    for (unsigned int i=column+1;i<tp_size;i++)
    {
        this->at(i,column)=-A.get(i,column)/m;
    }
}

template <class T,unsigned int tp_size> 
template<class U,unsigned int s>
SNgaussianMatrix<T,tp_size>::SNgaussianMatrix(const SNgeneric<U,s>& A , const unsigned int& c):
    column(c)
{
    populate_from(A);
}

// UTILITIES  ---------------------------------------

template <class T,unsigned int tp_size>
SpecialValue<T> SNgaussianMatrix<T,tp_size>::checkForSpecialElements(const unsigned int& i,const unsigned int& j) const
{
    if (i==j)
    {
        return SpecialValue<T>(1,true);
    }
    if (j!=column)
    {
        return SpecialValue<T>(0,true);
    }
    if (i<j)
    {
        return SpecialValue<T>(0,true);
    }
    return SpecialValue<T>(0,false);
}

// _AT AND _GET METHODS ---------------------------------------


// checkForSpecialElements(i,j)
// checks for element (i,j). If this is a special element 
// (a one whose value is fixed by the fact that we are a gaussian matrix)
// then it returns a 'SpecialValue<T>' with its boolean part set to true,
// meaning that this is a special value.
// In that case :
// - _get returns the corresponding value.
// - _at throws SNchangeNotAllowedException.
//
// If this is not a special value, it returns a 'SpecialValue<T>' with boolean part
// set to false. 
// In that case :
// - both _get and _at have to search in the stored values.


template <class T,unsigned int tp_size>
T SNgaussianMatrix<T,tp_size>::_get(unsigned int i,unsigned int j) const
{
    SpecialValue<T> sv=checkForSpecialElements(i,j);
    if (sv.special)
    {
        return sv.value;
    }
    return data.at(i-column-1);
}

template <class T,unsigned int tp_size>
T& SNgaussianMatrix<T,tp_size>::_at(unsigned int i,unsigned int j) 

    // The elements are stored in 
    //   std::array<T,tp_size> data
    // while the matrix only contains non fixed values on one column, under
    // the diagonal.
    //
    // First element of 'data' stores the first non fixed value.
    // Example for a gaussian matrix on column 1 :
    //
    //  1  0   0  0
    //  0  1   0  0
    //  0  d0  1  0
    //  0  d1  0  1
    //

{
    SpecialValue<T> sv=checkForSpecialElements(i,j);
    if (sv.special)
    {
        throw SNchangeNotAllowedException(i,j);
    }
    return data.at(i-column-1);
}

#endif
