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

// THE CLASS HEADER -----------------------------------------

template <class T,unsigned int tp_size>
class SNgaussianMatrix : public SNgeneric<T,tp_size>
{

    private:
        std::array<T,tp_size> data;     // the first 'c' are unused and remain uninitialized

        // populate the matrix from the elements of the given matrix
        template <class U,unsigned int s>
        void populate_from(const SNgeneric<U,s>&);
    public :
        const unsigned int column;

        template <class U,unsigned int s>
        SNgaussianMatrix(const SNgeneric<U,s>& , const unsigned int&);

        unsigned int getSize() const;
        T get(unsigned int,unsigned int) const;
        
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
    for (unsigned int i=column+1;i<tp_size;i++)
    {
        data.at(i)=-A.get(i,column)/A.get(column,column);
    }
}

template <class T,unsigned int tp_size> 
template<class U,unsigned int s>
SNgaussianMatrix<T,tp_size>::SNgaussianMatrix(const SNgeneric<U,s>& A , const unsigned int& c):
    column(c)
{
    populate_from(A);
}

// GETTER METHODS ---------------------------------------

template <class T,unsigned int tp_size>
unsigned int SNgaussianMatrix<T,tp_size>::getSize() const
{
    return tp_size;
}

template <class T,unsigned int tp_size>
T SNgaussianMatrix<T,tp_size>::get(unsigned int i,unsigned int j) const
{
    if (i>tp_size or j>tp_size)
    {
        throw SNoutOfRangeException(i,j,tp_size);
    }
    if (i==j)
    {
        return 1;
    }
    if (j!=column)
    {
        return 0;
    }
    if (i<j)
    {
        return 0;
    }
    return data.at(j);   
}

#endif
