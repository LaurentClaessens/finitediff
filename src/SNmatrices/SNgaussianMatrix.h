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
class SNgaussianMatrix
{

    private:
        std::array<T,tp_size> data;     // the first 'c' are unused and remain uninitialized
        unsigned int column;
    public :
        SNgaussianMatrix(const SNmatrix<T,tp_size>& , const unsigned int&);

        T get(unsigned int) const;
        
}

// CONSTRUCTOR  ---------------------------------------

template <class T,unsigned int tp_size>
SNgaussianMatrix::SNgaussianMatrix(const SNmatrix<T,tp_size>& A , const unsigned int& c):
    column(c)
{
    for (unsigned int i=c+1;i<tp_size;i++)
    {
        data.at(i)=-A.get(i,c)/A.get(c,c);
    }
}

// GETTER METHODS ---------------------------------------

template <class T,unsigned int tp_size>
T SNvector<>T,tp_size>::get(unsigned int i,unsigned int j) const
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
