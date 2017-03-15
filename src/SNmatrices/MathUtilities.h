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


/*
 Define some mathematical intermediate functions.
*/

#ifndef __MATHUTILITIES_H__055253__
#define __MATHUTILITIES_H__055253__

#include "../SNexceptions.cpp"

template <class T,class U>
bool componentWiseeEquality(const T& A,const U& B)
    // test the equality componentwise between to matrix objects.
{
    unsigned s=A.getSize();
    if (B.getSize()!=s)
    {
        return false;       // maybe an exception would be more appropriate.
    }
    for (unsigned int l=0;l<s;l++)
    {
        for (unsigned int c=0;c<s;c++)
        {
            if (A.get(l,c)!=B.get(l,c))
            {
                return false;
            }
        }
    }
    return true;
}

template <class T,class U,unsigned int s,unsigned int t>
T matrixProductComponent(const SNgeneric<T,s>& A,const SNgeneric<U,t>& B,unsigned int i,unsigned int j)
    // compute the component (i,j) of the matrix product A*B
    // - A,B  : matrices
    // - i,j  : the element of A*B to be computed
    //
    // The size compatibility is not checked since this functions is supposed 
    // to be called from something like an operator * that has already checked
    // the compatibility.
    //
    // As stated in the README.md, the type of the elements in the returned matrix
    // is the one of the first matrix of the product.
{
    T acc=0;
    for (unsigned int k=0;k<s;k++)
    {
        acc+=A.get(i,k)*B.get(k,j);
    }
    return acc;
}

#endif
