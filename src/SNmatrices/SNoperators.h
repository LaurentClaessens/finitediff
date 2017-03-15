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
 This file contains operators *,== between different types of matrices.
*/

#ifndef  OPERATORS_H__064802_
#define  OPERATORS_H__064802_

#include "SNgaussianMatrix.h"
#include "SNlowerTriangularMatrix.h"
#include "MathUtilities.h"
#include "../SNexceptions.cpp"



// PRODUCTS ------------------------------------------


// SNmatrix*SNmatrix.
// cannot do better than compute everything.

template <class U,class V,unsigned int s,unsigned int t>
SNmatrix<U,s> operator*(const SNmatrix<U,s>& A, const SNmatrix<V,t>& B)
{
    if (s!=t)
    {
        throw IncompatibleMatrixSizeException(s,t);
    }
    SNmatrix<U,s> ans;
    for (unsigned int i=0;i<s;i++)
    {
        for (unsigned int j=0;j<s;j++)
        {
            ans.at(i,j)=matrixProductComponent(A,B,i,j,s);
        }
    }
    return ans;   //relies on RVO.
}

// SNgaussianMatrix * SNmatrix
// can copy the first 'c' lines.

template <class U,class V,unsigned int s,unsigned int t>
SNmatrix<U,s> operator*(const SNgaussianMatrix<U,s>& A, const SNmatrix<V,t>& B)
{
    if (s!=t)
    {
        throw IncompatibleMatrixSizeException(s,t);
    }
    SNmatrix<U,s> ans;
    for (unsigned int i=c;i<c+1;i++)
    {
        for (unsigned int j=0;j<s;j++)
        {
            ans.at(i,j)=B.get(i,j);
        }
    }
    for (unsigned int i=c+1;i<s;i++)
    {
        for (unsigned int j=0;j<s;j++)
        {
            ans.at(i,j)=matrixProductComponent(A,B,i,j,s);
        }
    }
    return ans;  
}



// EQUALITIES ---------------------------------------

template <class U,class V,unsigned int s,unsigned int t>
bool operator==(const SNmatrix<U,s>& A,const SNmatrix<V,t>& B)
{
    return A.data==B.data;
}

template <class U,class V,unsigned int s,unsigned int t>
bool operator==(const SNgaussianMatrix<U,s>& G,const SNmatrix<V,t>& A)
{
    return componentWiseeEquality(G,A);
}

template <class U,class V,unsigned int s,unsigned int t>
bool operator==(const SNlowerTriangularMatrix<U,s>& G,const SNmatrix<V,t>& A)
{
    return componentWiseeEquality(G,A);
}


#endif
