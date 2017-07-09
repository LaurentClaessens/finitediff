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

#include "../SNpermutation.h"
#include "../SNgaussian.h"
#include "../SNmultiGaussian.h"
#include "../SNidentity.h"
#include "../SNlowerTriangular.h"
#include "../SNupperTriangular.h"
#include "../MathUtilities.h"
#include "../../exceptions/SNexceptions.cpp"


// PRODUCTS ------------------------------------------


// SUM ---------------------------------------

/** 
 * Sum of two SNmatrix.
 *
 * The return type is the one of the left argument.
 * THUS : this is *not* totally commutative. You may have
 * \f$ A+B\neq B+A \f$
 * */

template <class U,unsigned int s,class V,unsigned int t>
SNmatrix<U,s> operator+(const SNmatrix<U,s>& A,const SNmatrix<V,t>& B)
{

    // TODO : since this sum is not commutative, maybe one has to implement it
    // as member function.

    SNmatrix<U,s> new_matrix(A);
    for (unsigned int k=0;k<s*s;k++)
    {
        new_matrix.data.at(k)+=B.data.at(k);
    }
    return new_matrix;
}

// DIFFERENCE ---------------------------------------

template <class U,unsigned int s,class V,unsigned int t>
SNmatrix<U,s> operator-(const SNgeneric<U,s>& A,const SNidentity<V,t>& B)
{
    checkSizeCompatibility(A,B);
    SNmatrix<U,s> new_matrix(A);
    for (m_num k=0;k<s;k++)
    {
        new_matrix.at(k,k)-=1;
    }
    return new_matrix;
}


// EQUALITIES ---------------------------------------

template <class U,unsigned int s,class V,unsigned int t>
bool operator==(const SNgeneric<U,s>& A,const SNgeneric<V,t>& B)
{
    return componentWiseeEquality(A,B);
}

template <class U,unsigned int s,class V,unsigned int t>
bool operator==(const SNmatrix<U,s>& A,const SNmatrix<V,t>& B)
{
    checkSizeCompatibility(A,B);
    return A.data==B.data;
}

template <class U,unsigned int s,class V,unsigned int t>
bool operator==(const SNmatrix<U,s>& A,const SNupperTriangular<V,t>& B)
{
    return B==A;
}


/** 
 *\brief Return true if the two permutations are equal.
 * */
template <unsigned int tp_size>
bool operator==(const Mpermutation<tp_size>& A, const Mpermutation<tp_size>& B) 
{
    return A.data==B.data;
}

/** 
 *\brief Return true if the two permutations are equal.
 * */
template <unsigned int tp_size>
bool operator==(const MgenericPermutation<tp_size>& A, const MgenericPermutation<tp_size>& B) 
{
    for (unsigned int k=0;k<tp_size;++k)
    {
        if (A(k)!=B(k)   )
        {
            return false;
        }
    }
    return true;
}

#endif
