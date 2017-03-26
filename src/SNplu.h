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

#ifndef __SNPLU_H__142039__
#define __SNPLU_H__142039__


#include "SNmatrices/SNmatrix.h"
#include "SNmatrices/SNupperTriangular.h"
#include "SNmatrices/SNpermutation.h"
#include "SNmatrices/Mpermutation.h"


// THE CLASS HEADER -----------------------------------------

template <class T,unsigned int tp_size>
class SNplu
{

    friend SNplu<T,tp_size> SNmatrix<T,tp_size>::getPLU() const;

    private :
        Mpermutation<tp_size> m_P;              // we store the "mathematical" permutation, not the matrix.
        SNlowerTriangular<T,tp_size> m_L;
        SNupperTriangular<T,tp_size> m_U;

        /**  '_setU' takes a matrix and says that this is the 'U' one.
           - do not check that A is actually upper diagonale
           - erase the old U matrix.
           - only to be used at the end of the building process of the PLU
             decomposition, when the initial matrix is transformed into
             an upper diagonal one.
           */
        void _setU(const SNmatrix<T,tp_size>& A);
    public:

        SNupperTriangular<T,tp_size> getU() const;

        Mpermutation<tp_size> getMpermutation() const;
        SNpermutation<T,tp_size> getSNpermutation() const;
};

// CONSTRUCTORS -----------------------


template <class T,unsigned int tp_size>
void SNplu<T,tp_size>::_setU(const SNmatrix<T,tp_size>& A) 
{
    // non optimal because one could try to copy only the useful part.
    m_U=A;
}

// GETTER METHODS -----------------------

template <class T,unsigned int tp_size>
SNupperTriangular<T,tp_size> SNplu<T,tp_size>::getU() const
{
    return m_U;
}

template <class T,unsigned int tp_size>
Mpermutation<tp_size> SNplu<T,tp_size>::getMpermutation() const
{
    return  m_P;
}

#endif
