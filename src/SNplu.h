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
#include "SNmatrices/SNupperTriangularMatrix.h"
#include "SNmatrices/SNpermutation.h"
#include "SNmatrices/Mpermutation.h"


// THE CLASS HEADER -----------------------------------------

template <class T,unsigned int tp_size>
class SNplu
{

    friend SNplu<T,tp_size> SNmatrix<T,tp_size>::getPLU();

    private :
        SNpermutation<T,tp_size> m_P;
        SNlowerTriangularMatrix<T,tp_size> m_L;
        SNupperTriangularMatrix<T,tp_size> m_U;
    public:
        SNupperTriangularMatrix<T,tp_size> getU() const;

        Mpermutation<tp_size> getMpermutation() const;
        SNpermutation<T,tp_size> getSNpermutation() const;
};

// GETTER METHODS -----------------------

template <class T,unsigned int tp_size>
SNupperTriangularMatrix<T,tp_size> SNplu<T,tp_size>::getU() const
{
    return m_U;
}

template <class T,unsigned int tp_size>
Mpermutation<tp_size> SNplu<T,tp_size>::getMpermutation() const
{
    return  m_P.getMpermutation();
}

#endif
