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
        const Mpermutation<tp_size> data_P; 
        const SNlowerTriangular<T,tp_size> data_L;
        const SNupperTriangular<T,tp_size> data_U;
    public:

        /** @brief constructor from the already computed P,L and U.
         * */
        SNplu(const Mpermutation<tp_size>& mP,const SNlowerTriangular<T,tp_size>& mL,const SNupperTriangular<T,tp_size>& mU);

        const SNpermutation<T,tp_size> getP() const;
        const SNlowerTriangular<T,tp_size> getL() const;
        const SNupperTriangular<T,tp_size> getU() const;
        const Mpermutation<tp_size> getMpermutation() const;
};

// CONSTRUCTORS -----------------------

template <class T,unsigned int tp_size>
SNplu<T,tp_size>::SNplu(const Mpermutation<tp_size>& mP,const SNlowerTriangular<T,tp_size>& mL,const SNupperTriangular<T,tp_size>& mU):
    data_P(mP),
    data_L(mL),
    data_U(mU)
{}

// GETTER METHODS -----------------------

template <class T,unsigned int tp_size>
const SNpermutation<T,tp_size> SNplu<T,tp_size>::getP() const
{
    return SNpermutation<T,tp_size>(data_P);
}

template <class T,unsigned int tp_size>
const SNlowerTriangular<T,tp_size> SNplu<T,tp_size>::getL() const
{
    return data_L;
}


template <class T,unsigned int tp_size>
const SNupperTriangular<T,tp_size> SNplu<T,tp_size>::getU() const
{
    return data_U;
}


template <class T,unsigned int tp_size>
const Mpermutation<tp_size> SNplu<T,tp_size>::getMpermutation() const
{
    return  data_P;
}

#endif
