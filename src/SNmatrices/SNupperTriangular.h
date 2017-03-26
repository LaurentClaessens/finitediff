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

#ifndef __SNupperTriangular_H__103047
#define __SNupperTriangular_H__103047

#include <array>

#include "SNgeneric.h"
#include "../SNexceptions.cpp"

/*
   This represents a upper triangular matrix (the diagonal can be non zero).
*/

// THE CLASS HEADER -----------------------------------------

template <class T,unsigned int tp_size>
class SNupperTriangular : public SNgeneric<T,tp_size>
{
    private:
        std::array<T,tp_size*tp_size> data;     // many remain uninitialized

        T _get(const m_num,const m_num) const override;
        T& _at(m_num,m_num) override;
    public :
        /**
         * return `data` of the requested matrix. 
         * Since I cannot declare the constructor of SNupperTriangular
         * being a friend of SNmatrix (templated constructor ...), 
         * the constructor
         *  SNupperTriangular(const SNmatrix<T,tp_size>& A);
         * will use `_get_other_data` to get the `data` member of `A`.
         * */
        std::array<T,tp_size*tp_size> _get_other_data(const SNmatrix<T,tp_size>&) const;

        SNupperTriangular();
        SNupperTriangular(const SNmatrix<T,tp_size>& A);
};

// CONSTRUCTOR  ---------------------------------------

template <class T,unsigned int tp_size>
SNupperTriangular<T,tp_size>::SNupperTriangular(): 
    data() 
{ };


template <class T,unsigned int tp_size>
std::array<T,tp_size*tp_size> SNupperTriangular<T,tp_size>::_get_other_data(const SNmatrix<T,tp_size>& A) const
{
    return A.data;
}

template <class T,unsigned int tp_size>
SNupperTriangular<T,tp_size>::SNupperTriangular(const SNmatrix<T,tp_size>& A): 
    data(_get_other_data(A)) 
{};

// _GET AND _AT METHODS ---------------------------------------

template <class T,unsigned int tp_size>
T SNupperTriangular<T,tp_size>::_get(m_num l,m_num c) const
{
    if (l>c)
    {
        return 0;
    }
    return data.at(c*tp_size+l);
}

template <class T,unsigned int tp_size>
T& SNupperTriangular<T,tp_size>::_at(m_num l,m_num c) 
{
    if (l>c)
    {
        throw SNchangeNotAllowedException(l,c);
    }
    return data.at(c*tp_size+l);
}

#endif
