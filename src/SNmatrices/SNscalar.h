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

#ifndef __SNscalar_H__9561__
#define __SNscalar_H__9561__

#include <array>

#include "SNgeneric.h"

/** 
 * @brief A scalar matrix.
 *
 * A scalar matrix is a multiple of identity, that is a diagonal matrix
 * with the diagonal full of the same number.
 * */
template <class T,unsigned int tp_size>
class SNscalar : public SNgeneric<T,tp_size>
{
    private :
        T data; // the multiple.

        T _get(const m_num,const m_num) const override;

        /** 
         * @brief Return a reference to the scalar defining the matrix.
         *
         * Warning : modifying the element \f$ (3,3) \f$ in fact modify
         * in the same time all the elements \f$ (k,k) \f$.
         * */
        T& _at(m_num,m_num) override;
    public:
        SNscalar();
        SNscalar(const T& x);
};

// CONSTRUCTORS --------------------------------------------

// from nothing
template <class T,unsigned int tp_size>
SNscalar<T,tp_size>::SNscalar(): 
    data()
{}

// from a number
template <class T,unsigned int tp_size>
SNscalar<T,tp_size>::SNscalar(const T& x): 
    data(x)
{}



// _GET AND _AT METHODS ---------------------------------------

template <class T,unsigned int tp_size>
T SNscalar<T,tp_size>::_get(m_num l,m_num c) const
{
    if (l==c)
    {
        return data;
    }
    return 0;
}

template <class T,unsigned int tp_size>
T& SNscalar<T,tp_size>::_at(m_num l,m_num c) 
{
    if (l!=c)
    {
        throw SNchangeNotAllowedException(l,c);
    }
    return data;
}


#endif
