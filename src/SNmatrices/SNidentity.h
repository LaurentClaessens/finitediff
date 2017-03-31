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

#ifndef __SNIDENTITY_H__092622__
#define __SNIDENTITY_H__092622__

#include "SNgeneric.h"

/*
   This class represent the identity matrix.
*/
template <class T,unsigned int tp_size>
class SNidentity  : public SNgeneric<T,tp_size>
{

    private:
        T& _at(const m_num,const m_num) override;
        T _get(const m_num,const m_num) const override;

    public:
        SNidentity();
};

// CONSTRUCTORS, OPERATORS, ...  -------------------------------------------

template <class T,unsigned int tp_size>
SNidentity<T,tp_size>::SNidentity() { };

// _GET AND _AT METHODS ---------------------------

template <class T,unsigned int tp_size>
T& SNidentity<T,tp_size>::_at(const m_num i,const m_num j) 
{
        throw SNchangeNotAllowedException(i,j);
};

template <class T,unsigned int tp_size>
T SNidentity<T,tp_size>::_get(const m_num i,const m_num j) const
{
    if (i==j)
    {
        return 1;
    }
    return 0;
};

#endif
