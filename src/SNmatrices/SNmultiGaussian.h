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

#ifndef __SNMULTIGAUSSIAN_H__105525
#define __SNMULTIGAUSSIAN_H__105525


#include "SNlowerTriangular.h"
#include "m_num.h"

/** 
This class represent matrices that are product of gaussian matrices. They
are "partially" the lower triangular part of the PLU decomposition
during the process.

Some of the properties of these matrices are
- the inverse is easy
- they commute easily with the elementary permutation matrices when the
  indices of the latter are strictly larger than the number of non
  trivial columns here.
- The diagonal is filled by 1.

 */

template <class T,unsigned int tp_size>
class SNmultiGaussian : public SNgeneric<T,tp_size>
{
    private :
        SNlowerTriangular<T,tp_size> data_L;
        m_num data_last_column;       // the last non trivial column
    public:
        /** Construct gaussian matrix of the argument `A` */
        SNmultiGaussian(const SNgeneric<T,tp_size>& A);
};

// CONSTRUCTORS -------------------------------------------------

template <class T,unsigned int tp_size>
SNmultiGaussian<T,tp_size>::SNmultiGaussian(const SNgeneric<T,tp_size>& A):
    data_L(A.getGaussian(0)),
    data_last_column(0)
{ }

// UTILITIES  ---------------------------------------

template <class T,unsigned int tp_size>
SpecialValue<T> SNmultiGaussian<T,tp_size>::checkForSpecialElements(const m_num& i,const m_num& j) const
{
    if (i==j)
    {
        return SpecialValue<T>(1,true);
    }
    if (j>data_last_column)
    {
        return SpecialValue<T>(0,true);
    }
    if (i<j)
    {
        return SpecialValue<T>(0,true);
    }
    return SpecialValue<T>(0,false);      // In this case, the value "0" is dummy.
}


// _GET AND _AT METHODS ---------------------------------------

template <class T,unsigned int tp_size>
T SNmultiGaussian<T,tp_size>::_get(m_num i,m_num j) const
{
    SpecialValue<T> sv=checkForSpecialElements(i,j);
    if (sv.special)
    {
        return sv.value;
    }
    return data_L.at(i,j);  //if you change here, you have to change _at
}

template <class T,unsigned int tp_size>
T& SNmultiGaussian<T,tp_size>::_at(m_num i,m_num j) 
{
    SpecialValue<T> sv=checkForSpecialElements(i,j);
    if (sv.special)
    {
        throw SNchangeNotAllowedException(i,j);
    }
    return data_L.at(i,j);  //if you change here, you have to change _get
}



#endif
