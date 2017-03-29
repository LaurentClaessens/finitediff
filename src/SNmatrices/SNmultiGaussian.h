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


#include "SNlower.h"
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
class SNmultiGausian : public SNgeneric<T,tp_size>
{
    private :
        SNlowerTriangular<T,tp_size> data_L;
        m_num data_c;       // the last non trivial column

    public:
        // Initialize with a matrix A as argument produce the gaussian
        // matrix of A for column 0

}

#endif


