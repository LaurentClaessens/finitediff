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

#ifndef __SNlowerTriangularMatrix_H__103047
#define __SNlowerTriangularMatrix_H__103047

#include <array>

#include "SNgeneric.h"
#include "../SNexceptions.cpp"

/*
   This represents a lower triangular matrix (the diagonal can be non zero).
*/

// THE CLASS HEADER -----------------------------------------

template <class T,unsigned int tp_size>
class SNlowerTriangularMatrix : public SNgeneric<T,tp_size>
{

    private:
        std::array<T,tp_size*tp_size> data;     // many remain uninitialized
        T _get(const unsigned int,const unsigned int) const;
    public :
        SNlowerTriangularMatrix();

        unsigned int getSize() const;

};

// CONSTRUCTOR  ---------------------------------------

template <class T,unsigned int tp_size>
SNlowerTriangularMatrix<T,tp_size>::SNlowerTriangularMatrix(): data() { };

// GETTER METHODS ---------------------------------------

template <class T,unsigned int tp_size>
unsigned int SNlowerTriangularMatrix<T,tp_size>::getSize() const
{
    return tp_size;
};

template <class T,unsigned int tp_size>
T SNlowerTriangularMatrix<T,tp_size>::_get(unsigned int l,unsigned int c) const
{
    if (l<c)
    {
        return 0;
    }
    return data.at(c*tp_size+l);
}

template <class T,unsigned int tp_size>
T& SNlowerTriangularMatrix<T,tp_size>::at(unsigned int i,unsigned int j) 
{
    if (i<j)
    {
        throw SNchangeNotAllowedException(i,j);
    }
    return get(i,j);
}


#endif
