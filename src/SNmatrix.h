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

#ifndef __SNMATRIX_H__153113__
#define __SNMATRIX_H__153113__

#include <array>
#include <iostream>

#include "SNelement.h"

/*
This is my matrix type, designed for numerical computation. It represents a 
square matrix.

TEMPLATE ARGUMENTS

- class T : the type of the matrix entries (typically 'double')
- int tp_size : the size of the matrix (the number of entries is the square of that).

CREATING A MATRIX

Internally, SNmatrix uses a std::array to store the matrix entries, but you should
not try to take advantage of that.

The correct way to populate your matrix is :

            SNmatrix<int,2> sn;
            sn.at(0,0)=1;
            sn.at(0,1)=2;
            sn.at(1,0)=3;
            sn.at(1,1)=4;


Notice that the elements are numbered from 0 to tp_size-1. Not from 1.

DECOMPOSITION

The mathematics can be found (in French) here :
http://laurent.claessens-donadello.eu/pdf/mazhe.pdf


*/


// THE CLASS HEADER -----------------------------------------

template <class T,unsigned int tp_size>
class SNmatrix
{
    private:
        std::array<T,tp_size*tp_size> data;
        unsigned int size=tp_size;
    public:
        SNmatrix();
        unsigned int getSize() const;

        // return a reference to the value of the requested matrix entry.
        T& at(unsigned int,unsigned int);

        // return the matrix element on given (line,column).
        SNelement<T,tp_size> getElement(unsigned int line, unsigned int col);

        // return the larger element (in absolute value) on the given column
        SNelement<T,tp_size> getMaxOnColumn(unsigned int col) const;
};

// IMPLEMENTATIONS  -------------------------------------------


template <class T,unsigned int tp_size>
SNmatrix<T,tp_size>::SNmatrix(): data(){};



template <class T,unsigned int tp_size>
unsigned int SNmatrix<T,tp_size>::getSize() const
{
    return size;
};


template <class T,unsigned int tp_size>
T& SNmatrix<T,tp_size>::at(const unsigned int i,const unsigned int j)
{
    if (i>tp_size or j>tp_size)
    {
        std::cout<<"This SNmatrix has size "<<tp_size<<". Attempt to access element "<<i<<" , "<<j<<std::endl;
    }
    return data.at(j*tp_size+i);
};

template <class T,unsigned int tp_size>
SNelement<T,tp_size> SNmatrix<T,tp_size>::getMaxOnColumn(unsigned int col) const
{
    return getElement(1,col);
}

template <class T,unsigned int tp_size>
SNelement<T,tp_size> SNmatrix<T,tp_size>::getElement(unsigned int line, unsigned int col)
{
    return SNelement<T,tp_size>(line,col,*this);
}

#endif
