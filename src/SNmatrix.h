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
#include <cmath>

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


CONST ? 

Many methods are not 'const'. Here is the rationale.

- When you extract an element from the matrix with the method 'getElement', 
  the element contains a reference to the original matrix.

- The matrix itself is not constant because, after a 'PLU' decomposition, the
  matrix itself becomes the 'U' matrix and th original is lost.


 Note : the first point is maybe useless and some 'const' could be added in the
        future.

 Note : you cannot instantiate  
            SNmatrix<const double,7> sn;

When you extract an element, the latter has a reference (in fact a pointer)
to the original matrix.

EQUALITY OF MATRICES

The data of a matrix is encoded in an object
        std::array<T,tp_size*tp_size>
Two matrices are equal when these are equal. The copy constructor is based
on that idea.


PLU DECOMPOSITION

P : permutation matrix
L : lower triangular matrix
U : upper triangular matrix

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
        // the larger element on column 'col' under (or on) the line 'f_line'. 
        SNelement<T,tp_size> getLargerUnder(unsigned int f_line, unsigned int col);
    public:
        SNmatrix();
        SNmatrix(const SNmatrix<T,tp_size>&);

        template <class U,class V,unsigned int s>
        friend bool operator==(const SNmatrix<U,s>&,const SNmatrix<V,s>&);

        unsigned int getSize() const;

        // return a reference to the value of the requested matrix entry.
        T& at(unsigned int,unsigned int);

        // return the matrix element on given (line,column).
        SNelement<T,tp_size> getElement(unsigned int line, unsigned int col);

        // return the larger element (in absolute value) on the given column
        // In case of equality, return the last one (the larger line).
        //   The template type T has to accept arithmetic manipulations
        //   like abs, comparison.
        SNelement<T,tp_size> getLargerOnColumn(unsigned int col);
        // return the largest (absolute value) element under the diagonal
        // on the given column.
        SNelement<T,tp_size> getLargerUnderDiagonal(unsigned int col);

        // swap the lines l1 and l2. This is in-place replacement.
        // The matrix is changed.
        void swapLines(unsigned int l1,unsigned int l2);
};

// CONSTRUCTORS, OPERATORS, ...  -------------------------------------------

template <class T,unsigned int tp_size>
SNmatrix<T,tp_size>::SNmatrix(): data() { };

template <class T,unsigned int tp_size>
SNmatrix<T,tp_size>::SNmatrix(const SNmatrix<T,tp_size>& snm) : data(snm.data)  {};

template <class U,class V,unsigned int s>
bool operator==(const SNmatrix<U,s>& A,const SNmatrix<V,s>& B)
{
    return A.data==B.data;
}

// GETTER METHODS  -------------------------------------------

template <class T,unsigned int tp_size>
unsigned int SNmatrix<T,tp_size>::getSize() const
{
    return size;
};

template <class T,unsigned int tp_size>
SNelement<T,tp_size> SNmatrix<T,tp_size>::getElement(unsigned int line, unsigned int col)
{
    return SNelement<T,tp_size>(line,col,*this);
}

template <class T,unsigned int tp_size>
T& SNmatrix<T,tp_size>::at(const unsigned int i,const unsigned int j)
{
    if (i>tp_size or j>tp_size)
    {
        std::cout<<"This SNmatrix has size "<<tp_size<<". Attempt to access element "<<i<<" , "<<j<<std::endl;
    }
    return data.at(j*tp_size+i);
};


// GAUSS'S PIVOT METHODS

template <class T,unsigned int tp_size>
SNelement<T,tp_size> SNmatrix<T,tp_size>::getLargerUnder(unsigned int f_line, unsigned int col)
{
    T max_val=0;
    unsigned int max_line=0;

    for (unsigned int line=f_line;line<tp_size;line++)
    {
        if (std::abs(at(line,col))>max_val)
        {
            max_val=std::abs(at(line,col));
            max_line=line;
        };
    };
    return getElement(max_line,col);

}

template <class T,unsigned int tp_size>
SNelement<T,tp_size> SNmatrix<T,tp_size>::getLargerOnColumn(unsigned int col) 
{
    return getLargerUnder(0,col);
}

template <class T,unsigned int tp_size>
SNelement<T,tp_size> SNmatrix<T,tp_size>::getLargerUnderDiagonal(unsigned int col) 
{
    return getLargerUnder(col,col);
}


template <class T,unsigned int tp_size>
void SNmatrix<T,tp_size>::swapLines(unsigned int l1, unsigned int l2)
{
    if (l1!=l2)
    {
        for (unsigned int col=0;col<tp_size;col++)
        {
            T tmp = at(l1,col);
            at(l1,col)=at(l2,col);
            at(l2,col)=tmp;
        }
    }
}

#endif
