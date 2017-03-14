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
#include "SNmatrices/SNpermutation.h"


// THE CLASS HEADER -----------------------------------------

template <class T,unsigned int tp_size>
class SNplu
{

    friend SNplu<T,tp_size> SNmatrix<T,tp_size>::getPLU();

    private :
        SNpermutation<tp_size> permutations;
        SNmatrix<T,tp_size> m_L;
        SNmatrix<T,tp_size>& m_U;
    public:
        SNplu(SNmatrix<T,tp_size>& original);

        SNmatrix<T,tp_size>& getU();
        SNpermutation<tp_size> getPermutation() const;

};

// CONSTRUCTORS ---------------------------------------------- 

template <class T,unsigned int tp_size>
SNplu<T,tp_size>::SNplu(SNmatrix<T,tp_size>& original) : 
    permutations(),
    m_L(),
    m_U(original)
{}

// GETTER METHODS -----------------------

template <class T,unsigned int tp_size>
SNmatrix<T,tp_size>& SNplu<T,tp_size>::getU()
{
    return m_U;
}

template <class T,unsigned int tp_size>
SNpermutation<tp_size> SNplu<T,tp_size>::getPermutation() const
{
    return  permutations;
}

#endif
