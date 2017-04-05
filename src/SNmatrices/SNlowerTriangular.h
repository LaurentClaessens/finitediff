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

#ifndef __SNlowerTriangular_H__103047
#define __SNlowerTriangular_H__103047

#include <array>

#include "SNgeneric.h"
#include "../SNexceptions.cpp"

/*
   This represents a lower triangular matrix (the diagonal can be non zero).
*/

// THE CLASS HEADER -----------------------------------------

template <class T,unsigned int tp_size>
class SNlowerTriangular : public SNgeneric<T,tp_size>
{

    private:
        std::array<T,tp_size*tp_size> data;     // many remain uninitialized
        T _get(const m_num,const m_num) const override;
        T& _at(m_num,m_num) override;
    public :
        std::array<T,tp_size*tp_size> _get_other_data(const SNmatrix<T,tp_size>&) const;
        SNlowerTriangular();
        /** Construct a lower diagonal matrix from a SNmatrix.
         *
         * Whatever was over the diagonal is forgotten.
         * */
        SNlowerTriangular(const SNmatrix<T,tp_size>& A);
        /** Construct a lower diagonal matrix from a generic matrix.
         *
         * Due to the way the data is recorded in `SNmatrix` and 
         * `SNlowerTriangular`, it could be faster to initialize from a
         * `SNmatrix` than from a generic one (if you have the choice ...).
         * */
        SNlowerTriangular(const SNgeneric<T,tp_size>& A);
        /** Construct a lower triangular from a gaussian matrix */
        SNlowerTriangular(const SNgaussian<T,tp_size>& A);
};

// CONSTRUCTOR  ---------------------------------------

template <class T,unsigned int tp_size>
SNlowerTriangular<T,tp_size>::SNlowerTriangular(): data() { };

template <class T,unsigned int tp_size>
SNlowerTriangular<T,tp_size>::SNlowerTriangular(const SNgeneric<T,tp_size>& A)
{
    for (m_num l=0;l<tp_size;++l)
    {
        for (m_num c=0;c<l+1;++c)
        {
            this->at(l,c)=A.get(l,c);
        }
    }
}

// TODO : one has to factorize this function between here and 
// the same in SNupperTriangular.
template <class T,unsigned int tp_size>
std::array<T,tp_size*tp_size> SNlowerTriangular<T,tp_size>::_get_other_data(const SNmatrix<T,tp_size>& A) const
{
    return A.data;
}

template <class T,unsigned int tp_size>
SNlowerTriangular<T,tp_size>::SNlowerTriangular(const SNmatrix<T,tp_size>& A):
    data(_get_other_data(A))
{};

template <class T,unsigned int tp_size>
SNlowerTriangular<T,tp_size>::SNlowerTriangular(const SNgaussian<T,tp_size>& A)
{
    for (m_num i=0;i<tp_size;i++)
    {
        this->at(i,i)=1;
    }
    for (m_num c=0;c<tp_size;++c)
    {
        for (m_num l=c+1;l<tp_size;++l)
        {
            this->at(l,c)=0;
        }
    }
    for (m_num l= A.getColumn()+1;l<tp_size;++l )
    {
        this->at(l,A.getColumn())=A.get(l,A.getColumn());
    }
}

// _GET AND _AT METHODS ---------------------------------------

template <class T,unsigned int tp_size>
T SNlowerTriangular<T,tp_size>::_get(m_num l,m_num c) const
{
    if (l<c)
    {
        return 0;
    }
    return data.at(c*tp_size+l);
}

template <class T,unsigned int tp_size>
T& SNlowerTriangular<T,tp_size>::_at(m_num l,m_num c) 
{
    if (l<c)
    {
        throw SNchangeNotAllowedException(l,c);
    }
    return data.at(c*tp_size+l);
}

#endif
