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
#include "../exceptions/SNexceptions.cpp"

// forward definition
template <class T,unsigned int tp_size>
class SNmultiGaussian;

// THE CLASS HEADER -----------------------------------------

/**
* \brief Represents a lower triangular matrix (the diagonal can be non zero).
*/
template <class T,unsigned int tp_size>
class SNlowerTriangular : public SNgeneric<T,tp_size>
{

    private:
        std::array<T,tp_size*tp_size> data;     // many remain uninitialized
        T _get(const m_num&, const m_num&) const override;

        /** 
         * \brief Return by reference the content of element (`l`,`c`) of 
         * the matrix.
         *
         * If `l<c`, throws `SNnonAllowedChange` because  `_at` is destinated
         * to populate the matrix.
         * */
        T& _at(m_num l,m_num c) override;
    public :
        std::array<T,tp_size*tp_size> _get_other_data(const SNmatrix<T,tp_size>&) const;

        /** 
         * @brief Construct a lower triangular matrix with non initialized entries.
         * */
        SNlowerTriangular();

        /** 
         * @brief Initializes as a diagonal matrix full of `x`.
         * */
        explicit SNlowerTriangular(const T& x);

        /** \brief Construct a lower diagonal matrix from a SNmatrix.
         *
         * Whatever was over the diagonal is forgotten.
         * */
        explicit SNlowerTriangular(const SNmatrix<T,tp_size>& A);

        /** \brief Construct a lower diagonal matrix from a generic matrix.
         *
         * Due to the way the data is recorded in `SNmatrix` and 
         * `SNlowerTriangular`, it could be faster to initialize from a
         * `SNmatrix` than from a generic one (if you have the choice ...).
         * */
        explicit SNlowerTriangular(const SNgeneric<T,tp_size>& A);

        /** Construct a lower triangular from a gaussian matrix */
        //cppcheck-suppress noExplicitConstructor
        SNlowerTriangular(const SNgaussian<T,tp_size>& A);
        /** Construct a lower triangular from a multi-gaussian matrix */
        //cppcheck-suppress noExplicitConstructor
        SNlowerTriangular(const SNmultiGaussian<T,tp_size>& A);

        void swap(SNlowerTriangular<T,tp_size>& other);
};

// CONSTRUCTOR  ---------------------------------------

// from nothing
template <class T,unsigned int tp_size>
SNlowerTriangular<T,tp_size>::SNlowerTriangular():
    data()
{}

// from SNlowerTriangular
template <class T,unsigned int tp_size>
SNlowerTriangular<T,tp_size>::SNlowerTriangular(const T& x):
    data{}     // initialize full of zeroes.
{
    for (m_num i=0;i<tp_size;++i)
    {
        this->at(i,i)=x;
    }
}

// from SNgeneric
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
SNlowerTriangular<T,tp_size>::SNlowerTriangular(const SNgaussian<T,tp_size>& A):
    data{} // initialized full of zeroes.
{
    for (m_num i=0;i<tp_size;++i)
    {
        this->at(i,i)=1;
    }
    for (m_num l= A.getColumn()+1;l<tp_size;++l )
    {
        this->at(l,A.getColumn())=A.get(l,A.getColumn());
    }
}

// from multi-gaussian
template <class T,unsigned int tp_size>
SNlowerTriangular<T,tp_size>::SNlowerTriangular(const SNmultiGaussian<T,tp_size>& A):
    data{} // initialized full of zeroes.
{
    for (m_num i=0;i<tp_size;++i)
    {
        this->at(i,i)=1;
    }
    m_num last_column(A.getLastColumn());
    for (m_num col=0;col<last_column;++col)
    {
        for (m_num l= col+1;l<tp_size;++l )
        {
            this->at(l,col)=A.get(l,col);
        }
    }
}

template <class T,unsigned int tp_size>
void SNlowerTriangular<T,tp_size>::swap(SNlowerTriangular<T,tp_size>& other)
{
    std::swap(data,other.data);
}

// _GET AND _AT METHODS ---------------------------------------

template <class T,unsigned int tp_size>
T SNlowerTriangular<T,tp_size>::_get(const m_num& l,const m_num& c) const
{
    if (l<c)
    {
        return 0;
    }
    return data.at(c*tp_size+l);
}


template <class T,unsigned int tp_size>
T& SNlowerTriangular<T,tp_size>::_at(const m_num& l,const m_num& c) 
{
    if (l<c)
    {
        throw SNchangeNotAllowedException(l,c);
    }
    return data.at(c*tp_size+l);
}

#endif
