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
        SpecialValue<T> checkForSpecialElements(const m_num&,const m_num&) const;
        T _get(m_num,m_num) const override;
        T& _at(m_num,m_num) override;
    public:
        /** Construct gaussian matrix of the argument `A` */
        SNmultiGaussian(const SNgeneric<T,tp_size>& A);

        /** return the last non trivial column */
        m_num lastColumn() const;

        /**
         * The product \f$ AB \f$ is easy when \f$ A \f$ is
         * multigaussian and \f$ B \f$ is a gaussian matrix
         * for the next line.
         *
         * `A*=B` does
         * - `A=A*B` if A is multigaussian with last non trivial column 
         *   \f$ l_l \f$ and if B is gaussian for the column \f$ l_c+1 \f$.
         * - throw `IncompatibleMatrixSizeException` exception if the size are
         *   not the same.
         * - throw `ProbablyNotWhatYouWantException` if the requirements about
         *   the column are not fulfilled.
         * */
        void operator *=(const SNgaussian<T,tp_size>& other);
};

// CONSTRUCTORS -------------------------------------------------

template <class T,unsigned int tp_size>
SNmultiGaussian<T,tp_size>::SNmultiGaussian(const SNgeneric<T,tp_size>& A):
    data_last_column(0),
    data_L(A.getGaussian(0))
{ }

// GETTER METHODS  ---------------------------------------


template <class T,unsigned int tp_size>
m_num SNmultiGaussian<T,tp_size>::lastColumn() const
{
    return data_last_column;
}

// OPERATORS  ---------------------------------------

template <class T,unsigned int tp_size>
void SNmultiGaussian<T,tp_size>::operator *=(const SNgaussian<T,tp_size>& other)
{
    checkSizeCompatibility(*this,other);
    if (other.column!=data_last_column+1)
    {
        throw ProbablyNotWhatYouWantException("You are trying to multiply a multi-Gaussian matrix by a gaussian matrix whose column is not the next one. This is mathematically possible, but probably not what you want. However; this situation is not yet implemented.");
    }
    data_last_column++;
    for (m_num l=other.column+1;l<tp_size;++l)
    {
        this->at(l,other.column)+=other.get(l,other.column);
    }
}

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
    return data_L.get(i,j);  //if you change here, you have to change _at
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
