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
#include "SNidentity.h"
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
        /** 
         * The no-parameter constructor initializes `data_last_column` to
         * `tp_size+1`, and this is detected by `_at` so that you cannot
         * populate the matrix before to have initialized the maximal
         * non trivial column number.
         *
         * See `setLastColumn`.
         * */
        SNmultiGaussian();
       /**
        * Construct gaussian matrix of the first line of the argument `A` 
        *
        * A gaussian matrix is a particular case of multi-gaussian matrix.
        * */
        SNmultiGaussian(const SNgeneric<T,tp_size>& A);
        SNmultiGaussian(const SNgaussian<T,tp_size>& A);
        SNmultiGaussian(const SNmultiGaussian<T,tp_size>& A);

        void swap(SNmultiGaussian<T,tp_size>&);
        SNmultiGaussian<T,tp_size>& operator=(const SNmultiGaussian<T,tp_size>);

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

        /** return the number of the last non trivial column */
        m_num getLastColumn() const;
        /** Set the number of the last non trivial column
         *
         * Needed for optimization purpose : during the PLU decomposition,
         * one need to modify a multi-gaussian matrix at each step. The point
         * is to *modify* it, not re-creating a new one each time.
         * */
        void setLastColumn(const m_num& lc);
};

// CONSTRUCTORS -------------------------------------------------
//
// from nothing
template <class T,unsigned int tp_size>
SNmultiGaussian<T,tp_size>::SNmultiGaussian():
    data_L(SNidentity<T,tp_size>()),
    data_last_column(tp_size+1)     //force the user to initialize (see `_at`)
{ }

// from generic
template <class T,unsigned int tp_size>
SNmultiGaussian<T,tp_size>::SNmultiGaussian(const SNgeneric<T,tp_size>& A):
    data_L(A.getGaussian(0)),
    data_last_column(0)
{ }

// from multigaussian
template <class T,unsigned int tp_size>
SNmultiGaussian<T,tp_size>::SNmultiGaussian(const SNmultiGaussian<T,tp_size>& A):
    data_L(A.data_L),
    data_last_column(A.getLastColumn())
{ 
    debug_print<<"constructeur SNmultiGaussian"<<std::endl;
    debug_print<<"this -> last : "<<this->getLastColumn()<<std::endl;
    debug_print<<"A-> last : "<<A.getLastColumn()<<std::endl;
}

// from gaussian
template <class T,unsigned int tp_size>
SNmultiGaussian<T,tp_size>::SNmultiGaussian(const SNgaussian<T,tp_size>& A):
    data_last_column(A.getColumn())
{ 
    for (m_num c=0;c<A.getColumn();++c)
    {
        for (m_num l=c+1;l<tp_size;++l)
        {
            this->at(l,c)=A.get(l,c);
        }
    }
    for (m_num l=A.getColumn()+1;l<tp_size;++l)
    {
        this->at(l,A.getColumn())=A.get(l,A.getColumn());
    }
}

template <class T,unsigned int tp_size>
void SNmultiGaussian<T,tp_size>::swap(SNmultiGaussian<T,tp_size>& other)
{
    data_L.swap(other.data_L);
    data_last_column.swap(other.data_L);
}

// assignation
template <class T,unsigned int tp_size>
SNmultiGaussian<T,tp_size>& SNmultiGaussian<T,tp_size>::operator=(const SNmultiGaussian<T,tp_size> other)
{
    debug_print<<"je passe par operator= dans SNmultiGaussian"<<std::endl;
    swap(other);
    return *this;
}

// GETTER/SETTER METHODS  ---------------------------------------


template <class T,unsigned int tp_size>
m_num SNmultiGaussian<T,tp_size>::getLastColumn() const
{
    return data_last_column;
}

template <class T,unsigned int tp_size>
void SNmultiGaussian<T,tp_size>::setLastColumn(const m_num& lc)
{
    if (lc>tp_size-1)   // makes no sense to have a gaussian 
                        // behaviour on the last line.
    {
        throw OutOfRangeColumnNumber("The specified column number is larger than the size of the matrix.");
    }
    data_last_column=lc;
}

// OPERATORS  ---------------------------------------

template <class T,unsigned int tp_size>
void SNmultiGaussian<T,tp_size>::operator *=(const SNgaussian<T,tp_size>& other)
{
    checkSizeCompatibility(*this,other);
    if (other.getColumn()!=data_last_column+1)
    {
        throw ProbablyNotWhatYouWantException("You are trying to multiply a multi-Gaussian matrix by a gaussian matrix whose column is not the next one. This is mathematically possible, but probably not what you want. However; this situation is not yet implemented.");
    }
    data_last_column++;
    for (m_num l=other.getColumn()+1;l<tp_size;++l)
    {
        this->at(l,other.getColumn())+=other.get(l,other.getColumn());
    }
}

// UTILITIES  ---------------------------------------


/**
 * A `SpecialValue` represents a value in the matrix that is fixed by
 * the type of the matrix (like "0" at position (1,4) for a lower triangular
 * matrix).
 *
 * This is merely a pair "value,boolean". 
 * The boolean says if the requested element is special. 
 * - if true, the value is the special value
 * - if false, the value is dummy because it means that the requested element
 *   is not special in this kind of matrix.
 * */

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
    if (data_last_column==tp_size+1)
    {
        throw NotInitializedMemberException("You are trying to populate a 'SNmultiGaussian' before to initialize the member 'data_last_column'. Use setLastColumn().");
    }
    SpecialValue<T> sv=checkForSpecialElements(i,j);
    if (sv.special)
    {
        throw SNchangeNotAllowedException(i,j);
    }
    return data_L.at(i,j);  //if you change here, you have to change _get
}

#endif
