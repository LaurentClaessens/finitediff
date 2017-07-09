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
* This class represent matrices that are product of gaussian matrices. They
* are "partially" the lower triangular part of the PLU decomposition
* during the process.
*
* Some of the properties of these matrices are
*
* - the inverse is easy
* - they commute easily with the elementary permutation matrices when the
*  indices of the latter are strictly larger than the number of non
*  trivial columns here.
*- The diagonal is filled by 1.
*/
template <class T,unsigned int tp_size>
class SNmultiGaussian : public SNgeneric<T,tp_size>
{
    
    private :
        SNlowerTriangular<T,tp_size> data_L;
        m_num data_last_column;       // the last non trivial column


        //cppcheck-suppress unusedPrivateFunction
        SpecialValue<T> checkForSpecialElements(const m_num&,const m_num&) const;
        T _get(const m_num&, const m_num&) const override;
        T& _at(const m_num&, const m_num&) override;
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
        explicit SNmultiGaussian(const SNgeneric<T,tp_size>& A);
        
        /** 
         * @brief Construct a multi-gaussian from a gaussian one.
         *
         * Since a gaussian matrix is (mathematically) a multi one,
         * the implicit conversion is allowed.
         * */
        //cppcheck-suppress noExplicitConstructor
        SNmultiGaussian(const SNgaussian<T,tp_size>& A);

        /**
         * @brief Initialises as the diagonal matrix full of `x`
         * */
        explicit SNmultiGaussian(const T& x);
        
        //cppcheck-suppress noExplicitConstructor
        SNmultiGaussian(const SNmultiGaussian<T,tp_size>& A);

        void swap(SNmultiGaussian<T,tp_size>&);
        SNmultiGaussian<T,tp_size>& operator=(SNmultiGaussian<T,tp_size>);

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

        /** return the inverse matrix */
        SNmultiGaussian<T,tp_size> inverse() const;

        /*\brief Swap the lines \f$ i \f$ and \f$ j \f$.
         *
         * This only works when \f$ i \f$ and \f$ j \f$ are larger than
         * the last non trivial column.
         *
         * If \f$ i \f$ or \f$ j \f$ is lower or equal to the last non trivial line,
         * throw a `ProbablyNotWhatYouWantException`.
         *
         * In place replacement.
         */
        void swapLines(const m_num& i,const m_num& j);

        /** return the number of the last non trivial column */
        m_num getLastColumn() const;
        /** 
         * \brief Set the number of the last non trivial column
         *
         * Needed for optimization purpose : during the PLU decomposition,
         * one need to modify a multi-gaussian matrix at each step. The point
         * is to *modify* it, not re-creating a new one each time.
         * */
        void setLastColumn(const m_num& lc);

        /** 
         * \brief copies the first `max_l` lines from `other` to `this`.
         *
         * The first `max_l` lines are copied from `other` taking into account :
         * - the elements that cannot be changed here will not be changed.
         *   Typically, the element (0,0) of `other` will be neglected. Also every
         *   diagonal elements are ignored, and all the elements that are beyond
         *   `this->getLastColumn()`.
         * - The size of `other` must be the same.
         * - The template type of `other` must be convertible to the template
         *   type of `this`.
         * - The copy copies each line from column 0 to the diagonal (exclusively).
         *   It does not restrict itself to something like the minimum between
         *   the diagonal and `this->getLastColumn()`. Thus : **the real
         *   last column of the resulting matrix could be incorrect**.
         *
         *   The last point needs attention. One *could* just take care of that
         *   issue during the copying process by looking at non zero elements 
         *   located beyond the `data_last_column` column.
         *   But if these elements are "fake zeroes" like they should be zero but
         *   are, for numerical reasons, non vanishing; in this case we loose a
         *   good occasion to make then real zeroes.
         *
         *   You have to set by hand the last column attribute *before* to 
         *   call `setLastColumn`.
         *
         * \see `setLastColumn`
         *
         * \see copyFirstLines(SNgeneric<U,t>& ans, const SNgeneric<T,tp_size>& A,const m_num& max_l)
         */
        template <class U,unsigned int s>
        void setFirstLines(const SNgeneric<U,s>& other,const m_num& max_l);
};

// CONSTRUCTORS -------------------------------------------------

// from one number
template <class T,unsigned int tp_size>
SNmultiGaussian<T,tp_size>::SNmultiGaussian(const T& x):
    data_L{},
    data_last_column(0) 
{
    if (x!=1)
    {
        throw SNchangeNotAllowedException(0,0,"The one parameter constructor of 'SNmultiGaussian' only works with 1 as agrument, because the other diagonal matrices are not multigaussian.");
    }
}

// from nothing
template <class T,unsigned int tp_size>
SNmultiGaussian<T,tp_size>::SNmultiGaussian():
    data_L(),
    data_last_column(tp_size+1)     //force the user to initialize (see `_at`)
{ }

// from generic
template <class T,unsigned int tp_size>
SNmultiGaussian<T,tp_size>::SNmultiGaussian(const SNgeneric<T,tp_size>& A):
    data_L(A.getGaussian(0)),
    data_last_column(0)
{  }

// from multigaussian
template <class T,unsigned int tp_size>
SNmultiGaussian<T,tp_size>::SNmultiGaussian(const SNmultiGaussian<T,tp_size>& A):
    data_L(A.data_L),
    data_last_column(A.getLastColumn())
{ }

// from gaussian
template <class T,unsigned int tp_size>
SNmultiGaussian<T,tp_size>::SNmultiGaussian(const SNgaussian<T,tp_size>& A):
    data_L(1),       // initiate as the unit matrix
    data_last_column(A.getColumn())
{ 
    for (m_num l=A.getColumn()+1;l<tp_size;++l)
    {
        this->at(l,A.getColumn())=A.get(l,A.getColumn());
    }
}

template <class T,unsigned int tp_size>
void SNmultiGaussian<T,tp_size>::swap(SNmultiGaussian<T,tp_size>& other)
{
    data_L.swap(other.data_L);
    data_last_column.swap(other.data_last_column);
}

// assignation
template <class T,unsigned int tp_size>
SNmultiGaussian<T,tp_size>& SNmultiGaussian<T,tp_size>::operator=(SNmultiGaussian<T,tp_size> other)
{
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

template <class T,unsigned int tp_size>
template <class U,unsigned int s>
void SNmultiGaussian<T,tp_size>::setFirstLines
        (const SNgeneric<U,s>& other,const m_num& max_l)
{
    for (m_num line=0;line<max_l+1;++line)
    {
        for (m_num col=0; col < line ;++col)
        {
            this->at(line,col)=other.get(line,col);
        }
    }
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
    ++data_last_column;
    for (m_num l=other.getColumn()+1;l<tp_size;++l)
    {
        this->at(l,other.getColumn())+=other.get(l,other.getColumn());
    }
}

// MATHEMATICS -------------------------------------------


template <class T,unsigned int tp_size>
SNmultiGaussian<T,tp_size> SNmultiGaussian<T,tp_size>::inverse() const
{
    auto ans=*this;
    for (m_num col=0;col<getLastColumn()+1;++col)
    {
        // one can parallelize with respect to the columns,
        // not with respect to the lines.
        for (m_num line=col+1;line<tp_size;++line)
        {
            T acc=0;
            for (m_num k=col;k<line;++k)
            {
                acc+=this->get(line,k)*ans.get(k,col);
            }
            ans.at(line,col)= - acc;
        }
    }
    return ans;
}


template <class T,unsigned int tp_size>
void SNmultiGaussian<T,tp_size>::swapLines(const m_num& i,const m_num& j)
{
    if (i<=getLastColumn() or j<=getLastColumn())
    {
        throw ProbablyNotWhatYouWantException("You are trying to swap the lines "+std::to_string(i)+" and "+std::to_string(j)+" while the last non trivial column is "+std::to_string(getLastColumn()) );
    }
    for (m_num col=0;col<=getLastColumn();++col)
    {
        std::swap( this->at(i,col),this->at(j,col)  );
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
T SNmultiGaussian<T,tp_size>::_get(const m_num& i,const m_num& j) const
{
    SpecialValue<T> sv=checkForSpecialElements(i,j);
    if (sv.special)
    {
        return sv.value;
    }
    return data_L.get(i,j);  //if you change here, you have to change _at
}

template <class T,unsigned int tp_size>
T& SNmultiGaussian<T,tp_size>::_at(const m_num& i, const m_num& j) 
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
