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

#ifndef __SNgaussian_H__174236
#define __SNgaussian_H__174236

#include <array>

#include "SNgeneric.h"
#include "m_num.h"
#include "../SNexceptions.cpp"


// forward definition
template <class T,unsigned int tp_size>
class SNgeneric;
template <class T>
class SpecialValue;


// THE CLASS HEADER -----------------------------------------

/**
   This represents a gaussian matrix.

   A gaussian matrix is the multiplicative matrix which serves to eliminate
   a column under the diagonal.

   Let \f$A\f$ be a matrix and \f$ c \f$ be an index. Let \f$ m=A_{c,c}\f$. The gaussian matrix for \f$ A \f$ on the index \f$ c \f$ is full of zeros except for
   - \f$ 1\f$ on the diagonal
   - \f$ A_{i,c}/m\f$  on position \f$ (i,c)\f$ when \f$ i>c\f$. 

   That is there is something non trivial only on the column \f$ c \f$, under the diagonal.
 

 TODO : giving 'c' as a template argument, one can in average spare half of the memory by allowing 'c' as size in 'data' instead of 'tp_size'.

*/

template <class T,unsigned int tp_size>
class SNgaussian : public SNgeneric<T,tp_size>
{


    private:
        std::array<T,tp_size> data;     // see implementation of "_at"

    
    /** 

     checkForSpecialElements(i,j)
     checks for element (i,j). If this is a special element 
     (a one whose value is fixed by the fact that we are a gaussian matrix)
     then it returns a 'SpecialValue<T>' with its boolean part set to true,
     meaning that this is a special value.
     In that case :
     - _get returns the corresponding value.
     - _at throws SNchangeNotAllowedException.
    
    If this is not a special value, it returns a 'SpecialValue<T>' with boolean part
    set to false. 
    In that case :
    - both _get and _at have to search in the stored values.

     */

        SpecialValue<T> checkForSpecialElements(const m_num& i,const m_num& j) const;

        //** populate the matrix from the elements of the given matrix */
        template <class U,unsigned int s>
        void populate_from(const SNgeneric<U,s>&);

        /** Construct a matrix from its data. See the implementation of '_at' */
        SNgaussian(const std::array<T,tp_size>& d, const m_num& c);

        T _get(m_num,m_num) const override;
        T& _at(m_num,m_num) override;
    public :
        const m_num column;

        /** Construct a gaussian matrix from a generic one by 
         * - setting 1 on the diagonal (whatever there is in 'A'),
         * - keeping what is below the diagonal on column 'c' 
         * - setting 0 everywhere else  
         *   */
        template <class U,unsigned int s>
        SNgaussian(const SNgeneric<U,s>& A, const m_num& c);

        SNgaussian<T,tp_size> inverse() const;
        
};

// CONSTRUCTOR  ---------------------------------------

template <class T,unsigned int tp_size> 
template <class U, unsigned int s>
void SNgaussian<T,tp_size>::populate_from(const SNgeneric<U,s>& A)
{
    if (s!=tp_size)
    {
        throw IncompatibleMatrixSizeException(tp_size,s);
    }

    T m = A.get(column,column);
    for (m_num i=column+1;i<tp_size;i++)
    {
        this->at(i,column)=-A.get(i,column)/m;
    }
}

template <class T,unsigned int tp_size> 
template<class U,unsigned int s>
SNgaussian<T,tp_size>::SNgaussian(const SNgeneric<U,s>& A , const m_num& c):
    column(c)
{
    populate_from(A);
}

template <class T,unsigned int tp_size> 
SNgaussian<T,tp_size>::SNgaussian(const std::array<T,tp_size>& d, const m_num& c):
    data(d),
    column(c)
{}

// UTILITIES  ---------------------------------------

template <class T,unsigned int tp_size>
SpecialValue<T> SNgaussian<T,tp_size>::checkForSpecialElements(const m_num& i,const m_num& j) const
{
    if (i==j)
    {
        return SpecialValue<T>(1,true);
    }
    if (j!=column)
    {
        return SpecialValue<T>(0,true);
    }
    if (i<j)
    {
        return SpecialValue<T>(0,true);
    }
    return SpecialValue<T>(0,false);
}

// _AT AND _GET METHODS ---------------------------------------



template <class T,unsigned int tp_size>
T SNgaussian<T,tp_size>::_get(m_num i,m_num j) const
{
    SpecialValue<T> sv=checkForSpecialElements(i,j);
    if (sv.special)
    {
        return sv.value;
    }
    return data.at(i-column-1);  //if you change here, you have to change _at
}

template <class T,unsigned int tp_size>
T& SNgaussian<T,tp_size>::_at(m_num i,m_num j) 

    // The elements are stored in 
    //   std::array<T,tp_size> data
    // while the matrix only contains non fixed values on one column, under
    // the diagonal.
    //
    // First element of 'data' stores the first non fixed value.
    // Example for a gaussian matrix on column 1 :
    //
    //  1  0   0  0
    //  0  1   0  0
    //  0  d0  1  0
    //  0  d1  0  1
    //
    //
    // Only the first (tp_size-c-1) elements of 'data' are used.

{
    SpecialValue<T> sv=checkForSpecialElements(i,j);
    if (sv.special)
    {
        throw SNchangeNotAllowedException(i,j);
    }
    return data.at(i-column-1);  //if you change here, you have to change _get
}

// MATHEMATICS  ---------------------------------------


template <class T,unsigned int tp_size>
SNgaussian<T,tp_size> SNgaussian<T,tp_size>::inverse() const
{
    std::array<T,tp_size> new_data(data);
    for (unsigned int k=0;k<tp_size-column-1;++k)
    {
        new_data.at(k)=-new_data.at(k);
    }
    return SNgaussian(new_data,column);
}

#endif
