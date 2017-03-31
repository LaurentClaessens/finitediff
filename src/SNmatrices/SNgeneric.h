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

#ifndef __SNGENERIC_H__142708_
#define __SNGENERIC_H__142708_

#include <cmath>

#include "SNgaussian.h"
#include "SNline.h"
#include "m_num.h"

#include "MathUtilities.h"
#include "../DebugPrint.h"


// forward
template <class T,unsigned int tp_size>
class SNgaussian;

/**
 
   This is the base class for the other matrices types. 
   A matrix can be
   - diagonal
   - upper/lower triangular
   - a gaussian transformation
   - symmetric
   - etc.

   Each of these types have different implementations for storing the elements, 
   perform equality tests (`opertor==`), multiplication (`operator*`) and so on.

 
   What you need in your subclass :

   - override  `T  _get`
     - `_get(i,j)` has to return the value of the element (i,j) of your matrix.

   - override  `T& _at`
    - `_at(i,j)` has to return a _reference_ to the element (i,j) of your matrix.

    - In consequence, `_at(i,j)` should throw 'SNchangeNotAllowedException' if the requested element cannot be modified. For example requesting the element (1,4) of a lower triangular matrix.
            
    - `_at` is intended to populate the matrix.
            
    
As an example, on a lower diagonal matrix,
- `_get(1,3)` returns 0 (by value)
- `_at(1,3)` throws SNchangeNotAllowedException


*/

template <class T,unsigned int tp_size>
class SNgeneric
{
    private :

        /** Return by value the requested element of the matrix */
        virtual T _get(const m_num,const m_num) const=0;

        /** Return by reference the requested element of the matrix */
        virtual T& _at(const m_num,const m_num)=0;


        /**
         throws `SNoutOfRangeException` if the requested element is out of
         range (larger than 'tp_size').
         */
        virtual void checkRangeCorectness(const m_num&,const m_num&) const final;
    public:
        virtual unsigned int getSize() const final;

        virtual T& at(const m_num,const m_num) final;
        virtual T get(const m_num&,const m_num&) const final;

        virtual SNline<T,tp_size> getSNline(m_num l) const;

        /** 
         subtract the given matrix from this matrix.
         This is in-place replacement. Thus the least const in the 
         world.
          */
        template <class V,unsigned int s>
        void subtract(const SNgeneric<V,s>&);
        template <class V,unsigned int s>
        void subtract(const SNgaussian<V,s>&);

        /** 
         * return the gaussian matrix for the requested column 'c'
         * */
        SNgaussian<T,tp_size> getGaussian(const m_num c) const;

        /** 
        numerical equality test 'up to epsilon'.
        If the max norm of "this-A" is strictly larger than epsilon,
        return false.
        Else return true.

        For this reason, the template parameter `T` has to support `std::abs`
         */
        template <class V,unsigned int s>
        bool isNumericallyEqual(const SNgeneric<V,s>& A,const double& epsilon) const;
};


// UTILITIES --------------------------

template <class T>
class SpecialValue
{
    public :
        T value;
        bool special;
        SpecialValue(T v,bool s) :
            value(v),
            special(s)
        {}
};

template <class T,unsigned int tp_size>
void SNgeneric<T,tp_size>::checkRangeCorectness(const m_num& l,const m_num& c) const 
{
    if (l>tp_size or c>tp_size)
    {
        throw SNoutOfRangeException(l,c,getSize());
    }
}

template <class T,unsigned int tp_size>
SNline<T,tp_size> SNgeneric<T,tp_size>::getSNline(m_num l) const
{
    std::array<T,tp_size> al;
    for (m_num c=0;c<tp_size;c++)
    {
        al.at(c)=this->get(l,c);
    }
    return SNline<T,tp_size>(al);
}

// OPERATORS ------------------------------

template <class V,unsigned int s>
std::ostream& operator<<(std::ostream& stream,const SNgeneric<V,s>& snm)
{
    for (m_num l=0;l<s;l++)
    {
        stream<<snm.getSNline(l)<<std::endl;
    }
    return stream;
}

// GET SIZE ------------------------------

template <class T,unsigned int tp_size>
unsigned int SNgeneric<T,tp_size>::getSize() const
{
    return tp_size;
}

// GET AND AT METHODS ------------------------------

template <class T,unsigned int tp_size>
T SNgeneric<T,tp_size>::get(const m_num& i,const m_num& j) const
{
    checkRangeCorectness(i,j);
    return _get(i,j);
}

template <class T,unsigned int tp_size>
T& SNgeneric<T,tp_size>::at(const m_num i,const m_num j)
{
    checkRangeCorectness(i,j);
    return _at(i,j);
}


// MATHEMATICAL FUNCTIONALITIES ------------------------


template <class T,unsigned int tp_size>
SNgaussian<T,tp_size> SNgeneric<T,tp_size>::getGaussian(const m_num c) const
{
    return SNgaussian<T,tp_size>(*this,c);
}


template <class T,unsigned int tp_size>
template <class V,unsigned int s>
void SNgeneric<T,tp_size>::subtract(const SNgeneric<V,s>& S)
{
    checkSizeCompatibility(*this,S);
    for (m_num i=0;i<tp_size;i++)
    {
        for (m_num j=0;j<tp_size;j++)
        {
            this->at(i,j)-=S.get(i,j);
        }
    }
}

template <class T,unsigned int tp_size>
template <class V,unsigned int s>
void SNgeneric<T,tp_size>::subtract(const SNgaussian<V,s>& G)
{
    checkSizeCompatibility(*this,G);
    m_num c=G.column;

    // subtract the non trivial "half column"
    for (m_num i=c+1;i<tp_size;i++)
    {
        this->at(i,c)-=G.get(i,c);
    }
    // subtract the diagonal
    for (m_num i=0;i<tp_size;i++)
    {
        this->at(i,i)-=1;
    }
}

template <class T,unsigned int tp_size>
template <class V,unsigned int s>
bool SNgeneric<T,tp_size>::isNumericallyEqual(const SNgeneric<V,s>& A,const double& epsilon) const
{
    checkSizeCompatibility(*this,A);
    T abs_diff;
    for (m_num i=0;i<tp_size;i++)
    {
        for (m_num j=0;j<tp_size;j++)
        {
            abs_diff=std::abs(  this->get(i,j)-A.get(i,j)  );
            if (abs_diff>epsilon)
            {
                return false;
            }
        }
    }
    return true;
}

#endif
