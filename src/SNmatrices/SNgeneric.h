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

/*
 
   This is the base class for the other matrices types.
   A matrix can be
   - diagonal
   - upper/lower triangular
   - a gaussian transformation
   - symmetric
   - etc.

   Each of these types have different implementations for storing the elements, 
   perform equality tests (opertor==), multiplication (operator*) and so on.

   Much of these operation, on the least ressort, depend only on the 'get' method
   that returns an element of the matrix.
 
 
   What you need in your subclass :

   - override  T  _get   (mandatory)
       This one has to finish with 
       return SNgeneric<T,tp_size>::get(i,j);
   - override  T& _at    (mandatory)
        This one has to finish with
       return SNgeneric<T,tp_size>::at(i,j);

    _get(i,j) has to return the value of the element (i,j) of your matrix.

    _at(i,j) has to return a _reference_ to the element (i,j) of your matrix.
            So 'at' is intended to modify the matrix.
            In consequence, _at(i,j) should throw 'SNchangeNotAllowedException' if
            the requested element cannot be modified. 
            
            
    As an example, on a lower diagonal matrix,
    _get(1,3) returns 0 (by value)
    _at(1,3) throws SNchangeNotAllowedException


*/

#ifndef __SNGENERIC_H__142708_
#define __SNGENERIC_H__142708_

#include "SNgaussianMatrix.h"
#include "SNline.h"

#include "MathUtilities.h"

// THE CLASS HEADER -----------------------------------------

template <class T,unsigned int tp_size>
class SNgeneric
{
    //template <class V,unsigned int s>
    //friend std::ostream& operator<<(std::ostream&,const SNgeneric<V,s>&);

    private :
        virtual T& _at(const unsigned int,const unsigned int)=0;
        virtual T _get(const unsigned int,const unsigned int) const=0;

        // throws 'SNoutOfRangeException' if the requested element is out of
        // range (larger than 'tp_size').
        virtual void checkRangeCorectness(const unsigned int&,const unsigned int&) const final;
    public:
        virtual unsigned int getSize() const final;

        virtual T& at(unsigned int,unsigned int) final;
        virtual T get(const unsigned int&,const unsigned int&) const final;

        virtual SNline<T,tp_size> getSNline(unsigned int l) const;

        // subtract the given matrix from this matrix.
        // This is in-place replacement. Thus the least const in the 
        // world.
        template <class V,unsigned int s>
        void subtract(const SNgeneric<V,s>&);
        template <class V,unsigned int s>
        void subtract(const SNgaussianMatrix<V,s>&);

        // return the gaussian matrix for the requested column
        SNgaussianMatrix<T,tp_size> getGaussian(const unsigned int c) const;

        // numerical equality test 'up to epsilon'
        // If if max norm of "this-A" is strictly larger than epsilon,
        // return false.
        // Else return true.
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
void SNgeneric<T,tp_size>::checkRangeCorectness(const unsigned int& l,const unsigned int& c) const 
{
    if (l>tp_size or c>tp_size)
    {
        throw SNoutOfRangeException(l,c,getSize());
    }
}

template <class T,unsigned int tp_size>
SNline<T,tp_size> SNgeneric<T,tp_size>::getSNline(unsigned int l) const
{
    std::array<T,tp_size> al;
    for (unsigned int c=0;c<tp_size;c++)
    {
        al.at(c)=this->get(l,c);
    }
    return SNline<T,tp_size>(al);
}

// OPERATORS ------------------------------

template <class V,unsigned int s>
std::ostream& operator<<(std::ostream& stream,const SNgeneric<V,s>& snm)
{
    for (unsigned int l=0;l<s;l++)
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
T SNgeneric<T,tp_size>::get(const unsigned int& i,const unsigned int& j) const
{
    checkRangeCorectness(i,j);
    return _get(i,j);
}

template <class T,unsigned int tp_size>
T& SNgeneric<T,tp_size>::at(const unsigned int i,const unsigned int j)
{
    checkRangeCorectness(i,j);
    return _at(i,j);
}


// MATHEMATICAL FUNCTIONALITIES ------------------------


template <class T,unsigned int tp_size>
SNgaussianMatrix<T,tp_size> SNgeneric<T,tp_size>::getGaussian(const unsigned int c) const
{
    return SNgaussianMatrix<T,tp_size>(*this,c);
}


template <class T,unsigned int tp_size>
template <class V,unsigned int s>
void SNgeneric<T,tp_size>::subtract(const SNgeneric<V,s>& S)
{
    checkSizeCompatibility(*this,S);
    for (unsigned int i=0;i<tp_size;i++)
    {
        for (unsigned int j=0;j<tp_size;j++)
        {
            this->at(i,j)-=S.get(i,j);
        }
    }
}

template <class T,unsigned int tp_size>
template <class V,unsigned int s>
void SNgeneric<T,tp_size>::subtract(const SNgaussianMatrix<V,s>& G)
{
    checkSizeCompatibility(*this,G);
    unsigned int c=G.column;

    // subtract the non trivial "half column"
    for (unsigned int i=c+1;i<tp_size;i++)
    {
        this->at(i,c)-=G.get(i,c);
    }
    // subtract the diagonal
    for (unsigned int i=0;i<tp_size;i++)
    {
        this->at(i,i)-=1;
    }
}

template <class T,unsigned int tp_size>
template <class V,unsigned int s>
bool SNgeneric<T,tp_size>::isNumericallyEqual(const SNgeneric<V,s>& A,const double& epsilon) const
{
    checkSizeCompatibility(*this,A);
    T diff;
    for (unsigned int i=0;i<tp_size;i++)
    {
        for (unsigned int j=0;j<tp_size;j++)
        {
            diff=std::abs(  this->get(i,j)-A.get(i,j)  );
            if (diff>epsilon)
            {
                return false;
            }
        }
    }
    return true;
}

#endif
