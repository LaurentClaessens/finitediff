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

// THE CLASS HEADER -----------------------------------------

template <class T,unsigned int tp_size>
class SNgeneric
{
    private :
        virtual T _get(const unsigned int,const unsigned int) const=0;
        virtual T _at(const unsigned int,const unsigned int)=0;

        // throws 'SNoutOfRangeException' if the requested element is out of
        // range (larger than 'tp_size').
        virtual void checkRangeCorectness(const unsigned int&,const unsigned int&) const final;
    public:
        virtual unsigned int getSize() const final;

        virtual T& at(unsigned int,unsigned int) final;
        virtual T get(const unsigned int&,const unsigned int&) const final;
};


// UTILITIES --------------------------


template <class T,unsigned int tp_size>
void SNgeneric<T,tp_size>::checkRangeCorectness(const unsigned int& l,const unsigned int& c) const 
{
    if (l>tp_size or c>tp_size)
    {
        throw SNoutOfRangeException(l,c,getSize());
    }
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

#endif
