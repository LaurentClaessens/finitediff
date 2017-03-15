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
   - override  T& at  
       this one has to finish with 'return get(i,j);'
       So it will throw the 'SNoutOfRangeException' if needed and call T _get


*/

#ifndef __SNGENERIC_H__142708_
#define __SNGENERIC_H__142708_

// THE CLASS HEADER -----------------------------------------

template <class T,unsigned int tp_size>
class SNgeneric
{
    private :
        virtual T _get(const unsigned int,const unsigned int) const=0;
    public:
        unsigned int getSize() const final;

        // 'at' is not final because some derived matrices have to
        // throw 'SNchangeNotAllowedException' for some elements.
        // Your implementation of 'at' should finish with 
        // return get(i,j);
        T& at(unsigned int,unsigned int);

        // Do not override 'get'.
        // What you want to override is the private '_get' that has to return the
        // element. Here 'get' check for the matrix size (throws if needed) and then
        // return what '_get' provides.
        T get(const unsigned int&,const unsigned int&) const final;
};

template <class T,unsigned int tp_size>
T SNgeneric<T,tp_size>::get(const unsigned int& i,const unsigned int& j) const
{
    if (i>tp_size or j>tp_size)
    {
        throw SNoutOfRangeException(i,j,getSize());
    }
    return _get(i,j);
}

template <class T,unsigned int tp_size>
unsigned int SNgeneric<T,tp_size>::getSize() const
{
    return tp_size;
}

template <class T,unsigned int tp_size>
T& SNgeneric<T,tp_size>::at(const unsigned int i,const unsigned int j)
{
    return get(i,j);
}

#endif
