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


#ifndef __SNELEMENT_H__145452__
#define __SNELEMENT_H__145452__

/*
This class describes a matrix element from a matrix of type 'SNmatrix'.

An element contains
- line and column
- its value
*/

// THE CLASS HEADER -----------------------------------------


/** 
 * \brief This is a template class that represent an element of a `SNmatrix`.
 *
 * The template arguments `T` and `tp_size` have the same meaning. 
 * */
template <class T,unsigned int tp_size>
class SNelement
{
    private:
        const T value;
    public :
        SNelement(const unsigned int line,const unsigned int column,const T v);

        SNelement<T,tp_size> operator=(const SNelement<T,tp_size>& other);

        const unsigned int line;
        const unsigned int column;

        // return the value of the matrix element
        T getValue() const;
};

// CONSTRUCTOR, ASSIGNATION, ...  -------------------------------------------

template <class T,unsigned int tp_size>
SNelement<T,tp_size>::SNelement(const unsigned int l,const unsigned int c, const T v) : 
    value(v),
    line(l),
    column(c)
{}

// OTHER FUNCTIONALITIES -------------------------------------------

template <class T,unsigned int tp_size>
T SNelement<T,tp_size>::getValue() const
{
   return value;
}

#endif
