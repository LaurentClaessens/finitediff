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
- a reference to the matrix which it belongs to.


*/

template <class T,int tp_size>
class SNelement
{
    private :
        int line;
        int column;
        SNmatrix<T,tp_size>& snmatrix;
    public :
        SNelement(int line,int column,T value,const SNmatrix<T,tp_size>& snmatrix);
        T getValue();
};

SNelement::SNelement(int l,int c,T val,const SNmatrix<T,tp_size>& snm) : 
    line(l),
    column(c),
    smnatrix(snm)
{}

#endif
