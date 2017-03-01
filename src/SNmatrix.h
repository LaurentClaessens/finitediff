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

#ifndef __SNMATRIX_H__153113__
#define __SNMATRIX_H__153113__

#include <array>

/*
This is my matrix type, designed for numerical computation.
*/

template <class T,int tp_size>
class SNmatrix
{
    private:
        std::array<T,tp_size> data;
        int size=tp_size;
    public:
        SNmatrix();
        int getSize() const;
};

template <class T,int tp_size>
SNmatrix<T,tp_size>::SNmatrix(): data(){};

template <class T,int tp_size>
int SNmatrix<T,tp_size>::getSize() const
{
    return size;
};

#endif
