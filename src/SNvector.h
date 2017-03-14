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

#ifndef __SNVECTOR_H__170345
#define __SNVECTOR_H__170345_

#include <array>

/*
This is my vector type, designed for numerical computation. 
*/

// THE CLASS HEADER -----------------------------------------

template <class T,unsigned int tp_size>
class SNvector
{

    private:
        std::array<T,tp_size> data;
    public :
        T get(unsigned int) const;
        T& at(unsigned int);
}


template <class T,unsigned int tp_size>
T SNvector<>T,tp_size>::get(unsigned int i) const
{
    return data.at(i);
}
template <class T,unsigned int tp_size>
T& SNvector<>T,tp_size>::at(unsigned int i) const
{
    return data.at(i);
}


#endif
