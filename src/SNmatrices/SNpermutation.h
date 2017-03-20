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

#ifndef __SNPERMUTATION_H_060431__
#define __SNPERMUTATION_H_060431__

#include <array>

#include "Mpermutation.h"

/*
 This class represents a permutation matrix. 

 The permutation is recorded in a `Mpermutation` object.

*/

// THE CLASS HEADER -----------------------------------------

template <class T,unsigned int tp_size>
class SNpermutation : public SNgeneric<T,tp_size>
{
    template <class V,unsigned int s>
    friend std::ostream& operator<<(std::ostream&, SNpermutation<V,s>);
    
    private:
        Mpermutation<tp_size> data;
        T& _at(const m_num i,const m_num j) override;
        T _get(const m_num i,const m_num j) const override;
    public :
        SNpermutation(Mpermutation<tp_size> d); 
};


// OPERATORS -------------------------------

template <class T,unsigned int s>
std::ostream& operator<<(std::ostream& stream, SNpermutation<T,s> perm)
{
    for (unsigned int l=0;l<s;l++)
    {
        stream<<perm.data.at(l);
    }
    return stream;
}

// _AT AND _GET -------------------------------

template <class T,unsigned int tp_size>
T& SNpermutation<T,tp_size>::_at(const m_num i,const m_num j)
{
    std::cout<<"Not yet implemented _at for "<<i<<","<<j<<std::endl;
    return data.at(0);
};

#endif
