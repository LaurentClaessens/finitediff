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
        SNpermutation(const Mpermutation<tp_size>& d); 

        /** The constructor with zero arguments initializes to the identity. */
        SNpermutation(); 

        Mpermutation<tp_size> getMpermutation() const;

        /** returns the inverse matrix */
        SNpermutation<T,tp_size> inverse() const;
};

// CONSTRUCTORS -------------------------------


template <class T,unsigned int tp_size>
SNpermutation<T,tp_size>::SNpermutation(const Mpermutation<tp_size>& d) :
    data(d)
{}

template <class T,unsigned int tp_size>
SNpermutation<T,tp_size>::SNpermutation() {}


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
    throw SNchangeNotAllowedException(i,j);
};

template <class T,unsigned int tp_size>
T SNpermutation<T,tp_size>::_get(const m_num i,const m_num j) const 
{
    if (getMpermutation().image(j)==i)
    {
        return 1;
    }
    return 0;
};

template <class T,unsigned int tp_size>
Mpermutation<tp_size> SNpermutation<T,tp_size>::getMpermutation() const
{
    return data;
}

// MATHEMATICS -------------------------------

template <class T,unsigned int tp_size>
SNpermutation<T,tp_size> SNpermutation<T,tp_size>::inverse() const
{
    SNpermutation<T,tp_size> inv(this->getMpermutation().inverse());
    return inv;
}

#endif
