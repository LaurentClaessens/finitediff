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

/*
 This class represents a permutation matrix. It could be a 'SNmatrix', but since it is
 for sure a sparse matrix, we use a simpler implementation.


 The permutation is recorded in the array
        std::array<unsigned int,tp_size> permutations;
 The array
 (1,4,4,5)
 represents the permutation
 (1,1)(2,4)(3,4)(4,5).
*/

// THE CLASS HEADER -----------------------------------------

template <unsigned int tp_size>
class SNpermutation
{
    template <unsigned int s>
    friend std::ostream& operator<<(std::ostream&, SNpermutation<s>);
    
    private:
        std::array<unsigned int,tp_size> data;

    public:
        unsigned int& at(const unsigned int i);
};


// OPERATORS -------------------------------

template <unsigned int s>
std::ostream& operator<<(std::ostream& stream, SNpermutation<s> perm)
{
    for (unsigned int l=0;l<s;l++)
    {
        stream<<perm.data.at(l)<<std::endl;
    }
    return stream;
}

template <unsigned int tp_size>
unsigned int& SNpermutation<tp_size>::at(const unsigned int i)
{
    if (i>tp_size)
    {
        std::cout<<"This SNpermutation has size "<<tp_size<<". Attempt to access element "<<i<<std::endl;
    }
    return data.at(i);
};

#endif
