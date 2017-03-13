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
*/

// THE CLASS HEADER -----------------------------------------

template <unsigned int tp_size>
class SNpermutation
{
    template <unsigned int s>
    friend std::ostream& operator<<(std::ostream&, SNpermutation<s>&);
    
    private:
        std::array<tp_size> permutations;
}


// OPERATORS -------------------------------

template <unsigned int s>
std::ostream& operator<<(std::ostream& stream, SNpermutation<s>& perm)
{
    for (unsigned int l=0;l<s;l++)
    {
        stream<<perm.at(l)<<std::endl;
    }
    return stream;
}

