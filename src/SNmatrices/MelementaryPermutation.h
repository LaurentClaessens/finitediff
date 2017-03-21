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

#ifndef __MELEMENTARYPERMUTATION_H_095019_
#define __MELEMENTARYPERMUTATION_H_095019_

#include "../DebugPrint.h"

/*

    This class represents an elementary permutation, that is a 
    permutation that permutes only two elements.

*/

template <unsigned int tp_size>
class MelementaryPermutation
{
    template <unsigned int s>
    friend std::ostream& operator<<(std::ostream&, MelementaryPermutation<s>);

    private:
        unsigned int elA;
        unsigned int elB;
    public :
        MelementaryPermutation(unsigned int A,unsigned int B); 

        // return by value the image of 'k' by the permutation.
        // There are two ways to do that.
        unsigned int operator()(const unsigned int k) const;
        unsigned int get(const unsigned int k) const;

        Mpermutation<tp_size> operator*(const Mpermutation<tp_size>& ) const;
};


// OPERATORS ---------------------------------

template <unsigned int tp_size>
Mpermutation<tp_size> MelementaryPermutation<tp_size>::operator*(const Mpermutation<tp_size>& perm) const
{
    Mpermutation<tp_size> tmp;
    debug_print<<"Le this est : "<<std::endl;
    debug_print<<*this<<std::endl;

    debug_print<<"À composer avec "<<std::endl;
    debug_print<<perm<<std::endl;

    debug_print<<"Le tmp généré est : "<<std::endl<<tmp<<std::endl;
    for (unsigned int k=0;k<tp_size;++k)
    {
        debug_print<<"k est : "<<k<<std::endl;
        debug_print<<"perm.get(k) est : "<<perm.get(k)<<std::endl;
        tmp.at(k)=this->get( perm.get(k)  );
    }
    return tmp;
}

template <unsigned int s>
std::ostream& operator<<(std::ostream& stream, MelementaryPermutation<s> perm)
{
    stream<<perm.elA<<" <--> "<<perm.elB;
    return stream;
}

// ACTION ON THE INTEGERS ---------------------------------

template <unsigned int tp_size>
MelementaryPermutation<tp_size>::MelementaryPermutation(unsigned int A,unsigned int B) :
    elA(A),
    elB(B)
{
    if (elA > tp_size or elB>tp_size)
    {
        throw OutOfRangeConstructionElementaryPermutationException(A,B,tp_size);
    }
}
    
template <unsigned int tp_size>
unsigned int MelementaryPermutation<tp_size>::get(const unsigned int k) const
{
    if (k>tp_size)
    {
        throw PermutationIdexoutOfRangeException(k,tp_size);
    }
    if (k==elA)
    {
        return elB;
    }
    if (k==elB)
    {
        return elA;
    }
    return k;
}

template <unsigned int tp_size>
unsigned int MelementaryPermutation<tp_size>::operator()(const unsigned int k) const
{
    return get(k);
}


#endif


