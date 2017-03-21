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


/*

    This class represents an elementary permutation, that is a 
    permutation that permutes only two elements.

*/

template <unsigned int tp_size>
class MelementaryPermutation
{
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

template <unsigned int tp_size>
Mpermutation<tp_size> MelementaryPermutation<tp_size>::operator*(const Mpermutation<tp_size>& perm) const
{
    Mpermutation<tp_size> tmp(perm);
    for (unsigned int k=0;k<tp_size;++k)
    {
        tmp.at(k)=this->get( perm.get(k)  );
    }
    return tmp;
}

template <unsigned int tp_size>
MelementaryPermutation<tp_size>::MelementaryPermutation(unsigned int A,unsigned int B) :
    elA(A),
    elB(B)
{}
    
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


