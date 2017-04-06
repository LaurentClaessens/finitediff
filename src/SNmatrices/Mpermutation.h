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

#ifndef __MPERMUTATION_H_121119__
#define __MPERMUTATION_H_121119__

#include <array>

#include "MelementaryPermutation.h"

// THE CLASS HEADER -----------------------------------------

/**
 This class represents a permutation. (not a matrix)


 The permutation is recorded in the array
        `std::array<unsigned int,tp_size> permutations;`
 The array
 (a,b,c,d)
 represents the permutation that maps
 1 to a
 2 to b
 3 to c
 4 to d

 The numbers a,b,c,d must be different and in [0,tp_size].

!! THIS IS NOT THE USUAL CONVENTION IN MATH !!

In math, the usual convention when writing
(a,b,c,d)
is
a->b
b->c
c->d
d->a

*/

template <unsigned int tp_size>
class Mpermutation
{

    template <unsigned int s>
    friend std::ostream& operator<<(std::ostream&, Mpermutation<s>);
    
    private:
        std::array<unsigned int,tp_size> data;
    public :
        Mpermutation(const std::array<unsigned int,tp_size>& d); 
        Mpermutation(const MelementaryPermutation<tp_size>& ); 

        /** The no-argument constructors initializes to identity */
        Mpermutation(); 


        /** return by value the image of `k` */
        unsigned int operator()(const unsigned int k) const;
        /** return by value the image of `k` */
        unsigned int image(const unsigned int k) const;

        // return by reference the image of 'k' by the permutation
        // allows to populate.
        unsigned int& at(const unsigned int k);

        bool operator==(const Mpermutation<tp_size>& other) const;

        /** Return the inverse permutation */
        Mpermutation<tp_size> inverse() const;
};


// GETTER METHODS ----------------------------


template <unsigned int tp_size>
unsigned int& Mpermutation<tp_size>::at(const unsigned int k)
{
    return data.at(k);
}

// CONSTRUCTOR ----------------------------

template <unsigned int tp_size>
Mpermutation<tp_size>::Mpermutation(const std::array<unsigned int,tp_size>& d) :
    data(d)
{
    for (unsigned int k=0;k<tp_size;++k)
    {
        if (d.at(k)>tp_size-1) // Mpermutation<4> permutes the set {0,1,2,3}.
        {
            throw PermutationIdexoutOfRangeException(k,tp_size);
        }
    }
}

template <unsigned int tp_size>
Mpermutation<tp_size>::Mpermutation(const MelementaryPermutation<tp_size>& p )
{
    this->at(  p.getA()  )=p.getB();
    this->at(  p.getB()  )=p.getA();
}

template <unsigned int tp_size>
Mpermutation<tp_size>::Mpermutation() 
{
    for (unsigned int k=0;k<tp_size;++k)
    {
        at(k)=k;
    }
}

// OPERATORS -------------------------------

template <unsigned int tp_size>
bool Mpermutation<tp_size>::operator==(const Mpermutation<tp_size>& other) const
{
    return data==other.data;
}

template <unsigned int tp_size>
unsigned int Mpermutation<tp_size>::image(const unsigned int k) const
{
    if (k>tp_size)
    {
        throw PermutationIdexoutOfRangeException(k,tp_size);
    }
    return data.at(k);
}

template <unsigned int tp_size>
unsigned int Mpermutation<tp_size>::operator()(const unsigned int k) const
{
    return image(k);
}

template <unsigned int s>
std::ostream& operator<<(std::ostream& stream, Mpermutation<s> perm)
{
    for (unsigned int l=0;l<s;l++)
    {
        stream<<l<<"->"<<perm.data.at(l)<<std::endl;
    }
    return stream;
}


// MATHEMATICS -------------------------------


template <unsigned int tp_size>
Mpermutation<tp_size> Mpermutation<tp_size>::inverse () const
{
    Mpermutation<tp_size> inv;
    for (unsigned int k=0;k<tp_size;++k)
    {
        inv.at(  this->image(k)  )=k;
    }
    return inv;
}

#endif
