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

#include "../Utilities.h"

template <unsigned int tp_size>
class Mpermutation;

/**
    This class represents an elementary permutation, that is a 
    permutation that permutes only two elements.
*/

template <unsigned int tp_size>
class MelementaryPermutation
{
    template <unsigned int s>
    friend std::ostream& operator<<(std::ostream&, MelementaryPermutation<s>);

    private:
        /** The permutation permutes `data_elA` with `data_elB`  */
        const unsigned int data_elA; 
        /** The permutation permutes `data_elA` with `data_elB`  */
        const unsigned int data_elB;

    public :
        /** 
         * \param A an unsigned integer between `0` and `tp_size`
         * \param B an unsigned integer between `0` and `tp_size`
         *
         * Construct the permutation of \f$ A \f$ and \f$ B \f$.
         * */
        MelementaryPermutation(const unsigned int A,const unsigned int B); 

         /**
         *\brief Return the first of the two elements that are permuted.
         * */
        unsigned int getA() const;
         /**
         *\brief Return the second of the two elements that are permuted.
         * */
        unsigned int getB() const;

        /** 
        return by value the image of 'k' by the permutation.
         */
        unsigned int operator()(const unsigned int k) const;
        /** 
        return by value the image of 'k' by the permutation.
         */
        unsigned int image(const unsigned int k) const;

        /** 
         * The product is the composition of permutations.
         *
         * `perm1 * perm2` is the composition of the elementary permutations `perm1` and `perm2`.
         * The result is not an *elementary* permutation anymore.
         *
         * \see `Mpermutation`.
         * */
        Mpermutation<tp_size> operator*(const Mpermutation<tp_size>& ) const;
};


// GETTER/SETTER METHODS ---------------------------------

template <unsigned int tp_size>
unsigned int MelementaryPermutation<tp_size>::getA() const
{
    return data_elA;
}

template <unsigned int tp_size>
unsigned int MelementaryPermutation<tp_size>::getB() const
{
    return data_elB;
}

// OPERATORS ---------------------------------

template <unsigned int tp_size>
Mpermutation<tp_size> MelementaryPermutation<tp_size>::operator*(const Mpermutation<tp_size>& perm) const
{
    Mpermutation<tp_size> tmp;
    for (unsigned int k=0;k<tp_size;++k)
    {
        tmp.at(k)=this->image( perm.image(k)  );
    }
    return tmp;
}

template <unsigned int s>
std::ostream& operator<<(std::ostream& stream, MelementaryPermutation<s> perm)
{
    stream<<perm.getA()<<" <--> "<<perm.getB();
    return stream;
}

// ACTION ON THE INTEGERS ---------------------------------

template <unsigned int tp_size>
MelementaryPermutation<tp_size>::MelementaryPermutation(const unsigned int A,const unsigned int B) :
    data_elA(A),
    data_elB(B)
{
    if (getA() > tp_size or getB() > tp_size)
    {
        throw OutOfRangeConstructionElementaryPermutationException(A,B,tp_size);
    }
}
    
template <unsigned int tp_size>
unsigned int MelementaryPermutation<tp_size>::image(const unsigned int k) const
{
    if (k>tp_size)
    {
        throw PermutationIdexoutOfRangeException(k,tp_size);
    }
    if (k==getA())
    {
        return getB();
    }
    if (k==getB())
    {
        return getA();
    }
    return k;
}

template <unsigned int tp_size>
unsigned int MelementaryPermutation<tp_size>::operator()(const unsigned int k) const
{
    return image(k);
}


#endif


