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

#ifndef __MGENERICPERMUTATION_H__238823794
#define __MGENERICPERMUTATION_H__238823794



//forward
template <unsigned int>
class MelementaryPermutation;

/** 
 * \brief A permutation is a bijection of a finite subset of N
 *
 * A permutation is a bijection of a set of the form \f$ \{0,1,...,n\} \f$
 * for some \f$ n \f$.
 *
 * The template parameter `<unsigned int tp_size>` is the \f$ n \f$.
 *
 * This virtual base class provides the common operations of
 * `MelementaryPermutation` and `Mpermutation`.
 *
 * - Mainly the product `operator*` have the same implementation for the 4
 * possibilities of products `MelementaryPermutation` by `Mpermutation`.
 *
 * - Also rationalize the definition of `operator()` from the definition of
 *   the member function `image`.
 */


template <unsigned int tp_size>
class MgenericPermutation
{
    public:

        // For computing the image :
        // - each subclass define its own `image` member function.
        //
        // - The operator () is defined here
        //
        // - The member function `at` has not to be defined for each
        //   subclass because it makes no sense for an elementary permutation.
        //
        //   TODO : the latter is not really true. We could define `at` for
        //   everyone using the `SpecialValue` trick.

        /** return by value the image of `k` */
        virtual unsigned int operator()(const unsigned int k) const final;
        /** return by value the image of `k` */
        virtual unsigned int image(const unsigned int k) const=0;

};

template <unsigned int tp_size>
unsigned int MgenericPermutation<tp_size>::operator()(const unsigned int k) const
{
    return image(k);
}

#endif
