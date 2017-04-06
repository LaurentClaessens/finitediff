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

/** 
 * This file consists in a bunch of small functions.
*/


#ifndef __UTILITIES_H_062410__
#define __UTILITIES_H_062410__

#include <sstream>

/**
 * return the length of a value supposed to be a numerical type.
 *
 * Note :
 * ```
 * double a=4;
 * std::string s=std::to_string(a);
 * ```
 * in that case, `s` is something like "4.000000", and its size will be like 8 or 9.
 *
 * This function *does not* act like that and return the length of `value` as 
 * printed by `std::cout`.
 * */
template <class T>
unsigned int value_length(const T& value);


#endif
