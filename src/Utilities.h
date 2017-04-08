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
#include <iostream>

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
unsigned int value_length(const T& value)
{
    std::stringstream s;
    s<<value;
    return s.tellp();
}

/** The test functions use this function to announce (output) that they are 
 * starting */
void echo_function_test(const std::string& s);

/** Use this function to announce (output) that a `CPPUNIT_ASSERT`
 * is about to be launched.
 * */
void echo_single_test(const std::string& s);

/**
 So I can put `debug_print<<blahblah` everywhere in my code (for debugging purpose) and found them back with a simple `ack`.
 */ 
extern std::ostream& debug_print;

/** 
 * Print a small message and a matrix, for debug purpose.
 *
 *\param name A small descrition of the matrix that will be printed.
 *\param mte the matrix to print.
 * */
template <class M>
void debug_matrix_print(const std::string& name, const M& mtr)
{
    std::cout<<name<<std::endl;
    std::cout<<mtr<<std::endl;
}


/**
 * \brief Display a small warning.
 *
 * \param message The message to be displayed.
 *
 * We have several type of matrices and several overloading of the operators like
 * `==` and `*`. These overloads are optimized for the particular structure of
 * our matrices.
 *
 * However generic overloads are also implemented. We use `tooGenericWarning` to 
 * print a warning when a generic overload is used because it could be a sign that
 * we are lacking an overload (or that the matrices types are not the ones
 * we believe).
 *
 *
 * \see operator*(const SNgeneric<U,s>& A, const SNgeneric<V,t>& B)
 * */
void tooGenericWarning(const std::string& message);

#endif
