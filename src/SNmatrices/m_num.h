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


#ifndef __MNUM_H__094427__
#define __MNUM_H__094427__

#import <iostream>

/**
    This class is a wrapper for (a priori) `unsigned int`.

    It represents a number of line or column in a matrix.

    For the moment, the template parameter for the matrix size itself
    remains 'unsigned int'.
*/

class m_num
{
    private :
        unsigned int num;
    public :
        //cppcheck-suppress noExplicitConstructor
        m_num(const unsigned int n);  
        explicit m_num(const int n); 

        m_num operator++();  // ++i
        m_num operator++(int);  // i++

        bool operator >(const unsigned int& b) const;
        bool operator >(const m_num& b) const;
        bool operator >(const int& b) const;
        
        bool operator <(const unsigned int& b) const;
        bool operator <(const m_num& b) const;
        bool operator <(const int& b) const;

        /** Allows conversion to `unsigned int` */
        operator unsigned int() const;

        void swap(m_num& other);
};


#endif
