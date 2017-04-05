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

#include "m_num.h"
#include "../SNexceptions.cpp"

// CONSTRUCTOR --------------------------------


m_num::m_num(const unsigned int n) : 
    num(n)
{}

m_num::m_num(const int n) 
{
    if (n<0)
    {
        throw NegativeMatrixElementNumberException(n);
    }
    num=n;
}

// CONVERSIONS   ----------------------------------

m_num::operator unsigned int() const
{
    return num;
}

// INCREMENT  ----------------------------------

m_num m_num::operator++() 
{
    ++num;
    return *this;
}
m_num m_num::operator++(int) 
{
    m_num tmp(*this);
    ++num;
    return tmp;
}

// COMPARISON -------------------------- 

bool m_num::operator >(const unsigned int& b) const { return num>b; }
bool m_num::operator >(const m_num& b) const { return num>b.num; }
bool m_num::operator >(const int& b) const 
{ 
    return int(num)>b; 
}
bool m_num::operator <(const unsigned int& b) const 
{
    return num<b;
}
bool m_num::operator <(const m_num& b) const { return num<b.num; }
bool m_num::operator <(const int& b) const 
{ 
    return int(num)<b; 
}
