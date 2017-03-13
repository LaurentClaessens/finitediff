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

#ifndef __DEBUGPRINT_H__085438__
#define __DEBUGPRINT_H__085438__


/*
   I often need to pack my code with 'std::cout<<' for debugging purposes.o

   This class serves to make
   'debug_print<<'
   and to remove them easier with a small grep/ack.
*/

#include <iostream>

class DebugPrint
{
    template <class V>
    friend DebugPrint& operator<<(DebugPrint&, V);

    public:
        void endl() const;

};

void DebugPrint::endl() const
{
    std::cout<<std::endl;
}

template <class V>
DebugPrint& operator<<(DebugPrint& stream, V s) 
{
    std::cout<<s;
    return stream;
}

#endif
