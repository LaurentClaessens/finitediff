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

#ifndef __EXCEPTIONS_H__095622__
#define __EXCEPTIONS_H__095622__

#include <exception>
#include <string>

class SNOutOfRangeException : public std::exception
{
    private :
        unsigned int line;
        unsigned int column;
        unsigned int tp_size;

    public: 
        SNOutOfRangeException(const unsigned int i, const unsigned int j,const unsigned int s): 
            line(i),
            column(j),
            tp_size(s)
    {}
        virtual const char* what() const throw()
        {
            std::string s_line=std::to_string(line);
            std::string s_col=std::to_string(column);
            std::string s_size=std::to_string(tp_size);

            std::string text= "Attempt to access element ("+s_line+" , "+s_col+" ) while the matrix has size "+s_size;
            return text.c_str();
        }
};

#endif
