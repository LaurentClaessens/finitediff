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

class IncompatibleMatrixSizeException : public std::exception
{
    private :
        unsigned int size1;
        unsigned int size2;

    public: 
        SNchangeNotAllowedException(const unsigned int s1, const unsigned int s2): 
            size1(s1),
            size2(s2)
    {}
        virtual const char* what() const throw()
        {
            std::string s_line=std::to_string(line);
            std::string s_col=std::to_string(column);

            std::string text= "You cannot modify element ("+s_line+" , "+s_col+" ) with this kind of matrix. You should maybe use 'get' instead of 'at'.";
            return text.c_str();
        }
};


class SNchangeNotAllowedException : public std::exception
{
    private :
        unsigned int line;
        unsigned int column;

    public: 
        SNchangeNotAllowedException(const unsigned int i, const unsigned int j): 
            line(i),
            column(j)
    {}
        virtual const char* what() const throw()
        {
            std::string s_line=std::to_string(line);
            std::string s_col=std::to_string(column);

            std::string text= "You cannot modify element ("+s_line+" , "+s_col+" ) with this kind of matrix. You should maybe use 'get' instead of 'at'.";
            return text.c_str();
        }
};


class SNoutOfRangeException : public std::exception
{
    private :
        unsigned int line;
        unsigned int column;
        unsigned int tp_size;

    public: 
        SNoutOfRangeException(const unsigned int i, const unsigned int j,const unsigned int s): 
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
