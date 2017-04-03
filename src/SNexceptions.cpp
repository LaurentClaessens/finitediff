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
        IncompatibleMatrixSizeException(const unsigned int s1, const unsigned int s2): 
            size1(s1),
            size2(s2)
    {}
        virtual const char* what() const throw()
        {
            std::string s_s1=std::to_string(size1);
            std::string s_s2=std::to_string(size2);

            std::string text= "First matrix has size "+s_s1+" while second matrix has size "+s_s2;
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

class PermutationIdexoutOfRangeException : public std::exception
{

    /** This exception is raised when one try to access to an element out of range of a permutation
     *  ```
     *  Mpermutation<4> perm;
     *  // populate 'perm'
     *  unsigned int k=perm(5)   // throws
     *  ```
     *
     *  The constructor of `Mpermutation` from   `std::array<unsigned int,tp_size>`
     *  also throws when one element of the array is large than tp_size.
     * */

    private :
        unsigned int index;
        unsigned int tp_size;

    public: 
        PermutationIdexoutOfRangeException(const unsigned int i, const unsigned int t): 
            index(i),
            tp_size(t)
    {}
        virtual const char* what() const throw()
        {
            std::string s_index=std::to_string(index);
            std::string s_size=std::to_string(tp_size);

            std::string text= "Attempt to access element ("+s_index+" while the I am a permutation of integers from 0 to "+s_size;
            return text.c_str();
        }
};

class OutOfRangeConstructionElementaryPermutationException : public std::exception
{
    private :
        unsigned int elA;
        unsigned int elB;
        unsigned int tp_size;

    public: 
        OutOfRangeConstructionElementaryPermutationException(const unsigned int A, const unsigned int B,const unsigned int s): 
            elA(A),
            elB(B),
            tp_size(s)
    {}
        virtual const char* what() const throw()
        {
            std::string s_elA=std::to_string(elA);
            std::string s_elB=std::to_string(elB);
            std::string s_size=std::to_string(tp_size);

            std::string text= "Cannot create the elementary permutation "+s_elA+" <--> "+s_elB+" when tp_size is "+s_size+".";
            return text.c_str();
        }
};

class ProbablyNotWhatYouWantException : public std::exception
{
    private :
        std::string message;

    public: 
        ProbablyNotWhatYouWantException(const std::string m): 
            message(m)
    {}
        virtual const char* what() const throw()
        {
            return message.c_str();
        }
};

class NegativeMatrixElementNumberException : public std::exception
{
    private :
        int num;

    public: 
        NegativeMatrixElementNumberException(const int n): 
            num(n)
    {}
        virtual const char* what() const throw()
        {
            std::string s_num=std::to_string(num);

            std::string text= "Trying to access line or column with negative number : "+s_num;
            return text.c_str();
        }
};


#endif
