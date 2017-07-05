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

/** 
 * This exception is raised when one tries to multiply or compare matrices that have not the same size.
 * */
class IncompatibleMatrixSizeException : public std::exception
{
    private :
        unsigned int size1;
        unsigned int size2;

        std::string _msg;

        std::string message(const unsigned int size1, const unsigned int size2) const
        {
            std::string s_s1=std::to_string(size1);
            std::string s_s2=std::to_string(size2);

            return "First matrix has size "+s_s1+" while second matrix has size "+s_s2;
        };

    public: 
        IncompatibleMatrixSizeException(const unsigned int s1, const unsigned int s2): 
            size1(s1),
            size2(s2),
            _msg( message(size1,size2) )
        {}
        virtual const char* what() const throw()
        {
            return _msg.c_str();
        }
};

class SNchangeNotAllowedException : public std::exception
{
    private :
        std::string _msg;

        std::string message(const unsigned int i,const unsigned int j,const std::string& text) const
        {
            std::string s_line=std::to_string(i);
            std::string s_col=std::to_string(j);

            return "You cannot modify element ("+s_line+" , "+s_col+" ) with this kind of matrix. You should maybe use 'get' instead of 'at'."+" "+text;
        };

    public: 
        SNchangeNotAllowedException(const unsigned int line, const unsigned int column,const std::string& text=""): 
            _msg( message(line,column,text) )
        {}
        virtual const char* what() const throw()
        {
            return _msg.c_str();
        }
};

class NotInitializedMemberException : public std::exception
{
    private :
        std::string text;

    public: 
        explicit NotInitializedMemberException(std::string t): 
            text(t)
    {}
        virtual const char* what() const throw()
        {
            return text.c_str();
        }
};

class OutOfRangeColumnNumber : public std::exception
{
    private :
        std::string text;

    public: 
        explicit OutOfRangeColumnNumber(std::string t): 
            text(t)
    {}
        virtual const char* what() const throw()
        {
            return text.c_str();
        }
};


/**
* When one ask by `get` or `at` an element that is out
* of the range of the matrix. 
*
* ```
* SNmatrix<int,4> A;
* A.get(5,1);     // raises
* ```
**/
class SNoutOfRangeException : public std::exception
{
    private :
        std::string _msg;
        std::string message(const unsigned int i,const unsigned int j,const unsigned int size) const 
        {
            std::string s_line=std::to_string(i);
            std::string s_col=std::to_string(j);
            std::string s_size=std::to_string(size);

            return "Attempt to access element ("+s_line+" , "+s_col+" ) while the matrix has size "+s_size;
        };

    public: 
        SNoutOfRangeException(const unsigned int i, const unsigned int j,const unsigned int s): 
            _msg(message(i,j,s))
    {}
        virtual const char* what() const throw()
        {
            return _msg.c_str();
        }
};

/** This exception is raised when one try to access to an element out of range of a permutation
*  ```
*  Mpermutation<4> perm;
*  // populate 'perm'
*  unsigned int k=perm(5)   // throws
*  ```
*
*  The constructor of `Mpermutation` from   `std::array<unsigned int,tp_size>`
*  also throws when one element of the array is larger than tp_size.
* */
class PermutationIdexoutOfRangeException : public std::exception
{

    private :
        std::string _msg;

        std::string message(const unsigned int index, const unsigned int size) const
        {
            std::string s_index=std::to_string(index);
            std::string s_size=std::to_string(size);

            return "Attempt to access element ("+s_index+" while the I am a permutation of integers from 0 to "+s_size;
        };

    public: 
        PermutationIdexoutOfRangeException(const unsigned int index, const unsigned int size): 
            _msg(message(index,size))
    {}
        virtual const char* what() const throw()
        {
            return _msg.c_str();
        }
};

class OutOfRangeConstructionElementaryPermutationException : public std::exception
{
    private :
        std::string _msg;
        std::string message(const unsigned int elA,const unsigned int elB,const unsigned int size) const
        {
            std::string s_elA=std::to_string(elA);
            std::string s_elB=std::to_string(elB);
            std::string s_size=std::to_string(size);

            return "Cannot create the elementary permutation "+s_elA+" <--> "+s_elB+" when tp_size is "+s_size+".";
        }

    public: 
        OutOfRangeConstructionElementaryPermutationException(const unsigned int A, const unsigned int B,const unsigned int s): 
            _msg(message(A,B,s))
        {}
        virtual const char* what() const throw()
        {
            return _msg.c_str();
        }
};

class ProbablyNotWhatYouWantException : public std::exception
{
    private :
        std::string message;

    public: 
        explicit ProbablyNotWhatYouWantException(const std::string m): 
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
        std::string _msg;

        std::string message(const int n)
        {
            std::string s_num=std::to_string(n);
            return "Trying to access line or column with negative number : "+s_num;
        };

    public: 
        explicit NegativeMatrixElementNumberException(const int n): 
            _msg(message(n))
    {}
        virtual const char* what() const throw()
        {
            return _msg.c_str();
        }
};

#endif
