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


#ifndef __SNLINE_H__64621__
#define __SNLINE_H__64621__


/*
This class describes a matrix line from a matrix of type 'SNmatrix'.

An element contains
- its line number
- a pointer to the matrix which it belongs to.
  The latter it only accessible to the user trough the method 'getSNmatrix' that returns a reference
*/


// FORWARD DECLARATION ------------------------ 

template <class T,unsigned int tp_size>
class SNmatrix;

// THE CLASS HEADER -----------------------------------------

template <class T,unsigned int tp_size>
class SNline
{
    private :
        std::array<T,tp_size> data;
        unsigned int line;

        // Safe to use a raw pointer because I do not request any ownership.
        // The matrix is owned by someone else and will be released by that guy.
        SNmatrix<T,tp_size>* snmatrix;


    public :

        SNline(unsigned int line,SNmatrix<T,tp_size>& snmatrix);
        SNline(const std::array<T,tp_size>&);       // do not use this constructors otherwise than for testing purpose.

        template <class U,class V,unsigned int s>
        friend bool operator==(const SNline<U,s>&,const SNline<V,s>&);
        template <class V,unsigned int s>
        friend std::ostream& operator<<(std::ostream&, const SNline<V,s>&);

        // return the value of the ith element on the line.
        //  - at returns a reference (! non const)
        //  - get returns a value.
        T& at(const unsigned int i);
        T get(const unsigned int i) const;

        // return the matrix from which the element is extracted
        SNmatrix<T,tp_size>& getSNmatrix() const;

        // return the number of the first non-zero element in the line.
        unsigned int firstNonZeroColumn() const;

        // divide all the elements by the first non-zero element. So the line becomes
        // (0,0,1, blah )
        // the original matrix is not updated.
        void makeUnit();
};

// CONSTRUCTOR, ASSIGNATION, ...  -------------------------------------------

template <class T,unsigned int tp_size>
SNline<T,tp_size>::SNline(unsigned int l,SNmatrix<T,tp_size>& snm) : 
    line(l),
    snmatrix(&snm)
{ }

// the following constructor does not initialize the referenced matrix.
template <class T,unsigned int tp_size>
SNline<T,tp_size>::SNline(const std::array<T,tp_size>& ar) : data(ar) {};


// OPERATORS -------------------------------------------

template <class U,class V,unsigned int s>
bool operator==(const SNline<U,s>& A,const SNline<V,s>& B)
{
    return A.data==B.data;
}

template <class V,unsigned int s>
std::ostream& operator<<(std::ostream& stream, const SNline<V,s>& line)
{
    for (unsigned int c=0;c<s;c++)
    {
        stream<<line.get(c)<<",";
    }
    return stream;
}

// GETTER METHODS -------------------------------------------

template <class T,unsigned int tp_size>
T& SNline<T,tp_size>::at(const unsigned int i)
{
    if (i>tp_size)
    {
        std::cout<<"This SNline has size "<<tp_size<<". Attempt to access element "<<i<<" ."<<std::endl;
    }
    return data.at(i);
};

template <class T,unsigned int tp_size>
T SNline<T,tp_size>::get(const unsigned int i) const
{
    if (i>tp_size)
    {
        std::cout<<"This SNline has size "<<tp_size<<". Attempt to access element "<<i<<" ."<<std::endl;
    }
    return data.at(i);
};

template <class T,unsigned int tp_size>
SNmatrix<T,tp_size>& SNline<T,tp_size>::getSNmatrix() const
{
    return *snmatrix;
}

// OTHER FUNCTIONALITIES -------------------------------------------


template <class T,unsigned int tp_size>
unsigned int SNline<T,tp_size>::firstNonZeroColumn() const
    // return the number of the first non-zero element.
    // If they are all non-zero, return the size+1.
{
    for (unsigned int col=0;col<tp_size;col++)
    {
        if (get(col)!=0)
        {
            return col;
        }
    }
   return tp_size+1;
}

template <class T,unsigned int tp_size>
void SNline<T,tp_size>::makeUnit()
{
    unsigned int col=firstNonZeroColumn();
    if (col!=tp_size+1)
    {
        const T m = get(col);
        at(col)=1;      // the first one is by hand 1 (because we know it).
        for (unsigned int c=col+1;c<tp_size;c++)
        {
            at(c)=get(c)/m;
        }
    }
}

#endif
