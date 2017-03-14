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

#ifndef __SNMATRIX_H__153113__
#define __SNMATRIX_H__153113__

#include <array>
#include <iostream>
#include <cmath>

#include "SNelement.h"
#include "SNline.h"
#include "SNgaussianMatrix.h"
#include "../SNvector.h"

#include "../DebugPrint.h"
#include "../SNexceptions.cpp"
DebugPrint debug_print;

/*
This is my matrix type, designed for numerical computation. It represents a 
square matrix.
*/

// forward definition

template <class T,unsigned int tp_size>
class SNplu;

// THE CLASS HEADER -----------------------------------------

template <class T,unsigned int tp_size>
class SNmatrix
{
    friend class SNmatrixTest;
    friend class GaussTest;

    template <class V,unsigned int s>
    friend std::ostream& operator<<(std::ostream&, SNmatrix<V,s>&);

    template <class U,class V,unsigned int s>
    friend bool operator==(const SNmatrix<U,s>&,const SNmatrix<V,s>&);

    template <class U,class V,unsigned int s>
    friend bool operator==(const SNgaussianMatrix<U,s>&,const SNmatrix<V,s>&);

    private:
        std::array<T,tp_size*tp_size> data;
        unsigned int size=tp_size;
        // the larger element on column 'col' under (or on) the line 'f_line'. 
        SNelement<T,tp_size> getLargerUnder(unsigned int f_line, unsigned int col);

        // Substrat the given vector (line) from the line 'line'
        // Using the Gauss's elimination one need to do many differences like
        //  L_i -> L_i - k*L_1/m
        //  where 'k' is the first element of L_i and 'm' is the pivot.
        //  The function 'lineMinusLine' serves to not compute L_1/m
        //  as many times as the number of substitutions to do.
        void lineMinusLine(unsigned int line,SNline<T,tp_size> v);


        // return the larger element (in absolute value) on the given column
        // In case of equality, return the last one (the larger line).
        //   The template type T has to accept arithmetic manipulations
        //   like abs, comparison.
        SNelement<T,tp_size> getLargerOnColumn(unsigned int col);
        // return the largest (absolute value) element under the diagonal
        // on the given column.
        SNelement<T,tp_size> getLargerUnderDiagonal(unsigned int col);

        // From the line number "line", return a line normalized
        // in such a way that the first (non zero) element is 1.
        SNline<T,tp_size> gaussEliminationLine(unsigned int line);

    public:
        SNmatrix();
        SNmatrix(const SNmatrix<T,tp_size>&);


        unsigned int getSize() const;

        // return a reference to the value of the requested matrix entry.
        T& at(unsigned int,unsigned int);
        T get(unsigned int,unsigned int) const;

        // return the max of the absolute values of all the matrix elements
        T max_norm() const;

        // return the matrix element on given (line,column).
        SNelement<T,tp_size> getElement(unsigned int line, unsigned int col);

        // return the requested line
        SNline<T,tp_size> getSNline(unsigned int l);

        // Use the Gauss'elimination to transform the SNmatrix
        // to an upper triangular matrix.
        // In-place transformation !!
        void makeUpperTriangular();

        // swap the lines l1 and l2. This is in-place replacement.
        // The matrix is changed.
        void swapLines(unsigned int l1,unsigned int l2);

        // subtract the given matrix to 'this'.
        // this is a in-place replacement.
        template <class V,unsigned int s>
        void subtract(const SNmatrix<V,s>& B);

        // return the PLU decomposition.
        // This is a heavily non-const method because 'this' is transformed
        // into the U matrix.
        SNplu<T,tp_size> getPLU();

        // return the gaussian matrix for the requested column
        // we return a reference in order to avoid a copy.
        SNgaussianMatrix<T,tp_size>& getGaussian(unsigned int c) const;
};

// CONSTRUCTORS, OPERATORS, ...  -------------------------------------------

template <class T,unsigned int tp_size>
SNmatrix<T,tp_size>::SNmatrix(): data() { };

template <class T,unsigned int tp_size>
SNmatrix<T,tp_size>::SNmatrix(const SNmatrix<T,tp_size>& snm) : data(snm.data)  {};

template <class U,class V,unsigned int s>
bool operator==(const SNmatrix<U,s>& A,const SNmatrix<V,s>& B)
{
    return A.data==B.data;
}

template <class U,class V,unsigned int s>
bool operator==(const SNgaussianMatrix<U,s>& G,const SNmatrix<V,s>& A)
{
    for (unsigned int l=0;l<s;l++)
    {
        for (unsigned int c=0;c<s;c++)
        {
            if (G.get(l,c)!=A.get(l,c))
            {
                return false;
            }
        }
    }
    return true;
}

template <class V,unsigned int s>
std::ostream& operator<<(std::ostream& stream, SNmatrix<V,s>& snm)
{
    for (unsigned int l=0;l<s;l++)
    {
        stream<<snm.getSNline(l)<<std::endl;
    }
    return stream;
}

template <class T,unsigned int t> template<class V,unsigned int s>
void SNmatrix<T,t>::subtract(const SNmatrix<V,s>& B)
{
    for (unsigned int i=0;i<data.size();i++)
    {
        data.at(i)-=B.data.at(i);
    }
}

// GETTER METHODS  -------------------------------------------

template <class T,unsigned int tp_size>
unsigned int SNmatrix<T,tp_size>::getSize() const
{
    return size;
};

template <class T,unsigned int tp_size>
SNelement<T,tp_size> SNmatrix<T,tp_size>::getElement(unsigned int line, unsigned int col)
{
    return SNelement<T,tp_size>(line,col,*this);
}

template <class T,unsigned int tp_size>
SNline<T,tp_size> SNmatrix<T,tp_size>::getSNline(unsigned int l)
{
    std::array<T,tp_size> al;
    for (unsigned int c=0;c<tp_size;c++)
    {
        al.at(c)=get(l,c);
    }
    return SNline<T,tp_size>(al);
}

template <class T,unsigned int tp_size>
T& SNmatrix<T,tp_size>::at(const unsigned int i,const unsigned int j)
{
    if (i>tp_size or j>tp_size)
    {
        throw SNoutOfRangeException(i,j,tp_size);
    }
    return data.at(j*tp_size+i);
};

template <class T,unsigned int tp_size>
T SNmatrix<T,tp_size>::get(const unsigned int i,const unsigned int j) const
{
    if (i>tp_size or j>tp_size)
    {
        throw SNoutOfRangeException(i,j,tp_size);
    }
    return data.at(j*tp_size+i);
};


// GAUSS'S ELIMINATION METHODS


template <class T,unsigned int tp_size>
SNgaussianMatrix<T,tp_size>& SNmatrix<T,tp_size>::getGaussian(unsigned int c) const
{
    return SNgaussianMatrix<T,tp_size>(*this,c);
}

template <class T,unsigned int tp_size>
T SNmatrix<T,tp_size>::max_norm() const
{
    T m(0);
    for (T v:data)
    {
        T s=std::abs(v);
        if (s>m)
        {
            m=s;
        }
    }
    return m;
}

template <class T,unsigned int tp_size>
SNelement<T,tp_size> SNmatrix<T,tp_size>::getLargerUnder(unsigned int f_line, unsigned int col) 
{
    T max_val=0;
    unsigned int max_line=0;

    for (unsigned int line=f_line;line<tp_size;line++)
    {
        if (std::abs(get(line,col))>max_val)
        {
            max_val=std::abs(get(line,col));
            max_line=line;
        };
    };
    return getElement(max_line,col);
}

template <class T,unsigned int tp_size>
SNelement<T,tp_size> SNmatrix<T,tp_size>::getLargerOnColumn(unsigned int col) 
{
    return getLargerUnder(0,col);
}

template <class T,unsigned int tp_size>
SNelement<T,tp_size> SNmatrix<T,tp_size>::getLargerUnderDiagonal(unsigned int col) 
{
    return getLargerUnder(col,col);
}

template <class T,unsigned int tp_size>
void SNmatrix<T,tp_size>::swapLines(unsigned int l1, unsigned int l2)
{
    if (l1!=l2)
    {
        for (unsigned int col=0;col<tp_size;col++)
        {
            T tmp = get(l1,col);
            at(l1,col)=get(l2,col);
            at(l2,col)=tmp;
        }
    }
}


template <class T,unsigned int tp_size>
void SNmatrix<T,tp_size>::lineMinusLine(unsigned int line,SNline<T,tp_size> v)
{
    for (unsigned int c=0;c<tp_size;c++)
    {
        at(line,c)=get(line,c)-v.get(c);
    }
}

template <class T,unsigned int tp_size>
SNline<T,tp_size> SNmatrix<T,tp_size>::gaussEliminationLine(unsigned int line)
{
    SNline<T,tp_size> l=getSNline(line);
    l.makeUnit();
    return l;
}

template <class T,unsigned int tp_size>
SNplu<T,tp_size> SNmatrix<T,tp_size>::getPLU()

    // for each column :
    // - get the larger entry under the diagonal
    // - swap the 'larger' line with the current line
    // - the 'killing line' is the swapped line multiplied by the right number
    //   in such a way that the 'diagonal' entry is 1.
    // - use that 'killing line' to eliminate the column (under the diagonal)
    //
    // All the mathematics is explained with some details here :
    // http://laurent.claessens-donadello.eu/pdf/lefrido.pdf

{
    SNplu<T,tp_size> plu(*this);
    SNmatrix<T,tp_size> L;
    
    for (unsigned int c=0;c<tp_size;c++)
    {
        auto max_el = getLargerUnderDiagonal(c);

        if (max_el.getValue()!=0)   // not a column full of zero's
        {

            plu.permutations.at(c)=max_el.line;
            swapLines(c,max_el.line);
            auto killing_line=gaussEliminationLine(c);



            for (unsigned int l=c+1;l<tp_size;l++)
            {
                T m = get(l,c);  // the value to be eliminated

                // TODO : this is not optimal because
                // we already know the first 'c' differences are 0.
                lineMinusLine(l,m*killing_line);
            }
        }
        else
        {
            // if no permutations is done we record the trivial one.
            plu.permutations.at(c)=c;  
        }
    }
    return plu;
}

#endif
