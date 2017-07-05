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

#include "SNgeneric.h"
#include "SNelement.h"
#include "SNline.h"
#include "SNgaussian.h"
#include "SNupperTriangular.h"
#include "Mpermutation.h"
#include "SNpermutation.h"
#include "MelementaryPermutation.h"
#include "MathUtilities.h"
#include "SNoperators.h"
#include "../SNvector.h"
#include "../exceptions/SNexceptions.cpp"

#include "../Utilities.h"


// forward definition
template <class T,unsigned int tp_size>
class SNplu;


/*
This is my matrix type, designed for numerical computation. It represents a 
square matrix.

NOTE : if you want the identity matrix, there is the `SNidentity` class.
*/
template <class T,unsigned int tp_size>
class SNmatrix  : public SNgeneric<T,tp_size>
{

    friend class SNmatrixTest;
    friend class GaussTest;

    template <class U,unsigned int s,class V,unsigned int t>
    friend bool operator==(const SNmatrix<U,s>&,const SNmatrix<V,t>&);
    template <class U,unsigned int s,class V,unsigned int t>
    friend SNmatrix<U,s> operator+(const SNmatrix<U,s>& A,const SNmatrix<V,t>& B);
    
    friend std::array<T,tp_size*tp_size> SNupperTriangular<T,tp_size>::_get_other_data(const SNmatrix<T,tp_size>&) const;
    friend std::array<T,tp_size*tp_size> SNlowerTriangular<T,tp_size>::_get_other_data(const SNmatrix<T,tp_size>&) const;


    private:
        std::array<T,tp_size*tp_size> data;
        unsigned int size=tp_size;

        /**  the larger element on column 'col' under (or on) the line 'f_line'.*/
        SNelement<T,tp_size> getLargerUnder(m_num f_line, m_num col) const;

        // Substrat the given vector (line) from the line 'line'
        // Using the Gauss's elimination one need to do many differences like
        //  L_i -> L_i - k*L_1/m
        //  where 'k' is the first element of L_i and 'm' is the pivot.
        //  The function 'lineMinusLine' serves to not compute L_1/m
        //  as many times as the number of substitutions to do.
        void lineMinusLine(m_num line,SNline<T,tp_size> v);


        // return the larger element (in absolute value) on the given column
        // In case of equality, return the last one (the larger line).
        //   The template type T has to accept arithmetic manipulations
        //   like abs, comparison.
        SNelement<T,tp_size> getLargerOnColumn(m_num col) const;
        // return the largest (absolute value) element under the diagonal
        // on the given column.
        SNelement<T,tp_size> getLargerUnderDiagonal(m_num col) const;

        // From the line number "line", return a line normalized
        // in such a way that the first (non zero) element is 1.
        SNline<T,tp_size> gaussEliminationLine(m_num line);
        T& _at(const m_num,const m_num) override;
        T _get(const m_num,const m_num) const override;

        /** Set the matrix from another one */
        void _set_from(const SNgeneric<T,tp_size>&);
        void set_identity();
    public:
        SNmatrix();
        SNmatrix(const SNmatrix<T,tp_size>&);
        /**
         * Construct a SNmatrix as copy of a generic matrix.
         * Here we copy every elements.
         * */
        SNmatrix(const SNgeneric<T,tp_size>&);
        /** Creates a matrix full of x */
        SNmatrix(const T& x);


        // return the max of the absolute values of all the matrix elements
        T max_norm() const;

        // return the matrix element on given (line,column).
        SNelement<T,tp_size> getElement(m_num line, m_num col) const;

        // swap the lines l1 and l2. This is in-place replacement.
        // The matrix is changed.
        void swapLines(m_num l1,m_num l2);


        /** 
         * return the PLU decomposition as a `SNplu` object.
         */ 
        SNplu<T,tp_size> getPLU() const;

};

// CONSTRUCTORS  -------------------------------------------

template <class T,unsigned int tp_size>
SNmatrix<T,tp_size>::SNmatrix(): data() { };

template <class T,unsigned int tp_size>
SNmatrix<T,tp_size>::SNmatrix(const SNmatrix<T,tp_size>& snm) : data(snm.data)  {};

template <class T,unsigned int tp_size>
SNmatrix<T,tp_size>::SNmatrix(const T& v): 
    data()
{
    for (unsigned int k=0;k<tp_size*tp_size;++k)
    {
            data.at(k)=v;
    }
};

template <class T,unsigned int tp_size>
SNmatrix<T,tp_size>::SNmatrix(const SNgeneric<T,tp_size>& A)
{
    this->_set_from(A);
}

//  SOME ILLEGITIMATE(?) WAYS TO SET THE VALUES OF A MATRIX -----------------

template <class T,unsigned int tp_size>
void SNmatrix<T,tp_size>::_set_from(const SNgeneric<T,tp_size>& A)
{
    for (m_num i=0;i<tp_size;i++)
    {
        for (m_num j=0;j<tp_size;j++)
        {
            this->at(i,j)=A.get(i,j);
        }
    }
}
template <class T,unsigned int tp_size>
void SNmatrix<T,tp_size>::set_identity()
{
    for (m_num i=0;i<tp_size;i++)
    {
        for (m_num j=0;j<tp_size;j++)
        {
            this->at(i,j)=0;
        }
        this->at(i,i)=1;
    }
}

// GETTER METHODS  -------------------------------------------


template <class T,unsigned int tp_size>
SNelement<T,tp_size> SNmatrix<T,tp_size>::getElement(m_num line, m_num col) const
{
    return SNelement<T,tp_size>(line,col,this->get(line,col));
}


// _GET AND _AT METHODS ---------------------------

template <class T,unsigned int tp_size>
T& SNmatrix<T,tp_size>::_at(const m_num i,const m_num j) 
{
    return data.at(j*tp_size+i);
};

template <class T,unsigned int tp_size>
T SNmatrix<T,tp_size>::_get(const m_num i,const m_num j) const
{
    return data.at(j*tp_size+i);
};


// GAUSS'S ELIMINATION METHODS



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
SNelement<T,tp_size> SNmatrix<T,tp_size>::getLargerUnder(m_num f_line, m_num col) const
{
    T max_val=0;
    m_num max_line=0;

    for (m_num line=f_line;line<tp_size;line++)
    {
        if (std::abs(this->get(line,col))>max_val)
        {
            max_val=std::abs(this->get(line,col));
            max_line=line;
        };
    };
    return getElement(max_line,col);
}

template <class T,unsigned int tp_size>
SNelement<T,tp_size> SNmatrix<T,tp_size>::getLargerOnColumn(m_num col) const
{
    return getLargerUnder(0,col);
}

template <class T,unsigned int tp_size>
SNelement<T,tp_size> SNmatrix<T,tp_size>::getLargerUnderDiagonal(m_num col)  const
{
    return getLargerUnder(col,col);
}

template <class T,unsigned int tp_size>
void SNmatrix<T,tp_size>::swapLines(m_num l1, m_num l2)
{
    if (l1!=l2)
    {
        for (m_num col=0;col<tp_size;col++)
        {
            T tmp = this->get(l1,col);
            this->at(l1,col)=this->get(l2,col);
            this->at(l2,col)=tmp;
        }
    }
}


template <class T,unsigned int tp_size>
void SNmatrix<T,tp_size>::lineMinusLine(m_num line,SNline<T,tp_size> v)
{
    for (m_num c=0;c<tp_size;c++)
    {
        this->at(line,c)=this->get(line,c)-v.get(c);
    }
}

template <class T,unsigned int tp_size>
SNline<T,tp_size> SNmatrix<T,tp_size>::gaussEliminationLine(m_num line)
{
    SNline<T,tp_size> l=this->getSNline(line);
    l.makeUnit();
    return l;
}

template <class T,unsigned int tp_size>
SNplu<T,tp_size> SNmatrix<T,tp_size>::getPLU() const

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
    SNplu<T,tp_size> plu;

    Mpermutation<tp_size>& permutation=plu.data_P;

    // progressively become L
    SNmultiGaussian<T,tp_size> mM=SNidentity<T,tp_size>();  
    // this will progressively become U
    SNmatrix<T,tp_size> mU=*this;  


    for (m_num c=0;c<tp_size;c++)
    {
        auto max_el = mU.getLargerUnderDiagonal(c);

        if (max_el.getValue()!=0)   // not a column full of zero's
        {

            // We swap the line 'c' with max_el.line
            MelementaryPermutation<tp_size> el_perm(c,max_el.line);
            permutation=permutation*el_perm; 
            mU.swapLines(c,max_el.line);

            auto G=mU.getGaussian(c);

            // On the first iteration, there is nothing to swap.
            if (c!=0)
            {
                mM.swapLines(c,max_el.line);
                mM*=G.inverse();
            }
            else
            {
                mM=G.inverse();
            }
            auto killing_line=mU.gaussEliminationLine(c);
            for (m_num l=c+1;l<tp_size;l++)
            {
                T m = mU.get(l,c);  // the value to be eliminated

                // TODO : this is not optimal because
                // we already know the first 'c' differences are 0.
                mU.lineMinusLine(l,m*killing_line);
            }
        }
    }
    // at this point, the matrix mU should be the correct one.
    plu._setU(mU);
    plu._setL(mM);
    return plu;
}

#endif
