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


/*
* This file contains the multiplication operators between the different types
* of matrices.
*/

#ifndef  MULTIPLICATIONS_H__14063_
#define  MULTIPLICATIONS_H__14063_

#include "../SNpermutation.h"
#include "../SNgaussian.h"
#include "../SNmultiGaussian.h"
#include "../SNidentity.h"
#include "../SNlowerTriangular.h"
#include "../SNupperTriangular.h"
#include "../SNscalar.h"
#include "../MathUtilities.h"
#include "../../exceptions/SNexceptions.cpp"

/**
* \brief `SNgeneric` * `SNgeneric`.
*
* In general, I cannot do better than compute everything. 
*  
* The very point of making many different classes is NOT to use this 'default' implementation for the product.
*
* When this product is used, a warning is displayed.
*
* \see `tooGenericWarning`
*/
template <class U,class V,unsigned int s,unsigned int t>
SNmatrix<U,s> operator*(const SNgeneric<U,s>& A, const SNgeneric<V,t>& B)
{
    tooGenericWarning("Warning : using a very generic product 'SNgeneric * SNgeneric'. Can't you be more specific ?");
    checkSizeCompatibility(A,B);
    SNmatrix<U,s> ans;
    for (unsigned int i=0;i<s;i++)
    {
        for (unsigned int j=0;j<s;j++)
        {
            ans.at(i,j)=matrixProductComponent(A,B,i,j);
        }
    }
    return ans;   //relies on RVO.
}

// number * identity

/** 
 * @brief multiply a number by the identity matrix.
 *
 * This produces a scalar matrix.
 *
 * The resulting template type is the one of the number.
 * */
template <class U,class V,unsigned int s>
SNscalar<U,s> operator*
(const U& x, const SNidentity<V,s>& M)
{
    return SNscalar<U,s>(x);
}

// SNgaussian * SNgeneric

/** 
 * \brief Product `SNgaussian` * `SNgeneric`
 *
 * As far as the template parameters are concerned,
 * the answer is `SNmatrix<T,tp_size>` with
 *
 * - `T` is the type of the gaussian (the left operand)
 * - `tp_size` is the common size of the two matrices.
 *
 * Let \f$ G \f$ be gaussian with non trivial column \f$ c \f$ and \f$ E \f$
 * be generic. For the product \f$ GE \f$ one can copy the first \f$ c \f$ lines.
 *
 * For the other lines (\f$ i>c \f$) the sum \f$ (GE)_{ij}=\sum_kG_{ik}E_{kj} \f$
 * is non vanishing only with \f$ k=i \f$ and \f$ k=c \f$.
 *
 * */
template <class U,class V,unsigned int s,unsigned int t>
SNmatrix<U,s> operator*
(const SNgaussian<U,s>& A, const SNgeneric<V,t>& B)
{

    checkSizeCompatibility(A,B);

    SNmatrix<U,s> ans;

    const unsigned int c=A.getColumn();

    copyFirstLines(ans,B,c);

    // for the other lines, one has only two non vanishing elements
    // in the gaussian.
    for (m_num line=c+1;line<s;++line)
    {
        for (m_num col=0;col<s;++col)
        {
            ans.at(line,col)=B.get(line,col)+A.get(line,c)*B.get(c,col);
        }
    }
    return ans;
}


// SNgaussian * SNgaussian

template <class U,class V,unsigned int s,unsigned int t>
SNmultiGaussian<U,s> operator*
(const SNgaussian<U,s>& A, const SNgaussian<V,t>& B)
{
    checkSizeCompatibility(A,B);
    unsigned int tp_size=s;

    SNmultiGaussian<U,s> ans;
    if (A.getColumn()==B.getColumn())
    {
        ans.setLastColumn(A.getColumn());
        m_num col=A.getColumn();
        for (m_num line = col+1 ; line< tp_size ;++line )
        {
            ans.at(line,col)=A.get(line,col)+B.get(line,col);
        }
    }
    else if (A.getColumn()>B.getColumn())
    {
        ans.setLastColumn(A.getColumn());

            // The `_at` function in in `SNmultigauss` automatically 
            // returns '0' when the column number is larger than 
            // `getLastColumn`. In other words, these are 
            // "special values" and attempting to access them with `_at`
            // throws a `SNchangeNotAllowedException`. 
        for (m_num col=0;col <= ans.getLastColumn() ;++col)
        {
            if (col==A.getColumn())
            {
                for (m_num line=col+1;line<s;++line)
                {
                    ans.at(line,col)=A.get(line,col);
                }
            }
            else if (col==B.getColumn())
            {
                for (m_num line=B.getColumn()+1;line<A.getColumn()+1;++line)
                {
                    ans.at(line,col)=B.get(line,col);
                }
                for (m_num line=A.getColumn()+1;line<s;++line)
                {
                    ans.at(line,col)=A.get(line,A.getColumn())*B.get(A.getColumn(),col)+B.get(line,B.getColumn());
                }
            }
            else
            {
                for (m_num line=col+1;line<s;++line)
                {
                    ans.at(line,col)=0;
                }
            }
        }
    }
    else
    {
        ans.setLastColumn(B.getColumn());

        for (m_num c=0;c <= ans.getLastColumn();c++)
        {
            for (m_num l=c+1;l<s;l++)
            {
                ans.at(l,c)=A.get(l,c)+B.get(l,c);
            }
        }
    }
    return ans;
}

// SNmultiGaussian * SNgaussian

template <class U,class V,unsigned int s,unsigned int t>
SNmultiGaussian<U,s> operator*
(const SNmultiGaussian<U,s>& M, const SNgaussian<V,t>& G)
{
    checkSizeCompatibility(M,G);
    SNmultiGaussian<U,s> ans(M);

    if (M.getLastColumn()  >=  G.getColumn())
    {
        ans.setLastColumn(M.getLastColumn());

        m_num col=G.getColumn();
        for (m_num line=col+1;line<s;++line)
        {
            ans.at(line,col)=matrixProductComponent(M,G,line,col);
        }
    }
    else
    {
        ans.setLastColumn(G.getColumn());

        for (m_num l=G.getColumn()+1;l<s;l++)
        {
            ans.at(l,G.getColumn())+=G.get(l,G.getColumn());
        }
    }
    return ans;
}

// SNmultiGaussian * SNgeneric

/** 
 *\brief Product `SNmultiGaussian` * `SNgeneric`
 *
 * As far as the template parameters are concerned,
 * the answer is `SNmatrix<T,tp_size>` with
 *
 * - `T` is the type of the multi-gaussian (the left operand)
 * - `tp_size` is the common size of the two matrices.
 *
 * We populate the answer line by line. Let \f$ M \f$ be a multi-gaussian
 * and \f$ E \f$ a generic matrix. We denote by \f$ c \f$ the last non-trivial
 * column of \f$ M \f$.
 *
 * For computing \f$ (ME)_{ij}=\sum_kM_{ik}E_{kj} \f$ :
 *
 * - the first line is copied.
 * - if \f$ i<c \f$ we compute the product with the first \f$ i-1 \f$ elements, 
 *   and add the \f$ i \f$th.
 * - if \f$ i\geq c \f$ we compute the first \f$ c-1 \f$ products and 
 *   add the \f$ i \f$th.
 *
 * */
template <class U,class V,unsigned int s,unsigned int t>
SNmatrix<U,s> operator*
(const SNmultiGaussian<U,s>& M, const SNgeneric<V,t>& E)
{
    
    checkSizeCompatibility(M,E);
    const unsigned int tp_size=M.getSize(); // for homogeneous notations.
    const m_num last_col=M.getLastColumn();

    SNmatrix<U,s> ans;

    // copy the first line

    for (m_num col=0;col<tp_size;++col)
    {
        ans.at(0,col)=E.get(0,col);
    }

    // when the line number is smaller than 'last_col'
    for (m_num line=1;line <= last_col;++line)
    {
        for (m_num col=0;col<tp_size;++col)
        {
            U acc=0;
            for (m_num k=0; k < line;++k)
            {
                acc+=M.get(line,k)*E.get(k,col);
            }
            ans.at(line,col)=acc+E.get(line,col);
        }
    }

    // for the last lines
    for (m_num line=last_col+1;line<tp_size;++line)
    {
        for (m_num col=0;col<tp_size;++col)
        {
            U acc=0;
            for (m_num k=0;k <= last_col;++k)
            {
                acc+=M.get(line,k)*E.get(k,col);
            }
            ans.at(line,col)=acc+E.get(line,col);
        }
    }
    return ans;
}


// SNmultiGaussian * SNmultigaussian

template <class U,class V,unsigned int s,unsigned int t>
SNmultiGaussian<U,s> operator*
(const SNmultiGaussian<U,s>& A, const SNmultiGaussian<V,t>& B)
{
    checkSizeCompatibility(A,B);
    SNmultiGaussian<U,s> ans;
    ans.setLastColumn(std::max(A.getLastColumn(),B.getLastColumn()));

    for (m_num col=0;col<ans.getLastColumn();++col)
    {
        for (m_num line=col+1;line<s;++line)
        {
            U acc=0;
            // TODO : non optimal because the first and last products are 1*something and something*1.
            for (m_num k=col;k <= line;++k)
            {
                acc+=(A.get(line,k)*B.get(k,col));
            }
            ans.at(line,col)=acc;
        }
    }
    return ans;
}

// SNgaussian * SNlowerTriangular

template <class U,class V,unsigned int s,unsigned int t>
SNlowerTriangular<U,s> operator*
(const SNgaussian<U,s>& A, const SNlowerTriangular<V,t>& B)
{

    // gaussian * lower trig -> lower trig
    //
    // Copy the first 'c' lines (only under the diagonal)
    // Then compute the product for the other lines (only under the diagonal)

    checkSizeCompatibility(A,B);
    unsigned int size=A.getSize();
    unsigned int c=A.getColumn();
    SNlowerTriangular<U,s> ans;

    for (unsigned int i=0;i<c+1;i++)
    {
        for (unsigned int j=0;j<i+1;j++)
        {
            ans.at(i,j)=B.get(i,j);
        }
    }
    for (unsigned int i=c+1;i<size;i++)
    {
        for (unsigned int j=0;j<i+1;j++)
        {
            ans.at(i,j)=matrixProductComponent(A,B,i,j);
        }
    }
    return ans;
}

// SNgaussian * SNmultiGaussian

/** 
 * \brief Product `SNgaussian` * `SNmultiGaussian`
 *
 * The answer is `SNmultiGaussian<T,tp_size>` with
 *
 * - `T` is the type of the gaussian (the left operand)
 * - `tp_size` is the common size of the two matrices.
 * - the last non trivial column is the max of the column of the first argument and
 *   of the last non trivial of the second argument.
 *
 *
 * Let \f$ G \f$ be a gaussian matrix for the column \f$ col \f$ and \f$ M \f$ 
 * a multi-gaussian with max column \f$ l_c \f$.
 *
 * For computing the element \f$ (GM)_{ij}=\sum_kG_{ik}M{kj} \f$ we use the 
 * structure of \f$ G \f$.
 *
 * - If \f$ i\leq col \f$ we have \f$ (GM)_{ij}=M_{ij} \f$. So we copy the first 
 *   \f$ col \f$ lines.
 * - The sum has only two non vanishing terms and we have
 *   \f$ (GM)_{ij}=M_{ij}+G_{i,col}M_{col,j} \f$.
 *
 * */
template <class U,class V,unsigned int s,unsigned int t>
SNmultiGaussian<U,s> operator*
(const SNgaussian<U,s>& G, const SNmultiGaussian<V,t>& M)
{

    checkSizeCompatibility(G,M);

    const unsigned int tp_size=G.getSize(); // for homogeneity
    const m_num col=G.getColumn();
    const m_num last_col=M.getLastColumn();

    SNmultiGaussian<U,s> ans;
    ans.setLastColumn(  std::max(col,last_col)  );

    ans.setFirstLines(M,col);

    for (m_num i=col+1;i<tp_size;++i)   // loop over the next lines
    {
        for (m_num j=0;j<i;++j)
        {
            ans.at(i,j)=M.get(i,j)+G.get(i,col)*M.get(col,j);
        }
    }
    return ans;
}

// Mpermutation * Mpermutation

/** 
 * The multiplication "permutation1 * permutation2" 
 * is the composition. 
*/
template <unsigned int tp_size>
Mpermutation<tp_size> operator*
(const Mpermutation<tp_size>& p1, const Mpermutation<tp_size>& p2)
{
    Mpermutation<tp_size> new_perm;
    for (unsigned int i=0;i<tp_size;++i)
    {
        new_perm.at(i)=p1.image( p2.image(i) );
    }
    return new_perm;
}

/** 
 *\brief Product `MgenericPermutation` * `MgenericPermutation`
 *
 * The product is the composition.
 *
 * This product defines the products of `Mpermutation` and `MelementaryPermutation` (there are 4 possibilities).
 * 
 * */
template <unsigned int tp_size>
Mpermutation<tp_size> operator*(const MgenericPermutation<tp_size>& A, const MgenericPermutation<tp_size>& B) 
{
    Mpermutation<tp_size> tmp;
    for (unsigned int k=0;k<tp_size;++k)
    {
        tmp.at(k)=A.image( B.image(k)  );
    }
    return tmp;
}

#endif
