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
 Define some mathematical intermediate functions.
*/

#ifndef __MATHUTILITIES_H__055253__
#define __MATHUTILITIES_H__055253__

#include "../exceptions/SNexceptions.cpp"


// forward
template <class T,unsigned int tp_size>
class SNgeneric;

template <class U,unsigned int s,class V,unsigned int t>
void checkSizeCompatibility(const SNgeneric<U,s>&A, const SNgeneric<V,t>&B )
{
    // the 'getSize' test is in order to avoid 'unused parameter warning'.
    if (s!=t or A.getSize()!=B.getSize() )
    {
        throw IncompatibleMatrixSizeException(s,t);
    }
}

template <class T,class U>
bool componentWiseeEquality(const T& A,const U& B)
    // test the equality componentwise between to matrix objects.
{
    checkSizeCompatibility(A,B);
    unsigned int s=A.getSize();
    for (unsigned int l=0;l<s;l++)
    {
        for (unsigned int c=0;c<s;c++)
        {
            if (A.get(l,c)!=B.get(l,c))
            {
                return false;
            }
        }
    }
    return true;
}
/** 
   \brief compute the component (i,j) of the matrix product A*B

     \param A a matrix
     \param B a matrix
     \param i the line of the element of 'A*B' to be computed
     \param j the column of the element of 'A*B' to be computed
    
    The size compatibility is not checked since this functions is supposed 
     to be called from something like an operator * that has already checked
     the compatibility.
    
     As stated in the README.md, the type of the elements in the returned matrix
     is the one of the first matrix of the product.
 */
template <class T,class U,unsigned int s,unsigned int t>
T matrixProductComponent(const SNgeneric<T,s>& A,const SNgeneric<U,t>& B,unsigned int i,unsigned int j)
{
    T acc=0;
    for (unsigned int k=0;k<s;k++)
    {
        acc+=A.get(i,k)*B.get(k,j);
    }
    return acc;
}

/** 
 * \brief copy the first lines of a matrix in an other one.
 * 
 * \param ans the matrix in which one copies (non-const).
 * \param A the matrix from which one copies (const).
 * \param max_l the number of lines that have to be copied (const).
 *
 * ```
 * copyFirstLines(A,B,5);
 * ```
 * Now the 5 first lines in `A` are the same of the ones of `B`. The previously present elements
 * in the first 5 lines of `A` are lost.
 *
 * This function is useful to populate the "answer" matrix when multiplying a gaussian matrix
 * by something else.
 *
 * Remarks :
 * 
 * - the matrices `ans` and `A` must have same size;
 * - for the template type, `ans` is `SNgeneric<U,s>` and `A` 
 *   is SNgeneric<T,tp_size>. Obviously
 *   the conversion from `T` to `U` must be possible.
 *
 * \see SNmultiGaussian<T,tp_size>::setFirstLines(const SNgeneric<U,s>& other,const m_num& max_l);
 */
template <class U, unsigned int t,class T,unsigned int tp_size>
void copyFirstLines(SNgeneric<U,t>& ans, const SNgeneric<T,tp_size>& A,const m_num& max_l)
{
    checkSizeCompatibility(ans,A);
    for (m_num line=0;line<max_l+1;++line)
    {
        for (m_num col=0; col < tp_size;++col)
        {
            ans.at(line,col)=A.get(line,col);
        }
    }
}

#endif
