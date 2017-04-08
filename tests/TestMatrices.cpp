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

#include "../src/SNmatrices/SNmatrix.h"


/*
 This file contains some matrices ready to be tested.
*/

// Copy-paste for creating own matrix
    //FOOOO.at(0,0)=<++>; FOO.at(0,1)=<++>; FOO.at(0,2)=<++>; FOO.at(0,3)=<++>;
    //FOOOO.at(1,0)=<++>; FOO.at(1,1)=<++>; FOO.at(1,2)=<++>; FOO.at(1,3)=<++>;
    //FOOOO.at(2,0)=<++>; FOO.at(2,1)=<++>; FOO.at(2,2)=<++>; FOO.at(2,3)=<++>;
    //FOOOO.at(3,0)=<++>; FOO.at(3,1)=<++>; FOO.at(3,2)=<++>; FOO.at(3,3)=<++>;

auto testMatrixA()
/*
 0  3  6
 0  0  2
 0  3  6

 This matrix is not invertible. So it does not makes really sense to
 compute its PLU decomposition.
*/
{
    SNmatrix<double,3> B; 
    B.at(0,0)=0; B.at(0,1)=3; B.at(0,2)=6; 
    B.at(1,0)=0; B.at(1,1)=0; B.at(1,2)=2;
    B.at(2,0)=0; B.at(2,1)=3; B.at(2,2)=6;
    return B;
}

auto testMatrixA_U()
/*
 0  3  6
 0  3  6
 0  0  2
*/
{
    SNmatrix<double,3> B; 
    B.at(0,0)=0; B.at(0,1)=3; B.at(0,2)=6; 
    B.at(1,0)=0; B.at(1,1)=3; B.at(1,2)=6;
    B.at(2,0)=0; B.at(2,1)=0; B.at(2,2)=2;
    return B;
}

auto testMatrixB()
/*
 1  2  3
 2  5  0
 3  8  0
*/
{
    SNmatrix<double,3> B; 
    B.at(0,0)=1; B.at(0,1)=2; B.at(0,2)=3; 
    B.at(1,0)=2; B.at(1,1)=5; B.at(1,2)=0;
    B.at(2,0)=3; B.at(2,1)=8; B.at(2,2)=0;
    return B;
}

auto testMatrixB_U()
    // this one is the upper triangular of testMatrixB
/*
 1    2     3
 0  -2/3    3
 0    0    -3/2
*/
{
    SNmatrix<double,3> B; 
    B.at(0,0)=3; B.at(0,1)=8; B.at(0,2)=0; 
    B.at(1,0)=0; B.at(1,1)=-2.0/3; B.at(1,2)=3;
    B.at(2,0)=0; B.at(2,1)=0; B.at(2,2)=-3.0/2;
    return B;
}

auto testMatrixC()
    /*
     1  -3  5
     2  1   6
     1  2   6.1
    */
{
    SNmatrix<double,3> A;
    A.at(0,0)=1;
    A.at(0,1)=-3;       // Test negative and first line
    A.at(0,2)=5;
    A.at(1,0)=2;
    A.at(1,1)=1;
    A.at(1,2)=6;
    A.at(2,0)=1;
    A.at(2,1)=2;
    A.at(2,2)=6.1;      // test last line and decimal value
    return A;
}

auto testMatrixD()
    /*
    1 2 3
    4 5 6
    7 8 9
    */
{
    SNmatrix<double,3> A;
    A.at(0,0)=1; A.at(0,1)=2; A.at(0,2)=3;
    A.at(1,0)=4; A.at(1,1)=5; A.at(1,2)=6;
    A.at(2,0)=7; A.at(2,1)=8; A.at(2,2)=9;   
    return A;
}

auto testMatrixD_U()
    /*
    7    8      9
    0   6/7   12/7
    0    0      0
    */
{
    SNmatrix<double,3> A;
    A.at(0,0)=7; A.at(0,1)=8; A.at(0,2)=9;
    A.at(1,0)=0; A.at(1,1)=6./7; A.at(1,2)=12./7;
    A.at(2,0)=0; A.at(2,1)=0; A.at(2,2)=0;   
    return A;
}

auto testMatrixE()
    /*
    4 6 8 9
    5 1 7 1
    3 2 3 4
    2 5 6 7
    */
{
    SNmatrix<double,4> A;
    A.at(0,0)=4; A.at(0,1)=6; A.at(0,2)=8; A.at(0,3)=9;
    A.at(1,0)=5; A.at(1,1)=1; A.at(1,2)=7; A.at(1,3)=1;
    A.at(2,0)=3; A.at(2,1)=2; A.at(2,2)=3; A.at(2,3)=4;  
    A.at(3,0)=2; A.at(3,1)=5; A.at(3,2)=6; A.at(3,3)=7;  
    return A;
}

auto testMatrixE_U()
{
    SNmatrix<double,4> A;

 // The answer given by Sage :
 //sage: A=matrix(  [   [4,6,8,9],[5,1,7,1],[3,2,3,4],[2,5,6,7]  ]  )
 //sage: print(  A.LU()[2] )
 // 
 //[     5      1      7      1]
 //[     0   26/5   12/5   41/5]
 //[     0      0 -24/13  31/26]
 //[     0      0      0   1/24]


    A.at(0,0)=5; A.at(0,1)=1; A.at(0,2)=7; A.at(0,3)=1;
    A.at(1,0)=0; A.at(1,1)=26./5; A.at(1,2)=12./5; A.at(1,3)=41./5;
    A.at(2,0)=0; A.at(2,1)=0; A.at(2,2)=-24./13; A.at(2,3)=31./26;  
    A.at(3,0)=0; A.at(3,1)=0; A.at(3,2)=0; A.at(3,3)=1./24;  
    return A;
}

auto testMatrixF()
{
    SNmatrix<double,4> F;

    F.at(0,0)=1; F.at(0,1)=0; F.at(0,2)=3; F.at(0,3)=9;
    F.at(1,0)=6; F.at(1,1)=2; F.at(1,2)=3; F.at(1,3)=5;
    F.at(2,0)=7; F.at(2,1)=8; F.at(2,2)=1; F.at(2,3)=3;
    F.at(3,0)=7; F.at(3,1)=7; F.at(3,2)=4; F.at(3,3)=6;
    return F;
}

auto testMatrixG()
{

    // 3 0 0 0
    // 1 2 0 0
    // 3 4 5 0
    // 6 7 8 1

    SNlowerTriangular<double,4> A;
    A.at(0,0)=3; 
    A.at(1,0)=1; A.at(1,1)=2; 
    A.at(2,0)=3; A.at(2,1)=4; A.at(2,2)=5;
    A.at(3,0)=6; A.at(3,1)=7; A.at(3,2)=8; A.at(3,3)=1;
    
    return A;
}

auto testMatrixH()
{
    SNmatrix<double,4> F;
    F.at(0,0)=0.1; F.at(0,1)=-4; F.at(0,2)=10; F.at(0,3)=3;
    F.at(1,0)=-1.1; F.at(1,1)=-6; F.at(1,2)=1./3; F.at(1,3)=4./9;
    F.at(2,0)=10./7; F.at(2,1)=2; F.at(2,2)=3; F.at(2,3)=5;
    F.at(3,0)=4; F.at(3,1)=7; F.at(3,2)=12.34; F.at(3,3)=50;
    
    return F;
}

auto testMatrixI()
    /*
    2 5 6 7
    5 1 7 1
    3 2 3 4
    4 6 8 9
    */
{
    SNmatrix<double,4> A;
    A.at(0,0)=2; A.at(0,1)=5; A.at(0,2)=6; A.at(0,3)=7;  
    A.at(1,0)=5; A.at(1,1)=1; A.at(1,2)=7; A.at(1,3)=1;
    A.at(2,0)=3; A.at(2,1)=2; A.at(2,2)=3; A.at(2,3)=4;  
    A.at(3,0)=4; A.at(3,1)=6; A.at(3,2)=8; A.at(3,3)=9;
    return A;
}


template <class T,unsigned int tp_size>
class AutoTestMatrix
{
    public :
        SNmatrix<T,tp_size> A;
        SNmatrix<T,tp_size> ans_P;
        SNlowerTriangular<T,tp_size> ans_L;
        SNupperTriangular<T,tp_size> ans_U;

        AutoTestMatrix(const SNmatrix<T,tp_size>& mA,const SNmatrix<T,tp_size>& mP, const SNlowerTriangular<T,tp_size>& mL,const SNupperTriangular<T,tp_size>& mU) :
            A(mA),
            ans_P(mP),
            ans_L(mL),
            ans_U(mU)
        { }
};

