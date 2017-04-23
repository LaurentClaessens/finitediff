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


/** 
 * This is an example of PLU decomposition.
 *
 * How to use it for your own program ? See [here](examples/README.md).
 * */

#include <iostream>
#include "../src/SNplu.h"

int main ()
{

    // Here is how to create the following matrix :

    // 1 3 4 5 6
    // 2 5 7 8 9
    // 3 4 6 2 0
    // 6 7 8 7 1
    // 9 2 3 0 1

    SNmatrix<double,5> A;       // creating

    // populating with 'at'

    A.at(0,0)=1; A.at(0,1)=3; A.at(0,2)=4; A.at(0,3)=5;A.at(0,4)=6;
    A.at(1,0)=2; A.at(1,1)=5; A.at(1,2)=7; A.at(1,3)=8;A.at(1,4)=9;
    A.at(2,0)=3; A.at(2,1)=3; A.at(2,2)=6; A.at(2,3)=2;A.at(2,4)=0;
    A.at(3,0)=6; A.at(3,1)=7; A.at(3,2)=8; A.at(3,3)=7;A.at(3,4)=1;
    A.at(4,0)=9; A.at(4,1)=2; A.at(4,2)=3; A.at(4,3)=0;A.at(4,4)=1;
    
    std::cout<<"The matrix A is "<<std::endl;
    std::cout<<A<<std::endl;            // printing the matrix


    // compute the PLU decomposition :
    auto plu=A.getPLU();

    // get the matrices :

    auto cP=plu.getP();     // permutation
    auto cL=plu.getL();     // lower triangular
    auto cU=plu.getU();     // upper triangular

    // print :
    std::cout<<"P"<<std::endl;
    std::cout<<cP<<std::endl;
    std::cout<<"L"<<std::endl;
    std::cout<<cL<<std::endl;
    std::cout<<"U"<<std::endl;
    std::cout<<cU<<std::endl;

}
