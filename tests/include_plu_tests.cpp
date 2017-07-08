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
 * This test file tests the fact to be able to compute and print a PLU decompisition
 * with only one `include`.
 *
 * It does not tests the correctness of the PLU itself.
 * */

#include "../src/SNplu.h"

int main ()
{
    SNmatrix<double,5> A;       // creating
    A.at(0,0)=1; A.at(0,1)=3; A.at(0,2)=4; A.at(0,3)=5;A.at(0,4)=6;
    A.at(1,0)=2; A.at(1,1)=5; A.at(1,2)=7; A.at(1,3)=8;A.at(1,4)=9;
    A.at(2,0)=3; A.at(2,1)=3; A.at(2,2)=6; A.at(2,3)=2;A.at(2,4)=0;
    A.at(3,0)=6; A.at(3,1)=7; A.at(3,2)=8; A.at(3,3)=7;A.at(3,4)=1;
    A.at(4,0)=9; A.at(4,1)=2; A.at(4,2)=3; A.at(4,3)=0;A.at(4,4)=1;
    
    std::cout<<"The matrix A is "<<std::endl;
    std::cout<<A<<std::endl;            // printing the matrix

    auto plu=A.getPLU();
    auto cP=plu.getP();     // permutation
    auto cL=plu.getL();     // lower triangular
    auto cU=plu.getU();     // upper triangular

    std::cout<<"P"<<std::endl;
    std::cout<<cP<<std::endl;
    std::cout<<"L"<<std::endl;
    std::cout<<cL<<std::endl;
    std::cout<<"U"<<std::endl;
    std::cout<<cU<<std::endl;
}
