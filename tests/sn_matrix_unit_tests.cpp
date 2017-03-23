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

#include <cppunit/TestCase.h>
#include <cppunit/extensions/TypeInfoHelper.h>
#include <cppunit/TestAssert.h>

#include "../src/SNmatrices/SNmatrix.h"

#include "TestMatrices.cpp"

class SNmatrixTest : public CppUnit::TestCase
{
    private :
        void test_instantiate()
        {
            SNmatrix<double,7> sn1;
            CPPUNIT_ASSERT(sn1.getSize()==7);
        }
        void test_element_reference()
            // test if the reference to the original matrix in the matrix element works
        {
            echo_function_test("test_element_reference");
            SNmatrix<double,3> A;
            A.at(0,0)=1;
            A.at(0,1)=3/5;
            A.at(0,2)=3;
            A.at(1,0)=2;
            A.at(1,1)=5;
            A.at(1,2)=0;
            A.at(2,0)=3;
            A.at(2,1)=8;
            A.at(2,2)=0;

            auto el1=A.getElement(0,0);
            CPPUNIT_ASSERT(el1.getValue()==1);

            auto el2=A.getElement(0,1);
            CPPUNIT_ASSERT(el2.getValue()==3/5);
        }
        void test_populate()
        {
            echo_function_test("test_populate");
            SNmatrix<int,2> sn2;
            sn2.at(0,0)=1;
            sn2.at(0,1)=20;
            sn2.at(1,0)=3;
            sn2.at(1,1)=4;
            sn2.at(0,1)=2;
            CPPUNIT_ASSERT(sn2.getSize()==2);
            CPPUNIT_ASSERT(sn2.at(0,0)==1);
            CPPUNIT_ASSERT(sn2.at(0,1)==2);
        }
        void test_max_on_column()
            // test if one get the right maximal element on a column 
            // (for the Gauss's pivot)
        {
            echo_function_test("test_max_on_column");
            auto A=testMatrixC();

    /*
     1  -3  5
     2  1   6
     1  2   6.1
    */

            echo_single_test("test the max value on a column");
            auto max0=A.getLargerOnColumn(0);
            auto max1=A.getLargerOnColumn(1);
            auto max2=A.getLargerOnColumn(2);
            CPPUNIT_ASSERT(max0.getValue()==2);
            CPPUNIT_ASSERT(max1.getValue()==-3);
            CPPUNIT_ASSERT(max2.getValue()==6.1);


            A.at(2,2)=-2; 
            auto bmax0=A.getLargerUnderDiagonal(0);           // By the way, test the assignation operator.
            auto bmax1=A.getLargerUnderDiagonal(1);
            auto bmax2=A.getLargerUnderDiagonal(2);

            echo_single_test("test the max value under the diagonal");
            CPPUNIT_ASSERT(bmax0.getValue()==2);
            CPPUNIT_ASSERT(bmax1.getValue()==2);
            CPPUNIT_ASSERT(bmax2.getValue()==-2);
        }
        void test_copy_constructor()
            // test the copy constructor and the equality operator.
        {

            echo_function_test("test_copy_constructor");
            SNmatrix<int,2> A;
            A.at(0,0)=1;
            A.at(0,1)=2;
            A.at(1,0)=3;
            A.at(1,1)=4;

            auto B(A);
            echo_single_test("equality between the matrix and the copied one");
            CPPUNIT_ASSERT(A==B);
        }
        void test_swap_line()
        {
        auto A=testMatrixC();
    /*
     1  -3  5
     2  1   6
     1  2   6.1
    */

            auto B(A);

            A.swapLines(0,1);
            CPPUNIT_ASSERT(A.at(0,1)==1);
            A.swapLines(0,1);
            CPPUNIT_ASSERT(A==B);
        }
    void test_max_norm()
    {
        echo_function_test("test_max_norm");
        auto A=testMatrixC();
    /*
     1  -3  5
     2  1   6
     1  2   6.1
    */

        CPPUNIT_ASSERT(A.max_norm()==6.1);
        A.at(2,2)=-6.1;
        CPPUNIT_ASSERT(A.max_norm()==6.1);
        A.at(2,2)=-1;
        A.at(1,2)=1;
        CPPUNIT_ASSERT(A.max_norm()==5);
    }

    public :
        void runTest()
        {
            test_instantiate();
            test_populate();
            test_element_reference();
            test_max_on_column();
            test_copy_constructor();
            test_swap_line();
            test_max_norm();
        }
};

int main ()
{
    std::cout<<"SNmatrixTest"<<std::endl;
    SNmatrixTest sn_test;
    sn_test.runTest();
}
