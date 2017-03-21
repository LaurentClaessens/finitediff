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

class MultiplicationTest : public CppUnit::TestCase
{
    private :
        void test_gauss_times_matrix()
        {
            echo_function_test("test_gauss_times_matrix");
            auto F=testMatrixF();
            auto G=F.getGaussian(1);    // the gaussian matrix of F

// computations from Sage :
//  
// sage: F=matrix(  [ [1,0,3,9],[6,2,3,5],[7,8,1,3],[7,7,4,6] ]  )
// sage: G=matrix(   [     [1,0,0,0],[0,1,0,0],[0,-4,1,0],[0,-7/2,0,1]   ]   )
// sage: G*F
//[    1     0     3     9]
//[    6     2     3     5]
//[  -17     0   -11   -17]
//[  -14     0 -13/2 -23/2]

            SNmatrix<double,4> ans_G;
            ans_G.at(0,0)=1; ans_G.at(0,1)=0; ans_G.at(0,2)=0; ans_G.at(0,3)=0;
            ans_G.at(1,0)=0; ans_G.at(1,1)=1; ans_G.at(1,2)=0; ans_G.at(1,3)=0;
            ans_G.at(2,0)=0; ans_G.at(2,1)=-4; ans_G.at(2,2)=1; ans_G.at(2,3)=0;
            ans_G.at(3,0)=0; ans_G.at(3,1)=-7./2; ans_G.at(3,2)=0; ans_G.at(3,3)=1;

            double epsilon(0.00001);
            ans_G.subtract(G); 
            CPPUNIT_ASSERT(ans_G.max_norm()<epsilon);

            SNmatrix<double,4> GstarF;
            GstarF.at(0,0)=1; GstarF.at(0,1)=0; GstarF.at(0,2)=3; GstarF.at(0,3)=9;
            GstarF.at(1,0)=6; GstarF.at(1,1)=2; GstarF.at(1,2)=3; GstarF.at(1,3)=5;
            GstarF.at(2,0)=-17; GstarF.at(2,1)=0; GstarF.at(2,2)=-11; GstarF.at(2,3)=-17;
            GstarF.at(3,0)=-14; GstarF.at(3,1)=0; GstarF.at(3,2)=-13./2; GstarF.at(3,3)=-23./2;


            CPPUNIT_ASSERT(GstarF.get(2,1)<epsilon);            // test that the gaussian matrix did the work.
            CPPUNIT_ASSERT(GstarF.get(3,1)<epsilon);

            GstarF.subtract(G*F);
            CPPUNIT_ASSERT(GstarF.max_norm()<epsilon);
        }
        void test_gauss_times_lower_trig()
        {
            echo_function_test("test_gauss_times_lower_trig");
            auto A=testMatrixG();           // lower triangular
    // A=
    // 3 0 0 0
    // 1 2 0 0
    // 3 4 5 0
    // 6 7 8 1

            auto G=A.getGaussian(0); 

            SNmatrix<double,4> ans_G;
            ans_G.at(0,0)=1; ans_G.at(0,1)=0; ans_G.at(0,2)=0; ans_G.at(0,3)=0;
            ans_G.at(1,0)=-1./3; ans_G.at(1,1)=1; ans_G.at(1,2)=0; ans_G.at(1,3)=0;
            ans_G.at(2,0)=-1; ans_G.at(2,1)=0; ans_G.at(2,2)=1; ans_G.at(2,3)=0;
            ans_G.at(3,0)=-2; ans_G.at(3,1)=0; ans_G.at(3,2)=0; ans_G.at(3,3)=1;

            double epsilon(0.00001);
            echo_single_test("G numerically equal to ans_G");
            CPPUNIT_ASSERT(  G.isNumericallyEqual(ans_G,epsilon)  );

            SNmatrix<double,4> GstarA;
//[3 0 0 0]
//[0 2 0 0]
//[0 4 5 0]
//[0 7 8 1]
            GstarA.at(0,0)=3; GstarA.at(0,1)=0; GstarA.at(0,2)=0; GstarA.at(0,3)=0;
            GstarA.at(1,0)=0; GstarA.at(1,1)=2; GstarA.at(1,2)=0; GstarA.at(1,3)=0;
            GstarA.at(2,0)=0; GstarA.at(2,1)=4; GstarA.at(2,2)=5; GstarA.at(2,3)=0;
            GstarA.at(3,0)=0; GstarA.at(3,1)=7; GstarA.at(3,2)=8; GstarA.at(3,3)=1;

            echo_single_test("GstarA numerically equal to G*A");

            CPPUNIT_ASSERT(  GstarA.isNumericallyEqual(G*A,epsilon)  );
        }


    public :
        void runTest()
        {
            test_gauss_times_matrix();
            test_gauss_times_lower_trig();
        }
};

int main ()
{
    std::cout<<"Matrix multiplication tests"<<std::endl;
    MultiplicationTest mul_test;
    mul_test.runTest();
}
