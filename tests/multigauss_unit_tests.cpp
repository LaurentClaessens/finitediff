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

#include "../src/SNmatrices/SNmultiGaussian.h"
#include "../src/SNmatrices/SNidentity.h"
#include "TestMatrices.cpp"

class multigaussTests : public CppUnit::TestCase
{
    private :
        void working_tests()
        {
            echo_function_test("working tests");
            auto A=testMatrixI();
    /*
    2 5 6 7
    5 1 7 1
    3 2 3 4
    4 6 8 9
    */

            SNmultiGaussian<double,4> mg(A);
            SNmatrix<double,4> prod=mg*A;

            CPPUNIT_ASSERT(prod.get(1,0)==0);
            CPPUNIT_ASSERT(prod.get(2,0)==0);
            CPPUNIT_ASSERT(prod.get(3,0)==0);

        }
        void get_at_tests()
        {
            echo_function_test("get_at_tests");
            auto A=testMatrixI();
            SNmultiGaussian<double,4> mg(A);

            echo_single_test("A simple get");
            CPPUNIT_ASSERT(mg.get(3,0)==-2);

            echo_single_test("Modifying with 'at'");
            mg.at(2,0)=12;
            CPPUNIT_ASSERT(mg.get(2,0)==12);

        }
        void constructor_tests()
        {
            echo_function_test("constructor_tests");
            auto E=testMatrixE();
            double epsilon(0.0000001);

            echo_single_test("Constructor as gaussian matrix for zeroth column");
            
            SNmultiGaussian<double,4> mg(E);
            CPPUNIT_ASSERT(mg.isNumericallyEqual(E.getGaussian(0),epsilon));

            auto M0=E.getGaussian(0);
            auto M1=E.getGaussian(1);
            auto M2=E.getGaussian(2);

            auto i_mg=mg;       // we recall how was out multigaussian before the multiplication by M1.
            mg*=M1;

            // The point in passing through SNmatrix and testing
            // sn_mg+sn_M1-SNidentity
            // instead of testing directly
            // mg+M1-SNidentity
            // is that I do not want to define the sum of gaussian matrices.
            // As far as I know, there are no legitimate situations in which you want to 
            // sum such matrices.
            echo_single_test("SNmatrix from gaussian");
            SNmatrix<double,4> sn_i_mg(i_mg);
            SNmatrix<double,4> sn_M1(M1);
            CPPUNIT_ASSERT(i_mg==sn_i_mg);
            CPPUNIT_ASSERT(M1==sn_M1);

            echo_single_test("Product with a gaussian");
            CPPUNIT_ASSERT(mg.isNumericallyEqual( sn_i_mg+sn_M1-SNidentity<double,4>(),epsilon  ));
        }
    public:
        void runTest()
        {
            working_tests();
            get_at_tests();
            constructor_tests();
        }
};

int main ()
{
    std::cout<<"multigaussTests"<<std::endl;
    multigaussTests mg_test;
    mg_test.runTest();
}
