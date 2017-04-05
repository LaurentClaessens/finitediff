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
#include "../src/SNmatrices/SNoperators.h"
#include "TestMatrices.cpp"

class multigaussTests : public CppUnit::TestCase
{
    private :
        void wrong_order_works()
        {
            echo_function_test("wrong_order_works");
            auto E=testMatrixE();

            double epsilon(0.0001);
   /* E=
    4 6 8 9
    5 1 7 1
    3 2 3 4
    2 5 6 7
    */
            debug_print<<"E"<<std::endl;
            debug_print<<E<<std::endl;

            auto G1=E.getGaussian(1);
            auto G2=E.getGaussian(2);

            debug_matrix_print("G2",G2);
            debug_matrix_print("G1",G1);

            debug_matrix_print("G2*E",G2*E);

            SNmultiGaussian<double,4> ans12;
            ans12.setLastColumn(2);
            ans12.at(2,1)=-2;
            ans12.at(3,1)=-5;
            ans12.at(3,2)=-2;

            SNmultiGaussian<double,4> ans21;
            ans21.setLastColumn(2);
            ans21.at(2,1)=-2;
            ans21.at(3,1)=-1;
            ans21.at(3,2)=-2;


            debug_print<<"G1*G2"<<std::endl;
            debug_print<<G1*G2<<std::endl;
            debug_print<<"ans12"<<std::endl;
            debug_print<<ans12<<std::endl;

            echo_single_test("Product of two gaussians -- right order");
            CPPUNIT_ASSERT(ans12.isNumericallyEqual(G1*G2,epsilon) );

            echo_single_test("Product of two gaussians -- wrong order");
            CPPUNIT_ASSERT(ans21.isNumericallyEqual(G2*G1,epsilon) );

        }
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
        void associativity_check()
        {
            echo_function_test("associativity_check");

            auto E=testMatrixE();       //SNmatrix
            double epsilon(0.0000001);

            SNgaussian<double,4> G0(E);

            auto G1=(G0*E).getGaussian(1);
            auto G2=(G1*G0*E).getGaussian(2);

            SNmultiGaussian<double,4> mg=G0*(G1*G2);

            
            echo_single_test("G0*G1*G2");
            CPPUNIT_ASSERT(mg.isNumericallyEqual(G0*G1*G2,epsilon));
            echo_single_test("(G0*G1)*G2");
            CPPUNIT_ASSERT(mg==(G0*G1)*G2);
        }
        void multi_working_tests()
        {
            echo_function_test("multi working tests");

            auto E=testMatrixE();
            double epsilon(0.0000001);

            SNgaussian<double,4> G0(E);

            auto G1=(G0*E).getGaussian(1);
            auto G2=(G1*G0*E).getGaussian(2);


            SNmultiGaussian<double,4> mg=G2*(G1*G0);

            auto prod=mg*E;


            CPPUNIT_ASSERT(prod.get(1,0)==0);
            CPPUNIT_ASSERT(prod.get(2,0)==0);
            CPPUNIT_ASSERT(prod.get(3,0)==0);
            CPPUNIT_ASSERT(prod.get(2,1)==0);
            CPPUNIT_ASSERT(prod.get(3,1)==0);
            CPPUNIT_ASSERT(prod.get(3,2)==0);

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

            auto M1=E.getGaussian(1);

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
        void inverse_tests()
        {
            CPPUNIT_ASSERT(false);      // to be implemented.
        }
        void non_initialized_tests()
        {
            SNmultiGaussian<int,4> A;
            CPPUNIT_ASSERT_THROW(A.at(0,2)=3,NotInitializedMemberException);
            SNgaussian<int,4> B;
            CPPUNIT_ASSERT_THROW(B.at(0,2)=3,NotInitializedMemberException);

            SNmultiGaussian<int,4> C;
            CPPUNIT_ASSERT_THROW(C.setLastColumn(4),OutOfRangeColumnNumber);
            SNmultiGaussian<int,4> D;
            CPPUNIT_ASSERT_THROW(D.setLastColumn(4),OutOfRangeColumnNumber);

            SNmultiGaussian<int,4> E;
            CPPUNIT_ASSERT_THROW(E.setLastColumn(7),OutOfRangeColumnNumber);
            SNmultiGaussian<int,4> F;
            CPPUNIT_ASSERT_THROW(F.setLastColumn(7),OutOfRangeColumnNumber);
        }
    public:
        void runTest()
        {
            non_initialized_tests();
            wrong_order_works();
            associativity_check();
            multi_working_tests();
            working_tests();
            get_at_tests();
            constructor_tests();
            inverse_tests();
        }
};

int main ()
{
    std::cout<<"multigaussTests"<<std::endl;
    multigaussTests mg_test;
    mg_test.runTest();
}
