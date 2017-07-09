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
#include "../src/SNmatrices/operators/SNoperators.h"
#include "../src/SNmatrices/operators/multiplication.h"
#include "TestMatrices.cpp"
#include "../src/Utilities.h"

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

            auto G1=E.getGaussian(1);
            auto G2=E.getGaussian(2);

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
            
            SNmultiGaussian<double,4> ans;
            ans.setLastColumn(2);

            ans.at(1,0)=-1.25; 
            ans.at(2,0)=-0.75; ans.at(2,1)=-0.38461538461538464; 
            ans.at(3,0)=-0.5; ans.at(3,1)=0.3076923076923077; ans.at(3,2)=0.5833333333333335;


            // the following tests are essentially there to check that the
            // products are well defined. Mainly, they have not to throw 
            // `SNchangeNotAllowedException` due to bad line/column numbering
            // in the loops that are computing the products.

            auto c=G1*G0;
            auto da=G0*G1; 
            auto aa=G1*G2;
            auto ba=G0*E;

            auto a=G1*ba;
            auto b=c*E;
            auto ca=da*G2; // (G0*G1)*G2
            auto cb=G0*aa; // G0*(G1*G2)

            // avoid warning about unused variable.
            CPPUNIT_ASSERT(a.isNumericallyEqual(b,epsilon));
            CPPUNIT_ASSERT(ca.isNumericallyEqual(cb,epsilon));

            echo_single_test("(G0*G1)*G2");
            CPPUNIT_ASSERT(ans.isNumericallyEqual((G0*G1)*G2,epsilon));
            echo_single_test("G0*(G1*G2)");
            CPPUNIT_ASSERT(ans.isNumericallyEqual(G0*(G1*G2),epsilon));
            echo_single_test("G0*G1*G2");
            CPPUNIT_ASSERT(ans.isNumericallyEqual(G0*G1*G2,epsilon));
        }
        void product_tests()
        {
            echo_function_test("product_tests");

            auto E=testMatrixE();       //SNmatrix<double,4>
            auto F=testMatrixF();       //SNmatrix<double,4>
            auto H=testMatrixH();       //SNmatrix<double,4>

            // I cannot ask epsilon to be smaller because for the
            // answers, I copied by hand the print of Sage's answer.
            double epsilon(0.0001);

            auto G0=E.getGaussian(0);
            auto G1=F.getGaussian(1);
            auto G2=H.getGaussian(2);

            auto t21=G2*G1;
            auto t10=G1*G0;
            auto t21_0=t21*G0;  //(G2*G1)*G0
            auto t2_10=G2*t10;  //G2*(G1*G0)
            auto t210=G2*G1*G0;

            SNmultiGaussian<double,4> ans_21;
            ans_21.setLastColumn(2);

    
//[ 1.0    0.0                0.0          0.0]
//[ 0.0    1.0                0.0          0.0]
//[ 0.0   -4.0                1.0          0.0]
//[ 0.0   12.953320000000001  -4.11333     1.0]


            ans_21.at(2,1)=-4;
            ans_21.at(3,1)=12.95332;
            ans_21.at(3,2)=-4.11333;

            SNmultiGaussian<double,4> ans_10;
            ans_10.setLastColumn(1);

//[  1.0   0.0   0.0   0.0]
//[-1.25   1.0   0.0   0.0]
//[ 4.25  -4.0   1.0   0.0]
//[3.875  -3.5   0.0   1.0]

            ans_10.at(1,0)=-1.25;
            ans_10.at(2,0)=4.25;
            ans_10.at(3,0)=3.875;
            ans_10.at(2,1)=-4;
            ans_10.at(3,1)=-3.5;


            SNmultiGaussian<double,4> ans_210;
            ans_210.setLastColumn(2);

//[                1.0                 0.0           0.0    0.0]
//[              -1.25                 1.0           0.0    0.0]
//[               4.25                -4.0           1.0    0.0]
//[-13.606652500000001  12.953320000000001      -4.11333    1.0]

            ans_210.at(1,0)=-1.25;
            ans_210.at(2,0)=4.25;
            ans_210.at(3,0)=-13.6066525;
            ans_210.at(2,1)=-4;
            ans_210.at(3,1)=12.95332;
            ans_210.at(3,2)=-4.11333;

            CPPUNIT_ASSERT(t21.isNumericallyEqual(ans_21,epsilon));
            CPPUNIT_ASSERT(t10.isNumericallyEqual(ans_10,epsilon));
            CPPUNIT_ASSERT(t21_0.isNumericallyEqual(ans_210,epsilon));
            CPPUNIT_ASSERT(t2_10.isNumericallyEqual(ans_210,epsilon));
            CPPUNIT_ASSERT(t210.isNumericallyEqual(ans_210,epsilon));

            auto u_0E=G0*E;


            auto a=G1*(G0*E);
            auto b=(G1*G0)*E;
            SNmatrix<double,4> ans_0E;
            SNmatrix<double,4> ans_10E;

    

//[   4.0    6.0    8.0    9.0]
//[   0.0   -6.5   -3.0 -10.25]
//[   0.0   -2.5   -3.0  -2.75]
//[   0.0    2.0    2.0    2.5]
ans_0E.at(0,0)=4; ans_0E.at(0,1)=6; ans_0E.at(0,2)=8; ans_0E.at(0,3)=9;
ans_0E.at(1,0)=0; ans_0E.at(1,1)=-6.5; ans_0E.at(1,2)=-3; ans_0E.at(1,3)=-10.25;
ans_0E.at(2,0)=0; ans_0E.at(2,1)=-2.5; ans_0E.at(2,2)=-3; ans_0E.at(2,3)=-2.75;
ans_0E.at(3,0)=0; ans_0E.at(3,1)=2; ans_0E.at(3,2)=2; ans_0E.at(3,3)=2.5;


//[   4.0    6.0    8.0    9.0]
//[   0.0   -6.5   -3.0 -10.25]
//[   0.0   23.5    9.0  38.25]
//[   0.0  24.75   12.5 38.375]

ans_10E.at(0,0)=4; ans_10E.at(0,1)=6; ans_10E.at(0,2)=8; ans_10E.at(0,3)=9;
ans_10E.at(1,0)=0; ans_10E.at(1,1)=-6.5; ans_10E.at(1,2)=-3; ans_10E.at(1,3)=-10.25;
ans_10E.at(2,0)=0; ans_10E.at(2,1)=23.5; ans_10E.at(2,2)=9; ans_10E.at(2,3)=38.25;
ans_10E.at(3,0)=0; ans_10E.at(3,1)=24.75; ans_10E.at(3,2)=12.5; ans_10E.at(3,3)=38.375;

            CPPUNIT_ASSERT(u_0E.isNumericallyEqual(ans_0E,epsilon));
            CPPUNIT_ASSERT(a.isNumericallyEqual(ans_10E,epsilon));
            CPPUNIT_ASSERT(b.isNumericallyEqual(ans_10E,epsilon));

            // This is a test for the initialization of
            // data_last_column

            SNmultiGaussian<double,4> mag;
            mag=G2*(G1*G0);
            CPPUNIT_ASSERT(mag.isNumericallyEqual(t210,epsilon));
        }
        void multi_working_tests()
        {
            echo_function_test("multi working tests");

            auto E=testMatrixE();       //SNmatrix<double,4>

            SNgaussian<double,4> G0(E);

            auto G1=(G0*E).getGaussian(1);
            auto G2=(G1*G0*E).getGaussian(2);

            SNmultiGaussian<double,4> mg=G2*(G1*G0);
            auto prod=mg*E;

            echo_single_test("check vanishing components in G2*G1*G0*E");
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
            echo_function_test("inverse_tests");

            auto J=testMatrixJ();       //multigaussian
            auto iJ=J.inverse();
            auto ID=SNidentity<double,4>();
            double epsilon(0.0001);

            SNmultiGaussian<double,4> ans;
            ans.setLastColumn(1);
            ans.at(1,0)=-2; 
            ans.at(2,0)=5; ans.at(2,1)=-4;
            ans.at(3,0)=8; ans.at(3,1)=-7;

            echo_single_test("the inverse is correctly computed");
            CPPUNIT_ASSERT(iJ==ans);

            echo_single_test("J*iJ()");
            CPPUNIT_ASSERT(ID.isNumericallyEqual(J*iJ,epsilon));

            echo_single_test("J*J.inverse()");
            CPPUNIT_ASSERT(ID.isNumericallyEqual(J*J.inverse(),epsilon));

            echo_single_test("iJ*J");
            CPPUNIT_ASSERT(ID.isNumericallyEqual(iJ*J,epsilon));

            echo_single_test("J.inverse()*J");
            CPPUNIT_ASSERT(ID.isNumericallyEqual(J.inverse()*J,epsilon));

            echo_single_test("The inverse of the inverse");
            CPPUNIT_ASSERT(J.isNumericallyEqual(iJ.inverse(),epsilon));
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
        void swap_lines_tests()
        {
            echo_function_test("swap_lines_tests");
            SNmultiGaussian<double,5> K=testMatrixK();

            SNmultiGaussian<double,5> ans;
            ans.setLastColumn(2);
            ans.at(1,0)=2; 
            ans.at(2,0)=3; ans.at(2,1)=4;
            ans.at(4,0)=6; ans.at(4,1)=7; ans.at(4,2)=8;
            ans.at(3,0)=9; ans.at(3,1)=2; ans.at(3,2)=3;


            echo_single_test("basic swap");
            K.swapLines(3,4);
            CPPUNIT_ASSERT(K==ans);
    
            echo_single_test("throw when swapping wrong lines");
            K.swapLines(3,4);
            CPPUNIT_ASSERT_THROW(K.swapLines(4,2),ProbablyNotWhatYouWantException);
        }

        /** 
         * There is an implicit conversion from `SNmultiGaussian` to
         * `SNlowerTriangular`.
         * */
        void construct_lower_triangular()
        {
            SNmultiGaussian<double,5> M;
            SNlowerTriangular<double,5> L(M);
            CPPUNIT_ASSERT(M==L);
        }
    public:
        void runTest()
        {
            construct_lower_triangular();
            swap_lines_tests();
            inverse_tests();
            product_tests();
            multi_working_tests();
            non_initialized_tests();
            wrong_order_works();
            associativity_check();
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
