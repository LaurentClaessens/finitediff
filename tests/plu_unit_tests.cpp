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

#include "../src/SNmatrices/SNline.h"
#include "../src/SNplu.h"
#include "../src/SNmatrices/SNmatrix.h"
#include "TestMatrices.cpp"
#include "auto_tests_matrices.h"

#include "../src/Utilities.h"

template <class T,unsigned int tp_size>
void auto_test( const AutoTestMatrix<T,tp_size>& atm  )
{
    echo_function_test("One more matrix");
    auto plu=atm.A.getPLU();
    double epsilon(0.0000001);

    // The variables with 'c' prefix are the ones computed.

    auto cP=plu.getP();
    auto cL=plu.getL();
    auto cU=plu.getU();

    echo_single_test("product ans_PLU");
    auto ans_prod=atm.ans_P*atm.ans_L*atm.ans_U;
    CPPUNIT_ASSERT(ans_prod.isNumericallyEqual(atm.A,epsilon));


    echo_single_test("U factor");
    CPPUNIT_ASSERT(cU.isNumericallyEqual(atm.ans_U,epsilon));
    echo_single_test("P factor");
    CPPUNIT_ASSERT(cP.isNumericallyEqual(atm.ans_P,epsilon));
    echo_single_test("L factor");
    CPPUNIT_ASSERT(cL.isNumericallyEqual(atm.ans_L,epsilon));

    echo_single_test("product cPLU");
    auto c_prod=cP*cL*cU;
    CPPUNIT_ASSERT(c_prod.isNumericallyEqual(atm.A,epsilon));
}

class pluTest : public CppUnit::TestCase
{
    private :
        template <unsigned int s>
        void test_A(  SNmatrix<double,s> A,SNmatrix<double,s> A_U  )
        {
            echo_function_test("Testing a pair A,A_U");

            auto plu=A.getPLU();
            auto cP(plu.getP());
            auto cL(plu.getL());
            auto cU(plu.getU());

            double epsilon(0.0000001);
            echo_single_test("Testing a U");
            CPPUNIT_ASSERT(cU.isNumericallyEqual(A_U,epsilon));
            echo_single_test("Testing the product");

            debug_matrix_print("A",A);
            debug_matrix_print("cP*cL*cU",cP*cL*cU);

            CPPUNIT_ASSERT(A.isNumericallyEqual( cP*cL*cU,epsilon));
        }
        /** 
         * for adding a test here,
         * - add the matrix in `create_tests_matrices.sage`
         * - launch `sage create_tests_matrices` in a terminal
         * - add here a line `auto_test(atm_<chosen name>);`
         * */
        void launch_auto_tests_sage()
        {
            auto_test(atm_ooTBBRooWpZkfJ);
            auto_test(atm_FOO);
            auto_test(atm_BAR);
            auto_test(atm_ooIJXAooDhmylq);
        }
        void plu_from_PLU_tests()
        {
            echo_function_test("PLU from PLU");

            // We take the PLU of A and then try to get the PLU
            // of P,L and U.
            // The result of PLU of P should be (P,id,id).
            // The result of PLU of L should be (id,L,id).
            // The result of PLU of U should be (id,id,U).

            auto A=testMatrixL();
            auto plu=A.getPLU();

            // Not 'auto' because I want to compute the 
            // PLU decomposition from scratch.
            // If 'auto', the first will be a "SNpermutation" object
            // and is eligible to have a very optimized PLU decomposition.
            SNmatrix<double,5> mP(plu.getP());
            SNmatrix<double,5> mL(plu.getL());
            SNmatrix<double,5> mU(plu.getU());

            auto Pplu=mP.getPLU();
            auto Lplu=mL.getPLU();
            auto Uplu=mU.getPLU();

            SNmatrix<double,5> PP=Pplu.getP();
            SNmatrix<double,5> PL=Pplu.getL();
            SNmatrix<double,5> PU=Pplu.getU();

            SNmatrix<double,5> LP=Lplu.getP();
            SNmatrix<double,5> LL=Lplu.getL();
            SNmatrix<double,5> LU=Lplu.getU();

            SNmatrix<double,5> UP=Uplu.getP();
            SNmatrix<double,5> UL=Uplu.getL();
            SNmatrix<double,5> UU=Uplu.getU();

            auto ID=SNidentity<double,5>();

            CPPUNIT_ASSERT(PP==mP);
            CPPUNIT_ASSERT(PL==ID);
            CPPUNIT_ASSERT(PU==ID);

            CPPUNIT_ASSERT(LP==ID);
            CPPUNIT_ASSERT(LL==mL);
            CPPUNIT_ASSERT(LU==ID);

            CPPUNIT_ASSERT(UP==ID);
            CPPUNIT_ASSERT(UL==ID);
            CPPUNIT_ASSERT(UU==mU);
        }
    public:
        void runTest()
        {
            launch_auto_tests_sage();
            test_A(testMatrixA(),testMatrixA_U());
            test_A(testMatrixB(),testMatrixB_U());
            test_A(testMatrixD(),testMatrixD_U());
            test_A(testMatrixE(),testMatrixE_U());
            plu_from_PLU_tests();
        }
};

int main ()
{
    std::cout<<"pluTest"<<std::endl;
    pluTest plu_test;
    plu_test.runTest();
}
