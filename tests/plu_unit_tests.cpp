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

#include "../src/DebugPrint.h"

template <class T,unsigned int tp_size>
void auto_test(const SNmatrix<T,tp_size>& A, const SNpermutation<T,tp_size>& ans_P,const SNlowerTriangular<T,tp_size>& ans_L, SNupperTriangular<T,tp_size>& ans_U)
{
    echo_function_test("One more matrix");
    auto plu=A.getPLU();
    double epsilon(0.0000001);

    auto cP=plu.getP();
    auto cL=plu.getL();
    auto cU=plu.getU();

    // TODO : check that the product PLU is equal to A.

    echo_single_test("P factor");
    CPPUNIT_ASSERT(cP.isNumericallyEqual(ans_P,epsilon));
    echo_single_test("L factor");
    //CPPUNIT_ASSERT(cL.isNumericallyEqual(ans_L,epsilon));
    echo_single_test("U factor");
    CPPUNIT_ASSERT(cU.isNumericallyEqual(ans_U,epsilon));
}

class pluTest : public CppUnit::TestCase
{
    private :
        template <unsigned int s>
        void test_A(  SNmatrix<double,s> A,SNmatrix<double,s> A_U  )
        {
            auto plu=A.getPLU();
            auto mU(plu.getU());

            double epsilon(0.0000001);
            CPPUNIT_ASSERT(mU.isNumericallyEqual(A_U,epsilon));
        }
        void launch_auto_tests_sage()
        {
            auto A=testMatrix_FOO_A();
            auto ans_P=testMatrix_FOO_A_P();
            auto ans_L=testMatrix_FOO_A_L();
            auto ans_U=testMatrix_FOO_A_U();
            auto_test(A,ans_P,ans_L,ans_U);
        }
    public:
        void runTest()
        {
            test_A(testMatrixA(),testMatrixA_U());
            test_A(testMatrixB(),testMatrixB_U());
            test_A(testMatrixD(),testMatrixD_U());
            test_A(testMatrixE(),testMatrixE_U());
            launch_auto_tests_sage();
        }
};

int main ()
{
    std::cout<<"pluTest"<<std::endl;
    pluTest plu_test;
    plu_test.runTest();
}
