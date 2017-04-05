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

#include <array>

#include <cppunit/TestCase.h>
#include <cppunit/extensions/TypeInfoHelper.h>
#include <cppunit/TestAssert.h>

#include "../src/SNmatrices/SNgaussian.h"

#include "../src/SNmatrices/SNmatrix.h"
#include "../src/SNmatrices/SNpermutation.h"
#include "TestMatrices.cpp"

class SNgaussianTest : public CppUnit::TestCase
{
    private :
        void test_inverse()
        {
            echo_function_test("test_inverse");

            auto H=testMatrixH();
            auto Hg=H.getGaussian(2);
            auto iHg=Hg.inverse();
            auto ID=SNpermutation<double,4>();

            double epsilon(0.000001);
            CPPUNIT_ASSERT(ID.isNumericallyEqual(Hg*iHg,epsilon)  );
            CPPUNIT_ASSERT(ID.isNumericallyEqual(iHg*Hg,epsilon)  );

        }
    void test_assignation_to_lower_triangular()
    {
        echo_function_test("test_assignation_to_lower_triangular");
        SNlowerTriangular<double,4> L;
        auto H=testMatrixH();
        L=H.getGaussian(0);

        SNlowerTriangular<double,4> M(H.getGaussian(1));
    }
    public:
        void runTest()
        {
            test_inverse();
            test_assignation_to_lower_triangular();
        }
};

int main ()
{
    std::cout<<"sn_gaussian_unit_tests"<<std::endl;
    SNgaussianTest sn_gauss_tests;
    sn_gauss_tests.runTest();
}

