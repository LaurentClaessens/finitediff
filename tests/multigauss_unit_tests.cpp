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
#include "TestMatrices.cpp"

class multigaussTests : public CppUnit::TestCase
{
    private :
        void constructor_test()
        {
            auto E=testMatrixE();
            double epsilon(0.0000001);

            echo_single_test("Constructor from matrix");
            
            SNmultiGaussian mg(E);
            CPPUNIT_ASSERT(mg.isNumericallyEqual(E.getGaussian(0),epsilon));

            auto M0=E.getGaussian(0);
            auto M1=E.getGaussian(1);
            auto M2=E.getGaussian(2);

            mg*=M1;
        }
    public:
        void runTest()
        {
            constructor_test();
        }
};

int main ()
{
    std::cout<<"multigaussTests"<<std::endl;
    multigaussTests mg_test;
    mg_test.runTest();
}
