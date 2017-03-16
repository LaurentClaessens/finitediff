
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
#include "../src/SNexceptions.cpp"
#include "TestMatrices.cpp"

#include "../src/DebugPrint.h"
DebugPrint debug_print;

class ExceptionsTests : public CppUnit::TestCase
{
    private:
        void out_of_range_test()
        {
            SNmatrix<int,4> A;
            debug_print<<"test numéro 1";
            CPPUNIT_ASSERT_THROW(A.get(5,1),SNoutOfRangeException);
            SNlowerTriangularMatrix<int,4> B;
            debug_print<<"test numéro 2";
            CPPUNIT_ASSERT_THROW(B.at(1,5),SNoutOfRangeException);
        }
        void change_not_allowed_test()
        {
            SNlowerTriangularMatrix<int,4> A;
            debug_print<<"test numéro 3";
            CPPUNIT_ASSERT_THROW(A.at(1,2),SNchangeNotAllowedException);
        }

        void incompatible_matrix_size_test()
        {
            SNlowerTriangularMatrix<int,4> A;
            SNmatrix<double,2> B;
            debug_print<<"test numéro 4";
            debug_print.endl();
            CPPUNIT_ASSERT_THROW(A==B,IncompatibleMatrixSizeException);

            debug_print<<"test numéro 5";
            debug_print.endl();
            SNmatrix<double,2> C;
            SNmatrix<double,3> D;
            debug_print<<"test numéro 6";
            debug_print.endl();
            CPPUNIT_ASSERT_THROW(auto K=C*D,IncompatibleMatrixSizeException);
        }

    public:
        void runTest()
        {
            out_of_range_test();
            change_not_allowed_test();
            incompatible_matrix_size_test();
        }
};

int main ()
{
    std::cout<<"Exceptions tests"<<std::endl;
    ExceptionsTests ex_test;
    ex_test.runTest();
}
