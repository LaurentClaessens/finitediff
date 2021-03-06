
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
#include "../src/exceptions/SNexceptions.cpp"
#include "TestMatrices.cpp"

class ExceptionsTests : public CppUnit::TestCase
{
    private:
        void out_of_range_test()
        {
            SNmatrix<int,4> A;
            CPPUNIT_ASSERT_THROW(A.get(5,1),SNoutOfRangeException);
            SNlowerTriangular<int,4> B;
            CPPUNIT_ASSERT_THROW(B.at(1,5),SNoutOfRangeException);
        }
        void change_not_allowed_test()
        {
            SNlowerTriangular<int,4> A;
            CPPUNIT_ASSERT_THROW(A.at(1,2),SNchangeNotAllowedException);
        }

        void incompatible_matrix_size_test()
        {
            SNlowerTriangular<int,4> A;
            SNmatrix<double,2> B;

            // std::cout<<  is in order to avoid the warning about 
            // "equality comparison result unused"
            // I do no know how to do better ...
            CPPUNIT_ASSERT_THROW(std::cout<<(A==B),IncompatibleMatrixSizeException);

            SNmatrix<double,2> C;
            SNmatrix<double,3> D;
            CPPUNIT_ASSERT_THROW(C*D,IncompatibleMatrixSizeException);
        }

        void error_message_test()
        {
            echo_function_test("error_message_test");
            SNmatrix<double,2> C;
            SNmatrix<double,3> D;
            try
            {
                std::cout<<C*D;
            }
            catch (IncompatibleMatrixSizeException& e)
            {
                std::string ans("First matrix has size 2 while second matrix has size 3");
                CPPUNIT_ASSERT(e.what()==ans);
            }
        };

    public:
        void runTest()
        {
            error_message_test();
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
