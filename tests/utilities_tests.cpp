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

#include "../src/Utilities.h"

class UtilitiesTest : public CppUnit::TestCase
{
    private :
        void test_value_length()
        {
            echo_function_test("test_value_length");
            double x1(4);
            double x2(4.34);
            double x3(10./3);
            double x4(-12./7);
            double x5( 0.000000001  );

            CPPUNIT_ASSERT(value_length(x1)==1);
            CPPUNIT_ASSERT(value_length(x2)==4);
            CPPUNIT_ASSERT(value_length(x3)==7);
            CPPUNIT_ASSERT(value_length(x4)==8);
            CPPUNIT_ASSERT(value_length(x5)==5);
        }
    public :
        void runTest()
        {
            test_value_length();
        }
};
int main ()
{
    UtilitiesTest util_test;
    util_test.runTest();
}
