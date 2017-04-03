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

#include "../src/SNmatrices/m_num.h"
#include "../src/DebugPrint.h"

class mNumTest : public CppUnit::TestCase
{
    private:
        void increment_tests()
        {
            unsigned int k(4);
            m_num m(3);
            m++;
            CPPUNIT_ASSERT(m==4);
            CPPUNIT_ASSERT(m==k);
            ++m;
            CPPUNIT_ASSERT(m==5);
            CPPUNIT_ASSERT(m==k+1);

            m_num n=++m;
            CPPUNIT_ASSERT(n==6);
            CPPUNIT_ASSERT(m==6);

            m_num s=m++;
            CPPUNIT_ASSERT(s==6);
            CPPUNIT_ASSERT(m==7);
        }
    public:
        void runTest()
        {
            increment_tests();
        }
};
 
int main ()
{
    std::cout<<"m_num_unit_test"<<std::endl;
    mNumTest m_num_test;
    m_num_test.runTest();
}
