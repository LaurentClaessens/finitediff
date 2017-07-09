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
#include "TestMatrices.cpp"
#include "../src/Utilities.h"

class mNumTest : public CppUnit::TestCase
{
    private:
        void increment_tests()
        {
            echo_function_test("increment_tests");
            unsigned int k(4);
            m_num m(3);
            ++m;
            CPPUNIT_ASSERT(m==4);
            CPPUNIT_ASSERT(m==k);

            // cppcheck-suppress postfixOperator
            m++;
            CPPUNIT_ASSERT(m==5);
            CPPUNIT_ASSERT(m==k+1);

            m_num n=++m;
            CPPUNIT_ASSERT(n==6);
            CPPUNIT_ASSERT(m==6);

            m_num s=m++;
            CPPUNIT_ASSERT(s==6);
            CPPUNIT_ASSERT(m==7);
        }
        void loop_tests()
        {
            echo_function_test("loop_tests");
            std::array<int,4> v;
            for (m_num k=0;k<4;++k)
            {
                v.at(k)=k;
            }
            std::array<int,4> ans_v;
            ans_v.at(0)=0;
            ans_v.at(1)=1;
            ans_v.at(2)=2;
            ans_v.at(3)=3;
            CPPUNIT_ASSERT(v==ans_v);
            std::array<int,4> w;
            for (m_num k=0;k<4;++k)
            {
                w.at(k)=k;
            }
            CPPUNIT_ASSERT(w==ans_v);
        }
    public:
        void runTest()
        {
            increment_tests();
            loop_tests();
        }
};
 
int main ()
{
    std::cout<<"m_num_unit_test"<<std::endl;
    mNumTest m_num_test;
    m_num_test.runTest();
}
