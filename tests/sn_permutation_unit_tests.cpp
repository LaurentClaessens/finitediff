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

#include "../src/SNmatrices/SNline.h"
#include "../src/SNplu.h"
#include "../src/SNmatrices/SNmatrix.h"
#include "TestMatrices.cpp"

#include <iostream>
std::ostream& debug_print(std::cout);

class SNpermutationsTest : public CppUnit::TestCase
{
    private :
        void test_permutation()
        {
            echo_function_test("test_permutation");

            std::array<unsigned int, 4> a1{ {1, 2, 3,4} };
            Mpermutation<4> perm1(a1);         // there should be some errors here


            std::array<unsigned int, 4> a2{ {0,1, 2, 3} };
            Mpermutation<4> permID(a2);
            CPPUNIT_ASSERT(permID(0)==0);
            CPPUNIT_ASSERT(permID(1)==1);
            CPPUNIT_ASSERT(permID(2)==2);
            CPPUNIT_ASSERT(permID(3)==3);

            CPPUNIT_ASSERT(permID(12)==3);      // should do something about that

            std::array<unsigned int, 4> a3{ {1,2, 0, 3} };
            Mpermutation<4> perm2(a3);
            CPPUNIT_ASSERT(perm2(0)==1);
            CPPUNIT_ASSERT(perm2(1)==2);
            CPPUNIT_ASSERT(perm2(2)==0);
            CPPUNIT_ASSERT(perm2(3)==3);
         }
        void test_identity()
        {
            echo_function_test("test_identiry");

            std::array<unsigned int, 4> aID{ {0,1, 2, 3} };
            Mpermutation<4> permID(aID);

            std::array<unsigned int, 4> a3{ {1,2, 0, 3} };
            Mpermutation<4> perm3(a3);

            std::array<unsigned int, 4> a4{ {1,2, 0, 3} };
            Mpermutation<4> perm4(a4);

            std::array<unsigned int, 4> aa{ {1,2, 0, 3} };
            Mpermutation<4> ans_p1(aa);

            std::array<unsigned int, 4> a5{ {3,2, 1, 0} };
            Mpermutation<4> perm5(a5);

            std::array<unsigned int, 4> ab{ {3,0, 2, 1} };
            Mpermutation<4> ans_p2(ab);

            CPPUNIT_ASSERT(perm3==permID*perm3);
            CPPUNIT_ASSERT(perm4==permID*perm4);
            CPPUNIT_ASSERT(perm5==permID*perm5);
            CPPUNIT_ASSERT(ans_p1==ans_p1*permID);
            CPPUNIT_ASSERT(ans_p2==ans_p2*permID);
        }
        void test_product()
        {
            echo_function_test("test_product");
            std::array<unsigned int, 4> a3{ {1,2, 0, 3} };
            Mpermutation<4> perm3(a3);

            std::array<unsigned int, 4> a4{ {1,2, 0, 3} };
            Mpermutation<4> perm4(a4);

            std::array<unsigned int, 4> aa{ {1,2, 0, 3} };
            Mpermutation<4> ans_p1(aa);
            CPPUNIT_ASSERT(ans_p1==perm3*perm4);

            std::array<unsigned int, 4> a5{ {3,2, 1, 0} };
            Mpermutation<4> perm5(a5);

            std::array<unsigned int, 4> ab{ {3,0, 2, 1} };
            Mpermutation<4> ans_p2(ab);
            CPPUNIT_ASSERT(ans_p2==perm3*perm5);
        }
    public:
        void runTest()
        {
            test_permutation();
            test_identity();
            test_product();
        }
};

int main ()
{
    std::cout<<"sn_permutation_unit_tests"<<std::endl;
    SNpermutationsTest sn_perm_tests;
    sn_perm_tests.runTest();
}
