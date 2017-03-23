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

#include "../src/DebugPrint.h"

class SNpermutationsTest : public CppUnit::TestCase
{
    private :
        void test_permutation()
        {
            echo_function_test("test_permutation");



            echo_single_test("Test the identity permutation");
            std::array<unsigned int, 4> a2{ {0,1, 2, 3} };
            Mpermutation<4> permID(a2);
            CPPUNIT_ASSERT(permID(0)==0);
            CPPUNIT_ASSERT(permID(1)==1);
            CPPUNIT_ASSERT(permID(2)==2);
            CPPUNIT_ASSERT(permID(3)==3);

            echo_single_test("Test a constructor with out of range arguments");
            std::array<unsigned int, 4> a1{ {1, 2, 3,4} };
            CPPUNIT_ASSERT_THROW(Mpermutation<4> perm1(a1),PermutationIdexoutOfRangeException);

            echo_single_test("Test throwing when asking a too large number (>tp_size)");
            CPPUNIT_ASSERT_THROW(std::cout<<permID(12),PermutationIdexoutOfRangeException);

            echo_single_test("Test an arbitrary permutation");
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

            echo_single_test("perm3 * perm4");
            std::array<unsigned int, 4> aa{ {2,0, 1, 3} };
            Mpermutation<4> ans_p1(aa);
            CPPUNIT_ASSERT(ans_p1==perm3*perm4);
            echo_single_test("This one is commutative");
            CPPUNIT_ASSERT(perm4*perm3==perm3*perm4);

            std::array<unsigned int, 4> a5{ {3,2, 1, 0} };
            Mpermutation<4> perm5(a5);

            echo_single_test("perm3 * perm5");
            std::array<unsigned int, 4> ab{ {3,0, 2, 1} };
            Mpermutation<4> ans_p2(ab);
            CPPUNIT_ASSERT(ans_p2==perm3*perm5);

            echo_single_test("This one is not commutative : perm5 * perm3");
            std::array<unsigned int, 4> ac{ {2,1,3,0} };
            Mpermutation<4> ans_p3(ac);
            CPPUNIT_ASSERT(ans_p3==perm5*perm3);
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
