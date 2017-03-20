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

#include <iostream>
std::ostream& debug_print(std::cout);

class SNpermutationsTest : public CppUnit::TestCase
{
    private :
        void test_permutation()
        {
            echo_function_test("test_permutation");
            debug_print<<"MON TEST ACTUEL"<<std::endl;

            auto A=testMatrixE();
            auto plu=A.getPLU();
            std::cout<<"La permutation :"<<std::endl;
            std::cout<<plu.getPermutation()<<std::endl;
        }
    public:
        void runTest()
        {
            test_permutation();
        }
};

int main ()
{
    std::cout<<"sn_permutation_unit_tests"<<std::endl;
    SNpermutationsTest sn_perm_tests;
    sn_perm_tests.runTest();
}
