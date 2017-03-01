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
#include "../src/SNmatrix.h"

class SNmatrixTest : public CppUnit::TestCase
{
    private :
        void test1()
        {
            SNmatrix<double,7> sn1;
            CPPUNIT_ASSERT(sn1.getSize()==7);
        }
        void test3()
        {
            SNmatrix<double,3> A;
            A.at(0,0)=1;
            A.at(0,1)=2;
            A.at(0,2)=3;
            A.at(1,0)=2;
            A.at(1,1)=5;
            A.at(1,2)=0;
            A.at(2,0)=3;
            A.at(2,1)=8;
            A.at(2,2)=0;
        }
        void test2()
        {
            SNmatrix<int,2> sn2;
            sn2.at(0,0)=1;
            sn2.at(0,1)=20;
            sn2.at(1,0)=3;
            sn2.at(1,1)=4;
            sn2.at(0,1)=2;
            CPPUNIT_ASSERT(sn2.getSize()==2);
            CPPUNIT_ASSERT(sn2.at(0,0)==1);
            CPPUNIT_ASSERT(sn2.at(0,1)==2);
        }
    public :
        void runTest()
        {
            test1();
            test2();
            test3();
        }
};

int main ()
{
    SNmatrixTest sn_test;
    sn_test.runTest();
}
