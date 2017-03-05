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
#include "../src/RepeatFunction.h"

double square(double x)
// return the square of the given number.
{
    return x*x;
};

class RepeatFunctionTest : public CppUnit::TestCase
{
    public :
        void runTest()
        {
            RepeatFunction rf(square);
            CPPUNIT_ASSERT(rf(4)==16);
            CPPUNIT_ASSERT(rf(0)==0);
            CPPUNIT_ASSERT(rf(-2)==4);
            CPPUNIT_ASSERT(rf(1.5)==2.25);
        }
};


class SNmatrixTest : public CppUnit::TestCase
{
    private :
        void test1()
        {
            SNmatrix<double,7> sn1;
            CPPUNIT_ASSERT(sn1.getSize()==7);
        }
        void test_element_reference()
            // test if the reference to the original matrix in the matrix element works
        {
            SNmatrix<double,3> A;
            A.at(0,0)=1;
            A.at(0,1)=3/5;
            A.at(0,2)=3;
            A.at(1,0)=2;
            A.at(1,1)=5;
            A.at(1,2)=0;
            A.at(2,0)=3;
            A.at(2,1)=8;
            A.at(2,2)=0;

            SNelement<double,3> el1=A.getElement(0,0);
            CPPUNIT_ASSERT(el1.getValue()==1);

            SNelement<double,3> el2=A.getElement(0,1);
            CPPUNIT_ASSERT(el2.getValue()==3/5);

            double back1=el2.getSNmatrix().at(2,1);
            CPPUNIT_ASSERT(back1==8);

            A.at(2,2)=12;
            double back2=el2.getSNmatrix().at(2,2);
            CPPUNIT_ASSERT(back2==12);
        }
        void test_populate()
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
            test_populate();
            test_element_reference();
        }
};

int main ()
{
    SNmatrixTest sn_test;
    sn_test.runTest();

    RepeatFunctionTest rf_test;
    rf_test.runTest();
}
