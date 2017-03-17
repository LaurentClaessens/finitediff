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

class SNlineTest : public CppUnit::TestCase
{
    private :

        void test_first_nonzero()
        {
            std::array<double,4> d1;
            d1.at(0)=1;
            d1.at(1)=2;
            d1.at(2)=3;
            d1.at(3)=4;
            SNline<double,4> L1(d1);

            CPPUNIT_ASSERT(L1.firstNonZeroColumn()==0);      // check the case in which the first is non zero.

            L1.at(0)=0;
            CPPUNIT_ASSERT(L1.get(0)==0);                     // check re-assignation
            CPPUNIT_ASSERT(L1.firstNonZeroColumn()==1);

            std::array<double,3> d2;
            d2.at(0)=0;
            d2.at(1)=0;
            d2.at(2)=0;
            SNline<double,3> L2(d2);
            CPPUNIT_ASSERT(L2.firstNonZeroColumn()==4);      // check when they are all zero.

            std::array<double,3> d3;
            d3.at(0)=0;
            d3.at(1)=0;
            d3.at(2)=-0.1;
            SNline<double,3> L3(d3);
            CPPUNIT_ASSERT(L3.firstNonZeroColumn()==2);      // check the case in which the first non-zero is the last one.
        }
        void test_make_unit()
        {
            std::array<double,4> d1;
            d1.at(0)=1;
            d1.at(1)=2;
            d1.at(2)=3;
            d1.at(3)=4;

            SNline<double,4> L1(d1);
            SNline<double,4> A1(L1);        // check copy constructor
            L1.makeUnit();
            CPPUNIT_ASSERT(L1==A1);         // check equality operator

            L1.at(0)=0;
            CPPUNIT_ASSERT(L1.get(0)==0);
            CPPUNIT_ASSERT(L1.at(0)==0);
            std::array<double,4> a2;
            a2.at(0)=0;
            a2.at(1)=1;
            a2.at(2)=1.5;
            a2.at(3)=2;
            SNline<double,4> A2(a2);

            L1.makeUnit();
            CPPUNIT_ASSERT(L1==A2);
            L1.makeUnit();
            CPPUNIT_ASSERT(L1==A2);     // twice makeUnit.

            std::array<double,3> d2;
            d2.at(0)=0;
            d2.at(1)=0;
            d2.at(2)=0;
            SNline<double,3> L2(d2);
            SNline<double,3> A3(L2);
            L2.makeUnit();   
            CPPUNIT_ASSERT(L2==A3);

            std::array<double,3> d3;
            d3.at(0)=0;
            d3.at(1)=0;
            d3.at(2)=-0.1;
            SNline<double,3> L3(d3);
            L3.makeUnit();
            CPPUNIT_ASSERT(L3.get(0)==0);
            CPPUNIT_ASSERT(L3.get(1)==0);
            CPPUNIT_ASSERT(L3.get(2)==1);
        }

    void test_get_snline()
        // check the function 'getSNline'
    {
        SNmatrix<double,3> A;
        A.at(0,0)=1;
        A.at(0,1)=-3;      
        A.at(0,2)=5;
        A.at(1,0)=2;
        A.at(1,1)=1;
        A.at(1,2)=6;
        A.at(2,0)=1;
        A.at(2,1)=2;
        A.at(2,2)=6.1;  

        auto l0=A.getSNline(0);
        auto l1=A.getSNline(1);
        auto l2=A.getSNline(2);
        CPPUNIT_ASSERT(l0.get(0)==1);
        CPPUNIT_ASSERT(l0.get(1)==-3);
        CPPUNIT_ASSERT(l0.get(2)==5);
        CPPUNIT_ASSERT(l1.get(1)==1);
        CPPUNIT_ASSERT(l1.get(2)==6);
        CPPUNIT_ASSERT(l2.get(0)==1);
        CPPUNIT_ASSERT(l2.get(2)==6.1);
    }

    public :
        void runTest()
        {
            test_first_nonzero();
            test_make_unit();
            test_get_snline();
        }
};

int main ()
{
    std::cout<<"SNLineTest"<<std::endl;
    SNlineTest sl_test;
    sl_test.runTest();
}