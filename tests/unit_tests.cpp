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
#include "../src/SNline.h"
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

class SNelementTest : public CppUnit::TestCase
{
    private :
        void test_no_destruct(SNmatrix<int,2>& A)
        // extract an SNelement from the matrix A.
        // When we exit the function, the element is deleted and A must not
        // be deleted in the same time.
        // See test_delete.
        {
            A.at(1,1)=1;
            SNelement<int,2> el=A.getElement(0,0);
            int a=el.getValue();
            CPPUNIT_ASSERT(a==1);
        }
        void test_delete()
        {
            SNmatrix<int,2> A;
            A.at(0,0)=1;
            A.at(0,1)=20;
            A.at(1,0)=3;
            A.at(1,1)=4;
            A.at(0,1)=2;

            test_no_destruct(A);
            SNelement<int,2> el=A.getElement(0,0);      // Should not crash.
            CPPUNIT_ASSERT(el.getValue()==1);

            // check that the matrix used in 'test_no_delete' was not a copy.
            CPPUNIT_ASSERT(A.at(1,1)==1);
        }
        void test_assignation()
        {
            SNmatrix<int,2> A;
            A.at(0,0)=1;
            A.at(0,1)=20;
            A.at(1,0)=3;
            A.at(1,1)=4;
            A.at(0,1)=2;

            SNmatrix<int,2> B;
            B.at(0,0)=2;
            B.at(1,1)=5.3;              // notice the implicit rounding to 5

            SNelement<int,2> el=A.getElement(0,0);
            CPPUNIT_ASSERT(el.getValue()==1);

            el=B.getElement(0,0);                   // re-assignation. Now the matrix of 'el' should be B, no A anymore.

            // New value is ok
            CPPUNIT_ASSERT(el.getValue()==2);
            // No more the values of the old matrix
            CPPUNIT_ASSERT(!(el.getSNmatrix().at(1,1)==4));

            // Correct values of the new matrix
            CPPUNIT_ASSERT(el.getSNmatrix().at(1,1)==5);
            CPPUNIT_ASSERT(el.getSNmatrix().getElement(1,1).getValue()==5);

            // check that A itself is not modified.
            CPPUNIT_ASSERT(!(A.at(1,1)==5));
            CPPUNIT_ASSERT(A.at(1,1)==4);
            CPPUNIT_ASSERT(A.at(0,0)==1);
        };
        void test_chain_matrix_element()
            // The matrix given from SNelement.getSNmatrix is the original one, not a copy.
            // So modifying the matrix returned by one element should change the original matrix.
        {
            SNmatrix<int,2> A;
            A.at(0,0)=1;
            A.at(0,1)=20;

            SNelement<int,2> el=A.getElement(0,0);
            el.getSNmatrix().at(0,1)=12;
            CPPUNIT_ASSERT(A.at(0,1)=12);
        }
    public :
        void runTest()
        {
            test_delete();
            test_assignation();
            test_chain_matrix_element();
        }
};

class SNmatrixTest : public CppUnit::TestCase
{
    private :
        void test_instantiate()
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

            auto el1=A.getElement(0,0);
            CPPUNIT_ASSERT(el1.getValue()==1);

            auto el2=A.getElement(0,1);
            CPPUNIT_ASSERT(el2.getValue()==3/5);

            auto back1=el2.getSNmatrix().at(2,1);
            CPPUNIT_ASSERT(back1==8);

            // test if the matrix is a reference 
            // (thus 'automatically' updates itself)
            A.at(2,2)=12;
            auto back2=el2.getSNmatrix().at(2,2);
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
        void test_max_on_column()
            // test if one get the right maximal element on a column 
            // (for the Gauss's pivot)
        {
            SNmatrix<double,3> A;
            A.at(0,0)=1;
            A.at(0,1)=-3;       // Test negative and first line
            A.at(0,2)=5;
            A.at(1,0)=2;
            A.at(1,1)=1;
            A.at(1,2)=6;
            A.at(2,0)=1;
            A.at(2,1)=2;
            A.at(2,2)=6.1;      // test last line and decimal value

            auto max0=A.getLargerOnColumn(0);
            auto max1=A.getLargerOnColumn(1);
            auto max2=A.getLargerOnColumn(2);
            CPPUNIT_ASSERT(max0.getValue()==2);
            CPPUNIT_ASSERT(max1.getValue()==-3);
            CPPUNIT_ASSERT(max2.getValue()==6.1);


            A.at(2,2)=-2; 
            max0=A.getLargerUnderDiagonal(0);           // By the way, test the assignation operator.
            max1=A.getLargerUnderDiagonal(1);
            max2=A.getLargerUnderDiagonal(2);
            CPPUNIT_ASSERT(max0.getValue()==2);
            CPPUNIT_ASSERT(max1.getValue()==2);
            CPPUNIT_ASSERT(max2.getValue()==-2);
        }
        void test_copy_constructor()
            // test the copy constructor and the equality operator.
        {

            SNmatrix<int,2> A;
            A.at(0,0)=1;
            A.at(0,1)=2;
            A.at(1,0)=3;
            A.at(1,1)=4;

            auto B(A);
            CPPUNIT_ASSERT(A==B);
        }
        void test_swap_line()
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

            auto B(A);

            A.swapLines(0,1);
            CPPUNIT_ASSERT(A.at(0,1)==1);
            A.swapLines(0,1);
            CPPUNIT_ASSERT(A==B);
        }
    public :
        void runTest()
        {
            test_instantiate();
            test_populate();
            test_element_reference();
            test_max_on_column();
            test_copy_constructor();
            test_swap_line();
        }
};

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
            CPPUNIT_ASSERT(L1.at(0)==0);                     // check re-assignation
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
            L2.makeUnit();      //should crash for the moment

            std::array<double,3> d3;
            d3.at(0)=0;
            d3.at(1)=0;
            d3.at(2)=-0.1;
            SNline<double,3> L3(d3);
        }

    public :
        void runTest()
        {
            test_first_nonzero();
            test_make_unit();
        }
};

int main ()
{
    SNmatrixTest sn_test;
    sn_test.runTest();

    SNlineTest sl_test;
    sl_test.runTest();

    SNelementTest se_test;
    se_test.runTest();

    RepeatFunctionTest rf_test;
    rf_test.runTest();
}
