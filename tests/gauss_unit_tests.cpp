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

#include "../src/DebugPrint.h"

class GaussTest : public CppUnit::TestCase
{
    private :
        void test_elimination_line()
        {
        SNmatrix<double,3> A;
        A.at(0,0)=1;
        A.at(0,1)=-3;      
        A.at(0,2)=5;
        A.at(1,0)=2;
        A.at(1,1)=1;
        A.at(1,2)=6;
        A.at(2,0)=3;
        A.at(2,1)=2;
        A.at(2,2)=6.1;  

        auto lu0=A.gaussEliminationLine(0);
        auto lu1=A.gaussEliminationLine(1);
        auto lu2=A.gaussEliminationLine(2);

        CPPUNIT_ASSERT(lu0.get(0)==1);
        CPPUNIT_ASSERT(lu0.get(1)==-3);
        CPPUNIT_ASSERT(lu0.get(2)==5);
        CPPUNIT_ASSERT(lu1.get(0)==1);
        CPPUNIT_ASSERT(lu1.get(1)==0.5);
        CPPUNIT_ASSERT(lu1.get(2)==3);
        CPPUNIT_ASSERT(lu2.get(0)==1);
        CPPUNIT_ASSERT(lu2.get(1)==2.0/3);      // notice that 2/3 does not work
        CPPUNIT_ASSERT(lu2.get(2)==6.1/3);

        SNmatrix<double,3> B;
        B.at(1,0)=0;
        B.at(1,1)=0;
        B.at(1,2)=2;
        B.at(2,0)=0;
        B.at(2,1)=3;
        B.at(2,2)=6;  

        auto lt1=B.gaussEliminationLine(1);
        auto lt2=B.gaussEliminationLine(2);

        CPPUNIT_ASSERT(lt1.get(0)==0);
        CPPUNIT_ASSERT(lt1.get(1)==0);
        CPPUNIT_ASSERT(lt1.get(2)==1);
        CPPUNIT_ASSERT(lt2.get(0)==0);
        CPPUNIT_ASSERT(lt2.get(1)==1);
        CPPUNIT_ASSERT(lt2.get(2)==2);
        }
        void test_LminusL()
            // test the matrix manipulation
            // L_i -> L_i- m*L_k
            // that "eliminates" the line i using the line k.
        {
            auto A=testMatrixA();
                /*
               0  3  6
               0  0  2
               0  3  6
                */
            auto l2=A.getSNline(2);
            A.lineMinusLine(0,l2);
            CPPUNIT_ASSERT(A.get(0,0)==0);
            CPPUNIT_ASSERT(A.get(0,1)==0);
            CPPUNIT_ASSERT(A.get(0,2)==0);

            A.lineMinusLine(1,l2);
            CPPUNIT_ASSERT(A.get(1,0)==0);
            CPPUNIT_ASSERT(A.get(1,1)==-3);
            CPPUNIT_ASSERT(A.get(1,2)==-4);

            auto B(A);
            A.lineMinusLine(2,A.getSNline(0));
            CPPUNIT_ASSERT( A.getSNline(2)==B.getSNline(2) );
        }
        void test_upper_triangular()
        {
            echo_function_test("test_upper_triangular");

            debug_print<<"  ----------- LE TETST "<<std::endl;
            debug_print<<""<<std::endl;

            double epsilon(0.0000001);
            auto A=testMatrixB();
            /*
             1  2  3
             2  5  0
             3  8  0
            */
            auto A_U=testMatrixB_U();
            auto plu_A= A.getPLU();

            echo_single_test("plu.getU==A_U");

            debug_print<<"Calculée :"<<std::endl;
            debug_print<<plu_A.getU()<<std::endl;
            debug_print<<"réponse :"<<std::endl;
            debug_print<<A_U<<std::endl;

            CPPUNIT_ASSERT(A_U.isNumericallyEqual(plu_A.getU(),epsilon));

            auto B=testMatrixA();
                /*
                 0  3  6
                 0  0  2
                 0  3  6
                */
            auto B_U=testMatrixA_U();

            auto plu_B=B.getPLU();  
            plu_B.getU().subtract(B_U);
            CPPUNIT_ASSERT(B.max_norm()<epsilon);
        }
    public :
        void runTest()
        {
            test_elimination_line();
            test_LminusL();
            test_upper_triangular();
        }
};

int main ()
{
    std::cout<<"GaussTest"<<std::endl;
    GaussTest gauss_test;
    gauss_test.runTest();
}
