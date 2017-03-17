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

#include "../src/SNmatrices/SNmatrix.h"
#include "../src/SNmatrices/SNpermutation.h"
#include "../src/SNmatrices/SNline.h"
#include "../src/SNmatrices/SNlowerTriangularMatrix.h"
#include "../src/SNmatrices/SNoperators.h"
#include "../src/SNplu.h"

#include <iostream>
std::ostream& debug_print(std::cout);

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
            double epsilon(0.0000001);
            auto A=testMatrixB();
            /*
             1  2  3
             2  5  0
             3  8  0
            */
            auto A_U=testMatrixB_U();
            auto plu_A= A.getPLU();

            // check that 'A' is turned into its 'U'
            CPPUNIT_ASSERT(plu_A.getU()==A);
            A.subtract(A_U); 
            CPPUNIT_ASSERT(A.max_norm()<epsilon);

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

class pluTest : public CppUnit::TestCase
{
    private :
        void test_reference()
        {
            auto A=testMatrixA();
            auto plu=A.getPLU();
            CPPUNIT_ASSERT(A==plu.getU());

            // check that the 'U' matrix in PLU is a reference.
            // this us however bad practice to modify A after having
            // the PLU done.
            A.at(1,1)=0.132;
            CPPUNIT_ASSERT(plu.getU().get(1,1)==0.132);
        };
        template <unsigned int s>
        void test_A(  SNmatrix<double,s> A,SNmatrix<double,s> A_U  )
        {
            auto plu=A.getPLU();

            auto mU(plu.getU());
            mU.subtract(A_U);
            double epsilon(0.0000001);
            CPPUNIT_ASSERT(mU.max_norm()<epsilon);
        }
        void test_permutation()
        {
            debug_print<<"MON TEST ACTUEL";
            debug_print.endl();

            auto A=testMatrixE();
            auto plu=A.getPLU();
            std::cout<<"La permutation :"<<std::endl;
            std::cout<<plu.getPermutation()<<std::endl;
        }
    public:
        void runTest()
        {
            test_reference();
            std::cout<<"pluTest A"<<std::endl;
            test_A(testMatrixA(),testMatrixA_U());
            std::cout<<"pluTest B"<<std::endl;
            test_A(testMatrixB(),testMatrixB_U());
            std::cout<<"pluTest D"<<std::endl;
            test_A(testMatrixD(),testMatrixD_U());
            std::cout<<"pluTest E"<<std::endl;
            test_A(testMatrixE(),testMatrixE_U());
            test_permutation();
        }
};

class MultiplicationTest : public CppUnit::TestCase
{
    private :
        void test_gauss_times_matrix()
        {
            auto F=testMatrixF();
            auto G=F.getGaussian(1);    // the gaussian matrix of F

// computations from Sage :
//  
// sage: F=matrix(  [ [1,0,3,9],[6,2,3,5],[7,8,1,3],[7,7,4,6] ]  )
// sage: G=matrix(   [     [1,0,0,0],[0,1,0,0],[0,-4,1,0],[0,-7/2,0,1]   ]   )
// sage: G*F
//[    1     0     3     9]
//[    6     2     3     5]
//[  -17     0   -11   -17]
//[  -14     0 -13/2 -23/2]

            SNmatrix<double,4> ans_G;
            ans_G.at(0,0)=1; ans_G.at(0,1)=0; ans_G.at(0,2)=0; ans_G.at(0,3)=0;
            ans_G.at(1,0)=0; ans_G.at(1,1)=1; ans_G.at(1,2)=0; ans_G.at(1,3)=0;
            ans_G.at(2,0)=0; ans_G.at(2,1)=4; ans_G.at(2,2)=1; ans_G.at(2,3)=0;
            ans_G.at(3,0)=0; ans_G.at(3,1)=7./2; ans_G.at(3,2)=0; ans_G.at(3,3)=1;

            CPPUNIT_ASSERT(G==ans_G);

            SNmatrix<double,4> GstarF;
            GstarF.at(0,0)=1; GstarF.at(0,1)=0; GstarF.at(0,2)=3; GstarF.at(0,3)=9;
            GstarF.at(1,0)=6; GstarF.at(1,1)=2; GstarF.at(1,2)=3; GstarF.at(1,3)=5;
            GstarF.at(2,0)=-17; GstarF.at(2,1)=0; GstarF.at(2,2)=-11; GstarF.at(2,3)=-17;
            GstarF.at(3,0)=-14; GstarF.at(3,1)=0; GstarF.at(3,2)=-13./2; GstarF.at(3,3)=-23./2;

            CPPUNIT_ASSERT(G*F==GstarF);
        }
        void test_gauss_times_lower_trig()
        {
            auto A=testMatrixG();
            auto G=A.getGaussian(0); 

            SNmatrix<double,4> ans_G;
            ans_G.at(0,0)=1; ans_G.at(0,1)=0; ans_G.at(0,2)=0; ans_G.at(0,3)=0;
            ans_G.at(1,0)=-1./3; ans_G.at(1,1)=1; ans_G.at(1,2)=0; ans_G.at(1,3)=0;
            ans_G.at(2,0)=-1; ans_G.at(2,1)=0; ans_G.at(2,2)=1; ans_G.at(2,3)=0;
            ans_G.at(3,0)=-2; ans_G.at(3,1)=0; ans_G.at(3,2)=0; ans_G.at(3,3)=1;

            CPPUNIT_ASSERT(G==ans_G);

            SNmatrix<double,4> GstarA;
//[3 0 0 0]
//[0 2 0 0]
//[0 4 5 0]
//[0 7 8 1]
            GstarA.at(0,0)=3; GstarA.at(0,1)=0; GstarA.at(0,2)=0; GstarA.at(0,3)=0;
            GstarA.at(1,0)=0; GstarA.at(1,1)=2; GstarA.at(1,2)=0; GstarA.at(1,3)=0;
            GstarA.at(2,0)=0; GstarA.at(2,1)=4; GstarA.at(2,2)=5; GstarA.at(2,3)=0;
            GstarA.at(3,0)=0; GstarA.at(3,1)=7; GstarA.at(3,2)=8; GstarA.at(3,3)=1;

            CPPUNIT_ASSERT(G*A==GstarA);
        }

    public :
        void runTest()
        {
            test_gauss_times_matrix();
            test_gauss_times_lower_trig();
        }
};

int main ()
{



    std::cout<<"SNelementTest"<<std::endl;
    SNelementTest se_test;
    se_test.runTest();

    std::cout<<"GaussTest"<<std::endl;
    GaussTest gauss_test;
    gauss_test.runTest();

    std::cout<<"pluTest"<<std::endl;
    pluTest plu_test;
    plu_test.runTest();
}
