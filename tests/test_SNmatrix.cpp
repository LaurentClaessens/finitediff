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
#include "../src/SNmatrix.h"

class SNmatrixTest : public CppUnit::TestCase
{
    public :
        void runTest()
        {
            SNmatrix<double,7> sn;
            CPPUNIT_ASSERT(sn.getSize()==7);
        }
};

int main ()
{
    SNmatrixTest sn_test;
    sn_test.runTest();
}