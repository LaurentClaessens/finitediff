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

#include "Utilities.h"

// cppcheck-suppress unusedFunction
void echo_function_test(const std::string& s)
{
    std::cout<<"   function : "<<s<<std::endl;
}

// cppcheck-suppress unusedFunction
void echo_single_test(const std::string& s)
{
    std::cout<<"      single test : "<<s<<std::endl;
}

std::ostream& debug_print(std::cout);


// cppcheck-suppress unusedFunction
void tooGenericWarning(const std::string& message)
{
    std::cout<<message<<std::endl;
}
