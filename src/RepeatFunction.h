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


#ifndef __REPEATFUNCTION_H_53825__
#define __REPEATFUNCTION_H_53825__

#include <functional>

class RepeatFunction
/*
 This class takes a function as constructor and then becomes a 
 functor that repeats that function.
*/
{
    private :
        std::function<double(double)> f;
    public :
        RepeatFunction(std::function<double(double)>);
        double operator() (double x) const;
};

#endif
