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

/*
 
   This is the base class for the other matrices types.
   A matrix can be
   - diagonal
   - upper/lower triangular
   - a gaussian transformation
   - symmetric
   - etc.

   Each of these types have different implementations for storing the elements, 
   perform equality tests (opertor==), multiplication (operator*) and so on.

   Much of these operation, on the least ressort, depend only on the 'get' method
   that returns an element of the matrix.
 
*/
