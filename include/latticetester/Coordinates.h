// This file is part of LatticeTester.
//
// Copyright (C) 2012-2022  The LatticeTester authors, under the occasional supervision
// of Pierre L'Ecuyer at Universit� de Montr�al.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License a t
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef LATTICETESTER__COORDINATES_H
#define LATTICETESTER__COORDINATES_H

#include <iterator>
#include <set>
#include <map>
#include <iostream>
#include <vector>

using namespace std;

namespace LatticeTester {

  /**
   * \class Coordinates
   *
   * An object type that contains a set of coordinate indices, used to specify a projection.
   *
   * This is a standard C++ set, `std::set<std::size_t>`, with additional input and output operators.
   * These types of objects are created and returned by the classes in the `CoordinateSets` namespace.
   * This class also offers some input and output operators.
   * Important: The coordinates are assumed to start at 1.
   *
   * An object `coords` that contains the coordinates {1, 4, 5} can be created via
   * `Coordinates coords({1, 4, 5});` for example.
   */
  class Coordinates : public std::set<std::size_t> {
    public:

      /**
       * Constructs an empty set of coordinates.
       */
      Coordinates():
        std::set<value_type>()
    { }

      /**
       * Copy-constructor.
       */
      Coordinates(const Coordinates& other):
        std::set<value_type>(other)
    { }

      /**
       * Builds a set of coordinates by adding all the elements of coord in it.
       * */
      Coordinates(const std::vector<std::size_t>& coord) {
        for (auto it = coord.begin(); it != coord.end(); it++) {
           insert(*it);
        }
      }

      /**
       * Constructs a coordinate set populated with the values from `first`
       * (inclusively) to `last` (exclusively).
       */
      template<typename InputIterator>
        Coordinates(InputIterator first, InputIterator last):
          std::set<value_type>(first, last)
       { }
  };

  /**
   * \relates Coordinates
   * Formats the coordinate set `coords` and outputs it to `os`.
   */
  static std::ostream& operator<< (std::ostream& os, const Coordinates& coords)  {
	os << "{";
	Coordinates::const_iterator it = coords.begin();
	if (it != coords.end()) {
	  os << *it;
	  while (++it != coords.end())
	    os << "," << *it;
	}
	os << "}";
	return os;
  }


  /**
   * \relates Coordinates
   * Reads a formatted coordinate set from `is`.
   *
   * The input must consist of positive integers separated by whitespace and/or by
   * commas, and optionally enclosed in braces.  The ordering is not important.
   * Repeated values are ignored.
   * For example, the following strings are valid input that would produce
   * equivalent Coordinates objects:
   * - <tt>1 2 5</tt>
   * - <tt>1, 2, 5</tt>
   * - <tt>{1 2 5}</tt>
   * - <tt>{1,2,5}</tt>
   * - <tt>{1, 2, 5}</tt>
   * - <tt>2 5 1</tt>
   * - <tt>2 1 5 1</tt>
   */
  static std::istream& operator>> (std::istream& is, Coordinates& coords) {
     coords.clear();
  
     string digits = "0123456789";
     string sep = " \t,";
  
     // check if coordinate set is enclosed in braces
     bool with_braces = false;
     if (is.peek() == '{') {
       is.get();
       with_braces = true;
     }
  
     while (true) {
       if (with_braces && is.peek() == '}') {
         is.get();
         break;
       }
       if (digits.find(is.peek()) != string::npos) {
         // digit found
         Coordinates::value_type val;
         is >> val;
         coords.insert(val);
         continue;
       }
       if (sep.find(is.peek()) != string::npos) {
         // discard separator character
         is.get();
         continue;
       }
       break;
     }
     return is;
   }

} //end namespace

#endif
