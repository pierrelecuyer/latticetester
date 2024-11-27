// This file is part of LatticeTester.
//
// Copyright (C) 2012-2022  The LatticeTester authors, under the supervision
// of Pierre L'Ecuyer at Universit� de Montr�al.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef LATTICETESTER__FLEXTYPES_H
#define LATTICETESTER__FLEXTYPES_H

#include <cstdint>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/xdouble.h>
#include <NTL/quad_float.h>
// #include "latticetester/NTLWrap.h"  // This one is needed for the vector and matrix types.

/**
 * There are five admissible combinations of types for `(Int, Real)`.
 * They are represented by the five codes given below.
 * For example, to use Int = NTL::ZZ, Real = double in a program, it suffices o put:
 *   #define TYPES_CODE  ZD
 * at the very beginning of the file, *before* the present file is read.
 *
 * Another (more flexible) way of specifying the flexible types
 * `(Int, Real)` is to pass the types we want to use in the class and
 * function templates. See the guide and the examples to see how to do that.
*/


/**
 * This template function returns a character string that gives the values
 * of `Int` and `Real`, usually for printing purposes.
 * For example, `strTypes<NTL::ZZ,double>(str)` will return the string
 * `Int = NTL::ZZ, Real = double` in `str`.
 */
template<typename Int, typename Real>
void strTypes(std::string &str);

template<typename Int, typename Real>
void strTypes(std::string &str) {
   str.clear();
   if (std::is_same<Int, NTL::ZZ>::value) str += "Int = NTL::ZZ, ";
   if (std::is_same<Int, long>::value) str += "Int = long, ";
   if (std::is_same<Real, double>::value) str += "Real = double";
   if (std::is_same<Real, NTL::xdouble>::value) str += "Real = xdouble";
   if (std::is_same<Real, NTL::quad_float>::value) str += "Real = quad_float";
   if (std::is_same<Real, NTL::RR>::value) str += "Real = NTL::RR";
   // std::cout << "Types: " << str << "\n\n";
}

#define   LD  1
#define   ZD  2
//#define   LX  3
#define   ZX  4
//#define   LQ  5
#define   ZQ  6
//#define   LR  7
#define   ZR  8

#define IntVec NTL::Vec<Int>
#define IntMat NTL::Mat<Int>
#define RealVec NTL::Vec<Real>
#define RealMat NTL::Mat<Real>

// std::string strFlexTypes0;

// For the following to be useful, TYPES_CODE must be defined before this file is read!
// In this case, the types `Int` and `Real` are defined here, and also string `strFlexTypes`
// used to print in outputs which flexible types we are using.
#if    TYPES_CODE == LD
	  typedef int64_t  Int;
     typedef double  Real;
     std::string strFlexTypes = "Int = int64_t, Real = double";
#elif  TYPES_CODE == ZD
	  typedef NTL::ZZ Int;
     typedef double Real;
     std::string strFlexTypes = "Int = NTL::ZZ, Real = double";
#elif  TYPES_CODE ==  ZX
     typedef NTL::ZZ Int;
     typedef NTL::xdouble Real;
     std::string strFlexTypes = "Int = NTL::ZZ, Real = xdouble";
#elif  TYPES_CODE ==  ZQ
     typedef NTL::ZZ Int;
     typedef NTL::quad_float Real;
     std::string strFlexTypes = "Int = NTL::ZZ, Real = quad_float";
#elif  TYPES_CODE ==  ZR
     typedef NTL::ZZ Int;
     typedef NTL::RR Real;
     std::string strFlexTypes = "Int = NTL::ZZ, Real = NTL::RR";
#endif

//#ifdef TYPES_CODE
//     typedef NTL::Vec<Int> IntVec;
//     typedef NTL::Mat<Int> IntMat;
//     typedef NTL::Vec<Real> RealVec;
//     typedef NTL::Mat<Real> RealMat;
//#endif

//     }
#endif
