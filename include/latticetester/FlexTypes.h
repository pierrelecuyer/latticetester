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
#include "latticetester/NTLWrap.h"  // This one is needed for the vector and matrix types.

/**
 * There are five admissible combinations of types for (Int, Real).
 * They are represented by the five codes given below.
 * For example, to use Int = NTL::ZZ, Real = double  in a program, one should put:
 *   #define TYPES_CODE  ZD
 * at the very beginning of the file.
*/

#define   LD  1
#define   ZD  2
//#define   LQ  3
#define   ZQ  4
//#define   LX  5
#define   ZX  6
//#define   LR  7
#define   ZR  8

// std::string strFlexTypes;

#if    TYPES_CODE == LD
	  typedef int64_t  Int;
     typedef double  Real;
     std::string strFlexTypes = "Int = int64_t, Real = double";
#elif  TYPES_CODE == ZD
	  typedef NTL::ZZ Int;
     typedef double Real;
     std::string strFlexTypes = "Int = NTL::ZZ, Real = double";
#elif  TYPES_CODE ==  ZQ
     typedef NTL::ZZ Int;
     typedef NTL::quad_float Real;
     std::string strFlexTypes = "Int = NTL::ZZ, Real = quad_float";
#elif  TYPES_CODE ==  ZX
     typedef NTL::ZZ Int;
     typedef NTL::xdouble Real;
     std::string strFlexTypes = "Int = NTL::ZZ, Real = xdouble";
#elif  TYPES_CODE ==  ZR
     typedef NTL::ZZ Int;
     typedef NTL::RR Real;
     std::string strFlexTypes = "Int = NTL::ZZ, Real = NTL::RR";
#endif

#ifdef TYPES_CODE
     typedef NTL::vector<Int> IntVec;
     typedef NTL::matrix<Int> IntMat;
     typedef NTL::vector<Real> RealVec;
     typedef NTL::matrix<Real> RealMat;
#endif

#endif
