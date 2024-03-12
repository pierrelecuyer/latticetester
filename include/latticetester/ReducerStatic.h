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

#ifndef LATTICETESTER_REDUCERSTATIC_H
#define LATTICETESTER_REDUCERSTATIC_H

#include "NTL/LLL.h"
#include "NTL/tools.h"
#include "NTL/vector.h"
#include "NTL/matrix.h"
#include "NTL/ZZ.h"
#include "NTL/RR.h"

#include "latticetester/EnumTypes.h"
#include "latticetester/Util.h"
#include "latticetester/IntLattice.h"
#include "latticetester/BasisConstruction.h"
#include "latticetester/NTLWrap.h"
#include "latticetester/LLL_FPInt.h"

#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <cstdint>
#include <iostream>
#include <ctime>
#include <string>
#include <iomanip>
#include <set>
#include <map>
#include <cmath>
#include <cstdlib>
#include <type_traits>

using namespace LatticeTester;

namespace LatticeTester {

/**
 * This file provides only static functions to reduce a lattice basis in some ways
 * (pairwise, LLL, BKZ \cite rDIE75a, \cite mLEN82a, \cite mSCH91a).
 * The LLL and BKZ ones are wrappers to slightly modified versions of the NTL
 * functions described at https://libntl.org/doc/LLL.cpp.html .
 * These functions do not require the creation of an object like in `ReducerBB`.
 * They take an `IntMat` object that contains a basis and return the reduced basis
 * in the same object. In our versions, the basis can occupy only part of the `IntMat`
 * object so we can use the same object for several bases of various sizes,
 * the basis entries can be of either ZZ or `int64_t` type,
 * the shortest basis vector is always placed in the first row,
 * and the squared vector lengths are also returned in a vector.
 * The norm to measure the vector lengths is always taken as the Euclidean one.
 * The computations can be done either in `double` or in `RR` for the real numbers.
 *
 * All these functions take as input (and return as output) an `IntMat` object
 * whose first `dim` rows and columns (a square matrix) are assumed to form
 * be a set of independent vectors that form a basis for the lattice.
 * The same `IntMat` object can then be used for several lattices of different sizes.
 */

// class ReducerStatic {
typedef NTL::vector<Int> IntVec;
typedef NTL::matrix<Int> IntMat;

//public:

/**
 * This function uses the NTL implementation of the LLL reduction algorithm
 * with factor `delta`, presented in \cite mSCH91a (see also  \cite iLEC22l).
 * The reduction is applied to the first `dim` basis vectors and coordinates
 * when `dim > 0`, and to the entire basis (all vectors) when `dim=0`.
 * In the former case, the transformations are not applied to all the
 * columns, so we will no longer have a consistent basis for the original lattice
 * if it had more than `dim` dimensions. To recover a basis for the full lattice
 * in this case, we may save it before calling this function, or rebuild it.
 *
 * This function always uses the Euclidean norm.
 * The factor `delta` must be between 1/2 and 1. The closer it is to 1,
 * the more the basis is reduced, in the sense that the LLL
 * algorithm will enforce tighter conditions on the basis.
 * The returned basis always has its shortest vector in first place.
 * The vector pointed by `sqlen` (if given) will contain the square lengths of the basis vectors.
 * To recover these values in a `Vec<double> v` one can pass `&v` to the function.
 * The parameter `precision` specifies the precision of the floating point numbers
 * that the algorithm will use. `EnumTypes.h` provides a list of the possible values,
 * and their description is done in the module `LLL` of NTL.
 */
template<typename IntMat>
static void redLLLNTL(IntMat &basis, double delta = 0.99999, long dim = 0,
      NTL::Vec<double> *sqlen = 0, PrecisionType precision = DOUBLE);

/**
 * This static function implements an exact algorithm from NTL to perform the original LLL reduction.
 * This is slower than `redLLLNTL`, but more accurate.
 * It does not take the `dim` and `sqlen` parameters (for now).
 */
template<typename IntMat>
static void redLLLNTLExact(IntMat &basis, double delta = 0.99999);

/**
 * This calls the NTL implementation of the floating point version of the
 * BKZ reduction algorithm presented in \cite mSCH91a,
 * with reduction factor `delta`, block size `blocksize`, pruning parameter
 * `prune'; see \cite iLEC22l. The other paraameters have the same meaning
 * as in `redLLLNTL`.
 * The parameter `blocksize` gives the size of the blocks in the BKZ
 * reduction. Roughly, larger blocks means a stronger condition.
 * A `blocksize` of 2 is equivalent to LLL reduction.
 */
template<typename IntMat>
static void redBKZ(IntMat &basis, double delta = 0.999999, int64_t blocksize =
      10, long prune = 0, long dim = 0, NTL::Vec<double> *sqlen = 0,
      PrecisionType prec = DOUBLE);

// }

//============================================================================
// Implementation

//===========================================================================

// namespace LatticeTester {

// General implementation.
template<typename IntMat>
void redLLLNTL(IntMat &basis, double delta, long dim, NTL::Vec<double> *sqlen,
      PrecisionType precision) {
   MyExit(1, "redLLLNTL: General version not implemented.\n");
}

// A specialization for the case where Int = int64_t and Real = double.
template<>
void redLLLNTL(NTL::matrix<int64_t> &basis, double delta, long dim,
      NTL::Vec<double> *sqlen, PrecisionType precision) {
   if (precision == DOUBLE) {
      NTL::LLL_FPInt(basis, delta, dim, dim, sqlen);
      return;
   }
   MyExit(1, "redLLLNTL: Int = int64_t with precision != DOUBLE.\n");
}

// A specialization for the case where Int = ZZ.
// See `https://github.com/u-u-h/NTL/blob/master/doc/LLL.txt` for details
// about the `PrecisionType` choices.
template<>
void redLLLNTL(NTL::matrix<NTL::ZZ> &basis, double delta, long dim,
      NTL::Vec<double> *sqlen, PrecisionType precision) {
   NTL::matrix<NTL::ZZ> cpbasis;
   switch (precision) {
   case DOUBLE:
      NTL::LLL_FPInt(basis, delta, dim, dim, sqlen);
      break;
   case QUADRUPLE:
      cpbasis.SetDims(dim, dim);
      copy(basis, cpbasis, dim, dim);  // From Util
      NTL::LLL_QP(cpbasis, delta);
      // NTL::LLL_QP_lt(basis, delta, dim, dim, sqlen);
      copy(cpbasis, basis, dim, dim);
      break;
   case XDOUBLE:
      cpbasis.SetDims(dim, dim);
      copy(basis, cpbasis, dim, dim);  // From Util
      NTL::LLL_XD(cpbasis, delta);
      // NTL::LLL_XD_lt(basis, delta, dim, dim, sqlen);
      copy(cpbasis, basis, dim, dim);
      break;
   case RR:
      NTL::LLL_RR_lt(basis, delta, dim, dim, sqlen);
      break;
   default:
      MyExit(1, "redLLLNTL: undefined PrecisionType.");
   }
}

//=========================================================================

// Exact version for the general case.
template<typename IntMat>
void redLLLNTLExact(IntMat &basis, double delta) {
   MyExit(1, "redLLLNTLExact works only for Int = ZZ");
}

// Exact version, for Int = ZZ.
template<>
void redLLLNTLExact(NTL::matrix<NTL::ZZ> &basis, double delta) {
   NTL::ZZ det(0);
   int64_t denum;
   denum = round(1.0 / (1.0 - delta)); // We want (denum-1)/denum \approx delta.
   NTL::LLL(det, basis, denum - 1, denum);
}

//=========================================================================

// BKZ, general case.
template<typename IntMat>
void redBKZ(IntMat &basis, double delta, long blocksize, long prune, long dim,
      NTL::Vec<double> *sqlen, PrecisionType prec);

// Specialization for Int = int64_t.
template<>
void redBKZ(NTL::matrix<int64_t> &basis, double delta, long blocksize,
      long prune, long dim, NTL::Vec<double> *sqlen, PrecisionType precision) {
   if (precision == DOUBLE) {
      NTL::BKZ_FPInt(basis, delta, blocksize, prune, dim, dim, sqlen);
      return;
   }
   MyExit(1, "redBKZ: Int = int64_t with precision != DOUBLE.\n");
}

// Specialization for Int = ZZ.
template<>
void redBKZ(NTL::matrix<NTL::ZZ> &basis, double delta, long blocksize,
      long prune, long dim, NTL::Vec<double> *sqlen, PrecisionType precision) {
   NTL::matrix<NTL::ZZ> cpbasis;
   switch (precision) {
   case DOUBLE:
      NTL::BKZ_FPInt(basis, delta, blocksize, prune, dim, dim, sqlen);
      return;
   case QUADRUPLE:
      cpbasis.SetDims(dim, dim);
      LatticeTester::copy(basis, cpbasis, dim, dim);  // From Util
      NTL::BKZ_QP(cpbasis, delta, blocksize);
      LatticeTester::copy(cpbasis, basis);
      break;
   case XDOUBLE:
      cpbasis.SetDims(dim, dim);
      LatticeTester::copy(basis, cpbasis, dim, dim);  // From Util
      NTL::BKZ_XD(cpbasis, delta, blocksize);
      LatticeTester::copy(cpbasis, basis);
      break;
   case RR:
      NTL::BKZ_RR_lt(basis, delta, blocksize, prune, dim, dim, sqlen);
      break;
   default:
      MyExit(1, "Undefined precision type for redBKZ");
   }
}

}   // namespace LatticeTester

#endif 
