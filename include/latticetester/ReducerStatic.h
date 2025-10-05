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

#include "NTL/tools.h"
#include "NTL/vector.h"
#include "NTL/matrix.h"
#include "NTL/ZZ.h"
#include "NTL/RR.h"
#include "NTL/LLL.h"

#include "latticetester/EnumTypes.h"
#include "latticetester/Util.h"
#include "latticetester/IntLattice.h"
#include "latticetester/BasisConstruction.h"
#include "latticetester/LLL_FP64.h"
#include "latticetester/LLL_lt.h"

using namespace LatticeTester;

namespace LatticeTester {

/**
 * \file latticetester/ReducerStatic.h
 *
 * Static functions to reduce a lattice basis via LLL and BKZ.
 *
 * These functions are wrappers to slightly modified versions of the NTL
 * functions described at https://libntl.org/doc/LLL.cpp.html .
 * They do not require the creation of an object like in `ReducerBB`.
 * They take an `IntMat` object that contains a basis and return the reduced basis
 * in the same object. In our versions, the basis can occupy only part of the `IntMat`
 * object so we can use the same object for several bases of various sizes,
 * the basis entries can be of either ZZ or `int64_t` type,
 * the shortest basis vector is always placed in the first row,
 * and the squared vector lengths are also returned in a vector.
 * The norm to measure the vector lengths is always the Euclidean norm.
 * The computations can be done either in `double` or in `RR` for the real numbers.
 */


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
 * The factor `delta` must be in [1/2, 1). The closer it is to 1,
 * the more the basis is reduced, in the sense that the LLL
 * algorithm will enforce tighter conditions on the basis.
 * The returned basis always has its shortest vector in first place.
 * The vector pointed by `sqlen` (if given) will contain the square lengths of the basis vectors.
 * To recover these values in a `Vec<double> v` one can pass `&v` to the function.
 * The parameter `precision` specifies the precision of the floating point numbers
 * that the algorithm will use. `EnumTypes.h` provides a list of the possible values,
 * and their description is done in the module `LLL` of NTL.
 */
template<typename Int, typename Real>
static Real redLLL(IntMat &basis, double delta = 0.99999, long dim = 0, RealVec *sqlen = 0);

/**
 * This static function implements an exact algorithm from NTL to perform the original LLL reduction.
 * This is slower than `redLLLNTL`, but more accurate.
 * It does not take the `dim` and `sqlen` parameters (for now).
 */
template<typename Int, typename Real>
static void redLLLExact(IntMat &basis, double delta = 0.99999);

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
template<typename Int, typename Real>
static Real redBKZ(IntMat &basis, double delta = 0.99999, int64_t blocksize = 10, long prune = 0,
      long dim = 0, RealVec *sqlen = 0);



//============================================================================
// Implementation

// General implementation.
template<typename Int, typename Real>
Real redLLL(IntMat &basis, double delta, long dim, RealVec *sqlen) {
   myExit("redLLL: General version not implemented.\n");
}

// A specialization for the case where Int = int64_t and Real = double.
template<>
double redLLL(NTL::Mat<int64_t> &basis, double delta, long dim, NTL::Vec<double> *sqlen) {
   return NTL::LLL_FP64(basis, delta, dim, dim, sqlen);
}

// A specialization for the case where Int = ZZ and Real = double.
template<>
double redLLL(NTL::Mat<NTL::ZZ> &basis, double delta, long dim, NTL::Vec<double> *sqlen) {
   //NTL::LLL_FPInt(basis, delta, dim, dim, sqlen);
   return NTL::LLL_FP_lt(basis, delta, dim, dim, sqlen);
}

// A specialization for the case where Int = ZZ and Real = xdouble.
template<>
xdouble redLLL(NTL::Mat<NTL::ZZ> &basis, double delta, long dim, NTL::Vec<xdouble> *sqlen) {
   return NTL::LLL_XD_lt(basis, delta, dim, dim, sqlen);
}

// A specialization for the case where Int = ZZ and Real = quad_float.
template<>
quad_float redLLL(NTL::Mat<NTL::ZZ> &basis, double delta, long dim, NTL::Vec<quad_float> *sqlen) {
   return NTL::LLL_QP_lt(basis, delta, dim, dim, sqlen);
}

// A specialization for the case where Int = ZZ and Real = RR.
template<>
NTL::RR redLLL(NTL::Mat<NTL::ZZ> &basis, double delta, long dim, NTL::Vec<NTL::RR> *sqlen) {
   return NTL::LLL_RR_lt(basis, delta, dim, dim, sqlen);
}

//=========================================================================

// Exact version for the general case.
template<typename Int>
void redLLLExact(IntMat &basis, double delta) {
   myExit("redLLLNTLExact works only for Int = ZZ");
}

// Exact version, for Int = ZZ.
template<>
void redLLLExact(NTL::Mat<NTL::ZZ> &basis, double delta) {
   NTL::ZZ det(0);
   int64_t denum;
   denum = round(1.0 / (1.0 - delta)); // We want (denum-1)/denum \approx delta.
   NTL::LLL(det, basis, denum - 1, denum);
}

//=========================================================================

// BKZ, general case.
template<typename Int, typename Real>
Real redBKZ(IntMat &basis, double delta, long blocksize, long prune, long dim, RealVec *sqlen);

// Specialization for Int = int64_t.
template<>
double redBKZ(NTL::Mat<int64_t> &basis, double delta, long blocksize, long prune, long dim,
      NTL::Vec<double> *sqlen) {
   return NTL::BKZ_FP64(basis, delta, blocksize, prune, dim, dim, sqlen);
}

// Specialization for Int = ZZ and Real = double.
template<>
double redBKZ(NTL::Mat<NTL::ZZ> &basis, double delta, long blocksize, long prune, long dim,
      NTL::Vec<double> *sqlen) {
   return NTL::BKZ_FP_lt(basis, delta, blocksize, prune, dim, dim, sqlen);
}

// Specialization for Int = ZZ and Real = xdouble.   `precision` is not used.
template<>
xdouble redBKZ(NTL::Mat<NTL::ZZ> &basis, double delta, long blocksize, long prune, long dim,
      NTL::Vec<xdouble> *sqlen) {
   return NTL::BKZ_XD_lt(basis, delta, blocksize, prune, dim, dim, sqlen);
}

// Specialization for Int = ZZ and Real = quad_float.   `precision` is not used.
template<>
quad_float redBKZ(NTL::Mat<NTL::ZZ> &basis, double delta, long blocksize, long prune, long dim,
      NTL::Vec<quad_float> *sqlen) {
   return NTL::BKZ_QP_lt(basis, delta, blocksize, prune, dim, dim, sqlen);
}

// Specialization for Int = ZZ and Real = RR.   `precision` is not used.
template<>
NTL::RR redBKZ(NTL::Mat<NTL::ZZ> &basis, double delta, long blocksize, long prune, long dim,
      NTL::Vec<NTL::RR> *sqlen) {
   return NTL::BKZ_RR_lt(basis, delta, blocksize, prune, dim, dim, sqlen);
}

}   // namespace LatticeTester

#endif 
