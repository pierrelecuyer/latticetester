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

#ifndef LATTICETESTER_BASISCONSTRUCTION_H
#define LATTICETESTER_BASISCONSTRUCTION_H

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <set>
#include <map>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <cassert>
#include <type_traits>

#include "NTL/tools.h"
#include "NTL/ZZ.h"
#include "NTL/vec_ZZ.h"
#include "NTL/mat_ZZ.h"
//#include "NTL/LLL.h"

#include <latticetester/FlexTypes.h>
#include "latticetester/Util.h"
#include "latticetester/EnumTypes.h"
//#include "latticetester/IntLattice.h"
#include "latticetester/Coordinates.h"
#include "latticetester/LLL_FP64.h"
#include "latticetester/LLL_lt.h"

using namespace LatticeTester;
using namespace NTL;

/**
 * \file latticetester/BasisConstruction.h
 *
 * Static functions to construct a basis from a set of generating
 * vectors that are not necessarily independent, to construct a triangular basis,
 * to construct the basis for a projection over a given subset of coordinates,
 * and to obtain the \f$m\f$-dual of a given basis.
 *
 * The algorithms are described in the Lattice Tester guide \cite iLEC25ltg.
 * The implementation relies on %NTL and uses %NTL matrices.
 * When the basis turns out to have fewer rows than columns, most of the functions here
 * add implicitly the rescaled unit vectors to the set of generating vectors.
 * In that case, the basis matrix is always square and all the vectors of the form
 * \f$m \mathbf{e}_i\f$ belong to the lattice.
 *
 * %NTL already offers an efficient procedure to construct an LLL-reduced basis from a set
 * of generating vectors.  It is encapsulated in the `LLLConstruction0` function given below.
 * This function does not assume that the rescaled unit vectors \f$m \mathbf{e}_i \f$ belong
 * to the lattices and it does not even know about \f$m\f$.
 * The function `LLLBasisConstruction` adds those vectors to the set of generating vectors,
 * so it always returns a square basis.
 *
 * We also offer an alternative functions that construct a triangular basis from a set of
 * generating vectors. These functions are usually faster.
 * They always add the rescaled unit vectors implicitly to the set.
 * The function `lowerTriangularBasis` constructs a lower-triangular basis, while
 * `upperTriangularBasis` constructs an upper-triangular basis.
 *
 * To compute the \f$m\f$-dual of a given basis, we have a general (but slow) function
 * implemented in `mDualBasis`, and much faster functions in `mDualLowerTriangular`
 * and `mDualUpperTriangular` that works only when the basis is triangular.
 * The functions `mDualLowerTriangularMod0` and `mDualUpperTriangularMod0`
 * compute these \f$m\f$-dual bases, and also reduce their non-diagonal entries
 * mod \f$m\f$ towards 0, giving bases with shorter vectors for the \f$m\f$-dual lattice.
 *
 * We also have functions to compute the basis of a projection of a given lattice over
 * a specified set of coordinates.  The function `projectionConstructionLLL` does this
 * by using LLL to construct the basis of the projection, while `projectionConstructionUpperTri`
 * constructs an upper-triangular basis for the projection.
 * The function `projectionConstruction` takes the construction method as a parameter.
 *
 * All functions take as input a `IntMat` object that contains either a basis or a set of
 * generating vectors. By default, all rows and columns of this matrix object are used.
 * But in most functions, the user can ask to use only a subset of these rows and columns,
 * via the optional parameters `r` and `c`.  This can permit one to use the same `IntMat`
 * object for several numbers of dimensions, to avoid doing many object creations or resizing.
 *
 * All functions in this file are static, so there is no notion of
 * `BasisConstruction` object. We also avoid to create new objects (such as vectors and
 * matrices) inside these functions.  These functions can be called thousands or millions
 * of times in a program, and we want the user to be able to re-use the same vectors and
 * matrices over and over again instead of creating new ones.
 *
 * Note that when one of these functions is used for an `IntMat` object that store the primal
 * or \f$m\f$-dual basis inside an `IntLattice` object, only the primal or the \f$m\f$-dual
 * basis is changed, and therefore the  \f$m\f$-dual relationship will no longer hold after
 * the call. The user must be aware of that.  Most of the time, there is no need to reinstate
 * this relationship, because for example if we want to compute a shortest nonzero vector
 * is the \f$m\f$-dual lattice, then once we have an \f$m\f$-dual basis we no longer need
 * a primal basis.
 *
 * The programs `TestBasisConstructSmall` and `TestBasisConstructSpeed` in the examples
 * illustrate how to use these functions and make speed comparisons.
 */

namespace LatticeTester {

/**
 * This function takes a set of generating vectors of a lattice in matrix `gen` and
 * finds a lattice basis by applying LLL reduction with the given value of `delta`,
 * using the NTL implementation specified by `prec`.
 * The basis is returned in the first rows of `gen`, in a number of rows equal to its rank.
 * See the file `EnumTypes` and the documentation of LLL in NTL for the meaning and
 * choices for `prec`. The default value of delta is not very close to 1 because the
 * goal here is just to compute a basis. It can be taken much closer to 1 if we
 * really prefer a highly-reduced basis.
 * The optional parameters `r` and `c` indicate the numbers of rows and columns
 * of the matrix `gen` that are actually used. When they are 0 (the default values),
 * then these numbers are taken to be the dimensions of the matrix `gen`.
 * When `sqlen` is not 0, the square lengths of the basis
 * vectors are returned in this array, exactly as in the `LLL_lt.h` module.
 * The function returns the square length of the shortest basis vector.
 * The basis `gen` is never resized.
 *
 * This function *does not* assume that all vectors \f$m e_i\f$ belong to the lattice, so
 * it may return a basis matrix that has fewer rows than columns!
 * To make sure that these vectors belong to the lattice, we can add them
 * explicitly beforehand to the set of generating vectors, or call the next function,
 * which does that.
 */
template<typename Int, typename Real>
static Real LLLConstruction0(IntMat &gen, const double delta = 0.9, long r = 0, long c = 0,
RealVec *sqlen = 0);

/**
 * Similar to `LLLConstruction0`, except that this function adds explicitly all the vectors
 * \f$m \mathbf{e}_i\f$ to the generating set, so there will be `r + c` rows initially and
 * the function always returns a square matrix.
 * The matrix `gen` is not resized by this function, so it can remain larger
 * than the lattice dimension.
 */
template<typename Int, typename Real>
static Real LLLBasisConstruction(IntMat &gen, const Int &m, const double delta = 0.9, long r = 0,
      long c = 0, RealVec *sqlen = 0);

/**
 * Takes a set of generating vectors in the matrix `gen` and iteratively
 * transforms it into a lower triangular lattice basis into the matrix `basis`.
 * This lattice is assumed to contain all the vectors of the form @f$m \mathbf{e}_j@f$,
 * so these vectors are added implicitly to the generating set.
 * Apart from that, all the entries of `gen` given as input are assumed to be
 * reduced modulo the scaling factor `m` and all the computations are done modulo `m`.
 * The matrix `basis` is assumed to be large enough to contain the new basis,
 * whose dimension should be the number of columns taken from `gen`.
 * After the execution, `gen` will contain irrelevant information (garbage)
 * and `basis` will contain an upper triangular basis.
 * The algorithm is explained in the lattice tester guide.
 * Important: `gen` and `basis` must be different `IntMat` objects.
 *
 * When `r` and/or `c` are strictly positive, they specify the numbers of
 * rows and columns of `gen` that are actually used, as in `LLLConstruction0`.
 * The matrix `gen` is never resized.  The new basis should be c x c.
 */
template<typename Int>
static void lowerTriangularBasis(IntMat &basis, IntMat &gen, const Int &m, long r = 0, long c = 0);

/**
 * Same as `lowerTriangularBasis`, except that the returned basis is upper triangular.
 */
template<typename Int>
static void upperTriangularBasis(IntMat &basis, IntMat &gen, const Int &m, long r = 0, long c = 0);

/**
 * The old version from \cite rCOU96a and \cite iLEC00l.
 */
template<typename Int>
static void upperTriangularBasisOld96(IntMat &basis, IntMat &gen, const Int &m, long r = 0, long c =
      0);

/**
 * Takes a lower-triangular basis matrix `basis` and computes the m-dual basis `basisDual`.
 * The function assumes that each coefficient on the diagonal of `basis` is nonzero and divides `m`.
 * That is, the basis matrix must be square and m-invertible.
 * Since the basis is lower triangular, its m-dual will be upper triangular.
 * When `dim > 0`, it must give the number of rows and columns of the matrix `basis`
 * that is actually used. Otherwise (by default) the function uses `basis.numCols()`.
 */
template<typename Int>
static void mDualLowerTriangular(IntMat &basisDual, const IntMat &basis, const Int &m,
      long dim = 0);

/**
 * Same as `mDualLowerTriangular` and then reduces all the non-diagonal `mod m` towards 0.
 */
template<typename Int>
static void mDualLowerTriangularMod0(IntMat &basisDual, const IntMat &basis, const Int &m,
      long dim = 0);

/**
 * Takes a upper triangular basis matrix `basis` and computes the m-dual basis `basisDual`.
 * This function is the equivalent of mDualLowerTriangular for upper-triangular matrices.
 */
template<typename Int>
static void mDualUpperTriangular(IntMat &basisDual, const IntMat &basis, const Int &m,
      long dim = 0);

/**
 * Same as `mDualUpperTriangular` and then reduces all the non-diagonal `mod m` towards 0.
 */
template<typename Int>
static void mDualUpperTriangularMod0(IntMat &basisDual, const IntMat &basis, const Int &m,
      long dim = 0);

/**
 * This function does essentially the same thing as `mDualUpperTriangular`, but the
 * algorithm is slightly different. It uses the method described in \cite rCOU96a.
 */
template<typename Int>
static void mDualUpperTriangularOld96(IntMat &basisDual, const IntMat &basis, const Int &m,
      long dim = 0);

/**
 * This function assumes that `basis` contains a basis of the primal lattice
 * scaled by the factor `m`, not necessarily triangular, and it returns in `basisDual`
 * the m-dual basis. It uses matrix inversion and is rather slow.
 * It is currently implemented only for `Int = ZZ` and it also assumes that the dimensions
 * of the two `IntMat` objects is exactly the same as the dimensions of the lattices.
 * The reason for this is that we use an NTL function that works only under these conditions.
 */
template<typename Int>
static void mDualBasis(IntMat &basisDual, const IntMat &basis, const Int &m);

/**
 * This function overwrites the first `r` rows of matrix 'out' by a matrix formed by the
 * first `r` rows of the c columns of matrix `in` that are specified by `proj`,
 * where `c = size(proj)' is the cardinality of the projection `proj`.
 * Only the first c columns of the first `r` rows are overwritten;
 * the other entries are left unchanged.
 * If `r = 0`, then all the rows of matrix `in` are taken.
 * When `in` is a basis, `r` will usually be the dimension of that basis and `c` will be smaller.
 * The dimensions of the matrices `in` and `out` are always left unchanged.
 * The dimensions of `out` are assumed to be large enough.   After the call,
 * the matrix `out` will then contain a set of generating vectors for the projection `proj`.
 * The matrices `in` and `out` must be different IntMat objects, otherwise the program
 * halts with an error message.
 */
template<typename Int>
static void projectMatrix(IntMat &out, const IntMat &in, const Coordinates &coordSet, long r = 0);

/**
 * Constructs a basis for the projection onto `coordSet` of the lattice with basis `inBasis`,
 * using `LLLBasisConstruction`, and returns it in `projBasis`.
 * This returned basis is not triangular in general.
 * Its dimension will be the number of coordinates in `coordSet`.
 * The matrix `projBasis` must have enough columns to hold it and at least as many
 * rows as the number of rows that we use from `inBasis`.
 * When `r > 0`, only the first `r` rows of the matrix `inBasis` are used,
 * otherwise we use all rows of that matrix.
 * This `r` should normally be the true dimension `dim` of that basis, which is often
 * smaller than the size `maxDim` of the `IntMat` object that contains the basis.
 * The square Euclidean lengths of the basis vectors are returned in the array `sqlen` when
 * the latter is given.
 */
template<typename Int, typename Real>
static void projectionConstructionLLL(IntMat &projBasis, const IntMat &inBasis,
      const Coordinates &coordSet, const Int &m, const double delta = 0.9, long r = 0, RealVec *sqlen =
            0);

/**
 * Same as `projectionConstructionLLL`, but the construction is made using
 * `upperTriangularBasis`, so the returned basis is upper triangular.
 * When `r > 0`, only the first `r` rows of the matrix `inBasis` are actually used,
 * otherwise we use all rows of that matrix.
 * In the first version, we pass a matrix `genTemp` that will be used to store the
 * generating vectors of the projection before making the triangularization.
 * The two matrices `projBasis` and `genTemp` must have enough columns to hold the
 * projection and at least as many rows as the number of rows that we use from `inBasis`.
 * We pass `genTemp` as a parameter to avoid the internal creation of a new matrix each time,
 * in case we call this function several times.  Its contents will be modified.
 * In the second version, this matrix is not passed and a temporary one is created internally,
 * which may add a bit of overhead.
 */
template<typename Int>
static void projectionConstructionUpperTri(IntMat &projBasis, const IntMat &inBasis,
IntMat &genTemp, const Coordinates &coordSet, const Int &m, long r = 0);

template<typename Int>
static void projectionConstructionUpperTri(IntMat &projBasis, const IntMat &inBasis,
      const Coordinates &coordSet, const Int &m, long r = 0);

/**
 * In this version, the construction method is passed as a parameter. The default is LLL.
 * In the triangular case, a temporary matrix is created internally.
 */
template<typename Int>
static void projectionConstruction(IntMat &projBasis, const IntMat &inBasis,
      const Coordinates &coordSet, const Int &m, const ProjConstructType projType = LLLPROJ,
      const double delta = 0.9);

//============================================================================
// Implementation

// General case, no implementation.
template<typename Int, typename Real>
static Real LLLConstruction0(IntMat &gen, const double delta, long r, long c, RealVec *sqlen) {
   std::cerr << "LLLConstruction0: general case is not implemented.\n";
   exit(1);
}

// The int64_t implementation.
// This one works only for `precision == DOUBLE` and Real == double.
template<>// <long, double>
double LLLConstruction0(NTL::Mat<long> &gen, const double delta, long r, long c,
      NTL::Vec<double> *sqlen) {
   return NTL::LLL_FP64(gen, delta, r, c, sqlen);
   // return NTL::LLL_FPInt<long>(gen, delta, r, c, sqlen);
}

// The ZZ + double implementation.
template<>// <NTL::ZZ, double>
double LLLConstruction0(NTL::Mat<NTL::ZZ> &gen, const double delta, long r, long c,
      NTL::Vec<double> *sqlen) {
   return NTL::LLL_FP_lt(gen, delta, r, c, sqlen);
}

// The ZZ + xdouble implementation.
template<>
xdouble LLLConstruction0(NTL::Mat<NTL::ZZ> &gen, const double delta, long r, long c,
      NTL::Vec<xdouble> *sqlen) {
   return NTL::LLL_XD_lt(gen, delta, r, c, sqlen);
}

// The ZZ + quad_float implementation.
template<>
quad_float LLLConstruction0(NTL::Mat<NTL::ZZ> &gen, const double delta, long r, long c,
      NTL::Vec<quad_float> *sqlen) {
   return NTL::LLL_QP_lt(gen, delta, r, c, sqlen);
}

// The ZZ + RR implementation.
template<>
RR LLLConstruction0(NTL::Mat<NTL::ZZ> &gen, const double delta, long r, long c,
      NTL::Vec<NTL::RR> *sqlen) {
   return NTL::LLL_RR_lt(gen, delta, r, c, sqlen);
}

//===========================================================================

template<typename Int, typename Real>
Real LLLBasisConstruction(IntMat &gen, const Int &m, double delta, long r, long c, RealVec *sqlen) {
   if (c == 0) c = gen.NumCols();
   if (r == 0) r = gen.NumRows();
   // We add the m e_i row vectors and we perform the LLL with that.
   // No change in the dimensions of gen.
   int64_t i, j;
   for (i = r; i < r + c; i++) {
      for (j = 0; j < c; j++) {
         if (i == j) gen[i][j] = m;
         else gen[i][j] = 0;
      }
   }
   // std::cout << "Warning for LLLBasisConstruction: we had to add some rows!\n";
   // std::cout << "  c = " << c << ", rank = " << rank << "\n";
   return LLLConstruction0<Int, Real>(gen, delta, r + c, c, sqlen);
}

//==============================================================================

template<typename Int>
void lowerTriangularBasis(IntMat &basis, IntMat &gen, const Int &m, long dim1, long dim2) {
   // `dim1` and `dim2` are `s` and `t` in the guide.
   if (dim1 == 0) dim1 = gen.NumRows();
   if (dim2 == 0) dim2 = gen.NumCols();
   assert(dim2 <= basis.NumRows() && dim2 <= basis.NumCols());

   long i, j, k, l;
   Int gcd, gcdCopy;  // The gcd c_j in the guide, computed incrementally.
   Int c, d;          // Coefficients returned by the XGCD function.
   NTL::Vec<Int> coeff_gcd, coeff_xj, xj;  // Here we create new vectors!
   coeff_gcd.SetLength(dim1);  // The coefficients a_{1,j},...,a_{s,j} in the guide.
   coeff_xj.SetLength(dim1);   // The coefficients a_{1,j},...,a_{s,j} that define xj.
   xj.SetLength(dim2);         // The new basis vector x_j computed at step j.

   for (j = dim2 - 1; j >= 0; j--) {  // column j.
      // Here we compute the submatrix whose upper left corner is (j,j) in the upper triangular basis.
      // Find c_j and the coefficients `a_{i,j}` by applying the Euclidean algorithm multiple times.
      for (i = 0; i < dim1; i++) {
         NTL::rem(gen[i][j], gen[i][j], m);  // Reduce modulo m
         coeff_gcd[i] = 0;
      }
      gcd = m;
      for (i = dim1 - 1; i >= 0; i--) {
         // NTL::rem(gen[i][j], gen[i][j], m);
         if (gen[i][j] != 0) {
            // XGCD (g, c, d, const a, const b) does g = gcd(a, b) = a*c + b*d.
            gcdCopy = gcd;  // We need a copy for the `const a` parameter.
            NTL::XGCD(gcd, c, d, gcdCopy, gen[i][j]);
            coeff_gcd[i] = d;
            for (l = i + 1; l < dim1; l++) {
               NTL::mul(coeff_gcd[l], coeff_gcd[l], c);
               // NTL::rem(coeff_gcd[l], coeff_gcd[l], m);
            }
         }
      }
      // If `gcd = m`, then this basis (row) vector will be `x_j = m e_j`.
      if (gcd == m) {
         for (k = 0; k < dim2; k++) {
            if (k != j) basis[j][k] = 0;
            else basis[j][k] = m;
         }
      } else {
         // We now compute `x_j` and modify the `v_{i,k}` for `k <= j`.
         for (k = 0; k < dim2; k++)
            xj[k] = 0;
         for (i = 0; i < dim1; i++) {
            if (coeff_gcd[i] != 0) {
               for (k = 0; k <= j; k++) {
                  NTL::MulAddTo(xj[k], gen[i][k], coeff_gcd[i]);
                  // std::cout << "  Before: xj[k] = " << xj[k] << "\n";
                  ModuloTowardZero(xj[k], m, xj[k]);
                  // NTL::rem(xj[k], xj[k], m);     // In case we want always a positive remainder (modulo).
                  // std::cout << "  After: xj[k] = " << xj[k] << "\n";
                  // if (xj[k] < 0) NTL::add(xj[k], xj[k], m);    // No!
               }
            }
         }
         //std::cout << "  j = " << j << "\n";
         //std::cout << "  vector x_j = " << xj << "\n";
         // Next we update the vectors v_i.
         // We first calculate the coefficients with which x_j needs to be multiplied.
         for (i = 0; i < dim1; i++) {
            NTL::div(coeff_xj[i], gen[i][j], gcd);
            // std::cout << "  coeff_xj[i] = " << coeff_xj[i] << "\n";
            // NTL::rem(coeff_xj[i], coeff_xj[i], m);
            // std::cout << "  coeff_xj[i] = " << coeff_xj[i] << "\n";
         }
         //std::cout << "  vector x_j = " << xj << "\n";
         //std::cout << "  coeff_x_j = " << coeff_xj << "\n";
         // Update the components of index <= j of the old vectors v_i.
         for (i = 0; i < dim1; i++) {
            if (coeff_xj[i] != 0) {
               for (k = j; k >= 0; k--) {
                  NTL::MulSubFrom(gen[i][k], coeff_xj[i], xj[k]);
                  // if (abs(gen[i][k]) > m * m)  std::cout << "  vector gen[i][k] = " << gen[i][k] << "\n";
                  NTL::rem(gen[i][k], gen[i][k], m);
               }
            }
         }
         // This is row j of the upper-triangular basis.
         for (k = 0; k < dim2; k++)
            basis[j][k] = xj[k];
      }
   }
}

//===================================================================

template<typename Int>
void upperTriangularBasis(IntMat &basis, IntMat &gen, const Int &m, long dim1, long dim2) {
   // `dim1` and `dim2` are `s` and `t` in the guide.
   if (dim1 == 0) dim1 = gen.NumRows();
   if (dim2 == 0) dim2 = gen.NumCols();
   assert(dim2 <= basis.NumRows() && dim2 <= basis.NumCols());

   long i, j, k, l;
   Int gcd, gcdCopy;  // The gcd c_j in the guide, computed incrementally.
   Int c, d;          // Coefficients returned by the XGCD function.
   NTL::Vec<Int> coeff_gcd, coeff_xj, xj;  // Here we create new vectors!
   coeff_gcd.SetLength(dim1);  // The coefficients a_{1,j},...,a_{s,j} in the guide.
   coeff_xj.SetLength(dim1);   // The coefficients a_{1,j},...,a_{s,j} that define xj.
   xj.SetLength(dim2);         // The new basis vector x_j computed at step j.

   for (j = 0; j < dim2; j++) {  // column j.
      // Here we compute the submatrix whose upper left corner is (j,j) in the upper triangular basis.
      // Find c_j and the coefficients `a_{i,j}` by applying the Euclidean algorithm multiple times.
      for (i = 0; i < dim1; i++) {
         NTL::rem(gen[i][j], gen[i][j], m);  // Reduce modulo m to a non-negative value.
         //  gen[i][j] = gen[i][j] % m;      // Doing it this way is slower.
         coeff_gcd[i] = 0;
      }
      gcd = m;
      for (i = 0; i < dim1; i++) {
         if (gen[i][j] != 0) {
            // XGCD (g, c, d, const a, const b) does g = gcd(a, b) = a*c + b*d.
            gcdCopy = gcd;  // We need a copy for the `const a` parameter.
            NTL::XGCD(gcd, c, d, gcdCopy, gen[i][j]);
            coeff_gcd[i] = d;
            for (l = 0; l < i; l++) {
               NTL::mul(coeff_gcd[l], coeff_gcd[l], c);
               //NTL::rem(coeff_gcd[l], coeff_gcd[l], m);
            }
         }
      }
      // std::cout << "j = " << j << ", gcd = " << gcd << "\n";
      // If `gcd = m`, then this basis (row) vector will be `x_j = m e_j`.
      if (gcd == m) {
         for (k = 0; k < dim2; k++) {
            if (k != j) basis[j][k] = 0;
            else basis[j][k] = m;
         }
      } else {
         // We now compute `x_j` and modify the `v_{i,k}` for `k >= j`.
         for (k = 0; k < dim2; k++)
            xj[k] = 0;
         for (i = 0; i < dim1; i++) {
            if (coeff_gcd[i] != 0) {
               for (k = j; k < dim2; k++) {
                  NTL::MulAddTo(xj[k], gen[i][k], coeff_gcd[i]);
                  //NTL::rem(xj[k], xj[k], m);
                  ModuloTowardZero(xj[k], m, xj[k]);
               }
            }
         }
         // std::cout << "  vector x_j = " << xj << "\n";
         // Next we update the vectors v_i.
         // We first calculate the coefficients with which x_j needs to be multiplied.
         for (i = 0; i < dim1; i++) {
            NTL::div(coeff_xj[i], gen[i][j], gcd);
            NTL::rem(coeff_xj[i], coeff_xj[i], m);
         }
         // std::cout << "  coeff_x_j = " << coeff_xj << "\n";
         // Update the components of index >= j of the old vectors v_i.
         for (i = 0; i < dim1; i++) {
            if (coeff_xj[i] != 0) {
               for (k = j; k < dim2; k++) {
                  NTL::MulSubFrom(gen[i][k], coeff_xj[i], xj[k]);
               }
            }
         }
         // This is row j of the upper-triangular basis.
         for (k = 0; k < dim2; k++)
            basis[j][k] = xj[k];
      }
   }
}

//===================================================

//  This is the old triangularization method that we had in Modula-2 in 1996.
template<typename Matr, typename Int>
void upperTriangularBasisOld96(Matr &V, Matr &W, const Int &m, int64_t lin, int64_t col) {
   Int T1, T2, T3, T4, T5, T6, T7, T8;
   for (int64_t j = 0; j < col; j++) {
      for (int64_t i = 0; i < lin; i++) {
         // std::cout << "Before Modulo, j = " << j << ", i = " << i << "\n ";
         Modulo(W[i][j], m, W[i][j]);
      }
      int64_t r = 0;
      while (r < lin - 1) {
         while (NTL::IsZero(W[r][j]) && r < lin - 1)
            ++r;
         if (r < lin - 1) {
            int64_t s = r + 1;
            while (NTL::IsZero(W[s][j]) && s < lin - 1)
               ++s;
            if (!NTL::IsZero(W[s][j])) {
               Int temp;
               Euclide(W[r][j], W[s][j], T1, T2, T3, T4, temp);
               W[s][j] = temp;

               NTL::clear(W[r][j]);

               for (int64_t j1 = j + 1; j1 < col; j1++) {
                  T5 = T1 * W[r][j1];
                  T6 = T2 * W[s][j1];
                  T7 = T3 * W[r][j1];
                  T8 = T4 * W[s][j1];
                  W[s][j1] = T5 + T6;
                  Modulo(W[s][j1], m, W[s][j1]);
                  W[r][j1] = T7 + T8;
                  Modulo(W[r][j1], m, W[r][j1]);
               }
            } else {
               for (int64_t j1 = j; j1 < col; j1++) {
                  std::swap(W[r][j1], W[s][j1]);
               }
            }
            r = s;
         }
      }
      if (NTL::IsZero(W[lin - 1][j])) {
         for (int64_t j1 = 0; j1 < col; j1++) {
            if (j1 != j) NTL::clear(V[j][j1]);
            else V[j][j1] = m;
         }
      } else {
         Euclide(W[lin - 1][j], m, T1, T2, T3, T4, V[j][j]);
         for (int64_t j1 = 0; j1 < j; j1++)
            NTL::clear(V[j][j1]);
         for (int64_t j1 = j + 1; j1 < col; j1++) {
            T2 = W[lin - 1][j1] * T1;
            Modulo(T2, m, V[j][j1]);
         }
         Quotient(m, V[j][j], T1);
         for (int64_t j1 = j + 1; j1 < col; j1++) {
            W[lin - 1][j1] *= T1;
            Modulo(W[lin - 1][j1], m, W[lin - 1][j1]);
         }
      }

   }
}

//===================================================

template<typename Int>
void mDualLowerTriangular(IntMat &B, const IntMat &A, const Int &m, long dim) {
   // Note:  A = basis,  B = basisDual is upper triangular
   if (dim == 0) dim = A.NumRows();
   assert(dim <= A.NumCols());
   assert(dim <= B.NumRows() && dim <= B.NumCols());
   for (int64_t i = 0; i < dim; i++) {
      // Put zeros under the diagonal.
      for (int64_t j = 0; j < i; j++)
         NTL::clear(B[i][j]);
      // Set diagonal elements.
      NTL::div(B[i][i], m, A[i][i]);
      // Compute the other ones.
      for (int64_t j = i + 1; j < dim; j++) {
         NTL::clear(B[i][j]);
         for (int64_t k = i; k < j; k++)
            NTL::MulSubFrom(B[i][j], A[j][k], B[i][k]);
         NTL::div(B[i][j], B[i][j], A[j][j]);
      }
   }
}

//===================================================

template<typename Int>
void mDualUpperTriangular(IntMat &B, const IntMat &A, const Int &m, long dim) {
   // Note:  A = basis,  B = basisDual is lower triangular
   if (dim == 0) dim = A.NumRows();
   assert(dim <= A.NumCols());
   assert(dim <= B.NumRows() && dim <= B.NumCols());
   for (int64_t i = 0; i < dim; i++) {
      // Put zeros above the diagonal.
      for (int64_t j = i + 1; j < dim; j++)
         NTL::clear(B[i][j]);
      // Set diagonal elements.
      NTL::div(B[i][i], m, A[i][i]);
      // Compute the other ones.
      for (int64_t j = i - 1; j >= 0; j--) {
         NTL::clear(B[i][j]);
         for (int64_t k = j + 1; k <= i; k++)
            NTL::MulSubFrom(B[i][j], A[j][k], B[i][k]);
         NTL::div(B[i][j], B[i][j], A[j][j]);
      }
   }
}

//===================================================

template<typename Int>
void mDualLowerTriangularMod0(IntMat &B, const IntMat &A, const Int &m, long dim) {
   mDualLowerTriangular<Int>(B, A, m, dim);
   for (int64_t i = 0; i < dim; i++)
      for (int64_t j = i + 1; j < dim; j++)
         ModuloTowardZero(B[i][j], m, B[i][j]);
}

//===================================================

template<typename Int>
void mDualUpperTriangularMod0(IntMat &B, const IntMat &A, const Int &m, long dim) {
   mDualUpperTriangular<Int>(B, A, m, dim);
   for (int64_t i = 0; i < dim; i++)
      for (int64_t j = i - 1; j >= 0; j--)
         ModuloTowardZero(B[i][j], m, B[i][j]);
}

//======================================================
// This is the old Modula-2 version from Couture and L'Ecuyer (1996).
template<typename Int>
void mDualUpperTriangularOld96(IntMat &basisDual, const IntMat &basis, const Int &m, long dim) {
   if (!dim) dim = basis.NumRows();
   assert(dim <= basisDual.NumRows() && dim <= basisDual.NumCols());
   Int gcd;
   Int mm = m;            // Local copy of m that can be changed.
   for (int64_t i = 0; i < dim; i++) {
      for (int64_t j = i + 1; j < dim; j++)
         NTL::clear(basisDual[i][j]);
      if (!NTL::IsZero(basisDual[i][i])) {
         gcd = NTL::GCD(mm, basis[i][i]);
         mm *= basis[i][i] / gcd;
         basisDual *= basis[i][i] / gcd;
      }
      DivideRound(mm, basis[i][i], basisDual[i][i]);
      for (int64_t j = i - 1; j >= 0; j--) {
         NTL::clear(basisDual[i][j]);
         for (int64_t k = j + 1; k <= i; k++)
            basisDual[i][j] += basis[j][k] * basisDual[i][k];
         if (basisDual[i][j] != 0) basisDual[i][j] = -basisDual[i][j];
         if (!NTL::IsZero(basisDual[i][j] % basis[j][j])) {
            gcd = NTL::GCD(basisDual[i][j], basis[j][j]);
            mm *= basis[j][j] / gcd;
            basisDual[i][j] *= basis[j][j] / gcd;
         }
         DivideRound(basisDual[i][j], basis[j][j], basisDual[i][j]);
      }
   }
}

// A specialized implementation for ZZ. Not much different from the Int version.
template<>
void mDualUpperTriangularOld96(NTL::Mat<NTL::ZZ> &basisDual, const NTL::Mat<NTL::ZZ> &basis,
      const NTL::ZZ &m, long dim) {
   if (!dim) dim = basis.NumRows();
   assert(dim <= basisDual.NumRows() && dim <= basisDual.NumCols());
   NTL::ZZ gcd, fac;
   NTL::ZZ mm = m;            // Local copy of m that can be changed.
   for (int64_t i = 0; i < dim; i++) {
      for (int64_t j = i + 1; j < dim; j++)
         NTL::clear(basisDual[i][j]);
      if (!NTL::IsZero(basisDual[i][i])) {
         gcd = NTL::GCD(mm, basis[i][i]);
         div(fac, basis[i][i], gcd);
         mul(mm, mm, fac);
         mul(basisDual[i][i], basisDual[i][i], fac);
      }
      DivideRound(mm, basis[i][i], basisDual[i][i]);
      for (int64_t j = i - 1; j >= 0; j--) {
         NTL::clear(basisDual(i, j));
         for (int64_t k = j + 1; k <= i; k++)
            basisDual[i][j] += basis[j][k] * basisDual[i][k];
         //if (basisDual(i, j) != 0)
         NTL::negate(basisDual[i][j], basisDual[i][j]);
         if (!NTL::IsZero(basisDual[i][j] % basis[j][j])) {
            gcd = NTL::GCD(basisDual[i][j], basis[j][j]);
            div(fac, basis[j][j], gcd);
            mul(mm, mm, fac);
            mul(basisDual[i][j], basisDual[i][j], fac);
         }
         div(basisDual[i][j], basisDual[i][j], basis[j][j]);
      }
   }
}

// ===================================================

template<typename Int>
void mDualBasis(NTL::Mat<Int> &basisDual, const NTL::Mat<Int> &basis, const Int &m) {
   std::cerr << "mDualBasis is implemented only for NTL::ZZ integers.\n";
   exit(1);
}

// The specialization for the case where `Int = ZZ`.
template<>
void mDualBasis(NTL::Mat<NTL::ZZ> &basisDual, const NTL::Mat<NTL::ZZ> &basis, const NTL::ZZ &m) {
   NTL::ZZ det, fac;
   long dim = basis.NumRows();
   if (dim != basis.NumCols()) {
      std::cerr << "mDualBasis: the given basis matrix must be square.\n";
      exit(1);
   }
   inv(det, basisDual, basis);
   NTL::Mat < NTL::ZZ > C = basisDual;
   div(fac, det, m);
   for (int64_t i = 0; i < dim; i++) {
      for (int64_t j = 0; j < dim; j++) {
         div(basisDual[j][i], C[i][j], fac);
      }
   }
}

//=================================================================================
template<typename Int>
void projectMatrix(IntMat &out, const IntMat &in, const Coordinates &coordSet, long r) {
   if (in == out) {
      std::cout << "\n***** Error: in and out must be different IntMat objects " << std::endl;
      exit(1);
   }
   if (!r) r = in.NumRows();   // In case r=0.
   // We assume without testing that `out` is large enough for coordSet.size().
   long j = 0;
   for (auto it = coordSet.begin(); it != coordSet.end(); it++, j++) {
      for (long i = 0; i < r; i++)
         out[i][j] = in[i][*it - 1];
   }
}

//===================================================
template<typename Int, typename Real>
void projectionConstructionLLL(IntMat &projBasis, const IntMat &inBasis, const Coordinates &coordSet,
      const Int &m, const double delta, long r, RealVec *sqlen) {
   projectMatrix(projBasis, inBasis, coordSet, r);
   LLLBasisConstruction(projBasis, m, delta, r, coordSet.size(), sqlen);
}

//===================================================
template<typename Int>
void projectionConstructionUpperTri(IntMat &projBasis, const IntMat &inBasis, IntMat &genTemp,
      const Coordinates &coordSet, const Int &m, long r) {
   projectMatrix(genTemp, inBasis, coordSet, r);
   upperTriangularBasis(projBasis, genTemp, m, r, coordSet.size());
}

template<typename Int>
void projectionConstructionUpperTri(IntMat &projBasis, const IntMat &inBasis,
      const Coordinates &coordSet, const Int &m, long r) {
   long dim = coordSet.size();      // Dimension of projection
   if (!r) r = inBasis.NumRows();
   IntMat genTemp;
   genTemp.SetDims(r, dim); // Here an internal object is created and resized!
   projectMatrix(genTemp, inBasis, coordSet, r);
   upperTriangularBasis(projBasis, genTemp, m, r, dim);
}

//===================================================

template<typename Int>
void projectionConstruction(IntMat &projBasis, const IntMat &inBasis, const Coordinates &coordSet,
      const Int &m, const ProjConstructType projType, const double delta) {
   if (projType == LLLPROJ) projectionConstructionLLL(projBasis, inBasis, coordSet, m, 0, delta);
   if (projType == UPPERTRIPROJ) projectionConstructionUpperTri(projBasis, inBasis, coordSet, m);
}

} // end namespace LatticeTester

#endif
