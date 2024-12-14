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
//#include "latticetester/NTLWrap.h"
#include "latticetester/Util.h"
#include "latticetester/EnumTypes.h"
#include "latticetester/IntLattice.h"
#include "latticetester/Coordinates.h"
#include "latticetester/LLL_FP64.h"
#include "latticetester/LLL_lt.h"

using namespace LatticeTester;
using namespace NTL;

/**
 * \file latticetester/BasisConstruction.h
 *
 * This file offers functions to construct a basis from a set of generating
 * vectors that are not necessarily independent, to construct a triangular basis,
 * to construct the basis for a projection over a given subset of coordinates,
 * and to obtain the \f$m\f$-dual of a given basis.
 *
 * The implementation relies on NTL and uses NTL matrices.
 * When the basis turns out to have fewer rows than columns, some of the functions
 * add implicitly the rescaled unit vectors to the set of generating vectors.
 * In that case, the basis matrix is always square and all the vectors of the form
 * \f$m \be_i\f$ belong to the lattice.
 *
 * NTL already offers an efficient procedure to construct an LLL-reduced basis from a set
 * of generating vectors.  This is the most effective way of constructing a basis
 * and it is encapsulated in the `LLLConstruction0` function given below.
 * This function does not assume that the rescaled unit vectors \f$m \be_i\f$ belong
 * to the lattices and it does not even know about \f$m\f$.
 * The function `LLLBasisConstruction` adds those vectors to the set of generating vectors,
 * so it always returns a square basis.
 *
 * We also offer an alternative functions that construct a triangular basis from a set of
 * generating vectors. They always add the rescaled unit vectors implicitly to the set.
 * The function `lowerTriangularBasis` constructs a lower-triangular basis, while
 * `upperTriangularBasis` constructs an upper-triangular basis.
 *
 * To compute the  \f$m\f$-dual of a given basis, we have a general (but slow) function
 * implemented in `mDualBasis`, and a much faster function in `mDualUpperTriangular`
 * that works only when the basis is upper-triangular.
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
 * All functions in this file are static, so there is no reason to create any
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
 * The programs `BasisManipulationVerbose` and `BasisManipulation` in the examples
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
 * vectors are returned in this array, exactly as in the `LLL_FPZZflex.h` module.
 * These optional parameters are allowed to take non-default values only when
 * `Int==ZZ` and `precision = DOUBLE`.
 * The function returns the dimension (number of rows) of the newly computed basis,
 * which may differ from the number of rows of the `gen` object.
 * The latter is never resized.
 *
 * This function *does not* assume that all vectors \f$m e_i\f$ belong to the lattice, so
 * it may return a basis matrix that has fewer rows than columns!
 * To make sure that these vectors belong to the lattice, we can add them
 * explicitly beforehand to the set of generating vectors, or call the next function.
 */
template<typename Int, typename Real>
static long LLLConstruction0(IntMat &gen, const double delta = 0.9, long r = 0, long c = 0,
      RealVec *sqlen = 0);

/**
 * Similar to `LLLConstruction0`, except that in case the set of generating
 * vectors do not generate a full-dimensional lattice, it adds the vectors
 * \f$m e_i\f$ to the generating set, so it always returns a square matrix.
 * The matrix `gen` is not resized by this function, so it can remain larger
 * than the lattice dimension.
 */
template<typename Int, typename Real>
static void LLLBasisConstruction(IntMat &gen, const Int &m, const double delta = 0.9, long r = 0,
      long c = 0, RealVec *sqlen = 0);

/**
 * Takes a set of generating vectors in the matrix `gen` and iteratively
 * transforms it into a lower triangular lattice basis into the matrix `basis`.
 * This lattice is assumed to contain all the vectors of the form @f$m e_j@f$,
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
static void lowerTriangularBasis(IntMat &gen, IntMat &basis, const Int &m, long r = 0, long c = 0);

// This is a previous version.
template<typename Int>
static void lowerTriangularBasisCW(IntMat &gen, IntMat &basis, const Int &m, long r = 0, long c = 0);

/**
 * Same as `lowerTriangularBasis`, except that the returned basis is upper triangular.
 */
template<typename Int>
static void upperTriangularBasis(IntMat &gen, IntMat &basis, const Int &m, long r = 0, long c = 0);

// This is a previous version.
template<typename Int>
static void upperTriangularBasisCW(IntMat &gen, IntMat &basis, const Int &m, long r = 0, long c = 0);

/**
 * Takes an upper triangular basis matrix `basis` and computes the m-dual basis `basisDual`.
 * The function assumes that each coefficient on the diagonal of `basis` is nonzero and divides `m`.
 * That is, the basis matrix must be square and invertible.
 * The algorithm is described in the Lattice Tester guide \cite iLEC22l.
 * Since the basis is upper triangular, its m-dual will be lower triangular.
 * When `dim > 0`, it gives the number of rows and columns of the matrix `basis`
 * that is actually used. Otherwise (by default) the function uses `basis.numCols()`.
 */
template<typename Int>
static void mDualUpperTriangular(const IntMat &basis, IntMat &basisDual, const Int &m,
      long dim = 0);

/**
 * This function does essentially the same thing as `mDualUpperTriangular`, but the
 * algorithm is slightly different. It uses the method described in \cite rCOU96a.
 */
template<typename Int>
static void mDualUpperTriangular96(IntMat &basis, IntMat &basisDual, const Int &m, long dim = 0);

// static void mDualUpperTriangular96ZZ(NTL::Mat<NTL::ZZ> &basis,
//            NTL::Mat<NTL::ZZ> &basisDual, const NTL::ZZ &m, long dim = 0);

/**
 * This function assumes that `basis` contains a basis of the primal lattice
 * scaled by the factor `m`, not necessarily triangular, and it returns in `basisDual`
 * the m-dual basis. It uses matrix inversion and is rather slow.
 * It is currently implemented only for `Int = ZZ` and it also assumes that the dimensions
 * of the two `IntMat` objects is exactly the same as the dimensions of the lattices.
 * The reason for this is that we use an NTL function that works only under these conditions.
 */
template<typename Int>
static void mDualBasis(const IntMat &basis, IntMat &basisDual, const Int &m);

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
static void projectMatrix(const IntMat &in, IntMat &out, const Coordinates &proj, long r = 0);

/**
 * Constructs a basis for the projection `proj` of the lattice with basis `inBasis`,
 * using `LLLBasisConstruction`, and returns it in `projBasis`.
 * This returned basis is not triangular in general.
 * Its dimension will be the number of coordinates in `proj`.
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
static void projectionConstructionLLL(const IntMat &inBasis, IntMat &projBasis,
      const Coordinates &proj, const Int &m, const double delta = 0.9, long r = 0, RealVec *sqlen =
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
static void projectionConstructionUpperTri(const IntMat &inBasis, IntMat &projBasis,
      IntMat &genTemp, const Coordinates &proj, const Int &m, long r = 0);

template<typename Int>
static void projectionConstructionUpperTri(const IntMat &inBasis, IntMat &projBasis,
      const Coordinates &proj, const Int &m, long r = 0);

/**
 * In this version, the construction method is passed as a parameter. The default is LLL.
 * In the triangular case, a temporary matrix is created internally.
 */
template<typename Int>
static void projectionConstruction(const IntMat &inBasis, IntMat &projBasis,
      const Coordinates &proj, const Int &m, const ProjConstructType projType = LLLPROJ,
      const double delta = 0.9);


//============================================================================
// Implementation

// General case, no implementation.
template<typename Int, typename Real>
static long LLLConstruction0(IntMat &gen, const double delta, long r, long c, RealVec *sqlen) {
   std::cerr << "LLLConstruction0: general case is not implemented.\n";
   exit(1);
}

// The int64_t implementation.
// This one works only for `precision == DOUBLE` and Real == double.
template<>   // <long, double>
long LLLConstruction0(NTL::Mat<long> &gen, const double delta, long r, long c,
      NTL::Vec<double> *sqlen) {
   return NTL::LLL_FP64(gen, delta, r, c, sqlen);
   // return NTL::LLL_FPInt<long>(gen, delta, r, c, sqlen);
}

// The ZZ + double implementation.
template<>   // <NTL::ZZ, double>
long LLLConstruction0(NTL::Mat<NTL::ZZ> &gen, const double delta, long r, long c,
      NTL::Vec<double> *sqlen) {
   return NTL::LLL_FP_lt(gen, delta, r, c, sqlen);
}

// The ZZ + xdouble implementation.
template<>
long LLLConstruction0(NTL::Mat<NTL::ZZ> &gen, const double delta, long r, long c,
      NTL::Vec<xdouble> *sqlen) {
   return NTL::LLL_XD_lt(gen, delta, r, c, sqlen);
}

// The ZZ + quad_float implementation.
template<>
long LLLConstruction0(NTL::Mat<NTL::ZZ> &gen, const double delta, long r, long c,
      NTL::Vec<quad_float> *sqlen) {
   return NTL::LLL_QP_lt(gen, delta, r, c, sqlen);
}

// The ZZ + RR implementation.
template<>
long LLLConstruction0(NTL::Mat<NTL::ZZ> &gen, const double delta, long r, long c,
      NTL::Vec<NTL::RR> *sqlen) {
   return NTL::LLL_RR_lt(gen, delta, r, c, sqlen);
}

//============================================================================

template<typename Int, typename Real>
void LLLBasisConstruction(IntMat &gen, const Int &m, double delta, long r, long c, RealVec *sqlen) {
   //std::cout << "LLLBasisConstruction, before LLL:  c = " << c << "\n";
   int64_t rank = LLLConstruction0<Int, Real>(gen, delta, r, c, sqlen);
   //std::cout << "LLLBasisConstruction, after LLL:  c = " << c << ", rank = " << rank << "\n";
   //std::cout << "Basis: \n" << gen << "\n";
   if (c == 0) c = gen.NumCols();
   if (rank == c) return;  // We are done!

   // We now add the m e_i row vectors, and we redo the LLL with that.
   // No change in the dimensions of gen.
   int64_t i, j;
   for (i = rank; i < rank + c; i++) {
      for (j = 0; j < c; j++) {
         if (i == j) gen[i][j] = m;
         else gen[i][j] = 0;
      }
   }
   std::cout << "Warning for LLLBasisConstruction: we had to add some rows!\n";
   std::cout << "  c = " << c << ", rank = " << rank << "\n";
   rank = LLLConstruction0<Int, Real>(gen, delta, rank + c, c, sqlen);
}

//==============================================================================

// This one does not seem to work... It sometimes gives wrong results.
template<typename Int>
void lowerTriangularBasisCW(IntMat &gen, IntMat &basis, const Int &m, long dim1, long dim2) {
   // Note:  dim1 = r, dim2 = c.   The new basis should be c x c.
   NTL::Vec<Int> coeff_gcd, coeff_xj, xj; // Vectors are created locally here.
   Int gcd, gcdCopy, C, D;
   long i, j, k, l;
   // In case r or c is zero:
   if (dim1 == 0) dim1 = gen.NumRows();
   if (dim2 == 0) dim2 = gen.NumCols();
   assert(dim2 <= basis.NumRows() && dim2 <= basis.NumCols());
   // Allocate space for the vectors:
   coeff_gcd.SetLength(dim1);
   coeff_xj.SetLength(dim1);
   xj.SetLength(dim2);
   // if (basis.NumRows() != dim2 || basis.NumCols() != dim2)
   //    basis.SetDims(dim2, dim2);

   for (i = dim2 - 1; i > -1; i--) {
      // Reset these vectors to 0, as they may contain nonzero values from the previous i.
      for (j = dim1 - 1; j > -1; j--)
         coeff_gcd[j] = 0;
      for (j = dim2 - 1; j > -1; j--)
         xj[j] = 0;
      // Search for the first non-zero element in the row.
      for (k = dim1 - 1; (k > -1 && gen[dim1 - 1 - k][i] == 0); k--) {
      }
      //			if (gen[k][i] != 0)	break;
      // Reduce the other generators as they are used often in what follows.
      for (j = dim1 - 1; j > dim1 - k - 1; j--)
         NTL::rem(gen[j][i], gen[j][i], m);
      // The `else` case adds m e_i to the basis matrix.
      if (k > -1) {
         gcd = m;    // Will be GCD(m, gen[k][i]);
         coeff_gcd[k] = 1;
         gcdCopy = gcd;

         // Find the other coefficients by applying the Euclidean algorithm multiple times
         for (j = dim1 - 1; j > dim1 - k - 1; j--) {
            if (gen[j][i] == 0) coeff_gcd[j] = 0;
            else {
               NTL::XGCD(gcd, C, D, gcdCopy, gen[j][i]);
               coeff_gcd[j] = D;
               for (l = dim1 - j - 1 - 1; l > -1; l--) {
                  NTL::mul(coeff_gcd[dim1 - 1 - l], coeff_gcd[dim1 - 1 - l], C);
               }
               gcdCopy = gcd;
            }
         }
         // If gcd = m, then this basis (row) vector will be `m e_i`.
         if (gcd == m) {
            for (j = dim2 - 1; j > -1; j--) {
               if (j != i) basis[i][j] = 0;
               else basis[i][j] = m;
            }
         } else {
            // Reduce the coefficients found during the Euclidean algorithm.
            for (j = 0; j < dim1; j++)
               NTL::rem(coeff_gcd[dim1 - 1 - j], coeff_gcd[dim1 - 1 - j], m);
            // We have now found all the coefficients and can compute the vector x_i.
            for (l = dim1 - 1; l > -1; l--) {
               if (coeff_gcd[l] != 0) {
                  for (j = dim2 - 1; j > dim1 - 1 - i - 1; j--) {
                     NTL::MulAddTo(xj[j], gen[l][j], coeff_gcd[l]);
                  }
               }
            }
            // Next we calculate the new vectors v_i.
            // We first calculate the coefficients with which x_i needs to be multiplied.
            for (j = dim1 - 1; j > -1; j--) {
               NTL::div(coeff_xj[j], gen[j][i], gcd);
               NTL::rem(coeff_xj[j], coeff_xj[j], m);
            }
            for (j = dim2 - 1; j > -1; j--)
               NTL::rem(xj[j], xj[j], m);
            // Update the v_i
            for (l = dim1 - 1; l > -1; l--) {
               if (coeff_xj[l] != 0) {
                  for (j = dim2 - 1; j > dim1 - 1 - i - 1; j--) {
                     NTL::MulSubFrom(gen[l][j], coeff_xj[l], xj[j]);
                  }
               }
            }
            // Set the `i`th base vector.
            for (j = 0; j < dim2; j++)
               basis[i][j] = xj[j];
         }
      } else {
         for (j = dim2 - 1; j > -1; j--) {
            if (j != i) basis[i][j] = 0;
            else basis[i][j] = m;
         }
      }
   }
}


//===================================================================

template<typename Int>
void upperTriangularBasisCW(IntMat &gen, IntMat &basis, const Int &m, long dim1, long dim2) {
   // `dim1` and `dim2` are `s` and `t` in the guide.
   if (dim1 == 0) dim1 = gen.NumRows();
   if (dim2 == 0) dim2 = gen.NumCols();
   assert(dim2 <= basis.NumRows() && dim2 <= basis.NumCols());

   long i, j, k, l;
   Int gcd, gcdCopy, C, D;
   NTL::Vec<Int> coeff_gcd, coeff_xj, xj;  // Here we create new vectors!
   coeff_gcd.SetLength(dim1);
   coeff_xj.SetLength(dim1);
   xj.SetLength(dim2);

   for (j = 0; j < dim2; j++) {
      for (i = 0; i < dim1; i++)
         coeff_gcd[i] = 0;
      for (i = 0; i < dim2; i++)
         xj[i] = 0;
      // Search for the first non-zero element in column j.
      for (k = 0; (k < dim1 && gen[k][j] == 0); k++) {}
      // The reduction of the generating vectors is done only here, for column j.
      // Reduce the other generators as they are used often in what follows.
      for (i = k; i < dim1; i++) {
         NTL::rem(gen[i][j], gen[i][j], m);  // This should be done *before* searching for the first nonzero element?
      }
      // The `else` case adds m e_j to the basis matrix.
      if (k < dim1) {
         gcd = m;    // Will be GCD(m, gen[k][j]);
         coeff_gcd[k] = 1;
         gcdCopy = gcd;

         // Find the other coefficients by applying the Euclidean algorithm multiple times
         for (i = k; i < dim1; i++) {
            if
               (gen[i][j] == 0) coeff_gcd[i] = 0;
            else {
               // XGCD (g, c, d, a, b) does g = gcd(a, b) = a*c + b*d.
               NTL::XGCD(gcd, C, D, gcdCopy, gen[i][j]);
               coeff_gcd[i] = D;
               for (l = 0; l < i; l++) {
                  NTL::mul(coeff_gcd[l], coeff_gcd[l], C);
               }
               gcdCopy = gcd;
            }
         }
         // If gcd = m, then this basis (row) vector will be `m e_j`.
         if (gcd == m) {
            for (i = 0; i < dim2; i++) {
               if (i != j) basis[j][i] = 0;
               else basis[j][i] = m;
            }
         } else {
            // Reduce the coefficients found during the Euclidean algorithm.
            for (i = 0; i < dim1; i++) {
               NTL::rem(coeff_gcd[i], coeff_gcd[i], m);
            }
            // We have now found all the coefficients and can compute the vector x_j.
            for (k = 0; k < dim1; k++) {
               if (coeff_gcd[k] != 0) {
                  for (i = j; i < dim2; i++) {
                     NTL::MulAddTo(xj[i], gen[k][i], coeff_gcd[k]);
                     NTL::rem(xj[i], xj[i], m);
                  }
               }
            }
            // Next we calculate the new vectors v_j.
            // We first calculate the coefficients with which x_j needs to be multiplied.
            for (i = 0; i < dim1; i++) {
               NTL::div(coeff_xj[i], gen[i][j], gcd);
               NTL::rem(coeff_xj[i], coeff_xj[i], m);
            }
            for (i = 0; i < dim2; i++)
               NTL::rem(xj[i], xj[i], m);
            // Update the v_i's
            for (k = 0; k < dim1; k++) {
               if (coeff_xj[k] != 0) {
                  for (i = j; i < dim2; i++) {
                     NTL::MulSubFrom(gen[k][i], coeff_xj[k], xj[i]);
                  }
               }
            }
            // Set the `j`th basis vector.
            for (i = 0; i < dim2; i++)
               basis[j][i] = xj[i];
         }
      } else {
         for (i = 0; i < dim2; i++) {
            if (j != i) basis[j][i] = 0;
            else basis[j][i] = m;
         }
      }
   }
}

//===================================================================

template<typename Int>
void lowerTriangularBasis(IntMat &gen, IntMat &basis, const Int &m, long dim1, long dim2) {
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

   for (j = dim2-1; j >= 0; j--) {  // column j.
      // Here we compute the submatrix whose upper left corner is (j,j) in the upper triangular basis.
      // Find c_j and the coefficients `a_{i,j}` by applying the Euclidean algorithm multiple times.
      for (i = 0; i < dim1; i++) {
         NTL::rem(gen[i][j], gen[i][j], m);  // Reduce modulo m
         coeff_gcd[i] = 0;
      }
      gcd = m;
      for (i = dim1-1; i >= 0; i--) {
         // NTL::rem(gen[i][j], gen[i][j], m);
         if (gen[i][j] != 0) {
            // XGCD (g, c, d, const a, const b) does g = gcd(a, b) = a*c + b*d.
            gcdCopy = gcd;  // We need a copy for the `const a` parameter.
            NTL::XGCD(gcd, c, d, gcdCopy, gen[i][j]);
            coeff_gcd[i] = d;
            for (l = i+1; l < dim1; l++) {
               NTL::mul(coeff_gcd[l], coeff_gcd[l], c);
               NTL::rem(coeff_gcd[l], coeff_gcd[l], m);
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
                  std::cout << "  Before: xj[k] = " << xj[k] << "\n";
                  NTL::rem(xj[k], xj[k], m);
                  std::cout << "  After: xj[k] = " << xj[k] << "\n";
                  // if (xj[k] < 0) NTL::add(xj[k], xj[k], m);
               }
            }
         }
         std::cout << "  j = " << j << "\n";
         std::cout << "  vector x_j = " << xj << "\n";
         // Next we update the vectors v_i.
         // We first calculate the coefficients with which x_j needs to be multiplied.
         for (i = 0; i < dim1; i++) {
            NTL::div(coeff_xj[i], gen[i][j], gcd);
            // std::cout << "  coeff_xj[i] = " << coeff_xj[i] << "\n";
            // NTL::rem(coeff_xj[i], coeff_xj[i], m);
            // std::cout << "  coeff_xj[i] = " << coeff_xj[i] << "\n";
         }
         std::cout << "  vector x_j = " << xj << "\n";
         std::cout << "  coeff_x_j = " << coeff_xj << "\n";
         // Update the components of index <= j of the old vectors v_i.
         for (i = 0; i < dim1; i++) {
            if (coeff_xj[i] != 0) {
               for (k = j; k >= 0 ; k--) {
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
void upperTriangularBasis(IntMat &gen, IntMat &basis, const Int &m, long dim1, long dim2) {
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
               NTL::rem(coeff_gcd[l], coeff_gcd[l], m);
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
                  NTL::rem(xj[k], xj[k], m);
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

/**
 * For `B` to be `m`-dual to `A`, we have to have that \f$AB^t = mI\f$.
 * Since `A` is upper triangular, `B` will be a lower triangular matrix
 * with `A(i,i)*B(i,i) = m` for all `i` and
 * \f$ A_i \cdot B_j = 0\f$ for \f$i\neq j\f$. To get the second condition,
 * we simply have to recursively take for each line
 * \f[B_{i,j} = -\frac{1}{A_{j,j}}\sum_{k=j+1}^i A_{j,k} B_{i,k}.\f]
 */
template<typename Int>
void mDualUpperTriangular(const IntMat &A, IntMat &B, const Int &m, long dim) {
   // Note:  A = basis,  B = basisDual
   if (dim == 0) dim = A.NumRows();
   assert(dim <= A.NumCols());
   assert(dim <= B.NumRows() && dim <= B.NumCols());
   for (int64_t i = 0; i < dim; i++) {
      // Put zeros under the diagonal.
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

/*
 template<>
 void mDualUpperTriangular(const NTL::Mat<long> &A, NTL::Mat<long> &B,
 const long &m, long dim) {
 // Note:  A = basis,  B = basisDual
 if (dim == 0)
 dim = A.NumRows();
 assert(dim <= A.NumCols());
 assert(dim <= B.NumRows() && dim <= B.NumCols());
 for (int64_t i = 0; i < dim; i++) {
 // Put zeros under the diagonal.
 for (int64_t j = i + 1; j < dim; j++)
 B[i][j] = 0;
 // Set diagonal elements.
 B[i][i] = m / A[i][i];
 // Compute the other ones.
 for (int64_t j = i - 1; j >= 0; j--) {
 B[i][j] = 0;
 for (int64_t k = j + 1; k <= i; k++)
 NTL::MulSubFrom(B[i][j], A[j][k], B[i][k]);
 NTL::div(B[i][j], B[i][j], A[j][j]);
 }
 }
 }
 */

//======================================================
// This is the old version from Couture and L'Ecuyer (1996).
template<typename Int>
void mDualUpperTriangular96(IntMat &basis, IntMat &basisDual, const Int &m, long dim) {
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
            basisDual *= basis[j][j] / gcd;
         }
         DivideRound(basisDual[i][j], basis[j][j], basisDual[i][j]);
      }
   }
}

// A specialized implementation for ZZ. What is different from Int version?   ****
template<>
void mDualUpperTriangular96(NTL::Mat<NTL::ZZ> &basis, NTL::Mat<NTL::ZZ> &basisDual,
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
            gcd = NTL::GCD(basisDual(i, j), basis[j][j]);
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
void mDualBasis(const NTL::Mat<Int> &basis, NTL::Mat<Int> &basisDual, const Int &m) {
   std::cerr << "mDualBasis is implemented only for NTL::ZZ integers.\n";
   exit(1);
}

// The specialization for the case where `Int = ZZ`.
template<>
// void mDualBasis(const IntMat &basis, IntMat &basisDual,
void mDualBasis(const NTL::Mat<NTL::ZZ> &basis, NTL::Mat<NTL::ZZ> &basisDual,
      const NTL::ZZ &m) {
   NTL::ZZ det, fac;
   long dim = basis.NumRows();
   if (dim != basis.NumCols()) {
      std::cerr << "mDualBasis: the given basis matrix must be square.\n";
      exit(1);
   }
   inv(det, basisDual, basis);
   NTL::Mat<NTL::ZZ> C = basisDual;
   div(fac, det, m);
   for (int64_t i = 0; i < dim; i++) {
      for (int64_t j = 0; j < dim; j++) {
         div(basisDual[j][i], C[i][j], fac);
      }
   }
}

/*
 template<>
 void BasisConstruction<NTL::ZZ>::mDualBasis(const NTL::Mat<NTL::ZZ> &basis,
 NTL::Mat<NTL::ZZ> &basisDual, const NTL::ZZ &m, long dim) {
 NTL::ZZ d, fac;
 if (!dim) {
 dim = basis.NumRows();
 basisDual.SetDims(dim, dim);
 }
 IntMat copyBasis, copyBasisDual;
 copyBasis.SetDims(dim, dim);
 // We added the dim param. just to avoid this type of creation of new objects,
 // so this creation and copying defeats the purpose and seems to be nonsense!!!
 // But the function `inv` is very slow anyway.

 // Here we copy only the part that we need.
 for (int64_t i = 0; i < dim; i++) {
 for (int64_t j = 0; j < dim; j++)
 copyBasis[i][j] = basis[i][j];
 }
 inv(d, copyBasisDual, copyBasis);
 NTL::Mat<NTL::ZZ> C = copyBasisDual;
 div(fac, d, m);
 for (int64_t i = 0; i < basis.NumRows(); i++) {
 for (int64_t j = 0; j < basis.NumCols(); j++) {
 if (i < dim && j < dim)
 div(basisDual[j][i], C[i][j], fac);
 else
 basisDual[j][i] = 0;
 }
 }
 }
 */

//=================================================================================
template<typename Int>
void projectMatrix(const IntMat &in, IntMat &out, const Coordinates &proj, long r) {
   if (in == out) {
      std::cout << "\n***** Error: in and out must be different IntMat objects " << std::endl;
      exit (1);
   }
   if (!r) r = in.NumRows();   // In case r=0.
   // We assume without testing that `out` is large enough for proj.size().
   long j = 0;
   for (auto it = proj.begin(); it != proj.end(); it++, j++) {
      for (long i = 0; i < r; i++)
         out[i][j] = in[i][*it - 1];
   }
}

//===================================================
template<typename Int, typename Real>
void projectionConstructionLLL(const IntMat &inBasis, IntMat &projBasis, const Coordinates &proj,
      const Int &m, const double delta, long r, RealVec *sqlen) {
   projectMatrix(inBasis, projBasis, proj, r);
   LLLBasisConstruction(projBasis, m, delta, r, proj.size(), sqlen);
}

//===================================================

template<typename Int>
void projectionConstructionUpperTri(const IntMat &inBasis, IntMat &projBasis, IntMat &genTemp,
      const Coordinates &proj, const Int &m, long r) {
   projectMatrix(inBasis, genTemp, proj, r);
   upperTriangularBasis(genTemp, projBasis, m, r, proj.size());
}

template<typename Int>
void projectionConstructionUpperTri(const IntMat &inBasis, IntMat &projBasis,
      const Coordinates &proj, const Int &m, long r) {
   long dim = proj.size();      // Dimension of projection
   if (!r) r = inBasis.NumRows();
   IntMat genTemp;
   genTemp.SetDims(r, dim); // Here an internal object is created and resized!
   projectMatrix(inBasis, genTemp, proj, r);
   upperTriangularBasis(genTemp, projBasis, m, r, dim);
}

//===================================================

template<typename Int>
void projectionConstruction(const IntMat &inBasis, IntMat &projBasis, const Coordinates &proj,
      const Int &m, const ProjConstructType projType, const double delta) {
   if (projType == LLLPROJ) projectionConstructionLLL(inBasis, projBasis, proj, m, 0, delta);
   if (projType == UPPERTRIPROJ) projectionConstructionUpperTri(inBasis, projBasis, proj, m);
}

} // end namespace LatticeTester

#endif
