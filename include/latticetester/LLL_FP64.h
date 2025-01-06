//
// This file is part of LatticeTester.
//
#ifndef NTL_LLL_FP64__H
#define NTL_LLL_FP64__H

#include <NTL/vector.h>
#include <NTL/matrix.h>
//#include <latticetester/Util.h>

using namespace NTL;

/**
 * \file LLL_FP64.h
 * This file contains a modified version of the `LLL_FP` module of NTL.
 * With the modified functions, the basis entries are in `long`, alias `int64_t`.
 * We can also apply LLL or BKZ only to a submatrix
 * (first `r` rows and `c` columns) of the matrix `B` that is passed in and returned.
 * The returned basis will have `c` columns and at most `max(r,c)` rows
 * (the rank of the basis matrix), so it may not occupy the entire space in `B`.
 * With this flexibility, we can reserve a large block of
 * memory for the matrix `B` and reuse this same block (the same object `B`)
 * for thousands or millions of lattices that we want to analyze, even if the bases
 * have different dimensions.
 *
 * Another addition is the possibility to recover an array `sqlen` that gives the square
 * Euclidean lengths of the basis vectors, in `double`.
 * This array is maintained in the `LLL_FP` functions of NTL,
 * but it is hidden in the implementation and not accessible from outside.
 * In `IntLattice`, these lengths are maintained in a `RealVec` object.
 * In this module, we assume that `Real = double`.
 *
 * Each function returns the dimension of the computed basis (number of independent rows).
 * This basis is always returned in the upper-left corner of the matrix `B`.
 * This differs from the `LLL_FP` functions, which returns the zero vectors at the top.
 */

/**
 * These functions are in the NTL namespace.
 */
namespace NTL {

/**
 * This function is similar to `LLL_FP` in NTL, but only the first `r` rows and
 * first `c` columns of the matrix `B` are considered, and a basis is built
 * for the lattice generated by these (partial) rows.  The other elements of
 * `B` are ignored. The basis is returned in the upper left corner of `B`,
 * with the shortest basis vector always in the first row.
 * If `r=0`, then all the rows of the `IntMat` object are taken.
 * If `c=0`, then all the columns are taken.
 * The square lengths of the returned basis vectors are also returned in the
 * `double` vector `sqlen`,  in `sqlen[0],..., sqlen[d-1]`, if this vector given.
 * The indices of `B` and `sqlen` start at 0.
 * The function returns the dimension of the computed basis (the number of independent rows).
 */
long LLL_FP64(NTL::Mat<long> &B, const double delta = 0.99,
      long r = 0, long c = 0, NTL::Vec<double> *sqlen = 0);

/**
 * This function is similar to `BKZ_FP` in NTL, with the same modifications
 * as in LLL_FP64 above.
 */
long BKZ_FP64(NTL::Mat<long> &BB, const double delta = 0.99, long blocksize = 10,
      long prune = 0, long r = 0, long c = 0, NTL::Vec<double> *sqlen = 0);

}

#endif
