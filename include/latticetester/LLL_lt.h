#ifndef NTL_LLL_lt__H
#define NTL_LLL_lt__H

#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/mat_ZZ.h>
#include <NTL/vec_double.h>
#include <NTL/vec_xdouble.h>
#include <NTL/vec_quad_float.h>
#include <NTL/vec_RR.h>

/**
 * \file LLL_lt.h
 *
 * These static functions are slightly modified versions of the `LLL` and `BKZ`
 * functions of %NTL.
 *
 * With the modifications, one can apply the reductions only to an
 * upper-left corner of the basis matrix `B`, the first `r` rows and
 * first `c` columns.  The other elements of `B` are ignored.
 * The basis is returned in the same upper left corner of `B`,
 * with the shortest basis vector always in the first row.
 * This differs from the `LLL_FP` functions, which returns the zero vectors at the top.
 * If `r=0`, then all the rows of the `IntMat` object are taken.
 * If `c=0`, then all the columns are taken.
 * With this added flexibility the same block of memory can be reused
 * for millions of lattices even if their bases have different dimensions.
 * The square lengths of the returned basis vectors are also returned in the
 * `double` vector `sqlen`,  in `sqlen[0],..., sqlen[d-1]`, if this vector given.
 * The indices of `B` and `sqlen` start at 0.
 * Each function returns the dimension of the computed basis (number of independent rows).
 */

using namespace NTL;

/**
 * These functions are in the %NTL namespace.
 */
namespace NTL {

/**
 * This function is similar to `LLL_FP` in %NTL, but only the first `r` rows and
 * first `c` columns of the matrix `B` are considered, a basis is built
 * for the lattice generated by these (partial) rows,
 * and the square lengths of the returned basis vectors are returned in `sqlen`.
 */
long LLL_FP_lt(mat_ZZ &B, double delta = 0.99,
      long r = 0, long c = 0, vec_double *sqlen = 0);

/**
 * This function is similar to `BKZ_FP` in %NTL, with the same modifications
 * as in LLL_FPInt above.
 */
long BKZ_FP_lt(mat_ZZ &BB, double delta = 0.99,
      long blocksize = 10, long prune = 0, long r = 0, long c = 0, vec_double *sqlen = 0);


long LLL_XD_lt(mat_ZZ &B, double delta = 0.99,
      long r = 0, long c = 0, vec_xdouble *sqlen = 0);

/**
 * This function is similar to `BKZ_XD` in %NTL, with the same modifications
 * as in LLL_FPInt above.
 */
long BKZ_XD_lt(mat_ZZ &BB, double delta = 0.99,
      long blocksize = 10, long prune = 0, long r = 0, long c = 0, vec_xdouble *sqlen = 0);


long LLL_QP_lt(mat_ZZ &B, double delta = 0.99,
      long r = 0, long c = 0, vec_quad_float *sqlen = 0);

/**
 * This function is similar to `BKZ_FP` in %NTL, with the same modifications
 * as in LLL_FPInt above.
 */
long BKZ_QP_lt(mat_ZZ &BB, double delta = 0.99,
      long blocksize = 10, long prune = 0, long r = 0, long c = 0, vec_quad_float *sqlen = 0);


long LLL_RR_lt(mat_ZZ& B, double delta = 0.99, long r = 0, long c = 0, vec_RR* sqlen = 0);

/**
 * These two functions are wrappers of `BKZ_RR` in %NTL, with the same modifications
 * as in `LLL_RR_lt` above.
 */
long BKZ_RR_lt(mat_ZZ& B, double delta = 0.99,
      long blocksize = 10, long prune = 0, long r = 0, long c = 0, vec_RR* sqlen = 0);

}

#endif
