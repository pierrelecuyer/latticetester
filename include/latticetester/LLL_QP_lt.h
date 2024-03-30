#ifndef NTL_LLL_QP_lt__H
#define NTL_LLL_QP_lt__H

#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/mat_ZZ.h>
#include <NTL/vec_double.h>
#include <NTL/vec_quad_float.h>
// #include <latticetester/FlexTypes.h>

NTL_OPEN_NNS

// classical Gramm-Schmidt versions. modified for Lattice Tester


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
long LLL_QP_lt(mat_ZZ &B, double delta = 0.99,
      long r = 0, long c = 0, vec_quad_float *sqlen = 0);

/**
 * This function is similar to `BKZ_FP` in NTL, with the same modifications
 * as in LLL_FPInt above.
 */
long BKZ_QP_lt(mat_ZZ &BB, double delta = 0.99,
      long blocksize = 10, long prune = 0, long r = 0, long c = 0, vec_quad_float *sqlen = 0);


NTL_CLOSE_NNS

#endif
