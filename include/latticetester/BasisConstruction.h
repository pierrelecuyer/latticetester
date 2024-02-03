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
#include <type_traits>

#include "NTL/tools.h"
#include "NTL/ZZ.h"
#include "NTL/vec_ZZ.h"
#include "NTL/mat_ZZ.h"
#include "NTL/LLL.h"

#include "latticetester/NTLWrap.h"
#include "latticetester/EnumTypes.h"
#include "latticetester/IntLattice.h"
#include "latticetester/Util.h"
#include "latticetester/Coordinates.h"
#include "latticetester/LLL64.h"
#include "latticetester/LLL_FPZZflex.h"

namespace LatticeTester {

/**
 * This static class offers methods (functions) to construct a basis from a set of generating
 * vectors that are not necessarily independent, to construct a triangular basis,
 * to construct the basis for a projection over a given subset of coordinates,
 * and to obtain the \f$m\f$-dual of a given basis.
 * The implementation relies on NTL and uses NTL matrices.
 * When the basis turns out to have fewer rows than columns, some of the methods
 * add implicitly the rescaled unit vectors to the set of generating vectors.
 * In that case, the basis matrix is always square and all the vectors of the form
 * \f$m \be_i\f$ belong to the lattice.
 *
 * NTL already offers an efficient method to construct an LLL-reduced basis from a set
 * of generating vectors.  This is the most effective way of constructing a basis
 * and it is encapsulated in the `LLLConstruction0` method given below.
 * This method does not assume that the rescaled unit vectors \f$m \be_i\f$ belong
 * to the lattices and it does not even know about \f$m\f$.
 * The method `LLLBasisConstruction` adds those vectors to the set of generating vectors,
 * so it always returns a square basis.
 * ****   We must make sure that it works in `dim < maxDim` dimensions as well!   ****
 * ****   Currently, it seems that this works only when `Int==ZZ` and `precision = DOUBLE`.
 *
 * We also offer an alternative methods that construct a triangular basis from a set of
 * generating vectors. They always add the rescaled unit vectors implicitly to the set.
 * The method `lowerTriangularBasis` constructs a lower-triangular basis, while
 * `upperTriangularBasis` constructs an upper-triangular basis.
 *
 * To compute the  \f$m\f$-dual of a given basis, we have a general method implemented in
 * `mDualBasis`, and a much faster method in `mDualUpperTriangular` that works only when
 * the basis is upper-triangular.
 *
 * We also have functions to compute the basis of a projection of a given lattice over
 * a specified set of coordinates.  The function `projectionConstructionLLL` does this
 * by using LLL to construct the basis of the projection, while `projectionConstructionUpperTri`
 * constructs an upper-triangular basis for the projection.
 * The function `projectionConstruction` takes the construction method as a parameter.
 *
 * All functions take as input a `IntMat` object that contains either a basis or a set of
 * generating vectors. By default, all rows and columns of this matrix object are used.
 * But in some functions, the user can ask to use only a subset of these rows and columns,
 * via the optional parameters `r` and `c`.  This can permit one to use the same `IntMat`
 * object for several numbers of dimensions, to avoid object creations.
 *
 * All functions in this class are static, so there is no reason to create any
 * `BasisConstruction` object. We also avoid to create new objects (such as vectors and
 * matrices) inside these functions.  These functions can be called thousands or millions
 * of times in a program, and we want the user to be able to re-use the same vectors and
 * matrices over and over again instead of creating new ones.
 *
 * The programs `BasisManipulationVerbose` and `BasisManipulation` in the examples
 * illustrate how to use these functions and make speed comparisons.
 *
 */

template<typename Int> class BasisConstruction {

private:
    typedef NTL::vector<Int> IntVec;
    typedef NTL::matrix<Int> IntMat;
    // typedef NTL::matrix<int64_t> mat_long;
    // typedef NTL::vector<int64_t> vec_long;
public:

    /**
     * This function takes a set of generating vectors of a lattice in matrix `gen` and
     * finds a lattice basis by applying LLL reduction with the given value of `delta`,
     * using the NTL implementation specified by `prec`.
     * The basis is returned in the first rows of `gen`, in a number of rows equal to its rank.
     * See the class `EnumTypes` and the documentation of LLL in NTL for the meaning and
     * choices for `prec`. The default value of delta is not very close to 1 because the
     * goal here is just to compute a basis. It can be taken much closer to 1 if we
     * really prefer a highly-reduced basis.
     * The optional parameters `r` and `c` indicate the numbers of rows and columns
     * of the matrix `gen` that are actually used. When they are 0 (the default values),
     * then these numbers are taken to be the dimensions of the matrix `gen`.
     * When `sqlen` is not 0, the square lengths of the basis
     * vectors are returned in this array, exactly as in the `LLL_FPZZflex.h` module.
     * These optional parameters are allowed only when `Int==ZZ` and `precision = DOUBLE`.
     * The function returns the dimension (number of rows) of the newly computed basis,
     * which may differ from the number of rows of the `gen` object.
     * The latter is never resized.
     *
     * This function *does not* assume that all vectors \f$m e_i\f$ belong to the lattice, so
     * it may return a basis matrix that has fewer rows than columns!
     * To make sure that these vectors belong to the lattice, we can add them
     * explicitly beforehand to the set of generating vectors, or call the next method.
     */
    static long LLLConstruction0(IntMat &gen, double delta = 0.9, long r = 0,
            long c = 0, double *sqlen = 0, PrecisionType precision = DOUBLE);

    /**
     * Similar to `LLLConstruction0`, except that in case the set of generating
     * vectors do not generate a full-dimensional lattice, it adds the vectors
     * \f$m e_i\f$ to the generating set, so it always returns a square matrix.
     * The matrix `gen` is not resized by this function, so it can remain larger
     * than the lattice dimension.
     */
    static void LLLBasisConstruction(IntMat &gen, const Int &m, double delta =
            0.9, long r = 0, long c = 0, double *sqlen = 0,
            PrecisionType precision = DOUBLE);

    /**
     * Takes a set of generating vectors in the matrix `gen` and iteratively
     * transforms it into a lower triangular lattice basis into the matrix `basis`.
     * This lattice is assumed to contain all the vectors of the form @f$m e_j@f$,
     * so these vectors are added implicitly to the generating set.
     * Apart from that, all the entries of `gen` given as input are assumed to be
     * reduced modulo the scaling factor `m` and all the computations are done modulo `m`.
     * The matrix `basis` is resized automatically (internally) to a square matrix with
     * the appropriate dimensions.
     * After the execution, `gen` will contain irrelevant information (garbage)
     * and `basis` will contain an upper triangular basis.
     * The algorithm is explained in the lattice tester guide.
     * Important: `gen` and `basis` must be different `IntMat` objects.
     *
     * When `r` and/or `c` are strictly positive, they specify the numbers of
     * rows and columns of `gen` that are actually used, as in `LLLConstruction0`.
     * The matrix `gen` is never resized.
     */
    static void lowerTriangularBasis(IntMat &gen, IntMat &basis, const Int &m,
            long r = 0, long c = 0);

    /**
     * Same as `lowerTriangularBasis`, except that the returned basis is upper triangular.
     */
    static void upperTriangularBasis(IntMat &gen, IntMat &basis, const Int &m,
            long r = 0, long c = 0);

    /**
     * Takes an upper triangular basis matrix `basis` and computes the m-dual basis `basisDual`.
     * The method assumes that each coefficient on the diagonal of `basis` is nonzero and divides `m`.
     * That is, the basis matrix must be square and invertible.
     * The algorithm is described in the Lattice Tester guide \cite iLEC22l.
     * Since the basis is upper triangular, its m-dual will be lower triangular.
     * When `dim > 0`, it gives the number of rows and columns of the matrix `basis`
     * that is actually used. Otherwise (by default) the method uses `basis.numCols()`.
     */
    static void mDualUpperTriangular(const IntMat &basis, IntMat &basisDual,
            const Int &m, long dim = 0);

    /**
     * This function does essentially the same thing as `mDualUpperTriangular`, but the
     * algorithm is slightly different. It uses the method described in \cite rCOU96a.
     */
    static void mDualUpperTriangular96(IntMat &basis, IntMat &basisDual,
            const Int &m, long dim = 0);

    static void mDualUpperTriangular96ZZ(NTL::matrix<NTL::ZZ> &basis,
            NTL::matrix<NTL::ZZ> &basisDual, const NTL::ZZ &m, long dim = 0);

    /**
     * This function assumes that `basis` contains a basis of the primal lattice
     * scaled by the factor `m`, not necessarily triangular, and it returns in `basisDual`
     * the m-dual basis.  It is currently implemented only for `Int = ZZ`.
     * It used matrix inversion and is rather slow.
     * ****   This function currently uses always `dim = basis.NumRows()`.  Not good!
     */
    static void mDualBasis(const IntMat &basis, IntMat &basisDual,
            const Int &m);

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
    static void projectMatrix(const IntMat &in, IntMat &out,
            const Coordinates &proj, long r = 0);

    /**
     * Same like projectMatrix but the number of rows is equal to the number of projection coordinates.
     *
     * *****  What is this used for?  In general, this will not generate a basis.
     */
    //static void projectMatrix_minRows(const IntMat &in, IntMat &out,
    //      const Coordinates &proj);
    /**
     * ******   If only for an LCG, does not belong here!   ******
     * This function directly calculates the projection of the dual matrix for an LCG    ******
     * if the basis itself is given in the standard form.
     */
    // static void projectMatrixDual(const IntMat &in, IntMat &out, const Coordinates &proj);
    /**
     * Constructs a basis for the projection `proj` of the lattice with basis `inBasis`,
     * using `LLLBasisConstruction`, and puts it in `projBasis`.
     * This basis is not triangular in general.
     * This will overwrite the matrix `projBasis`, or part of it.
     * It does not compute the dual.
     * When `r > 0`, only the first `r` rows of the matrix `inBasis` are actually used.
     * This `r` should normally be the true dimension `dim` of that basis, which is often
     * smaller than the size `maxDim` of the `IntMat` object that contains the basis.
     * The square Euclidean lengths of the basis vectors are returned in the array `sqlen` when
     * the latter is given.
     */
    static void projectionConstructionLLL(const IntMat &inBasis,
            IntMat &projBasis, const Coordinates &proj, const Int &m,
            const double delta = 0.9, long r = 0, double *sqlen = 0);

    /**
     * Same as `projectionConstructionLLL`, but the construction is made using
     * `upperTriangularBasis`, so the returned basis is upper triangular.
     * In the first version, we pass a matrix `genTemp` that will be used to store the
     * generating vectors of the projection before making the triangularization.
     * We pass it as a parameter to avoid the internal creation of a new matrix each time.
     * In the second version, this matrix is not passed and a temporary one is created internally.
     */
    static void projectionConstructionUpperTri(const IntMat &inBasis,
            IntMat &projBasis, IntMat &genTemp, const Coordinates &proj,
            const Int &m, long r = 0);

    static void projectionConstructionUpperTri(const IntMat &inBasis,
            IntMat &projBasis, const Coordinates &proj, const Int &m,
            long r = 0);

    /**
     * In this version, the construction method is passed as a parameter. The default is LLL.
     * In the triangular case, a temporary matrix is created internally.
     */
    static void projectionConstruction(const IntMat &inBasis, IntMat &projBasis,
            const Coordinates &proj, const Int &m,
            const ProjConstructType projType = LLLPROJ,
            const double delta = 0.9);

};

//============================================================================
// Implementation

// The int64_t implementation.
// This one works only for r = c = 0 (that is, dim == maxDim) for now!
template<>
long BasisConstruction<int64_t>::LLLConstruction0(NTL::matrix<int64_t> &gen,
        double delta, long r, long c, double *sqlen, PrecisionType precision) {
    long num = gen.NumRows();   // Number of generating vectors.
    int64_t rank = 0;
    if ((precision == DOUBLE) & (c == 0) & (r == 0) & (sqlen == 0))
        rank = LLL64_FP(gen, delta);
    // rank = LLL_FP64flex(gen, delta, r, c, sqlen);    //  TO DO !!!!
    else
        std::cerr
                << "LLLConstruction0 for int64_t: implemented only for precision=DOUBLE and c=r=0.\n";
    // Move the zero rows to the bottom.
    for (long i = 0; i < rank; i++) {
        NTL::swap(gen[i], gen[num - rank + i]);
    }
    // gen.SetDims(rank, gen.NumCols());  // Here, gen is resized.
    return rank;
}

// The ZZ implementation.  We can have `dim < maxDim` only if `precision == DOUBLE`.
template<>
long BasisConstruction<NTL::ZZ>::LLLConstruction0(NTL::matrix<NTL::ZZ> &gen,
        double delta, long r, long c, double *sqlen, PrecisionType precision) {
    long rank = 0;
    if (precision == DOUBLE)
        // This function returns the nonzero vectors in the top rows of `gen`.
        // It never resizes the matrix `gen`.                                    ***
        return rank = NTL::LLL_FPZZflex(gen, delta, r, c, sqlen);

    // If precision is not DOUBLE, we are not allowed to specify r or c
    // with positive values that differ from the dimensions of gen.
    if (((r > 0) & (r != gen.NumRows())) | (((c > 0) & (c != gen.NumCols()))))
        std::cerr
                << "LLLConstruction0: r and c params not allowed when precision != DOUBLE.\n";

    switch (precision) {
    case QUADRUPLE:
        rank = NTL::LLL_QP(gen, delta);
        break;
    case XDOUBLE:
        rank = NTL::LLL_XD(gen, delta);
        break;
    case RR:
        rank = NTL::LLL_RR(gen, delta);
        break;
    default:
        std::cerr << "LLLConstruction0: unknown precision type.\n";
    }
    // NTL puts the zero rows at the top of the matrix `gen`.
    // Here we move the nonzero vectors to the top (first `rank` rows).
    long num = gen.NumRows();
    for (long i = 0; i < rank; i++) {
        NTL::swap(gen[i], gen[num - rank + i]);
    }
    // gen.SetDims(rank, gen.NumCols());
    return rank;
}

//============================================================================

template<typename Int>
void BasisConstruction<Int>::LLLBasisConstruction(IntMat &gen, const Int &m,
        double delta, long r, long c, double *sqlen, PrecisionType precision) {
    int64_t rank = LLLConstruction0(gen, delta, r, c, sqlen, precision);
    if (!c)
        c = gen.NumCols();
    if (rank == c)
        return;  // We are done!

    // We now add the m e_i row vectors, and we redo the LLL with that.
    // No change in the dimensions of gen.
    int64_t i, j;
    for (i = rank; i < rank + c; i++) {
        for (j = 0; j < c; j++) {
            if (i == j)
                gen[i][j] = m;
            else
                gen[i][j] = 0;
        }
    }
    rank = LLLConstruction0(gen, delta, rank + c, c, sqlen);
    std::cout << "LLLBasisConstruction: we had to add some rows!\n";
}

//===================================================================

template<typename Int>
void BasisConstruction<Int>::upperTriangularBasis(IntMat &gen, IntMat &basis,
        const Int &m, long dim1, long dim2) {
    IntVec coeff_gcd, coeff_xi, xi;  // Here we create new vectors!
    Int gcd, gcd_tower, C, D;
    long i, j, k, l;
    // long dim1 = r;
    // long dim2 = c;
    // In case dim1 or dim2 is zero:
    if (dim1 == 0)
        dim1 = gen.NumRows();
    if (dim2 == 0)
        dim2 = gen.NumCols();

    // Allocate space for the vectors:
    coeff_gcd.SetLength(dim1);
    coeff_xi.SetLength(dim1);
    xi.SetLength(dim2);
    if (basis.NumRows() != dim2 || basis.NumCols() != dim2)
        basis.SetDims(dim2, dim2);

    for (i = 0; i < dim2; i++) {
        // Reset these vectors to 0, as they may contain nonzero values from the previous i.
        // xi.clear();   // This call causes a segmentation fault in the int64_t case!
        // coeff_gcd.clear();
        for (j = 0; j < dim1; j++) {
            coeff_gcd[j] = 0;
        }
        for (j = 0; j < dim2; j++) {
            xi[j] = 0;
        }
        // Search for the first non-zero element in the row.
        for (k = 0; (k < dim1 && gen[k][i] == 0); k++) {
        }
        //			if (gen[k][i] != 0)	break;
        // Reduce the other generators as they are used often in what follows.
        for (j = k; j < dim1; j++) {
            NTL::rem(gen[j][i], gen[j][i], m);
        }
        // The `else` case adds m e_i to the basis matrix.
        if (k < dim1) {
            gcd = m;    // Will be GCD(m, gen[k][i]);
            coeff_gcd[k] = 1;
            gcd_tower = gcd;

            // Find the other coefficients by applying the Euclidean algorithm multiple times
            for (j = k; j < dim1; j++) {
                if (gen[j][i] == 0)
                    coeff_gcd[j] = 0;
                else {
                    NTL::XGCD(gcd, C, D, gcd_tower, gen[j][i]);
                    coeff_gcd[j] = D;
                    for (l = 0; l < j; l++) {
                        NTL::mul(coeff_gcd[l], coeff_gcd[l], C);
                    }
                    gcd_tower = gcd;
                }
            }
            // If gcd = m, then this basis (row) vector will be `m e_i`.
            if (gcd == m) {
                for (j = 0; j < dim2; j++) {
                    if (j != i)
                        basis[i][j] = 0;
                    else
                        basis[i][j] = m;
                }
            } else {
                // Reduce the coefficients found during the Euclidean algorithm.
                for (j = 0; j < dim1; j++) {
                    NTL::rem(coeff_gcd[j], coeff_gcd[j], m);
                }
                // We have now found all the coefficients and can compute the vector x_i.
                for (k = 0; k < dim1; k++) {
                    if (coeff_gcd[k] != 0) {
                        for (j = i; j < dim2; j++) {
                            NTL::MulAddTo(xi[j], gen[k][j], coeff_gcd[k]);
                        }
                    }
                }
                // Next we calculate the new vectors v_i.
                // We first calculate the coefficients with which x_i needs to be multiplied.
                for (j = 0; j < dim1; j++) {
                    NTL::div(coeff_xi[j], gen[j][i], gcd);
                    NTL::rem(coeff_xi[j], coeff_xi[j], m);
                }
                for (j = 0; j < dim2; j++)
                    NTL::rem(xi[j], xi[j], m);
                // Update the v_i
                for (k = 0; k < dim1; k++) {
                    if (coeff_xi[k] != 0) {
                        for (j = i; j < dim2; j++) {
                            NTL::MulSubFrom(gen[k][j], coeff_xi[k], xi[j]);
                        }
                    }
                }
//				// Set the `i`th base vector.
                for (j = 0; j < dim2; j++)
                    basis[i][j] = xi[j];
            }
        } else {
            for (j = 0; j < dim2; j++) {
                if (j != i)
                    basis[i][j] = 0;
                else
                    basis[i][j] = m;
            }
        }
    }
}

//==============================================================================

template<typename Int>
void BasisConstruction<Int>::lowerTriangularBasis(IntMat &gen, IntMat &basis,
        const Int &m, long dim1, long dim2) {
    IntVec coeff_gcd, coeff_xi, xi; // Several vectors are created locally here.
    Int gcd, gcd_tower, C, D;
    long i, j, k, l;
    // long dim1 = r;
    // long dim2 = c;
    // In case r or c is zero:
    if (!dim1)
        dim1 = gen.NumRows();
    if (!dim2)
        dim2 = gen.NumCols();

    // Allocate space for the vectors:
    coeff_gcd.SetLength(dim1);
    coeff_xi.SetLength(dim1);
    xi.SetLength(dim2);
    if (basis.NumRows() != dim2 || basis.NumCols() != dim2)
        basis.SetDims(dim2, dim2);

    for (i = dim2 - 1; i > -1; i--) {
        // Reset these vectors to 0, as they may contain nonzero values from the previous i.
        // xi.clear();   // This call causes a segmentation fault in the int64_t case!
        // coeff_gcd.clear();
        for (j = dim1 - 1; j > -1; j--)
            coeff_gcd[j] = 0;
        for (j = dim2 - 1; j > -1; j--)
            xi[j] = 0;
        // Search for the first non-zero element in the row.
        for (k = dim1 - 1; (k > -1 && gen[dim1 - 1 - k][i] == 0); k--) {
        }
        //			if (gen[k][i] != 0)	break;
        // Reduce the other generators as they are used often in what follows.
        for (j = dim1 - 1; j > dim1 - k - 1; j--) {
            NTL::rem(gen[j][i], gen[j][i], m);
        }
        // The `else` case adds m e_i to the basis matrix.
        if (k > -1) {
            gcd = m;    // Will be GCD(m, gen[k][i]);
            coeff_gcd[k] = 1;
            gcd_tower = gcd;

            // Find the other coefficients by applying the Euclidean algorithm multiple times
            for (j = dim1 - 1; j > dim1 - k - 1; j--) {
                if (gen[j][i] == 0)
                    coeff_gcd[j] = 0;
                else {
                    NTL::XGCD(gcd, C, D, gcd_tower, gen[j][i]);
                    coeff_gcd[j] = D;
                    for (l = dim1 - j - 1 - 1; l > -1; l--) {
                        NTL::mul(coeff_gcd[dim1 - 1 - l],
                                coeff_gcd[dim1 - 1 - l], C);
                    }
                    gcd_tower = gcd;
                }
            }
            // If gcd = m, then this basis (row) vector will be `m e_i`.
            if (gcd == m) {
                for (j = dim2 - 1; j > -1; j--) {
                    if (j != i)
                        basis[i][j] = 0;
                    else
                        basis[i][j] = m;
                }
            } else {
                // Reduce the coefficients found during the Euclidean algorithm.
                for (j = 0; j < dim1; j++) {
                    NTL::rem(coeff_gcd[dim1 - 1 - j], coeff_gcd[dim1 - 1 - j],
                            m);
                }
                // We have now found all the coefficients and can compute the vector x_i.
                for (l = dim1 - 1; l > -1; l--) {
                    if (coeff_gcd[l] != 0) {
                        for (j = dim2 - 1; j > dim1 - 1 - i - 1; j--) {
                            NTL::MulAddTo(xi[j], gen[l][j], coeff_gcd[l]);
                        }
                    }
                }
                // Next we calculate the new vectors v_i.
                // We first calculate the coefficients with which x_i needs to be multiplied.
                for (j = dim1 - 1; j > -1; j--) {
                    NTL::div(coeff_xi[j], gen[j][i], gcd);
                    NTL::rem(coeff_xi[j], coeff_xi[j], m);
                }
                for (j = dim2 - 1; j > -1; j--)
                    NTL::rem(xi[j], xi[j], m);
                // Update the v_i
                for (l = dim1 - 1; l > -1; l--) {
                    if (coeff_xi[l] != 0) {
                        for (j = dim2 - 1; j > dim1 - 1 - i - 1; j--) {
                            NTL::MulSubFrom(gen[l][j], coeff_xi[l], xi[j]);
                        }
                    }
                }
                // Set the `i`th base vector.
                for (j = 0; j < dim2; j++)
                    basis[i][j] = xi[j];
            }
        } else {
            for (j = dim2 - 1; j > -1; j--) {
                if (j != i)
                    basis[i][j] = 0;
                else
                    basis[i][j] = m;
            }
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
void BasisConstruction<Int>::mDualUpperTriangular(const IntMat &A, IntMat &B,
        const Int &m, long dim) {
    if (dim == 0) {
        dim = A.NumRows(); // Assumes that the dimension is the number of rows of `IntMat`.
        // B.SetDims(dim, dim);  // B is resized to the dim of the basis.   *******
    }
    if (dim > A.NumCols()) {
        MyExit(1,
                ":mDualUpperTriangular: not enough columns, basis matrix must be square");
        return;
    }
    assert(dim <= B.NumRows());
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

//======================================================

// This is the old version from Couture and L'Ecuyer (1996).
template<typename Int>
void BasisConstruction<Int>::mDualUpperTriangular96(IntMat &basis,
        IntMat &basisDual, const Int &m, long dim) {
    if (m < 1) {
        MyExit(1, "mDualUpperTriangular96: m must be a positive integer");
        return;
    }
    if (!dim) {
        dim = basis.NumRows();
        basisDual.SetDims(dim, dim);   // Here B is resized.
    }
    if (dim != basis.NumCols()) {
        MyExit(1, ":mDualUpperTriangular96: basis matrix must be square");
        return;
    }
    // We must have a triangular basis matrix in the first place.
    // if (!CheckTriangular(basis, basis.NumRows(), Int(0)))
    //   MyExit (1, "mDualUpperTriangular96:  Basis not triangular");
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
            if (basisDual[i][j] != 0)
                basisDual[i][j] = -basisDual[i][j];
            if (!NTL::IsZero(basisDual[i][j] % basis[j][j])) {
                gcd = NTL::GCD(basisDual[i][j], basis[j][j]);
                mm *= basis[j][j] / gcd;
                basisDual *= basis[j][j] / gcd;
            }
            DivideRound(basisDual[i][j], basis[j][j], basisDual[i][j]);
        }
    }
}

// A specialized implementation for ZZ.   WHY THIS?   What is different from Int version?   ****
template<>
void BasisConstruction<NTL::ZZ>::mDualUpperTriangular96ZZ(
        NTL::matrix<NTL::ZZ> &basis, NTL::matrix<NTL::ZZ> &basisDual,
        const NTL::ZZ &m, long dim) {
    if (m < 1) {
        MyExit(1, "mDualUpperTriangular96: m must be a positive integer");
        return;
    }
    if (!dim) {
        dim = basis.NumRows();
        basisDual.SetDims(dim, dim);
    }
    if (dim != basis.NumCols()) {
        MyExit(1, ":mDualUpperTriangular96: basis matrix must be square");
        return;
    }
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
void BasisConstruction<Int>::mDualBasis(const IntMat &basis, IntMat &basisDual,
        const Int &m) {
    std::cerr << "mDualBasis is implemented only for NTL::ZZ integers.\n";
    exit(1);
}

// The specialization for the case where `Int = ZZ`.
template<>
void BasisConstruction<NTL::ZZ>::mDualBasis(const NTL::matrix<NTL::ZZ> &basis,
        NTL::matrix<NTL::ZZ> &basisDual, const NTL::ZZ &m) {
    NTL::ZZ det, fac;
    int64_t dim = basis.NumRows();
    if (dim != basis.NumCols()) {
        std::cerr << "mDualBasis: the given basis matrix must be square.\n";
        exit(1);
    }
    inv(det, basisDual, basis);
    NTL::matrix < NTL::ZZ > C = basisDual;
    div(fac, det, m);
    for (int64_t i = 0; i < dim; i++) {
        for (int64_t j = 0; j < dim; j++) {
            div(basisDual[j][i], C[i][j], fac);
        }
    }
}

/*
 template<>
 void BasisConstruction<NTL::ZZ>::mDualBasis(const NTL::matrix<NTL::ZZ> &basis,
 NTL::matrix<NTL::ZZ> &basisDual, const NTL::ZZ &m, long dim) {
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
 NTL::matrix<NTL::ZZ> C = copyBasisDual;
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
void BasisConstruction<Int>::projectMatrix(const IntMat &in, IntMat &out,
        const Coordinates &proj, long r) {
    if (in == out) MyExit(1, "in and out must be different IntMat objects.");
    if (!r) r = in.NumRows();   // In case r=0.
    // We assume without testing that out is large enough for proj.size().
    for (auto it = proj.begin(), long j = 0; it != proj.end(); it++, j++) {
        for (long i = 0; i < r; i++)  out[i][j] = in[i][*it - 1];
    }
}

//=================================================================================

/*

 template<>
 void BasisConstruction<Int>::projectMatrix_minRows (const IntMat &in,
 IntMat &out, const Coordinates &proj) {
 if (in == out) MyExit(1, "in and out must be different IntMat objects.");
 int64_t projSize = (int64_t) proj.size();
 if (out.NumRows() != projSize || out.NumCols() != projSize)
 out.SetDims(projSize, projSize);
 auto it = proj.cbegin();
 for (int64_t i = 0; i < projSize; i++) {
 out[0][i] = in[0][*it-1];
 for (int j = 1; j < projSize; j++) {
 out[j][i] = 0;
 if (i==j)
 out[j][i] = in[j][j];
 }
 it++;
 }
 // }
 // else
 //	MyExit(1, "A projection coordinate exceeds the dimension of the current basis.");

 };

 //=================================================================================

 template<>
 void BasisConstruction<Int>::projectMatrixDual (const IntMat &in,
 IntMat &out, const Coordinates &proj) {
 if (in == out) MyExit(1, "in and out must be different IntMat objects.");
 int64_t projSize = (int64_t) proj.size();
 if (out.NumRows() != projSize || out.NumCols() != projSize)
 out.SetDims(projSize, projSize);
 auto it = proj.cbegin();
 for (int64_t i = 0; i < projSize; i++) {
 out[i][0] = -in[0][*it-1];
 for (int j = 1; j < projSize; j++) {
 out[j][i] = 0;
 if (i==j)
 out[j][i] = 1;
 }
 it++;
 }
 out[0][0] = in[1][1];
 //}
 //else
 //	MyExit(1, "A projection coordinate exceeds the dimension of the current basis.");

 };

 */

//===================================================
/*
 template<typename Int>
 void BasisConstruction<Int>::projectMatrix0(const IntMat &in, IntMat &out,
 const Coordinates &proj) {
 if (in == out)
 MyExit(1, "in and out must be different IntMat objects.");
 int inDim = in.NumCols();
 uint64_t lat_dim = in.NumCols();
 std::size_t projSize = proj.size();
 out.SetDims(in.NumRows(), projSize);
 auto it = proj.cbegin();
 for (std::size_t i = 0; i < projSize; i++) {
 if (*it <= lat_dim) {
 for (int j = 0; j < inDim; j++) {
 out[j][i] = in[j][*it - 1];
 }
 } else
 MyExit(1,
 "A projection coordinate exceeds the dimension of the current basis.");
 it++;
 }
 }
 */

//===================================================
template<typename Int>
void BasisConstruction<Int>::projectionConstructionLLL(const IntMat &inBasis,
        IntMat &projBasis, const Coordinates &proj, const Int &m,
        const double delta, long r, double *sqlen) {
    projectMatrix(inBasis, projBasis, proj, r);
    LLLBasisConstruction(projBasis, m, delta, r, proj.size(), sqlen);
}

//===================================================

template<typename Int>
void BasisConstruction<Int>::projectionConstructionUpperTri(
        const IntMat &inBasis, IntMat &projBasis, IntMat &genTemp,
        const Coordinates &proj, const Int &m, long r) {
    //  if (!r) r = inBasis.NumRows();
    projectMatrix(inBasis, genTemp, proj, r);
    upperTriangularBasis(genTemp, projBasis, m, r, proj.size());
}

template<typename Int>
void BasisConstruction<Int>::projectionConstructionUpperTri(
        const IntMat &inBasis, IntMat &projBasis, const Coordinates &proj,
        const Int &m, long r) {
    long dim = proj.size();      // Dimension of projection
    IntMat genTemp;
    genTemp.SetDims(dim, dim); // Here an internal object is created and resized!
    projectMatrix(inBasis, genTemp, proj, r);
    upperTriangularBasis(genTemp, projBasis, m, r, dim);
}

//===================================================

template<typename Int>
void BasisConstruction<Int>::projectionConstruction(const IntMat &inBasis,
        IntMat &projBasis, const Coordinates &proj, const Int &m,
        const ProjConstructType projType, const double delta) {
    if (projType == LLLPROJ)
        projectionConstructionLLL(inBasis, projBasis, proj, m, delta);
    if (projType == UPPERTRIPROJ)
        projectionConstructionUpperTri(inBasis, projBasis, proj, m);
}

template class BasisConstruction<std::int64_t> ;
template class BasisConstruction<NTL::ZZ> ;

} // end namespace LatticeTester

#endif
