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
#include "latticetester/LLL_FPZZtest.h"

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
 *
 * We also offer an alternative methods that construct a triangular basis from a set of
 * generating vectors. They always add the rescaled unit vectors implicitly to the set.
 * The method `lowerTriangularBasis` constructs a lower-triangular basis, while
 * `upperTriangularBasis` constructs an upper-triangular basis.
 *
 * To compute the  \f$m\f$-dual of a given basis, we have a general method implemented in
 * `mDualBasis`, and a faster method in `mDualUpperTriangular` that works only when
 * the basis is upper-triangular.
 *
 * We also have functions to compute the basis of a projection of a given lattice over
 * a specified set of coordinates.  The method `projectionConstructionLLL` does this
 * by using LLL to construct the basis of the projection, while `projectionConstructionUpperTri`
 * constructs an upper-triangular basis for the projection.
 *
 * All functions in this class are static, no there is no reason to create any
 * `BasisConstruction` object. We also avoid to create new objects (such as vectors and
 * matrices) inside these functions.  These functions can be called thousands or millions
 * of times in a program, and we want the user to be able to re-use the same vectors and
 * matrices over and over again instead of creating new ones.
 *
 * The programs `BasisManipulationVerbose` and `BasisManipulation` in the examples
 * illustrate how to use these functions and make some speed comparisons.
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
	 * using the NTL implementation specified by `prec`. It returns this basis in `gen`.
	 * See the class `EnumTypes` and the documentation of LLL in NTL for the meaning and
	 * choices for `prec`.
	 * This function is implemented only for the `NTL::ZZ` type, because it uses NTL.
	 * It *does not* assume that all vectors \f$m e_i\f$ belong to the lattice, so
	 * it may return a basis matrix that has fewer rows than columns!
	 * To make sure that these vectors belong to the lattice, we can add them
	 * explicitly beforehand to the set of generating vectors, or call the next method.
	 */
	static void LLLConstruction0(IntMat &gen, double delta = 0.999999,
			PrecisionType precision = DOUBLE);

   /**
    * Similar to `LLLConstruction`, except that in case the set of generating
    * vectors do not generate a full-dimensional lattice, it adds the vectors
    * \f$m e_i\f$ to the generating set, so it always returns a square matrix.
    */
	static void LLLBasisConstruction(IntMat &gen, const Int &m, double delta = 0.999999,
			PrecisionType precision = DOUBLE);

	/**
	 * Takes a set of generating vectors in the matrix `gen` and iteratively
	 * transforms it into a lower triangular lattice basis into the matrix `basis`.
	 * `gen` and `basis` must have the same number of rows and the same number of columns.
	 * All the entries of `gen` given as input are assumed to be reduced modulo `m`.
	 * All the computations are done modulo the scaling factor `m`.
	 * After the execution, `gen` will contain irrelevant information (garbage)
	 * and `basis` will contain an upper triangular basis.
	 * Perhaps with zero rows at the end, in general, unless we assume implicitly that
	 * all vectors of the form m e_i are in the generating set.
	 * The algorithm is explained in the lattice tester guide.
	 * Important: `gen` and `basis` must be different `IntMat` objects.
	 */
	static void lowerTriangularBasis(IntMat &gen, IntMat &basis, const Int &m);

	/**
	 * Similar to `lowerTriangularBasis`, except that the returned basis is upper triangular.
	 */
	static void upperTriangularBasis(IntMat &gen, IntMat &basis, const Int &m);

	/**
	 * Takes an upper triangular basis matrix `basis` and computes the m-dual basis `basisDual`.
	 * The method assumes that each coefficient on the diagonal of `basis` is nonzero and divides `m`.
	 * That is, the basis matrix must be square and invertible.
	 * The algorithm is described in the Lattice Tester guide \cite iLEC22l.
	 * Since the basis is upper triangular, its m-dual will be lower triangular.
	 */
	static void mDualUpperTriangular(const IntMat &basis, IntMat &basisDual,
			const Int &m);

	/**
	 * This function does essentially the same thing as `mDualUpperTriangular`, but it is
	 * slightly different and slower. It uses the method described in \cite rCOU96a.
	 */
	static void mDualUpperTriangular96(IntMat &basis, IntMat &basisDual, Int &m);

	/**
	 * This function assumes that `basis` contains a basis of the primal lattice
	 * scaled by the factor `m`, not necessarily triangular, and it returns in `basisDual`
	 * the m-dual basis.
	 */
	static void mDualBasis(const IntMat &basis, IntMat &basisDual, Int &m);

	/**
	 * Overwrites the matrix 'out' by a matrix formed by the columns of matrix `in` that are
	 * specified by `proj`, and resizes the matrix `out` accordingly.
	 * The matrix `out` will then contain a set of generating vectors for the projection `proj`.
	 * The matrices `in` and `out` must be different IntMat objects, otherwise the program
	 * halts with an error message.
	 */
	static void projectMatrix(const IntMat &in, IntMat &out, const Coordinates &proj);
	
	/**
	 * Constructs a basis for the projection `proj` of the lattice with basis `inBasis`,
	 * using `LLLConstruction`, and puts it in `projBasis`. This basis is not triangular in general.
	 * This will overwrite the matrix `projBasis` and will change its dimensions if needed.
	 * It does not compute the dual.
	 */
	static void projectionConstructionLLL(IntMat &inBasis,
			IntMat &projBasis, const Coordinates &proj, const Int &m, const double delta = 0.9);
	
	/**
	 * Same as `projectionConstructionLLL`, but the construction is made using
	 * `upperTriangularBasis`, so the returned basis is upper triangular.
	 * The matrix `genTemp` will be used to store the generating vectors of the
	 * projection before making the triangularization. We pass it as a parameter
	 * to avoid the interval creation of a new matrix each time.
	 * If `genTemp` is not passed as input, then a temporary matrix will be created internally.
	 */
	static void projectionConstructionUpperTri(IntMat &inBasis,
			IntMat &projBasis, const Coordinates &proj, const Int &m, IntMat &genTemp);

	static void projectionConstructionUpperTri(IntMat &inBasis,
			IntMat &projBasis, const Coordinates &proj, const Int &m);	

	static void projectionConstruction(IntMat &inBasis,
			IntMat &projBasis, const Coordinates &proj, const Int &m, const ProjConstructType t_proj = LLLPROJ, const double delta = 0.9);
	
	
};

//============================================================================
// Implementation

template<typename Int>
void BasisConstruction<Int>::LLLConstruction0(IntMat &gen, double delta,
		PrecisionType prec) {
	std::cerr << "LLLConstruction0 works only for int64_t or NTL::ZZ integers.\n";
	std::cerr << "Aborting.\n";
	exit(1);
}

template<>
void BasisConstruction<int64_t>::LLLConstruction0(NTL::matrix<int64_t> &gen,
		double delta, PrecisionType prec) {
	long num = gen.NumRows();
	int64_t rank=num;
	if (prec == DOUBLE)
		rank = LLL64_FP (gen, delta);
	else
		std::cerr << "LLLConstruction0 for int64_t: implemented only for prec=DOUBLE.\n";
	for (long i = 0; i < rank; i++) {
		NTL::swap(gen[i], gen[num - rank + i]);
	}
	gen.SetDims(rank, gen.NumCols());
}

template<>
void BasisConstruction<NTL::ZZ>::LLLConstruction0(NTL::matrix<NTL::ZZ> &gen,
		double delta, PrecisionType prec) {
	long rank;
	switch (prec) {
	case DOUBLE:
		rank = NTL::LLL_FP(gen, delta);
	    // rank = LLL_FPZZtest(gen, delta);
		break;
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
	long num = gen.NumRows();
	for (long i = 0; i < rank; i++) {
		NTL::swap(gen[i], gen[num - rank + i]);
	}
	gen.SetDims(rank, gen.NumCols());
}

//============================================================================

template<typename Int>
void BasisConstruction<Int>::LLLBasisConstruction(IntMat &gen, const Int &m, double delta,
		PrecisionType prec) {
	LLLConstruction0 (gen, delta, prec);
	int64_t rank = gen.NumRows();
	int64_t dim = gen.NumCols();
	if (rank == dim)
		return;
    gen.SetDims(rank + dim, dim);
    // We now add the m e_i vectors, and we redo the LLL.
    int64_t i, j;
    for (i=rank; i < rank+dim; i++) {
        for (j=0; j < dim; j++) {
        	if (i==j) gen[i][j] = m; else gen[i][j] = 0;
        }
    }
	LLLConstruction0 (gen, delta, prec);
	std::cout << "LLLBasisConstruction: we had to add some rows!\n";
}

//===================================================================

template<typename Int>
void BasisConstruction<Int>::upperTriangularBasis
         (IntMat &gen, IntMat &basis, const Int &m) {
	IntVec coeff_gcd, coeff_xi, xi;  // Here we create new vectors!
	Int gcd, gcd_tower, C, D;
	long dim1 = gen.NumRows();
	long dim2 = gen.NumCols();
	long i, j, k, l;
	
	// Define dimensions of vectors
	coeff_gcd.SetLength(dim1);
	coeff_xi.SetLength(dim1);
	xi.SetLength(dim2);
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
		for (k = 0; (k < dim1 && gen[k][i] == 0); k++) {}
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
			if (gcd==m) {
				for (j = 0; j < dim2; j++) {
				  if (j != i)
					  basis[i][j] = 0;
				  else
				  	basis[i][j] = m;
				}
			}
			else {
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
				// Set the `i`th base vector.
				basis[i] = xi;
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
		const Int &m) {
	IntVec coeff_gcd, coeff_xi, xi;  // Vectors are created locally here.
	Int gcd, gcd_tower, C, D;
	long dim1 = gen.NumRows();
	long dim2 = gen.NumCols();
	long i, j, k, l;

	//Define dimensions of vectors
	coeff_gcd.SetLength(dim1);
	coeff_xi.SetLength(dim1);
	xi.SetLength(dim2);
    basis.SetDims(dim2, dim2);
    
	for (i = dim2-1; i > -1; i--) {
		// Reset these vectors to 0, as they may contain nonzero values from the previous i.
		// xi.clear();   // This call causes a segmentation fault in the int64_t case!
		// coeff_gcd.clear();
		for (j = dim1-1; j >-1; j--) 
		    coeff_gcd[j] = 0;
		for (j = dim2-1; j >-1; j--) 
		    xi[j] = 0;
		// Search for the first non-zero element in the row.
		for (k = dim1-1; (k > -1 && gen[dim1-1-k][i] == 0); k--) {}
		//			if (gen[k][i] != 0)	break;
		// Reduce the other generators as they are used often in what follows.
		for (j = dim1-1; j > dim1-k-1; j--) {
		    NTL::rem(gen[j][i], gen[j][i], m);
		}
		// The `else` case adds m e_i to the basis matrix.
		if (k > -1) {
			gcd = m;    // Will be GCD(m, gen[k][i]);
			coeff_gcd[k] = 1;
			gcd_tower = gcd;

			// Find the other coefficients by applying the Euclidean algorithm multiple times
			for (j = dim1-1; j >dim1-k-1; j--) {
				if (gen[j][i] == 0)
					coeff_gcd[j] = 0;
				else {
					NTL::XGCD(gcd, C, D, gcd_tower, gen[j][i]);
					coeff_gcd[j] = D;
					for (l = dim1-j-1-1; l > -1; l--) {
						NTL::mul(coeff_gcd[dim1-1-l], coeff_gcd[dim1-1-l], C);
					}
					gcd_tower = gcd;
				}
			}
			// If gcd = m, then this basis (row) vector will be `m e_i`.
			if (gcd==m) {
				for (j = dim2-1; j > -1; j--) {
				  if (j != i)
					  basis[i][j] = 0;
				  else
					  basis[i][j] = m;
				}
			}
			else {
				// Reduce the coefficients found during the Euclidean algorithm.
				for (j = 0; j < dim1; j++) {
				  NTL::rem(coeff_gcd[dim1-1-j], coeff_gcd[dim1-1-j], m);
				}
				// We have now found all the coefficients and can compute the vector x_i.
				for (l = dim1-1; l > -1; l--) {
					if (coeff_gcd[l] != 0) {
						for (j = dim2-1; j > dim1-1-i-1; j--) {
							NTL::MulAddTo(xi[j], gen[l][j], coeff_gcd[l]);
						}
					}
				}
				// Next we calculate the new vectors v_i.
				// We first calculate the coefficients with which x_i needs to be multiplied.
				for (j = dim1-1; j > -1; j--) {
					NTL::div(coeff_xi[j], gen[j][i], gcd);
					NTL::rem(coeff_xi[j], coeff_xi[j], m);
				}
				for (j = dim2-1; j> -1; j--) 
					NTL::rem(xi[j], xi[j], m);
				// Update the v_i
				for (l = dim1-1; l > -1; l--) {
					if (coeff_xi[l] != 0) {
						for (j = dim2-1; j > dim1-1-i-1; j--) {
							NTL::MulSubFrom(gen[l][j], coeff_xi[l], xi[j]);
						}
					}
				}
				// Set the `i`th base vector.
				basis[i] = xi;
			}
		} else {
			for (j = dim2-1; j > -1; j--) {
				if (j != i)
					basis[i][j] = 0;
				else
					basis[i][j] = m;
			}
		}
	}
}

//======================================================

// This is the old version from Couture and L'Ecuyer (1996).
template<typename Int>
void BasisConstruction<Int>::mDualUpperTriangular96(IntMat &basis,
		IntMat &basisDual, Int &m) {
	if (m < 1) {
		std::cerr << "m has to be a positive integer.\n";
		exit(1);  return;
	}
	// We must have a triangular basis matrix in the first place.
	// if (!CheckTriangular(basis, basis.NumRows(), Int(0)))
	//	  MyExit (1, "mDualUpperTriangular96:  Basis not triangular");
	long dim = basis.NumRows();
	if (dim != basis.NumCols()) {
		std::cout << ":mDualUpperTriangular96: basis matrix must be square.\n";
		return;
	}
	basisDual.SetDims(dim, dim);
	Int mm = m;            // Local copy of m that can be changed.
	for (int64_t i = 0; i < dim; i++) {
		for (int64_t j = i + 1; j < dim; j++)
			NTL::clear(basisDual(i, j));
		if (!NTL::IsZero(basisDual(i, i))) {
			Int gcd = NTL::GCD(mm, basis(i, i));
			mm *= basis(i, i) / gcd;
			basisDual *= basis(i, i) / gcd;
		}

		DivideRound(mm, basis(i, i), basisDual(i, i));
		for (int64_t j = i - 1; j >= 0; j--) {
			NTL::clear(basisDual(i, j));
			for (int64_t k = j + 1; k <= i; k++)
				basisDual(i, j) += basis(j, k) * basisDual(i, k);
			if (basisDual(i, j) != 0)
				basisDual(i, j) = -basisDual(i, j);
			if (!NTL::IsZero(basisDual(i, j) % basis(j, j))) {
				Int gcd = NTL::GCD(basisDual(i, j), basis(j, j));
				mm *= basis(j, j) / gcd;
				basisDual *= basis(j, j) / gcd;
			}
			DivideRound(basisDual(i, j), basis(j, j), basisDual(i, j));
		}
	}
}

template<>
void BasisConstruction<NTL::ZZ>::mDualUpperTriangular96( NTL::matrix<NTL::ZZ>  &basis,
		NTL::matrix<NTL::ZZ>  &basisDual, NTL::ZZ &m) {
	if (m < 1) {
		std::cerr << "m has to be a positive integer.\n";
		exit(1);  return;
	}
	NTL::ZZ gcd, fac;
	NTL::ZZ mm = m;            // Local copy of m that can be changed.
	int64_t dim;
	dim = basis.NumRows();
	basisDual.SetDims(dim, dim);
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
				negate(basisDual[i][j], basisDual[i][j]);
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
		const Int &m) {
	int64_t dim = A.NumRows();
	B.SetDims(dim, dim);
	for (int64_t i = 0; i < dim; i++) {
		for (int64_t j = i + 1; j < dim; j++)
			NTL::clear(B(i, j));
		DivideRound(m, A(i, i), B(i, i));
		for (int64_t j = i - 1; j >= 0; j--) {
			NTL::clear(B(i, j));
			for (int64_t k = j + 1; k <= i; k++)
				B(i, j) += A(j, k) * B(i, k);
			if (B(i, j) < 0)
				B(i, j) = -B(i, j);
			DivideRound(B(i, j), A(j, j), B(i, j));
		}
	}
}

template<>
void BasisConstruction<int64_t>::mDualUpperTriangular(
		const NTL::matrix<int64_t> &A, NTL::matrix<int64_t> &B, const long &m) {
	long dim = A.NumRows();
    long i, j, k;
	for (i = 0; i < dim; i++) {
		for (j = i + 1; j < dim; j++)
			NTL::clear(B[i][j]);
		NTL::div(B[i][i], m, A[i][i]);
		for (j = i - 1; j >= 0; j--) {
			NTL::clear(B[i][j]);
			for (k = j + 1; k <= i; k++)
				NTL::MulSubFrom(B[i][j], A[j][k], B[i][k]);
			NTL::div(B[i][j], B[i][j], A[j][j]);
		}
	}
}


template<>
void BasisConstruction<NTL::ZZ>::mDualUpperTriangular(
		const NTL::matrix<NTL::ZZ> &A, NTL::matrix<NTL::ZZ> &B,	const NTL::ZZ &m) {
	int64_t dim = A.NumRows();
	for (int64_t i = 0; i < dim; i++) {
		for (int64_t j = i + 1; j < dim; j++)
			NTL::clear(B[i][j]);
		div(B[i][i], m, A[i][i]); 
		for (int64_t j = i - 1; j >= 0; j--) {
			NTL::clear(B[i][j]);
			for (int64_t k = j + 1; k <= i; k++)
				MulSubFrom(B[i][j], A[j][k], B[i][k]);
			div(B[i][j], B[i][j], A[j][j]);		
		}
	}
}

//===================================================

template<typename Int>
void BasisConstruction<Int>::mDualBasis(const IntMat &basis, IntMat &basisDual,
		Int &m) {
	std::cerr << "mDualBasis is implemented only for NTL::ZZ integers.\n";
	std::cerr << "Aborting.\n";
	exit(1);
}

// The specialization for the case where `Int = ZZ`.
template<>
void BasisConstruction<NTL::ZZ>::mDualBasis(
		const NTL::matrix<NTL::ZZ> &basis, NTL::matrix<NTL::ZZ> &basisDual, NTL::ZZ &m) {
	NTL::ZZ d, fac;
	int64_t dim = basis.NumRows();
	if (dim != basis.NumCols()) {
		std::cerr << "mDualBasis: the given basis matrix must be square.\n";
		exit(1);
	}
	inv(d, basisDual, basis);
	NTL::matrix<NTL::ZZ> C = basisDual;
	div(fac, d, m);
	for (int64_t i = 0; i < dim; i++) {
		for (int64_t j = 0; j < dim; j++) {
			div(basisDual[j][i], C[i][j], fac);
		}
	}
}

//=================================================================================

// The matrix out will always be resized to a number of columns equal to the
// number of coordinates in `proj`.
template<>
void BasisConstruction<NTL::ZZ>::projectMatrix (const IntMat &in,
		IntMat &out, const Coordinates &proj) {
	if (in == out) MyExit(1, "in and out must be different IntMat objects.");
    int inDim = in.NumCols();
    uint64_t lat_dim = in.NumCols();   
	std::size_t projSize = proj.size();
	out.SetDims(in.NumRows(), projSize);   // here we resize the matrix each time!
	auto it = proj.cbegin();
	for (std::size_t i = 0; i < projSize; i++) {
	    if (*it <= lat_dim) {
			for (int j = 0; j < inDim; j++) {
			   out[j][i] = in[j][*it];
			}
		}
		else
			MyExit(1, "A projection coordinate exceeds the dimension of the current basis.");
		it++;
	}
};

template<>
void BasisConstruction<Int>::projectionConstructionLLL(
		IntMat &inBasis, IntMat &projBasis, const Coordinates &proj, const Int &m, double delta) {
	projectMatrix(inBasis, projBasis, proj);
	LLLBasisConstruction(projBasis, m, delta);
}

template<>
void BasisConstruction<Int>::projectionConstructionUpperTri(
		IntMat &inBasis, IntMat &projBasis, const Coordinates &proj, const Int &m, IntMat &genTemp) {
	projectMatrix(inBasis, genTemp, proj);
	upperTriangularBasis(genTemp, projBasis, m);
}

template<>
void BasisConstruction<Int>::projectionConstructionUpperTri(
		IntMat &inBasis, IntMat &projBasis, const Coordinates &proj, const Int &m) {
	IntMat genTemp;
	projectMatrix(inBasis, genTemp, proj);
	upperTriangularBasis(genTemp, projBasis, m);
}


template<>
void BasisConstruction<Int>::projectionConstruction(IntMat &inBasis,
			IntMat &projBasis, const Coordinates &proj, const Int &m, const ProjConstructType t_proj, const double delta) {
	if (t_proj == LLLPROJ)
		projectionConstructionLLL(inBasis, projBasis, proj, m, delta);
	if (t_proj == UPPERTRIPROJ)
		projectionConstructionUpperTri(inBasis, projBasis, proj, m);
}
	
//void BasisConstruction<Int>::projectionConstruction(IntMat &inBasis, IntMat &projBasis, const Coordinates &proj, const Int &m, double delta)

template class BasisConstruction<std::int64_t>;
template class BasisConstruction<NTL::ZZ>;

} // end namespace LatticeTester

#endif
