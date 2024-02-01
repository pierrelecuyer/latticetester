// This file is part of LatticeTester.
//
// Copyright (C) 2012-2022  The LatticeTester authors, under the occasional supervision
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

#ifndef LATTICETESTER_INTLATTICE_H
#define LATTICETESTER_INTLATTICE_H

#include "latticetester/NTLWrap.h"
#include "latticetester/EnumTypes.h"
#include "latticetester/Util.h"
#include "latticetester/Coordinates.h"
#include "latticetester/BasisConstruction.h"

#include <string>
#include <sstream>
#include <iomanip>
#include <cassert>

namespace LatticeTester {

/**
 * An `IntLattice` object is an integral lattice, with its basis or its `m`-dual basis, or both.
 * There are tools to perform simple manipulations on those lattice bases.
 * The value of `m` must be chosen in a way that all coordinates of the basis and
 * of its `m`-dual are integers, so they can be represented exactly.
 * The basis or its `m`-dual is rescaled by `m`, which is typically the smallest integer with this property.
 *
 * The dimension `dim` of the lattice is the number of independent vectors that form a basis.
 * Usually, these vectors also have `dim` coordinates, but in general they may have more.
 * The basis and/or the `m`-dual basis are stored in `IntMat` arrays (from NTL) of sizes `maxDim x maxDim`,
 * where `maxDim is usually fixed to a value as large as the largest `dim` that we want to handle.
 * These arrays are then allocated only once and never have to be resized, which improves speed.
 * A boolean variable `withPrimal` indicates if we maintain the primal basis and a variable
 * `withDual` indicates if we maintain the dual basis. At least one of them (or both) should be true.
 *
 * A norm is also chosen in `NormType` to measure the vector lengths; by default it is the
 * Euclidean norm.
 * Methods and attributes are offered to compute and store the norms of the basis and/or
 * the m-dual basis vectors, to permute these vectors, sort them by length, etc.
 * An `IntLattice` object contains protected variables to store all these quantities.
 * For better efficiency, we should avoid creating many `IntLattice` objects, for example when
 * making searches for good lattices. We should try to reuse the same one as much as we can.
 *
 * Often, we want to assess the quality of projections of the lattice over subsets of coordinates.
 * The method `buildProjection` computes a basis for the lattice defined as the
 * projection of the full lattice on a subset \f$I\f$  of coordinates indices
 * defined by a `Coordinates` object.  When computing figures of merit, one may want to
 * recompute a basis for a large number of different projections, which can be specified
 * by a `CoordinateSets` object.
 *
 * NOTE: There are no methods to copy or overwrite only the dual lattice!!!  Maybe not needed?
 *
 * The class `IntLatticeExt` extends this class and contains virtual methods that must
 * be defined in its subclasses.
 */

template<typename Int, typename Real>
class IntLattice {

private:
	typedef NTL::vector<Int> IntVec;
	typedef NTL::matrix<Int> IntMat;
	typedef NTL::vector<Real> RealVec;
	typedef NTL::matrix<Real> RealMat;

public:

	/**
	 * Constructs a lattice whose basis is the identity, in `maxDim` dimensions,
	 * with the specified norm type, and the scaling factor `m`.
	 * The parameter `withPrimal` and `withDual` indicate if the primal basis and the m-dual
	 * basis are maintained and updated (or not).  At least one of them must be `true`.
	 */
	IntLattice(const Int m, const int64_t maxDim, bool withPrimal = true, bool withDual = false,
			   NormType norm = L2NORM);

	/**
	 * Similar to the previous constructor, except that the primal basis is given in
	 * `basis`, which must be an `IntMat` object of size `maxDim` by `maxDim`.
     * Here, `withPrimal` is set to `true`, and the m-dual basis will be computed
     * only if `withDual == true`.
	 */
	IntLattice(const IntMat basis, const Int m, const int64_t maxDim,
			   bool withDual = false, NormType norm = L2NORM);

	/**
	 * Constructs a lattice with the given basis and given m-dual basis for the given `m`,
	 * in `maxDim` dimensions, and with the specified norm type.
	 * In this case, by default, both the primal and m-dual basis will be maintained.
	 * The two `IntMat` objects must be of size `maxDim` by `maxDim`.
	 */
	IntLattice(const IntMat primalbasis, const IntMat dualbasis,
			   const Int m, const int64_t maxDim, NormType norm = L2NORM);

	/**
	 * Copy constructor. Makes a deep copy of `lat` into `*this` new object.
	 */
	IntLattice(const IntLattice<Int, Real> &lat);

	/**
	 * Destructor.
	 */
	~IntLattice();

	/**
	 * Makes a deep copy of the lattice `lat` into this (existing) object.
	 * New matrix and vector objects are constructed to store the bases and norms.
	 */
	void copyLattice(const IntLattice<Int, Real> &lat);

	/*
	 * *** Previously named `copyLattice`.
	 * Overwrites the first `dim` elements of the basis of the lattice `lat` over the elements
	 * of the basis of the current object, in the upper left corner of the basis matrix.
	 * The vector norms and the dual basis (if available) are also overwritten.
	 * The difference with `copyLattice` is that here, no new matrix or vector is constructed;
	 * the previous ones are re-used. Requirement: `dim <= maxDim`.
	 */
	void overwriteLattice(const IntLattice<Int, Real> &lat, long d);
	 
	/**
	 * **** FROM PIERRE:  Not clear to me if we want this general method here.   ****
	 * **** It is already in BasisConstruction, no?
	 * ****  Can we just get rid of this one?   *****
	 *
	 * Builds an upper triangular basis for the projection `proj` for this lattice
	 * and replaces the current `lattice` object by this projection.
	 * If the latter maintains a dual basis, then the (triangular) m-dual basis is also updated.
	 * Note that the same `lattice` objects can be used when calling this method several
	 * times to examine different projections.
	 */
	void buildProjection(IntLattice<Int, Real> *lattice, const Coordinates &proj);

	/**
	 * Initializes a vector containing the norms of the basis vectors to -1
	 * for all components.  It means the norms are no longer up to date.
	 * Bad name !!!!
	 */
	// void initVecNorm();

	/**
	 * Returns the `IntMat` object that contains the basis of this lattice.
	 */
	IntMat& getBasis() {
		return m_basis;
	}
	
	/**
	 * Changes the basis to 'basis', the modulus to 'm', and the current dimension to 'dim'.
	 * Does not change `maxDim`.
	 */
	void setBasis(const IntMat basis, const Int m, const int64_t dim,
			bool withDual = false, NormType norm = L2NORM) {
	   this->m_modulo=m;
	   this->m_withDual=withDual;
	   this->m_norm=norm;
  	   this->m_dim=dim;
       this->m_basis=basis;
	   // this->m_vecNorm.resize(dim);
	   setNegativeNorm();     
	}
 
    /**
     * Changes only the 'basis' and the current dimension to 'dim'.
     * Does not change `m` and `maxDim`.
     */
 	void setBasis(const IntMat basis, const int64_t dim) {
	  this->m_dim=dim;
	  this->m_basis=basis;
	  // this->m_vecNorm.resize(dim);
	  setNegativeNorm();     
 	}

	/**
	 * Returns the m-dual basis represented in a matrix.
	 */
	IntMat& getDualBasis() {
		return m_dualbasis;
	}
	
	/* 
	 * This function calculates a projection 'projBasis' of the basis.
	 * The coordinates of the projection are given by 'coord'.
	*/
	//void getProjBasis (const Coordinates & coord, IntMat projBasis);
	
	/* 
	 * This function calculates a projection 'projBasisDual' of the dual basis.
	 * The coordinates of the projection are given by 'coord'.
	*/
	//void getProjBasisDual (const Coordinates & coord, IntMat projBasisDual);
	
	/* 
	 * This function calculates a projection 'projBasis' of the basis and of the dual basis.
	 * The coordinates of the projection are given by 'coord'.
	*/
	//void getProjBasisPrimalDual (const Coordinates & coord, IntMat projBasis, IntMat projDualBasis);

	/**
	 * Returns the current dimension of the lattice, which is the dimension of the basis vectors,
	 * and also usually the number of independent vectors in the basis.
	 */
	int64_t getDim() const {
		return m_dim;
	}

	/**
	 * Returns the `NormType` used by this lattice.
	 */
	NormType getNormType() const {
		return m_norm;
	}

	/**
	 * Returns the norm (squared in case of the L^2 norm) of the i-th vector of the basis,
	 * with the index i starting at 0.
	 */
	Real getVecNorm(const int64_t &i) {
		return m_vecNorm[i];
	}

	/**
	 * Returns the norm (squared in case of the L^2 norm) of each basis vector, in a vector.
	 */
	RealVec getVecNorm() const {
		return m_vecNorm;
	}

	/**
	 * Returns the norm (squared in case of the L^2 norm) of the i-th vector of the m-dual basis.
	 * The index i starts at 0.
	 */
	Real getDualVecNorm(const int64_t &i) {
		return m_dualvecNorm[i];
	}

	/**
	 * Returns the norm (squared in case of the L^2 norm) of each vector of the m-dual basis, in a vector.
	 */
	RealVec getDualVecNorm() const {
		return m_dualvecNorm;
	}

	/**
	 * Returns the scaling factor `m`.
	 */
	Int getModulo() const {
		return m_modulo;
	}

	/**
	 * Sets the dimension of the basis to `dim`. This does not change `maxDim` nor any of the
	 * basis vectors, but only the dimension variable.
	 */
	void setDim(const int64_t &dim) {
		if (dim > 0)
			m_dim = dim;
	}

	/**
	 * Sets the `NormType` used by this lattice to `norm`.
	 */
	void setNormType(const NormType &norm) {
		m_norm = norm;
	}

	/**
	 * Sets the norm of the `i`-th component of the basis to `value`, which is assumed to
	 * be the correct value.  To recompute the norm, use `updateVecNorm(const int64_t&)` instead.
	 */
	void setVecNorm(const Real &value, const int64_t &i) {
		m_vecNorm[i] = value;
	}

	/**
	 * Sets the norm of the `i`-th component of the m-dual basis to `value`,  which is assumed to
	 * be the correct value.  To recompute the norm, use `updateDualVecNorm(const int64_t&)` instead.
	 */
	void setDualVecNorm(const Real &value, const int64_t &i) {
		m_dualvecNorm[i] = value;
	}

	/**
	 * Returns `true` iff we maintain a primal basis.
	 */
	bool withPrimal() {
		return m_withPrimal;
	}

	/**
	 * Sets the `withPrimal` flag to `flag`. This flag indicates whether or
	 * not we maintain a primal basis for this `IntLattice`.  It is the flag
	 * returned by `withPrimal()`.
	 */
	void setPrimalFlag(bool flag) {
		m_withPrimal = flag;
	}

	/**
	 * Returns `true` iff we maintain an m-dual basis.
	 */
	bool withDual() {
		return m_withDual;
	}

	/**
	 * Sets the `withDual` flag to `flag`. This flag indicates whether or
	 * not we maintain an `m`-dual basis for this `IntLattice`. It is the flag
	 * returned by `withDual()`.
	 */
	void setDualFlag(bool flag) {
		m_withDual = flag;
	}

	/**
	 * Sets to -1 the first `dim` values in the array containing the norms of the basis vectors,
	 * where `dim` is the current dimension of the lattice.
	 * This means that these norms are no longer up to date.
	 */
	void setNegativeNorm();

	/**
	 * Sets the value of the `i`-th component in the array containing the
	 * norms of the basis vectors to -1.
	 */
	void setNegativeNorm(const int64_t &i) {
		m_vecNorm[i] = -1;
	}

	/**
	 * Sets all the values in the array containing the norms of the dual basis
	 * vectors to -1, to indicate that these norms are no longer up to date.
	 */
	void setDualNegativeNorm();

	/**
	 * Sets the value of the `i`-th component in the array containing the
	 * norms of the m-dual basis vectors to -1.
	 */
	void setDualNegativeNorm(const int64_t &i) {
		m_dualvecNorm[i] = -1;
	}

	/**
	 * ****   FROM PIERRE: We may ask that this function returns the shortest one!   *****
	 * Updates the array containing the norms of the basis vectors by recomputing them.
	 * Returns the length (squared in case of the L^2 norm) of the shortest vector in the current basis.
	 *
	 * **** We could either return the length as a Real, or perhaps just the index of
	 * **** the vector that is shortest!  Then we can easily look at its length if needed.
	 */
	Real updateVecNorm();

	/**
	 * Updates the array containing the norm of the basis vectors from the `d`-th
	 * component to the last, by recomputing them.
	 * Putting `d=0` recomputes all the norms.
	 */
	Real updateVecNorm(const int64_t &d);

	/**
	 * Updates the array containing the m-dual basis vectors norms by recomputing them.
	 * Assumes that the dual basis is available.
	 * Returns the length (squared in case of the L^2 norm) of the shortest vector in the
	 * current m-dual basis.
	 */
	Real updateDualVecNorm();
	
	/**
	 * Updates the 'd'-th entry of the array containing the m-dual basis vectors norms.
	 * Only the first c components are used for calculating the norm.
	 * */
	void updateSingleDualVecNorm(const int64_t &d, const int64_t & c);

	/**
	 * Updates the array containing the m-dual basis vectors norms from the `d`-th
	 * component to the last by recomputing them.
	 * */
	Real updateDualVecNorm(const int64_t &d);

	/**
	 * Updates the array containing the m-dual basis vectors norms from the `d`-th
	 * component to the last by recomputing them. Only the first c components are
	 * used for calculating the norm
	 * */
	void updateDualVecNorm(const int64_t &d, const int64_t & c);
	

	/**
	 * Updates the `i`-th value of the array containing the square norms of the
	 * basis vectors by recomputing it using the `L2NORM`.
	 */
	void updateScalL2Norm(const int64_t i);

	/**
	 * Updates the `k1`-th to the `k2-1`-th values of the array containing
	 * the square norms of the basis vectors by recomputing them using the `L2NORM`.
	 */
	void updateScalL2Norm(const int64_t k1, const int64_t k2);

	/**
	 * Updates the `i`-th value of the array containing the square norms of the
	 * m-dual basis vectors by recomputing it using the `L2NORM`.
	 */
	void updateDualScalL2Norm(const int64_t i);

	/**
	 * Updates the `k1`-th to the `k2-1`-th values of the array containing
	 * the square norms of the m-dual basis vectors by recomputing them using the `L2NORM`.
	 */
	void updateDualScalL2Norm(const int64_t k1, const int64_t k2);

	/*
	 * ****  PIERRE: The two functions below are doing a lot of work, perhaps more than needed!
	 * ****  They recompute the norms even when the norms are available!
	 * ****  Also, We really want to avoid computing square roots, in all cases!    ****
	 * ****  Finally, in realistic cases, `double` will be too small!!!  We need Real.
	 * See `updateVecNorm` above.
	 *
	 * Returns the length (squared in case of the L^2 norm) of the shortest vector in the current basis.
	 */
	double getShortestLengthBasis();

	/*
	 * Returns the length (squared in case of the L^2 norm) of the shortest vector in the current m-dual basis.
	 */
	double getShortestLengthDualBasis();

	/**
	 * Exchanges vectors `i` and `j` in the basis. This also changes the
	 * m-dual basis vectors and the arrays containing secondary information
	 * about the two bases (like the norms) accordingly.
	 */
	void permute(int64_t i, int64_t j);

	/**
	 * Exchanges vectors `i` and `j` in the basis without changing the m-dual.
	 */
	void permutePrimal(int64_t i, int64_t j);

	/**
	 * Exchanges vectors `i` and `j` in the `m`-dual basis without changing the primal.
	 */
	void permuteDual(int64_t i, int64_t j);

	/**
	 * Exchanges the primal and m-dual bases and vector norms, and the indicator variables
	 * `withPrimal` and `withDual`.
	 */
	 void dualize ();

	/**
	 * Returns `true` iff the m-dual basis contained in the object really is
	 * the m-dual of the current primal basis. Otherwise, or if either the primal or
	 * m-dual basis is not maintained, it returns false.
	 */
	bool checkDuality();

	/**
	 * Sorts the primal basis vectors with indices greater of equal to `d` by
	 * increasing length. The m-dual vectors are permuted accordingly. Assumes
	 * that the lengths (norms) of the corresponding basis vectors are up to date.
	 */
	void sortBasis (int64_t d);
	
	/**
	 * Sorts the primal basis vectors with indices greater of equal to `d` by
	 * increasing length. The m-dual vectors (if maintained) are **not** permuted.
	 */
	void sortPrimalBasis (int64_t d);
	
	/**
	 * Sorts the `m`-dual basis vectors with indices greater of equal to `d` by
	 * increasing length. The primal basis vectors are **not** permuted.
	 */
	void sortDualBasis (int64_t d);

	/**
	 * Returns a string that contains the primal basis vectors and their norms.
	 */
	std::string toStringBasis() const;

	/**
	 * Returns a string with the `m`-dual basis vectors and their norms.
	 */
	std::string toStringDualBasis() const;

	/**
	 * Returns a string that represents this lattice and its parameters.
	 * It contains the dimension, the norm used, the basis and m-dual basis vectors and
	 * the basis and dual basis vector norms.
	 */
	std::string toString() const;

	/**
	 * Writes on standard output the string returned by `toString`.
	 */
	void write() const;

protected:

    /**
	 * The scaling factor `m` used for rescaling the lattice.
	 */
	Int m_modulo;

	/**
	 * The maximal dimension for this lattice.  It should correspond to the size of the `IntMat`
	 * objects that contain the basis and/or its m-dual.
	 */
	int64_t m_maxDim;

	/**
	 * The current dimension of the lattice, which is the number of (independent) vectors
	 * in the basis. It cannot exceed the number of coordinates in those vectors.
	 * It also cannot exceed m_maxDim.
	 */
	int64_t m_dim;

	/**
	 * The rows of the m_dim x m_dim upper left corner of this matrix are the primal basis vectors.
	 */
	IntMat m_basis;

	/**
	 * The rows of the m_dim x m_dim upper left corner this matrix are the m-dual basis vectors.
	 * May not be initialized.  When m_withDual = true, it must be initialized.
	 */
	IntMat m_dualbasis;

	/**
	 * The NormType used to measure the vector lengths for this lattice.
	 * It is used for the basis reduction and compute a shortest vector, for example.
	 */
	NormType m_norm;

	/**
	 * A vector that stores the norm of each basis vector.
	 * In case of the L_2 norm, it contains the square norm instead.
	 * A value of -1 means that the norm is not up to date.
	 */
	RealVec m_vecNorm;

	/**
	 * Similar to vecNorm, but for the m-dual basis.
	 */
	RealVec m_dualvecNorm;

	/**
	 * This variable is `true` iff a primal basis is available.
	 *
	 * NOTE: We thought about adding this variable for the case where we DO NOT
	 * maintain the primal, but only the dual. We decided to leave it out for now,
	 * to avoid adding more parameters to the functions.
	 */
	// bool m_withPrimal;

	/**
	 * This variable is `true` iff an m-dual basis is available.
	 */
	bool m_withDual;

	/**
	 * `true` iff the current basis is triangular (in case we want to check).
	 */
	// bool m_triangularBasis;

	/**
	 * The dimension of the current projection. It should not exceed m_dim.
	 */
	// int64_t m_dimProj;

	/**
	 * The primal basis of the current projection.
	 */
	//IntMat m_basisProj;

	/**
	 * The m-dual basis of the current projection.   NEEDED ?   **********
	 */
	// IntMat m_dualbasisProj;

	/**
	 * Allocates space to the vectors m_basisProj and m_dualbasisProj used internally to store
	 * the bases for the current projection.
	 * This should not be called directly by the user.
	 * It is called by the constructors and the copy method.
	 */
	//void initProj();

};

// class IntLattice

//===========================================================================

template<typename Int, typename Real>
IntLattice<Int, Real>::IntLattice(const Int m, const int64_t maxDim, bool withPrimal,
		 bool withDual, NormType norm) {
//		: m_modulo(m), m_dim(dim), m_withDual(withDual), m_norm(norm) {
    this->m_modulo=m;
    this->m_maxDim=maxDim;
    this->m_dim=maxDim;
    this->m_withPrimal=withPrimal;
    this->m_withDual=withDual;
    this->m_norm=norm;
	this->m_basis.resize(maxDim, maxDim);
	this->m_vecNorm.resize(maxDim);
	setNegativeNorm();
}

//===========================================================================

template<typename Int, typename Real>
IntLattice<Int, Real>::IntLattice (const IntMat basis, const Int m,
		const int64_t maxDim, bool withDual, NormType norm) {
//		: m_basis(basis), m_modulo(m), m_dim(dim), m_withDual(withDual), m_norm(norm) {
    this->m_modulo=m;
    this->m_maxDim=maxDim;
    this->m_dim=maxDim;
    this->m_withPrimal=true;
    this->m_withDual=withDual;
    this->m_norm=norm;
    assert (basis.NumRows() == maxDim);
    this->m_basis=basis;
    this->m_vecNorm.resize(maxDim);
	setNegativeNorm();
}

/*=========================================================================*/

template<typename Int, typename Real>
IntLattice<Int, Real>::IntLattice(const IntMat primalbasis,
		const IntMat dualbasis, const Int m, const int64_t maxDim, NormType norm) :
		IntLattice<Int, Real>(primalbasis, m, maxDim, true, norm) {
    assert (dualbasis.NumRows() == maxDim);
    this->m_withPrimal=true;
    this->m_withDual=true;
    this->m_dualbasis = dualbasis;
	this->m_dualvecNorm.resize(maxDim);
	setDualNegativeNorm();
}

/*=========================================================================*/

template<typename Int, typename Real>
IntLattice<Int, Real>::IntLattice(const IntLattice<Int, Real> &lat) {
	copyLattice(lat);
}

/*=========================================================================*/

template<typename Int, typename Real>
IntLattice<Int, Real>::~IntLattice() {
	// kill();
	this->m_basis.IntMat::kill();              // Ok ?
	this->m_dualbasis.IntMat::kill();
	this->m_vecNorm.kill();
	this->m_dualvecNorm.kill();
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::copyLattice(
		const IntLattice<Int, Real> &lat) {
	this->m_modulo = lat.m_modulo;
	this->m_maxDim = lat.m_maxDim;
	this->m_dim = lat.m_dim;
	this->m_basis = IntMat(lat.m_basis);
	this->m_norm = lat.m_norm;
	this->m_vecNorm = RealVec(lat.m_vecNorm);
	this->m_withDual = lat.m_withDual;
	this->m_dualbasis = IntMat(lat.m_dualbasis);
	this->m_dualvecNorm = RealVec(lat.m_dualvecNorm);
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::overwriteLattice (
		const IntLattice<Int, Real> &lat, long dim) {
	if (dim <= this->m_maxDim) {
		CopyMatr(this->m_basis, lat.m_basis, dim);
		CopyVect(this->m_vecNorm, lat.m_vecNorm, dim);
		this->m_withDual = lat.m_withDual;
		if (this->m_withDual) {
			// We want to avoid resizing!
			// this->m_dualbasis.resize(this->m_basis.size1(),
			//		this->m_basis.size1());
			// this->m_dualvecNorm.resize(this->m_basis.size1());
			CopyMatr(this->m_dualbasis, lat.m_dualbasis, dim);
			CopyVect(this->m_dualvecNorm, lat.m_dualvecNorm, dim);
		}
		this->m_modulo = lat.m_modulo;
	}
	else
		std::cout << "IntLattice::overwriteLattice: dim > m_maxDim"
				<< std::endl;
	}


//===========================================================================

/*
template<typename Int, typename Real>
void IntLattice<Int, Real>::init() {
	int64_t dim = m_dim;
	this->setNegativeNorm();
}
*/

//===========================================================================

//template<typename Int, typename Real>
//void IntLattice<Int, Real>::initProj() {
//	// Reserves space for the projections in up to the dimension dim of the full lattice.
//	int64_t dim = m_dim;
//	this->setNegativeNorm();
//	this->m_basisProj.resize(dim, dim);   // Basis of current projection.
//	if (this->m_withDual) {
//		this->m_dualbasisProj.resize(dim, dim);
//		// double temp;   // Used only for m_lgVolDual2.
//		// NTL::conv (temp, this->m_modulo);
//		// m_lgVolDual2 = new double[dim+1];
//		// m_lgm2 = 2.0 * Lg (temp);
//		// m_lgVolDual2[1] = m_lgm2;
//	}
//}

//===========================================================================

template<typename Int>
class BasisConstruction;

/*  We definitely do not want this function as it is.    ****** */
//  Takes a projection of `*lattice` and puts it in the present object.

template<typename Int, typename Real>
void IntLattice<Int, Real>::buildProjection(
		IntLattice<Int, Real> *lattice, const Coordinates &proj) {
	const int64_t dim = this->getDim();
	int64_t i = 0;
	// We create two new matrices each time we build a projection!!!  Not good.   ******
	// We do not use the matrices initialized by `initProj`  ???   **********
	IntMat temp, temp2;
	temp.SetDims(dim, dim);  // dim of current lattice. We resize two objects here!  BAD!  *****
	temp2.SetDims(dim, dim);
	for (auto iter = proj.begin(); iter != proj.end(); ++iter) {
		// iter runs over the retained columns for the projection.
		//  What if a column number  *iter  exceeds dim-1  ????
		for (int64_t j = 0; j < dim; j++) {
			temp(j, i) = this->m_basis(j, (*iter));
		}
		++i;
	}
	// The generating vectors of proj are now in temp.
	// We construct a triangular basis for the projection `lattice` and put it in temp2.
	// The dimension of this projection is assumed to be the projection size,
	// so `temp2` will be a square invertible matrix.
	lattice->setDim(static_cast<int64_t>(proj.size()));  // Changes the lattice dimension!
	// lattice->m_order = m_order;
	// BasisConstruction<Int> bc;
	BasisConstruction<Int>::upperTriangularBasis(temp, temp2, this->m_modulo);
	temp2.SetDims(lattice->getDim(), lattice->getDim());
	lattice->setNegativeNorm();
	lattice->m_basis = temp2;       //  Current basis is replaced by temp2.
	lattice->m_withDual = this->m_withDual;
	if (this->m_withDual) {
		BasisConstruction<Int>::mDualUpperTriangular(lattice->m_basis, lattice->m_dualbasis,
				this->m_modulo);
		lattice->setDualNegativeNorm();
	}
}


//=========================================================================

/* **** The following is not needed since we can just call `projectMatrix` directly.  Right?

template<typename Int, typename Real>
void IntLattice<Int, Real>::getProjBasis (const Coordinates & coord, IntMat projBasis) {	
	 BasisConstruction<Int>::projectMatrix(this->getBasis(), projBasis, coord);
}

//=========================================================================

//  **** Warning: We cannot just project the m-dual basis on a set of coordinates!!!!!
//       See the guide.
template<typename Int, typename Real>
void IntLattice<Int, Real>::getProjBasisDual (const Coordinates & coord, IntMat projDualBasis) {
	 BasisConstruction<Int>::projectMatrix(this->getDualBasis(), projDualBasis, coord);
}

//=========================================================================

template<typename Int, typename Real>
void IntLattice<Int, Real>::getProjBasisPrimalDual (const Coordinates & coord, IntMat projBasis, 
		IntMat projDualBasis) {
	BasisConstruction<Int>::projectMatrix(this->getBasis(), projBasis, coord);
	BasisConstruction<Int>::projectMatrix(this->getDualBasis(), projDualBasis, coord);
}
*/


/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::setNegativeNorm() {
	for (int64_t i = 0; i < this->m_dim; i++) {
		this->m_vecNorm[i] = -1;
	}
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::setDualNegativeNorm() {
	for (int64_t i = 0; i < this->m_dim; i++) {
		this->m_dualvecNorm[i] = -1;
	}
}

/*=========================================================================*/

template<typename Int, typename Real>
IntLattice<Int, Real>::updateVecNorm() {
	updateVecNorm(0);
}

/*=========================================================================*/

template<typename Int, typename Real>
IntLattice<Int, Real>::updateVecNorm(const int64_t &d) {
	assert(d >= 0);
	for (int64_t i = d; i < this->m_dim; i++) {
		NTL::matrix_row<IntMat> row(this->m_basis, i);
		if (this->m_norm == L2NORM) {
			ProdScal<Int>(row, row, this->m_dim, this->m_vecNorm[i]);
		} else {
			CalcNorm<IntVec, Real>(row, this->m_dim, this->m_vecNorm[i],
					this->m_norm);
		}
	}
}

/*=========================================================================*/

template<typename Int, typename Real>
double IntLattice<Int, Real>::updateDualVecNorm() {
	return updateDualVecNorm(0);
}

/*=========================================================================*/

template<typename Int, typename Real>
double IntLattice<Int, Real>::updateDualVecNorm(const int64_t &d) {
	assert(d >= 0);
	assert(this->m_withDual);
	for (int64_t i = d; i < this->m_dim; i++) {
		NTL::matrix_row<IntMat> row(this->m_dualbasis, i);
		if (this->m_norm == L2NORM) {
			ProdScal<Int>(row, row, this->m_dim, this->m_dualvecNorm[i]);
		} else {
			CalcNorm<IntVec, Real>(row, this->m_dim, this->m_dualvecNorm[i],
					this->m_norm);
		}
	}
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::updateDualVecNorm(const int64_t &d, const int64_t &c) {
	assert(d >= 0);
	assert(this->m_withDual);
	for (int64_t i = 0; i < d+1; i++) {
		NTL::matrix_row<IntMat> row(this->m_dualbasis, i);
		if (this->m_norm == L2NORM) {
			ProdScal<Int>(row, row, c, this->m_dualvecNorm[i]);
		} else {
			CalcNorm<IntVec, Real>(row, c, this->m_dualvecNorm[i], this->m_norm);
		}
	}
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::updateSingleDualVecNorm(const int64_t &d, const int64_t &c) {
	assert(d >= 0);
	assert(this->m_withDual);
	NTL::matrix_row<IntMat> row(this->m_dualbasis, d);
	if (this->m_norm == L2NORM) {
		ProdScal<Int>(row, row, c, this->m_dualvecNorm[d]);
	} else {
		CalcNorm<IntVec, Real>(row, c, this->m_dualvecNorm[d], this->m_norm);
	}
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::updateScalL2Norm(const int64_t i) {
	NTL::matrix_row<IntMat> row(this->m_basis, i);
	ProdScal<Int>(row, row, this->m_dim, this->m_vecNorm[i]);
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::updateScalL2Norm(const int64_t k1,
		const int64_t k2) {
	for (int64_t i = k1; i < k2; i++) {
		updateScalL2Norm(i);
	}
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::updateDualScalL2Norm(const int64_t i) {
	NTL::matrix_row<IntMat> row(this->m_dualbasis, i);
	ProdScal<Int>(row, row, this->m_dim, this->m_dualvecNorm[i]);
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::updateDualScalL2Norm(const int64_t k1,
		const int64_t k2) {
	for (int64_t i = k1; i < k2; i++) {
		updateDualScalL2Norm(i);
	}
}


//=========================================================================

// ****  We cannot return a `double` here, this is usually too small!

template<typename Int, typename Real>
double IntLattice<Int, Real>::getShortestLengthBasis() {
   double out;
   Real temp;
   this->updateVecNorm(0);   // This function could return the index of the smallest!   *****
   temp = this->getVecNorm(0);
   for (int i = 1; i < this->getBasis().NumRows(); i++) {
      if (this->getVecNorm(i) < temp) { temp = this->getVecNorm(i);}
   }
   NTL::conv(out, temp);
   // if (this->getNormType()==L2NORM) out = sqrt(out);
   return out;
}

//=========================================================================

template<typename Int, typename Real>
double IntLattice<Int, Real>::getShortestLengthDualBasis() {
   double out;
   Real temp;
   this->updateVecNorm(0);
   temp = this->getDualVecNorm(0);
   for (int i = 1; i < this->getDualBasis().NumRows(); i++) {
      if (this->getDualVecNorm(i) < temp) { temp = this->getDualVecNorm(i);}
   }
   NTL::conv(out, temp);
   // if (this->getNormType()==L2NORM) out = sqrt(out);
   return out;
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::permute(int64_t i, int64_t j) {
	if (i == j)
		return;
	for (int64_t k = 0; k < this->m_dim; k++) {
		swap9(this->m_basis(j, k), this->m_basis(i, k));
		if (this->m_withDual) {
			swap9(this->m_dualbasis(j, k), this->m_dualbasis(i, k));
		}
	}
	swap9(this->m_vecNorm[i], this->m_vecNorm[j]);
	if (this->m_withDual) {
		swap9(this->m_dualvecNorm[i], this->m_dualvecNorm[j]);
	}
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::permutePrimal(int64_t i, int64_t j) {
	if (i == j)
		return;
	for (int64_t k = 0; k < this->m_dim; k++) {
		swap9(this->m_basis(j, k), this->m_basis(i, k));
	}
	swap9(this->m_vecNorm[i], this->m_vecNorm[j]);
}

//===========================================================================

template<typename Int, typename Real>
    void IntLattice<Int, Real>::dualize () {
    std::swap (this->m_basis, this->m_dualbasis);
    std::swap (this->m_vecNorm, this->m_dualvecNorm);
    // std::swap (this->m_withPrimal, this->m_withDual);
  }

/*=========================================================================*/

template<typename Int, typename Real>
bool IntLattice<Int, Real>::checkDuality() {
	if (!this->m_withDual) {
		std::cout << "Calling IntLattice::checkDuality with undefined m-dual"
				<< std::endl;
		return false;
	}
	Int S;
	int64_t dim = getDim();
	for (int64_t i = 0; i < dim; i++) {
		for (int64_t j = 0; j < dim; j++) {
			NTL::matrix_row<const IntMat> row1(this->m_basis, i);
			NTL::matrix_row<const IntMat> row2(this->m_dualbasis, j);
			ProdScal<Int>(row1, row2, dim, S);
			if (j != i) {
				if (S != 0) {
					std::cout << "******  checkDuality failed for V[" << i
							<< "] and W[" << j << "]" << std::endl;
					return false;
				}
			} else if (S != this->m_modulo) {
				std::cout << "******  checkDuality failed for i, j = " << i
						<< " , " << j << std::endl;
				return false;
			}
		}
	}
	return true;

}



/*=========================================================================*/

/*
 * We assume that the (square) lengths are already updated.
 * This gives flexibility to the user to put something else than
 * the square Euclidean length in vecNorm.
 */
template<typename Int, typename Real>
void IntLattice<Int, Real>::sortBasis (int64_t d) {
	int64_t dim = getDim();
	for (int64_t i = 0; i < dim; i++) {
		if (getVecNorm(i) < 0) {
			std::cout << "\n***** ERROR: in sort, Negative norm for i = " << i
					<< ",  dim = " << dim << std::endl;
		}
	}
	// This is a rather inefficient sort, in O(d^2) operations! 
	for (int64_t i = d; i < dim; i++) {
		int64_t k = i;
		for (int64_t j = i + 1; j < dim; j++) {
			if (getVecNorm(j) < getVecNorm(k))
				k = j;
		}
		if (i != k)
			permute(i, k);
	}
}

/*=========================================================================*/

/*
 * We assume that the (square) lengths are already updated.
 * This gives flexibility to the user to use something else than
 * the square Euclidean length.
 */
template<typename Int, typename Real>
void IntLattice<Int, Real>::sortBasisPrimal (int64_t d) {
	int64_t dim = getDim();
	for (int64_t i = 0; i < dim; i++) {
		if (getVecNorm(i) < 0) {
			std::cout << "\n***** ERROR: sort   Negative norm for i = " << i
					<< ",  dim = " << dim << std::endl;
		}
	}
	// This sort takes O(d^2) operations! 
	for (int64_t i = d; i < dim; i++) {
		int64_t k = i;
		for (int64_t j = i + 1; j < dim; j++) {
			if (getVecNorm(j) < getVecNorm(k))
				k = j;
		}
		if (i != k)
			permuteNoDual(i, k);
	}
}

/*=========================================================================*/

template<typename Int, typename Real>
std::string IntLattice<Int, Real>::toString() const {
	std::ostringstream os;
	os << "Dim = " << this->m_dim << " \n \n";
	os << std::setprecision(10) << "Primal basis vectors:\n";
	for (int64_t i = 0; i < this->m_dim; i++) {
		os << this->m_basis[i];
		//for (int64_t j = 0; j < this->m_dim; j++) {
		//  os <<  this->m_basis(i,j);
		//}
		os << "\n";
	}
	os << "\nm-Dual basis vectors:\n";
	for (int64_t i = 0; i < this->m_dim; i++) {
		if (this->m_withDual) {
			os << this->m_dualbasis[i];
			//for (int64_t j = 0; j < this->m_dim; j++) {
			//  os << this->m_dualbasis(i,j);
			//}
			os << "\n";
		}
	}
	os << "\n";
	os << "Norm used: " << toStringNorm(this->m_norm) << "\n"
			<< std::endl;
	os << "Norm of each Basis vector: \n";
	os << "Primal";
	if (this->m_withDual)
		os << "\t\tDual\n";
	os << "\n";

	for (int64_t i = 0; i < this->m_dim; i++) {
		if (this->m_vecNorm[i] < 0) {
			os << "NaN OR Not computed";
		} else {
			if (this->m_norm == L2NORM) {
				os << NTL::sqrt(this->m_vecNorm[i]);
			} else {
				os << this->m_vecNorm[i];
			}
		}
		os << "\t";
		if (this->m_withDual) {
			if (this->m_dualvecNorm[i] < 0)
				os << "NaN OR Not computed";
			else {
				if (this->m_norm == L2NORM) {
					os << NTL::sqrt(this->m_dualvecNorm[i]);
				} else {
					os << this->m_dualvecNorm[i];
				}
			}
		}
		os << "\n";
	}
	os << std::endl;
	return os.str();
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::write() const {
	std::cout << this->toString() << "\n";
}

/*=========================================================================*/

template<typename Int, typename Real>
std::string IntLattice<Int, Real>::toStringBasis() const {
	std::ostringstream os;
	os << "Primal Basis:\n";
	os << "  Dim = " << this->m_dim << " \n";
	for (int64_t i = 0; i < this->m_dim; i++) {
		os << "    [";
		for (int64_t j = 0; j < this->m_dim; j++)
			os << " " << std::setprecision(15) << this->m_basis(i, j);
		os << " ]\n";
	}
	os << "  Norms:\n";
	os << "    [";
	for (int64_t i = 0; i < this->m_dim; i++) {
		if (this->m_vecNorm[i] < 0) {
			os << "-1" << " ";
		} else {
			os << this->m_vecNorm[i] << " ";
		}
	}
	os << "]" << std::endl;
	return os.str();
}

/*=========================================================================*/

template<typename Int, typename Real>
std::string IntLattice<Int, Real>::toStringDualBasis() const {
	std::ostringstream os;
	os << "m-Dual Basis:\n";
	os << "  Dim = " << this->m_dim << " \n";
	for (int64_t i = 0; i < this->m_dim; i++) {
		os << "    [";
		for (int64_t j = 0; j < this->m_dim; j++)
			os << " " << std::setprecision(15) << this->m_dualbasis(i, j);
		os << " ]\n";
	}
	os << "  Norms:\n";
	os << "    [";
	for (int64_t i = 0; i < this->m_dim; i++) {
		if (this->m_dualvecNorm[i] < 0) {
			os << "-1" << " ";
		} else {
			os << this->m_dualvecNorm[i] << " ";
		}
	}
	os << "]" << std::endl;
	return os.str();
}

template class IntLattice<std::int64_t, double> ;
template class IntLattice<NTL::ZZ, double> ;
template class IntLattice<NTL::ZZ, NTL::RR> ;

} // namespace LatticeTester

#endif
