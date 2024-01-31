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

#ifndef LATTICETESTER_INTLATTICEEXT_H
#define LATTICETESTER_INTLATTICEEXT_H

#include "latticetester/IntLattice.h"
/*
 #include "latticetester/NormaBestLat.h"
 #include "latticetester/NormaBestBound.h"
 #include "latticetester/NormaLaminated.h"
 #include "latticetester/NormaRogers.h"
 #include "latticetester/NormaMinkL1.h"
 #include "latticetester/NormaPalpha.h"
 #include "latticetester/NormaMinkL2.h"
 #include "latticetester/Normalizer.h"
 */
#include "latticetester/Coordinates.h"
#include "latticetester/Lacunary.h"
#include "latticetester/Util.h"
#include "latticetester/BasisConstruction.h"
#include <cassert>

namespace LatticeTester {

/**
 * This abstract class extends `IntLattice` and is a skeleton for the
 * specialized subclasses that define specific types of lattices.
 * It is not intended to be used directly, but only via subclasses.
 * An `IntLatticeExt` object is an `IntLattice` with additional (virtual) methods
 * that must be implemented in subclasses.
 *
 * There are virtual methods to construct a basis and/or an m-dual basis of the full lattice,
 * to construct a basis and/or an m-dual basis for the projection of the full lattice on a subset
 * \f$I\f$  of coordinates indices specified by a `Coordinates` object,
 * and to extend the current basis and/or its m-dual by one coordinate.
 *
 * ***   Is the following still used?   ***
 * An `IntLatticeExt` object keeps a current projection of the whole lattice on a subset
 * of coordinates, as well as a basis (and perhaps its m-dual basis) for this current projection.
 * When computing figures of merit, this projection is frequently changed and the basis must
 * be updated.  The method `buildProjection` takes care of that.
 *
 * **REMOVE ?**
 * (I do not think we need to have the following here, but only where we compute the FOMs.)
 *
 * The lattices considered here are assumed to have a special structure, which is used
 * for the computation of the lattice density and the normalization constants in the
 * figures of merit.  It is assumed that the lattice has rank \f$k\f$ and that the
 * rescaling was done by multiplying the primal basis vectors by \f$m\f$.
 * All the lattices considered in the LatMRG and LatNet Builder software tools
 * have this property.
 * A lattice of rank \f$k\f$ with integer vectors modulo \f$m\f$ contains
 * \f$m^k\f$ distinct vectors (modulo $m$). If we divide the basis vectors by \f$m\f$,
 * this gives \f$m^k\f$ vectors per unit of volume, so \f$m^k\f$ is the density of the
 * original (non-scaled) lattice. This number is used to obtain bounds on the shortest vector length,
 * which are used to normalize the shortest vector length in the spectral test.
 * This class offers methods to compute and store the constants
 * \f$ \log_2(m^{2i}) \f$ for \f$ 1 \leq i \leq k \f$,
 * which are used elsewhere to compute the normalization constants.
 */

template<typename Int, typename Real>
class IntLatticeExt: public IntLattice<Int, Real> {
private:
	typedef NTL::vector<Int> IntVec;
	typedef NTL::matrix<Int> IntMat;
	typedef NTL::vector<Real> RealVec;
	typedef NTL::matrix<Real> RealMat;

public:

	/**
	 * A constructor that initializes the primal and dual bases with the
	 * identity matrix. The dimension of the lattice is set to `maxDim`
	 * and the norm type is set to `norm`.
	 * @param m The scaling factor `m` for the integer coordinates
	 * @param maxDim The maximal dimension for which this lattice can be
	 * expanded/tested
	 * @param withDual Specifies whether this object contains a dual or not
	 * @param norm  The norm type to measure the vector lengths.
	 */
	IntLatticeExt(Int m, int64_t maxDim, bool withDual=false, NormType norm = L2NORM);

	/**
	 * Copy constructor that makes a copy of `lat`. The maximal dimension
	 * of the created basis is set equal to the current dimension in `lat`.
	 */
	IntLatticeExt(const IntLatticeExt<Int, Real> &lat);

	/**
	 * Makes a copy of `lattice` into this object. It copies the
	 * internal vectors and matrices using the NTL assign operator =
	 * (see https://libntl.org/doc/matrix.cpp.html).
	 */
	void copy(const IntLatticeExt<Int, Real> &lat);

	/**
	 * Destructor. Depends on the specific subclass.
	 */
	virtual ~IntLatticeExt();

	/**
	 * This returns the rank (order) of the lattice.
	 */
	// int64_t getOrder() const { return m_order; }

	/**
	 * This virtual method builds a basis for the lattice in `dim` dimensions.
	 * This `dim` must not exceed `maxDim`.
	 *
	 * WHERE IS THE BASIS STORED?  DOES THIS CREATE A NEW OBJECT?   ******
	 */
	virtual void buildBasis(int64_t dim);

 	/** THIS IS FOR TESTING ONLY (CW)
	 * This virtual method builds a basis for the lattice in `dim` dimensions.
	 * This `dim` must not exceed `maxDim`. In contrast to buildBasis above,
	 * the `IntMat` object that holds the basis will have
	 * dimensions 'maxDim' x 'maxDim' and the entries that exceed 'dim' are set to 0.
	 * This function must be implemented in subclasses.
	 */
	virtual void buildBasisFullMatrix(int64_t dim);

	/**
	 * This virtual method builds only the dual basis for the lattice in `dim` dimensions.
	 * This `dim` must not exceed `maxDim`. buildDualBasis(d) does nothing, it must be
	 * implemented in subclasses.
	 */
	virtual void buildDualBasis(int64_t dim);

	/** THIS IS FOR TESTING ONLY (CW)
	 * This virtual method builds a basis for the dual lattice in `dim` dimensions.
	 * This `dim` must not exceed `maxDim`. In contrast to buildDualBasis, the basis matrix
	 * has dimension 'maxDim' x 'maxDim' and the entries which exceed 'd' are set to 0.<
	 * The fucntion must be implemented in subclasses.
	 */
	virtual void buildDualBasisFullMatrix(int64_t dim);

        /**
         * This virtual method builds only the dual basis in 'dim' dimensions while setting
         * the number of columns to fixed a value 'c'. It must be implemented in subclasses.
         */
    virtual void buildDualBasis (int64_t d, int64_t c);

	/**
	 * Builds a basis for the projection of the lattice over the coordinates in `proj`
	 * and returns it in `projBasis`.
	 */
	virtual void buildBasisProj (IntMat &projBasis, const Coordinates &proj);

	/**
	 * Builds a basis for the m-dual of the projection of the lattice over the coordinates in `proj`
	 * and returns it in `projBasis`.
	 */
	virtual void buildDualBasisProj (IntMat &projBasis, const Coordinates &proj);

    /**
	 * Increments the dimension of the basis and dual basis vectors by one.
	 * This implementation initializes the added components to `0` and does not
	 * compute the value taken by the added components and vector. It also
	 * resets the norm vectors to -1. The implementation in this
	 * class is meant to be overridden by subclasses.
	 *
	 * I THINK THIS SHOULD HAVE NO IMPLEMENTATION AT ALL!
	 * The current implementation does not increase the dimension at all.
	 * It is really used somewhere?
	 */
	virtual void incDimBasis();
 
    /** THIS IS FOR TESTING ONLY (CW)
         * Increases the current dimension of only the (primal) lattice basis by 1
         * under the assumption the dual basis matrix has dimension 'maxDim' x 'maxDim'. 
         * while the dimension of the basis is 'd'-1. This implementation is meant to be overridden 
	 * by subclasses as well.
     */
    virtual void incDimBasisFullMatrix (int64_t d);

	/**
	 * Increments the dimension of only the dual basis vectors by one.  
	 * This implementation works as incDim and is meant to be overridden 
	 * by subclasses as well.
	 */
	virtual void incDimDualBasis();
	
        /**
         * Increases the current dimension of only the dual lattice basis by 1
         * while fixing the number of columns to 'c'. Note that dim + 1 <= c <= maxDim
         * must hold. This implementation is meant to be overridden 
	     * by subclasses as well.
         */
    virtual void incDimDualBasis (int64_t c);
    
        /** THIS IS FOR TESTING ONLY (CW)
         * Increases the current dimension of only the dual lattice basis by 1
         * under the assumption the dual basis matrix has dimension 'maxDim' x 'maxDim'. 
         * while the dimension of the basis is 'd'-1. This implementation is meant to be overridden 
	 * by subclasses as well.
         */
    virtual void incDimDualBasisFullMatrix (int64_t d);

	/**
	 * Computes and stores the logarithm of the normalization factors
	 * (<tt>m_lgVolDual2</tt>) in all dimensions up to `MaxDim`, for this
	 * lattice. Here, `lgm2` must be \f$\log_g m^2\f$ and the computed values are
	 * those returned by `getLgVolDual2` below.
	 * These normalization contants are for the Euclidean norm only.
	 */
	// void calcLgVolDual2 (double lgm2);
	/**
	 * Returns \f$\log m^{2i}\f = i \log m^2$  for \f$1\le i \le k\f$,
	 * and \f$\log m^{2k}\f$ otherwise, where \f$k\f$ is the lattice rank (or order).
	 */
	// double getLgVolDual2 (int64_t i) const { return m_lgVolDual2[i]; }
	/**
	 * REMOVE: This method precomputes the log of the lattice density (or of its dual),
	 * as a function of the dimension.  These values are part of the normalization constants used to get
	 * the normalized merit from the shortest vectors in the lattice. If
	 * `dualF` is `true`, the values are computed for the m-dual
	 * lattice, otherwise they are computed for the primal lattice.
	 *
	 * Currently, this only computes the log of m^(k/dim) or its inverse.
	 * ** Done only once in a search? **
	 */
	// void computeNormalConstants(bool dualF);
	//    void fixLatticeNormalization (bool dualF);


	/**
	 * REMOVE: This depends on the lattice only via the density.
	 * This method is used nowhere else inside lattice tester.
	 * In LatMRG, it is used only once in progs/Test.h
	 * In Latnet builder, it seems to be used nowhere.
	 *
	 * Creates and returns a Normalizer corresponding to the normalization
	 * type `norma` and the density of the current lattice.
	 * The argument `alpha` = \f$\alpha\f$ is used only for the
	 * \f$P_{\alpha}\f$ measure. For all other cases, it is unused.
	 * The returned Normalizer is returned but not stored in this object.
	 * It contains the complete normalization constants for the number of dimensions
	 * of this lattice.
	 */
	// LatticeTester::Normalizer * getNormalizer (NormaType norma,
	//    int64_t alpha, bool dualF);


	/**
	 * Selects and stores a vector of indices with lacunary values.
	 */
//	virtual void setLac(const Lacunary<Int>& lac);

	/**
	 * Returns a string describing the lattice.
	 */
	virtual std::string toString() const;

protected:

	/**
	 * \copydoc LatticeTester::IntLattice::kill()
	 *  ** USEFUL ? **
	 */
	virtual void kill();

	/**
	 * REMOVE?  The order (rank) of the basis. Only defined in certain subclasses.
	 */
	// int64_t m_order;

	/**
	 * A vector of normalization constants.  See `calcLgVolDual2`.
	 */
	// double *m_lgVolDual2;
	/**
	 * \f$\log_2 (m^2)\f$.
	 */
	// double m_lgm2;

};
// Class IntLatticeExt

//===========================================================================

template<typename Int, typename Real>
IntLatticeExt<Int, Real>::IntLatticeExt(Int m, int64_t maxDim, bool withDual, NormType norm) :
		IntLattice<Int, Real>(m, maxDim, withDual, norm) {
	this->m_maxDim = maxDim;
	// this->m_withDual = withDual;
	// m_order = k;
	// Reserves space for the lattice and the projections (via initProj()) in up to maxDim dimensions.
	//this->initProj();
	this->m_basis.resize(this->m_dim, this->m_dim);
	this->m_vecNorm.resize(this->m_dim);
	this->setNegativeNorm();
	if (withDual) {
		this->m_dualbasis.resize(this->m_dim, this->m_dim);
		this->m_dualvecNorm.resize(this->m_dim);
		this->setDualNegativeNorm();
	}
}

//===========================================================================

template<typename Int, typename Real>
IntLatticeExt<Int, Real>::IntLatticeExt(const IntLatticeExt<Int, Real> &lat) :
		IntLattice<Int, Real>(lat) {
	this->m_withDual = lat.m_withDual;
	// m_order = lat.m_order;
	//this->initProj();
	//this->m_basisProj = lat.m_basisProj;
	if (this->m_withDual) {
		this->setDualNegativeNorm();
		//this->m_dualbasisProj = lat.m_dualbasisProj;
	}
}

//===========================================================================

template<typename Int, typename Real>
void IntLatticeExt<Int, Real>::kill() {
	// IntLattice<Int, Real>::kill();
}

//===========================================================================

template<typename Int, typename Real>
IntLatticeExt<Int, Real>::~IntLatticeExt() {
	//this->m_basisProj.kill();
	//this->m_dualbasisProj.kill();
	// IntLatticeExt<Int, Real>::~IntLatticeExt();
}

//===========================================================================

template<typename Int, typename Real>
void IntLatticeExt<Int, Real>::copy(const IntLatticeExt<Int, Real> &lat) {
	// Uses the NTL assignment operator = to make a copy of the bases.
	// m_order = lat.getOrder();
	this->m_modulo = lat.m_modulo;
	//m_m2 = lat.m_m2;
	this->m_basis = lat.m_basis;
	if (lat.m_withDual)
		this->m_dualbasis = lat.m_dualbasis;
	//this->initProj();  // Resizes the bases to the dimensions of the current object.
}

//===========================================================================

template<typename Int, typename Real>
void IntLatticeExt<Int, Real>::buildBasis(int64_t d) {
	// To be re-implemented in subclasses.
	MyExit(1, " buildBasis(d) does nothing, it must be implemented in subclass");
	d++;  // eliminates compiler warning
}

//===========================================================================

template<typename Int, typename Real>
void IntLatticeExt<Int, Real>::buildBasisFullMatrix(int64_t d) {
	// To be re-implemented in subclasses.
	MyExit(1, " buildBasisFullMatrix(d) does nothing, it must be implemented in subclass");
	d++;  // eliminates compiler warning
}

//===========================================================================

template<typename Int, typename Real>
void IntLatticeExt<Int, Real>::buildDualBasis(int64_t d) {
	// To be re-implemented in subclasses.
	MyExit(1, " buildDualBasis(d) does nothing, it must be implemented in subclass");
	d++;  // eliminates compiler warning
}

//===========================================================================

template<typename Int, typename Real>
void IntLatticeExt<Int, Real>::buildDualBasis(int64_t d, int64_t c) {
	// To be re-implemented in subclasses.
	MyExit(1, " buildDualBasis(d) does nothing, it must be implemented in subclass");
	d++;  // eliminates compiler warning
}

//===========================================================================

template<typename Int, typename Real>
void IntLatticeExt<Int, Real>::buildDualBasisFullMatrix(int64_t d) {
	// To be re-implemented in subclasses.
	MyExit(1, " buildDualBasisFullMatrix(d) does nothing, it must be implemented in subclass");
	d++;  // eliminates compiler warning
}

//===========================================================================

/*
 *
 * I THINK THIS SHOULD HAVE NO IMPLEMENTATION AT ALL!
 * The current implementation does not increase the dimension at all.
 * It is really used somewhere?
 */
template<typename Int, typename Real> void
IntLatticeExt<Int, Real>::incDimBasis() {
	IntLatticeExt<Int, Real> lattmp(*this);
	int64_t dim = this->getDim();
	this->m_basis.resize(dim + 1, dim + 1);
	this->m_vecNorm.resize(dim + 1);
	if (this->m_withDual) {
		this->m_dualbasis.resize(dim + 1, dim + 1);
		this->m_dualvecNorm.resize(dim + 1);
	}
	for (int64_t i = 0; i < dim; i++) {
		for (int64_t j = 0; j < dim; j++) {
			this->m_basis(i, j) = lattmp.m_basis(i, j);
			if (this->m_withDual)
				this->m_dualbasis(i, j) = lattmp.m_dualbasis(i, j);
		}
		this->m_vecNorm(i) = lattmp.m_vecNorm(i);
		if (this->m_withDual)
			this->m_dualvecNorm(i) = lattmp.m_dualvecNorm(i);
	}
	this->setNegativeNorm(dim);
	if (this->m_withDual)
		this->setDualNegativeNorm(dim);
	this->setDim(dim + 1);
	return;
}

//===========================================================================

template<typename Int, typename Real>
void IntLatticeExt<Int, Real>::incDimBasisFullMatrix(int64_t d) {
	MyExit(1, " incDimBasisFullMatrix(d) does nothing, it must be implemented in subclass");
	d++;  // eliminates compiler warning
}

//===========================================================================

template<typename Int, typename Real>
void IntLatticeExt<Int, Real>::incDimDualBasis() {
	IntLatticeExt<Int, Real> lattmp(*this);
	int64_t dim = this->getDim();
	this->m_dualbasis.resize(dim + 1, dim + 1);
	this->m_dualvecNorm.resize(dim + 1);

	for (int64_t i = 0; i < dim; i++) {
		for (int64_t j = 0; j < dim; j++) {
				this->m_dualbasis(i, j) = lattmp.m_dualbasis(i, j);
		}
		this->m_dualvecNorm(i) = lattmp.m_dualvecNorm(i);
	}
	this->setDualNegativeNorm(dim);
	this->setDim(dim + 1);
	return;
}

//===========================================================================

template<typename Int, typename Real>
void IntLatticeExt<Int, Real>::incDimDualBasis(int64_t c) {
	MyExit(1, " buildDualBasis(d) does nothing, it must be implemented in subclass");
	c++;  // eliminates compiler warning
}

//===========================================================================

template<typename Int, typename Real>
void IntLatticeExt<Int, Real>::incDimDualBasisFullMatrix(int64_t d) {
	MyExit(1, " incDimDualBasisFullMatrix(d) does nothing, it must be implemented in subclass");
	d++;  // eliminates compiler warning
}

//===========================================================================
/**
 template<typename Int, typename Real>
 void IntLatticeExt<Int, Real>::calcLgVolDual2 (double lgm2)
 {
 if(!(this->m_withDual)) return;
 int64_t dim = this->getDim();
 int64_t rmax = std::min(m_order, dim);

 m_lgVolDual2[1] = lgm2;
 for (int64_t r = 2; r <= rmax; r++)
 m_lgVolDual2[r] = m_lgVolDual2[r - 1] + lgm2;
 // WARNING [David]: one version had `m_order` instead of `rmax`.
 for (int64_t r = rmax + 1; r <= dim; r++)
 m_lgVolDual2[r] = m_lgVolDual2[r - 1];
 }
 */

//===========================================================================
/**
 template<typename Int, typename Real>
 void IntLatticeExt<Int, Real>::computeNormalConstants(
 bool dualF)
 {
 // Normalization factor: dual to primal : m^(k/dim) -> 1/m^(k/dim)
 // This is the part of the normalization that depends on the lattice density.
 if (( dualF && m_lgVolDual2[1] < 0.0) ||
 (!dualF && m_lgVolDual2[1] > 0.0)) {
 for (int64_t i = 0; i < this->getDim(); i++)
 m_lgVolDual2[i] = -m_lgVolDual2[i];
 }
 //   for (int64_t i = 1; i <= getMaxDim(); i++)
 //      std::cout << " fix  " << m_lgVolDual2[i] << endl;
 }ss
 */

//===========================================================================

/**
 template<typename Int, typename Real>
 Normalizer<Real> * IntLatticeExt<Int, Real>::getNormalizer(
 NormaType norma, int64_t alpha, bool dualF)
 {
 int64_t dim = this->getDim();
 Normalizer<Real> *normal;

 Real logDensity;

 // The primal lattice density is assumed to be m^k, and m^{-k} for the dual.
 if (dualF) // dual basis
 logDensity = - m_order * NTL::log(this->m_modulo);
 else // primal basis
 logDensity = m_order * NTL::log(this->m_modulo);

 // We create a normalizer normal for the given density, and return it.
 // This normalizer is not stored in this object.
 switch (norma) {
 case BESTLAT:
 normal = new NormaBestLat<Real> (logDensity, dim);
 break;
 case BESTBOUND:
 normal = new NormaBestBound<Real> (logDensity, dim);
 break;
 case LAMINATED:
 normal = new NormaLaminated<Real> (logDensity, dim);
 break;
 case ROGERS:
 normal = new NormaRogers<Real> (logDensity, dim);
 break;
 case MINKL1:
 normal = new NormaMinkL1<Real> (logDensity, dim);
 break;
 case MINK:
 normal = new NormaMinkL2<Real> (logDensity, dim);
 break;
 case NONE:
 normal = new Normalizer<Real> (logDensity, dim, "Norma_generic");
 break;
 default:
 std::cout << "normalizer:   no such case";
 exit (2);
 }
 return normal;
 }
 */

//===========================================================================
template<typename Int, typename Real>
std::string IntLatticeExt<Int, Real>::toString() const {
	// To be re-implemented in subclasses.
	assert(0);
	return std::string();
}

//===========================================================================

template class IntLatticeExt<std::int64_t, double> ;
template class IntLatticeExt<NTL::ZZ, double> ;
template class IntLatticeExt<NTL::ZZ, NTL::RR> ;

} // End namespace LatticeTester

#endif
