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
 * There are virtual methods to construct a basis and/or an m-dual basis of the full lattice
 * and to extend the current basis and/or its m-dual by one coordinate.
 * These methods must be implemented in the subclasses.
 * The `IntLattice` base class already implements methods to construct a basis and/or
 * an m-dual basis for the projection of the full lattice on a subset
 * of coordinates indices specified by a `Coordinates` object.  These general default
 * implementation can often be overridden by faster specialized implementations in the subclasses.
 *
 * Recall that the lattices in `IntLattice` (and here) have `dim < maxDim` dimensions and
 * the bases are stored in `IntMat` objects of dimensions `maxDim x maxDim`.
 *
 * ***  The following has been changed.   ***
 * An `IntLatticeExt` object keeps a current projection of the whole lattice on a subset
 * of coordinates, as well as a basis (and perhaps its m-dual basis) for this current projection.
 * When computing figures of merit, this projection is frequently changed and the basis must
 * be updated.  The method `buildProjection` takes care of that.
 * ****   This has been changed:  The projection must now be kept in a separate `IntLattice` object.
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
	IntLatticeExt(Int m, int64_t maxDim, bool withPrimal=false, bool withDual=false,
			      NormType norm = L2NORM);

	/**
	 * Copy constructor that makes a copy of `lat`. The maximal dimension
	 * of the created basis is set equal to the current maximal dimension in `lat`.
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
	 * This returns the rank (order) of the lattice.   Needed?
	 */
	// int64_t getOrder() const { return m_order; }

	/**
	 * This virtual method builds a basis for the lattice in `dim` dimensions,
	 * and store it in the upper-left part of the internal `m_basis` variable,
	 * which is assumed to be an `IntMat` object of dimensions maxDim x maxDim.
	 * The parameter `dim` must not exceed `maxDim`.
	 * If `withDual` is true, it also builds an m-dual basis.
	 */
	virtual void buildBasis(int64_t dim) {};

 	/** THIS IS FOR TESTING ONLY (CW)
	 * This virtual method builds a basis for the lattice in `dim` dimensions.
	 * This `dim` must not exceed `maxDim`. In contrast to buildBasis above,
	 * the `IntMat` object that holds the basis will have
	 * dimensions 'maxDim' x 'maxDim' and the entries that exceed 'dim' are set to 0.
	 * This function must be implemented in subclasses.
	 */
	// virtual void buildBasisFullMatrix(int64_t dim) {};

	/**
	 * Similar to `buildBasis`, but for the m-dual only.
	 * This virtual method builds only the m-dual basis for the lattice in `dim` dimensions.
	 * This `dim` must not exceed `maxDim`. The flag  `withDual` is assumed to be true.
	 */
	virtual void buildDualBasis(int64_t dim) {};

	/** THIS IS FOR TESTING ONLY (CW)
	 * This virtual method builds a basis for the dual lattice in `dim` dimensions.
	 * This `dim` must not exceed `maxDim`. In contrast to buildDualBasis, the basis matrix
	 * has dimension 'maxDim' x 'maxDim' and the entries which exceed 'd' are set to 0.<
	 * The fucntion must be implemented in subclasses.
	 */
	// virtual void buildDualBasisFullMatrix(int64_t dim) {};

    /**
     * **** This seems to be the same as the `buildDualBasis` above.
     *
     * This virtual method builds only the dual basis in 'dim' dimensions while setting
     * the number of columns to fixed a value 'c'. It must be implemented in subclasses.
     */
    // virtual void buildDualBasis (int64_t d, int64_t c) {};

    /**
	 * Increments the dimension `dim` of the basis by one. One coordinate is added
	 * to each basis vector and one new basis vector is added as well.
	 * Usually, the other coordinates are left unchanged.
	 * If `withDual` is true, then the dimension of the m-dual is also increased by 1.
	 */
	virtual void incDimBasis() {};
 
    /** THIS IS FOR TESTING ONLY (CW)
         * Increases the current dimension of only the (primal) lattice basis by 1
         * under the assumption the dual basis matrix has dimension 'maxDim' x 'maxDim'. 
         * while the dimension of the basis is 'd'-1. This implementation is meant to be overridden 
	 * by subclasses as well.
     */
    // virtual void incDimBasisFullMatrix (int64_t d) {};

	/**
	 * Similar to `incDimBasis`. Increments the dimension `dim` by 1, but only for the dual basis.
	 */
	virtual void incDimDualBasis() {};
	
        /**
         * **** Same as the other one, it seems.
         *
         * Increases the current dimension of only the dual lattice basis by 1
         * while fixing the number of columns to 'c'. Note that dim + 1 <= c <= maxDim
         * must hold. This implementation is meant to be overridden 
	     * by subclasses as well.
         */
    // virtual void incDimDualBasis (int64_t c) {};
    
        /** THIS IS FOR TESTING ONLY (CW)
         * Increases the current dimension of only the dual lattice basis by 1
         * under the assumption the dual basis matrix has dimension 'maxDim' x 'maxDim'. 
         * while the dimension of the basis is 'd'-1. This implementation is meant to be overridden 
	 * by subclasses as well.
     */
    // virtual void incDimDualBasisFullMatrix (int64_t d) {};

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
	 * Returns a string describing this lattice.
	 */
	virtual std::string toString() const { };

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
IntLatticeExt<Int, Real>::IntLatticeExt(Int m, int64_t maxDim, bool withPrimal, bool withDual, NormType norm) :
		IntLattice<Int, Real>(m, maxDim, withPrimal, withDual, norm) {
	// this->m_maxDim = maxDim;
	// m_order = k;
	// Reserves space for the lattice and the projections (via initProj()) in up to maxDim dimensions.
	//this->initProj();
	if (withPrimal) {
	    this->m_basis.resize(this->m_maxDim, this->m_maxDim);
	    this->m_vecNorm.resize(this->m_maxDim);
	    this->setNegativeNorm();
	if (withDual) {
		this->m_dualbasis.resize(this->m_maxDim, this->m_maxDim);
		this->m_dualvecNorm.resize(this->m_maxDim);
		this->setDualNegativeNorm();
	}
}

//===========================================================================

template<typename Int, typename Real>
IntLatticeExt<Int, Real>::IntLatticeExt(const IntLatticeExt<Int, Real> &lat) :
		IntLattice<Int, Real>(lat) {
	// this->m_withDual = lat.m_withDual;
	// m_order = lat.m_order;
	if (this->m_withPrimal)  this->setPrimalNegativeNorm();
	if (this->m_withDual)  this->setDualNegativeNorm();
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
	if (lat.m_withPrimal) this->m_basis = lat.m_basis;
	if (lat.m_withDual) this->m_dualbasis = lat.m_dualbasis;
}

//===========================================================================

template class IntLatticeExt<std::int64_t, double> ;
template class IntLatticeExt<NTL::ZZ, double> ;
template class IntLatticeExt<NTL::ZZ, NTL::RR> ;

} // End namespace LatticeTester

#endif
