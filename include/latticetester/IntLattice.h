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
 * We use only the upper-left `dim x dim` corner of each array to store the actual basis.
 * These arrays can then be allocated only once and never have to be resized, which improves speed.
 * A boolean variable `withPrimal` indicates if we maintain the primal basis and a variable
 * `withDual` indicates if we maintain the dual basis. At least one of them (or both) should be true.
 *
 * A norm is also chosen in `NormType` to measure the vector lengths; by default it is the
 * Euclidean norm.
 * Methods and attributes are offered to compute and store the norms of the basis and/or
 * the m-dual basis vectors, to permute these vectors, sort them by length, etc.
 * The norms are for the vectors of the rescaled lattice and its m-dual.
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
    IntLattice(const Int m, const int64_t maxDim, bool withPrimal = true,
            bool withDual = false, NormType norm = L2NORM);

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
    IntLattice(const IntMat primalbasis, const IntMat dualbasis, const Int m,
            const int64_t maxDim, NormType norm = L2NORM);

    /**
     * Copy constructor. Makes a deep copy of `lat` into `*this` new object.
     */
    IntLattice(const IntLattice<Int, Real> &lat);

    /**
     * Destructor.
     */
    virtual ~IntLattice();

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
     * Takes the projection of the lattice represented by the present object
     * over the set of coordinates determined by `proj`, finds a set of generating
     * vectors for that projection, and builds a basis and perhaps its m-dual basis
     * (if maintained) for the projection. After that, the `IntLattice` object
     * given by `projLattice` will represent the projection. It is assumed that the basis
     * of the present lattice is already available and contains all the coordinates in `proj`.
     * If `projLattice.withDual() == true`, the primal basis of the projection is obtained by
     * a triangularization method and a lower-triangular basis of the m-dual is also computed.
     * Otherwise, only the primal basis is computed, using LLL with the parameter `delta`.
     * The modulus `m` is assumed to be the same for `projLattice` and for the current object.
     * The `maxDim` of the `projLattice` object must be large enough so it can holds the projection.
     * The variables associated with the projection (dimension, norms, etc.) are also updated.
     * This method can be called several times with the same `projLattice` object
     * to examine several different projections.
     * Note that representing each projection as an `IntLattice` object is required when
     * we want to call `Reducer::shortestVector` for several projections.
     */
    virtual void buildProjection(IntLattice<Int, Real> *projLattice,
            const Coordinates &proj, double delta = 0.99);

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
        assert(this->m_maxDim == basis.NumRows());
        this->m_modulo = m;
        this->m_withPrimal = true;
        this->m_withDual = withDual;
        this->m_norm = norm;
        this->m_dim = dim;
        this->m_basis = basis;
        setNegativeNorm();
    }

    /**
     * Changes only the 'basis' and the current dimension to 'dim'.
     * Does not change `m` and `maxDim`.
     */
    void setBasis(const IntMat basis, const int64_t dim) {
        assert(this->m_maxDim == basis.NumRows());
        this->m_dim = dim;
        this->m_basis = basis;
        setNegativeNorm();
    }

    /**
     * Returns the m-dual basis represented in a matrix.
     */
    IntMat& getDualBasis() {
        return m_dualbasis;
    }

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
        assert(dim <= this->m_maxDim);
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
     * Updates the array containing the norms of the basis vectors by recomputing them.
     * We could also return the index of the vector that is shortest!
     */
    void updateVecNorm();

    /**
     * Updates the array containing the norm of the basis vectors from the `d`-th
     * component to the last, by recomputing them.
     * Putting `d=0` recomputes all the norms.
     */
    void updateVecNorm(const int64_t &d);

    /**
     * Updates the 'd'-th entry of the array containing the basis vectors norms.
     * Only the first c components are used for calculating the norm.
     * */
    void updateSingleVecNorm(const int64_t &d, const int64_t &c);

    
    /**
     * Updates the array containing the m-dual basis vectors norms by recomputing them.
     * Assumes that the dual basis is available.
     */
    void updateDualVecNorm();

    /**
     * Updates the 'd'-th entry of the array containing the m-dual basis vectors norms.
     * Only the first c components are used for calculating the norm.
     * */
    void updateSingleDualVecNorm(const int64_t &d, const int64_t &c);

    /**
     * Updates the array containing the m-dual basis vectors norms from the `d`-th
     * component to the last by recomputing them.
     * */
    void updateDualVecNorm(const int64_t &d);

    /**
     * Updates the array containing the m-dual basis vectors norms from the `d`-th
     * component to the last by recomputing them. Only the first c components are
     * used for calculating the norm
     * */
    void updateDualVecNorm(const int64_t &d, const int64_t &c);

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
    void dualize();

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
    void sortBasis(int64_t d);

    /**
     * Sorts the primal basis vectors with indices greater of equal to `d` by
     * increasing length. The m-dual vectors (if maintained) are **not** permuted.
     */
    void sortPrimalBasis(int64_t d);

    /**
     * Sorts the `m`-dual basis vectors with indices greater of equal to `d` by
     * increasing length. The primal basis vectors are **not** permuted.
     */
    void sortDualBasis(int64_t d);

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
     * This variable is `true` iff a primal basis is available.
     */
    bool m_withPrimal;

    /**
     * This variable is `true` iff an m-dual basis is available.
     */
    bool m_withDual;

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

};

// class IntLattice

//===========================================================================

template<typename Int, typename Real>
IntLattice<Int, Real>::IntLattice(const Int m, const int64_t maxDim,
        bool withPrimal, bool withDual, NormType norm) {
//		: m_modulo(m), m_dim(dim), m_withDual(withDual), m_norm(norm) {
    this->m_modulo = m;
    this->m_maxDim = maxDim;
    this->m_dim = maxDim;
    this->m_withPrimal = withPrimal;
    this->m_withDual = withDual;
    this->m_norm = norm;
    this->m_basis.resize(maxDim, maxDim);
    this->m_vecNorm.resize(maxDim);
    if (withDual) {
    	m_dualbasis.resize(maxDim, maxDim);
    	m_dualvecNorm.resize(maxDim);
    }
    setNegativeNorm();
}

//===========================================================================

template<typename Int, typename Real>
IntLattice<Int, Real>::IntLattice(const IntMat basis, const Int m,
        const int64_t maxDim, bool withDual, NormType norm) {
//		: m_basis(basis), m_modulo(m), m_dim(dim), m_withDual(withDual), m_norm(norm) {
    this->m_modulo = m;
    this->m_maxDim = maxDim;
    this->m_dim = maxDim;
    this->m_withPrimal = true;
    this->m_withDual = withDual;
    this->m_norm = norm;
    assert(basis.NumRows() == maxDim);
    this->m_basis = basis;
    this->m_vecNorm.resize(maxDim);
    if (withDual) {
    	m_dualbasis.resize(maxDim, maxDim);
    	m_dualvecNorm.resize(maxDim);
    }
    setNegativeNorm();
}

/*=========================================================================*/

template<typename Int, typename Real>
IntLattice<Int, Real>::IntLattice(const IntMat primalbasis,
        const IntMat dualbasis, const Int m, const int64_t maxDim,
        NormType norm) :
        IntLattice<Int, Real>(primalbasis, m, maxDim, true, norm) {
    assert(dualbasis.NumRows() == maxDim);
    this->m_withPrimal = true;
    this->m_withDual = true;
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
void IntLattice<Int, Real>::copyLattice(const IntLattice<Int, Real> &lat) {
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
void IntLattice<Int, Real>::overwriteLattice(const IntLattice<Int, Real> &lat,
        long dim) {
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
    } else
        std::cout << "IntLattice::overwriteLattice: dim > m_maxDim"
                << std::endl;
}

//===========================================================================
//		this->m_dualbasisProj.resize(dim, dim);
//		// double temp;   // Used only for m_lgVolDual2.
//		// NTL::conv (temp, this->m_modulo);
//		// m_lgVolDual2 = new double[dim+1];
//		// m_lgm2 = 2.0 * Lg (temp);
//		// m_lgVolDual2[1] = m_lgm2;
//	}
//}

//===========================================================================

template<typename Int>class BasisConstruction;  // Needed? CW: Yes! Otherwise it does not compile.

template<typename Int, typename Real>
void IntLattice<Int, Real>::buildProjection(IntLattice<Int, Real> *projLattice,
        const Coordinates &proj, double delta) {
    // We assume here that this and lattice have the same m.
    projLattice->setDim (proj.size());  // Number of coordinates in the projection.
    if (!projLattice->m_withDual) { // This builds only the primal basis.
        BasisConstruction<Int>::projectionConstructionLLL(this->m_basis,
                projLattice->m_basis, proj, this->m_modulo, delta, proj.size());//, CW
                //projLattice->m_vecNorm);
    } else { // This builds both the primal and the m-dual bases.
        BasisConstruction<Int>::projectionConstructionUpperTri(this->m_basis,
                projLattice->m_basis, proj, this->m_modulo, this->m_dim);
        BasisConstruction<Int>::mDualUpperTriangular(projLattice->m_basis,
                projLattice->m_dualbasis, this->m_modulo, projLattice->m_dim);
        this->setNegativeNorm();
        this->setDualNegativeNorm();
    }
}

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
void IntLattice<Int, Real>::updateVecNorm() {
    updateVecNorm(0);
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::updateVecNorm(const int64_t &d) {
    assert(d >= 0);
    for (int64_t i = d; i < this->m_dim; i++) {
        NTL::matrix_row<IntMat> row(this->m_basis, i);  // Is this making a copy of the row?  ******
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
void IntLattice<Int, Real>::updateDualVecNorm() {
    updateDualVecNorm(0);
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::updateSingleVecNorm(const int64_t &d,
        const int64_t &c) {
    assert(d >= 0);
    NTL::matrix_row<IntMat> row(this->m_basis, d);
    if (this->m_norm == L2NORM) {
        ProdScal<Int>(row, row, c, this->m_vecNorm[d]);
    } else {
        CalcNorm<IntVec, Real>(row, c, this->m_vecNorm[d], this->m_norm);
    }
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::updateDualVecNorm(const int64_t &d) {
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
void IntLattice<Int, Real>::updateDualVecNorm(const int64_t &d,
        const int64_t &c) {
    assert(d >= 0);
    assert(this->m_withDual);
    for (int64_t i = 0; i < d + 1; i++) {
        NTL::matrix_row<IntMat> row(this->m_dualbasis, i);
        if (this->m_norm == L2NORM) {
            ProdScal<Int>(row, row, c, this->m_dualvecNorm[i]);
        } else {
            CalcNorm<IntVec, Real>(row, c, this->m_dualvecNorm[i],
                    this->m_norm);
        }
    }
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::updateSingleDualVecNorm(const int64_t &d,
        const int64_t &c) {
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
void IntLattice<Int, Real>::dualize() {
    std::swap(this->m_basis, this->m_dualbasis);
    std::swap(this->m_vecNorm, this->m_dualvecNorm);
    bool temp = this->m_withPrimal;
    this->m_withPrimal = this->m_withDual;
    this->m_withDual = temp;
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
void IntLattice<Int, Real>::sortBasis(int64_t d) {
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
void IntLattice<Int, Real>::sortPrimalBasis(int64_t d) {
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
            permutePrimal(i, k);
    }
}

/*=========================================================================*/

inline double sqrtReal(const double &a) { return std::sqrt(a); }

inline NTL::RR sqrtReal(const NTL::RR &a) { return NTL::sqrt(a); }

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
    os << "Norm used: " << toStringNorm(this->m_norm) << "\n" << std::endl;
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
                os << sqrtReal(this->m_vecNorm[i]);
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
                    os << sqrtReal(this->m_dualvecNorm[i]);
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
