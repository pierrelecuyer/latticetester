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
    IntLatticeExt(Int m, int64_t maxDim, bool withPrimal = false,
            bool withDual = false, NormType norm = L2NORM);

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
    virtual void buildBasis(int64_t dim) {
    }
    ;

    /**
     * Similar to `buildBasis`, but builds only the m-dual basis, in `dim` dimensions.
     * This `dim` must not exceed `maxDim`. The flag  `withDual` is assumed to be true.
     */
    virtual void buildDualBasis(int64_t dim) {
    }
    ;

    /**
     * Increments the dimension `dim` of the basis by one. One coordinate is added
     * to each basis vector and one new basis vector is added as well.
     * Usually, the other coordinates are left unchanged.
     * If `withDual` is true, then the dimension of the m-dual is also increased by 1.
     */
    virtual void incDimBasis() {
    }
    ;

    /**
     * Similar to `incDimBasis`. Increments the dimension `dim` by 1, but only for the dual basis.
     */
    virtual void incDimDualBasis() {
    }
    ;

    /**
     * Returns a string describing this lattice.
     */
    virtual std::string toString() const {
        return "";
    }
    ;

protected:

    /**
     * \copydoc LatticeTester::IntLattice::kill()
     *  ** USEFUL ? **
     */
    virtual void kill();

};
// Class IntLatticeExt

//===========================================================================

template<typename Int, typename Real>
IntLatticeExt<Int, Real>::IntLatticeExt(Int m, int64_t maxDim, bool withPrimal,
        bool withDual, NormType norm) :
        IntLattice<Int, Real>(m, maxDim, withPrimal, withDual, norm) {
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
}
//===========================================================================

template<typename Int, typename Real>
IntLatticeExt<Int, Real>::IntLatticeExt(const IntLatticeExt<Int, Real> &lat) :
        IntLattice<Int, Real>(lat) {
    // this->m_withDual = lat.m_withDual;
    if (this->m_withPrimal)
        this->setNegativeNorm();
    if (this->m_withDual)
        this->setDualNegativeNorm();
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
    this->m_modulo = lat.m_modulo;
    if (lat.m_withPrimal)
        this->m_basis = lat.m_basis;
    if (lat.m_withDual)
        this->m_dualbasis = lat.m_dualbasis;
}

//===========================================================================

template class IntLatticeExt<std::int64_t, double> ;
template class IntLatticeExt<NTL::ZZ, double> ;
template class IntLatticeExt<NTL::ZZ, NTL::RR> ;

} // End namespace LatticeTester

#endif
