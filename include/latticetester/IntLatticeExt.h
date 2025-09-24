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

#include "latticetester/Util.h"
#include "latticetester/IntLattice.h"
#include "latticetester/Coordinates.h"
#include "latticetester/BasisConstruction.h"
#include <cassert>

namespace LatticeTester {

/**
 * \class IntLatticeExt
 *
 * This abstract class extends `IntLattice` with additional (virtual) methods
 * that must be implemented in subclasses that define specific types of lattices.
 *
 * These virtual methods permit one to construct a basis and/or an m-dual basis of the lattice
 * and to extend the current basis and/or its m-dual by one coordinate.
 * The `IntLattice` base class already implements methods to construct a basis and/or
 * an m-dual basis for the projection of the full lattice on a subset
 * of coordinates indices specified by a `Coordinates` object.  These general default
 * implementations are often overridden in subclasses by faster specialized implementations.
 *
 * Recall that the lattices in `IntLattice` (and here) have `dim <= maxDim` dimensions and
 * the bases are stored in `IntMat` objects of dimensions `maxDim x maxDim`.
 */

template<typename Int, typename Real>
class IntLatticeExt: public IntLattice<Int, Real> {

public:

   /**
    * A constructor that reserves the space for the primal and m-dual bases and vector lengths,
    * but does not initialize them. The maximal dimension of the lattice is set to `maxDim`,
    * and the norm type is set to `norm`.
    * In the second version, the scaling factor for the integer coordinates is set to `m`.
    */
   IntLatticeExt(int64_t maxDim, NormType norm = L2NORM);

   IntLatticeExt(Int m, int64_t maxDim, NormType norm = L2NORM);

   /**
    * Destructor. Depends on the specific subclass.
    */
   virtual ~IntLatticeExt();

   /**
    * This virtual method builds a basis for the lattice in `dim` dimensions,
    * and store it in the upper-left part of the internal `m_basis` variable,
    * which is assumed to be an `IntMat` object of dimensions maxDim x maxDim.
    * The parameter `dim` must not exceed `maxDim`.
    */
   virtual void buildBasis(int64_t dim) {
   }
   ;

   /**
    * Similar to `buildBasis`, but builds only the m-dual basis, in `dim` dimensions.
    * This `dim` must not exceed `maxDim`.
    */
   virtual void buildDualBasis(int64_t dim) {
   }
   ;

   /**
    * Increments the dimension `dim` of the basis by one. One coordinate is added
    * to each basis vector and one new basis vector is added as well.
    * Usually, the other coordinates are left unchanged.
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

protected:

   /**
    * \copydoc LatticeTester::IntLattice::kill()
    */
   virtual void kill();

};
// Class IntLatticeExt

//============================================================================
// IMPLEMENTTION

//===========================================================================

template<typename Int, typename Real>
IntLatticeExt<Int, Real>::IntLatticeExt(int64_t maxDim, NormType norm) :
      IntLattice<Int, Real>(maxDim, norm) {
   this->m_basis.SetDims(this->m_maxDim, this->m_maxDim);
   this->m_vecNorm.SetLength(this->m_maxDim);
   this->setNegativeNorm();
   this->m_dualbasis.SetDims(this->m_maxDim, this->m_maxDim);
   this->m_dualvecNorm.SetLength(this->m_maxDim);
   this->setDualNegativeNorm();
}

template<typename Int, typename Real>
IntLatticeExt<Int, Real>::IntLatticeExt(Int m, int64_t maxDim, NormType norm) :
      IntLatticeExt<Int, Real>(maxDim, norm) {
   this->m_modulo = m;
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

template class IntLatticeExt<std::int64_t, double> ;
template class IntLatticeExt<NTL::ZZ, double> ;
template class IntLatticeExt<NTL::ZZ, xdouble> ;
template class IntLatticeExt<NTL::ZZ, quad_float> ;
template class IntLatticeExt<NTL::ZZ, NTL::RR> ;

} // End namespace LatticeTester

#endif
