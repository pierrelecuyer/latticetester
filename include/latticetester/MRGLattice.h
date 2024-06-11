// This file is part of LatticeTester.
//
// Copyright (C) 2012-2022  The LatticeTester authors, under the occasional supervision
// of Pierre L'Ecuyer at Universit? de Montr?al.
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

#ifndef LATTICETESTER_MRGLATTICE_H
#define LATTICETESTER_MRGLATTICE_H

#include <cassert>
#include "latticetester/Util.h"
#include "latticetester/EnumTypes.h"
#include "latticetester/IntLatticeExt.h"
#include "latticetester/BasisConstruction.h"

namespace LatticeTester {

/**
 * This subclass of `IntLatticeExt` defines an MRG lattce and is similar to Rank1Lattice,
 * but constructs lattices associated with multiple recursive generators (MRGs) with 
 * modulus m, order k, and vector of multipliers a = (a_1, . . . , a_k).
 *
 * The functions `buildBasis` and `incDimBasis` always build and update the primal basis.
 * To build the m-dual, use `buildDualBasis`.
 */

template<typename Int, typename Real>
class MRGLattice: public IntLatticeExt<Int, Real> {

private:
   typedef NTL::vector<Int> IntVec;
   typedef NTL::matrix<Int> IntMat;
   typedef NTL::vector<Real> RealVec;

public:

   /**
    * This constructor takes as input the modulus `m`, the vector of multipliers `aa`,
    * and the norm used to measure the vector lengths.
    * The maximal dimension `maxDim` will be the maximal dimension of the basis.
    * This constructor does not build the basis, to leave
    * more flexibility in the dimension when doing so.
    */
   MRGLattice(const Int &m, const IntVec &aa, int64_t maxDim, NormType norm = L2NORM);

   /*
    * Copy constructor.
    */
   MRGLattice(const MRGLattice<Int, Real> &Lat);

   /**
    * Assigns `Lat` to this object.
    */
   MRGLattice& operator=(const MRGLattice<Int, Real> &Lat);

   /**
    * Destructor.
    */
   ~MRGLattice();

   /**
    * Sets the vector of multipliers. The order of the lattice is set equal to
    * the length of this vector.
    */
   void setaa(const IntVec &aa);

   /**
    * Builds a basis in `dim` dimensions. This `dim` must not exceed `this->maxDim()`.
    * This initial primal basis will be upper triangular.
    */
   void buildBasis(int64_t dim);

   /**
    * Builds both the primal and an m-dual lower triangular basis directly
    * in `dim` dimensions.  This `dim` must not exceed `maxDim`.
    */
   void buildDualBasis(int64_t dim);

   /**
    * Increases the current dimension of the primal basis by 1 and updates the basis.
    * The new increased dimension must not exceed `maxDim`.
    */
   void incDimBasis();

   /**
    * Increases the current dimension of both the primal and m-dual basis by 1.
    * The new increased dimension must not exceed `maxDim`.
    * This function uses the simplified method for MRG lattices given in the lattice tester guide.
    * It requires the original primal basis to update the m-dual.
    */
   void incDimDualBasis();

   /**
    * This method overrides its namesake in `IntLattice`. A basis for the projection of this
    * `MRGLattice` over the coordinates in `proj` is returned in `projLattice`.
    * The implementation used here exploits the rank-k lattice structure and it
    * is faster than the general one. See Section 5.5 of the LatMRG guide.
    */
   void buildProjection(IntLattice<Int, Real> &projLattice, const Coordinates &proj) override;

   /**
    * Overrides the same function from `IntLattice`.
    * When all the first k coordinates {1,...,k} belong to the projection, both the primal
    * and m-dual constructions are direct, just by selecting the rows and columns
    * whose indices are in `proj`. Otherwise an upper triangular basis is constructed
    * for the basis and a lower-triangular basis for the m-dual.
    */
   void buildProjectionDual(IntLattice<Int, Real> &projLattice, const Coordinates &proj) override;

   /**
    * Returns the first `dim` components of the generating vector \f$\ba\f$ as a string,
    * where `dim` is the current lattice dimension.
    */
   std::string toStringCoef() const;

   // Order of this MRG.
   int m_order;

   /**
    * The coefficients used for the MRG recurrence.
    */
   IntVec m_aCoeff;

protected:

   /**
    * These variants take the basis as a parameter for more flexibility.
    * They are used inside buildBais, buildBasisDual, incDimBasis, etc., with either `m_basis` or `m_basis0`.
    */
   void buildBasis0(IntMat &basis, int64_t d);

   void incDimBasis0(IntMat &basis, int64_t d);

   void buildProjection0(IntMat &basis, IntMat &pbasis, const Coordinates &proj);

   /**
    * This matrix contains the initial basis V for the MRG lattice as defined in
    * Section 4.1 of the LatMRG guide. It is used only to calculate or update the m-dual
    * basis. Its dimension is enlarged whenever needed and never reduced.
    */
   IntMat m_basis0;

   /**
    * The current dimension of `m_basis0`.
    */
   int64_t m_dim0 = 0;

   /**
    * This auxiliary matrix is used to store generating vectors of projections.
    */
   IntMat m_genTemp;

   /**
    * Hidden variable used across the `buildProjection` functions, to avoid recomputing it.
    * This variable is `true` if the first m_order coordinates are all in `proj`.
    */
   bool projCase1 = true; // Tif the first m_order coordinates are all in current `proj`.

};

//============================================================================

template<typename Int, typename Real>
MRGLattice<Int, Real>::MRGLattice(const Int &m, const IntVec &aa, int64_t maxDim, NormType norm) :
      IntLatticeExt<Int, Real>(m, maxDim, norm) {
   this->m_maxDim = maxDim;
   setaa(aa);
   m_order = aa.length();
   m_dim0 = 0;
   m_basis0.resize(maxDim, maxDim);
   m_genTemp.resize(maxDim, maxDim);
}

//============================================================================

template<typename Int, typename Real>
MRGLattice<Int, Real>::~MRGLattice() {
   m_aCoeff.kill();
   // m_aux_basis.kill();
}

//============================================================================
template<typename Int, typename Real>
MRGLattice<Int, Real>& MRGLattice<Int, Real>::operator=(const MRGLattice<Int, Real> &lat) {
   if (this == &lat) return *this;
   this->copy(lat);
   m_aCoeff = lat.m_aCoeff;
   m_order = lat.m_order;
   m_basis0 = lat.m_basis0;
   return *this;
}

//============================================================================

template<typename Int, typename Real>
MRGLattice<Int, Real>::MRGLattice(const MRGLattice<Int, Real> &lat) :
      IntLatticeExt<Int, Real>(lat.m_modulo, lat.getDim(),
            lat.getNormType()) {
   m_aCoeff = lat.m_aCoeff;
   m_order = lat.m_order;
   m_basis0 = lat.m_basis0;
   // Should also copy the basis and all other variables?
}

//============================================================================

/**
 * Sets the generating vector to `aa`.
 */
template<typename Int, typename Real>
void MRGLattice<Int, Real>::setaa(const IntVec &aa) {
   m_aCoeff = aa;
   m_order = aa.length();
   this->m_dim = 0;  // Current basis is now invalid.
   this->m_dimdual = 0;
   this->m_dim0 = 0;
}

//============================================================================

// Builds an upper-triangular basis directly in `d` dimensions, as explained in Section 4.1 of
// the guide of LatMRG, puts it in `basis`.  Must have d < m_maxdim.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildBasis0(IntMat &basis, int64_t d) {
   assert(d <= this->m_maxDim);
   int64_t dk = min(d, m_order);
   int64_t i, j, k;
   // Put the identity matrix in the upper left corner
   for (i = 0; i < dk; i++) {
      for (j = 0; j < dk; j++)
         basis[i][j] = (i == j);  // Avoid "if" statements.
   }
   if (d > m_order) {
      // Put m times the identity matrix to the lower right part and the zero matrix to the lower left part.
      for (i = m_order; i < d; i++)
         for (j = 0; j < d; j++)
            basis[i][j] = this->m_modulo * (i == j);
      // Fill the rest of the first m_order rows
      for (i = 0; i < m_order; i++) {
         for (j = m_order; j < d; j++) {
            basis[i][j] = 0;
            // Calculate the components of v_{i,j}. The first component of the coefficient is m_aCoeff[0] here
            for (k = 0; k < m_order; k++)
               basis[i][j] += m_aCoeff[k] * basis[i][j - (k + 1)] % this->m_modulo;
         }
      }
   }
}

//============================================================================

// An upper-triangular basis is built directly, as explained in the guide of LatMRG.
// The dimension `maxDim` of the `IntMat` array is unchanged, but the basis dimension isset to `d`.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildBasis(int64_t d) {
   this->setDim(d);
   this->buildBasis0(this->m_basis, d);
}

//============================================================================

// Builds the m-dual basis in a direct way in d dimensions. This uses m_basis0.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildDualBasis(int64_t d) {
   if (d > m_dim0) {  // Compute the initial primal basis up to d dimensions.
      this->buildBasis0(this->m_basis0, d);
      m_dim0 = d;
   }
   this->setDimDual(d);
   int64_t dk = min(d, m_order);
   int64_t i, j;
   for (i = 0; i < dk; i++)
      for (j = 0; j < d; j++)
         this->m_dualbasis[i][j] = (i == j) * this->m_modulo;
   for (i = dk; i < d; i++) {
      for (j = 0; j < dk; j++)
         this->m_dualbasis[i][j] = -this->m_basis[j][i];
      for (j = dk; j < d; j++)
         this->m_dualbasis[i][j] = (i == j);
   }
}

//============================================================================

// Increases the dimension of given basis from d-1 to d dimensions.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::incDimBasis0(IntMat &basis, int64_t d) {
   // int64_t d = 1 + this->getDim();  // New current dimension.
   assert(d <= this->m_maxDim);
   int64_t i, j, k;
   // Add new row and new column of the primal basis.
   for (j = 0; j < d - 1; j++)
      basis[d - 1][j] = 0;
   if (d - 1 >= m_order) basis[d - 1][d - 1] = this->m_modulo;
   else basis[d - 1][d - 1] = 1;
   for (i = 0; i < d - 1; i++) {
      basis[i][d - 1] = 0;
      if (d - 1 >= m_order) {
         for (k = 0; k < m_order; k++)
            basis[i][d - 1] += m_aCoeff[k] * basis[i][d - 1 - (k + 1)] % this->m_modulo;
      }
   }
}

//============================================================================

template<typename Int, typename Real>
void MRGLattice<Int, Real>::incDimBasis() {
   int64_t d = 1 + this->getDim();  // New current dimension.
   this->setDim(d);
   this->incDimBasis0(this->m_basis, d);
}

//============================================================================

template<typename Int, typename Real>
void MRGLattice<Int, Real>::incDimDualBasis() {
   int64_t d = 1 + this->getDimDual();  // New current dimension.
   this->setDimDual(d);
   while (m_dim0 < d) {  // Increase m_basis0 if needed.
      m_dim0++;
      this->incDimBasis0(m_basis0, m_dim0);
   }
   int64_t i;
   // Add one extra coordinate to each vector.
   for (i = 0; i < d - 1; i++) {
      this->m_dualbasis[i][d - 1] = 0;
      this->m_dualbasis[d - 1][i] = 0;
   }
   // Add the new vector the dual basis
   for (i = 0; i < m_order; i++)
      this->m_dualbasis[d - 1][i] = -this->m_basis0[i][d - 1];
   if (d - 1 < m_order) this->m_dualbasis[d - 1][d - 1] = this->m_modulo;
   else this->m_dualbasis[d - 1][d - 1] = 1;
}

//============================================================================

// We must be able to do this with either m_basis or m_basis0.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildProjection0(IntMat &basis, IntMat &pbasis,
      const Coordinates &proj) {
   int64_t d = proj.size();
   // assert (proj.end() <= unsigned(this->m_dim));
   // IntMat &pbasis = projLattice.getBasis(); // We want to construct this basis of the projection.
   int64_t i, j;
   // Check if we are in case 1.
   //  *** This seems to assume that the coordinates of the projection are always in increasing order, right?   ******
   projCase1 = true; // This holds if the first m_order coordinates are all in `proj`.
   if (d < (unsigned) m_order) projCase1 = false;
   else {
      j = 0;
      for (auto it = proj.begin(); it != proj.end(); it++, j++) {
         if (j < m_order) {
            if (*it != unsigned(j + 1)) projCase1 = false;
         } else break;
      }
   }
   if (projCase1) { // We first compute the first m_order rows of the projection basis.
      for (i = 0; i < m_order; i++) {
         j = 0;
         for (auto it = proj.begin(); it != proj.end(); it++, j++)
            pbasis[i][j] = this->m_basis[i][*it - 1];
      }
      // Then the other rows.
      for (i = m_order; i < d; i++)
         for (j = 0; j < d; j++)
            pbasis[i][i] = this->m_modulo * (i == j);
   } else {
      // In this case we need to use the generic algorithm
      j = 0;
      for (auto it = proj.begin(); it != proj.end(); it++, j++) {
         for (i = 0; i < unsigned(d); i++)
            m_genTemp[i][j] = this->m_basis[i][*it - 1];
      }
      upperTriangularBasis(m_genTemp, pbasis, this->m_modulo, d, d);
   }
}

//============================================================================

// This cannot use m_basis0, only m_basis, because basis0 may not exist.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildProjection(IntLattice<Int, Real> &projLattice,
      const Coordinates &proj) {
   // int64_t d = proj.size();
   // assert (proj.end() <= unsigned(this->m_dim));
   // IntMat &pbasis = projLattice.getBasis(); // We want to construct this basis of the projection.
   this->buildProjection0 (this->m_basis, projLattice.getBasis(), proj);
   projLattice.setDim(proj.size());
}

//============================================================================
template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildProjectionDual(IntLattice<Int, Real> &projLattice,
      const Coordinates &proj) {
   int64_t d = proj.size();
   // assert (proj.end() <= unsigned(this->m_dim0));
   IntMat &pbasis = projLattice.getBasis();  // Reference to basis of projection.
   IntMat &dualBasis = projLattice.getDualBasis();   // And its m-dual.

   // We first build a basis for the primal projection.
   this->buildProjection0 (this->m_basis0, pbasis, proj);
   projLattice.setDim(d);

   // Then we compute its m-dual.
/*
   int64_t i, j, k;
   if (projCase1) { // Compute the m-dual basis directly.    ***  Is it really faster and worthwhile?   ****
      for (j = 0; j < m_order; j++) {
         dualBasis[j][j] = this->m_modulo;
         i = m_order;
         auto it = proj.begin();
         for (k = 0; k < m_order - 1; k++)
            it++;
         for (it++; it != proj.end(); ++it, ++i) {
            dualBasis[i][j] = -this->m_basis0[j][*it - 1]; // First column.
         }
      }
      for (i = 0; i < d; i++)
         for (j = m_order; j < d; j++)
            dualBasis[i][j] = (i == j);  // ???
   } else {
*/
      mDualUpperTriangular(pbasis, dualBasis, this->m_modulo, d);
      projLattice.setDimDual(d);
//   }
}

//============================================================================
template<typename Int, typename Real>
std::string MRGLattice<Int, Real>::toStringCoef() const {
   return toString(m_aCoeff, 0, this->getDim());
}

//============================================================================

template class MRGLattice<std::int64_t, double> ;
template class MRGLattice<NTL::ZZ, double> ;
template class MRGLattice<NTL::ZZ, xdouble> ;
template class MRGLattice<NTL::ZZ, quad_float> ;
template class MRGLattice<NTL::ZZ, NTL::RR> ;

} // End namespace LatticeTester

#endif
