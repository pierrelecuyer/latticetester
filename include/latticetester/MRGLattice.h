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

#include "latticetester/Util.h"
#include "latticetester/EnumTypes.h"
#include "latticetester/IntLatticeExt.h"
#include "latticetester/BasisConstruction.h"

namespace LatticeTester {

/**
 * This subclass of `IntLatticeExt` defines an MRG lattce and is similar to Rank1Lattice,
 * but constructs lattices associated with multiple recursive generators (MRGs) with 
 * modulus m, order k, and vector of multipliers a = (a1, . . . , ak).
 *
 * The functions `buildBasis` and `incDimBasis` always build and update the primal
 * basis, because it is required to update the m-dual. To build only the m-dual,
 * use `buildDualBasis`.
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
    * The variable `withPrimal` indicates if the primal basis will be maintained or not,
    * and `withDual` indicates if the dual basis will be maintained or not.
    * This constructor does not build the basis, to leave
    * more flexibility in the dimension when doing so.
    */
   MRGLattice(const Int &m, const IntVec &aa, int64_t maxDim, bool withPrimal =
         false, bool withDual = false, NormType norm = L2NORM);

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
    * This function can be called only when `withPrimal` is set to true.
    * If `withDual` is true, it also builds an m-dual lower-triangular basis.
    */
   void buildBasis(int64_t dim);

   /**
    * Builds only an m-dual lower triangular basis (and not the primal) directly
    * in `dim` dimensions.  This `dim` must not exceed `maxDim`.
    */
   void buildDualBasis(int64_t dim);

   /**
    * Increases the current dimension of the primal basis by 1 and updates the basis.
    * This function can be called only when `withPrimal` is set to true.
    * If `withDual`, it also increases the m-dual basis and makes it the m-dual of the primal basis.
    * The new increased dimension must not exceed `maxDim`.
    */
   void incDimBasis();

   /**
    * Increases the current dimension of only the m-dual basis by 1.
    * The primal basis is left unchanged (not updated).
    * The new increased dimension must not exceed `maxDim`.
    * This method uses the simplified method for MRG lattices given in the lattice tester guide.
    */
   void incDimDualBasis();

   /**
    * This method overrides its namesake in `IntLattice`. The projection of this
    * `MRGLattice` over the coordinates in `proj` is returned in `projLattice`.
    * The implementation used here exploits the rank-k lattice structure and it
    * is faster than the general one. See Section 5.5 of the guide.
    */
   void buildProjection(IntLattice<Int, Real> &projLattice,
         const Coordinates &proj) override;

   /**
    * Overrides the same function from `IntLattice`.
    * When all the first k coordinates {1,...,k} belong to the projection, both the primal
    * and m-dual constructions are direct, just by selecting the rows and columns
    * whose indices are in `proj`. Otherwise an upper triangular basis is constructed.
    */
   void buildProjectionDual(IntLattice<Int, Real> &projLattice,
         const Coordinates &proj) override;

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
    * Vector of multipliers (generating vector) of the rank 1 lattice rule.
    * They are stored for up to `maxDim()` dimensions.
    * The first coordinate has index 0.
    */
   // IntVec m_a;   // Not used!
   /**
    * This auxillary matrix contains the standard basis for the chosen maxDim.
    * It is only used to calculate the dual matrix in case that the primal
    * basis is not maintained.
    */
   // IntMat m_aux_basis;

   /**
    * This auxillary matrix is used to store generating vectors of projections.
    */
   IntMat m_genTemp;

   /**
    * Hidden variable used across the `buildProjection` functions, to avoid recomputing it.
    * This variable is `true` if the first m_order coordinates are all in `proj`.
    */
   bool projCase1 = true; // This holds if the first m_order coordinates are all in `proj`.

};

//============================================================================

template<typename Int, typename Real>
MRGLattice<Int, Real>::MRGLattice(const Int &m, const IntVec &aa,
      int64_t maxDim, bool withPrimal, bool withDual, NormType norm) :
      IntLatticeExt<Int, Real>(m, maxDim, withPrimal, withDual, norm) {
   this->m_maxDim = maxDim;
   setaa(aa);
   m_order = aa.length();
   //m_aux_basis.resize(maxDim, maxDim);
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
MRGLattice<Int, Real>& MRGLattice<Int, Real>::operator=(
      const MRGLattice<Int, Real> &lat) {
   if (this == &lat)
      return *this;
   this->copy(lat);
   m_aCoeff = lat.m_aCoeff;
   m_order = lat.m_order;
   return *this;
}

//============================================================================

template<typename Int, typename Real>
MRGLattice<Int, Real>::MRGLattice(const MRGLattice<Int, Real> &lat) :
      IntLatticeExt<Int, Real>(lat.m_modulo, lat.getDim(), lat.m_withPrimal,
            lat.m_withDual, lat.getNormType()) {
   m_aCoeff = lat.m_aCoeff;
   m_order = lat.m_order;
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
}

//============================================================================

// An upper-triangular basis is built directly, as explained in the guide of LatMRG.
// The dimension `maxDim` of the `IntMat` array is unchanged, but the basis dimension isset to `d`.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildBasis(int64_t d) {
   this->m_dim = d;
   int64_t dk = min(d, m_order);
   int64_t i, j, k;
   // Put the identity matrix in the upper left corner
   for (i = 0; i < dk; i++) {
      for (j = 0; j < dk; j++)
         this->m_basis[i][j] = (i == j);  // Avoid "if" statements.
   }
   if (d > m_order) {
      // Put m times the identity matrix to the lower right part and the zero matrix to the lower left part.
      for (i = m_order; i < d; i++)
         for (j = 0; j < d; j++)
            this->m_basis[i][j] = this->m_modulo * (i == j);
      // Fill the rest of the first m_order rows
      for (i = 0; i < m_order; i++) {
         for (j = m_order; j < d; j++) {
            this->m_basis[i][j] = 0;
            // Calculate the components of v_{i,j}. The first component of the coefficient is m_aCoeff[0] here
            for (k = 0; k < m_order; k++)
               this->m_basis[i][j] += m_aCoeff[k]
                     * this->m_basis[i][j - (k + 1)] % this->m_modulo;
         }
      }
   }
   this->setNegativeNorm();
}

//============================================================================

// This one builds the m-dual basis in a direct way, in d dimensions.
template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildDualBasis(int64_t d) {
   this->m_dim = d; //       ****  Attention: same dim variable for both primal and dual...  ****
   int64_t dk = min(d, m_order);
   int64_t i, j, k;

   for (i = 0; i < dk; i++) {
      for (j = 0; j < dk; j++) {
         this->m_dualbasis[i][j] = (i == j);
      }
   }
   if (d > m_order) {
      // Fill the upper right with zeros
      for (i = 0; i < m_order; i++) {
         for (j = m_order; j < d; j++) {
            this->m_dualbasis[i][j] = 0;
         }
      }
      for (i = m_order; i < d; i++) {
         for (j = 0; j < d; j++) {
            this->m_dualbasis[i][j] = (i == j);
            if (j < m_order)
               // Calculate the components of w_{i,j} = -v_{j,i}.
               // The first component of the coefficient is m_aCoeff[0] here
               for (k = 0; k < m_order; k++) {
                  this->m_dualbasis[i][j] -= m_aCoeff[k]
                        * this->m_dualbasis[i - (k + 1)][j] % this->m_modulo;
               }
         }
      }
      // Change these values to their negatives.   **** We now compute the negative directly.  *****
      //for (i = m_order; i < d; i++) {
      //   for (j = 0; j < m_order; j++)
      //      this->m_dualbasis[i][j] = -this->m_dualbasis[i][j];
      //}
   }
   for (i = 0; i < dk; i++) {
      this->m_dualbasis[i][i] = this->m_modulo;
   }
   this->setDualNegativeNorm();
}

//============================================================================

template<typename Int, typename Real>
void MRGLattice<Int, Real>::incDimBasis() {
   // assert(this->m_withPrimal);      // m_withPrimal must be true.
   int64_t d = 1 + this->getDim();  // New current dimension.
   assert(d <= this->m_maxDim);
   this->m_dim = d;
   int64_t i, j, k;
   // Add new row and new column of the primal basis.
   for (j = 0; j < d - 1; j++)
      this->m_basis[d - 1][j] = 0;
   if (d - 1 >= m_order) {
      this->m_basis[d - 1][d - 1] = this->m_modulo;
   } else
      this->m_basis[d - 1][d - 1] = 1;
   for (i = 0; i < d - 1; i++) {
      this->m_basis[i][d - 1] = 0;
      if (d - 1 >= m_order) {
         for (k = 0; k < m_order; k++)
            this->m_basis[i][d - 1] += m_aCoeff[k]
                  * this->m_basis[i][d - 1 - (k + 1)] % this->m_modulo;
      }
   }
   this->setNegativeNorm();
}

//============================================================================

template<typename Int, typename Real>
void MRGLattice<Int, Real>::incDimDualBasis() {
   this->incDimBasis();  // We first increase the primal basis by one dimension.
   int64_t d = this->m_dim;
   // assert(d <= m_maxDim);
   // m_dim = d;
   int64_t i;
   // Add one extra coordinate to each vector.
   for (i = 0; i < d - 1; i++) {
      this->m_dualbasis[i][d - 1] = 0;
      this->m_dualbasis[d - 1][i] = 0;
   }
   // Add the new vector the dual basis
   for (i = 0; i < m_order; i++)
      this->m_dualbasis[d - 1][i] = -this->m_basis[i][d - 1];
   if (d - 1 < m_order) {
      this->m_dualbasis[d - 1][d - 1] = this->m_modulo;
   } else
      this->m_dualbasis[d - 1][d - 1] = 1;
}

//============================================================================
template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildProjection(IntLattice<Int, Real> &projLattice,
      const Coordinates &proj) {
   int64_t d = proj.size();
   IntMat &basis = projLattice.getBasis(); // We want to construct this basis of the projection.
   int64_t i, j;
   // Check if we are in case 1.
   //  *** Do we assume that the coordinates of the projection are always in increasing order?   ******
   projCase1 = true; // This holds if the first m_order coordinates are all in `proj`.
   if (d < (unsigned) m_order)
      projCase1 = false;
   else {
      j = 0;
      for (auto it = proj.begin(); it != proj.end(); it++, j++) {
         if (j < m_order) {
            if (*it != unsigned(j + 1))
               projCase1 = false;
         } else
            break;
      }
   }
   if (projCase1) { // We first compute the first m_order rows of the projection basis.
      for (i = 0; i < m_order; i++) {
         j = 0;
         for (auto it = proj.begin(); it != proj.end(); it++, j++)
            basis[i][j] = this->m_basis[i][*it - 1];
      }
      // Then the other rows.
      for (i = m_order; i < d; i++)
         for (j = 0; j < d; j++)
            basis[i][i] = this->m_modulo * (i == j);
   } else {
      // In this case we need to use the generic algorithm
      j = 0;
      for (auto it = proj.begin(); it != proj.end(); it++, j++) {
         for (i = 0; i < unsigned(d); i++)
            m_genTemp[i][j] = this->m_basis[i][*it - 1]; // ****  Check if m_genTemp  is large enough.
      }
      upperTriangularBasis(m_genTemp, basis, this->m_modulo, d, d);
   }
}

//============================================================================
template<typename Int, typename Real>
void MRGLattice<Int, Real>::buildProjectionDual(
      IntLattice<Int, Real> &projLattice, const Coordinates &proj) {
   int64_t i, j, k;
   long d = proj.size();
   IntMat &basis = projLattice.getBasis();  // Reference to basis.
   IntMat &dualBasis = projLattice.getDualBasis();
   // We first build a basis for the projection.
   this->buildProjection (projLattice, proj);
   // Then we compute its m-dual.
   if (projCase1) { // Compute the m-dual basis directly.    ***  Is it really faster and worthwhile?   ****
      for (j = 0; j < m_order; j++) {
         dualBasis[j][j] = this->m_modulo;
         i = m_order;
         auto it = proj.begin();
         for (k = 0; k < m_order - 1; k++)
            it++;
         for (it++; it != proj.end(); ++it, ++i) {
            dualBasis[i][j] = -this->m_basis[j][*it - 1]; // First column.
         }
      }
      for (i = 0; i < d; i++)
         for (j = m_order; j < d; j++)
            dualBasis[i][j] = (i == j);
   } else {
      mDualUpperTriangular(basis, dualBasis, this->m_modulo, proj.size());
      projLattice.setDim(proj.size());
   }
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
