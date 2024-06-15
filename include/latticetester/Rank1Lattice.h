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

#ifndef LATTICETESTER_RANK1LATTICE_H
#define LATTICETESTER_RANK1LATTICE_H

#include "latticetester/Util.h"
#include "latticetester/EnumTypes.h"
#include "latticetester/IntLatticeExt.h"

namespace LatticeTester {

/**
 * This subclass of `IntLatticeExt` defines a general rank 1 lattice rule in \f$d\f$ dimensions,
 * whose points \f$\mathbb{u}_i\f$ are defined by
 * \f{equation}{
 *    \mathbf{u}_i = (i \mathbf{a} \mod m)/m,
 * \f}
 * where \f$\mathbf{a} = (a_1,a_2,\dots,a_d) \in \mathbb{Z}_m^d\f$ is called the generating vector.
 * The lattice is rescaled simply by removing the division by \f$m\f$.
 * Given \f$\mathbf{a}\f$, a basis of the rescaled (integer) lattice is given by
 * \f{align*}{
 *   \mathbf{v}_1 & = \mathbf{a} \\
 *   \mathbf{v}_2 & = m \mathbf{e}_2 \\
 *   \vdots & \\
 *   \mathbf{v}_d & = m \mathbf{e}_d�\\
 * \f}
 * where \f$\mathbf{e}_i\f$ is the \f$i^\text{th}\f$ unit vector.
 *
 * A condition often required when building a rank 1 lattice is that
 * \f$\gcd(a_j, m) = 1\f$ for \f$j =1, \dots,d\f$.  When this condition is
 * verified, each lower-dimensional projection of the lattice contains the same number
 * of points (has the same density \f$m\f$) as the full lattice.
 * When searching for lattices that satisfy this condition, one may assume
 * without loss of generality generality that \f$a_1 = 1\f$.
 * Under that condition, it is straightforward to construct a basis for a projection that
 * contains the first coordinate, and also its m-dual basis.
 * In this class, we assume that \f$a_1 = 1\f$ and \f$\gcd(a_j, m) = 1\f$ for all \f$j\f$.
 * We exploit these conditions when building bases for the lattice, the projections,
 * and their m-duals.
 *
 * The functions `buildBasis` and `incDimBasis` always build and update the primal
 * basis, because it is required to update the m-dual.  To build only the m-dual,
 * use `buildDualBasis`.
 */

template<typename Int, typename Real>
class Rank1Lattice: public IntLatticeExt<Int, Real> {

private:
   typedef NTL::vector<Int> IntVec;
   typedef NTL::matrix<Int> IntMat;
   typedef NTL::vector<Real> RealVec;

public:

   /**
    * This constructor takes as input the modulus `m`, the generating vector `aa`,
    * and the norm used to measure the vector lengths.
    * The maximal dimension `maxDim` will be the length `d` of the vector `aa`.
    * The coefficient @f$a_j@f$ will be `aa[j-1]`. One must have `aa[0] = 1`.
    * The variable `withPrimal` indicates if the primal basis will be maintained or not,
    * and `withDual` indicates if the dual basis will be maintained or not.
    * This constructor does not build the basis, to leave
    * more flexibility in the dimension when doing so.
    */
   Rank1Lattice(const Int &m, const IntVec &aa, NormType norm = L2NORM);

   /**
    * Constructor for the special case of a Korobov lattice.
    * Here the generating vector has the form @f$\ba = (1, a, a^2 mod m, a^3 mod m, ...)@f$
    * where @f$a@f$ is an integer such that @f$1 < a < m@f$.
    */
   Rank1Lattice(const Int &m, const Int &a, int64_t maxDim, NormType norm = L2NORM);

   /**
    * This constructor does not specify the generating vector.
    */
   Rank1Lattice(const Int &m, int64_t maxDim, NormType norm = L2NORM);

   /**
    * Copy constructor.
    */
   Rank1Lattice(const Rank1Lattice<Int, Real> &Lat);

   /**
    * Assigns `Lat` to this object.
    */
   Rank1Lattice& operator=(const Rank1Lattice<Int, Real> &Lat);

   /**
    * Destructor.
    */
   ~Rank1Lattice();

   /**
    * Sets the generating vector to `aa`.  The maximum dimension will be set to its dimension.
    */
   void setaa(const IntVec &aa);

   /**
    * Sets this lattice to a Korobov lattice with multiplier `a`.
    * The maximum dimension is unchanged.
    */
   void seta(const Int &a);

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
    * This method uses the simplified method given in the lattice tester guide:
    * the new m-dual basis vector is simply  @f$\bw_d = (-a_d, 0, ..., 0, 1)@f$.
    */
   void incDimDualBasis();

   /**
    * This method overrides its namesake in `IntLattice`. The projection of this
    * `Rank1Lattice` over the coordinates in `proj` is returned in `projLattice`.
    * The implementation used here exploits the rank-1 lattice structure and it
    * is simpler and faster than the general one. See Section 5.5 of the guide.
    * When the first coordinate is 1 and belongs to the projection, both the primal
    * and m-dual constructions are direct, just by selecting the rows and columns
    * whose indices are in `proj`. Otherwise, those rows and columns plus the
    * first row form a set of `proj.size()+1` generating vectors for the primal.
    * This number must not exceed the `maxDim` of `projLattice`.
    */
   void buildProjection(IntLattice<Int, Real> &projLattice, const Coordinates &proj) override;

   /**
    * Similar to `buildProjection`, but builds a basis for the m-dual of the projection.
    */
   void buildProjectionDual(IntLattice<Int, Real> &projLattice, const Coordinates &proj) override;

   /**
    * Returns the first `dim` components of the generating vector \f$\ba\f$ as a string,
    * where `dim` is the current lattice dimension.
    */
   std::string toStringCoef() const;

protected:

   /**
    * Vector of multipliers (generating vector) of the rank 1 lattice rule.
    * They are stored for up to `maxDim()` dimensions.
    * The first coordinate has index 0.
    */
   IntVec m_a;
};

//============================================================================

template<typename Int, typename Real>
Rank1Lattice<Int, Real>::Rank1Lattice(const Int &m, int64_t maxDim, NormType norm) :
      IntLatticeExt<Int, Real>(m, maxDim, norm) {
   this->m_a.SetLength(maxDim);
}

//============================================================================

template<typename Int, typename Real>
Rank1Lattice<Int, Real>::Rank1Lattice(const Int &m, const IntVec &aa, NormType norm) :
      IntLatticeExt<Int, Real>(m, aa.length(), norm) {
   this->m_maxDim = aa.length();
   this->setaa(aa);
}

//============================================================================

template<typename Int, typename Real>
Rank1Lattice<Int, Real>::Rank1Lattice(const Int &m, const Int &a, int64_t maxDim, NormType norm) :
      IntLatticeExt<Int, Real>(m, maxDim, norm) {
   this->m_a.SetLength(maxDim);
   this->seta(a);    // Set to a Korobov lattice.
}

//============================================================================

template<typename Int, typename Real>
Rank1Lattice<Int, Real>::~Rank1Lattice() {
   this->m_a.kill();
}

//============================================================================
template<typename Int, typename Real>
Rank1Lattice<Int, Real>& Rank1Lattice<Int, Real>::operator=(const Rank1Lattice<Int, Real> &lat) {
   if (this == &lat) return *this;
   this->copy(lat);
   this->m_a = lat.m_a;
   return *this;
}

//============================================================================

template<typename Int, typename Real>
Rank1Lattice<Int, Real>::Rank1Lattice(const Rank1Lattice<Int, Real> &lat) :
      IntLatticeExt<Int, Real>(lat.m_modulo, lat.m_maxDim, lat.getNormType()) {
   this->m_a = lat.m_a;
   // Should also copy the basis and all other variables?
}

//============================================================================

/**
 * Sets the generating vector to `aa`.
 */
template<typename Int, typename Real>
void Rank1Lattice<Int, Real>::setaa(const IntVec &aa) {
   if (aa[0] != 1) MyExit(1, "Rank1Lattice::setaa: must have aa[0] == 1.");
   this->m_a = aa;
}

//============================================================================

/**
 * Sets this lattice to a Korobov lattice with multiplier `a`, with dim = maxDim.
 */
template<typename Int, typename Real>
void Rank1Lattice<Int, Real>::seta(const Int &a) {
   this->m_a[0] = 1;
   for (int64_t i = 1; i < this->m_maxDim; i++) {
      this->m_a[i] = (a * this->m_a[i - 1]) % this->m_modulo;
   }
}

//============================================================================

// An upper-triangular basis is built directly, as explained in the guide of Lattice Tester.
// The dimension `maxDim` of the `IntMat` array is unchanged.
template<typename Int, typename Real>
void Rank1Lattice<Int, Real>::buildBasis(int64_t d) {
   // assert(d <= this->m_maxDim);
   this->setDim(d);
   int64_t i, j;
   // This builds an upper-triangular primal basis in a standard direct way.
   for (j = 0; j < d; j++)
      this->m_basis[0][j] = this->m_a[j];
   for (i = 1; i < d; i++)
      for (j = 0; j < d; j++)
         this->m_basis[i][j] = this->m_modulo * (i == j);
   this->setNegativeNorm();
}

//============================================================================

// This builds the m-dual basis, also in a direct way.
template<typename Int, typename Real>
void Rank1Lattice<Int, Real>::buildDualBasis(int64_t d) {
   // assert(d <= this->m_maxDim);
   this->setDimDual(d);
   int64_t i, j;
   this->m_dualbasis[0][0] = this->m_modulo;
   for (j = 1; j < d; j++)
      this->m_dualbasis[0][j] = 0;
   for (i = 1; i < d; i++) {
      this->m_dualbasis[i][0] = -this->m_a[i];
      for (j = 1; j < d; j++)
         this->m_dualbasis[i][j] = (i == j);
   }
   this->setDualNegativeNorm();
}

//============================================================================

template<typename Int, typename Real>
void Rank1Lattice<Int, Real>::incDimBasis() {
   int64_t d = 1 + this->getDim();  // New current dimension.
   this->setDim(d);
   int64_t i, j;
   // Add new row and new column of the primal basis.
   for (j = 0; j < d - 1; j++)
      this->m_basis[d - 1][j] = 0;
   this->m_basis[d - 1][d - 1] = this->m_modulo;
   for (i = 0; i < d - 1; i++) {
      this->m_basis[i][d - 1] = (this->m_a[d - 1] * this->m_basis[i][0]) % this->m_modulo;
   }
   // this->setNegativeNorm();
}

//============================================================================

template<typename Int, typename Real>
void Rank1Lattice<Int, Real>::incDimDualBasis() {
   int64_t d = 1 + this->getDim();
   this->setDimDual(d);
   // Add one extra coordinate to each vector.
   for (int64_t i = 0; i < d; i++) {
      this->m_dualbasis[i][d - 1] = 0;
      this->m_dualbasis[d - 1][i] = 0;
   }
   this->m_dualbasis[d - 1][0] = -m_a[d - 1];
   this->m_dualbasis[d - 1][d - 1] = 1;
   // this->setDualNegativeNorm();
}

//============================================================================
template<typename Int, typename Real>
void Rank1Lattice<Int, Real>::buildProjection(IntLattice<Int, Real> &projLattice,
      const Coordinates &proj) {
   // Builds only the primal basis for the projection.
   // We use the method described in the Lattice Tester guide, section 5.5.
   // It *does not* assume that a basis for `this` has been computed before.
   long i, j;
   long d = proj.size();     // Number of coordinates in the projection.
   projLattice.setDim(d);
   IntMat &basis = projLattice.getBasis();  // Reference to the projection basis.
   bool case1 = m_a[*proj.begin() - 1] == 1; // We have a_{i_1} = 1.
   // We first compute the first row of the basis matrix.
   if (case1) {  // The easy (usual) case.
      j = 0;
      for (auto it = proj.begin(); it != proj.end(); it++, j++)
         basis[0][j] = m_a[*it - 1];
   } else {
      Int c1, b1, b2;   // c1 will be gcd(a_{i_1}, m). We should always have c1=1.
      // Recall:  XGCD (g, c, d, a, b) does g = gcd(a, b) = a*c + b*d.
      NTL::XGCD(c1, b1, b2, m_a[*proj.begin() - 1], this->m_modulo);
      if (c1 > 1) MyExit(1, "Rank1Lattice::buildProjection: c1 > 1.");
      basis[0][0] = 1;
      auto it = proj.begin();
      j = 1;
      for (it++; it != proj.end(); it++, j++)
         basis[0][j] = (m_a[*it - 1] * b1) % this->m_modulo;
   }
   // Then the other rows.
   for (i = 1; i < d; i++)
      for (j = 0; j < d; j++)
         basis[i][j] = this->m_modulo * (i == j);
}

//============================================================================
template<typename Int, typename Real>
void Rank1Lattice<Int, Real>::buildProjectionDual(IntLattice<Int, Real> &projLattice,
      const Coordinates &proj) {
   // The dual basis is computed directly, without using the primal.
   // We use the method described in the Lattice Tester guide, section 5.5.
   // It *does not* assume that a basis for `this` has been computed before.
   long i, j;
   long d = proj.size();     // Number of coordinates in the projection.
   projLattice.setDim(d);
   projLattice.setDimDual(d);
   IntMat &dualBasis = projLattice.getDualBasis();
   dualBasis[0][0] = this->m_modulo;
   i = 1;
   auto it = proj.begin();
   bool case1 = m_a[*proj.begin() - 1] == 1; // We have a_{i_1} = 1.
   // We first compute the first column of the m-dual basis.
   if (case1) {
      for (it++; it != proj.end(); ++it, ++i)
         dualBasis[i][0] = m_a[*it - 1]; // First column.
   } else {
      Int c1, b1, b2; // c1 will be gcd(a_{i_1}, m). We should always have c1=1.
      // Recall:  XGCD (g, c, d, a, b) does g = gcd(a, b) = a*c + b*d.
      NTL::XGCD(c1, b1, b2, m_a[*proj.begin() - 1], this->m_modulo);
      if (c1 > 1) MyExit(1, "Rank1Lattice::buildProjection: c1 > 1.");
      for (it++; it != proj.end(); it++, i++)
         dualBasis[i][0] = -(m_a[*it - 1] * b1 % this->m_modulo);
   }
   // Then we initialize all the other columns.
   for (i = 0; i < d; i++)
      for (j = 1; j < d; j++)
         dualBasis[i][j] = (i == j);
}

//============================================================================
template<typename Int, typename Real>
std::string Rank1Lattice<Int, Real>::toStringCoef() const {
   return toString(this->m_a, 0, this->getDim());
}

//============================================================================

template class Rank1Lattice<std::int64_t, double> ;
template class Rank1Lattice<NTL::ZZ, double> ;
template class Rank1Lattice<NTL::ZZ, xdouble> ;
template class Rank1Lattice<NTL::ZZ, quad_float> ;
template class Rank1Lattice<NTL::ZZ, NTL::RR> ;

} // End namespace LatticeTester

#endif
