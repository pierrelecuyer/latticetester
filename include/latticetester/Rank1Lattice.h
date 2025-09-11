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
 * \class Rank1Lattice
 *
 * This subclass of `IntLatticeExt` defines a general rank-1 lattice rule in \f$t\f$ dimensions,
 * whose \f$m\f$ points are \f$ \mathbf{u}_i = (i \mathbf{a} \bmod m)/m \f$ for \f$i = 0,\dots,m-1\f$,
 * where \f$\mathbf{a} = (a_1,a_2,\dots,a_t) \in \mathbb{Z}_m^t\f$ is the generating vector,
 * \f$a_1 = 1\f$, and \f$\gcd(a_j, m) = 1\f$ for \f$j =2, \dots,t\f$.
 *
 * The lattice is rescaled simply by removing the division by \f$m\f$.
 * A simple upper-triangular basis for the rescaled lattice is given by
 * \f$\mathbf{v}_1 = \mathbf{a},\ \mathbf{v}_2 = m \mathbf{e}_2, \dots,
 * \mathbf{v}_t = m \mathbf{e}_t \f$,
 * where \f$\mathbf{e}_j\f$ is the \f$j^\text{th}\f$ unit vector.
 * The above conditions ensure that each lower-dimensional projection of
 * \f$\Lambda_t \cap [0,m)^t\f$ contains exactly \f$m\f$ distinct points.
 * Under these conditions, it is also straightforward to construct a basis for any projection
 * of \f$\Lambda_t\f$ and for its \f$m\f$-dual, as explained in Section 5.4 of the guide.
 * We exploit these properties.
 * The functions `buildBasis` and `incDimBasis` always build and update only the rescaled primal basis.
 * To build only an \f$m\f$-dual basis, use `buildDualBasis`.
 *
 * The dimension \f$t\f$ of the generating vector \f$\mathbf{a}\f$ may differ from the maximal dimension `maxDim`
 * of a lattice basis. It is specified separately either via the parameter `dimaa` or as the length
 * of the given vector ``aa`.
 */

template<typename Int, typename Real>
class Rank1Lattice: public IntLatticeExt<Int, Real> {

public:

   /**
    * This constructor takes as input the modulus `m`, the generating vector \f$\mathbf{a} = \f$`aa`,
    * and the norm used to measure the vector lengths.
    * The parameter `maxDim` gives the maximal dimension of a basis.
    * The constructor allocates space for that dimension, but it does not build the basis.
    * The current basis dimension is initially 0.
    * The maximal length of the generating vector will be the length `t` of the given vector `aa`.
    * The coefficient @f$a_j@f$ will be `aa[j-1]`. One must have `aa[0] = 1`.
    */
   Rank1Lattice(const Int &m, const IntVec &aa, int64_t maxDim, NormType norm = L2NORM);

   /**
   * In this version, the maximal dimension `maxDim` is set to the length of the vector `aa`.
   */
   Rank1Lattice(const Int &m, const IntVec &aa, NormType norm = L2NORM);

   /**
    * Constructor for the special case of a Korobov lattice.
    * Here the generating vector has the form @f$\mathbf{a} = (1, a, a^2 \bmod m, a^3 \bmod m, ...)@f$
    * where @f$a@f$ is an integer such that @f$1 < a < m@f$ and \f$\gcd(a, m) = 1\f$.
    * The maximal dimension of a basis will be `maxDim` and the maximal length of the
    * generating vector will be `dimaa`.
    * The basis is not initialized, so its current dimension is initially set to 0.
    */
   Rank1Lattice(const Int &m, const Int &a, int64_t maxDim, int64_t dimaa, NormType norm = L2NORM);

   /**
    * This version assumes that `dimaa = maxDim`.
    */
   Rank1Lattice(const Int &m, const Int &a, int64_t maxDim, NormType norm = L2NORM);

   /**
    * This constructor does not specify the generating vector but reserves space for it.
    */
   Rank1Lattice(const Int &m, int64_t maxDim, int64_t dimaa, NormType norm = L2NORM);

   /**
    * This version assumes that `dimaa = maxDim`.
    */
   Rank1Lattice(const Int &m, int64_t maxDim, NormType norm = L2NORM);

   /*
    * Copy constructor.
    */
   // Rank1Lattice(const Rank1Lattice<Int, Real> &Lat);

   /*
    * Assigns `Lat` to this object.
    */
   // Rank1Lattice& operator=(const Rank1Lattice<Int, Real> &Lat);

   /**
    * Destructor.
    */
   ~Rank1Lattice();

   /**
    * Sets the generating vector to `aa`. Its maximum dimension (value of `t`) will be set to the length of `aa`.
    */
   void setaa(const IntVec &aa);

   /**
    * Sets this lattice to a Korobov lattice with multiplier `a`
    * and sets the length of the generating vector to `dimaa`.
    */
   void seta(const Int &a, int64_t dimaa);

   /**
    * Sets this lattice to a Korobov lattice with multiplier `a`.
    * The maximum dimension of the generating vector is unchanged.
    */
   void seta(const Int &a);

   /**
    * Returns the generating vector `aa`.
    */
   IntVec getaa();

   /**
    * Builds an upper-triangular primal basis in `dim` dimensions. This `dim` must not exceed `this->maxDim()`.
    */
   void buildBasis(int64_t dim);

   /**
    * Builds a lower-triangular m-dual basis (and not the primal) directly
    * in `dim` dimensions.  This `dim` must not exceed `this->maxDim()`.
    */
   void buildDualBasis(int64_t dim);

   /**
    * Increases the current dimension of the primal basis by 1 and updates the basis.
    * The new increased dimension must not exceed `maxDim()`.
    */
   void incDimBasis();

   /**
    * Increases the current dimension of only the m-dual basis by 1.
    * The primal basis is left unchanged (not updated).
    * The new increased dimension must not exceed `maxDim()`.
    * This method uses the simplified method given in the lattice tester guide:
    * the new m-dual basis vector is simply  @f$\mathbf{w}_s = (-a_s, 0, ..., 0, 1)@f$.
    */
   void incDimDualBasis();

   /**
    * Builds a basis for the projection of this `Rank1Lattice` onto the coordinates
    * in `coordSet` and puts it as the `m_basis` in `projLattice`.
    * The construction is direct, just by selecting the rows and columns
    * whose indices are in `coordSet`, and using the vector `aa`. See Section 5.4 of the guide.
    * The implementation is simpler and faster than the general one.
    * The dimension `coordSet.size()` must not exceed the `maxDim` of `projLattice`.
    * Note that `projLattice` is allowed to be the same as the current `IntLattice` object,
    * in which case the projection basis will be built into the current `m_basis`.
    */
   void buildProjection(IntLattice<Int, Real> &projLattice, const Coordinates &coordSet) override;

   /**
    * Similar to `buildProjection`, but builds a basis for the m-dual of the projection.
    */
   void buildProjectionDual(IntLattice<Int, Real> &projLattice, const Coordinates &coordSet) override;

   /**
    * Returns the first `dim` components of the generating vector \f$\mathbf{a}\f$ as a string,
    * where `dim` is the current lattice basis dimension.
    */
   std::string toStringCoef() const;

protected:

   /**
    * Vector of multipliers (generating vector) of the rank 1 lattice rule.
    * They are stored for up to `m_dimaa` dimensions.
    * The first coordinate has index 0.
    */
   IntVec m_aa;

   /**
    * Maximal (and current) dimension of the vector of coefficients `m_aa`.
    */
   int64_t m_dimaa;

};

//============================================================================
// IMPLEMENTTION

//============================================================================

template<typename Int, typename Real>
Rank1Lattice<Int, Real>::Rank1Lattice(const Int &m, const IntVec &aa, int64_t maxDim, NormType norm) :
      IntLatticeExt<Int, Real>(m, maxDim, norm) {
   this->m_dimaa = aa.length();
   this->setaa(aa);
}

//============================================================================

template<typename Int, typename Real>
Rank1Lattice<Int, Real>::Rank1Lattice(const Int &m, const IntVec &aa, NormType norm) :
      IntLatticeExt<Int, Real>(m, aa.length(), norm) {
   this->setaa(aa);
}

//============================================================================

template<typename Int, typename Real>
Rank1Lattice<Int, Real>::Rank1Lattice(const Int &m, const Int &a, int64_t maxDim,  int64_t dimaa, NormType norm) :
      IntLatticeExt<Int, Real>(m, maxDim, norm) {
   this->m_dimaa = dimaa;
   this->m_aa.SetLength(dimaa);
   this->seta(a);    // Set to a Korobov lattice.
}

//============================================================================

template<typename Int, typename Real>
Rank1Lattice<Int, Real>::Rank1Lattice(const Int &m, const Int &a, int64_t maxDim, NormType norm) :
      IntLatticeExt<Int, Real>(m, maxDim, norm) {
   this->m_dimaa = maxDim;
   this->m_aa.SetLength(maxDim);
   this->seta(a);
}

//============================================================================

template<typename Int, typename Real>
Rank1Lattice<Int, Real>::Rank1Lattice(const Int &m, int64_t maxDim, int64_t dimaa, NormType norm) :
      IntLatticeExt<Int, Real>(m, maxDim, norm) {
   this->m_dimaa = dimaa;
   this->m_aa.SetLength(dimaa);
}

//============================================================================

template<typename Int, typename Real>
Rank1Lattice<Int, Real>::Rank1Lattice(const Int &m, int64_t maxDim, NormType norm) :
      IntLatticeExt<Int, Real>(m, maxDim, norm) {
   this->m_dimaa = maxDim;
   this->m_aa.SetLength(maxDim);
}

//============================================================================

template<typename Int, typename Real>
Rank1Lattice<Int, Real>::~Rank1Lattice() {
   this->m_aa.kill();
}

/*
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
*/

//============================================================================

/*
 * Sets the generating vector to `aa`.
 */
template<typename Int, typename Real>
void Rank1Lattice<Int, Real>::setaa(const IntVec &aa) {
   if (aa[0] != 1) myExit("Rank1Lattice::setaa: must have aa[0] == 1.");
   this->m_dimaa = aa.length();
   this->m_aa = aa;
}

//============================================================================

/*
 * Sets this lattice to a Korobov lattice with multiplier `a` and
 * generating vector of length `dimaa`.
 */
template<typename Int, typename Real>
void Rank1Lattice<Int, Real>::seta(const Int &a, int64_t dimaa) {
   this->m_dimaa = dimaa;
   this->m_aa.SetLength(dimaa);
   this->seta(a);
}

//============================================================================

template<typename Int, typename Real>
void Rank1Lattice<Int, Real>::seta(const Int &a) {
   this->m_aa[0] = 1;
   for (int64_t j = 1; j < this->m_dimaa; j++) {
      this->m_aa[j] = (a * this->m_aa[j - 1]) % this->m_modulo;
   }
}

//============================================================================

template<typename Int, typename Real>
IntVec Rank1Lattice<Int, Real>::getaa() {
   return m_aa;
}

//============================================================================

// An upper-triangular basis is built directly, as explained in the guide of Lattice Tester.
// The dimension `maxDim` of the `IntMat` array is unchanged.
template<typename Int, typename Real>
void Rank1Lattice<Int, Real>::buildBasis(int64_t d) {
   assert(d <= this->m_maxDim);
   this->setDim(d);
   int64_t i, j;
   for (j = 0; j < d; j++)
      this->m_basis[0][j] = this->m_aa[j];
   for (i = 1; i < d; i++)
      for (j = 0; j < d; j++)
         this->m_basis[i][j] = this->m_modulo * (i == j);
   // this->setNegativeNorm();
}

//============================================================================

// This builds the m-dual basis, also in a direct way.
template<typename Int, typename Real>
void Rank1Lattice<Int, Real>::buildDualBasis(int64_t d) {
   assert(d <= this->m_maxDim);
   this->setDimDual(d);
   int64_t i, j;
   this->m_dualbasis[0][0] = this->m_modulo;
   for (j = 1; j < d; j++)
      this->m_dualbasis[0][j] = 0;
   for (i = 1; i < d; i++) {
      this->m_dualbasis[i][0] = -this->m_aa[i];
      for (j = 1; j < d; j++)
         this->m_dualbasis[i][j] = (i == j);
   }
   // this->setDualNegativeNorm();
}

//============================================================================

template<typename Int, typename Real>
void Rank1Lattice<Int, Real>::incDimBasis() {
   int64_t d = 1 + this->getDim();  // New current dimension.
   this->setDim(d);
   int64_t i, j;
   // Add new row and new column in the primal basis.
   for (j = 0; j < d - 1; j++)
      this->m_basis[d - 1][j] = 0;
   this->m_basis[d - 1][d - 1] = this->m_modulo;
   for (i = 0; i < d - 1; i++) {
      this->m_basis[i][d - 1] = (this->m_aa[d - 1] * this->m_basis[i][0]) % this->m_modulo;
   }
   // this->setNegativeNorm();   // These norms are not used in LLL and BKZ.
}

//============================================================================

template<typename Int, typename Real>
void Rank1Lattice<Int, Real>::incDimDualBasis() {
   int64_t d = 1 + this->getDimDual();
   this->setDimDual(d);
   // Add one extra zero coordinate to each vector.
   for (int64_t i = 0; i < d; i++) {
      this->m_dualbasis[i][d - 1] = 0;
      this->m_dualbasis[d - 1][i] = 0;
   }
   this->m_dualbasis[d - 1][0] = -m_aa[d - 1];
   this->m_dualbasis[d - 1][d - 1] = 1;
   // this->setDualNegativeNorm();
}

//============================================================================
template<typename Int, typename Real>
void Rank1Lattice<Int, Real>::buildProjection(IntLattice<Int, Real> &projLattice, const Coordinates &coordSet) {
   // Builds only the primal basis for the projection.
   // We use the method described in the Lattice Tester guide, section 5.5.
   // It *does not* assume that a basis for `this` has been computed before.
   long i, j;
   long d = coordSet.size();     // Number of coordinates in the projection.
   projLattice.setDim(d);
   IntMat &basis = projLattice.getBasis();  // Reference to the projection basis.

   // We first compute the first row of the basis matrix.
   bool case1 = m_aa[*coordSet.begin() - 1] == 1; // We have a_{i_1} = 1.
   if (case1) {  // The easy (usual) case.
      j = 0;
      for (auto it = coordSet.begin(); it != coordSet.end(); it++, j++)
         basis[0][j] = m_aa[*it - 1];
   } else {
      Int c1, b1, b2;   // c1 will be gcd(a_{i_1}, m). We should always have c1=1.
      // Recall:  XGCD (c1, b1, b2, a, m) does c1 = gcd(a, m) = a*b1 + m*b2.
      NTL::XGCD(c1, b1, b2, m_aa[*coordSet.begin() - 1], this->m_modulo);
      if (c1 > 1) myExit("Rank1Lattice::buildProjection: c1 > 1.");
      // We multiply the first row by b1 to recover case 1.
      basis[0][0] = 1;
      auto it = coordSet.begin();
      j = 1;
      for (it++; it != coordSet.end(); it++, j++)
         basis[0][j] = (m_aa[*it - 1] * b1) % this->m_modulo;
   }

   // Then we put the vectors m e_i in the other rows.
   for (i = 1; i < d; i++)
      for (j = 0; j < d; j++)
         basis[i][j] = this->m_modulo * (i == j);
}

//============================================================================
template<typename Int, typename Real>
void Rank1Lattice<Int, Real>::buildProjectionDual(IntLattice<Int, Real> &projLattice,
      const Coordinates &coordSet) {
   // The dual basis is also computed directly, without using the primal.
   // We use the method described in the Lattice Tester guide.
   // It *does not* assume that a primal basis for `this` has been computed before.
   long i, j;
   long d = coordSet.size();     // Number of coordinates in the projection.
   projLattice.setDimDual(d);
   IntMat &dualBasis = projLattice.getDualBasis();
   dualBasis[0][0] = this->m_modulo;
   i = 1;
   auto it = coordSet.begin();
   // We first compute the first column of the m-dual basis.
   bool case1 = m_aa[*coordSet.begin() - 1] == 1; // We have a_{i_1} = 1.
   if (case1) {
      for (it++; it != coordSet.end(); ++it, ++i)
         dualBasis[i][0] = m_aa[*it - 1]; // First column.
   } else {
      std::cout << "buildProjDual: What are we doing in Case 2? \n";
      Int c1, b1, b2; // c1 will be gcd(a_{i_1}, m). We should always have c1=1.
      NTL::XGCD(c1, b1, b2, m_aa[*coordSet.begin() - 1], this->m_modulo);
      if (c1 > 1) myExit("Rank1Lattice::buildProjection: c1 > 1.");
      // We multiply the first column by b1 to recover case 1.
      for (it++; it != coordSet.end(); it++, i++)
         dualBasis[i][0] = -(m_aa[*it - 1] * b1 % this->m_modulo);
   }
   // We put the columns of the identity matrix in the other columns.
   for (i = 0; i < d; i++)
      for (j = 1; j < d; j++)
         dualBasis[i][j] = (i == j);
}

//============================================================================
template<typename Int, typename Real>
std::string Rank1Lattice<Int, Real>::toStringCoef() const {
   return toString(this->m_aa, 0, this->getDim());
}

//============================================================================

template class Rank1Lattice<std::int64_t, double> ;
template class Rank1Lattice<NTL::ZZ, double> ;
template class Rank1Lattice<NTL::ZZ, xdouble> ;
template class Rank1Lattice<NTL::ZZ, quad_float> ;
template class Rank1Lattice<NTL::ZZ, NTL::RR> ;

} // End namespace LatticeTester

#endif
