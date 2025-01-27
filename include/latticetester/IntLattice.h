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

#include <string>
#include <sstream>
#include <iomanip>
#include <cassert>

#include <NTL/xdouble.h>
#include <NTL/RR.h>

#include "latticetester/FlexTypes.h"
#include "latticetester/EnumTypes.h"
#include "latticetester/Util.h"
#include "latticetester/Coordinates.h"
#include "latticetester/BasisConstruction.h"
#include "latticetester/LLL_lt.h"

using namespace NTL;
using namespace LatticeTester;

namespace LatticeTester {

/**
 * \class latticetester/IntLattice.h
 *
 * An `IntLattice` object is an integral lattice, with its basis or its `m`-dual basis, or both.
 * There are tools to perform simple manipulations on those lattice bases.
 * The value of `m` must be chosen in a way that all coordinates of the basis and
 * of its `m`-dual are integers, so they can be represented exactly.
 * The basis or its `m`-dual is rescaled by `m`, which is typically the smallest integer with this property.
 *
 * The dimension `dim` of the lattice is the number of independent vectors that form a basis.
 * Usually, these vectors also have `dim` coordinates, but in general they may have more.
 * The basis and/or the `m`-dual basis are stored in `IntMat` arrays (from %NTL) of sizes `maxDim x maxDim`,
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
 * The class `IntLatticeExt` extends this class and contains virtual methods that must
 * be defined in its subclasses.
 */

template<typename Int, typename Real>
class IntLattice {

public:

   /**
    * Constructs a lattice whose basis is the identity, in `maxDim` dimensions,
    * with the specified norm type, and the scaling factor `m`.
    */
   IntLattice(const Int m, const int64_t maxDim, NormType norm = L2NORM);

   /**
    * Similar to the previous constructor, except that the primal basis is given in
    * `basis`, which must be an `IntMat` object of size `maxDim` by `maxDim`.
    */
   IntLattice(const IntMat basis, const Int m, const int64_t maxDim, NormType norm = L2NORM);

   /**
    * Constructs a lattice with the given basis and given m-dual basis for the given `m`,
    * in `maxDim` dimensions, and with the specified norm type.
    * In this case, by default, both the primal and m-dual basis will be maintained.
    * The two `IntMat` objects must be of size `maxDim` by `maxDim`.
    */
   IntLattice(const IntMat primalbasis, const IntMat dualbasis, const Int m, const int64_t maxDim,
         NormType norm = L2NORM);

   /*
    * Copy constructor. Makes a deep copy of `lat` into `*this` new object.   ****  Remove?
    */
   // IntLattice(const IntLattice<Int, Real> &lat);

   /**
    * Destructor.
    */
   virtual ~IntLattice();

   /*
    * Makes a deep copy of the lattice `lat` into this (existing) object.
    * New matrix and vector objects are constructed to store the bases and norms.
    */
   // void copyLattice(const IntLattice<Int, Real> &lat);

   /*
    * *** Previously named `copyLattice`.
    * Overwrites the first `dim` elements of the basis of the lattice `lat` over the elements
    * of the basis of the current object, in the upper left corner of the basis matrix.
    * Also first `dimdual` elements for the dual.
    * The vector norms are also overwritten.
    * The difference with `copyLattice` is that here, no new matrix or vector is constructed;
    * the previous ones are re-used. Requirement: `dim <= maxDim`.
    */
   // void overwriteLattice(const IntLattice<Int, Real> &lat, long dim, long dimdual);

   /**
    * Constructs a set of generating vectors for the projection of the present lattice,
    * over the set of coordinates determined by `proj`, and builds a basis for that
    * projection using an LLL procedure from `BasisConstruction`.
    * It is assumed that a basis for the present lattice is already available and contains
    * all the coordinates in `proj`.
    * The `maxDim` of the `projLattice` object must be large enough so it can holds the projection.
    * The modulus `m` is assumed to be the same for `projLattice` and for the current object.
    * This function can be called several times with the same `projLattice` object
    * to examine several different projections.
    * Note that representing each projection as an `IntLattice` object is required when
    * we want to call `Reducer::shortestVector` for several projections.
    * This function can be overridden by more efficient ones in certain classes.
    */
   virtual void buildProjectionLLL(IntLattice<Int, Real> &projLattice, const Coordinates &proj,
         double delta = 0.5);

   /**
    * Builds a basis for the projection of the present lattice over the set of coordinates
    * determined by `proj`. This becomes the basis in `projLattice`.
    * By default, it uses an upper-triangular construction.
    * It is assumed that a basis for the present lattice is already available and contains
    * all the coordinates in `proj`.
    */
   virtual void buildProjection(IntLattice<Int, Real> &projLattice, const Coordinates &proj);

   /**
    * Similar to `buildProjection`, except that it builds an upper-triangular basis
    * for the primal and a lower-triangular basis for the $m$-dual of the projection,
    * both in `projLattice`. The dimensions in the latter are updated.
    * We should use this function when we need a basis for the $m$-dual of the projection.
    */
   virtual void buildProjectionDual(IntLattice<Int, Real> &projLattice, const Coordinates &proj);

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
   void setBasis(const IntMat basis, const Int m, const int64_t dim, NormType norm = L2NORM) {
      assert(this->m_maxDim == basis.NumRows());
      m_modulo = m;
      m_norm = norm;
      m_dim = dim;
      m_basis = basis;
      setNegativeNorm();
   }

   /**
    * Changes only the 'basis' and the current dimension to 'dim'.
    * Does not change `m` and `maxDim`.
    */
   void setBasis(const IntMat basis, const int64_t dim) {
      assert(this->m_maxDim == basis.NumRows());
      m_dim = dim;
      m_basis = basis;
      setNegativeNorm();
   }

   /**
    * Returns the m-dual basis represented in a matrix.
    */
   IntMat& getDualBasis() {
      return m_dualbasis;
   }

   /**
    * Returns the maximal dimension for this lattice.
    */
   int64_t getMaxDim() const {
      return m_maxDim;
   }

   /**
    * Returns the current dimension of the primal lattice, which is the dimension of the basis vectors,
    * and also usually the number of independent vectors in the basis.
    */
   int64_t getDim() const {
      return m_dim;
   }

   /**
    * Returns the current dimension of the m-dual lattice.
    */
   int64_t getDimDual() const {
      return m_dimdual;
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
   Int getModulus() const {
      return m_modulo;
   }

   /**
    * Sets the dimension of the primal basis to `dim`. This does not change `maxDim` nor any of the
    * basis vectors, but only the dimension variable.
    */
   void setDim(const int64_t dim) {
      assert(dim <= this->m_maxDim);
      // if (dim > 0)
      m_dim = dim;
   }

   /**
    * Sets the dimension of the primal basis to `dim`. This does not change `maxDim` nor any of the
    * basis vectors, but only the dimension variable.
    */
   void setDimDual(const int64_t dim) {
      // std::cout << "IntLattice::setDimDual: dim = " << dim << ", maxdim = " << this->m_maxDim << std::endl;
      assert(dim <= this->m_maxDim);
      // if (dim > 0)
      m_dimdual = dim;
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
    * Exchanges vectors `i` and `j` in the basis. Also exchanges secondary information
    * about these two vectors (like the norms).
    */
   void permute(int64_t i, int64_t j);

   /**
    * Exchanges vectors `i` and `j` in the `m`-dual basis without changing the primal.
    */
   void permuteDual(int64_t i, int64_t j);

   /**
    * Exchanges the primal and m-dual bases, their dimensions, and vector norms.
    */
   void dualize();

   /**
    * Returns `true` iff the m-dual basis contained in the object really is
    * the m-dual of the current primal basis. Otherwise, it returns false.
    */
   bool checkDuality();

   /**
    * Sorts the primal basis vectors with indices greater of equal to `d` by
    * increasing length. The m-dual vectors are permuted accordingly. Assumes
    * that the lengths (norms) of the corresponding basis vectors are up to date.
    * The dual basis vectors are **not** permuted.
    */
   void sortBasis(int64_t d);

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
    * The current dimension of the primal basis, which is the number of (independent) vectors
    * in it. It cannot exceed the number of coordinates in those vectors.
    * It also cannot exceed `m_maxDim`.
    */
   int64_t m_dim;

   /**
    * The current dimension of the m-dual basis, which is the number of (independent) vectors
    * in its basis. It cannot exceed the number of coordinates in those vectors.
    * It also cannot exceed `m_maxDim`.
    */
   int64_t m_dimdual;

   /**
    * The rows of the m_dim x m_dim upper left corner of this matrix are the primal basis vectors.
    */
   IntMat m_basis;

   /**
    * The rows of the m_dim x m_dim upper left corner this matrix are the m-dual basis vectors.
    * May not be initialized.
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
IntLattice<Int, Real>::IntLattice(const Int m, const int64_t maxDim, NormType norm) {
   m_modulo = m;
   m_maxDim = maxDim;
   m_dim = 0;
   m_dimdual = 0;
   m_norm = norm;
   m_basis.SetDims(maxDim, maxDim);
   m_vecNorm.SetLength(maxDim);
   m_dualbasis.SetDims(maxDim, maxDim);
   m_dualvecNorm.SetLength(maxDim);
   setNegativeNorm();
}

//===========================================================================

/*  // Maybe we remove this.   The IntMat basis does not tell the dimension!
 template<typename Int, typename Real>
 IntLattice<Int, Real>::IntLattice(const IntMat basis, const Int m,
 const int64_t maxDim, NormType norm) {
 //		: m_basis(basis), m_modulo(m), m_dim(dim), m_norm(norm) {
 m_modulo = m;
 m_maxDim = maxDim;
 m_dim = ?????
 // m_dimdual = 0;
 m_norm = norm;
 assert(basis.NumRows() == maxDim);
 m_basis = basis;
 m_vecNorm.SetLength(maxDim);
 // m_dualbasis.SetDims(maxDim, maxDim);
 m_dualvecNorm.SetLength(maxDim);
 setNegativeNorm();
 }
 */
/*=========================================================================*/

/*
 template<typename Int, typename Real>
 IntLattice<Int, Real>::IntLattice(const IntMat primalbasis,
 const IntMat dualbasis, const Int m, const int64_t maxDim,
 NormType norm) :
 IntLattice<Int, Real>(primalbasis, m, maxDim, norm) {
 assert(dualbasis.NumRows() == maxDim);
 m_dualbasis = dualbasis;
 m_dimdual = ??????
 m_dualvecNorm.SetLength(maxDim);
 setDualNegativeNorm();
 }
 */

/*=========================================================================

template<typename Int, typename Real>
IntLattice<Int, Real>::IntLattice(const IntLattice<Int, Real> &lat) {
   copyLattice(lat);
}
*/

/*=========================================================================*/

template<typename Int, typename Real>
IntLattice<Int, Real>::~IntLattice() {
   // kill();
   this->m_basis.IntMat::kill();              // Ok ?
   this->m_dualbasis.IntMat::kill();
   this->m_vecNorm.kill();
   this->m_dualvecNorm.kill();
}

/*=========================================================================
 *
template<typename Int, typename Real>
void IntLattice<Int, Real>::copyLattice(const IntLattice<Int, Real> &lat) {
   this->m_modulo = lat.m_modulo;
   this->m_maxDim = lat.m_maxDim;
   this->m_dim = lat.m_dim;
   this->m_dimdual = lat.m_dimdual;
   this->m_basis = IntMat(lat.m_basis);
   this->m_norm = lat.m_norm;
   this->m_vecNorm = RealVec(lat.m_vecNorm);
   this->m_dualbasis = IntMat(lat.m_dualbasis);
   this->m_dualvecNorm = RealVec(lat.m_dualvecNorm);
}
*/
/*=========================================================================

template<typename Int, typename Real>
void IntLattice<Int, Real>::overwriteLattice(const IntLattice<Int, Real> &lat, long dim,
      long dimdual) {
   if (dim <= this->m_maxDim) {
      CopyMatr(this->m_basis, lat.m_basis, dim);
      CopyVect(this->m_vecNorm, lat.m_vecNorm, dim);
      CopyMatr(this->m_dualbasis, lat.m_dualbasis, dimdual);
      CopyVect(this->m_dualvecNorm, lat.m_dualvecNorm, dimdual);
      this->m_modulo = lat.m_modulo;
      this->m_dim = lat.m_dim;
      this->m_dimdual = lat.m_dimdual;
   } else std::cout << "IntLattice::overwriteLattice: dim > m_maxDim" << std::endl;
}
*/

//===========================================================================

template<typename Int, typename Real>
void IntLattice<Int, Real>::buildProjectionLLL(IntLattice<Int, Real> &projLattice,
      const Coordinates &proj, double delta) {
   // We assume here that this and lattice have the same scaling factor m.
   projLattice.setDim(proj.size());  // Number of coordinates in the projection.
   projectionConstructionLLL<Int, Real>(projLattice.m_basis, this->m_basis,
         proj, this->m_modulo, delta, proj.size());
}

//===========================================================================

// We use the first m_dim rows of the current basis as generating vectors.
template<typename Int, typename Real>
void IntLattice<Int, Real>::buildProjection(IntLattice<Int, Real> &projLattice,
      const Coordinates &proj) {
   // We assume here that this and lattice have the same m.
   projLattice.setDim(proj.size());  // Number of coordinates in the projection.
   projectionConstructionUpperTri<Int>(projLattice.m_basis, this->m_basis, proj, this->m_modulo,
         this->m_dim);
}

//===========================================================================

// This one builds both the primal and the m-dual bases of the projection.
template<typename Int, typename Real>
void IntLattice<Int, Real>::buildProjectionDual(IntLattice<Int, Real> &projLattice,
      const Coordinates &proj) {
   // We assume here that this and lattice have the same m.
   projLattice.setDim(proj.size());  // Number of coordinates in the projection.
   projLattice.setDimDual(proj.size());
   projectionConstructionUpperTri<Int>(projLattice.m_basis, this->m_basis, proj, this->m_modulo,
         this->m_dim);
   mDualUpperTriangular(projLattice.m_dualbasis, projLattice.m_basis, this->m_modulo,
         projLattice.m_dim);
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
   for (int64_t i = 0; i < this->m_dimdual; i++) {
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
      IntVec &row = this->m_basis[i];
      if (this->m_norm == L2NORM) {
         ProdScal<Int>(row, row, this->m_dim, this->m_vecNorm[i]);
      } else {
         CalcNorm<Int, Real>(row, this->m_dim, this->m_vecNorm[i], this->m_norm);
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
void IntLattice<Int, Real>::updateSingleVecNorm(const int64_t &d, const int64_t &c) {
   assert(d >= 0);
   IntVec &row = this->m_basis[d];
   //NTL::matrix_row<Int> row(this->m_basis, d);
   if (this->m_norm == L2NORM) {
      ProdScal<Int>(row, row, c, this->m_vecNorm[d]);
   } else {
      CalcNorm<Int, Real>(row, c, this->m_vecNorm[d], this->m_norm);
   }
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::updateDualVecNorm(const int64_t &d) {
   assert(d >= 0);
   for (int64_t i = d; i < this->m_dimdual; i++) {
      IntVec &row = this->m_dualbasis[i];
      //NTL::matrix_row<Int> row(this->m_dualbasis, i);
      if (this->m_norm == L2NORM) {
         ProdScal<Int>(row, row, this->m_dimdual, this->m_dualvecNorm[i]);
      } else {
         CalcNorm<Int, Real>(row, this->m_dimdual, this->m_dualvecNorm[i], this->m_norm);
      }
   }
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::updateDualVecNorm(const int64_t &d, const int64_t &c) {
   assert(d >= 0);
   for (int64_t i = 0; i < d + 1; i++) {
      IntVec &row = this->m_dualbasis[i];
      //NTL::matrix_row<Int> row(this->m_dualbasis, i);
      if (this->m_norm == L2NORM) {
         ProdScal<Int>(row, row, c, this->m_dualvecNorm[i]);
      } else {
         CalcNorm<Int, Real>(row, c, this->m_dualvecNorm[i], this->m_norm);
      }
   }
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::updateSingleDualVecNorm(const int64_t &d, const int64_t &c) {
   assert(d >= 0);
   IntVec &row = this->m_dualbasis[d];
   //NTL::matrix_row<Int> row(this->m_dualbasis, d);
   if (this->m_norm == L2NORM) {
      ProdScal<Int>(row, row, c, this->m_dualvecNorm[d]);
   } else {
      CalcNorm<Int, Real>(row, c, this->m_dualvecNorm[d], this->m_norm);
   }
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::updateScalL2Norm(const int64_t i) {
   IntVec &row = this->m_basis[i];
   // NTL::matrix_row<Int> row(this->m_basis, i);
   ProdScal<Int>(row, row, this->m_dim, this->m_vecNorm[i]);
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::updateScalL2Norm(const int64_t k1, const int64_t k2) {
   for (int64_t i = k1; i < k2; i++) {
      updateScalL2Norm(i);
   }
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::updateDualScalL2Norm(const int64_t i) {
   IntVec &row = this->m_dualbasis[i];
   // NTL::matrix_row<Int> row(this->m_dualbasis, i);
   ProdScal<Int>(row, row, this->m_dimdual, this->m_dualvecNorm[i]);
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::updateDualScalL2Norm(const int64_t k1, const int64_t k2) {
   for (int64_t i = k1; i < k2; i++) {
      updateDualScalL2Norm(i);
   }
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::permute(int64_t i, int64_t j) {
   if (i == j) return;
   for (int64_t k = 0; k < this->m_dim; k++) {
//      std::cout << " IntLattice::permute, m_basis(j, k), m_basis(i, k) = " <<
//            this->m_basis(j, k) << "  " << this->m_basis(i, k) << "\n";
      swap9(this->m_basis[j][k], this->m_basis[i][k]);
   }
   swap9(this->m_vecNorm[i], this->m_vecNorm[j]);
}

/*=========================================================================*/

template<typename Int, typename Real>
void IntLattice<Int, Real>::permuteDual(int64_t i, int64_t j) {
   if (i == j) return;
   for (int64_t k = 0; k < this->m_dimdual; k++) {
      swap9(this->m_dualbasis[j][k], this->m_dualbasis[i][k]);
   }
   swap9(this->m_dualvecNorm[i], this->m_dualvecNorm[j]);
}

//===========================================================================

template<typename Int, typename Real>
void IntLattice<Int, Real>::dualize() {
   std::swap(this->m_basis, this->m_dualbasis);
   std::swap(this->m_dim, this->m_dimdual);
   std::swap(this->m_vecNorm, this->m_dualvecNorm);
}

/*=========================================================================*/

template<typename Int, typename Real>
bool IntLattice<Int, Real>::checkDuality() {
   Int S;
   int64_t dim = getDim();
   for (int64_t i = 0; i < dim; i++) {
      for (int64_t j = 0; j < dim; j++) {
         IntVec row1 = this->m_basis[i];
         IntVec row2 = this->m_dualbasis[j];
         //NTL::matrix_row<const IntMat> row1(this->m_basis, i);
         //NTL::matrix_row<const IntMat> row2(this->m_dualbasis, j);
         ProdScal<Int>(row1, row2, dim, S);
         if (j != i) {
            if (S != 0) {
               std::cout << "******  checkDuality failed for V[" << i << "] and W[" << j << "]"
                     << std::endl;
               return false;
            }
         } else if (S != this->m_modulo) {
            std::cout << "******  checkDuality failed for i, j = " << i << " , " << j << std::endl;
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
         std::cout << "\n***** ERROR: in sort, Negative norm for i = " << i << ",  dim = " << dim
               << std::endl;
      }
   }
   // This is a rather inefficient sort, in O(d^2) operations!
   for (int64_t i = d; i < dim; i++) {
      int64_t k = i;
      for (int64_t j = i + 1; j < dim; j++) {
         if (getVecNorm(j) < getVecNorm(k)) k = j;
      }
      if (i != k) permute(i, k);
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
   for (int64_t i = 0; i < this->m_dimdual; i++) {
      // if (this->m_withDual) {
      os << this->m_dualbasis[i];
      //for (int64_t j = 0; j < this->m_dimdual; j++) {
      //  os << this->m_dualbasis(i,j);
      //}
      os << "\n";
      // }
   }
   os << "\n";
   os << "Norm used: " << toStringNorm(this->m_norm) << "\n" << std::endl;
   os << "Norm of each Basis vector: \n";
   os << "Primal";
   // if (this->m_withDual)
   os << "\t\tDual\n";
   os << "\n";

   for (int64_t i = 0; i < this->m_dim; i++) {
      if (this->m_vecNorm[i] < 0) {
         os << "NaN OR Not computed";
      } else {
         if (this->m_norm == L2NORM) {
            os << std::sqrt(conv<double>(this->m_vecNorm[i]));
            // os << sqrtReal(this->m_vecNorm[i]);
         } else {
            os << this->m_vecNorm[i];
         }
      }
      os << "\t";
      // if (this->m_withDual) {
      if (this->m_dualvecNorm[i] < 0) os << "NaN OR Not computed";
      else {
         if (this->m_norm == L2NORM) {
            os << std::sqrt(conv<double>(this->m_dualvecNorm[i]));
            // os << sqrtReal(this->m_dualvecNorm[i]);
         } else {
            os << this->m_dualvecNorm[i];
         }
      }
      // }
      os << "\n";
   }
   os << std::endl;
   return os.str();
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
         os << " " << std::setprecision(15) << this->m_basis[i][j];
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
   os << "  Dim = " << this->m_dimdual << " \n";
   for (int64_t i = 0; i < this->m_dimdual; i++) {
      os << "    [";
      for (int64_t j = 0; j < this->m_dimdual; j++)
         os << " " << std::setprecision(15) << this->m_dualbasis[i][j];
      os << " ]\n";
   }
   os << "  Norms:\n";
   os << "    [";
   for (int64_t i = 0; i < this->m_dimdual; i++) {
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
template class IntLattice<NTL::ZZ, xdouble> ;
template class IntLattice<NTL::ZZ, quad_float> ;
template class IntLattice<NTL::ZZ, NTL::RR> ;

} // namespace LatticeTester

#endif
