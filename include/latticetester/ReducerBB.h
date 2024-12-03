// This file is part of LatticeTester.
//
// Copyright (C) 2012-2022  The LatticeTester authors, under the supervision
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

#ifndef LATTICETESTER_REDUCERBB_H
#define LATTICETESTER_REDUCERBB_H

#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <cstdint>
#include <iostream>
#include <ctime>
#include <string>
#include <iomanip>
#include <set>
#include <map>
#include <cmath>
#include <cstdlib>
#include <type_traits>

#include "NTL/LLL.h"
#include "NTL/tools.h"
#include "NTL/ZZ.h"
#include "NTL/RR.h"

#include "latticetester/EnumTypes.h"
#include "latticetester/Util.h"
#include "latticetester/IntLattice.h"
#include "latticetester/BasisConstruction.h"
//#include "latticetester/LLL_FPInt.h"
//#include "latticetester/LLL_RR_lt.h"

namespace LatticeTester {

/**
 * This `ReducerBB` class provides facilities to reduce the basis of a lattice
 * (an `IntLattice` object) in various ways (pairwise, LLL, BKZ, Minkowski
 * \cite rDIE75a, \cite mLEN82a, \cite mSCH91a),
 * and to find a shortest nonzero vector in the lattice using a BB algorithm \cite rFIN85a.
 * Each `Reducer` must have an internal `IntLattice` object which is given upon construction
 * and can also be changed later via `setIntLattice`.
 * The reduction methods are applied to this internal object.
 *
 * The `shortestVector` and `reductMinkowski` methods do not apply any pre-reduction by themselves
 * Before calling them, one should always reduce the basis separately beforehand
 * with an LLL or BKZ reduction, because it drastically reduces the size of the BB search.
 * These two methods have no static version.
 * To use them one must create `Reducer` object, which maintains several
 * internal variables, vectors, and matrices.  The recommended way is to create
 * a single `Reducer` object with a maximal dimension large enough,
 * and then call the `shortestVector` and `reductMinkowski` methods with the relevant
 * `IntLattice` object as a parameter.
 * *****  In fact, it can be always the same IntLattice object all the time!!!  *****
 * The norm type, dimension, basis, vector lengths, etc.
 * will be taken from this `IntLattice` object.  The dimensions of the internal vectors
 * and matrices can be larger than required; the methods will use only the
 * entries that are needed for the given `IntLattice` basis.
 * Creating a new `Reducer` object for each `IntLattice` that we want to handle is very
 * inefficient and should be avoided, especially when we want to examine several
 * projections for several lattices.
 *
 * Most of the methods do not use or change the m-dual lattice.
 * To reduce the m-dual basis or find a shortest nonzero vector in it,
 * one should first dualize the lattice (`IntLattice::dualize` does that)
 * and then apply the desired methods.
 */

template<typename Int, typename Real>
class ReducerBB {

// using namespace LatticeTester;

public:

   /**
    * Constructor that initializes the reducer to work on the lattice `lat`.
    * The maximal dimension will be that of `lat`.
    */
   ReducerBB(IntLattice<Int, Real> &lat);

   /**
    * Constructor of a `ReducerBB` than can handle up to `maxDim` dimensions.
    * Space will be reserved once for all for the internal `IntMat` objects such as
    * the Gram and Cholesky matrices. These same matrices will be re-used over and over
    * to avoid repeated space reallocation.  It is recommended to create a reducer with
    * a large enough max dim in the first place, using this constructor.
    */
   ReducerBB(int64_t maxDim);

   /**
    * Copy constructor.
    */
   ReducerBB(const ReducerBB<Int, Real> &red);

   /**
    * Destructor.
    */
   ~ReducerBB();

   /**
    * Assignment operator that makes a deep copy of `red`
    * into the current object, using `copy`.
    */
   ReducerBB<Int, Real>& operator=(const ReducerBB<Int, Real> &red);

   /**
    * Initializes all matrices used in the following.
    */
   void init(int64_t maxDim);

   /**
    * Copies `red` into the current object.
    */
   void copy(const ReducerBB<Int, Real> &red);

   /**
    * Computes a shortest non-zero vector for the `IntLattice` stored in this `ReducerBB` object,
    * with respect to the norm in this `IntLattice`,
    * using the BB algorithm described in \cite rLEC97c and \cite iLEC22l.
    * The admissible norm types here are `L1NORM` and `L2NORM` (see `EnumTypes.h`).
    * If the constant `ReducerBB::MaxNodesBB` (see below) is exceeded
    * during the branch-and-bound, the method aborts and returns
    * `false`. Otherwise, it returns `true`. If the reduction was
    * successful, the new reduced basis can be accessed via `getIntLattice()`.
    * The shortest vector is placed in first position and its norm is updated.
    * The other vectors are not sorted.
    *
    * This function uses only the basis of the internal lattice, its vector
    * lengths, and scalar products. It never uses its m-dual.
    * To compute a shortest vector in the m-dual, one must first call `dualize`
    * on the target `IntLattice` object.
    * It is strongly recommended to use `redBKZ` or `redLLLNTL` to pre-reduce
    * the basis before invoking this method; this is not done automatically.
    * It is assumed that `getDim` returns the correct dimension of the internal lattice.
    */
   bool shortestVector();

   /**
    * In this version, the lattice is passed as a parameter.
    * It will become the new `IntLattice` of this `ReducerBB` object.
    * This method calls `setIntLattice(lat)`, so if the max dimension for the
    * ReducerBB is not large enough for `lat`, all the internal variables of this
    * reducer will be reset and the vectors and matrices will be automatically enlarged.
    * In particular, the bounds set by `setBoundL2` have to be reset.   *****
    * If the max dimension is large enough, only a pointer is changed.
    */
   bool shortestVector(IntLattice<Int, Real> &lat);

   /**
    * Reduces the current basis to a Minkowski-reduced basis with respect
    * to the Euclidean norm, assuming that the first \f$d\f$ vectors are
    * already reduced and sorted. If `MaxNodesBB` is exceeded during one
    * of the branch-and-bound step, the method aborts and returns `false`.
    * Otherwise it returns `true`, the basis is reduced and sorted by
    * increasing vector lengths. For a full reduction, just omit the `d` parameter.
    */
   bool reductMinkowski(int64_t d = 0);

   /**
    * In this version, the lattice is passed as a parameter.
    * It will become the new `IntLattice` of this `ReducerBB` object,
    * exactly as in `shortestVector(IntLattice)`.
    */
   bool reductMinkowski(IntLattice<Int, Real> &lat, int64_t d = 0);

   /**
    * Returns the current shortest vector found in this lattice.  It may not be in the basis.
    */
   IntVec getShortVec() {
      return m_bv;
   }

   /**
    * Returns the *square Euclidean length* of the current shortest basis vector in the lattice.
    */
   Real getMinLength2() {
      return m_lMin2;
   }

   /**
    * Returns the length of the current shortest basis vector in the lattice,
    * which is stored in a local variable.
    * This depends only on the lattice, but this length is stored and used in this class.
    */
   Real getMinLength() {
      if (m_lat->getNormType() == L2NORM) return sqrt(m_lMin2);
      else return m_lMin;
   }

   /**
    * Returns the length of the current *last* basis vector in the lattice.
    */
   Real getMaxLength() {
      if (m_lat->getNormType() == L2NORM) return sqrt(m_lat->getVecNorm(m_lat->getDim() - 1));
      else return m_lat->getVecNorm(m_lat->getDim() - 1);
   }

   /**
    * Sets a vector of bounds on the square of the acceptable shortest
    * vector lengths (for the Euclidean norm),
    * in dimensions from `dim1+1` to `dim2`. `thresholds[i]` must
    * contain a lower bound on the square of the length of the shortest
    * vector in the lattice in dimension `i+1`.
    * This bound will be used during during the Branch-and-Bound step
    * when computing a shortest lattice vector.  As soon as a vector shorter
    * than the bound is found, the BB algorithm will stop. This is useful
    * when we search for a good lattice with the spectral test, since
    * it can reduce the work considerably.
    * If these bounds are not set, the default values of 0 are used.
    * It is recommended to set these bounds before calling `shortestVector`
    * for the first time when making searches for good lattices.
    */
   void setBoundL2(const RealVec &thresholds, int64_t dim1, int64_t dim2);

   /**
    * Sets to `lat` the `IntLattice` object on which this reducer object will be working.
    * If it is already `lat`, nothing is done.
    * If the dimension of `lat` is larger than the max dimension for
    * this `ReducerBB`, the latter is increased and the `ReducerBB` is re-initialized
    * with new internal variables having the appropriate dimensions.
    * Otherwise, the internal variables are left unchanged.
    * It is recommended to create the reducer with a large enough max dim in the first place.
    */
   void setIntLattice(IntLattice<Int, Real> &lat) {
      m_lat = &lat;
      if (m_lat->getDim() > m_maxDim) init(m_lat->getDim());
   }

   /**
    * Returns the IntLattice object that this object is working on.
    */
   IntLattice<Int, Real>* getIntLattice() {
      return m_lat;
   }

   /**
    * Sets the decomposition method used by this reduced to compute bounds
    * in the BB procedures.
    */
   void setDecompTypeBB(DecompTypeBB decomp) {
      m_decomp = decomp;
   }

   /**
    * The maximum number of nodes in the branch-and-bound tree when
    * calling `shortestVector` or `reductMinkowski`. When this number is
    * exceeded, the method aborts and returns `false`.
    */
   int64_t maxNodesBB = 10000000;


private:

   /**
    * Finds a shortest vector of the primal basis using branch-and-bound (BB).
    * This is called by `shortestVector`.
    */
   bool redBBShortVec();

   /**
    * Recursive procedure that tries to find a shorter vector using BB.
    * It is called by `redBBShortVec`. The parameter j indicates the level of the
    * BB tree. The first call is for level i=0.
    */
   bool tryZShortVec(int64_t j, bool &smaller, NormType norm);

   /**
    * Used by `reductMinkowski` to try shorten the vectors of the primal basis using
    * branch-and-bound.
    */
   bool redBBMink(int64_t i, int64_t d, int64_t Stage, bool &smaller, bool taboo[] = NULL);

   /**
    * Recursive procedure used in `redBBMink` to try find shorter vectors.
    */
   bool tryZMink(int64_t j, int64_t i, int64_t Stage, bool &smaller, const IntMat &WTemp);

   /**
    * Performs pairwise reductions. This method tries to reduce each basis
    * vector with index larger than \f$d\f$ and distinct from \f$i\f$ by
    * adding to it a multiple of the \f$i\f$-th vector. Always uses the
    * Euclidean norm.
    */
   void pairwiseRedPrimal(int64_t i, int64_t d, bool taboo[] = NULL);

   /**
    * Performs pairwise reductions, trying to reduce every other vector of
    * the *dual* basis by adding multiples of the \f$i\f$-th vector. That
    * may change the \f$i\f$-th vector in the primal basis. Each such dual
    * reduction is actually performed only if that does not increase the
    * length of vector \f$i\f$ in the primal basis. Always uses the
    * Euclidean norm.
    */
   // void pairwiseRedDual(int64_t i, bool taboo[] = NULL);

   /**
    * Computes a Cholesky decomposition of the basis. Returns in `C0` the
    * elements of the upper triangular matrix of the Cholesky
    * decomposition that are above the diagonal. Returns in `DC2` the
    * squared elements of the diagonal.
    */
   bool calculCholesky(RealVec &DC2, RealMat &C0);

   /**
    * Recalculates the first \f$n\f$ entries of the \f$j^{th}\f$ column of
    * the Cholesky matrix of order 2.
    */
   void calculCholesky2LLL(int64_t n, int64_t j);

   /**
    * Recalculates the entry (\f$i\f$, \f$j\f$) of the Cholesky matrix of order 2.
    */
   void calculCholesky2Ele(int64_t i, int64_t j);

   /**
    * Permutes the \f$i^{th}\f$ and the \f$j^{th}\f$ line, and the
    * \f$i^{th}\f$ and the \f$j^{th}\f$ column of the scalar products matrix.
    */
   void permuteGramVD(int64_t i, int64_t j, int64_t n);

   /*
    * Recomputes the element in row `j` and column `j` of the matrix of scalar products.
    */
   void miseAJourGramVD(int64_t j);

   /*
    * Computes and stores in `m_gramVD` the matrix of scalar products `m_lat->V[i]*m_lat->V[j]`.
    * Equivalent to computing `m_lat->V * transpose(m_lat->V)`.  Used only by redLLLOld.
    */
   void calculGramVD();

   /**
    * Reduces the Cholesky matrix by adding a multiple of the i-th vector
    * to the j-th vector. It updates the Gram-Schmidt matrix.
    */
   void reductionFaible(int64_t i, int64_t j);

   /**
    * This is an old implementation of LLL translated from an old Modula-2 code.
    * It is considerably slower than the NTL versions when Int = ZZ.
    * We leave it here mostly to enable comparisons.
    * This function performs a LLL basis reduction with factor `delta` \cite iLEC22l.
    * The reduction is applied to the first `dim` basis vectors and coordinates
    * when `dim > 0`, and to the entire basis (all vectors) when `dim=0`.
    * Note that in the former case, the transformations are not applied to all the
    * columns, so we no longer have a consistent basis for the original lattice,
    * so we should make internal copies just in case dim > 0.
    * If we want a reduced basis for a subset of coordinates, we should first make
    * a copy to get a basis for these coordinates, before invoking this LLL.
    *
    * This function always uses the Euclidean norm.
    * The factor `delta` must be between 1/2 and 1. The closer it is to 1,
    * the more the basis is reduced, in the sense that the LLL
    * algorithm will enforce tighter conditions on the basis,
    * but this will require more work. The reduction algorithm is
    * applied until `maxcpt` successful transformations have been done,
    * or until the basis is correctly LLL reduced with factor `delta`.
    */
   void redLLLOld(double delta = 0.999999, std::int64_t maxcpt = 1000000000, int64_t dim = 0);

   /*
    * In this function, we assume that we have found a new shorter vector
    * \f$ \bv = \sum_{i=1}^t z_i \bv_i\f$ and we want to insert it in the basis.
    * If \f$z_j = \pm 1\f$ for some \f$j\f$, we can simply exchange \f$\bv\f$ with  \f$\bv_j\f$
    * in the basis. Otherwise, we transform the basis and the \f$z_j\f$'s so that one of the
    * nonzero \f$z_j\f$'s becomes equal to  \f$\pm 1\f$, as explained in the guide, and then
    * we make the exchange. After that, we exchange  \f$\bv\f$ with the vector in first position.
    * Note that some of the basis vectors may become very large.
    * This corresponds to the ``transformation of stage 3'' described in \cite rAFF85a.
    */
   void insertBasisVector(std::vector<std::int64_t> &z);

   /*
    * This function provides an alternative to `insertBasisVector` when we use the L2 norm.
    * It adds the new vector \f$\bv = \sum_{i=1}^t z_i \bv_i\f$ to the basis to form a set
    * of `dim+1` generating vectors, and it applies LLL to recover a new basis that contains \f$\bv\f$.
    */
   void insertBasisVectorLLL(std::vector<std::int64_t> &z);

   /**
    * Debug function that sorts and prints the primal and dual bases
    * to standard output, using the `write` function.
    */
   void tracePrintBases(char *mess);

   /**
    * The lattice that this object is working on.
    */
   IntLattice<Int, Real> *m_lat;

   /**
    * A vector that contains a lower bound on the acceptable (squared) length
    * of the shortest vector, for each number of dimensions.
    * If any vector of the lattice is shorter than this bound,
    * we stop the reduction immediately and reject this lattice since
    * we already know that its shortest vector is too small.
    * This is used in `RedBBShortVec`, and useful for the seek procedure in LatMRG.
    * We could also just check this after the BKZ reduction and after the BB,
    * in the seek procedures.
    */
   RealVec m_BoundL2;

   /**
    * Whenever the number of nodes in the branch-and-bound tree exceeds
    * <tt>MINK_LLL</tt> in the method <tt>reductMinkowski</tt>,
    * `PreRedLLLMink` is automatically set to `true` for the next call;
    * otherwise it is set to `false`.
    */
   std::int64_t MINK_LLL = 500000;

   /**
    * Pre-reduction flag for `reductMinkowski`.
    * When true, LLL is performed automatically at certain steps of the reduction.
    */
   bool PreRedLLLMink = true;

   /**
    * Decomposition method used for computing bounds in the BB procedures.
    * can be ``CHOLESKY` (the default) or `TRIANGULAR`.
    */
   DecompTypeBB m_decomp = CHOLESKY;

   /*
    * Local working variables for this class.
    * They are used inside the basis reduction and short vector methods, and
    * are declared here to avoid passing them as parameters across the methods.
    * The matrices and vectors are sized to some maximum dimensions in init(),
    * which must be large enough for all computations handled by this ReducerBB object.
    */
   int64_t m_maxDim, m_dim; // maximum dimension and current dimension.
   IntVec m_bv;   // Saves current shortest vector in primal basis
   // IntVec m_bw;   // Saves current shortest vector in dual basis (for Mink) ??
   Real m_lMin;   // The norm of the shortest vector in the primal basis
                  // according to the norm considered
   Real m_lMin2;  // Squared L2-norm of the shortest vector in the primal basis.

   Int m_bs;      // Used in pairwiseRedPrimal and pairwiseRedDual.
   Real m_ns;     // Used in pairwiseRedPrimal and pairwiseRedDual.
   // RealVec m_nv;  // Used in pairwiseRedDual.

   RealVec m_n2, m_dc2; // Vectors used in BB algorithms, in tryZShortVec and tryZMink.
   RealMat m_c0, m_c2;  // Matrices used in the Cholesky decomposition.
                        // We must avoid resizing them, because this is expensive!
   int64_t *m_IC;       // Indexes used in Cholesky

   std::vector<std::int64_t> m_zLI;   // Vector of values of z_i.
   RealVec m_zLR;     // Same values of z_i in floating point;
                      // we need them in FP when calculating bounds.
   std::vector<std::int64_t> m_zShort;  // Values of z_i for shortest vector.
   std::int64_t m_countNodes = 0;  // Number of visited nodes in the BB tree
   std::int64_t m_cpt;  // Number of successes in pre-reduction transformations
   bool m_foundZero;    // = true -> the zero vector has been handled

   RealMat m_cho2, m_gramVD;  // Matrices used by redLLLOld.

};
// End of class ReducerBB


//=========================================================================
// IMPLEMENTATION

// Constructor for a given `IntLattice` object `lat`.
template<typename Int, typename Real>
ReducerBB<Int, Real>::ReducerBB(IntLattice<Int, Real> &lat) {
   m_lat = &lat;
   init(m_lat->getDim());
}

//=========================================================================

// Constructor.
template<typename Int, typename Real>
ReducerBB<Int, Real>::ReducerBB(int64_t maxDim) {
   init(maxDim);
}

//=========================================================================

// Initialization.
template<typename Int, typename Real>
void ReducerBB<Int, Real>::init(int64_t maxDim) {
   m_maxDim = maxDim;
   int64_t dim1 = maxDim;
   int64_t dim2 = maxDim;
   if (dim2 <= 2) dim2++;
   m_c0.SetDims(dim1, dim1);
   m_c2.SetDims(dim1, dim1);
   m_cho2.SetDims(dim2, dim2);
   m_gramVD.SetDims(dim2, dim2);
   // Indices in m_IC can go as high as maxDim + 2.
   m_IC = new int64_t[2 + dim1];

   //m_nv.SetLength(dim1);
   m_bv.SetLength(dim1);
   // m_bw.SetLength(dim1);
   m_n2.SetLength(dim1);
   m_zLR.SetLength(dim1);
   m_zLI.resize(dim1);
   m_zShort.resize(dim1);
   m_dc2.SetLength(dim1);
   m_BoundL2.SetLength(dim1);

   m_lMin = std::numeric_limits<double>::max();
   m_lMin2 = m_lMin;
   for (int64_t i = 0; i < dim1; i++) {
      m_zLI[i] = -1;
      m_zShort[i] = -1;
      m_BoundL2[i] = -1;
      m_IC[i] = -1;
   }
   m_countNodes = 0;
   m_foundZero = false;
   m_cpt = 0;
   PreRedLLLMink = false;
   maxNodesBB = 1000000000;
}

//=========================================================================

template<typename Int, typename Real>
ReducerBB<Int, Real>::ReducerBB(const ReducerBB<Int, Real> &red) {
   copy(red);
}

//=========================================================================

template<typename Int, typename Real>
ReducerBB<Int, Real>&
ReducerBB<Int, Real>::operator=(const ReducerBB<Int, Real> &red) {
   if (this != &red) copy(red);
   return *this;
}

//=========================================================================

template<typename Int, typename Real>
void ReducerBB<Int, Real>::copy(const ReducerBB<Int, Real> &red) {
   m_lat = red.m_lat;
   m_maxDim = red.m_maxDim;
   m_c0 = red.m_c0;
   m_c2 = red.m_c2;
   m_dc2 = red.m_dc2;
   // m_nv = red.m_nv;
   m_bv = red.m_bv;
   // m_bw = red.m_bw;
   m_n2 = red.m_n2;
   m_zLR = red.m_zLR;
   m_zLI = red.m_zLI;
   m_zShort = red.m_zShort;
   m_cho2 = red.m_cho2;
   m_gramVD = red.m_gramVD;
   m_lMin = red.m_lMin;
   m_lMin2 = red.m_lMin2;
   m_BoundL2 = red.m_BoundL2;
   m_countNodes = 0;
   m_foundZero = false;
   m_cpt = 0;
   PreRedLLLMink = false;
   maxNodesBB = 1000000000;
   if (m_IC != 0) delete[] m_IC;
   m_IC = new int64_t[3 + m_lat->getDim()];
   for (int64_t i = 0; i < 2 + m_lat->getDim(); i++)
      m_IC[i] = red.m_IC[i];
}

//=========================================================================

template<typename Int, typename Real>
ReducerBB<Int, Real>::~ReducerBB() {
   m_c0.kill();
   m_c2.kill();
   m_cho2.kill();
   m_gramVD.kill();
   // m_nv.kill();
   m_bv.kill();
   // m_bw.kill();
   m_n2.kill();
   m_zLR.kill();
   m_zLI.resize(0);
   m_zShort.resize(0);
   m_dc2.kill();
   m_BoundL2.kill();
   delete[] m_IC;
}

//=========================================================================

template<typename Int, typename Real>
void ReducerBB<Int, Real>::setBoundL2(const RealVec &thresholds, int64_t dim1, int64_t dim2) {
   m_BoundL2.SetLength(dim2);
   for (int64_t i = dim1; i < dim2; i++)
      m_BoundL2[i] = thresholds[i];
}

//=========================================================================

void negativeCholesky() {
   std::cout << "\n***** Negative diagonal element in Cholesky Decomposition\n" << std::endl;
}

template<typename Int, typename Real>
bool ReducerBB<Int, Real>::calculCholesky(RealVec &DC2, RealMat &C0) {
   /*
    * Returns in C0 the elements of the upper triangular matrix of the
    * Cholesky decomposition that are above the diagonal. Returns in DC2 the
    * squared elements of the diagonal. These elements are rescaled by EchVV
    * when SISquares= true.
    */
   const int64_t dim = m_lat->getDim();
   int64_t k, j, i;
   Real m2;
   // C2[i][j] = C0[i][j] * C2[i][i] if i != j.
   // C2[i][i] = DC2[i].
   NTL::conv(m2, m_lat->getModulus());
   m2 = m2 * m2;
   int64_t d = dim;
   //if(m_lat->withDual())
   //  d = dim / 2; // If we use the Dual, we compute Cholesky
   // with the Dual

   // Compute the d first lines of C0 with the primal Basis.
   for (i = 0; i < d; i++) {
      m_lat->updateScalL2Norm(i);
      IntVec &row1 = m_lat->getBasis()[i];
      // NTL::Mat_row<Int> row1(m_lat->getBasis(), i);
      for (j = i; j < dim; j++) {
         if (j == i) NTL::conv(m_c2[i][i], m_lat->getVecNorm(i));
         else {
            IntVec &row2 = m_lat->getBasis()[j];
            // NTL::Mat_row<Int> row2(m_lat->getBasis(), j);
            ProdScal<Int>(row1, row2, dim, m_c2[i][j]);
         }
         for (k = 0; k < i; k++)
            m_c2[i][j] -= C0[k][i] * m_c2[k][j];
         if (i == j) {
            DC2[i] = m_c2[i][i];
            if (DC2[i] < 0.0) {
               negativeCholesky();
               return false;
            }
         } else C0[i][j] = m_c2[i][j] / DC2[i];
         // add for test
         /**  if(i!=j && i<j)
          std::cout<< "C0("<<i<<","<<j<<")="<<C0[i][j]<<" ";
          else if (i==j)
          std::cout<< "C0("<<i<<","<<j<<")="<<DC2[i]<<" ";
          else
          std::cout<< "C0("<<i<<","<<j<<")="<<"0"<<" ";	*/
      }
      // std::cout<<""<<std::endl;
   }
   return true;
}

// =========================================================================

template<typename Int, typename Real>
void ReducerBB<Int, Real>::insertBasisVector(std::vector<std::int64_t> &z) {
   int64_t i, j, h;
   std::int64_t q;
   const int64_t dim = m_lat->getDim();
   j = dim - 1;
   while (z[j] == 0)
      --j;
   while (abs(z[j]) > 1) {
      i = j - 1;
      while (z[i] == 0)
         --i;
      // On a 2 indices i < j tels que |z_j| > 1 et z_i != 0.
      while (z[j]) {
         // Truncate the quotient towards 0.
         q = z[i] / z[j];
         if (q) {
            // On ajoute q * v[i] au vecteur v[j]
            z[i] -= q * z[j];
            IntVec &row1 = m_lat->getBasis()[j];
            IntVec &row2 = m_lat->getBasis()[i];
            ModifVect(row1, row2, q, dim);
            // ModifVectModulo(row1, row2, q, m_lat->getModulus(), dim);  // In Util.h
            m_lat->setNegativeNorm(j);
         }
         // Permute the two vectors.
         std::swap(z[i], z[j]);
         m_lat->permute(i, j);
      }
      j = i;
   }
   // Put the new vector v in place of v_j, then in first place.
   IntVec &rowj = m_lat->getBasis()[j];
   for (h = 0; h < dim; h++)
      rowj[h] = m_bv[h];
   m_lat->permute(0, j);
   m_lat->setVecNorm(m_lMin2, 0);
   // m_lat->setNegativeNorm(0);
}

//=========================================================================

template<typename Int, typename Real>
void ReducerBB<Int, Real>::insertBasisVectorLLL(std::vector<std::int64_t> &z) {
   // Add the new vector `v` as a new row to the basis.
   const int64_t dim = m_lat->getDim();
   IntVec &newrow = m_lat->getBasis()[dim];
   for (int64_t h = 0; h < dim; h++)
      newrow[h] = m_bv[h];
   LLLConstruction0<Int, Real>(m_lat->getBasis(), 0.5, dim+1, dim);  //, m_lat->getVecNorm()*);
}

//=========================================================================

// If m_countNodes > MaxNodesBB, returns false, otherwise returns true.
template<typename Int, typename Real>
bool ReducerBB<Int, Real>::tryZMink(int64_t j, int64_t i, int64_t Stage, bool &smaller,
      const IntMat &WTemp)  {
   std::int64_t max0, min0;
   Real x, dc;
   Real center;
   std::int64_t zhigh, zlow, h;
   bool high;
   int64_t k;
   Real S1, S2, S3, S4, mR;

   const int64_t dim = m_lat->getDim();
   NTL::conv(mR, m_lat->getModulus());

   ++m_countNodes;
   if (m_countNodes > maxNodesBB) {
      std::cerr << "-------- m_countNodes > maxNodesBB = " << maxNodesBB << std::endl;
      return false;
   }

   // Calcul d'un intervalle contenant les valeurs admissibles de z_j.
   center = 0.0;
   if (j < dim - 1) {
      // Calcul du centre de l'intervalle.
      for (k = j + 1; k < dim; k++)
         center = center - m_c0(j, k) * m_zLR[k];

      // Distance du centre aux extremites de l'intervalle.
      // We use the L2 norm for this, since other norms are not allowed.
      dc = sqrt((m_lMin2 - m_n2[j]) / m_dc2[j]);

      /* Calcul de deux entiers ayant la propriete qu'un entier */
      /* non-compris entre (i.e. strictement exterieur `a) ceux-ci */
      /* n'appartient pas a l'intervalle.  */
      x = center - dc;
      NTL::conv(min0, trunc(x));
      if (x > 0.0) ++min0;

      x = center + dc;
      NTL::conv(max0, trunc(x));
      if (x < 0.0) --max0;

      // En vue du choix initial de z_j. On determine zlow et zhigh.
      if (min0 > max0) {
         return true;
      }
      if (min0 == max0) {
         zlow = min0;
         zhigh = max0 + 1;
         high = false;
      } else if (center >= 0.0) {
         NTL::conv(zlow, trunc(center));
         zhigh = zlow + 1;
         NTL::conv(h, trunc(2.0 * center));
         high = h & 1;
      } else {
         NTL::conv(zhigh, trunc(center));
         zlow = zhigh - 1;
         NTL::conv(h, -trunc(2.0 * center));
         high = (h & 1) == 0;
      }
   } else {  // j = dim-1
      zlow = 0;
      high = true;
      if (Stage == 2) {
         min0 = 1;
         max0 = 1;
         zhigh = 1;
      } else {
         min0 = 2;
         zhigh = 2;
         NTL::conv(max0, trunc(sqrt((m_lMin2 - m_n2[j]) / m_dc2[j])));
      }
   }

   Real temp;
   /* On essaie maintenant chacun des z[j] dans l'intervalle, en      */
   /* commencant par le centre puis en alternant d'un cote a l'autre. */
   while (zlow >= min0 || zhigh <= max0) {

      if (high) {
         m_zLI[j] = zhigh;
      } else {
         m_zLI[j] = zlow;
      }
      m_zLR[j] = m_zLI[j];

      // Calcul de m_n2[j-1].
      x = m_zLR[j] - center;

      if (j == 0) {
         Real tmps_n2 = m_n2[0] + x * x * m_dc2[0];
         if (tmps_n2 < m_lMin2) {
            // On verifie si on a vraiment trouve un vecteur plus court
            // NTL::Mat_row<const Int> row1(m_lat->getBasis(), dim - 1);
            IntVec &row1 = m_lat->getBasis()[dim-1];
            m_bv = row1;
            for (k = 0; k < dim - 1; k++) {
               if (m_zLI[k] != 0) {
                  IntVec &row1 = m_lat->getBasis()[k];
                  //NTL::Mat_row<Int> row1(m_lat->getBasis(), k);
                  ModifVect(m_bv, row1, m_zLI[k], dim);
               }
            }
            if (Stage == 3) {
               IntVec &row1 = m_lat->getBasis()[dim-1];
               //NTL::Mat_row<Int> row1(m_lat->getBasis(), dim - 1);
               ModifVect(m_bv, row1, m_zLR[dim - 1] - 1.0, dim);
            }

            ProdScal<Int>(m_bv, m_bv, dim, S1);
            NTL::conv(S4, m_lat->getVecNorm(dim - 1));
            if (S1 < S4) {
               if (Stage == 2) {
                  smaller = true;
                  if (!PreRedLLLMink) m_zShort = m_zLI;
                  else {
                     for (k = 1; k < dim; k++) {
                        // NTL::Mat_row<Int> row1(WTemp, k);
                        ProdScal<Int>(m_bv, WTemp[k], dim, S2);
                        Quotient(S2, mR, S3);
                        NTL::conv(m_zShort[k], S3);
                     }
                     m_zShort[dim - 1] = 1;
                  }
               } else if (Stage == 3 && !PreRedLLLMink) {
                  if (GCD2vect(m_zLI, i, dim) == 1) {
                     m_zShort = m_zLI;
                     smaller = true;
                  }
               } else {
                  for (k = 0; k < dim; k++) {
                     // NTL::Mat_row<Int> const row1(WTemp, k);
                     ProdScal<Int>(m_bv, WTemp[k], dim, S2);
                     Quotient(S2, mR, S3);
                     NTL::conv(m_zShort[k], S3);
                  }
                  if (GCD2vect(m_zShort, i, dim) == 1) {
                     smaller = true;
                  }
               }
               if (smaller) {
                  NTL::conv(temp, S1);
                  m_lat->setVecNorm(temp, dim - 1);
                  return true;
               }
            }
         }
      } else { // j > 0
         m_n2[j - 1] = m_n2[j] + x * x * m_dc2[j];
         if (m_lMin2 >= m_n2[j - 1]) {
            if (!tryZMink(j - 1, i, Stage, smaller, WTemp)) return false;
            // Des qu'on a trouve quelque chose, on sort de la recursion
            // et on retourne dans reductMinkowski.
            if (smaller) return true;
         }
      }
      if (high) {
         ++zhigh;
         if (zlow >= min0) high = false;
      } else {
         --zlow;
         if (zhigh <= max0) high = true;
      }
   }
   return true;
}

//=========================================================================

template<typename Int, typename Real>
bool ReducerBB<Int, Real>::redBBMink(int64_t i, int64_t d, int64_t Stage, bool &smaller,
      bool taboo[]) {
   /*
    * Tries to shorten m_lat->getBasis()[i] using branch-and-bound.
    * Used in Minkowski Reduction.
    * Stage is 2 or 3.
    * z[i] = 1 if Stage = 2, z[i] >= 2 if Stage = 3.
    * Stops and returns false if not finished after examining MaxNodesBB
    * nodes in the branch-and-bound tree.  When succeeds, returns true.
    * Assumes that the norm is Euclidean.
    */
   // bool withDual = false;  // m_lat->withDual();
   const int64_t dim = m_lat->getDim();
   IntMat VTemp(NTL::INIT_SIZE, dim, dim);
   IntMat WTemp(NTL::INIT_SIZE, dim, dim);
   bool TabooTemp[dim];
   Real tmp;
   smaller = false;

   // Approximation du carre de la longueur de Vi.
   if (m_lat->getVecNorm(i) < 0) {
      IntVec &row1 = m_lat->getBasis()[i];
      ProdScal<Int>(row1, row1, dim, tmp);
      //  ProdScal<Int> (m_lat->getBasis()[i], m_lat->getBasis()[i],
      //            dim, tmp);
      m_lat->setVecNorm(tmp, i);
   }
   NTL::conv(m_lMin2, m_lat->getVecNorm(i));
   m_lat->updateVecNorm();
   m_lat->permute(i, dim - 1);

   int64_t h;
   if (PreRedLLLMink) {
      // On memorise la base courante.
      VTemp = m_lat->getBasis();
      for (h = 0; h < dim; h++)
         TabooTemp[h] = true;
      redLLL<Int, Real>(m_lat->getBasis(), 0.99999, dim);
      // redLLLOld(1.0, 1000000, dim - 1);
      m_lat->updateVecNorm();
   }
   if (!calculCholesky(m_dc2, m_c0)) return false;
   m_countNodes = 0;
   m_n2[dim - 1] = 0.0;
   if (!tryZMink(dim - 1, i, Stage, smaller, WTemp)) return false;

   if (PreRedLLLMink) {
      /* On remet l'anciennne base, celle d'avant LLL, avant de considerer
       la prochaine m_lat->dimension.  */
      m_lat->getBasis() = VTemp;
      m_lat->updateVecNorm();
      for (h = 0; h < dim; h++)
         taboo[h] = TabooTemp[h];
   }
   if (smaller) {
      /* On a trouve un plus court vecteur.  On ameliore
       m_lat->getBasis()[k].  */
      int64_t k = 0;
      if (Stage == 2) k = dim - 1;
      else insertBasisVector(m_zShort);
      for (h = 0; h < dim; h++) {
         if ((m_zShort[h] != 0) && (h != k)) {
            if (Stage == 2) {
               if (h >= d) taboo[h] = false;
            }
         }
      }
   } else if (Stage == 2) taboo[dim - 1] = true;

   m_lat->permute(i, dim - 1);
   // trace( "APRES redBBMink");
   return true;
}

//=========================================================================

template<typename Int, typename Real>
bool ReducerBB<Int, Real>::tryZShortVec(int64_t j, bool &smaller, NormType norm) {
   /*
    * This recursive procedure uses a lower-triangular matrix obtained either via the Cholesky
    * decomposition or by taking a lower-triangular basis. It works for the
    * L1NORM and L2NORM.  It is called initially with j=dim-1 by redBBShortVec
    * to find a shortest vector via the BB procedure.
    * If `m_countNodes > MaxNodesBB`, it returns `false`, otherwise it returns `true`,
    * meaning that a shorter vector has been found and placed in `m_bv`,
    * its square length was put in `m_lMin2`, and the vector of the `z_j`'s in in `m_zShort`.
    * The basis itself is not modified, this is done afterwards in `insertBasisVector`.
    */

   // For a non-recursive implementation, these variables should be arrays indexed by j.
   Real dc, x, center, mn_xsquare_md;
   std::int64_t min0, max0;     // Bornes de l'intervalle pour les z_j.
   std::int64_t zlow, zhigh; // Valeur courante a gauche et a droite du centre.
   bool high; // Indicates if we are on the right (true) or the left of the center.
   int64_t k;
   std::int64_t temp;
   const int64_t dim = m_lat->getDim();
   ++m_countNodes;
   if (m_countNodes > maxNodesBB) {
      std::cerr << "*****   m_countNodes > maxNodesBB = " << maxNodesBB << std::endl;
      return false;
   }

   /* Compute an interval that contains the admissible values of zj. */
   /* This computation is for the L2 norm, but also works for the L1 norm. */
   /* 1. Compute the center of the interval.  */
   center = 0.0;
   dc = 0.0;
   if (m_decomp == CHOLESKY) {
      for (k = j + 1; k < dim; ++k)
         center -= m_c0[j][k] * m_zLR[k];
      // This dc is the distance from the center to the boundaries.
      // m_lMin2 contains the square length of current shortest vector with the selected norm.
      dc = sqrt((m_lMin2 - m_n2[j]) / m_dc2[j]);
   }
   if (m_decomp == TRIANGULAR && norm == L2NORM) {
      for (k = 0; k < j; ++k)
         center -= m_c0[j][k] * m_zLR[k];
      dc = sqrt((m_lMin2 - m_n2[j]) / m_dc2[j]);
   }
   if (m_decomp == TRIANGULAR && norm == L1NORM) {
      for (k = j + 1; k < dim; ++k)
         center -= m_c0[k][j] * m_zLR[k];
      dc = (m_lMin - m_n2[j]) / m_c0[j][j];
   }
   // Compute two integers min0 and max0 that are the min and max integers in the interval.
   if (!m_foundZero) min0 = 0;     // We are at the beginning, min will be zero.
   else {
      x = center - dc;
      NTL::conv(min0, trunc(x));
      if (x > 0.0) ++min0;
   }
   x = center + dc;
   NTL::conv(max0, trunc(x));
   if (x < 0.0) --max0;

   //   std::cout << "Borne Min="<<min0<<"     Borne Max="<<max0<<"     j="<<j<<std::endl;
   // Compute initial values of zlow and zhigh, the search pointers on each side of the interval.
   if (min0 > max0) return true;
   if (min0 == max0) {
      zlow = min0;
      zhigh = max0 + 1;
      high = false;
   } else if (center >= 0.0) {
      NTL::conv(zlow, trunc(center));
      zhigh = zlow + 1;
      NTL::conv(temp, trunc(2.0 * center));
      high = (temp & 1);
   } else {
      NTL::conv(zhigh, trunc(center));
      zlow = zhigh - 1;
      NTL::conv(temp, -trunc(2.0 * center));
      high = (temp & 1) == 0;
   }

   // We try each zj in the interval, starting in the center and alternating between left and right.
   while (zlow >= min0 || zhigh <= max0) {
      if (high) m_zLI[j] = zhigh;
      else m_zLI[j] = zlow;
      m_zLR[j] = m_zLI[j];

      // Computing m_n2[j-1].
      x = m_zLR[j] - center;
      mn_xsquare_md = m_n2[j] + x * x * m_dc2[j]; // We pre-compute this to save time.
      if (j == 0) {
         // All the zj have been selected: we now have a candidate vector to test!
         if (m_lMin2 > mn_xsquare_md) {
            /* Check if we have a shorter nonzero vector. */
            if (!m_foundZero) {
               // The first vector found will always be zero, we discard it.
               m_foundZero = true;
            } else {
               SetZero(m_bv, dim);
               for (k = 0; k < dim; k++) {
                  if (m_zLI[k] != 0) {
                     IntVec &row1 = m_lat->getBasis()[k];
                     // NTL::Mat_row<Int> row1(m_lat->getBasis(), k);
                     ModifVect(m_bv, row1, m_zLI[k], dim);
                  }
               }
               // The new shortest vector is now in `m_bv`.
               if (m_lat->getNormType() == L2NORM) {
                  ProdScal<Int>(m_bv, m_bv, dim, x);
               } else {
                  // Compute the square length for the L1 norm.
                  CalcNorm<IntVec, Real>(m_bv, dim, x, m_lat->getNormType());
                  x = x * x;
               }
               if (x < m_lMin2) {
                  // The new vector is shorter.
                  smaller = true;
                  NTL::conv(m_lMin2, x);
                  m_zShort = m_zLI;
                  // m_bw = m_bv;
               }
            }
         }
      } else if (m_lMin2 > mn_xsquare_md) {
         // There is still hope; we continue the recursion.
         m_n2[j - 1] = mn_xsquare_md;
         if (!tryZShortVec(j - 1, smaller, norm)) return false;
      } else m_n2[j - 1] = mn_xsquare_md;
      if (high) {
         ++zhigh;
         if (zlow >= min0) high = false;
      } else {
         --zlow;
         if (zhigh <= max0) high = true;
      }
   }
   return true;
}

//=========================================================================

template<typename Int, typename Real>
bool ReducerBB<Int, Real>::redBBShortVec() {
   /*
    * Finds shortest non-zero vector, using branch-and-bound, with L1 or L2 norm.
    *   ****  It may work only for the L2 norm, for now.  ****
    * Stops and returns false if not finished after examining MaxNodesBB nodes in the
    * branch-and-bound tree.  When succeeds, returns true, and the squared shortest
    * vector length will be in m_lMin2.
    *
    * This function uses (directly or indirectly) the following class variables:
    *    m_lMin, m_lMin2, m_decomp, m_boundL2, m_n2, m_countNodes, m_foundZero,
    *    m_bv, m_zShort, m_c0, m_zLR, m_zLI, m_dc2, ....  and more.
    * From m_lat (current lattice object):
    *    norm, sortBasisNoDual, updateScalL2Norm, getBasis,
    */
   // std::cout << " redBBShortVec, entering \n";
   NormType norm = m_lat->getNormType();
   if ((norm != L1NORM) & (norm != L2NORM)) {
      std::cerr << "RedBBShortVec: only L1 and L2 norms are supported";
      return false;
   }
   const int64_t dim = m_lat->getDim();  // Lattice dimension
   //std::cout << " Start redBBShortVec, dim  = " << dim << "\n";
   //std::cout << " Start redBBShortVec, basis = \n" << m_lat->getBasis() << "\n";

   bool smaller = false;  // Will change when we find a smaller vector.
   int64_t k;
   Real x(0.0);
   // if (norm != L2NORM)  m_lat->setNegativeNorm();
   // Here we sort the basis by L2 lengths, otherwise Cholesky will fail more rapidly
   // due to floating-point errors.
   m_lat->updateScalL2Norm(0, dim);
   m_lat->sortBasis(0);
   // std::cout << " redBBShortVec, check norm, norm = L_" << norm << "\n";

   // Approximate the square norm of the current shortest vector.
   if (norm == L2NORM) {
      NTL::conv(m_lMin2, m_lat->getVecNorm(0));
   } else {
      // Looking for the shortest vector in basis according to the considered norm
      // NTL::Mat_row<Int> row1(m_lat->getBasis(), 0);
      IntVec &row1 = m_lat->getBasis()[0];
      CalcNorm<IntVec, Real>(row1, dim, m_lMin, norm);
      for (k = 1; k < dim; k++) {
         // NTL::Mat_row<Int> row2(m_lat->getBasis(), k);
         IntVec &row2 = m_lat->getBasis()[k];
         CalcNorm<IntVec, Real>(row2, dim, x, norm);
         if (x < m_lMin) m_lMin = x;
      }
      m_lMin2 = m_lMin * m_lMin;  // Squared shortest length with L1 norm.
   }
   //std::cout << "\n\n ==============================================\n";
   //std::cout << "Entering redBBShortVec, square length of shortest = " <<  m_lMin2 << "\n";
   //std::cout << " Initial basis after LLL = \n" << m_lat->getBasis() << "\n";

   // If we already have a shorter vector than the minimum threshold, we stop right away.
   // This is useful for the seek programs in LatMRG.
   if (m_lMin2 <= m_BoundL2[dim - 1]) return false;

   if (m_decomp == CHOLESKY) {
      /* Perform the Cholesky decomposition; if it fails we exit. */
      if (!calculCholesky(m_dc2, m_c0)) return false;
   } else if (m_decomp == TRIANGULAR) {  // Just for testing; this is very slow!
      /* Perform a triangular decomposition:
       * Instead of doing the following, perhaps we may just assume that the basis
       * is already lower triangular?
       */
      IntMat m_v, m_v2;   // Here we create new matrices each time!!!!
      m_v.SetDims(dim, dim);
      m_v2.SetDims(dim, dim);
      Int mod = m_lat->getModulus();
      CopyMatr(m_v, m_lat->getBasis(), dim, dim);
      lowerTriangularBasis(m_v, m_v2, mod);
      // CopyMatr(m_lat->getBasis(), m_v2, dim, dim);
      for (int64_t i = 0; i < dim; i++) {
         for (int64_t j = 0; j < dim; j++) {
            if (i != j) m_c0[i][j] = NTL::conv < Real > (m_v2[i][j]) / NTL::conv < Real
                  > (m_v2[i][i]);
            else m_c0[i][j] = NTL::conv < Real > (m_v2[i][j]);
         }
      }
      for (int64_t i = 0; i < dim; i++) {
         m_dc2[i] = m_c0[i][i] * m_c0[i][i];
      }
   } else {
      std::cerr << "RedBBShortVec:decomp value not supported";
      return false;
   }

   // std::cout << " redBBShortVec, after Cholesky decomp \n";
   /* Perform the branch and bound.  */
   /* m_n2[j] will be the sum of terms |z*k|^2 ||v*k||^2 for k > j.  */
   m_n2[dim - 1] = 0.0;
   m_countNodes = 0;
   smaller = false;
   m_foundZero = false;
   // std::cout << " redBBShortVec, basis before tryZ = \n" << m_lat->getBasis() << "\n";
   if (!tryZShortVec(dim - 1, smaller, norm)) // We search for a shortest vector.
      return false;
   // std::cout << " redBBShortVec, after tryShortVec \n";
   if (smaller) {
      //std::cout << " redBBShortVec, found a shorter vector, square length: " << m_lMin2 << "\n";
      //std::cout << " redBBShortVec, m_bv = " << m_bv << "\n";
      //std::cout << " redBBShortVec, before transform Stage 3, basis = \n" << m_lat->getBasis() << "\n";

      // We found a shorter vector. It is in m_bv and its square length is in m_lMin2.
      // The short vector m_bv is not yet put in the basis.
      // The following does that and is useful only if we want to continue working with this basis.

      insertBasisVector(m_zShort); // Is this useful and OK for L1 ???
      // insertBasisVectorLLL(m_zShort); // Is this OK for L1.  This one is just a bit slower.
      //std::cout << " redBBShortVec, after transform Stage 3, basis = \n" << m_lat->getBasis() << "\n";
      //std::cout << " redBBShortVec, m_bv = " << m_bv << "\n";
      // std::cout << " redBBShortVec, m_bv[0] = " << m_bv[0] << "\n";
      // std::cout << " redBBShortVec, basis after update norm = \n" << m_lat->getBasis() << "\n";

      // After the following, the new current shortest vector will be in m_lat->getBasis()(0).
      // In the case of L1NORM, we must check if it is really smaller.
/*
      if (norm == L2NORM) {
         m_lat->permute(k, 0);
      }
      else {
         if (x < m_lMin) {
            m_lMin = x;
            m_lMin2 = m_lMin * m_lMin;
            m_lat->permute(k, 0);
         }
      }
*/
      // std::cout << " redBBShortVec, after permute, basis = \n" << m_lat->getBasis() << "\n";
   }
   // m_lat->updateVecNorm();
   // m_lat->sortBasis(0);
   //std::cout << " redBBShortVec, basis at exit= \n" << m_lat->getBasis() << "\n";
   //std::cout << "Exiting redBBShortVec, square length of shortest = " <<  m_lMin2 << "\n";
   return true;
}

//=========================================================================

// This works only for the L2 norm.
template<typename Int, typename Real>
bool ReducerBB<Int, Real>::reductMinkowski(int64_t d) {
   // bool withDual = false;  //  m_lat->withDual();
   const int64_t dim = m_lat->getDim();
   int64_t i;
   std::int64_t totalNodes = 0;
   bool found;
   bool smaller;      // A smaller vector has been found
   bool taboo[dim];   // taboo[i]=true means that vector i should not be modified.

   do {
      // The first d vectors should not be modified.
      for (i = 0; i < d; i++)
         taboo[i] = true;
      for (i = d; i < dim; i++)
         taboo[i] = false;
      found = false;
      do {
         // redDieter(d, taboo);  // Not very efficient...  Change to BKZ.
         m_lat->setNegativeNorm(d);
         m_lat->updateVecNorm(d);
         m_lat->sortBasis(d);
         found = false;
         for (i = 0; i < dim; i++) {
            if (!taboo[i]) {
               // On essaie de reduire le i-eme vecteur.
               if (!redBBMink(i, d, 2, smaller, taboo)) return false;
               totalNodes += m_countNodes;
               if (smaller) {
                  found = true;
               }
            }
         }
      } while (found);
      // Stage 3
      if (dim > 7) {
         for (i = d; i < dim; i++) {
            if (!redBBMink(i, d, 3, smaller, taboo)) return false;
            totalNodes += m_countNodes;
            if (smaller) found = true;
         }
      }
   } while (found);
   m_lat->setNegativeNorm();
   m_lat->updateScalL2Norm(0, dim);
   m_lMin2 = m_lat->getVecNorm(0);
   return true;
}

//=========================================================================

template<typename Int, typename Real>
void ReducerBB<Int, Real>::tracePrintBases(char *message) {
   std::cout << std::endl << "================================= " << message << std::endl;
   //std::cout << "dim = " << m_lat->getDim () << std::endl;
   m_lat->setNegativeNorm();
   m_lat->updateVecNorm();
   m_lat->sortBasis(0);
   m_lat->write();
}

//=========================================================================

template<typename Int, typename Real>
bool ReducerBB<Int, Real>::reductMinkowski(IntLattice<Int, Real> &lat, int64_t d) {
   setIntLattice(lat);
   return ReducerBB<Int, Real>::reductMinkowski(d);
}

//=========================================================================

template<typename Int, typename Real>
bool ReducerBB<Int, Real>::shortestVector() {
   return ReducerBB<Int, Real>::redBBShortVec();
}

template<typename Int, typename Real>
bool ReducerBB<Int, Real>::shortestVector(IntLattice<Int, Real> &lat) {
   setIntLattice(lat);
   return ReducerBB<Int, Real>::redBBShortVec();
}

//============================================================================

template class ReducerBB<std::int64_t, double> ;
template class ReducerBB<NTL::ZZ, double> ;
template class ReducerBB<NTL::ZZ, xdouble> ;
template class ReducerBB<NTL::ZZ, quad_float> ;
template class ReducerBB<NTL::ZZ, NTL::RR> ;

}     // namespace LatticeTester

#endif 
