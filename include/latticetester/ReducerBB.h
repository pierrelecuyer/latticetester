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
 * \file latticetester/ReducerBB.h
 *
 * This `ReducerBB` class provides functions to find a shortest nonzero vector in the lattice
 * using a BB algorithm as in \cite mFIN85a,
 * and to compute a Minkowski basis reduction as in \cite rAFF85a.
 *
 * Each `Reducer` must have an internal `IntLattice` object which is given upon construction
 * and can be changed via `setIntLattice`.
 * The `shortestVector` and `reductMinkowski` functions are applied to this internal object
 * and will use the norm associated with that `IntLattice` object.
 * These functions do not apply any pre-reduction by themselves.
 * Before calling them, one should always pre-reduce the basis via LLL or BKZ,
 * because it drastically reduces the size of the BB search tree.
 *
 * The `Reducer` object maintains several internal variables, vectors, and matrices
 * used by the `shortestVector` and `reductMinkowski` functions.
 * It is recommended to create a single `Reducer` object with a large enough maximal dimension
 * and then call the main functions with the relevant `IntLattice` object as a parameter.
 * One may also re-use the same `IntLattice` objects for many different lattices,
 * for example when performing a search for a good lattice.
 * The norm type, dimension, basis, vector lengths, etc.
 * will be taken from this internal `IntLattice` object.  The dimensions of the internal vectors
 * and matrices can be larger than required; the methods will use only the
 * entries that are needed for the given `IntLattice` basis.
 * Creating a new `Reducer` object for each `IntLattice` that we want to handle is very
 * inefficient and should be avoided, especially when we want to examine several
 * projections for several lattices.
 *
 * The functions in this class do not use or change the m-dual lattice.
 * To find a shortest nonzero vector in the m-dual lattice,
 * one should dualize the lattice (`IntLattice::dualize` does that)
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

   /*
    * Copy constructor.
    */
   //ReducerBB(const ReducerBB<Int, Real> &red);

   /**
    * Destructor.
    */
   ~ReducerBB();

   /*
    * Assignment operator that makes a deep copy of `red`
    * into the current object, using `copy`.
    */
   // ReducerBB<Int, Real>& operator=(const ReducerBB<Int, Real> &red);

   /*
    * Copies `red` into the current object.
    */
   //void copy(const ReducerBB<Int, Real> &red);

   /**
    * Computes a shortest non-zero vector for the `IntLattice` stored in this `ReducerBB` object,
    * for the norm stored in this `IntLattice`,
    * using the BB algorithm described in \cite rLEC97c and \cite iLEC22l.
    * The admissible norm types are `L1NORM` and `L2NORM` (see `EnumTypes.h`).
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
    * In particular, the bounds set by `setBounds2` have to be reset.
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
    * Returns the *square length* of the current shortest basis vector in the lattice,
    * usually computed with the current norm.
    */
   Real getMinLength2() {
      return m_lMin2;
   }

   /**
    * Returns the *length* (not squared) of the current shortest basis vector in the lattice,
    * which is stored in a local variable.
    * This depends only on the lattice, but this length is stored and used in this class.
    */
   Real getMinLength() {
      if (m_lat->getNormType() == L2NORM) return sqrt(m_lMin2);
      else return m_lMin1;
   }

   /**
    * Returns the length of the current *last* basis vector in the lattice.
    */
   Real getMaxLength() {
      if (m_lat->getNormType() == L2NORM) return sqrt(m_lat->getVecNorm(m_lat->getDim() - 1));
      else return m_lat->getVecNorm(m_lat->getDim() - 1);
   }

   /**
    * Returns the number of nodes that have been visited in the BB tree.
    */
    int64_t getCountNodes() {
      return m_countNodes;
   }

   /**
    * Returns the number of leaves visited in the BB tree.
    */
    int64_t getCountLeaves() {
       return m_countLeaves;
   }

   /**
    * Returns the maximal absolute value of a `z_j` in the latest BB.
    */
    int64_t getMaxZj() {
       return m_maxZj;
    }

   /**
    * Sets a vector of bounds on the square of the acceptable shortest
    * vector lengths, for each dimension from `dim1+1` to `dim2`.
    * For each `i`, `thresholds[i]` must contain a lower bound on the square length
    * of the shortest lattice vector in dimension `i+1`, for the current norm.
    * This bound will be used during during the Branch-and-Bound step
    * when computing a shortest lattice vector.  As soon as a vector shorter
    * than the bound is found, the BB algorithm will stop. This is useful to reduce
    * work when we search for a good lattice with the spectral test.
    * If these bounds are not set, the default values of 0 are used.
    * It is recommended to set these bounds before calling `shortestVector`
    * for the first time when making searches for good lattices.
    */
   void setBounds2(const RealVec &thresholds, int64_t dim1, int64_t dim2);

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
    * in the BB procedures. It can be `CHOLESKY` or `TRIANGULAR`.
    */
   void setDecompTypeBB(DecompTypeBB decomp) {
      m_decomp = decomp;
   }

   /**
    * Sets the level of verbosity in the terminal output.
    * The default value is 0 (minimal output).
    * Values from 1 to 4 give increasingly more details.
    */
   void setVerbosity(int64_t verbose) {
      m_verbose = verbose;
   }

   /**
    * The maximum number of nodes in the branch-and-bound tree when
    * calling `shortestVector` or `reductMinkowski`. When this number is
    * exceeded, the method aborts and returns `false`.
    */
   int64_t maxNodesBB = 10000000;


private:

   /**
    * Initializes all matrices and vectors used in this object.
    */
   void init(int64_t maxDim);

   /**
    * Recursive procedure that tries to find a shorter vector using BB.
    * It is called by `redBBShortVec`. The parameter j indicates the level of the
    * BB tree. The first call is for level i=0.
    */
   bool tryZShortVec(int64_t j, bool &smaller, NormType norm);

   bool tryZShortVecOld(int64_t j, bool &smaller, NormType norm);

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
    * Computes the LDL Cholesky decomposition of the basis. Returns in `m_L` the
    * lower-triangular matrix of the Cholesky decomposition that are below the diagonal.
    * Returns in `m_dc2` the squared elements of the diagonal.
    * All these elements are in `Real`.
    * The upper-triangular part of `m_L` (including the diagonal) is not initialized.
    */
   bool calculCholeskyLDL();

   /**
    * Computes a lower-triangular basis `L` with elements `ell_{i,j}`.
    * Then put in `m_L` the elements `\tilde\ell_{i,j} = \ell_{i,j}/\ell_{i,i}`
    * below the diagonal, the elements `\ell_{j,j}` on the diagonal,
    * and puts in `m_dc2` the squared elements of the diagonal of `L`.
    * All these elements are in `Real`.
    */
   bool calculTriangularL();

   /**
    * In this function, we assume that we have found a new shorter vector
    * \f$ \bv = \sum_{i=1}^t z_i \bv_i\f$ and we want to insert it in the basis.
    * If \f$z_j = \pm 1\f$ for some \f$j\f$, we can simply exchange \f$\bv\f$ with  \f$\bv_j\f$
    * in the basis. Otherwise, we transform the basis and the \f$z_j\f$'s so that one of the
    * nonzero \f$z_j\f$'s becomes equal to  \f$\pm 1\f$, as explained in the guide, and then
    * we make the exchange. After that, we exchange \f$\bv\f$ with the vector in first position.
    * Note that the basis transformation can make some of the old basis vectors very large.
    * This function works even if v is not shorter than the previous basis vectors.
    * It does not look at the vector norms, or at which norm we use.
    * It implements (roughly) the ``transformation of stage 3'' described in \cite rAFF85a.
    */
   void insertBasisVector(std::vector<std::int64_t> &z);

   /**
    * This function provides an alternative to `insertBasisVector` when we use the \f$L^2\f$ norm
    * and we are sure that \f$\bv\f$ is actually a shortest vector.
    * It adds the new vector \f$\bv = \sum_{i=1}^t z_i \bv_i\f$ to the basis to form a set
    * of `dim+1` generating vectors, and it applies LLL to recover a new basis that contains \f$\bv\f$.
    * Note that if \f$\bv\f$ is not a shortest vector for the \f$L^2\f$ norm, it may be removed
    * from the basis by LLL.
    */
   void insertBasisVectorLLL(std::vector<std::int64_t> &z);

   /*
    * Debug function that sorts and prints the primal and dual bases
    * to standard output, using the `write` function.
    */
   // void tracePrintBasis(char *mess);

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
   RealVec m_Bounds2;

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
    * Private working variables for this class.
    * They are used inside the basis reduction and short vector methods, and
    * are declared here to avoid passing them as parameters across the methods.
    * The matrices and vectors are sized to some maximum dimensions in init(),
    * which must be large enough for all computations handled by this ReducerBB object.
    */
   int64_t m_maxDim; // maximum dimension.
   int64_t m_dim; // Current dimension.
   Real m_lMin1;  // The norm of the shortest vector in the primal basis
                  // according to the norm considered
   Real m_lMin2;  // Squared L2-norm of the shortest vector in the primal basis.
   IntVec m_bv;   // Saves current shortest vector in primal basis
   IntVec m_vtest;  // A new shortest vector candidate to be tested.

   // Vectors used in BB algorithms, in tryZShortVec and tryZMink.
   // m_sjp[j] will be the sum of terms |z*k|^2 ||v*k||^2 for k > j, updated in `tryZShortVec`.
   // Exception: for triangular method, it is the sum of |z*k|^p ||v*k||^p.
   // m_dc2[j] corresponds to d_j in the description of BB in the user guide.
   // For Cholesky, these are the elements of the diagonal matrix D.
   // These vectors are updated only one coordinate at a time in the BB tree.
   RealVec m_sjp, m_dc2;

   // Matrices used in the Cholesky and triangular decomposition, for the BB bounds.
   // For Cholesky, m_L represents `tilde L` from the guide and `m_c2` is used internally.
   // For Triangular, m_L contains the lower-triangular basis.
   // We try to avoid resizing them!
   RealMat m_L, m_c2;

   std::vector<std::int64_t> m_z;   // Vector of (integer) values of z_i.   *************  TOO SMALL!
   RealVec m_zLR;   // Same vector in floating point (Real), needed when calculating bounds.
   std::vector<std::int64_t> m_zShort;  // Values of z_i for shortest vector.
   int64_t m_maxZj = 0;    // Largest absolute value of a `z_j` in latest BB.

   int64_t m_countNodes = 0;  // Number of visited nodes in the BB tree
   int64_t m_countLeaves = 0;  // Number of visited leaves in the BB tree
   bool m_foundZero;    // = true -> the zero vector has been handled

   int64_t m_verbose = 0;  // Indicates how much details are printed.
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
   m_L.SetDims(dim1, dim1);
   m_c2.SetDims(dim1, dim1);
   // Indices in m_IC can go as high as maxDim + 2.
   // m_IC = new int64_t[2 + dim1];

   m_vtest.SetLength(dim1);
   m_bv.SetLength(dim1);
   m_sjp.SetLength(dim1);
   m_zLR.SetLength(dim1);
   m_z.resize(dim1);
   m_zShort.resize(dim1);
   m_dc2.SetLength(dim1);
   m_Bounds2.SetLength(dim1);

   m_lMin1 = std::numeric_limits<double>::max();
   m_lMin2 = m_lMin1;
   for (int64_t i = 0; i < dim1; i++) {
      m_z[i] = -1;
      m_zShort[i] = -1;
      m_Bounds2[i] = -1;
   }
   m_countNodes = 0;
   m_countLeaves = 0;
   m_maxZj = 0;
   m_foundZero = false;
   PreRedLLLMink = false;
   maxNodesBB = 1000000000;
}

//=========================================================================

/*
template<typename Int, typename Real>
ReducerBB<Int, Real>::ReducerBB(const ReducerBB<Int, Real> &red) {
   copy(red);
}
*/

//=========================================================================

/*
template<typename Int, typename Real>
ReducerBB<Int, Real>&
ReducerBB<Int, Real>::operator=(const ReducerBB<Int, Real> &red) {
   if (this != &red) copy(red);
   return *this;
}
*/

//=========================================================================

/*
template<typename Int, typename Real>
void ReducerBB<Int, Real>::copy(const ReducerBB<Int, Real> &red) {
   m_lat = red.m_lat;
   m_maxDim = red.m_maxDim;
   m_L = red.m_L;
   m_c2 = red.m_c2;
   m_dc2 = red.m_dc2;
   m_vtest = red.m_vtest;
   m_bv = red.m_bv;
   // m_bw = red.m_bw;
   m_sjp = red.m_sjp;
   m_zLR = red.m_zLR;
   m_z = red.m_z;
   m_zShort = red.m_zShort;
   // m_cho2 = red.m_cho2;
   // m_gramVD = red.m_gramVD;
   m_lMin1 = red.m_lMin1;
   m_lMin2 = red.m_lMin2;
   m_Bounds2 = red.m_Bounds2;
   m_countNodes = 0;
   m_countLeaves = 0;
   m_foundZero = false;
   PreRedLLLMink = false;
   maxNodesBB = 1000000000;
}
*/

//=========================================================================

template<typename Int, typename Real>
ReducerBB<Int, Real>::~ReducerBB() {
   m_L.kill();
   m_c2.kill();
   m_bv.kill();
   m_vtest.kill();
   m_sjp.kill();
   m_zLR.kill();
   m_z.resize(0);
   m_zShort.resize(0);
   m_dc2.kill();
   m_Bounds2.kill();
}

//=========================================================================

template<typename Int, typename Real>
void ReducerBB<Int, Real>::setBounds2(const RealVec &thresholds, int64_t dim1, int64_t dim2) {
   m_Bounds2.SetLength(dim2);
   for (int64_t i = dim1; i < dim2; i++)
      m_Bounds2[i] = thresholds[i];
}

//=========================================================================

template<typename Int, typename Real>
bool ReducerBB<Int, Real>::calculCholeskyLDL() {
   const int64_t dim = m_lat->getDim();
   // Compute first dim rows of lower-triangular part of L.
   for (int64_t j = 0; j < dim; j++) {  // Column j.
      for (int64_t i = j; i < dim; i++) {  // Row i >= j.
         // m_c2[i][j] will first contain a_{i,j}, and finally `d_j \tilde\ell_{i,j}`
         // m_L[i][j] will contain `\tilde\ell_{i,j}`
         ProdScal<Int>(m_lat->getBasis()[i], m_lat->getBasis()[j], dim, m_c2[i][j]);
         for (int64_t k = 0; k < j; k++)
            m_c2[i][j] -= m_L[j][k] * m_c2[i][k];
         if (i == j) {
            // m_L[i][i] = sqrt(m_c2[i][i]);  // Not really needed, just for printing m_L.
            m_dc2[i] = m_c2[i][i];
            if (m_dc2[i] < 0.0) {
               std::cout << "\n*** Negative diag. element in Cholesky Decomp.\n" << std::endl;
               return false;
            }
         } else {
            m_L[i][j] = m_c2[i][j] / m_dc2[j];
            m_L[j][i] = 0;   // Just for when we print this matrix.
         }
      }
   }
   return true;
}

//=========================================================================

template<typename Int, typename Real>
bool ReducerBB<Int, Real>::calculTriangularL() {
   // If the basis is already lower triangular, this will be fast.
   const int64_t dim = m_lat->getDim();
   // IntMat &basis = m_lat->getBasis();
   IntMat copybasis;  // Here we create two new matrices each time!
   IntMat &tribasis =  m_lat->getBasis();  // The current basis, will become triangular.
   copybasis.SetDims(dim, dim); // A copy of the current basis.
   // tribasis.SetDims(dim, dim);  // Will be a lower-triangular basis. We need it
   // Int mod = m_lat->getModulus();
   CopyPartMat<IntMat>(copybasis, m_lat->getBasis(), dim, dim);  // Copy current basis into `copybasis`.
   lowerTriangularBasis(tribasis, copybasis, m_lat->getModulus());  // Here `copybasis` may be destroyed.
   if (m_verbose > 2)
      std::cout << " triangularL, lower triangular basis = \n" << m_lat->getBasis() << "\n";
   for (int64_t i = 0; i < dim; i++) {
      for (int64_t j = 0; j < dim; j++) {
         if (i != j) m_L[i][j] = NTL::conv <Real> (tribasis[i][j]) / NTL::conv<Real> (tribasis[j][j]);
         else m_L[j][j] = NTL::conv <Real> (tribasis[j][j]);
      }
   }
   // std::cout << " triangularL, lower triangular basis m_L = \n" << m_L << "\n";
   for (int64_t i = 0; i < dim; i++) {
      m_dc2[i] = m_L[i][i] * m_L[i][i];  // These are the square diagonal elements.
   }
   return true;
}

// =========================================================================

// This function puts the new vector v = m_bv in the basis without looking at its length.
// It does not have to be shorter than the current basis vectors.
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
      // We have 2 indices i < j such that |z_j| > 1 and z_i != 0.
      while (z[j]) {
         // Truncate the quotient towards 0.
         q = z[i] / z[j];
         if (q) {
            // We add q * v[i] to v[j]
            z[i] -= q * z[j];
            IntVec &row1 = m_lat->getBasis()[j];
            IntVec &row2 = m_lat->getBasis()[i];
            ModifVect(row1, row2, q, dim);
            // ModifVectModulo(row1, row2, q, m_lat->getModulus(), dim);  // In Util.h
            m_lat->setNegativeNorm(j);
         }
         // Permute the two vectors v[i] and v[j].
         std::swap(z[i], z[j]);
         m_lat->permute(i, j);
      }
      j = i;
   }
   // Put the new vector v in place of v_j, then in first place, and update its L2 norm.
   IntVec &rowj = m_lat->getBasis()[j];
   for (h = 0; h < dim; h++)
      rowj[h] = m_bv[h];
   m_lat->permute(0, j);
   //  if (norm == L2NORM) m_lat->setVecNorm(m_lMin2, 0); // Ok only if v is shortest.
   m_lat->setNegativeNorm(0);
}

//=========================================================================

// This one works only if v is a shortest vector for the L2 norm!
template<typename Int, typename Real>
void ReducerBB<Int, Real>::insertBasisVectorLLL(std::vector<std::int64_t> &z) {
   // Add the new vector `v` as a new row to the basis.
   const int64_t dim = m_lat->getDim();
   IntVec &newrow = m_lat->getBasis()[dim];
   for (int64_t h = 0; h < dim; h++)
      newrow[h] = m_bv[h];
   LLLConstruction0<Int, Real>(m_lat->getBasis(), 0.5, dim+1, dim);
}


//=========================================================================

template<typename Int, typename Real>
bool ReducerBB<Int, Real>::tryZShortVecOld(int64_t j, bool &smaller, NormType norm) {
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
   Real dc, center, x, m_spjm1;
   std::int64_t min0, max0;     // Interval boundaries for the z_j.
   std::int64_t zlow, zhigh;    // Current pointers on the left and right of the center.
   bool high;      // Indicates if we are on the right (true) or the left of the center.
   bool stillHope;
   int64_t i, k;
   std::int64_t temp;
   const int64_t dim = m_lat->getDim();
   IntVec vtest;  // A new shortest vector candidate to be tested.

   ++m_countNodes;
   if (m_countNodes > maxNodesBB) {
      std::cerr << "*****   m_countNodes > maxNodesBB = " << maxNodesBB << std::endl;
      return false;
   }

   // Compute an interval that contains the admissible values of z_j.
   // This computation works for the L2 or L1 norm.
   // m_sjp[j] is s_j(p) in the guide.
   center = 0.0;   dc = 0;
   for (i = j+1;  i < dim; ++i)
       center -= m_L[i][j] * m_zLR[i];  // This is `c_j` in the guide.
   // m_L contains the \tilde\ell_{i,j} of the guide.
   // This dc is the distance from the center to the boundaries.
   // m_lMin2 contains the square length of current shortest vector with the selected norm.
   if (m_decomp == CHOLESKY) {
      dc = sqrt((m_lMin2 - m_sjp[j]) / m_dc2[j]);
      // std::cout << " With Cholesky, j = " << j << ", center = " << center << ",  dc = " << dc << "\n";
   }
   if (m_decomp == TRIANGULAR && norm == L2NORM) {
       dc = sqrt(m_lMin2 - m_sjp[j]) / m_L[j][j];
   }
   if (m_decomp == TRIANGULAR && norm == L1NORM) {
      dc = (m_lMin1 - m_sjp[j]) / m_L[j][j];
      // if (dc < 0.0) dc = 0.0;   // Not needed and bad.
      //std::cout << " Triangular + L1, j = " << j << ", center = " << center << ",  dc = " << dc << "\n";
   }

   // Compute two integers min0 and max0 that are the min and max integers in the interval.
   if (!m_foundZero) min0 = 0;     // We are at the beginning for j=dim-1, min will be zero.
   else {
      x = center - dc;
      NTL::conv(min0, trunc(x));
      if (x > 0.0) ++min0;
   }
   x = center + dc;
   NTL::conv(max0, trunc(x));
   if (x < 0.0) --max0;

   if (m_verbose > 3) {
      for (i = 0;  i <= j; ++i) m_z[i] = 0;
      std::cout << "tryZ: vector z: " << m_z << "\n";
      std::cout << "  j="<<j << ",  center[j]= " << center << ", low bound= "<< min0 << ", up bound= "<< max0 << std::endl;
   }
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
      if (high) m_z[j] = zhigh;
      else m_z[j] = zlow;    // For j = dim-1, this will be 0.
      NTL::conv (m_zLR[j], m_z[j]);

      // Computing m_sjp[j-1].
      x = m_zLR[j] - center;
      if (m_decomp == TRIANGULAR && norm == L1NORM)
         m_spjm1 = m_sjp[j] + abs(x) * m_L[j][j]; // This is s_{j-1}(1) in the guide.
      else
         m_spjm1 = m_sjp[j] + x * x * m_dc2[j]; // This is s_{j-1}(2).

      if (j == 0) {
         // All the zj have been selected: we now have a candidate vector to test!
         if (m_decomp == TRIANGULAR && norm == L1NORM)
            stillHope =  (m_lMin1 > m_spjm1);
         else
            stillHope =  (m_lMin2 > m_spjm1);
         if (stillHope) {
            ++m_countLeaves;
            // Length of shortest is not too yet too small. Check if we have a shorter nonzero vector.
            if (!m_foundZero) {
               // The first vector found will always be zero, we discard it.
               m_foundZero = true;
               //std::cout << " Zero vector z: " << m_z << "\n";
            } else {
               vtest.SetLength(dim);
               for (k = 0; k < dim; k++)
                  vtest[k] = 0;
               for (k = 0; k < dim; k++) {
                  if (m_z[k] != 0)
                     ModifVect(vtest, m_lat->getBasis()[k], m_z[k], dim);
               }
               // The new shortest vector is now in `vtest`. We compute its norm.
               if (norm == L2NORM) {
                  ProdScal<Int>(vtest, vtest, dim, x);
                  smaller = (x < m_lMin2);
               } else {
                  // Compute the square length for the L1 norm.
                  CalcNorm<Int, Real>(vtest, dim, x, L1NORM);
                  //std::cout << " vector z: " << m_z << "\n";
                  //std::cout << " vector vtest: " << vtest << ",  L1 norm x = " << x << "\n";
                  smaller = (x < m_lMin1);
                  if (smaller) NTL::conv(m_lMin1, x);
                  x = x * x;
               }
               // This must work for either L1 or L2 norm, m_lMin2 is the square norm in both cases.
               if (smaller) {
                  // The new vector is shorter.
                  NTL::conv(m_lMin2, x);  // Update the smallest square norm.
                  m_bv = vtest;
                  m_zShort = m_z;
                  //std::cout << "tryZ: Testing shorter vector candidate vtest: " << vtest << ",  x = " << x << "\n";
                  // std::cout << " Shorter vector m_bv: " << m_bv << "\n";
                  //std::cout << "tryZ: found shorter vector vtest, Shortest L1 norm: " << m_lMin1 << ",  m_lMin2: " << m_lMin2 << "\n";
                  if (m_verbose > 2) {
                       std::cout << "tryZ: found shorter vector candidate m_vtest: " << m_vtest << "\n";
                       std::cout << "      vector z: " << m_zShort << "\n";
                       std::cout << "      square norm: " << m_lMin2 << "\n";
                       if (norm == L1NORM) std::cout << "      L1 norm: " << m_lMin1 << "\n";
                  }
               }
            }
         }
      } else {
         if (m_decomp == TRIANGULAR && norm == L1NORM)
            stillHope = (m_lMin1 > m_spjm1);
         else
            stillHope = (m_lMin2 > m_spjm1);
         if (stillHope) {
            // There is still hope; we continue the recursion.
            m_sjp[j - 1] = m_spjm1;
            if (!tryZShortVecOld(j - 1, smaller, norm)) return false;
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

   // To avoid creating them in recursive calls, some variables below would have to be arrays indexed by j.
   Real center, dc, x, spjm1;
   int64_t min0, max0;     // Interval boundaries for the z_j.
   int64_t zlow, zhigh;    // Current pointers on the left and right of the center.
   int64_t temp;
   bool high;      // Indicates if we are on the right (true) or the left of the center.
   bool stillHope;
   int64_t i, k;
   const int64_t dim = m_lat->getDim();

   ++m_countNodes;
   if (m_countNodes > maxNodesBB) {
      std::cerr << "*****   m_countNodes > maxNodesBB = " << maxNodesBB << std::endl;
      return false;
   }

   // Compute an interval that contains the admissible values of z_j.
   // This computation works for the L2 or L1 norm.
   // m_sjp[j] is s_j(p) in the guide, was computed one level higher in the tree.
   center = 0.0;
   for (i = j+1;  i < dim; ++i)
       center -= m_L[i][j] * m_zLR[i];  // This is `c_j` in the guide.

   // m_L contains the \tilde\ell_{i,j} of the guide.
   // This dc is the distance from the center to the boundaries.
   // m_lMin2 contains the square length of current shortest vector with the selected norm.
   dc = 0.0;
   if (m_decomp == CHOLESKY) {
      dc = sqrt((m_lMin2 - m_sjp[j]) / m_dc2[j]);
      // std::cout << " With Cholesky, j = " << j << ", center[j] = " << center << ",  dc = " << dc << "\n";
   }
   else if (m_decomp == TRIANGULAR) {
      if (norm == L2NORM)
         dc = sqrt(m_lMin2 - m_sjp[j]) / m_L[j][j];
      else if (norm == L1NORM) {
         dc = (m_lMin1 - m_sjp[j]) / m_L[j][j];
         // if (dc < 0.0) dc = 0.0;   // Not needed and bad.
         //std::cout << " Triangular + L1, j = " << j << ", center[j] = " << center << ",  dc = " << dc << "\n";
      }
   }
   else
      std::cerr << "*** tryZ: m_decomp undefined\n";

   // Compute two integers min0 and max0 that are the min and max integers in the interval.
   if (!m_foundZero)
      min0 = 0;     // We are at the beginning for j=dim-1, min will be zero.
   else {
      x = center - dc;
      NTL::conv(min0, trunc(x));
      if (x > 0.0) ++min0;
   }
   x = center + dc;
   NTL::conv(max0, trunc(x));
   if (x < 0.0) --max0;

   if (m_verbose > 3) {
      for (i = 0;  i <= j; ++i) m_z[i] = 0;
      std::cout << "tryZ: vector z: " << m_z << "\n";
      std::cout << "  j="<<j << ",  center[j]= " << center <<
            ", low bound= "<< min0 << ", up bound= "<< max0 << std::endl;
   }
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
      if (high) m_z[j] = zhigh;
      else m_z[j] = zlow;    // For j = dim-1, this will be 0.
      NTL::conv (m_zLR[j], m_z[j]);
      m_maxZj = std::max(m_maxZj, std::abs(m_z[j]));  // Update largest absolute z_j.

      // Computing m_sjp[j-1].
      x = m_zLR[j] - center;
      if (m_decomp == TRIANGULAR && norm == L1NORM)
         spjm1 = m_sjp[j] + abs(x) * m_L[j][j]; // This is s_{j-1}(1) in the guide.
      else
         spjm1 = m_sjp[j] + x * x * m_dc2[j]; // This is s_{j-1}(2).

      if (j == 0) {
         // All the z_j have been selected: we now have a candidate vector to test!
         if (m_decomp == TRIANGULAR && norm == L1NORM)
            stillHope =  (m_lMin1 > spjm1);
         else
            stillHope =  (m_lMin2 > spjm1);
         if (stillHope) {
            ++m_countLeaves;
            // Length of shortest is not too yet too small. Check if we have a shorter nonzero vector.
            if (!m_foundZero) {
               // The first vector found will always be zero, we discard it.
               m_foundZero = true;
               //std::cout << " Zero vector z: " << m_z << "\n";
            } else {
               for (k = 0; k < dim; k++)
                  m_vtest[k] = 0;
               for (k = 0; k < dim; k++) {
                  if (m_z[k] != 0)
                     ModifVect(m_vtest, m_lat->getBasis()[k], m_z[k], dim);
               }
               // The new shortest vector is now in `m_vtest`. We compute its norm.
               if (norm == L2NORM) {
                  ProdScal<Int>(m_vtest, m_vtest, dim, x);
                  smaller = (x < m_lMin2);
               } else {
                  // Compute the square length for the L1 norm.
                  CalcNorm<Int, Real>(m_vtest, dim, x, L1NORM);
                  //std::cout << " vector z: " << m_z << "\n";
                  //std::cout << " vector m_vtest: " << m_vtest << ",  L1 norm x = " << x << "\n";
                  smaller = (x < m_lMin1);
                  if (smaller) NTL::conv(m_lMin1, x);  // Updateb new min1.
                  x = x * x;
               }
               // This must work for either L1 or L2 norm, m_lMin2 is the square norm in both cases.
               if (smaller) {
                  // The new vector is shorter.
                  NTL::conv(m_lMin2, x);  // Update the smallest square norm.
                  m_bv = m_vtest;
                  m_zShort = m_z;
                  if (m_verbose > 2) {
                     std::cout << "tryZ: found shorter vector candidate m_vtest: " << m_vtest << "\n";
                     std::cout << "      vector z: " << m_zShort << "\n";
                     std::cout << "      square norm: " << m_lMin2 << "\n";
                     if (norm == L1NORM) std::cout << "      L1 norm: " << m_lMin1 << "\n";
                  }
               }
            }
         }
      } else {
         if (m_decomp == TRIANGULAR && norm == L1NORM)
            stillHope = (m_lMin1 > spjm1);
         else
            stillHope = (m_lMin2 > spjm1);
         if (stillHope) {
            // There is still hope; we continue the recursion.
            m_sjp[j - 1] = spjm1;
            if (!tryZShortVec(j - 1, smaller, norm)) return false;  // Recursive call.
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
bool ReducerBB<Int, Real>::shortestVector() {
   /*
    * Finds a shortest non-zero vector, using branch-and-bound, with L1 or L2 norm.
    * When it succeeds, returns true, and the squared shortest
    * vector length will be in m_lMin2, regardless of the selected norm.
    *
    * This function uses (directly or indirectly) the following class variables:
    *    m_lMin1, m_lMin2, m_decomp, m_Bounds2, m_sjp, m_countNodes, m_foundZero,
    *    m_bv, m_zShort, m_L, m_zLR, m_z, m_dc2, ....  and more.
    * From m_lat (current lattice object):
    *    norm, sortBasisNoDual, updateScalL2Norm, getBasis,
    */
   NormType norm = m_lat->getNormType();
   if ((norm != L1NORM) & (norm != L2NORM)) {
      std::cerr << "shortestVector: only L1 and L2 norms are supported";
      return false;
   }
   const int64_t dim = m_lat->getDim();  // Lattice dimension
   //std::cout << " Start shortestVector, basis = \n" << m_lat->getBasis() << "\n";

   bool smaller = false;  // Will change if/when we find a smaller vector.
   Real x(0.0);
   // In all cases, we update the L2 norm and we sort the basis by L2 lengths.
   // This is better for Cholesky, otherwise the decomp will fail more rapidly
   // due to floating-point errors.
   if (norm != L2NORM)  m_lat->setNegativeNorm();
   m_lat->updateScalL2Norm(0, dim);  // `m_vecNorm` will now contain the square L2 norms.
   m_lat->sortBasis(0);              // Vectors are sorted by L2 norms.
   // for (int64_t k = 0; k < dim; k++)  m_bv[k] = m_lat->getBasis()[0][k];
   m_bv = m_lat->getBasis()[0];
   //std::cout << " shortestVector, after sort, m_bv = " << m_bv <<
   //          ",  squared L2 norm = " << m_lat->getVecNorm(0) << "\n";

   // We put in m_lMin2 the approximate square norm of the shortest vector in the basis,
   // for either the L1 or L2 norm.  For the L2 norm, we know it is the first vector.
   if (norm == L2NORM) {
      NTL::conv(m_lMin2, m_lat->getVecNorm(0));
   } else {
      // Looking for the shortest vector in basis according to the L1 norm.
      CalcNorm<Int, Real>(m_lat->getBasis()[0], dim, m_lMin1, norm);
      for (int64_t k = 1; k < dim; k++) {
         CalcNorm<Int, Real>(m_lat->getBasis()[k], dim, x, norm);
         if (x < m_lMin1) m_lMin1 = x;
      }
      m_lMin2 = m_lMin1 * m_lMin1;  // Squared shortest length with L1 norm.
      //std::cout << " In SV before BB, Shortest L1 norm: " << m_lMin1 << ",  m_lMin2: " << m_lMin2 << "\n";
   }

   // If we already have a shorter vector than the minimum threshold, we stop right away.
   // This is useful for the seek programs in LatMRG.
   if (m_lMin2 <= m_Bounds2[dim - 1]) return false;

   // std::cout << " shortestVector, basis before decomposition = \n" << m_lat->getBasis() << "\n";
   if (m_decomp == CHOLESKY) {
      if (!calculCholeskyLDL()) return false;
   } else if (m_decomp == TRIANGULAR) {  // Just for testing; this is very slow!
      if (!calculTriangularL()) return false;
   } else {
      std::cerr << "shortestVector: decomp value not supported";
      return false;
   }
   if (m_verbose > 2) {
      // std::cout << "shortestVector: basis after decomposition = \n" << m_lat->getBasis() << "\n";
      std::cout << "shortestVector: matrix m_L = \n" << m_L  << "\n";
      std::cout << " vector m_dc2 = " << m_dc2  << "\n";
   }
   // Perform the branch and bound.
   // The following variables are used and updated in `tryZShortVec`.
   m_sjp[dim - 1] = 0.0;
   m_countNodes = 0;
   m_countLeaves = 0;
   m_maxZj = 0;
   smaller = false;
   m_foundZero = false;
   for (long j=0; j < dim; j++) m_z[j] = 0;
   if (!tryZShortVec (dim - 1, smaller, norm)) // We search for a shortest vector.  ******
      return false;
   if (smaller) {
      // We found a shorter vector. It is in m_bv and its square length is in m_lMin2.
      // The short vector m_bv is not yet put in the basis.
      // The following does that and is useful only if we want to continue working with this basis.
      //std::cout << "shortestVector: found smaller vector in BB, shortest L1 norm: " << m_lMin1 << ",  m_lMin2: " << m_lMin2 << "\n";
      //std::cout << "shortest vector = " << m_bv << "\n";
      insertBasisVector(m_zShort);
      // std::cout << "After insertBasis, Shortest L1 norm: " << m_lMin1 << ",  m_lMin2: " << m_lMin2 << "\n";
   }
   if (m_verbose > 1) {
      std::cout << "shortestVector: total number of calls to tryZ: " << m_countNodes << "\n";
      std::cout << "  shortest vector = " << m_bv << "\n";
      std::cout << "  square length m_lMin2 = " << m_lMin2 << "\n";
   }
   return true;
}


//=========================================================================

// If m_countNodes > MaxNodesBB, returns false, otherwise returns true.
template<typename Int, typename Real>
bool ReducerBB<Int, Real>::tryZMink(int64_t j, int64_t i, int64_t Stage, bool &smaller,
      const IntMat &WTemp)  {
   std::int64_t max0, min0;
   Real x, dc, center;
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
   if (j < dim - 1) {
      center = 0.0;
      for (i = j+1;  i < dim; ++i)
          center -= m_L[i][j] * m_zLR[i];  // This is `c_j` in the guide.
      dc = sqrt((m_lMin2 - m_sjp[j]) / m_dc2[j]);

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
      center = 0;
      zlow = 0;
      high = true;
      if (Stage == 2) {
         min0 = 1;
         max0 = 1;
         zhigh = 1;
      } else {
         min0 = 2;
         zhigh = 2;
         NTL::conv(max0, trunc(sqrt((m_lMin2 - m_sjp[j]) / m_dc2[j])));
      }
   }

   Real temp;
   /* On essaie maintenant chacun des z[j] dans l'intervalle, en      */
   /* commencant par le centre puis en alternant d'un cote a l'autre. */
   while (zlow >= min0 || zhigh <= max0) {
      if (high) {
         m_z[j] = zhigh;
      } else {
         m_z[j] = zlow;
      }
      m_zLR[j] = m_z[j];
      m_maxZj = max (m_maxZj, abs(m_z[j]));

      // Calcul de m_sjp[j-1].
      x = m_zLR[j] - center;
      if (j == 0) {
         Real tmps_n2 = m_sjp[0] + x * x * m_dc2[0];
         if (tmps_n2 < m_lMin2) {
            // On verifie si on a vraiment trouve un vecteur plus court
            IntVec &row1 = m_lat->getBasis()[dim-1];
            m_vtest = row1;
            for (k = 0; k < dim - 1; k++) {
               if (m_z[k] != 0) {
                  IntVec &row1 = m_lat->getBasis()[k];
                  //NTL::Mat_row<Int> row1(m_lat->getBasis(), k);
                  ModifVect(m_vtest, row1, m_z[k], dim);
               }
            }
            if (Stage == 3) {
               IntVec &row1 = m_lat->getBasis()[dim-1];
               //NTL::Mat_row<Int> row1(m_lat->getBasis(), dim - 1);
               ModifVect(m_vtest, row1, m_zLR[dim - 1] - 1.0, dim);
            }

            ProdScal<Int>(m_vtest, m_vtest, dim, S1);
            NTL::conv(S4, m_lat->getVecNorm(dim - 1));
            if (S1 < S4) {
               if (Stage == 2) {
                  smaller = true;
                  if (!PreRedLLLMink) m_zShort = m_z;
                  else {
                     for (k = 1; k < dim; k++) {
                        // NTL::Mat_row<Int> row1(WTemp, k);
                        ProdScal<Int>(m_vtest, WTemp[k], dim, S2);
                        Quotient(S2, mR, S3);
                        NTL::conv(m_zShort[k], S3);
                     }
                     m_zShort[dim - 1] = 1;
                  }
               } else if (Stage == 3 && !PreRedLLLMink) {
                  if (GCD2vect(m_z, i, dim) == 1) {
                     m_zShort = m_z;
                     smaller = true;
                  }
               } else {
                  for (k = 0; k < dim; k++) {
                     // NTL::Mat_row<Int> const row1(WTemp, k);
                     ProdScal<Int>(m_vtest, WTemp[k], dim, S2);
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
                  m_bv = m_vtest;
                  return true;
               }
            }
         }
      } else { // j > 0
         m_sjp[j - 1] = m_sjp[j] + x * x * m_dc2[j];
         if (m_lMin2 >= m_sjp[j - 1]) {
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
    * Tries to shorten m_lat->getBasis()[i] using branch-and-bound, in Minkowski Reduction.
    * Stage is 2 or 3.  z[i] = 1 if Stage = 2, z[i] >= 2 if Stage = 3.
    * Stops and returns false if not finished after examining MaxNodesBB
    * nodes in the branch-and-bound tree.  When succeeds, returns true.
    * Assumes that the norm is Euclidean.
    */
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
      m_lat->updateVecNorm();
   }
   if (!calculCholeskyLDL()) return false;
   m_countNodes = 0;
   m_countLeaves = 0;
   m_maxZj = 0;
   m_sjp[dim - 1] = 0.0;
   if (!tryZMink(dim - 1, i, Stage, smaller, WTemp)) return false;

   if (PreRedLLLMink) {
      // On remet l'anciennne base, celle d'avant LLL, avant de considerer
      // la prochaine m_lat->dimension.
      m_lat->getBasis() = VTemp;
      m_lat->updateVecNorm();
      for (h = 0; h < dim; h++)
         taboo[h] = TabooTemp[h];
   }
   if (smaller) {
      /* On a trouve un plus court vecteur qui ameliore m_lat->getBasis()[k].  */
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
   } else if (Stage == 2)
      taboo[dim - 1] = true;
   m_lat->permute(i, dim - 1);
   return true;
}

//=========================================================================

// This works only for the L2 norm.
template<typename Int, typename Real>
bool ReducerBB<Int, Real>::reductMinkowski(int64_t d) {
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
bool ReducerBB<Int, Real>::reductMinkowski(IntLattice<Int, Real> &lat, int64_t d) {
   setIntLattice(lat);
   return ReducerBB<Int, Real>::reductMinkowski(d);
}

//=========================================================================

template<typename Int, typename Real>
bool ReducerBB<Int, Real>::shortestVector(IntLattice<Int, Real> &lat) {
   setIntLattice(lat);
   return ReducerBB<Int, Real>::shortestVector();
}

//============================================================================

template class ReducerBB<std::int64_t, double> ;
template class ReducerBB<NTL::ZZ, double> ;
template class ReducerBB<NTL::ZZ, xdouble> ;
template class ReducerBB<NTL::ZZ, quad_float> ;
template class ReducerBB<NTL::ZZ, NTL::RR> ;

}     // namespace LatticeTester

#endif 
