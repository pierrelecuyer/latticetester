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

#include "NTL/LLL.h"
#include "NTL/tools.h"
#include "NTL/ZZ.h"
#include "NTL/RR.h"

#include "latticetester/EnumTypes.h"
#include "latticetester/Util.h"
#include "latticetester/IntLattice.h"
// #include "latticetester/BasisConstruction.h"
#include "latticetester/NTLWrap.h"
//#include "latticetester/LLL_FPInt.h"
//#include "latticetester/LLL_RR_lt.h"

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

private:
    // Local typedefs for matrix and vector types needed in the class.
    typedef NTL::vector<Int> IntVec;
    typedef NTL::matrix<Int> IntMat;
    typedef NTL::vector<Real> RealVec;
    typedef NTL::matrix<Real> RealMat;

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
     * Initializes all matrices used in the following
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
     *
     * This function uses only the basis of the internal lattice, its vector
     * lengths, and scalar products. It never uses its m-dual.
     * To compute a shortest vector in the m-dual, one must first call `dualize`
     * on the target `IntLattice` object.
     * It is strongly recommended to use `redBKZ` or `redLLLNTL` to pre-reduce
     * the basis before invoking this method; this is not done automatically.
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
     * This method performs pairwise reduction sequentially on all vectors
     * of the basis whose indices are greater of equal to `dim >=0`,
     * as proposed in \cite rDIE75a.
     * For this function to work, both the primal and m-dual bases must be
     * maintained together.
     * The boolean vector `taboo[]` is used internally in Minkowski reduction only:
     * when this vector is not `NULL`, `taboo[j]=true' means that the j-th vector
     * should not be modified.
     */
    void redDieter(int64_t dim, bool taboo[] = NULL);

    /**
     * Same as `redDieter(dim)` but the choice of
     * vectors on which to perform pairwise reduction is randomized,
     * using a simple RNG from the standard C library, with the given integer seed.
     */
    void redDieterRandomized(int64_t dim, int64_t seed);

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
    void redLLLOld(double delta = 0.999999, std::int64_t maxcpt = 1000000000,
            int64_t dim = 0);

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
        if (m_lat->getNormType() == L2NORM)
            return sqrt(m_lMin2);
        else
            return m_lMin;
    }

    /**
     * Returns the length of the current *last* basis vector in the lattice.
     */
    Real getMaxLength() {
        if (m_lat->getNormType() == L2NORM)
            return sqrt(m_lat->getVecNorm(m_lat->getDim() - 1));
        else
            return m_lat->getVecNorm(m_lat->getDim() - 1);
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
        if (m_lat->getDim() > m_maxDim)
            init(m_lat->getDim());
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
    bool redBBMink(int64_t i, int64_t d, int64_t Stage, bool &smaller,
            bool taboo[] = NULL);

    /**
     * Recursive procedure used in `redBBMink` to try find shorter vectors.
     */
    bool tryZMink(int64_t j, int64_t i, int64_t Stage, bool &smaller,
            const IntMat &WTemp);

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
    void pairwiseRedDual(int64_t i, bool taboo[] = NULL);

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

    /*
     * Method used in `shortestVector` to perform a transformation of
     * stage 3 described in \cite rAFF85a.
     * We call this when we have found a shorter vector which is a linear combination of the
     * previous basis vectors with coefficients given in z. This procedure updates the basis
     * so that the new shortest vector is in first and the other vectors are adjusted accordingly.
     * With the L2 norm, an alternative might be to add the new shortest vector to
     * the basis to get dim+1 generating vectors, and apply LLL to recover a new basis.
     *  ***  We should test this and compare the speeds.
     *
     * Here, we assume that \f$\sum_{i=1}^t z_i V_i\f$ is a
     * short vector that will enter the basis. Then we try to reduce some vectors
     * by looking for indices \f$i < j\f$ such that \f$|z_j| > 1\f$ and
     * \f$q=\lfloor z_i/z_j\rfloor\not0\f$, and adding \f$q V_i\f$ to
     * \f$V_j\f$ when this happens. We return in \f$k\f$ the last index
     * \f$j\f$ such that \f$|z_j|=1\f$.
     * This method does not change the m-dual.
     */
    void transformStage3ShortVec(std::vector<std::int64_t> &z, int64_t &k);

    /**
     * Method used in `reductMinkowski` to perform a transformation of
     * stage 3 described in \cite rAFF85a. Assumes that
     * \f$\sum_{i=1}^t z_i V_i\f$ is a
     * short vector that will enter the basis. Tries to reduce some vectors
     * by looking for indices \f$i < j\f$ such that \f$|z_j| > 1\f$ and
     * \f$q=\lfloor z_i/z_j\rfloor\not0\f$, and adding \f$q V_i\f$ to
     * \f$V_j\f$ when this happens. Returns in \f$k\f$ the last index
     * \f$j\f$ such that \f$|z_j|=1\f$.    **This method updates the m-dual.**
     */
    void transformStage3Mink(std::vector<std::int64_t> &z, int64_t &k);

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
     * Maximum number of transformations in the method `PreRedDieter`.
     * After <tt>MAX_PRE_RED</tt> successful transformations have been
     * performed, the prereduction is stopped.
     *
     */
    std::int64_t MAX_PRE_RED = 1000000;

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
    IntVec m_bw;   // Saves current shortest vector in dual basis
    Real m_lMin;   // The norm of the shortest vector in the primal basis
                   // according to the norm considered
    Real m_lMin2;  // Squared L2-norm of the shortest vector in the primal basis.

    Int m_bs;      // Used in pairwiseRedPrimal and pairwiseRedDual.
    Real m_ns;     // Used in pairwiseRedPrimal and pairwiseRedDual.
    RealVec m_nv;  // Used in pairwiseRedDual.

    RealVec m_n2, m_dc2; // Vectors used in BB algorithms, in tryZShortVec and tryZMink.
    RealMat m_c0, m_c2;  // Matrices used in the Cholesky decomposition.
                         // We must avoid resizing them, because this is expensive!
    int64_t *m_IC;       // Indexes used in Cholesky

    std::vector<std::int64_t> m_zLI;   // Vector of values of z_i.
    RealVec m_zLR;     // Same values of z_i in floating point;
                       // we need them in FP when calculating bounds.
    std::vector<std::int64_t> m_zShort;  // Values of z_i for shortest vector.
    std::int64_t m_countNodes = 0;  // Number of visited nodes in the BB tree
    std::int64_t m_countDieter = 0; // Number of attempts since last successful
                                    // Dieter transformation
    std::int64_t m_cpt;  // Number of successes in pre-reduction transformations
    bool m_foundZero;    // = true -> the zero vector has been handled

    RealMat m_cho2, m_gramVD;  // Matrices used by redLLLOld.

};
// End class ReducerBB

//=========================================================================

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
    if (dim2 <= 2)
        dim2++;
    m_c0.resize(dim1, dim1);
    m_c2.resize(dim1, dim1);
    m_cho2.resize(dim2, dim2);
    m_gramVD.resize(dim2, dim2);
    // Indices in m_IC can go as high as maxDim + 2.
    m_IC = new int64_t[2 + dim1];

    m_nv.resize(dim1);
    m_bv.resize(dim1);
    m_bw.resize(dim1);
    m_n2.resize(dim1);
    m_zLR.resize(dim1);
    m_zLI.resize(dim1);
    m_zShort.resize(dim1);
    m_dc2.resize(dim1);
    m_BoundL2.resize(dim1);

    m_lMin = std::numeric_limits<double>::max();
    m_lMin2 = m_lMin;
    for (int64_t i = 0; i < dim1; i++) {
        m_zLI[i] = -1;
        m_zShort[i] = -1;
        m_BoundL2[i] = -1;
        m_IC[i] = -1;
    }
    m_countNodes = 0;
    m_countDieter = 0;
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
    if (this != &red)
        copy(red);
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
    m_nv = red.m_nv;
    m_bv = red.m_bv;
    m_bw = red.m_bw;
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
    m_countDieter = 0;
    m_foundZero = false;
    m_cpt = 0;
    PreRedLLLMink = false;
    maxNodesBB = 1000000000;
    if (m_IC != 0)
        delete[] m_IC;
    m_IC = new int64_t[3 + m_lat->getDim()];
    for (int64_t i = 0; i < 2 + m_lat->getDim(); i++)
        m_IC[i] = red.m_IC[i];
}

//=========================================================================

template<typename Int, typename Real>
ReducerBB<Int, Real>::~ReducerBB() {
    m_c0.clear();
    m_c2.clear();
    m_cho2.clear();
    m_gramVD.clear();
    m_nv.clear();
    m_bv.clear();
    m_bw.clear();
    m_n2.clear();
    m_zLR.clear();
    m_zLI.clear();
    m_zShort.clear();
    m_dc2.clear();
    m_BoundL2.clear();
    delete[] m_IC;
}

//=========================================================================

template<typename Int, typename Real>
void ReducerBB<Int, Real>::setBoundL2(const RealVec &thresholds, int64_t dim1,
        int64_t dim2) {
    m_BoundL2.resize(dim2);
    for (int64_t i = dim1; i < dim2; i++)
        m_BoundL2[i] = thresholds[i];
}

//=========================================================================

template<typename Int, typename Real>
inline void ReducerBB<Int, Real>::calculGramVD() {
    // Returns in m_gramVD the matrix of scalar products m_lat->V[i] * m_lat->V[j].
    // The vector m_lat->V.vecNorm contains only the values m_lat->V[i] * m_lat->V[i].
    // Used only by redLLLOld.
    // %%%  Could easily be made static.
    const int64_t dim = m_lat->getDim();
    Int temp;
    for (int64_t i = 0; i < dim; i++) {
        NTL::matrix_row<IntMat> row1(m_lat->getBasis(), i);
        temp = 0;
        for (int64_t h = 0; h < dim; h++)
            NTL::MulAddTo(temp, row1[h], row1[h]);
        NTL::conv(m_gramVD[i][i], temp);
        for (int64_t j = i + 1; j < dim; j++) {
            NTL::matrix_row<IntMat> row2(m_lat->getBasis(), j);
            temp = 0;
            for (int64_t h = 0; h < dim; h++)
                NTL::MulAddTo(temp, row1[h], row2[h]);
            NTL::conv(m_gramVD[i][j], temp);
            m_gramVD[j][i] = m_gramVD[i][j];
        }
    }
}

//=========================================================================

template<typename Int, typename Real>
inline void ReducerBB<Int, Real>::miseAJourGramVD(int64_t j) {
    // Recalcule la j-ieme ligne et la j-ieme colonne de la matrice des
    // produits scalaires.  Used only in redLLLOld.
    // %%%  Could easily be made static by passing the rows and gramVD as parameters.
    const int64_t dim = m_lat->getDim();
    for (int64_t i = 0; i < dim; i++) {
        NTL::matrix_row<IntMat> row1(m_lat->getBasis(), i);
        NTL::matrix_row<IntMat> row2(m_lat->getBasis(), j);
        ProdScal<Int>(row1, row2, dim, m_gramVD[i][j]);
        m_gramVD[j][i] = m_gramVD[i][j];
    }
}

//=========================================================================

template<typename Int, typename Real>
void ReducerBB<Int, Real>::permuteGramVD(int64_t i, int64_t j, int64_t n) {
    // %%%  Could swap the rows directly.
    int64_t k;
    for (k = 0; k < n; k++) {
        std::swap(m_gramVD[i][k], m_gramVD[j][k]);
    }
    for (k = 0; k < n; k++) {
        std::swap(m_gramVD[k][i], m_gramVD[k][j]);
    }
}

//=========================================================================

template<typename Int, typename Real>
void ReducerBB<Int, Real>::calculCholesky2LLL(int64_t n, int64_t j) {
    // Recalcule les n premieres entrees de la colonne j de la matrice de
    // Cholesky d'ordre 2.   Pour redLLLOld.
    m_cho2[0][j] = m_gramVD[0][j];
    for (int64_t i = 1; i <= n; i++) {
        m_cho2[i][j] = m_gramVD[i][j];
        for (int64_t k = 0; k < i; k++) {
            m_cho2[i][j] -= (m_cho2[k][j] / m_cho2[k][k]) * m_cho2[k][i];
        }
    }
}

//=========================================================================

template<typename Int, typename Real>
inline void ReducerBB<Int, Real>::calculCholesky2Ele(int64_t i, int64_t j) {
    // Recalcule l'entree [i][j] de la matrice de Cholesky d'ordre 2.
    // Used only by redLLLOld.
    m_cho2[i][j] = m_gramVD[i][j];
    for (int64_t k = 0; k < i; k++) {
        m_cho2[i][j] -= m_cho2[k][i] * (m_cho2[k][j] / m_cho2[k][k]);
    }
}

//=========================================================================

void negativeCholesky() {
    std::cout << "\n***** Negative diagonal element in Cholesky Decomposition\n"
            << std::endl;
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
    NTL::conv(m2, m_lat->getModulo());
    m2 = m2 * m2;
    int64_t d = dim;
    //if(m_lat->withDual())
    //  d = dim / 2; // If we use the Dual, we compute Cholesky
    // with the Dual

    // Compute the d first lines of C0 with the primal Basis.
    for (i = 0; i < d; i++) {
        m_lat->updateScalL2Norm(i);
        NTL::matrix_row<IntMat> row1(m_lat->getBasis(), i);
        for (j = i; j < dim; j++) {
            if (j == i)
                NTL::conv(m_c2[i][i], m_lat->getVecNorm(i));
            else {
                NTL::matrix_row<IntMat> row2(m_lat->getBasis(), j);
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
            } else
                C0[i][j] = m_c2[i][j] / DC2[i];
            //add for test
            /**  if(i!=j && i<j)
             std::cout<< "C0("<<i<<","<<j<<")="<<C0[i][j]<<" ";
             else if (i==j)
             std::cout<< "C0("<<i<<","<<j<<")="<<DC2[i]<<" ";
             else
             std::cout<< "C0("<<i<<","<<j<<")="<<"0"<<" ";	*/
        }
        // std::cout<<""<<std::endl;
    }

    // Compute the d last lines of C0 with the dual Basis.
    /* This operation with the dual is needed in case of high dimension
     * and large number (30 bits). The choleski decomposition
     * convert numbers in double which is not sufficient in that case.
     * You need to use RR of the NTL library for this calculation.
     */
    //if(m_lat->withDual()){
    //  for (i = dim-1; i >= d; i--)
    //  {
    //    m_lat->updateDualScalL2Norm (i);
    //    for (j = i; j >= 0; j--) {
    //      if (j == i)
    //        NTL::conv (m_c2[i][i], m_lat->getDualVecNorm (i));
    //      else {
    //        NTL::matrix_row<IntMat> row1(m_lat->getDualBasis(), i);
    //        NTL::matrix_row<IntMat> row2(m_lat->getDualBasis(), j);
    //        ProdScal<Int> (row1, row2, dim, m_c2[i][j]);
    //      }
    //      for (k = i + 1; k < dim; k++){
    //        m_c2[i][j] -= C0[k][i] * m_c2[k][j];
    //      }
    //      if (i != j)
    //        C0[i][j] = m_c2[i][j] / m_c2[i][i];
    //    }
    //    DC2[i] = m2 / m_c2[i][i];
    //    if (DC2[i] < 0.0) {
    //      negativeCholesky();
    //      return false;
    //    }
    //    for (j = i + 1; j < dim; j++) {
    //      C0[i][j] = -C0[j][i];
    //      for (k = i + 1; k < j; k++) {
    //        C0[i][j] -= C0[k][i] * C0[k][j];
    //      }
    //    }
    //  }
    //}
    return true;
}

//===========================================================================

template<typename Int, typename Real>
void ReducerBB<Int, Real>::pairwiseRedPrimal(int64_t i, int64_t d, bool taboo[]) {
    const int64_t dim = m_lat->getDim();
    ++m_countDieter;
    m_lat->updateScalL2Norm(i);
    bool modifFlag;
    std::cout << " redDieter, entering pairwiseRedPrimal, i = " << i << ", d = " << d << " **\n";
    std::cout << " vecNorm(i) = " << m_lat->getVecNorm(i) << "\n";

    for (int64_t j = d; j < dim; j++) {
        NTL::matrix_row<IntMat> row1(m_lat->getBasis(), i);
        std::cout << " row(i) = " << row1 << "\n";
        if (i == j)
            continue;
        modifFlag = false;
        NTL::matrix_row<IntMat> row2(m_lat->getBasis(), j);
        ProdScal<Int>(row1, row2, dim, m_ns);
        std::cout << " redPrimal, before divideRound, j = " << j << "\n";
        std::cout << " m_ns = " << m_ns << ",  vecNorm(i) = " << m_lat->getVecNorm(i) << "\n";
        DivideRound<Real>(m_ns, m_lat->getVecNorm(i), m_ns);
        std::cout << " redPrimal, after divideRound, j = " << j << "\n";
        if (m_ns == 0)
            continue;
        std::cout << " redPrimal, before conv, m_ns = " << m_ns << "\n";
        NTL::conv(m_bs, m_ns);
        std::cout << " redPrimal, before if (m_ns ...), m_ns = " << m_ns << "\n";
        if (m_ns < 1000 && m_ns > -1000) {
            m_lat->updateScalL2Norm(j);
                NTL::matrix_row<IntMat> row2(m_lat->getBasis(), j);
                //NTL::matrix_row<IntMat> row1(m_lat->getBasis(), i);
                ModifVect(row2, row1, -m_bs, dim);

            // Verify that m_lat->getBasis()[j] is really shorter
                NTL::matrix_row<IntMat> row1(m_lat->getBasis(), j);
                ProdScal<Int>(row1, row1, dim, m_ns);
            if (m_ns >= m_lat->getVecNorm(j)) {
                NTL::matrix_row<IntMat> row2(m_lat->getBasis(), j);
                //NTL::matrix_row<IntMat> row1(m_lat->getBasis(), i);
                ModifVect(row2, row1, m_bs, dim);
            } else {
                modifFlag = true;
                m_lat->setVecNorm(m_ns, j);
            }
            //std::cout << " redDieter, end of the if, j = " << j << "\n";
         } else {
            std::cout << " redDieter, in the else, j = " << j << "\n";
            NTL::matrix_row<IntMat> row2(m_lat->getBasis(), j);
            //NTL::matrix_row<IntMat> row1(m_lat->getBasis(), i);
            ModifVect(row2, row1, -m_bs, dim);
            m_lat->setNegativeNorm(j);
            modifFlag = true;
        }
        if (modifFlag) {
            std::cout << " redDieter, in modifFlag \n";
            m_countDieter = 0;
            ++m_cpt;
//            if (m_lat->withDual()) {
            if (false) {
                //NTL::matrix_row<IntMat> row1(m_lat->getDualBasis(), i);
                NTL::matrix_row<IntMat> row2(m_lat->getDualBasis(), j);
                ModifVect(row1, row2, m_bs, dim);
                m_lat->setDualNegativeNorm(i);

            }
            if (taboo) {
                taboo[i] = false;
                taboo[j] = false;
            }
        }
    }
    std::cout << " redDieter, end of pairwiseRedPrimal, i = " << i << ", d = " << d << " *****\n";
}

//=========================================================================

template<typename Int, typename Real>
void ReducerBB<Int, Real>::pairwiseRedDual(int64_t i, bool taboo[]) {
    int64_t j;
    const int64_t dim = m_lat->getDim();
    std::cout << " redDieter, entering pairwiseRedDual, i = " << i << "  **\n";

    ++m_countDieter;
    m_lat->updateDualScalL2Norm(i);
    NTL::matrix_row<IntMat> row9(m_lat->getBasis(), i);
    m_bv = row9;
    for (j = 0; j < dim; j++) {
        NTL::matrix_row<IntMat> row2(m_lat->getDualBasis(), j);
        if (i != j) {
            NTL::matrix_row<IntMat> row1(m_lat->getDualBasis(), i);
            ProdScal<Int>(row1, row2, dim, m_ns);
            // ProdScal<Int> (m_lat->getDualBasis ()[i], m_lat->getDualBasis ()[j],
            //           dim, m_ns);
            //std::cout << " redDual, before divideRound, j = " << j << "\n";
            //std::cout << " m_ns = " << m_ns << ",  DualVecNorm(i) = " << m_lat->getDualVecNorm(i) << "\n";
            DivideRound<Real>(m_ns, m_lat->getDualVecNorm(i), m_nv[j]);
            //std::cout << " after divideRound, j = " << j << ", m_nv[j] = " << m_nv[j] << "\n";
            if (m_nv[j] != 0) {
                NTL::conv(m_bs, m_nv[j]);
                std::cout << " after NTL::conv,  m_bs = " << m_bs << "\n";
                NTL::matrix_row<IntMat> row7(m_lat->getBasis(), j);
                ModifVect(m_bv, row7, m_bs, dim);
                std::cout << " after ModifVect \n";
            }
        }
    }
    m_lat->updateScalL2Norm(i);
    ProdScal<Int>(m_bv, m_bv, dim, m_ns);
    if (m_ns < m_lat->getVecNorm(i)) {
        ++m_cpt;
        m_countDieter = 0;
        NTL::matrix_row<IntMat> row6(m_lat->getBasis(), i);
        for (j = 0; j < dim; j++)
            row6[j] = m_bv[j];
        m_lat->setNegativeNorm(i);
        if (taboo)
            taboo[i] = false;
        m_lat->setVecNorm(m_ns, i);
        for (j = 0; j < dim; j++) {
            NTL::matrix_row<IntMat> row1(m_lat->getDualBasis(), j);
            if (i != j && m_nv[j] != 0) {
                NTL::conv(m_bs, -m_nv[j]);
                NTL::matrix_row<IntMat> row2(m_lat->getDualBasis(), i);
                std::cout << " redDual, before modifVect \n";
                ModifVect(row1, row2, m_bs, dim);
                //  ModifVect (m_lat->getDualBasis ()[j], m_lat->getDualBasis ()[i],
                //            m_bs, dim);
                m_lat->setDualNegativeNorm(j);
                if (taboo)
                    taboo[j] = false;
            }
        }
    }
    std::cout << " redDieter, end of pairwiseRedDual, i = " << i << "  *****\n";
}

//=========================================================================

template<typename Int, typename Real>
void ReducerBB<Int, Real>::redDieter(int64_t d, bool taboo[]) {
    std::int64_t BoundCount;
    const int64_t dim = m_lat->getDim();
    bool withDual = false;   // m_lat->withDual();
    //std::cout << " redDieter, withDual = " << withDual << "\n";

    m_lat->updateScalL2Norm(d, dim);
    m_lat->sortBasis(d);
    int64_t i = dim - 1;
    m_cpt = 0;
    m_countDieter = 0;
    BoundCount = 2 * dim - d;
    do {
        pairwiseRedPrimal(i, d, taboo);
        if (i > d && withDual) {
            pairwiseRedDual(i, taboo);
        }
        if (i < 1)
            i = dim - 1;
        else
            --i;
    } while (!(m_countDieter >= BoundCount || m_cpt > MAX_PRE_RED));
    std::cout << " End of redDieter \n";
}

//=========================================================================

template<typename Int, typename Real>
void ReducerBB<Int, Real>::redDieterRandomized(int64_t d, int64_t seed) {
    std::int64_t BoundCount;
    const int64_t dim = m_lat->getDim();
    bool withDual = false;  // m_lat->withDual();

    m_lat->updateScalL2Norm(d, dim);
    //m_lat->getDualBasis ().updateScalL2Norm (d, dim);
    m_lat->sortBasis(d);
    int64_t i = dim - 1;
    m_cpt = 0;
    m_countDieter = 0;
    BoundCount = 2 * dim - d;
    srand(seed);   // We use a simple RNG from the standard C library.rimal
    do {
        pairwiseRedPrimal(rand() % dim, d);
        if (i > d && withDual)
            pairwiseRedDual(i);
        if (i < 1)
            i = dim - 1;
        else
            --i;
    } while (!(m_countDieter >= BoundCount || m_cpt > MAX_PRE_RED));
}

//=========================================================================

/**
 * Reduce the Cholesky matrix with adding a multiple of the i-th vector
 * to the j-th vector. It updates the Gram-Schmidt matrix.
 */
template<typename Int, typename Real>
void ReducerBB<Int, Real>::reductionFaible(int64_t i, int64_t j) {
    /*
     * Reduit la matrice de Cholesky (d'ordre 2) en ajoutant un multiple du
     * vecteur i au vecteur j, si possible.  Modifie le vecteur dual W_i en
     * consequence et remet a jour la matrice des produits scalaires.
     * Utilise seulement par redLLLOld.
     */
    Real cte;
    std::int64_t cteLI;
    // bool withDual = false;  // m_lat->withDual();
    cte = m_cho2[i][j] / m_cho2[i][i];

    const int64_t dim = m_lat->getDim();

    if (abs(cte) < std::numeric_limits<double>::max()) {
        // On peut representer cte en int64_t
        if (abs(cte) > 0.5) {
            NTL::conv(cteLI, Round(cte));
            NTL::matrix_row<IntMat> row1(m_lat->getBasis(), j);
            NTL::matrix_row<IntMat> row2(m_lat->getBasis(), i);
            ModifVect(row1, row2, -cteLI, dim);
            //if(withDual){
            //  NTL::matrix_row<IntMat> row3(m_lat->getDualBasis(), i);
            //  NTL::matrix_row<IntMat> row4(m_lat->getDualBasis(), j);
            //  ModifVect (row3, row4, cteLI, dim);
            //}

        } else
            return;

    } else {
        // On represente cte en double.
        if (abs(cte) < std::numeric_limits<long double>::max())
            cte = Round(cte);
        NTL::matrix_row<IntMat> row1(m_lat->getBasis(), j);
        NTL::matrix_row<IntMat> row2(m_lat->getBasis(), i);
        ModifVect(row1, row2, -cte, dim);
        //if(withDual){
        //  NTL::matrix_row<IntMat> row3(m_lat->getDualBasis(), i);
        //  NTL::matrix_row<IntMat> row4(m_lat->getDualBasis(), j);
        //  ModifVect (row3, row4, cte, dim);
        //}

    }
    m_lat->setNegativeNorm(j);
    m_lat->updateVecNorm(j);
    //if(withDual){
    //  m_lat->setDualNegativeNorm (i);
    //  m_lat->updateDualVecNorm(i);
    //}
    miseAJourGramVD(j);
    calculCholesky2LLL(i, j);
}

// =========================================================================

template<typename Int, typename Real>
void ReducerBB<Int, Real>::transformStage3Mink(std::vector<std::int64_t> &z,
        int64_t &k) {
    int64_t j, i;
    std::int64_t q;
    const int64_t dim = m_lat->getDim();
    bool withDual = false;  // m_lat->withDual();

    j = dim - 1;
    while (z[j] == 0)
        --j;
    while (abs(z[j]) > 1) {
        i = j - 1;
        while (z[i] == 0)
            --i;
        // On a 2 indices i < j tels que |z_j| > 1 et z_i != 0.
        while (z[j]) {
            // Troncature du quotient vers 0
            q = z[i] / z[j];
            if (q) {
                // On ajoute q * v[i] au vecteur m_lat->getBasis()[j]
                z[i] -= q * z[j];
                NTL::matrix_row<IntMat> row2(m_lat->getBasis(), i);
                NTL::matrix_row<IntMat> row1(m_lat->getBasis(), j);
                //    ModifVect (m_lat->getBasis ()[j], m_lat->getBasis ()[i],
                //            q, dim);
                ModifVect(row1, row2, q, dim);
                m_lat->setNegativeNorm(j);

                if (withDual) {
                    NTL::matrix_row<IntMat> row3(m_lat->getDualBasis(), i);
                    NTL::matrix_row<IntMat> row4(m_lat->getDualBasis(), j);
                    //    ModifVect (m_lat->getDualBasis ()[i], m_lat->getDualBasis ()[j],
                    //             -q, dim);
                    ModifVect(row3, row4, -q, dim);
                    m_lat->setDualNegativeNorm(i);
                }
            }
            // Permutation.
            std::swap(z[i], z[j]);
            m_lat->permute(i, j);
        }
        j = i;
    }
    k = j;
}

//=========================================================================

/*
 * We call this when we have found a shorter vector which is a linear combination of the
 * previous basis vectors with coefficients given in z.  This procedure updates the basis
 * so that the new shortest vector in first and the other vectors are adjusted accordingly.
 * Note that an alternative for the L2 norm might be to add the new shortest vector to
 * the basis to get dim+1 generating vectors, and apply LLL to recover a new basis.
 */
template<typename Int, typename Real>
void ReducerBB<Int, Real>::transformStage3ShortVec(std::vector<std::int64_t> &z,
        int64_t &k) {
    int64_t j, i;
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
            // Troncature du quotient vers 0
            q = z[i] / z[j];
            if (q) {
                // On ajoute q * v[i] au vecteur m_lat->getBasis()[j]
                z[i] -= q * z[j];
                NTL::matrix_row<IntMat> row2(m_lat->getBasis(), i);
                NTL::matrix_row<IntMat> row1(m_lat->getBasis(), j);
                ModifVect(row1, row2, q, dim);
                m_lat->setNegativeNorm(j);
            }
            // Permutation.
            std::swap(z[i], z[j]);
            m_lat->permute(i, j);
        }
        j = i;
    }
    k = j;
}

//=========================================================================

template<typename Int, typename Real>
bool ReducerBB<Int, Real>::tryZMink(int64_t j, int64_t i, int64_t Stage,
        bool &smaller, const IntMat &WTemp)
        // If m_countNodes > MaxNodesBB, returns false, otherwise returns true.
        {
    std::int64_t max0, min0;
    Real x, dc;
    Real center;
    std::int64_t zhigh, zlow, h;
    bool high;
    int64_t k;
    Real S1, S2, S3, S4, mR;

    const int64_t dim = m_lat->getDim();
    NTL::conv(mR, m_lat->getModulo());

    ++m_countNodes;
    if (m_countNodes > maxNodesBB) {
        std::cerr << "-------- m_countNodes > maxNodesBB = " << maxNodesBB
                << std::endl;
        return false;
    }

    // Calcul d'un intervalle contenant les valeurs admissibles de zj.
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
        if (x > 0.0)
            ++min0;

        x = center + dc;
        NTL::conv(max0, trunc(x));
        if (x < 0.0)
            --max0;

        // En vue du choix initial de zj. On determine zlow et zhigh.
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
                NTL::matrix_row<const IntMat> row1(m_lat->getBasis(), dim - 1);
                m_bv = row1;
                for (k = 0; k < dim - 1; k++) {
                    if (m_zLI[k] != 0) {
                        NTL::matrix_row<const IntMat> row1(m_lat->getBasis(),
                                k);
                        ModifVect(m_bv, row1, m_zLI[k], dim);
                    }
                }
                if (Stage == 3) {
                    NTL::matrix_row<const IntMat> row1(m_lat->getBasis(),
                            dim - 1);
                    ModifVect(m_bv, row1, m_zLR[dim - 1] - 1.0, dim);
                }

                ProdScal<Int>(m_bv, m_bv, dim, S1);
                NTL::conv(S4, m_lat->getVecNorm(dim - 1));
                if (S1 < S4) {
                    if (Stage == 2) {
                        smaller = true;
                        if (!PreRedLLLMink)
                            m_zShort = m_zLI;
                        else {
                            for (k = 1; k < dim; k++) {
                                NTL::matrix_row<const IntMat> row1(WTemp, k);
                                ProdScal<Int>(m_bv, row1, dim, S2);
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
                            NTL::matrix_row<const IntMat> row1(WTemp, k);
                            ProdScal<Int>(m_bv, row1, dim, S2);
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
                if (!tryZMink(j - 1, i, Stage, smaller, WTemp))
                    return false;
                // Des qu'on a trouve quelque chose, on sort de la recursion
                // et on retourne dans reductMinkowski.
                if (smaller)
                    return true;
            }
        }
        if (high) {
            ++zhigh;
            if (zlow >= min0)
                high = false;
        } else {
            --zlow;
            if (zhigh <= max0)
                high = true;
        }
    }
    return true;
}

//=========================================================================

template<typename Int, typename Real>
bool ReducerBB<Int, Real>::redBBMink(int64_t i, int64_t d, int64_t Stage,
        bool &smaller, bool taboo[]) {
    /*
     * Tries to shorten m_lat->getBasis()[i] using branch-and-bound.
     * Used in Minkowski Reduction.
     * Stage is 2 or 3.
     * z[i] = 1 if Stage = 2, z[i] >= 2 if Stage = 3.
     * Stops and returns false if not finished after examining MaxNodesBB
     * nodes in the branch-and-bound tree.  When succeeds, returns true.
     * Assumes that the norm is Euclidean.
     */
    bool withDual = false;  // m_lat->withDual();
    const int64_t dim = m_lat->getDim();
    IntMat VTemp(dim, dim), WTemp(dim, dim);
    bool TabooTemp[dim];
    Real tmp;
    smaller = false;

    // Approximation du carre de la longueur de Vi.
    if (m_lat->getVecNorm(i) < 0) {
        NTL::matrix_row<IntMat> row1(m_lat->getBasis(), i);
        ProdScal<Int>(row1, row1, dim, tmp);
        //  ProdScal<Int> (m_lat->getBasis()[i], m_lat->getBasis()[i],
        //            dim, tmp);
        m_lat->setVecNorm(tmp, i);
    }
    NTL::conv(m_lMin2, m_lat->getVecNorm(i));

    if (Stage == 3 && withDual) {
        if (m_lat->getDualVecNorm(i) < 0) {
            NTL::matrix_row<IntMat> row1(m_lat->getDualBasis(), i);
            ProdScal<Int>(row1, row1, dim, tmp);
            //   ProdScal<Int> (m_lat->getDualBasis()[i], m_lat->getDualBasis()[i],
            //            dim, tmp);
            m_lat->setDualVecNorm(tmp, i);
        }
        Real m2;
        NTL::conv(m2, m_lat->getModulo());
        m2 = m2 * m2;
        if (m_lMin2 * m_lat->getDualVecNorm(i) < 4 * m2)
            return true; // if the angle between the basis vector i and the dual
        // basis vector i is between -Pi/3 and Pi/3
    }
    if (withDual) {
        m_lat->updateDualVecNorm();
    }
    m_lat->updateVecNorm();
    m_lat->permute(i, dim - 1);

    int64_t k, h;

    if (PreRedLLLMink) {
        // On memorise la base courante.
        VTemp = m_lat->getBasis();
        if (withDual) {
            WTemp = m_lat->getDualBasis();
        }
        for (h = 0; h < dim; h++)
            TabooTemp[h] = true;
        redLLLOld(1.0, 1000000, dim - 1);
        m_lat->updateVecNorm();
        if (withDual) {
            m_lat->updateDualVecNorm();
        }
    }
    if (!calculCholesky(m_dc2, m_c0))
        return false;
    m_countNodes = 0;
    m_n2[dim - 1] = 0.0;
    if (!tryZMink(dim - 1, i, Stage, smaller, WTemp))
        return false;

    if (PreRedLLLMink) {
        /* On remet l'anciennne base, celle d'avant LLL, avant de considerer
         la prochaine m_lat->dimension.  */
        m_lat->getBasis() = VTemp;
        m_lat->updateVecNorm();
        if (withDual) {
            m_lat->getDualBasis() = WTemp;
            m_lat->updateDualVecNorm();
        }
        for (h = 0; h < dim; h++)
            taboo[h] = TabooTemp[h];
    }
    if (smaller) {
        /* On a trouve un plus court vecteur.  On ameliore
         m_lat->getBasis()[k].  */
        if (Stage == 2)
            k = dim - 1;
        else
            transformStage3Mink(m_zShort, k);
        NTL::matrix_row<IntMat> row1(m_lat->getBasis(), k);
        for (h = 0; h < dim; h++)
            row1(h) = m_bv[h];
        //  m_lat->getBasis ()[k] = m_bv;
        m_lat->setNegativeNorm(k);
        if (m_zShort[k] < 0 && withDual) {
            NTL::matrix_row<IntMat> row2(m_lat->getDualBasis(), k);
            ChangeSign(row2, dim);
        }
        /* Mise a jour des vecteurs de la base duale selon le nouveau
         m_lat->getBasis()[k] */
        for (h = 0; h < dim; h++) {
            if ((m_zShort[h] != 0) && (h != k)) {
                if (withDual) {
                    NTL::matrix_row<IntMat> row1(m_lat->getDualBasis(), h);
                    NTL::matrix_row<IntMat> row2(m_lat->getDualBasis(), k);
                    ModifVect(row1, row2, -m_zShort[h], dim);
                    m_lat->setDualNegativeNorm(h);
                }
                if (Stage == 2) {
                    if (h >= d)
                        taboo[h] = false;
                }
            }
        }
    } else if (Stage == 2)
        taboo[dim - 1] = true;

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
     * meaning that a shortest vector has been found.
     * The square length of the current shortest vector is maintained in `m_lMin2`.
     * When j=1 and we find a shorter vector, we return the corresponding `z_j`'s in `m_zShort`.
     */

    /* Pour une implantation non recursive, ces variables devraient
     etre des tableaux indices par j. */
    Real dc, x, center, mn_xsquare_md;
    std::int64_t min0, max0;     // Bornes de l'intervalle pour les z_j.
    std::int64_t zlow, zhigh; // Valeur courante a gauche et a droite du centre.
    bool high; // Indicates if we are on the right (true) or the left of the center.
    int64_t k;
    std::int64_t temp;
    const int64_t dim = m_lat->getDim();
    ++m_countNodes;
    if (m_countNodes > maxNodesBB) {
        std::cerr << "*****   m_countNodes > maxNodesBB = " << maxNodesBB
                << std::endl;
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
    if (!m_foundZero)
        min0 = 0;     // We are at the beginning, min will be zero.
    else {
        x = center - dc;
        NTL::conv(min0, trunc(x));
        if (x > 0.0)
            ++min0;
    }
    x = center + dc;
    NTL::conv(max0, trunc(x));
    if (x < 0.0)
        --max0;

    //   std::cout << "Borne Min="<<min0<<"     Borne Max="<<max0<<"     j="<<j<<std::endl;
    // Compute initial values of zlow and zhigh, the search pointers on each side of the interval.
    if (min0 > max0)
        return true;
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
        if (high)
            m_zLI[j] = zhigh;
        else
            m_zLI[j] = zlow;
        m_zLR[j] = m_zLI[j];

        // Computing m_n2[j-1].
        x = m_zLR[j] - center;
        mn_xsquare_md = m_n2[j] + x * x * m_dc2[j]; // We pre-compute this to save time.
        if (j == 0) {
            // All the zj have been selected: we have a vector to test.
            if (m_lMin2 > mn_xsquare_md) {
                /* Check if we have a shorter nonzero vector. */
                if (!m_foundZero) {
                    // Le premier vecteur trouve sera zero.
                    m_foundZero = true;
                } else {
                    SetZero(m_bv, dim);
                    for (k = 0; k < dim; k++) {
                        if (m_zLI[k] != 0) {
                            NTL::matrix_row<IntMat> row1(m_lat->getBasis(), k);
                            ModifVect(m_bv, row1, m_zLI[k], dim);
                        }
                    }
                    // The new vector is now in m_bv.
                    if (m_lat->getNormType() == L2NORM) {
                        ProdScal<Int>(m_bv, m_bv, dim, x);
                    } else {
                        // Compute the square length for the L1 norm.
                        CalcNorm<IntVec, Real>(m_bv, dim, x,
                                m_lat->getNormType());
                        x = x * x;
                    }
                    if (x < m_lMin2) {
                        // The new vector is shorter.
                        smaller = true;
                        NTL::conv(m_lMin2, x);
                        m_zShort = m_zLI;
                        m_bw = m_bv;   // ??????
                    }
                }
            }
        } else if (m_lMin2 > mn_xsquare_md) {
            // There is still hope; we continue the recursion.
            m_n2[j - 1] = mn_xsquare_md;
            if (!tryZShortVec(j - 1, smaller, norm))
                return false;
        } else
            m_n2[j - 1] = mn_xsquare_md;
        if (high) {
            ++zhigh;
            if (zlow >= min0)
                high = false;
        } else {
            --zlow;
            if (zhigh <= max0)
                high = true;
        }
    }
    return true;
}

//=========================================================================

template<typename Int, typename Real>
bool ReducerBB<Int, Real>::redBBShortVec() {
    /*
     * Finds shortest non-zero vector, using branch-and-bound, with L1 or L2 norm.
     * Stops and returns false if not finished after examining MaxNodesBB nodes in the
     * branch-and-bound tree.  When succeeds, returns true, and the squared shortest
     * vector length will be in m_lMin2.
     *
     * This function uses (directly or indirectly) the following class variables:
     *    m_lMin, m_lMin2, m_decomp, m_boundL2, m_n2, m_countNodes, m_foundZero,
     *    m_bw, m_zShort, m_c0, m_zLR, m_zLI, m_dc2, ....  and more.
     * From m_lat (current lattice object):
     *    norm, sortBasisNoDual, updateScalL2Norm, getBasis,
     *
     *
     */
    NormType norm = m_lat->getNormType();
    if ((norm != L1NORM) & (norm != L2NORM)) {
        std::cerr << "RedBBShortVec: only L1 and L2 norms are supported";
        return false;
    }
    const int64_t dim = m_lat->getDim();  // Lattice dimension
    //std::cout << " Start redBBShortVec, dim  = " << dim << "\n";
    //std::cout << " Start redBBShortVec, basis = \n" << m_lat->getBasis() << "\n";

    bool smaller = false;  // Will change when we find a smaller vector.
    int64_t k, h;
    Real x;
    //  if (norm != L2NORM)  m_lat->setNegativeNorm();
    // Here we sort the basis by L2 lengths, otherwise Cholesky will fail more rapidly
    // due to floating-point errors.
    m_lat->updateScalL2Norm(0, dim);
    m_lat->sortBasis(0);

    // Approximate the square norm of the current shortest vector.
    if (norm == L2NORM) {
        NTL::conv(m_lMin2, m_lat->getVecNorm(0));
    } else {
        // Looking for the shortest vector in basis according to the considered norm
        NTL::matrix_row<IntMat> row1(m_lat->getBasis(), 0);
        CalcNorm<IntVec, Real>(row1, dim, m_lMin, norm);
        for (k = 1; k < dim; k++) {
            NTL::matrix_row<IntMat> row2(m_lat->getBasis(), k);
            CalcNorm<IntVec, Real>(row2, dim, x, norm);
            if (x < m_lMin)
                m_lMin = x;
        }
        m_lMin2 = m_lMin * m_lMin;  // Squared shortest length with L1 norm.
    }
    // std::cout << " redBBShortVec, Sq length shortest = " <<  m_lMin2 << "\n";

    // If we already have a shorter vector than the minimum threshold, we stop right away.
    // This is useful for the seek programs in LatMRG.
    if (m_lMin2 <= m_BoundL2[dim - 1])
        return false;

    if (m_decomp == CHOLESKY) {
        /* Perform the Cholesky decomposition; if it fails we exit. */
        if (!calculCholesky(m_dc2, m_c0))
            return false;
    } else if (m_decomp == TRIANGULAR) {  // Just for testing; this is very slow!
        /* Perform a triangular decomposition:
         * Instead of doing the following, perhaps we may just assume that the basis
         * is already lower triangular?
         */
        IntMat m_v, m_v2;   // Here we create new matrices each time!!!!
        m_v.resize(dim, dim);
        m_v2.resize(dim, dim);
        Int mod = m_lat->getModulo();
        CopyMatr(m_v, m_lat->getBasis(), dim, dim);
        lowerTriangularBasis (m_v, m_v2, mod);
        // CopyMatr(m_lat->getBasis(), m_v2, dim, dim);
        for (int64_t i = 0; i < dim; i++) {
            for (int64_t j = 0; j < dim; j++) {
                if (i != j)
                    m_c0[i][j] = NTL::conv < Real > (m_v2[i][j]) / NTL::conv
                            < Real > (m_v2[i][i]);
                else
                    m_c0[i][j] = NTL::conv < Real > (m_v2[i][j]);
            }
        }
        for (int64_t i = 0; i < dim; i++) {
            m_dc2[i] = m_c0[i][i] * m_c0[i][i];
        }
    } else {
        std::cerr << "RedBBShortVec:decomp value not supported";
        return false;
    }

    //std::cout << " redBBShortVec, after decomp \n";
    /* Perform the branch and bound.  */
    /* m_n2[j] will be the sum of terms |z*k|^2 ||v*k||^2 for k > j.  */
    m_n2[dim - 1] = 0.0;
    m_countNodes = 0;
    smaller = false;
    m_foundZero = false;
    if (!tryZShortVec(dim - 1, smaller, norm)) // We search for a shortest vector.
        return false;
    //std::cout << " redBBShortVec, after tryShortVec \n";
    if (smaller) {
        //std::cout << " redBBShortVec, found a shorter vector \n";
        // We found a shorter vector. Its square length is in m_lMin2.
        transformStage3ShortVec(m_zShort, k); // Is this useful and OK for L1 ???
        //std::cout << " redBBShortVec, after transformStage 3 \n";
        NTL::matrix_row<IntMat> row1(m_lat->getBasis(), k);
        //std::cout << " redBBShortVec, after row1, k = " << k << "\n";
        //std::cout << " redBBShortVec, row1 = " << row1 << "\n";
        //std::cout << " redBBShortVec, m_bw = " << m_bw << "\n";
        //std::cout << " redBBShortVec, m_bw[1] = " << m_bw[1] << "\n";
        for (h = 0; h < dim; h++) {
            //std::cout << " redBBShortVec, row1 = m_bw, h = " << h << ", m_bw[h] = " << m_bw[h] << ", m_bv[h] = " << m_bv[h] << "\n";
            // row1(h) = m_bw[h];   ??????????????
            row1[h] = m_bw[h];
        }
        //std::cout << " redBBShortVec, before set negative norm (k) \n";
        m_lat->setNegativeNorm(k);
        //std::cout << " redBBShortVec, after set negative norm (k) \n";

        /* The new current shortest vector will be in
         m_lat->getBasis()(0). */
        /* In the case of L1NORM, we must check if it is really smaller.  */
        if (norm == L2NORM)
            m_lat->permute(k, 0);
        else {
            NTL::matrix_row<IntMat> row5(m_lat->getBasis(), k);
            CalcNorm(row5, dim, x, norm);
            if (x < m_lMin) {
                m_lMin = x;
                m_lMin2 = m_lMin * m_lMin;
                m_lat->permute(k, 0);
            }
        }
    }
    m_lat->updateVecNorm();
    m_lat->sortBasis(0);
    //std::cout << " redBBShortVec, before exit \n";
    return true;
}

//=========================================================================

// This works only for the L2 norm.
template<typename Int, typename Real>
bool ReducerBB<Int, Real>::reductMinkowski(int64_t d) {
    bool withDual = false;  //  m_lat->withDual();
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
            redDieter(d, taboo);  // Not very efficient...  Change to BKZ.
            m_lat->setNegativeNorm(d);
            m_lat->updateVecNorm(d);
            if (withDual) {
                m_lat->setDualNegativeNorm(d);
                m_lat->updateDualVecNorm(d);
            }
            m_lat->sortBasis(d);
            found = false;

            for (i = 0; i < dim; i++) {
                if (!taboo[i]) {
                    // On essaie de reduire le i-eme vecteur.
                    if (!redBBMink(i, d, 2, smaller, taboo))
                        return false;
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
                if (!redBBMink(i, d, 3, smaller, taboo))
                    return false;
                totalNodes += m_countNodes;
                if (smaller)
                    found = true;
            }
        }
    } while (found);
    m_lat->setNegativeNorm();
    m_lat->updateScalL2Norm(0, dim);
    if (withDual) {
        m_lat->setDualNegativeNorm();
        m_lat->updateDualScalL2Norm(0, dim);
    }
    m_lMin2 = m_lat->getVecNorm(0);
    return true;
}

//=========================================================================

template<typename Int, typename Real>
void ReducerBB<Int, Real>::tracePrintBases(char *message) {
    std::cout << std::endl << "================================= " << message
            << std::endl;
    //std::cout << "dim = " << m_lat->getDim () << std::endl;
    m_lat->setNegativeNorm();
    //m_lat->setDualNegativeNorm();
    m_lat->updateVecNorm();
    //m_lat->updateDualVecNorm();
    m_lat->sortBasis(0);
    m_lat->write();
}

//=========================================================================

// Translation of our old LLL implementation from Modula-2 to C++. Rather slow.
template<typename Int, typename Real>
void ReducerBB<Int, Real>::redLLLOld(double delta, std::int64_t maxcpt,
        int64_t Max) {
    // This is a C translation of the old implementation by R. Couture.
    // Effectue la pre-reduction de B au sens de Lenstra-Lenstra-Lovasz. N'utilise
    // pas les vecteurs m_lat->getBasis().vecNorm et  Wm_lat->getDualBasis().

    // bool withDual = false;   // m_lat->withDual();
    const int64_t REDBAS_e = 40;
    int64_t i, j, k, h;
    Real Cho0ij;
    Real limite;
    std::int64_t cpt;

    const int64_t dim = m_lat->getDim();
    if (Max == 0)
        Max = dim;
    cpt = 0;
    calculGramVD();
    limite = 1.0;
    for (k = 1; k <= REDBAS_e; k++)
        limite *= 2.0;
    limite *= dim;
    m_cho2[0][0] = m_gramVD[0][0];
    m_cho2[0][1] = m_gramVD[0][1];
    m_IC[0] = 1;
    m_cho2[1][1] = m_gramVD[1][1]
            - m_cho2[0][1] * (m_cho2[0][1] / m_cho2[0][0]);
    m_IC[1] = 1;
    for (i = 2; i < dim; i++)
        m_IC[i] = -1;
    h = 0;

    while (h < Max - 1 && cpt < maxcpt) {
        if (m_gramVD[h + 1][h + 1] > limite) {
            for (i = h; i >= 0; i--)
                reductionFaible(i, h + 1);
        } else
            reductionFaible(h, h + 1);

        calculCholesky2Ele(h + 1, h + 1);
        if (m_IC[h + 1] == -1)
            m_IC[h + 1] = h + 1;
        if (m_cho2[h + 1][h + 1] / m_cho2[h][h]
                + (m_cho2[h][h + 1]) / m_cho2[h][h]
                        * (m_cho2[h][h + 1] / m_cho2[h][h]) < delta) {
            ++cpt;
            m_lat->permute(h, h + 1);
            permuteGramVD(h, h + 1, dim);
            m_cho2[h][h] = m_gramVD[h][h];
            for (i = 0; i < h; i++) {
                std::swap(m_cho2[i][h], m_cho2[i][h + 1]);
                m_cho2[h][h] -= m_cho2[i][h] * (m_cho2[i][h] / m_cho2[i][i]);
            }
            if (h == 0) {
                Cho0ij = m_cho2[0][1] / m_cho2[0][0];
                if (abs(Cho0ij) > 0.5) {
                    m_IC[0] = 1;
                    m_IC[1] = -1;
                    h = 0;
                } else {
                    m_cho2[1][1] = m_gramVD[1][1]
                            - m_cho2[0][1] * m_cho2[0][1] / m_cho2[0][0];
                    calculCholesky2LLL(2, 2);
                    m_IC[0] = 2;
                    m_IC[1] = 2;
                    m_IC[2] = 2;
                    h = 1;
                }
            } else {
                m_IC[h] = h + 1;
                m_IC[h + 1] = -1;
                --h;
            }
        } else {
            for (i = 0; i <= h + 2; i++) {
                if (h + 2 > m_IC[i]) {
                    if (h + 2 < dim)
                        calculCholesky2Ele(i, h + 2);
                    m_IC[i] = h + 2;
                }
            }
            ++h;
        }
    }

    if (cpt == maxcpt) {
        std::cout << "***** in redLLLOld cpt > maxcpt = " << maxcpt
                << std::endl;
    }
    for (j = 2; j < Max; j++) {
        for (i = j - 2; i >= 0; i--)
            reductionFaible(i, j);
    }
    m_lat->setNegativeNorm();
}

//=========================================================================

template<typename Int, typename Real>
bool ReducerBB<Int, Real>::reductMinkowski(IntLattice<Int, Real> &lat,
        int64_t d) {
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
