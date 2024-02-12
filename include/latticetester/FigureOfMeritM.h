// This file is part of LatticeTester.
//
// Copyright (C) 2012-2023  The LatticeTester authors, under the supervision
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

#ifndef LATTICETESTER_FIGUREOFMERITM_H
#define LATTICETESTER_FIGUREOFMERITM_H

#include <iostream>
#include <cstdint>

#include "NTL/tools.h"
#include "NTL/ZZ.h"
#include "NTL/vec_ZZ.h"
#include "NTL/mat_ZZ.h"
#include "NTL/LLL.h"
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/RR.h>

#include "latticetester/FlexTypes.h"
#include "latticetester/EnumTypes.h"
#include "latticetester/BasisConstruction.h"
#include "latticetester/Util.h"
#include "latticetester/IntLattice.h"
#include "latticetester/IntLatticeExt.h"
#include "latticetester/LLL_FPZZflex.h"
#include "latticetester/Rank1Lattice.h"
#include "latticetester/Reducer.h"
#include "latticetester/Normalizer.h"
#include "latticetester/CoordinateSets.h"

namespace LatticeTester {

/**
 * This class offers functions to calculate the figure of merit (FOM)
 * \f{equation}{
 *    M_{t_1,\dots,t_d} = \min\left[ \min_{I\in S(t_1)} \omega_I \ell_I/\ell_I^*(\eta_I),\;
 *    \min_{2\le s\le d}\, \min_{I\in S(s,t_s)} \omega_I \ell_I/\ell_I^*(\eta_I) \right],
 * \f}
 * defined in Section 10 of the guide, for any given `IntLatticeExt` object.
 * The FOM is computed only for the (rescaled) primal lattice, the m-dual is never used.
 * The projections in \f$S(t_1)$\f are those over successive coordinates in up to \f$t_1$\f
 * dimensions, while those is  @f$S(s,t_s)$@f are projections over \f$s$\f distinct coordinates
 * that are non necessarily successive and are all in the set \f$\{1,\dots,t_s\}$\f,
 * for each order @f$s > 1$@f.
 * There are two variants for the latter: the first (default) variant takes
 * @f$S(s,t_s)$@f as just defined (also defined in the guide), and the other considers
 * only the set @f$S^{(1)}(s,t_s)$@f of projections that contain coordinate 1.
 * The parameter `includeFirst` in the constructor determines which variant is taken.
 *
 * The lengths of the shortest vectors in the projections can be calculated exactly by using the
 * BB algorithm after applying some pre-reduction, or they can be just approximated
 * by the lengths of the shortest basis vector obtained after applying some pre-reduction
 * such as LLL or BKZ.  The latter is much faster but not exact.
 * The `ReductionType` parameter in the constructor selects the method.
 * Note that we need a `Reducer` object only when BB is applied.
 * I other cases, we just use static methods for the reduction.
 *
 * The constructor requires the vector @f$(t_1,\dots,t_d)$@f, the type of reduction
 * that will be used to compute or approximate the vector lengths,
 * a `Reducer` object used for the reduction in case the reduction method includes BB,
 * a `Normalizer` object used to normalize the merit values, and the optional
 * `includeFirst` parameter in case we want to put it to `true`.
 * The function `setTVector` permits one to set (or reset) the vector  @f$(t_1,\dots,t_d)$@f.
 * One can use `setReductionMethod` to reset the type of reduction or change its parameters
 * `delta` and `blocksize` from their default values of 0.99999 and 10.
 * The `Normalizer` can be changed via `setNormalizer`.
 *
 * The method `computeMerit` computes the FOM.
 * The computation is stopped (early exit) as soon as we know that the value of the FOM
 * will be outside the interval `[low, high]`.  By default,
 * these two bounds are 0 and 1, but they can be changed via the function `setBounds`.
 *
 * Note: The class works only for the case where  "PrecisionType == DOUBLE".   
 * This is a limitation.
 */

template<typename Int> class FigureOfMeritM {

private:
    typedef NTL::vector<Int> IntVec;
    typedef NTL::matrix<Int> IntMat;
public:

    /*
     * This constructor will call `setTVector` with the given vector `t`
     * and `includeFirst` variable,
     * then set the `ReductionType`, `Reducer`, and `Normalizer` to the given values.
     */
    FigureOfMeritM(const vector<int64_t> &t, ReductionType &meth,
            Reducer<Int, Real> &red, Normalizer &norma, bool includeFirst = false);

    /*
     * Sets the vector @f$(t_1,..., t_d)$@f in the FOM definition to the vector `t`.
     * Note that the values of @f$t_1,..., t_d$@f are taken from `t[0],...,t[d-1]`,
     * respectively.  When `includeFirst` is `true`, we consider only the non-successive
     * projections that contain coordinate 1.
     * See the doc of the class `FromRanges` in `CoordinateSets` for more details.
     */
    void setTVector(const vector<int64_t> &t, bool includeFirst = false);

    /*
     * Sets the reduction method and the parameters delta and blocksize used
     * for the LLL and BKZ algorithm.
     */
    void setReductionMethod(ReductionType &meth, double delta = 0.99999,
            int64_t blocksize = 10);

    /*
     * Sets the normalizer to `norma`.
     */
    void setNormalizer(Normalizer &norma) {
        m_norma = &norma;
    }

    /*
     * Sets the low and the high bound for the FOM.
     * The FOM computation is stopped as soon as we know it is outside these bounds.
     */
    void setBounds(double &low, double &high) {
        m_highbound = high;
        m_lowbound = low;
    }

    /*
     * This function computes and returns the value of the FOM for the given lattice 'lat'.
     * The function returns 0 if the computation was not completed for some reason
     * (early exit, error, etc.).
     * The parameter `proj` points to an `IntLattice` object that is used to store the
     * projections when computing the FOM.  The `maxDim` in this object must be large
     * enough so it can store any of the projections: it must be at least \f$d$\f and
     * at least \f$t_1$\f.
     * Re-using this object permits one to avoid creating new objects internally.
     */
    double computeMerit(IntLatticeExt<Int, Real> &lat, IntLattice<Int, Real> *proj);

    /*
     * This function computes and returns the FOM for all projections
     * over sets of successive coordinates of the form
     * {1, 2, ..., m_t.size()} to {1, 2, ..., m_t[0]}.
     * It returns 0 if the computation was not completed for some reason.
     * The parameter `proj` is like for `computeMerit`.
     */
    double computeMeritSucc(IntLatticeExt<Int, Real> &lat,
            IntLattice<Int, Real> *proj);

    /*
     * This function computes and returns the FOM for all projections
     * over sets on non-successive coordinates determined by
     *  @f$S(s,t_s)$@f or  @f$S^{(1)}(s,t_s)$@f.
     * It returns 0 if the computation was not completed for some reason.
     * The parameter `proj` is like for `computeMerit`.
     */
    double computeMeritNonSucc(IntLatticeExt<Int, Real> &lat,
            IntLattice<Int, Real> *proj);

protected:

    /*
     * The vector 'm_t' defines the set of projections for which the FOM is computed.
     * `m_t[s-1]` represents t_s, for s=1,...,d.
     */
    vector<int64_t> m_t;

    /*
     * The reduction type used to compute (or estimate) the FOM.
     */
    ReductionType m_reductionMethod;

    /*
     * The `delta` parameter for LLL and BKZ reductions.
     */
    double m_delta = 0.99999;

    /*
     * Blocksize for the BKZ algorithm.
     */
    int64_t m_blocksize = 10;

    /*
     * Internal `CoordinateSets` object used to store the set of projections
     * that are considered for the non-successive coordinates.
     * It is constructed and populated by the `setTVector` function.
     * This object will contain a set of projections of each order.
     */
    CoordinateSets::FromRanges *m_coordRange;

    /*
     * Internal normalizer object for storing normalizing values in FiguresOfMeritM class
     */
    Normalizer *m_norma;

    /*
     * Internal reducer object used for finding the shortest vector of a projection
     */
    Reducer<Int, Real> *m_red;

    /*
     * Pointer to a vector used to store the square Euclidean lengths of the basis vectors
     * after an LLL or BKZ pre-reduction via the static methods of `LLL_FPZZflex`.
     */
    double *m_sqlen;
    /*
     * Indicates if the first coordinate will always be included in all the projections
     * over the non-successive coordinates.  If true, we use @f$S^{(1)}(s,t_s)$@f  in
     * Section 10 of the user's guide, otherwise we use @f$S(s,t_s)$@f and we therefore
     * have a larger set of projections.
     */
    // bool m_firstCoordinateAlwaysIn = false;

    /*
     * As soon as we know the FOM is above this bound, its calculation stopped.
     */
    double m_highbound = 1.0;

    /*
     * As soon as we know the FOM is below this bound, its calculation stopped.
     */
    double m_lowbound = 0.0;

    /*
     * The next three variables below are introduced just to simplify the code.
     * This one indicates whether the BB algorithm will be performed or not.
     */
    bool m_doingBB;

    /*
     * Indicates whether the LLL algorithm is performed.
     */
    bool m_doingLLL;

    /*
     * Indicates whether the BKZ algorithm is performed.
     */
    bool m_doingBKZ;
    
};

//============================================================================
// Implementation

template<typename Int>
FigureOfMeritM<Int>::FigureOfMeritM(const vector<int64_t> &t,
        ReductionType &meth, Reducer<Int, Real> &red, Normalizer &norma, bool includeFirst) {
    setTVector(t, includeFirst);
    setReductionMethod(meth);
    setNormalizer(norma);
    m_red = &red;
}

//=========================================================================

template<typename Int>
void FigureOfMeritM<Int>::setTVector(const vector<int64_t> &t, bool includeFirst) {
    m_t = t;
    // m_firstCoordinateAlwaysIn = includeFirst;
    // Clear the CoordinateSets object.
    m_coordRange = new CoordinateSets::FromRanges;
    /* Defines the lower bound for the range of coordinates that belong to a projection.
     * It is 2 if the first coordinate belongs to all the projections, because we do
     * not have to consider coordinate 1. Otherwise it is 1.
     * In the first case, coordinate 1 is added to all the projections as an extra coord.
     * See the doc of the class `FromRanges` in `CoordinateSets`.
     */
    int64_t min_dim = 1;
    if (includeFirst) min_dim = 2;
    int64_t d = static_cast<int>(t.size()); // Number of orders for the projections.
    for (int64_t i = 1; i < d; i++)
        // Adds the set of projections of order i.
        m_coordRange->includeOrder(i+1, min_dim, t[i], includeFirst);  // Change here *****
}

//=========================================================================

template<typename Int>
void FigureOfMeritM<Int>::setReductionMethod(ReductionType &meth, double delta,
        int64_t blocksize) {
    m_reductionMethod = meth;
    if (m_reductionMethod == BKZBB || m_reductionMethod == LLLBB
            || m_reductionMethod == PAIRBB) {
        m_doingBB = true;
    } else {
        m_doingBB = false;
    }
    if (m_reductionMethod == BKZBB || m_reductionMethod == BKZ) {
        m_doingBKZ = true;
    } else {
        m_doingBKZ = false;
    }
    if (m_reductionMethod == LLLBB || m_reductionMethod == LLL) {
        m_doingLLL = true;
    } else {
        m_doingLLL = false;
    }
    m_delta = delta;
    m_blocksize = blocksize;
}

//=========================================================================

template<typename Int>
double FigureOfMeritM<Int>::computeMerit(IntLatticeExt<Int, Real> &lat,
        IntLattice<Int, Real> *proj) {
    double merit = 0;
    double minmerit = 1.0;

    this->m_sqlen = new double[this->m_t.size() + 1];
    merit = computeMeritNonSucc(lat, proj);
    if (merit < minmerit)
        minmerit = merit;
    // In any of these cases the calculation is stopped
    if (minmerit == 0)
        return 0;

    this->m_sqlen = new double[this->m_t[0]];
    merit = computeMeritSucc(lat, proj);
    if (merit < minmerit)
        minmerit = merit;
    // In any of these cases the calculation is stopped
    if (minmerit == 0 || merit > this->m_highbound)
        return 0;

    return minmerit;
}

//=========================================================================
template<typename Int>
double FigureOfMeritM<Int>::computeMeritSucc(
        IntLatticeExt<Int, Real> &lat, IntLattice<Int, Real> *proj) {
    double merit = 0;
    double minmerit = 1.0;
    int64_t lower_dim = static_cast<int64_t>(this->m_t.size());
    lat.buildBasis(lower_dim);
    for (int64_t j = lower_dim + 1; j < this->m_t[0] + 1; j++) {
        lat.incDimBasis();
        if (m_doingBKZ) {
            BKZ_FPZZflex(lat.getBasis(), this->m_delta, this->m_blocksize, 
                    j, j, this->m_sqlen);
        } else if (m_doingLLL) {
            LLL_FPZZflex(lat.getBasis(), this->m_delta, j, j, this->m_sqlen);
        } else if (this->m_reductionMethod == PAIRBB) {
            this->m_red->redDieter(0);
            this->m_sqlen[0] = lat.getVecNorm(0);
        }
        if (!m_doingBB) {
            if (lat.getNormType() == L2NORM) {
                NTL::conv(merit,
                        sqrt(this->m_sqlen[0]) / this->m_norma->getBound(j));
            } else {
                NTL::conv(merit,
                		this->m_sqlen[0] / this->m_norma->getBound(j));
            }
        } else {
            if (!m_red->shortestVector(lat))
                return 0;
            NTL::conv(merit,
                    m_red->getMinLength() / this->m_norma->getBound(j));
        }
        if (merit < minmerit)
            minmerit = merit;
        if (minmerit <= this->m_lowbound)
            return 0;
    }
    return minmerit;
}

//=========================================================================
template<typename Int>
double FigureOfMeritM<Int>::computeMeritNonSucc(
        IntLatticeExt<Int, Real> &lat, IntLattice<Int, Real> *proj) {
    double merit = 0;
    double minmerit = 1.0;
    Coordinates coord;

    for (auto it = m_coordRange->begin(); it != m_coordRange->end(); it++) {
        coord = *it;
        lat.buildProjection(proj, coord, this->m_delta);  // Must have withDual = false   ***
        if (m_doingBKZ) {
            BKZ_FPZZflex(proj->getBasis(), this->m_delta, this->m_blocksize, 
                    coord.size(), coord.size(), this->m_sqlen);
        } else if (m_doingLLL) {
            LLL_FPZZflex(proj->getBasis(), this->m_delta, coord.size(),
                    coord.size(), this->m_sqlen);
        } else if (this->m_reductionMethod == PAIRBB) {
            this->m_red->redDieter(0);
            this->m_sqlen[0] = lat.getVecNorm(0);
        }
        if (!m_doingBB) {
            if (lat.getNormType() == L2NORM) {
                NTL::conv(merit,
                        sqrt(this->m_sqlen[0])
                                / this->m_norma->getBound(coord.size()));
                std::cout << merit << "\n";
            } else {
                NTL::conv(merit,
                		this->m_sqlen[0]
                                / this->m_norma->getBound(coord.size()));
            }
        } else {
            if (!m_red->shortestVector(*proj))
                return 0;
            merit = NTL::conv<double>(
                    m_red->getMinLength()
                            / this->m_norma->getBound(coord.size()));
        }
        if (merit < minmerit)
            minmerit = merit;
        if (minmerit <= this->m_lowbound)
            return 0;
    }
    return minmerit;
}

template class FigureOfMeritM<NTL::ZZ> ;

} // end namespace LatticeTester

#endif
