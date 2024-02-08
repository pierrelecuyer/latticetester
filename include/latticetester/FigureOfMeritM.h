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

#include "latticetester/NormaBestLat.h"
#include "latticetester/NormaBestBound.h"
#include "latticetester/WeightsUniform.h"
#include "latticetester/CoordinateSets.h"

namespace LatticeTester {

/**
 * This class offers methods (functions) to calculate the figure of merit M for a given
 * IntLattice object. The FoM M is defined as the minimum of the normalized shortest
 * vector lengths for a set of projections which is determined by a vector (t_1,...,t_d).
 * These projections can consist of successive coordinates as well as non-successive coordinates.
 * The exact lengths of the shortest vectors can either be calculated exactly by using the 
 * BB algorithm or they can be approximated using the shortest base vector length after LLL, BKZ 
 * or pairwise pre-reduction has been applied, which is much faster but less exact. 
 * 
 * Moreoover, the bounds and the parameters of the reduction method also have default values.
 * These can be changed by the methods 'setReductionMethod' and 'setBounds'.
 * 
 * The main method of this class is 'computeMeritM' which calculates the actual figure of merit
 * M for a given lattice.
 *
 */

template<typename Int> class FigureOfMeritM {

private:
    typedef NTL::vector<Int> IntVec;
    typedef NTL::matrix<Int> IntMat;
public:

    /*
     * Constructor of a FigureOfMeritM object. Needs a reducer object
     * to be passed. The vector 't' defines the set of
     * dimensions for which the figure of merit is calculated. It contains
     * the values of t_1,..., t_d in the definition of the FoM.
     *
     * ****  I put `includeFirst` to have the same name as in CoordinateSets
     *       and the same default value, for consistency.
     */
    FigureOfMeritM(const vector<int64_t> &t, ReductionType &meth,
            Reducer<Int, Real> &red, Normalizer &norma, bool includeFirst = false);

    /*
     * Sets the vector 't' @f$=  (t_1,..., t_d)$@f in the definition of this FOM.
     * ******   Are the values in t[0],..., t[d-1]  or in t[1],...,t[d] ??????   *******
     *          It should be the first case.   This must be said.
     *          t[0] will contain $t_1$, etc.
     */
    void setTVector(const vector<int64_t> &t);

    /*
     * Sets the reduction method and the parameter delta for LLL / BKZ algorithm
     * as well as the blocksize used in the BKZ algorithm.
     */
    void setReductionMethod(ReductionType &meth, double delta = 0.99999,
            int64_t blocksize = 10);

    /*
     * Directly sets a new coordinate set.
     *   Why do we want this??  It may be inconsistent with the definition of M.  ******
     */
    void setCoordinateSets(CoordinateSets::FromRanges coordRange) {
        m_coordRange = coordRange;
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
     * Sets the normalizer to `norma`.
     */
    void setNormalizer(Normalizer &norma) {
        m_norma = &norma;
    }
    
    /* 
     * Sets a boolean variable that indicates if the first coordinate will always
     * be included in all the projections over the non-successive coordinates.
     * If true, the FOM will @f$S^{(1)}(s,t_s)$@f  in Section 10 of the user's guide,
     * otherwise it will use @f$S(s,t_s)$@f, which contains a larger number of projections.
     * The default value is `false`.
     *
     * ******   Problem:  If we call this after the constructor,
     *          m_coordRange  is no longer correct!   So this "set..." does not work.
     *          This variable must be set in the constructor.
     */
    void setFirstCoordinateAlwaysIn(bool &includeFirst) {
        m_firstCoordinateAlwaysIn = includeFirst;
    }

    /*
     * This function calculates the Figure of Merit M of a given lattice 'lat'
     * and should be called by the user. The vector 't' defines the set of
     * dimensions for which the figure of merit is calculated. It contains
     * the values of t_1,..., t_d in the definition of the FoM.
     * The IntLattice object 'proj' is needed for saving projections.
     * The value 0 is returned if an error occurs
     * while calculating the shortest vector.
     *
     * In practice, we may also want to recover which projection gave the worst `merit`.  ********
     * And in some cases, we will want to recover the merit for each of the projections.  ********
     *
     * I removed `M`.  No longer relevant.
     */
    double computeMerit(IntLatticeExt<Int, Real> &lat,
            IntLattice<Int, Real> *proj);

    /*
     * This function calculates the Figure of Merit for all projections
     * consisting of successive coordinates of the forms
     * {1, 2, ..., m_t.size()} to {1, 2, ..., m_t[0]}
     * for the primal lattice.
     * The value 0 is returned if an error occurs while calculating the shortest vector.
     *
     * I removed "Primal".   Was no longer relevant.
     */
    double computeMeritSucc(IntLatticeExt<Int, Real> &lat,
            IntLattice<Int, Real> *proj);

    /*
     * This functions calculates the figure of merit of the primal lattice
     * for all projections consisting of non-successive coordinates.
     * The value 0 is returned if an error occurs while calculating the
     * shortest vector.
     */
    double computeMeritNonSucc(IntLatticeExt<Int, Real> &lat,
            IntLattice<Int, Real> *proj);

    /*
     * Internal `CoordinateSets` object used to store the sets of projections
     * that are considered. It is populated by the `setTVector` function.
     * This object will contain a set of projections of each order.
     */
    CoordinateSets::FromRanges m_coordRange;

    /*
     * Internal normalizer object for storing normalizing values in FiguresOfMeritM class
     */
    Normalizer *m_norma;

    /*
     * Internal reducer object used for finding the shortest vector of a projection
     */
    Reducer<Int, Real> *m_red;

protected:

    /*
     * The vector 'm_t' defines the set of projections for which the FOM is
     * computed . It contains the values of t_1,..., t_d in the definition of the FOM.
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
     * Indicates if the first coordinate will always be included in all the projections
     * over the non-successive coordinates.  If true, we use @f$S^{(1)}(s,t_s)$@f  in
     * Section 10 of the user's guide, otherwise we use @f$S(s,t_s)$@f and we therefore
     * have a larger set of projections.
     */
    bool m_firstCoordinateAlwaysIn = true;

    /*
     * As soon as we know the FOM is above this bound, its calculation stopped.
     */
    double m_highbound = 1.0;

    /*
     * As soon as we know the FOM is below this bound, its calculation stopped.
     */
    double m_lowbound = 0.0;

    /*
     * The vector 'm_b' is used to store the length of the basis vectors after pre-reduction has
     * been performed.    *****  What is this?   The basis is in an array of double  ?????  *****
     */
    double *m_b;

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
    m_firstCoordinateAlwaysIn = includeFirst;  // This must be set first.    ********
    setTVector(t);
    setReductionMethod(meth);
    setNormalizer(norma);
    m_red = &red;
}

//=========================================================================

template<typename Int>
void FigureOfMeritM<Int>::setTVector(const vector<int64_t> &t) {
    m_t = t;
    // Clear the CoordinateSets object.
    for (Coordinates::size_type i = 0; i < m_t.size(); i++)
        m_coordRange.excludeOrder(i);
    /* Defines the lower bound for the range of coordinates that belong to a projection.
     * It is 2 if the first coordinate belongs to all the projections, because we do
     * not have to consider coordinate 1. Otherwise it is 1.
     * In the first case, coordinate 1 is added to all the projections as an extra coord.
     * See the doc of the class `FromRanges` in `CoordinateSets`.
     */
    int64_t min_dim = 1;
    if (m_firstCoordinateAlwaysIn) min_dim = 2;
    int64_t d = static_cast<int>(t.size()); // Number of orders for the projections.
    // ******   Are the values in t[0],..., t[d-1]  or in t[1],...,t[d] ??????   *******
    for (int64_t i = 1; i < d; i++)
        // Adds the set of projections of order i.
        // I think we want orders from 2 to d.  So the first i below should be i+1.  *******
        m_coordRange.includeOrder(i, min_dim, t[i], true);
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

    m_b = new double[m_t.size()];
    merit = computeMeritNonSucc(lat, proj);
    if (merit < minmerit)
        minmerit = merit;
    // In any of these cases the calculation is stopped
    if (minmerit == 0 || minmerit < m_lowbound || minmerit > m_highbound)
        return 0;

    m_b = new double[lat.getDim()];
    merit = computeMeritSucc(lat, proj);
    if (merit < minmerit)
        minmerit = merit;
    // In any of these cases the calculation is stopped
    if (minmerit == 0 || minmerit < m_lowbound || minmerit > m_highbound)
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
            //  ******   Why not use LLL_FPZZflex for this one as well?   *********
            this->m_red->redBKZ(lat.getBasis(), this->m_delta,
                    this->m_blocksize);
        } else if (m_doingLLL) {
            LLL_FPZZflex(lat.getBasis(), this->m_delta, j, j, this->m_b);
        } else if (this->m_reductionMethod == PAIRBB) {
            this->m_red->redDieter(0);
        }
        if (!m_doingBB) {
            lat.updateSingleVecNorm(0, j);
            if (lat.getNormType() == L2NORM) {
                NTL::conv(merit,
                        sqrt(lat.getVecNorm(0)) / this->m_norma->getBound(j));
            } else {
                NTL::conv(merit,
                        lat.getVecNorm(0) / this->m_norma->getBound(j));
            }
        } else {
            if (!m_red->shortestVector(lat))
                return 0;
            NTL::conv(merit,
                    m_red->getMinLength() / this->m_norma->getBound(j));
        }
        if (merit < minmerit)
            minmerit = merit;
        if (minmerit <= this->m_lowbound || minmerit > this->m_highbound)
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

    for (auto it = m_coordRange.begin(); it != m_coordRange.end(); it++) {
        coord = *it;
        lat.buildProjection(proj, *it, this->m_delta);  // Must have withDual = false   ***
        if (m_doingBKZ) {
            //  ******   Why not use LLL_FPZZflex for this one as well?   *********
            this->m_red->redBKZ(proj->getBasis(), this->m_delta,
                    this->m_blocksize);
        } else if (m_doingLLL) {
            LLL_FPZZflex(proj->getBasis(), this->m_delta, coord.size(),
                    coord.size(), this->m_b);
        } else if (this->m_reductionMethod == PAIRBB) {
            this->m_red->redDieter(0);
        }
        if (!m_doingBB) {
            proj->updateSingleVecNorm(0, coord.size());
            if (lat.getNormType() == L2NORM) {
                NTL::conv(merit,
                        sqrt(proj->getVecNorm(0))
                                / this->m_norma->getBound(proj->getDim()));
            } else {
                NTL::conv(merit,
                        proj->getVecNorm(0)
                                / this->m_norma->getBound(proj->getDim()));
            }
        } else {
            if (!m_red->shortestVector(*proj))
                return 0;
            merit = NTL::conv<double>(
                    m_red->getMinLength()
                            / this->m_norma->getBound(proj->getDim()));
        }
        if (merit < minmerit)
            minmerit = merit;
        if (minmerit <= this->m_lowbound || minmerit > this->m_highbound)
            return 0;
    }
    return minmerit;
}

template class FigureOfMeritM<NTL::ZZ> ;

} // end namespace LatticeTester

#endif
