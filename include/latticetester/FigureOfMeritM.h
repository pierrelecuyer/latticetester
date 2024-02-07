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
 * Essential ingredients to pass to the constructor:
 * (1) the vector (t_1,...,t_d);
 * (2) the reduction method to be used;
 * (3) a Reducer object;
 * (4) a Normalizer object.
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
     */
    FigureOfMeritM(const vector<int64_t> &t, ReductionType &meth,
            Reducer<Int, Real> &red, Normalizer & norma);

    /*
     * Sets a new vector 't' =  t_1,..., t_d in the definition of the FoM.
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
     */
    void setCoordinateSets(CoordinateSets::FromRanges coordRange) {
        m_coordRange = coordRange;
    }

    /*
     * Method to set the weights used during the calculation of the
     * figure of merit.
     */
    //void setWeights (Weights & w) { m_weights = &w; }
    /*
     * Method to set the low and the high bound used during the calculation of the
     * figure of merit.
     */
    void setBounds(double &low, double &high) {
        m_highbound = high;
        m_lowbound = low;
    }

    /*
     * A function which allows to set set the normalizer to `norma`.
     */
    void setNormalizer(Normalizer &norma) {
        m_norma = &norma;
    }
    
    /* 
     * Sets if first variable shall always be included for the non-successive coordinates
     */
    void setFirstCoordinateAlwaysIn (bool & first) {
        m_firstCoordinateAlwaysIn = first;
    }

    /*
     * This function calculates the Figure of Merit M of a given lattice 'lat'
     * and should be called by the user. The vector 't' defines the set of
     * dimensions for which the figure of merit is calculated. It contains
     * the values of t_1,..., t_d in the definition of the FoM.
     * The IntLattice object 'proj' is needed for saving projections.
     * The value 0 is returned if an error occurs
     * while calculating the shortest vector.
     */
    double computeMeritM(IntLatticeExt<Int, Real> &lat,
            IntLattice<Int, Real> *proj);

    /*
     * This function calculates the Figure of Merit for all projections
     * consisting of successive coordinates of the forms
     * {1, 2, ..., m_t.size()} to {1, 2, ..., m_t[0]}
     * for the primal lattice.
     * The value 0 is returned if an error occurs while calculating the shortest vector.
     */
    double computeMeritMSuccPrimal(IntLatticeExt<Int, Real> &lat,
            IntLattice<Int, Real> *proj);

    /*
     * This functions calculates the figure of merit of the primal lattice
     * for all projections consisting of non-successive coordinates.
     * The value 0 is returned if an error occurs while calculating the
     * shortest vector.
     */
    double computeMeritMNonSuccPrimal(IntLatticeExt<Int, Real> &lat,
            IntLattice<Int, Real> *proj);

    /*
     * If set to true successive coordinates a considered first when calculating FigureOfMerit
     */
    //bool m_succCoordFirst = false;

    /*
     * Internal CoordinateSets object is used to store sets of coordinates.
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
     * The vector 'm_t' defines the set of
     * dimensions for which the figure of merit is calculated. It contains
     * the values of t_1,..., t_d in the definition of the FOM M.
     */
    vector<int64_t> m_t;
    
    /*
     * The reduction type used to compute (or estimate) the FOM.
     */
    ReductionType m_reductionMethod = BKZBB;
    
    /*
     * Variable containing the Weights for the FoM
     */
    //Weights *m_weights;
    
    /*
     * The boolean variable indicates whether the BB algorithm needs to be performed or not.
     */
    bool m_doingBB;
    
    /*
     * The delta parameter for LLL and BKZ reduction
     */
    double m_delta = 0.99999;

    /*
     * Blocksize of BKZ algorithm
     */
    int64_t m_blocksize = 10;
    
    /*
     * Indicates if first variable shall always be included for the non-successive coordinates
     */
    bool m_firstCoordinateAlwaysIn = true;  
    
    /*
     * If FoM is above this bound, then calculation of FoM is stopped
     */
    double m_highbound = 1;

    /*
     * If FoM is below this bound, then calculation of FoM is stopped.
     */
    double m_lowbound = 0;
    
     /*
     * The vector 'm_b' is used to store the length of the basis after pre-reduction has
     * been performed.    *****  What is this?   The basis is in an array of double  ?????  *****
     */
    double *m_b;
    
};

//============================================================================
// Implementation

template<typename Int>
FigureOfMeritM<Int>::FigureOfMeritM(const vector<int64_t> &t,
        ReductionType &meth, Reducer<Int, Real> &red, Normalizer & norma) {
    setTVector(t);
    setReductionMethod(meth);
    setNormalizer(norma);
    m_red = &red;
}

//=========================================================================

template<typename Int>
void FigureOfMeritM<Int>::setTVector(const vector<int64_t> &t) {
    // Clear CoordinateSets object
    for (Coordinates::size_type i = 0; i < m_t.size(); i++)
        m_coordRange.excludeOrder(i);
    /* Defines the lower bound for the range of dimensions.
     * Is equal to 2 for stationary projection because the first coordinate is always included then.
     * For non-stationary projections it is equal to 1.
     */
    int64_t min_dim;
    if (m_firstCoordinateAlwaysIn) {
        min_dim = 2;
    } else
        min_dim = 1;
    int64_t lower_dim = static_cast<int>(t.size());
    for (int64_t i = 1; i < lower_dim; i++)
        m_coordRange.includeOrder(i, min_dim, t[i], true);
    m_t = t;
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
    m_delta = delta;
    m_blocksize = blocksize;
}

//=========================================================================

template<typename Int>
double FigureOfMeritM<Int>::computeMeritM(IntLatticeExt<Int, Real> &lat,
        IntLattice<Int, Real> *proj) {
    double merit = 0;
    double minmerit = 1.0;

    m_b = new double[m_t.size()];
    merit = computeMeritMNonSuccPrimal(lat, proj);
    if (merit < minmerit)
        minmerit = merit;
// In any of these cases the calcuation is stopped
    if (minmerit == 0 || minmerit < m_lowbound || minmerit > m_highbound)
        return 0;

    m_b = new double[lat.getDim()];
    merit = computeMeritMSuccPrimal(lat, proj);
    if (merit < minmerit)
        minmerit = merit;
// In any of these cases the calcuation is stopped
    if (minmerit == 0 || minmerit < m_lowbound || minmerit > m_highbound)
        return 0;

    return minmerit;
}

//=========================================================================
template<typename Int>
double FigureOfMeritM<Int>::computeMeritMSuccPrimal(
        IntLatticeExt<Int, Real> &lat, IntLattice<Int, Real> *proj) {
    double merit = 0;
    double minmerit = 1.0;
    int64_t lower_dim = static_cast<int64_t>(this->m_t.size());
    lat.buildBasis(lower_dim);
    for (int64_t j = lower_dim + 1; j < this->m_t[0] + 1; j++) {
        lat.incDimBasis();
        if (this->m_reductionMethod == BKZBB
                || this->m_reductionMethod == BKZ) {
            this->m_red->redBKZ(lat.getBasis(), this->m_delta,
                    this->m_blocksize);
        } else if (this->m_reductionMethod == LLLBB
                || this->m_reductionMethod == LLL) {
            LLL_FPZZflex(lat.getBasis(), this->m_delta, j, j, this->m_b);
        } else if (this->m_reductionMethod == PAIRBB) {
            this->m_red->redDieter(0);
        }
        if (!m_doingBB) {
            lat.updateSingleVecNorm(0, j);
            NTL::conv(merit,
                    sqrt(lat.getVecNorm(0)) / this->m_norma->getBound(j));
        } else {
            if (!m_red->shortestVector(lat))
                return 0;
            merit = NTL::conv<double>(
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
double FigureOfMeritM<Int>::computeMeritMNonSuccPrimal(
        IntLatticeExt<Int, Real> &lat, IntLattice<Int, Real> *proj) {
    double merit = 0;
    double minmerit = 1.0;
    Coordinates coord;

    for (auto it = m_coordRange.begin(); it != m_coordRange.end(); it++) {
        coord = *it;
        lat.buildProjection(proj, *it, this->m_delta);
        if (this->m_reductionMethod == BKZBB
                || this->m_reductionMethod == BKZ) {
            this->m_red->redBKZ(proj->getBasis(), this->m_delta,
                    this->m_blocksize);
        } else if (this->m_reductionMethod == LLLBB
                || this->m_reductionMethod == LLL) {
            LLL_FPZZflex(proj->getBasis(), this->m_delta, coord.size(),
                    coord.size(), this->m_b);
        } else if (this->m_reductionMethod == PAIRBB) {
            this->m_red->redDieter(0);
        }
        if (!m_doingBB) {
            proj->updateSingleVecNorm(0, coord.size());
            NTL::conv(merit,
                    sqrt(proj->getVecNorm(0))
                            / this->m_norma->getBound(proj->getDim()));
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
