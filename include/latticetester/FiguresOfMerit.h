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

#ifndef LATTICETESTER_FIGURESOFMERIT_H
#define LATTICETESTER_FIGURESOFMERIT_H

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
#include "latticetester/Rank1Lattice.h"
#include "latticetester/Reducer.h"


#include "latticetester/NormaBestLat.h"
#include "latticetester/NormaBestBound.h"
#include "latticetester/WeightsUniform.h"
#include "latticetester/CoordinateSets.h"

namespace LatticeTester {

/**
 * This class offers methods (functions) to calculate figure of merit for a given 
 * IntLattice object
 *
 */

template<typename Int> class FiguresOfMerit {

private:
    typedef NTL::vector<Int> IntVec;
    typedef NTL::matrix<Int> IntMat;
public:	
        
    /*
     * Constructor of a FiguresOfMerit object. Needs a reducer object
     * to be passed. The vector 't' defines the set of 
     * dimensions for which the figure of merit is calculated. It contains
     * the values of t_1,..., t_d in the definition of the FoM.
     */
    FiguresOfMerit (Reducer<Int, Real> & red, const vector<int64_t> & t);
    
    /*
     * Sets a new vector 't' =  t_1,..., t_d in the definition of the FoM.
     */
    void setTVector(const vector<int64_t> & t);
    
    /*
     * Directly sets a new coordinate set.
     */
    void setCoordinateSets (CoordinateSets::FromRanges coordRange) { m_coordRange = coordRange; }
    
    /*
     * Method to set the weights used during the calculation of the 
     * figure of merit. 
     */
    void setWeights (Weights & w) { m_weights = &w; }
	
    /*
     * This function calculates the Figure of Merit M of a given lattice 'lat'
     * and should be called by the user. The vector 't' defines the set of 
     * dimensions for which the figure of merit is calculated. It contains
     * the values of t_1,..., t_d in the definition of the FoM.
     * The IntLattice object 'proj' is needed for saving projections.
     *  The value 0 is returned if an error occurs
     * while calculating the shortest vector.
     */
    double computeMeritM (IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj);
    
    /*
     * This function calculates the Figure of Merit for a given lattice 'lat' 
     * without applying any projection. 
     */
    double computeMeritMNoProj (IntLatticeExt<Int, Real> & lat);
    
    /*
     * This functions calculates the Figure of Merit for all projections 
     * consisting of successive coordinates of the forms 
     * {1, 2, ..., m_t.size()} to {1, 2, ..., m_t[0]}
     * The value 0 is returned if an error occurs while calculating the shortest vector.
     */
    double computeMeritMSucc (IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj);
	
    /*
     * This functions calculates the Figure of Merit for all projections 
     * consisting of non-successive coordinates. The value 0 is returned if 
     * an error occurs while calculating the shortest vector.
     */
//    double computeMeritMNonSucc (IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj);
    
    /*
     * ToDo: For the primal lattice
     */
    double computeMeritMNonSuccPrimal (IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj);

    /*
     * ToDo: For the dual lattice
     */
    double computeMeritMNonSuccDual (IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj);


    /*
     * A function which allows to set set the normalizer to user-specified value.
    */
    void setNormalizer (Normalizer & norm) { m_norma  = &norm; }

    /* 
     * Defines the type of the figure of merit
     */	
    MeritType m_fom = MERITM; 
    
    /*
     * If FoM is above this bound, then calculation of FoM is stopped
     */	
    double m_highbound = 1; 
    
    /*
     * If FoM is below this bound, then calculation of FoM is stopped.
     */	
    double m_lowbound = 0; 
    
    /*
     * If true, the FOM is calculated for the dual lattice.
     */
    bool m_fomInDual = false; 
    
    /* 
     * The reduction type used to compute (or estimate) the FOM.
     */
     ReductionType m_reductionMethod = BKZBB; 
    
    /*
     * The delta parameter for LLL or BKZ reduction
     */
    double m_delta = 0.99999; 
    
    /*
     * Blocksize of BKZ algorithm
     */
    int64_t m_blocksize = 10; 
    
    /* An integer giving the maximum number of nodes to be examined in any given 
     * BB procedure. When that value is exceeded, the BB is stopped and the generator 
     * is discarded. The number of discarded generators is given in the results. A small 
     * value of this maxnodesBB parameter can make the program run faster (sometimes much 
     * faster), permitting to examine more generators, but will increase the chances of 
     * rejecting good generators. The default value is 10^8 .
     * 
     */
    int64_t m_maxNodesBB = 10000000;
    
    /*
     * If set to true successive coordinates a considered first when calculating FigureOfMerit
     */
    bool m_succCoordFirst = false;
    
    /*
     * Sets if first variable shall always be included for the non-successive coordinates
     */
    bool m_projectStationary = true;
    
    /*
     * incDualOnly decides if only the dual basis shall be updated 
     * when increasing the dimension (but not the primal basis) 
     */
    bool m_incDualOnly = false;
    	
    /*
     * Type of projection construction
     */	
    ProjConstructType m_pctype = UPPERTRIPROJ;
    
    /*
     * Variable containing the Weights for the FoM 
     */
    Weights *m_weights; 
    
    /* The vector 'm_t' defines the set of 
     * dimensions for which the figure of merit is calculated. It contains
     * the values of t_1,..., t_d in the definition of the FOM M.
     */ 
    vector<int64_t> m_t; 
    
    /*
     * Stores the matrix of the projection basis
     */
    IntMat m_projBasis; 
    
    /*
     * Matrix necessary for intermediate steps of calculation
     */
    IntMat m_temp;
	    
    /*
     * CoordinateSets object is used to store sets of coordinates
     */
    CoordinateSets::FromRanges m_coordRange;  
    
    /*
     * Normalizer for storing normalizing values in FiguresOfMerit class
     */
    Normalizer *m_norma;
    
    /*
     * Reducer object used for finding the shortest vector of a projection
    */
    Reducer<Int, Real> *m_red;     
    
    /*
     * test
     */
    vector<int64_t> ttest;
};

//============================================================================
// Implementation

template<typename Int>
FiguresOfMerit<Int>::FiguresOfMerit (Reducer<Int, Real> & red, const vector<int64_t> & t) {    
   setTVector(t);
   m_red = &red;
}

template<typename Int>
void FiguresOfMerit<Int>::setTVector (const vector<int64_t> & t) { 
   // Clear CoordinateSets object
   for (Coordinates::size_type i = 0; i < m_t.size(); i++)
      m_coordRange.excludeOrder(i);
   /* Defines the lower bound for the range of dimensions. 
   * Is equal to 2 for stationary projection because the first coordinate is always included then.
   * For non-stationary projections it is equal to 1.
   */   
   int64_t min_dim; 
   if (m_projectStationary) {
      min_dim = 2;
   } else min_dim = 1;
   int64_t lower_dim = static_cast<int>(t.size());
   for (int64_t i = 1; i < lower_dim; i++)
       m_coordRange.includeOrder(i,min_dim,t[i],true);
   m_t = t;
}
  

template<typename Int>
double FiguresOfMerit<Int>::computeMeritM (IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj) {
   double merit = 0;
   double minmerit = 1.0;
   
   if (m_succCoordFirst == true) { 
       minmerit = computeMeritMSucc(lat, proj);
       // In any of these cases the calcuation is stopped
       if (minmerit == 0 || minmerit < m_lowbound || minmerit > m_highbound) return 0;
   }   
   
   if (!m_fomInDual) merit = computeMeritMNonSuccPrimal(lat, proj);
   else  merit = computeMeritMNonSuccDual(lat, proj);
   if (merit < minmerit) minmerit = merit;
   // In any of these cases the calcuation is stopped
   if (minmerit == 0 || minmerit < m_lowbound || minmerit > m_highbound) return 0;

   if (m_succCoordFirst == false) {
       merit = computeMeritMSucc(lat, proj);
       if (merit < minmerit) minmerit = merit;
       // In any of these cases the calcuation is stopped
       if (minmerit == 0 || minmerit < m_lowbound || minmerit > m_highbound) return 0;
   } 
   return minmerit; 
}

//=========================================================================

template<typename Int>
double FiguresOfMerit<Int>::computeMeritMNoProj (IntLatticeExt<Int, Real> & lat) {
   double shortest, merit;
   merit = 0.0;
   if (m_fomInDual) lat.dualize();
   lat.updateVecNorm();
   lat.sort(0); 
   if (m_reductionMethod == BKZBB || m_reductionMethod == BKZ) {
      m_red->redBKZ(lat.getBasis(), m_delta, m_blocksize);  
   } else if (m_reductionMethod == LLLBB || m_reductionMethod == LLL) {
      m_red->redLLLNTL(lat.getBasis(), m_delta);  
   } else if (m_reductionMethod == PAIRBB) {
      m_red->redDieter(0);
   }
   if (m_reductionMethod == BKZBB || m_reductionMethod == LLLBB || m_reductionMethod == PAIRBB) {
          if (!m_red->shortestVector(lat)) return 0;
          shortest = NTL::conv<double>(m_red->getMinLength());
   } else
   {
       shortest = lat.getShortestLengthBasis();
   }
   merit = shortest / m_norma->getBound(lat.getDim());
   if (m_fomInDual) lat.dualize();
   return merit;
}

//=========================================================================

template<typename Int>
double FiguresOfMerit<Int>::computeMeritMSucc (IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj) {
   double merit = 0;
   double minmerit = 1.0;
   double dim = lat.getBasis().NumCols();
   int64_t low = static_cast<int64_t>(m_t.size());
   
   lat.buildBasis(low+1);
   merit = computeMeritMNoProj(lat);
   BasisConstruction<Int>::mDualBasis(lat.getBasis(), lat.getDualBasis(), lat.getModulo());
   if (merit == 0) return merit;
   for (int j = low +1; j < m_t[0]; j++)
   {
       if (m_incDualOnly == true) {
           lat.incDimDualBasis();
       } else {
           lat.incDimBasis();
       }
       merit = computeMeritMNoProj(lat);
       BasisConstruction<Int>::mDualBasis(lat.getBasis(), lat.getDualBasis(), lat.getModulo());
       if (merit < minmerit) minmerit = merit;
       if (merit == 0 || merit < m_lowbound || merit > m_highbound) return 0;
   }   
   if (m_t[0] < dim)
      lat.buildBasis(dim);
   return minmerit;
}


//=========================================================================
template<typename Int>
double FiguresOfMerit<Int>::computeMeritMNonSuccPrimal (IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj) {
    double merit = 0;
    double minmerit = 1.0;
    double shortest = 0.0;
    for (auto it = m_coordRange.begin(); it != m_coordRange.end(); it++){
    	//calculat LLL projection
       BasisConstruction<Int>::projectionConstructionLLL(lat.getBasis(), m_projBasis, *it, lat.getModulo(), m_delta); 
       // Define IntLattice based on projected basis
       proj.setBasis(m_projBasis, lat.getModulo(), m_projBasis.NumCols());
       // Calculate shortest vector
       if (m_reductionMethod == BKZBB || m_reductionMethod == LLLBB || m_reductionMethod == PAIRBB) {
           if (!m_red->shortestVector(proj)) return 0;
           shortest = NTL::conv<double>(m_red->getMinLength());
       } else {
           shortest = proj.getShortestLengthBasis();
       }
       merit = shortest / m_norma->getBound(proj.getDim());
       if (merit < minmerit) minmerit = merit;
       if (merit == 0 || merit < m_lowbound || merit > m_highbound) return 0;
    }
    return minmerit;  
}

//=========================================================================
template<typename Int>
double FiguresOfMerit<Int>::computeMeritMNonSuccDual (IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj) {
    double merit = 0;
    double minmerit = 1.0;
    double shortest = 0.0;
    for (auto it = m_coordRange.begin(); it != m_coordRange.end(); it++){
       // Calculate upper-triangular basis for projection
       BasisConstruction<Int>::projectionConstructionUpperTri(lat.getBasis(), m_projBasis, *it, lat.getModulo());
       // Calculate upper-triangular mdual basis
       BasisConstruction<Int>::mDualUpperTriangular(m_projBasis, m_temp, lat.getModulo());
       // Define IntLattice based on mdual basis
       proj.setBasis(m_temp, lat.getModulo(), m_projBasis.NumCols());
       // Apply selected reduction technique
       if (m_reductionMethod == BKZBB || m_reductionMethod == BKZ) {
           m_red->redBKZ(proj.getBasis(), m_delta, m_blocksize);  
       } else if (m_reductionMethod == LLLBB || m_reductionMethod == LLL) {
           m_red->redLLLNTL(proj.getBasis(), m_delta);  
       } else if (m_reductionMethod == PAIRBB) {
           m_red->redDieter(0);
       }
       // Calculate shortest vector
       if (m_reductionMethod == BKZBB || m_reductionMethod == LLLBB || m_reductionMethod == PAIRBB) {
           if (!m_red->shortestVector(proj)) return 0;
           shortest = NTL::conv<double>(m_red->getMinLength());
       } else {
           shortest = proj.getShortestLengthBasis();
       }
       merit = shortest / m_norma->getBound(proj.getDim());
       if (merit < minmerit) minmerit = merit;
       if (merit == 0 || merit < m_lowbound || merit > m_highbound) return 0;
    }
    return minmerit;  
}

//=========================================================================
	
template class FiguresOfMerit<NTL::ZZ>;
//template class FiguresOfMerit<std::int64_t>;

} // end namespace LatticeTester

#endif
