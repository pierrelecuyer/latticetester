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
 * This class offers methods (functions) to calculate figure of merit for a given 
 * IntLattice object
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
    FigureOfMeritM (const vector<int64_t> & t, ReductionType & meth, Reducer<Int, Real> & red);
    
    /*
     * Sets a new vector 't' =  t_1,..., t_d in the definition of the FoM.
     */
    void setTVector(const vector<int64_t> & t);
    
    /*
     * Sets the reduction method and the parameter delta for LLL / BKZ algorithm
     * as well as the blocksize used in the BKZ algorithm.
     */
    void setReductionMethod(ReductionType & meth, double delta = 0.99999, int64_t blocksize = 10);
    
    /*
     * Directly sets a new coordinate set.
     */
    void setCoordinateSets (CoordinateSets::FromRanges coordRange) { m_coordRange = coordRange; }
    
    /*
     * Method to set the weights used during the calculation of the 
     * figure of merit. 
     */
    //void setWeights (Weights & w) { m_weights = &w; }
    
    /*
     * Method to set the low and the high bound used during the calculation of the 
     * figure of merit. 
     */
	void setBounds(double & low, double & high) {m_highbound = high; m_lowbound = low;}
    
    /*
     * This function calculates the Figure of Merit M of a given lattice 'lat'
     * and should be called by the user. The vector 't' defines the set of 
     * dimensions for which the figure of merit is calculated. It contains
     * the values of t_1,..., t_d in the definition of the FoM.
     * The IntLattice object 'proj' is needed for saving projections.
     * The value 0 is returned if an error occurs
     * while calculating the shortest vector.
     */
    double computeMeritM (IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj);
    
    /*
     * WAS NECESSARY FOR TESTING PURPOSES (FoMCalc.h) ONLY AND MAY BE DELETED AFTER TESTING
     * This function calculates the Figure of Merit for a given lattice 'lat' 
     * without applying any projection. 
     */
    double computeMeritMNoProj (IntLatticeExt<Int, Real> & lat, bool & m_fomInDual);
    
    /*
     * This function calculates the Figure of Merit for all projections 
     * consisting of successive coordinates of the forms 
     * {1, 2, ..., m_t.size()} to {1, 2, ..., m_t[0]}
     * for the primal lattice.
     * The value 0 is returned if an error occurs while calculating the shortest vector.
     */
    double computeMeritMSuccPrimal (IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj);
    
    /*
     * This functions calculates the figure of merit of the primal lattice 
     * for all projections consisting of non-successive coordinates. 
     * The value 0 is returned if an error occurs while calculating the 
     * shortest vector.
     */
    double computeMeritMNonSuccPrimal (IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj);

    /*
     * A function which allows to set set the normalizer to user-specified value.
    */
    void setNormalizer (Normalizer & norm) { m_norma  = &norm; }
        
    /* 
     * The reduction type used to compute (or estimate) the FOM.
     */
     ReductionType m_reductionMethod = BKZBB; 
    
    /* 
     * The delta parameter for LLL and BKZ reduction
     */
    double m_delta = 0.99999; 
    
    /*
     * Blocksize of BKZ algorithm
     */
    int64_t m_blocksize = 10; 
    
    /* 
     * If FoM is above this bound, then calculation of FoM is stopped
     */	
    double m_highbound = 1; 
    
    /* 
     * If FoM is below this bound, then calculation of FoM is stopped.
     */	
    double m_lowbound = 0; 
    
    /*
     * The boolean variable indicates whether the BB algorithm needs to be performed or not.
     */
    bool m_doingBB;
    
    /* 
     * If set to true successive coordinates a considered first when calculating FigureOfMerit
     */
    //bool m_succCoordFirst = false;
    
    /* 
     * Sets if first variable shall always be included for the non-successive coordinates
     */
    bool m_firstCoordinateAlwaysIn= true;
    	
    /*
     * Type of projection construction
     */	
    ProjConstructType m_pctype = UPPERTRIPROJ;

    /* 
     * Variable containing the Weights for the FoM 
     */
    //Weights *m_weights; 
    
    /* 
     * The vector 'm_t' defines the set of 
     * dimensions for which the figure of merit is calculated. It contains
     * the values of t_1,..., t_d in the definition of the FOM M.
     */ 
    vector<int64_t> m_t; 
    
    /*
     * The vector 'm_b' is used to store the basis after pre-reduction has
     * been performed.           *****   The basis is in an array of double  ?????
    */
    double *m_b;
    
    /* 
     * Internal variable which stores the matrix of the projection basis.
     */
    IntMat m_projBasis; 
    
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
    
};

//============================================================================
// Implementation

template<typename Int>
FigureOfMeritM<Int>::FigureOfMeritM (const vector<int64_t> & t, ReductionType & meth, Reducer<Int, Real> & red) {    
   setTVector(t);
   setReductionMethod(meth);
   m_red = &red;
}

//=========================================================================

template<typename Int>
void FigureOfMeritM<Int>::setTVector (const vector<int64_t> & t) { 
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
   } else min_dim = 1;
   int64_t lower_dim = static_cast<int>(t.size());
   for (int64_t i = 1; i < lower_dim; i++)
       m_coordRange.includeOrder(i,min_dim,t[i],true);
   m_t = t;
}

//=========================================================================

template<typename Int>
void FigureOfMeritM<Int>::setReductionMethod(ReductionType & meth, double delta, int64_t blocksize) {
	m_reductionMethod = meth;
	if (m_reductionMethod == BKZBB || m_reductionMethod == LLLBB || m_reductionMethod == PAIRBB) {
		m_doingBB = true;
	}
	else {
		m_doingBB = false;
	}
	m_delta = delta;
	m_blocksize = blocksize;
}

//=========================================================================

template<typename Int>
double FigureOfMeritM<Int>::computeMeritM (IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj) {
   double merit = 0;
   double minmerit = 1.0;
   
//   if (m_succCoordFirst == true) {    
//	   m_b = new double[lat.getDim()];
//     minmerit = computeMeritMSuccPrimal(lat, proj);
//     // In any of these cases the calcuation is stopped
//     if (minmerit == 0 || minmerit < m_lowbound || minmerit > m_highbound) return 0;
//   }   

   m_b = new double[m_t.size()];
   merit = computeMeritMNonSuccPrimal(lat, proj);
   if (merit < minmerit) minmerit = merit;
   // In any of these cases the calcuation is stopped
   if (minmerit == 0 || minmerit < m_lowbound || minmerit > m_highbound) return 0;
   
// if (m_succCoordFirst == false) {    
	   m_b = new double[lat.getDim()];
       merit = computeMeritMSuccPrimal(lat, proj);
       if (merit < minmerit) minmerit = merit;
       // In any of these cases the calcuation is stopped
       if (minmerit == 0 || minmerit < m_lowbound || minmerit > m_highbound) return 0;
//   }    
   
   return minmerit; 
}

// Was necessary only for testing (FoMCalc.h) - may be deleted later
//=========================================================================
template<typename Int>
double FigureOfMeritM<Int>::computeMeritMNoProj (IntLatticeExt<Int, Real> & lat, bool & m_fomInDual) {
   double shortest, merit;
   merit = 0.0;
   // Switches primal and dual lattice if calculations shall be done for the dual
   if (m_fomInDual) lat.dualize();
   lat.updateVecNorm();
   lat.sortBasis(0); 
   // Apply pre-reduction
   if (m_reductionMethod == BKZBB || m_reductionMethod == BKZ) {
      m_red->redBKZ(lat.getBasis(), m_delta, m_blocksize);  
   } else if (m_reductionMethod == LLLBB || m_reductionMethod == LLL) {
      m_red->redLLLNTL(lat.getBasis(), m_delta);  
   } else if (m_reductionMethod == PAIRBB) {
      m_red->redDieter(0);
   }
   // Calculate the shortest vector or get it
   if (m_doingBB) {
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
// CODE IS UPDATED BUT HAS NOT BEEN TESTED (in particual the use of buildDualBasisFullMatrix and incDimBasisFullMatrix)
template<typename Int>
double FigureOfMeritM<Int>::computeMeritMSuccPrimal (IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj) {
   double merit = 0;
   double minmerit = 1.0;
   int64_t lower_dim = static_cast<int64_t>(this->m_t.size());
   lat.buildBasisFullMatrix(lower_dim+1);
   if (this->m_reductionMethod == BKZBB || this->m_reductionMethod == BKZ) {
      this->m_red->redBKZ(lat.getBasis(), this->m_delta, this->m_blocksize);  
   } else if (this->m_reductionMethod == LLLBB || this->m_reductionMethod == LLL) {
     LLL_FPZZflex(lat.getBasis(), this->m_delta, lower_dim+1, lower_dim+1, this->m_b);
   } else if (this->m_reductionMethod == PAIRBB) {
     this->m_red->redDieter(0);
   }
   if (!m_doingBB) {
	      NTL::conv(merit, sqrt(this->m_b[0]) / this->m_norma->getBound(lower_dim+1));
     } else {
       m_projBasis.SetDims(lower_dim+1,lower_dim+1); //m_projBasis resized -> Can this be avoided?
       for (int k = 0; k < lower_dim+1; k++) {
          for (int l = 0; l < lower_dim+1; l++) m_projBasis[k][l] = lat.getBasis()[k][l];
       }
       proj.setBasis(m_projBasis, m_projBasis.NumCols());
       if (!m_red->shortestVector(proj)) return 0;
       merit = NTL::conv<double>(m_red->getMinLength() / this->m_norma->getBound(lower_dim+1));
     }
   if (merit == 0) return merit;
   minmerit = merit;
   for (int64_t j = lower_dim+2; j < this->m_t[0] + 1; j++)
   {
	   lat.incDimBasisFullMatrix(j);  
	   if (this->m_reductionMethod == BKZBB || this->m_reductionMethod == BKZ) {
	       this->m_red->redBKZ(lat.getBasis(), this->m_delta, this->m_blocksize);
	    } else if (this->m_reductionMethod == LLLBB || this->m_reductionMethod == LLL) {
	    	LLL_FPZZflex(lat.getBasis(), this->m_delta, j, j, this->m_b);
	    } else if (this->m_reductionMethod == PAIRBB) {
	       this->m_red->redDieter(0);
	    }
     if (!m_doingBB) {
	      NTL::conv(merit, sqrt(this->m_b[0]) / this->m_norma->getBound(j));
     } else {
       m_projBasis.SetDims(j,j); //m_projBasis resized -> Can this be avoided?
       for (int k = 0; k < j; k++) {
          for (int l = 0; l < j; l++) m_projBasis[k][l] = lat.getBasis()[k][l];
       }
       proj.setBasis(m_projBasis, m_projBasis.NumCols());
       if (!m_red->shortestVector(proj)) return 0;
       merit = NTL::conv<double>(m_red->getMinLength() / this->m_norma->getBound(j));
     }
       if (merit < minmerit) minmerit = merit;
       if (minmerit <= this->m_lowbound || minmerit > this->m_highbound) return 0;
   }  
   return minmerit; 
}


//=========================================================================
template<typename Int>
double FigureOfMeritM<Int>::computeMeritMNonSuccPrimal (IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj) {
    double merit = 0;
    double minmerit = 1.0;
    //TODO: Question: Should we use LLL_FPZZFlex instead of the LLL implementation to have access to the length of the shortest vector?
    for (auto it = m_coordRange.begin(); it != m_coordRange.end(); it++){
    	//calculate LLL projection
    	//CHANGES ON 01/30
       //BasisConstruction<Int>::projectionConstructionLLL(lat.getBasis(), m_projBasis, *it, lat.getModulo(), m_delta);
       BasisConstruction<Int>::projectMatrix(lat.getBasis(), m_projBasis, *it, 0);
       // Define IntLattice based on projected basis
       proj.setBasis(m_projBasis, lat.getModulo(), m_projBasis.NumCols());
       // Apply selected reduction technique
       if (this->m_reductionMethod == BKZBB || this->m_reductionMethod == BKZ) {
           this->m_red->redBKZ(proj.getBasis(), this->m_delta, this->m_blocksize);  
       } else if (this->m_reductionMethod == LLLBB || this->m_reductionMethod == LLL) {
           LLL_FPZZflex(proj.getBasis(), this->m_delta, proj.getBasis().NumRows(), proj.getBasis().NumCols(), this->m_b);    
       } else if (this->m_reductionMethod == PAIRBB) {
           this->m_red->redDieter(0);
       }
       // TODO: CHECK IF IT IS NECESSARY TO ADD DIMENSIONS TO THE REDUCED BASIS
       if (!m_doingBB) {
         // Use first basis vector as proxy for shortest length
	       NTL::conv(merit, sqrt(this->m_b[0]) / this->m_norma->getBound(proj.getDim()));
       } else {
         if (!m_red->shortestVector(proj)) return 0;
         merit = NTL::conv<double>(m_red->getMinLength() / this->m_norma->getBound(proj.getDim()));
       }
       //** OLD ** //
       // Calculate shortest vector
      // if (m_reductionMethod == BKZBB || m_reductionMethod == LLLBB || m_reductionMethod == PAIRBB) {
      //     if (!m_red->shortestVector(proj)) return 0;
      //     shortest = NTL::conv<double>(m_red->getMinLength());
      // } else {
      //     //TODO: Better NOT use this function but use LLL_FPZZFlex output
      //     shortest = proj.getShortestLengthBasis();
      // }
      // merit = shortest / m_norma->getBound(proj.getDim());
       if (merit < minmerit) minmerit = merit;
       if (merit == 0 || merit < m_lowbound || merit > m_highbound) return 0;
    }
    return minmerit;  
}


//=========================================================================
//template<typename Int>
//double FigureOfMeritM<Int>::computeMeritMNonSuccPrimal (IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj) {
//    double merit = 0;
//    double minmerit = 1.0;
//    double shortest = 0.0;
//    //TODO: Question: Should we use LLL_FPZZFlex instead of the LLL implementation to have access to the length of the shortest vector?
//    for (auto it = m_coordRange.begin(); it != m_coordRange.end(); it++){
//    	//calculate LLL projection
//       BasisConstruction<Int>::projectionConstructionLLL(lat.getBasis(), m_projBasis, *it, lat.getModulo(), m_delta); 
//       // Define IntLattice based on projected basis
//       proj.setBasis(m_projBasis, lat.getModulo(), m_projBasis.NumCols());
//       // Calculate shortest vector
//       if (m_reductionMethod == BKZBB || m_reductionMethod == LLLBB || m_reductionMethod == PAIRBB) {
//           if (!m_red->shortestVector(proj)) return 0;
//           shortest = NTL::conv<double>(m_red->getMinLength());
//       } else {
//           //TODO: Better NOT use this function but use LLL_FPZZFlex output
//           shortest = proj.getShortestLengthBasis();
//       }
//       merit = shortest / m_norma->getBound(proj.getDim());
//       if (merit < minmerit) minmerit = merit;
//       if (merit == 0 || merit < m_lowbound || merit > m_highbound) return 0;
//    }
//    return minmerit;  
//}

//=========================================================================
	
template class FigureOfMeritM<NTL::ZZ>;

} // end namespace LatticeTester

#endif
