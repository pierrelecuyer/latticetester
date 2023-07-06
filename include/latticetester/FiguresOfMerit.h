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
	 * Defines the type of the figure of merit
	 */	
	MeritType m_fom = MERITM; 
	
	/*
	 * If FoM is above this bound, then calculation of FoM is stopped
	 */	
	double m_highbound = 1; 
	
	/*
	 * If FoM is below this bound, then calculation of FoM is stopped
	 */	
	double m_lowbound = 0; 
	
	/*
	 * Figure of merit is calculated for the dual lattice if dual is true
	 */
	bool m_dual = false; 
	
	/* 
	 * Defines the pre-reduction type. Indicates how the shortest vectors are computed 
	 * or estimated.
	 */
	PreReductionType m_reductionMethod = BKZ; 
	
	/*
	 * delta-Parameter for LLL or BKZ reduction
	 */
	double m_delta = 0.9; 
	
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
    * A function which allows to set maxNodesBB and passes the value to the reducer
     */
	void setWeights (int64_t & nodes) { m_maxNodesBB = nodes; m_red->maxNodesBB = nodes; }
	
	/*
	 * Variable decides BB algorithm is performed
	 */
	bool m_performBB = true; 
	
	/*
	 * If set to true successive corrdinates a considered first when calculating FigureOfMerit
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
	 * Variable used to store the m of the RNG
	 */
	NTL::ZZ m_m; 
	
	/*
	 * Decides if the primal lattice shall be read out directly because 
	 * it is not necessary to apply a projection construction
	 */
	bool m_ReadOutPrimal = false;
	
	/*
	 * Variable containing the Weights for the FoM 
	 */
	Weights *m_weights; 
	        
    /*
     * Initizalization of a FiguresOfMerit object. Puts a weight of 1 to everything. 
     * This way, the user is not forced to define weights himself or to pass
     * a corresponding object. Moreover, the parameter 'maxDim' determines the 
     * maximal dimension of the reducer object which is created (once) upon
     * initizalization.
     */
    FiguresOfMerit(Weights & w, Reducer<Int, Real> & red) {
    	m_weights = &w; 
    	m_red = &red;
	}
    
	
	/*
	 * This function calculates the Figure of Merit M or Q based 
	 * on the chosen MeritType of a given lattice 'lat'. The vector 't'
	 * defines the set of dimensions for which the figure of merit
	 * is calculated needs to be passed as second input variable.  
	 * The reducer object needs to be passed to avoid
	 * creating these objects inside the FiguresOfMerit class.
	 * The IntLattice object 'proj' is needed for saving projections.
	 * The value 0 is returned if an error occurs while calculating 
	 * the shortest vector.
	 */
	double computeMerit(IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj, const IntVec & t);
	
	/*
	 * This function calculates the Figure of Merit M and should be
	 * called by the user. The parameters have the same meaning as for
	 * computeMerit. The value 0 is returned if an error occurs
	 * while calculating the shortest vector.
	 */
	double computeMeritM(IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj, const IntVec & t);
	
	/*
	 * This function calculates the Figure of Merit Q and should be
	 * called by the user. The parameters have the same meaning as for
	 * computeMerit. The value 0 is returned if an error occurs
	 * while calculating the shortest vector.
	 */
	double computeMeritQ(IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj, const IntVec & t);	
	
	/*
	 * This function calculates the Figure of Merit for a single projection
	 * of a given lattice 'lat'. The variable 'Coord' sets the coordinates
	 * of the projection to use. The variable 'useLatBasis' indicates if the
	 * primal / dual basis stored in the Intlattice 'lat' can be directly used 
	 * because it is already upper triangular. The value 0 is returned if 
	 * an error occurs while calculating the shortest vector.
	 */
	double computeMeritProj(IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj, 
			const Coordinates & Coord, bool useLatBasis);
	
	/*
	 * TODO
	 */
	double computeMeritProjDual(IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj, 
			const Coordinates & Coord, bool useLatBasis);


	/*
	 * This functions calculates the Figure of Merit for all projections 
	 * consisting of successive coordinates of the forms 
	 * {1, 2, ..., low} to {1, 2, ..., upp}
	 * Note that upp > low must hold. The value 0 is returned if 
	 * an error occurs while calculating the shortest vector.
	 */
	double computeMeritSucc(IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj, 
			const int64_t & low, const Int & upp);
	
	/*
	 * This functions calculates the Figure of Merit for all projections 
	 * consisting of non-successive coordinates. The value 0 is returned if 
	 * an error occurs while calculating the shortest vector.
	 */
	double computeMeritVec(IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj, const IntVec & t);
	/*
	 * Stores the matrix of the projection basis
	 */
	IntMat m_projBasis; 
	
	/*
	 * CoordinateSets object is used to store sets of coordinates
	 */
    CoordinateSets::FromRanges *m_CoordRange;  
    
    /*
     * Normalizer for storing normalizing values in FiguresOfMerit class
     */
    Normalizer *m_norma;

    /*
     * A function which allows to set set the normalizer to user-specified value.
    */
    void setNormalizer (Normalizer & norm) { m_norma  = &norm; }
    
    /*
     * Reducer object used for finding the shortest vector of a projection
    */
    Reducer<Int, Real> *m_red;            
};

//============================================================================
// Implementation

template<typename Int>
double FiguresOfMerit<Int>::computeMerit(IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj, 
		const IntVec & t) {
   double merit = 0;
   double minmerit = 1.0;
   int64_t low_dim;
   
   if (m_succCoordFirst == true) { 
	   low_dim = t.length();
	   minmerit = computeMeritSucc(lat, proj, low_dim, t[0]);
	   if (minmerit == 0) return minmerit;
	   if (minmerit < m_lowbound) return minmerit;
	   if (minmerit > m_highbound) return minmerit;
   }   
   
   merit = computeMeritVec(lat, proj, t);
   if (merit < minmerit) minmerit = merit;
   if (merit == 0) return merit;
   if (minmerit < m_lowbound) return minmerit;
   if (minmerit > m_highbound) return minmerit;

   if (m_succCoordFirst == false) {
	   low_dim = t.length();
	   merit = computeMeritSucc(lat, proj, low_dim, t[0]);
	   if (merit < minmerit) minmerit = merit;
	   if (merit == 0) return merit;
	   if (minmerit < m_lowbound) return minmerit;
	   if (minmerit > m_highbound) return minmerit;
   } 
   return minmerit; 
}

//=========================================================================

template<typename Int>
double FiguresOfMerit<Int>::computeMeritM(IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj, 
		const IntVec & t) {
   double out;
   m_fom = MERITM;
   out = computeMerit(lat, proj, t);
   m_m = lat.getModulo();
   return out;
}

//=========================================================================

template<typename Int>
double FiguresOfMerit<Int>::computeMeritQ(IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj, 
		const IntVec & t) {
   double out;
   m_fom = MERITQ;
   out = computeMerit(lat, proj, t);
   m_m = lat.getModulo();
   return out;
}

//=========================================================================

template<typename Int>
double FiguresOfMerit<Int>::computeMeritProj(IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj, 
		const Coordinates & Coord, bool useLatBasis) {
   double shortest, merit;
   merit = 0.0;
   m_m = lat.getModulo();

   if (useLatBasis) {
       m_projBasis = lat.getBasis();
   } else {
       BasisConstruction<Int>::projectionConstruction(lat.getBasis(), m_projBasis, Coord, m_m, m_pctype, m_delta); 
       IntMat projBasisDual;
       if (m_pctype == UPPERTRIPROJ) {
          BasisConstruction<Int>::mDualUpperTriangular(m_projBasis, projBasisDual, m_m);
       } else {                   
          BasisConstruction<Int>::mDualBasis(m_projBasis, projBasisDual, m_m);
       }      
   }
   proj = IntLattice<Int, Real> (m_projBasis, m_m, m_projBasis.NumCols());
   //double log_density=(double)(-log(abs(NTL::determinant(proj->getBasis())))); 
   //Normalizer* norma = new NormaBestLat(log_density, max_dim);
   proj.updateVecNorm();
   proj.sort(0); 
   if (m_reductionMethod == BKZ) {
      m_red->redBKZ(proj.getBasis(), m_delta, m_blocksize);  
   } else {
      m_red->redLLLNTL(proj.getBasis(), m_delta);  
   } 
   if (m_fom == MERITQ) {
	   if (m_performBB == true) {if (!m_red->reductMinkowski(proj, 0)) return 0;}
         merit = NTL::conv<double>(m_red->getMinLength()) / NTL::conv<double>(m_red->getMaxLength());
   } else {
	  if (m_performBB == true) {
          if (!m_red->shortestVector(proj)) return 0;
	  }
      shortest = NTL::conv<double>(m_red->getMinLength());
      // std::cout << norma->getBound((Coord).size()) << "\n";
      merit = m_weights->getWeight(Coord) * shortest/m_norma->getBound((Coord).size());
   }
   return merit;
}

//=========================================================================

template<typename Int>
double FiguresOfMerit<Int>::computeMeritProjDual(IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj, 
		const Coordinates & Coord, bool useLatBasis) {
   double shortest, merit;
   merit = 0.0;
   m_m = lat.getModulo();

   if (useLatBasis) {
       m_projBasis = lat.getDualBasis(); 
   } else {
       IntMat projBasisDual;
       if (m_ReadOutPrimal == true) {
	      m_projBasis = lat.getBasis();
       } else {
		  BasisConstruction<Int>::projectionConstruction(lat.getBasis(), m_projBasis, Coord, m_m, m_pctype, m_delta); 
       }
       if (m_pctype == UPPERTRIPROJ) {
          BasisConstruction<Int>::mDualUpperTriangular(m_projBasis, projBasisDual, m_m);
       } else {                   
           BasisConstruction<Int>::mDualBasis(m_projBasis, projBasisDual, m_m);
       }                         
       m_projBasis = projBasisDual;
   }		   
   proj = IntLattice<Int, Real> (m_projBasis, m_m, m_projBasis.NumCols());
   //double log_density=(double)(-log(abs(NTL::determinant(proj->getBasis())))); 
   //Normalizer* norma = new NormaBestLat(log_density, max_dim);
   proj.updateVecNorm();
   proj.sort(0); 
   if (m_reductionMethod == BKZ) {
      m_red->redBKZ(proj.getBasis(), m_delta, m_blocksize);  
   } else {
      m_red->redLLLNTL(proj.getBasis(), m_delta);  
   } 
   if (m_fom == MERITQ) {
	   if (m_performBB == true) {if (!m_red->reductMinkowski(proj, 0)) return 0;}
         merit = NTL::conv<double>(m_red->getMinLength()) / NTL::conv<double>(m_red->getMaxLength());
   } else {
	  if (m_performBB == true) {
          if (!m_red->shortestVector(proj)) return 0;
	  }
      shortest = NTL::conv<double>(m_red->getMinLength());
      // std::cout << norma->getBound((Coord).size()) << "\n";
      merit = m_weights->getWeight(Coord) * shortest/m_norma->getBound((Coord).size());
   }
   return merit;
}

//=========================================================================

template<typename Int>
double FiguresOfMerit<Int>::computeMeritSucc(IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj, 
		const int64_t & low, const Int & upp) {
   Coordinates Coord;	
   double merit = 0;
   double minmerit = 1.0;
   double dim = lat.getBasis().NumCols();
   Coord.clear();
   for (int j = 0; j < low + 1; j++) Coord.insert(j+1);
   lat.buildBasis(low+1);
   if (m_dual) { merit = computeMeritProjDual(lat, proj, Coord, true);
   } else merit = computeMeritProj(lat, proj, Coord, true);
   if (merit == 0) return merit;
   for (int j = low +1; j < upp; j++)
   {
       Coord.insert(j+1);
       if (m_incDualOnly == true) {
    	   lat.incDimDual();
       } else {
    	   lat.incDim();
       }
     //lat.buildBasis(j+1);
       if (m_dual) { merit = computeMeritProjDual(lat, proj, Coord, true);
       } else merit = computeMeritProj(lat, proj, Coord, true);
       if (merit < minmerit) minmerit = merit;
	   if (merit == 0) return merit;
       if (merit < m_lowbound) return minmerit;
       if (merit > m_highbound) return minmerit;
   }   
   if (upp < dim)
      lat.buildBasis(dim);
   return minmerit;
}

//=========================================================================

template<typename Int>
double FiguresOfMerit<Int>::computeMeritVec(IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj, 
		const IntVec & t) {
	Coordinates Coord;	
	double merit = 0;
	double minmerit = 1.0;
	int64_t max_dim, min_dim;
	max_dim = lat.getDim();  
	if (m_projectStationary) {
       min_dim = 2;
	} else min_dim = 1;

	for (int i = 1; i < t.length(); i++) {
       NTL::conv(max_dim, t[i]);
	   CoordinateSets::FromRanges m_CoordRange(i, i, min_dim, max_dim, m_projectStationary);  
       for (auto it = m_CoordRange.begin(); it != m_CoordRange.end(); it++){
          Coord = *it;
          if (m_dual) { merit = computeMeritProjDual(lat, proj, Coord, false);
          } else merit = computeMeritProj(lat, proj, Coord, false);
          if (merit < minmerit) minmerit = merit;
   	      if (merit == 0) return merit;
	      if (merit < m_lowbound) return minmerit;
	      if (merit > m_highbound) return minmerit;
       }
    }
	return minmerit;  
}

//=========================================================================
	
template class FiguresOfMerit<NTL::ZZ>;
//template class FiguresOfMerit<std::int64_t>;

} // end namespace LatticeTester

#endif
