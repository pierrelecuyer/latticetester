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
    * A function which allows to set maxNodesBB and passes the value to the reducer
     */
	void setmaxNodesBB (int64_t & nodes) { m_maxNodesBB = nodes; m_red->maxNodesBB = nodes; }

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
	 * Variable used to store the modulo m of the RNG.
	 */
	Int m_m; 
	
	/*
	 * Decides if the primal lattice shall be read out directly because 
	 * it is not necessary to apply a projection construction.
	 */
	bool m_ReadOutPrimal = false;
	
	/*
	 * Variable containing the Weights for the FoM 
	 */
	Weights *m_weights; 
	        
    /*
     * Constructor of a FiguresOfMerit object. Puts a weight of 1 to everything. 
     * This way, the user is not forced to define weights himself or to pass
     * a corresponding object. Moreover, the parameter 'maxDim' determines the 
     * maximal dimension of the reducer object which is created (once) upon
     * initialization.
     */
    FiguresOfMerit(Weights & w, Reducer<Int, Real> & red) {
    	m_weights = &w;     // Should it be optional?     **********
    	m_red = &red;
	}
	
	/*
	 * This function calculates the Figure of Merit M of a given lattice 'lat'
	 * and should be called by the user. The vector 't' defines the set of 
	 * dimensions for which the figure of merit is calculated.
	 * The IntLattice object 'proj' is needed for saving projections.
	 *  The value 0 is returned if an error occurs
	 * while calculating the shortest vector.
	 */
	double computeMeritM(IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj, const IntVec & t);
	
	/*
	 * This function calculates the Figure of Merit Q and should be
	 * called by the user. The parameters have the same meaning as for
	 * computeMerit. The value 0 is returned if an error occurs
	 * while calculating the shortest vector.
	 * DOES CURRENTLY NOT WORK
	 */
	double computeMeritQ(IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj, const IntVec & t);	
	
	/*
	 * This function calculates the Figure of Merit for a single projection 'proj'. 
	 * The variable 'Coord' sets the coordinates of the projection to use. 
	 * The value 0 is returned if an error occurs during the calculation of the shortest vector.
	 */
	double computeMeritProj(IntLattice<Int, Real> & proj, const Coordinates & Coord);
	
	/*
	 * This function calculates the Figure of Merit for a given lattice 'lat'. 
	 */
	double computeMeritLat(IntLatticeExt<Int, Real> & lat);
	
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
	double computeMeritNonSucc(IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj, const IntVec & t);
	
	/*
	 * Stores the matrix of the projection basis
	 */
	IntMat m_projBasis; 
	
	/* 
	 * This function calculates the variable m_projBasis based on the 
	 * desired method to derive it. The coordinates of the projection are given by
	 * 'Coord;. The parameter 'useLatBasis" indicates whether m_projBasis shall
	 * be read out directly from the IntLattice 'lat'. 
	*/

	void getProjBasis(IntLatticeExt<Int, Real> & lat, const Coordinates & Coord, bool useLatBasis);
	
	/*
	 * Same as for getProjBasis but for the dual.
	 */
	void getProjBasisDual(IntLatticeExt<Int, Real> & lat, const Coordinates & Coord, bool useLatBasis);
	
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
    
    /*
     * This function calculates the Figure of Merit M using method A: 
     * This method uses incDimBasis to increase the basis at each step, 
     * as explained in my 1997 paper.  
    */
    double computeMeritMSucc_MethodA(IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj, const IntVec & t);
    
    /*
     * This function calculates the Figure of Merit M using method B:
     * if we do not maintain the dual, incDim will only give you a basis for 
     * the primal. Then you will need to do something else to get the dual basis, 
     * in case your FOM is in the dual. You may compute a triangular basis from 
     * the primal basis, then get the m-dual of it, which will of course be triangular, 
     * then reduce this triangular dual basis via LLL or BKZ, then apply BB.
    */
    double computeMeritMSucc_MethodB(IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj, const IntVec & t);

    /*
     * This function calculates the Figure of Merit M using method C:
     * With this one, you take the primal basis, then compute the m-dual of 
     * it using the general mDualBasis method, then reduce this m-dual basis 
     * via LLL or BKZ, then apply BB.
    */
    double computeMeritMSucc_MethodC(IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj, const IntVec & t);
        
    /*
    * This function calculates the Figure of Merit M using method D:
    * The simplest is to just reconstruct the basis for the m-dual lattice from 
    * scratch as in section 4.1.1 of the guide each time we increase the dimension t by 1.  
    * Then we apply LLL or BKZ + BB directly to this m-dual basis. 
    */
    double computeMeritMSucc_MethodD(IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj, const IntVec & t);
            
    /*
    * This function calculates the Figure of Merit M using method E:
    * That is, we just update the m-dual basis, then apply LLL or BKZ to it, then BB
    */
    double computeMeritMSucc_MethodE(IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj, const IntVec & t);
      
};

//============================================================================
// Implementation

//template<typename Int>
//double FiguresOfMerit<Int>::computeMerit(IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj, 
//		const IntVec & t) {
//   double merit = 0;
//   double minmerit = 1.0;
//   int64_t low_dim;
//   
//   if (m_succCoordFirst == true) { 
//	   low_dim = t.length();
//	   minmerit = computeMeritSucc(lat, proj, low_dim, t[0]);
//	   if (minmerit == 0) return minmerit;
//	   if (minmerit < m_lowbound) return minmerit;
//	   if (minmerit > m_highbound) return minmerit;
//   }   
//   
//   merit = computeMeritNonSucc(lat, proj, t);
//   if (merit < minmerit) minmerit = merit;
//   if (merit == 0) return merit;
//   if (minmerit < m_lowbound) return minmerit;
//   if (minmerit > m_highbound) return minmerit;
//
//   if (m_succCoordFirst == false) {
//	   low_dim = t.length();
//	   merit = computeMeritSucc(lat, proj, low_dim, t[0]);
//	   if (merit < minmerit) minmerit = merit;
//	   if (merit == 0) return merit;
//	   if (minmerit < m_lowbound) return minmerit;
//	   if (minmerit > m_highbound) return minmerit;
//   } 
//   return minmerit; 
//}

//=========================================================================

template<typename Int>
double FiguresOfMerit<Int>::computeMeritM(IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj, 
		const IntVec & t) {
   double out;
   m_fom = MERITM;
   m_m = lat.getModulo();
   //out = computeMerit(lat, proj, t);
   //return out;
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
   
   merit = computeMeritNonSucc(lat, proj, t);
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
double FiguresOfMerit<Int>::computeMeritQ(IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj, 
		const IntVec & t) {
   double out;
   m_fom = MERITQ;
   m_m = lat.getModulo();
   //out = computeMerit(lat, proj, t);
   //return out;
   return 0;
}

//=========================================================================

template<typename Int>
void FiguresOfMerit<Int>::getProjBasis(IntLatticeExt<Int, Real> & lat, const Coordinates & Coord, bool useLatBasis) {
	   
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
}

template<typename Int>
void FiguresOfMerit<Int>::getProjBasisDual(IntLatticeExt<Int, Real> & lat, const Coordinates & Coord, bool useLatBasis) {
	   
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
}


//=========================================================================

template<typename Int>
double FiguresOfMerit<Int>::computeMeritProj(IntLattice<Int, Real> & proj, const Coordinates & Coord) {
   double shortest, merit;
   merit = 0.0;
   
   proj = IntLattice<Int, Real> (m_projBasis, m_m, m_projBasis.NumCols());
   proj.updateVecNorm();
   proj.sort(0); 
   if (m_reductionMethod == BKZBB || m_reductionMethod == BKZ) {
      m_red->redBKZ(proj.getBasis(), m_delta, m_blocksize);  
   } else {
      m_red->redLLLNTL(proj.getBasis(), m_delta);  
   } 
   //std::cout << proj.getBasis();
   if (m_reductionMethod == BKZBB || m_reductionMethod == LLLBB || m_reductionMethod == PAIRBB) {
       if (!m_red->shortestVector(proj)) return 0;
   }
   shortest = NTL::conv<double>(m_red->getMinLength());
   //std::cout << shortest << "\n";
   // std::cout << norma->getBound((Coord).size()) << "\n";
   merit = m_weights->getWeight(Coord) * shortest/m_norma->getBound((Coord).size());
   return merit;
}

//=========================================================================

template<typename Int>
double FiguresOfMerit<Int>::computeMeritLat(IntLatticeExt<Int, Real> & lat) {
   double shortest, merit;
   merit = 0.0;
   m_m = lat.getModulo();
   Coordinates Coord;
   
   lat.updateVecNorm();
   lat.sort(0); 
   if (m_reductionMethod == BKZBB || m_reductionMethod == BKZ) {
      m_red->redBKZ(lat.getBasis(), m_delta, m_blocksize);  
   } else {
      m_red->redLLLNTL(lat.getBasis(), m_delta);  
   } 
   //std::cout << lat.getBasis();
   if (m_reductionMethod == BKZBB || m_reductionMethod == LLLBB || m_reductionMethod == PAIRBB) {
          if (!m_red->shortestVector(lat)) return 0;
   }
   shortest = NTL::conv<double>(m_red->getMinLength());
   //std::cout << shortest << "\n";
   // std::cout << norma->getBound((Coord).size()) << "\n";
   for (int j = 0; j < lat.getBasis().NumCols(); j++) Coord.insert(j+1);
   merit = m_weights->getWeight(Coord) * shortest/m_norma->getBound((Coord).size());
   return merit;
}

//=========================================================================

template<typename Int>
double FiguresOfMerit<Int>::computeMeritSucc(IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj, 
		const int64_t & low, const Int & upp) {
   Coordinates Coord;	
   Coord.clear();
   double merit = 0;
   double minmerit = 1.0;
   double dim = lat.getBasis().NumCols();
   for (int j = 0; j < low + 1; j++) Coord.insert(j+1);
   lat.buildBasis(low+1);
//   if (m_fomInDual) { getProjBasisDual(lat, Coord, true); 
//   } else getProjBasis(lat, Coord, true); 
//   merit = computeMeritProj(proj, Coord);
   if (m_fomInDual) lat.dualize();
   merit = computeMeritLat(lat);
   BasisConstruction<Int>::mDualBasis(lat.getBasis(), lat.getDualBasis(), m_m);
   if (m_fomInDual) lat.dualize();
   if (merit == 0) return merit;
   for (int j = low +1; j < upp; j++)
   {
	   // Covers methods B and E
       Coord.insert(j+1);
       if (m_incDualOnly == true) {
    	   lat.incDimDualBasis();
       } else {
    	   lat.incDimBasis();
       }
     //lat.buildBasis(j+1);
//       if (m_fomInDual) { getProjBasisDual(lat, Coord, true); 
//       } else getProjBasis(lat, Coord, true); 
//       merit = computeMeritProj(proj, Coord);
       if (m_fomInDual) lat.dualize();
       merit = computeMeritLat(lat);
       BasisConstruction<Int>::mDualBasis(lat.getBasis(), lat.getDualBasis(), m_m);
       if (m_fomInDual) lat.dualize();
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
double FiguresOfMerit<Int>::computeMeritNonSucc(IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj, 
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
          if (m_fomInDual) { getProjBasisDual(lat, Coord, false); 
          } else getProjBasis(lat, Coord, false); 
          merit = computeMeritProj(proj, Coord);
          if (merit < minmerit) minmerit = merit;
   	      if (merit == 0) return merit;
	      if (merit < m_lowbound) return minmerit;
	      if (merit > m_highbound) return minmerit;
       }
    }
	return minmerit;  
}

//=========================================================================

//template<typename Int>
//double FiguresOfMerit<Int>::computeMeritMSucc_MethodA(IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj, 
//		const IntVec & t) {
//   m_m = lat.getModulo();
//   Coordinates Coord;	
//   Coord.clear();
//   double merit = 0;
//   double minmerit = 1.0;
//   double shortest;
//   int64_t low;
//   
//   low = t.length();
//   double dim = lat.getBasis().NumCols();
//   for (int j = 0; j < low + 1; j++) Coord.insert(j+1);
//   lat.buildBasis(low+1);
//   if (m_fomInDual) { getProjBasisDual(lat, Coord, true); // Reads out dual basis directly from lattice
//   } else getProjBasis(lat, Coord, true); // Reads out basis directly from lattice
//   //merit = computeMeritLat(lat, Coord);
//   merit = computeMeritProj(lat, proj, Coord);
//
//   if (t[0] < dim)
//      lat.buildBasis(dim);	     
//
//   return minmerit; 
//}

template<typename Int>
double FiguresOfMerit<Int>::computeMeritMSucc_MethodA(IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj, 
		const IntVec & t) {
   m_m = lat.getModulo();
   Coordinates Coord;	
   Coord.clear();
   double merit = 0;
   double minmerit = 1.0;
   int64_t low;
   m_m = lat.getModulo();
   
   low = t.length();
   for (int j = 0; j < low + 1; j++) Coord.insert(j+1);
   lat.buildBasis(low+1);
//   if (m_fomInDual) { getProjBasisDual(lat, Coord, false); 
//   } else getProjBasis(lat, Coord, false); 
//   merit = computeMeritProj(proj, Coord);
   if (m_fomInDual) lat.dualize();
   merit = computeMeritLat(lat);
   BasisConstruction<Int>::mDualBasis(lat.getBasis(), lat.getDualBasis(), m_m);
   if (m_fomInDual) lat.dualize();
   if (merit == 0) return merit;
   for (int j = low +1; j < t[0]; j++)
   {
       Coord.insert(j+1);
       lat.incDimBasis(); // Increase the dimension of the lattice (and potentially the dual) 
//       if (m_fomInDual) { getProjBasisDual(lat, Coord, false); 
//       } else getProjBasis(lat, Coord, false); 
//       merit = computeMeritProj(proj, Coord);
       if (m_fomInDual) lat.dualize();
       merit = computeMeritLat(lat);
       BasisConstruction<Int>::mDualBasis(lat.getBasis(), lat.getDualBasis(), m_m);
       if (m_fomInDual) lat.dualize();
       if (merit < minmerit) minmerit = merit;
 	   if (merit == 0) return merit;
       if (merit < m_lowbound) return minmerit;
       if (merit > m_highbound) return minmerit;
   }   
//   if (t[0] < dim)
//      lat.buildBasis(dim);	     

   return minmerit; 
}

//=========================================================================

template<typename Int>
double FiguresOfMerit<Int>::computeMeritMSucc_MethodB(IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj, 
		const IntVec & t) {
   m_m = lat.getModulo();
   Coordinates Coord;	
   Coord.clear();
   double merit = 0;
   double minmerit = 1.0;
   int64_t low;
   m_m = lat.getModulo();
   
   low = t.length();
   double dim = lat.getBasis().NumCols();
   for (int j = 0; j < low + 1; j++) Coord.insert(j+1);
   lat.buildBasis(low+1);
   m_pctype = UPPERTRIPROJ; // Use upper triangular algorithm to calculate the dual
//   if (m_fomInDual) { getProjBasisDual(lat, Coord, false); 
//   } else getProjBasis(lat, Coord, false); 
//   merit = computeMeritProj(proj, Coord);
   if (m_fomInDual) lat.dualize();
   merit = computeMeritLat(lat);
   BasisConstruction<Int>::mDualBasis(lat.getBasis(), lat.getDualBasis(), m_m);
   if (m_fomInDual) lat.dualize();
   if (merit == 0) return merit;
   for (int j = low +1; j < t[0]; j++)
   {
       Coord.insert(j+1);
       lat.incDimBasis(); // Increase the dimension of the lattice (and potentially the dual) 
//       if (m_fomInDual) { getProjBasisDual(lat, Coord, false); 
//       } else getProjBasis(lat, Coord, false); 
//       merit = computeMeritProj(proj, Coord);
       if (m_fomInDual) lat.dualize();
       merit = computeMeritLat(lat);
       BasisConstruction<Int>::mDualBasis(lat.getBasis(), lat.getDualBasis(), m_m);
       if (m_fomInDual) lat.dualize();
       if (merit < minmerit) minmerit = merit;
 	   if (merit == 0) return merit;
       if (merit < m_lowbound) return minmerit;
       if (merit > m_highbound) return minmerit;
   }   
   if (t[0] < dim)
      lat.buildBasis(dim);	     

   return minmerit; 
}

//=========================================================================

template<typename Int>
double FiguresOfMerit<Int>::computeMeritMSucc_MethodC(IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj, 
		const IntVec & t) {
   m_m = lat.getModulo();
   Coordinates Coord;	
   Coord.clear();
   double merit = 0;
   double minmerit = 1.0;
   int64_t low;
   m_m = lat.getModulo();
   
   low = t.length();
   double dim = lat.getBasis().NumCols();
   for (int j = 0; j < low + 1; j++) Coord.insert(j+1);
   lat.buildBasis(low+1);
   m_pctype = LLLPROJ; // Use the general method to calculate the dual
//   if (m_fomInDual) { getProjBasisDual(lat, Coord, false); 
//   } else getProjBasis(lat, Coord, false); 
//   merit = computeMeritProj(proj, Coord);
   if (m_fomInDual) lat.dualize();
   merit = computeMeritLat(lat);
   BasisConstruction<Int>::mDualBasis(lat.getBasis(), lat.getDualBasis(), m_m);
   if (m_fomInDual) lat.dualize();
   if (merit == 0) return merit;
   for (int j = low +1; j < t[0]; j++)
   {
       Coord.insert(j+1);
       lat.incDimBasis(); // Increase the dimension of the lattice (and potentially the dual) 
//       if (m_fomInDual) { getProjBasisDual(lat, Coord, false); 
//       } else getProjBasis(lat, Coord, false); 
//       merit = computeMeritProj(proj, Coord);
       if (m_fomInDual) lat.dualize();
       merit = computeMeritLat(lat);
       BasisConstruction<Int>::mDualBasis(lat.getBasis(), lat.getDualBasis(), m_m);
       if (m_fomInDual) lat.dualize();
       if (merit < minmerit) minmerit = merit;
 	   if (merit == 0) return merit;
       if (merit < m_lowbound) return minmerit;
       if (merit > m_highbound) return minmerit;
   }   
   if (t[0] < dim)
      lat.buildBasis(dim);	     

   return minmerit; 
}

//=========================================================================

template<typename Int>
double FiguresOfMerit<Int>::computeMeritMSucc_MethodD(IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj, 
		const IntVec & t) {
   m_m = lat.getModulo();
   Coordinates Coord;	
   Coord.clear();
   double merit = 0;
   double minmerit = 1.0;
   int64_t low;
   m_m = lat.getModulo();
   
   low = t.length();
   double dim = lat.getBasis().NumCols();
   for (int j = 0; j < low + 1; j++) Coord.insert(j+1);
   lat.buildDualBasis(low+1);
//   if (m_fomInDual) { getProjBasisDual(lat, Coord, true); 
//   } else getProjBasis(lat, Coord, true); 
//   merit = computeMeritProj(proj, Coord);
   if (m_fomInDual) lat.dualize();
   merit = computeMeritLat(lat);
   if (m_fomInDual) lat.dualize();
   if (merit == 0) return merit;
   for (int j = low +1; j < t[0]; j++)
   {
       Coord.insert(j+1);
       lat.buildDualBasis(j+1); // Build (only) the dual basis from scratch 
//       if (m_fomInDual) { getProjBasisDual(lat, Coord, true); 
//       } else getProjBasis(lat, Coord, true); 
//       merit = computeMeritProj(proj, Coord);
       if (m_fomInDual) lat.dualize();
       merit = computeMeritLat(lat);
       if (m_fomInDual) lat.dualize();
       if (merit < minmerit) minmerit = merit;
 	   if (merit == 0) return merit;
       if (merit < m_lowbound) return minmerit;
       if (merit > m_highbound) return minmerit;
   }   
   if (t[0] < dim)
      lat.buildBasis(dim);	     

   return minmerit; 
}

//=========================================================================

template<typename Int>
double FiguresOfMerit<Int>::computeMeritMSucc_MethodE(IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> & proj, 
		const IntVec & t) {
   m_m = lat.getModulo();
   Coordinates Coord;	
   Coord.clear();
   double merit = 0;
   double minmerit = 1.0;
   int64_t low;
   m_m = lat.getModulo();
   
   low = t.length();
   double dim = lat.getBasis().NumCols();
   for (int j = 0; j < low + 1; j++) Coord.insert(j+1);
   lat.buildDualBasis(low+1);
//   if (m_fomInDual) { getProjBasisDual(lat, Coord, true); 
//   } else getProjBasis(lat, Coord, true); 
//   merit = computeMeritProj(proj, Coord);
   if (m_fomInDual) lat.dualize();
   merit = computeMeritLat(lat);
   if (m_fomInDual) lat.dualize();
   if (merit == 0) return merit;
   for (int j = low +1; j < t[0]; j++)
   {
       Coord.insert(j+1);
       lat.incDimDualBasis(); // Build (only) the dual basis from scratch 
//       if (m_fomInDual) { getProjBasisDual(lat, Coord, true); 
//       } else getProjBasis(lat, Coord, true); 
//       merit = computeMeritProj(proj, Coord);
       if (m_fomInDual) lat.dualize();
       merit = computeMeritLat(lat);
       if (m_fomInDual) lat.dualize();
       if (merit < minmerit) minmerit = merit;
 	   if (merit == 0) return merit;
       if (merit < m_lowbound) return minmerit;
       if (merit > m_highbound) return minmerit;
   }   
   if (t[0] < dim)
      lat.buildBasis(dim);	     
	     

   return minmerit; 
}

//=========================================================================
	
template class FiguresOfMerit<NTL::ZZ>;
//template class FiguresOfMerit<std::int64_t>;

} // end namespace LatticeTester

#endif
