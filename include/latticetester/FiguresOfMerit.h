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
	MeritType fom = MERITM; 
	
	/*
	 * If true (the default value), we maximize the FOM. 
	 * With false (which is unusual), we minimize the FOM. This can be used
	 * to find bad generators.
	 * !!! NOT YET IMPLEMENTED!!!
	 */
	bool maximize = true;
	
	/*
	 * If FoM is below this bound, then calculation of FoM is stopped
	 */	
	double lowbound = 0; // 
	
	/*
	 * Figure of merit is calculated for the dual lattice if dual is true
	 */
	bool dual = false; 
	
	/* 
	 * The norm used to measure vector lengths. Can be L2NORM (the default value) 
	 * or L1NORM.
	 * !!! NOT YET IMPLEMENTED !!!
	 */
	NormType norm = L2NORM;
	
	/*
	 * The choice of normalization for the lengths of the shortest nonzero vectors, 
	 * used in the computation of the FOMs based on the spectral test with the L2 norm. 
	 * In all cases except for the Roger’s bound, the normalization constants are 
	 * available only up to t = 48 dimensions. 
	 * !!! NOT YET IMPLEMENTED !!!
	 */
	CriterionType normalizer;
	
	/* 
	 * Defines the pre-reduction type. Indicates how the shortest vectors are computed 
	 * or estimated.
	 */
	PreReductionType reductionMethod = BKZ; 
	
	/*
	 * delta-Parameter for LLL or BKZ reduction
	 */
	double delta = 0.9; 
	
	/*
	 * Blocksize of BKZ algorithm
	 */
	int64_t blocksize = 10; 
	
	/* 
	 * An integer giving the maximum number of nodes to be examined in any given BB procedure. 
	 * When that value is exceeded, the BB is stopped and the generator is discarded. 
	 * The number of discarded generators is given in the results. A small value of this 
	 * maxnodesBB parameter can make the program run faster (sometimes much faster), 
	 * permitting to examine more generators, but will increase the chances of rejecting good 
	 * generators. The default value is 10^8.
	 * !!! NOT YET IMPLEMENTED !!!
	 */
	int64_t maxnodesBB = 100000000;
	
	/*
	 * Variable decides BB algorithm is performed
	 */
	bool performBB = true; 
	
	/*
	 * If set to true successive corrdinates a considered first when calculating FigureOfMerit
	 */
	bool succCoordFirst = false;
	
	/*
	 * Sets if first variable shall always be included for the non-successive coordinates
	 */
	bool projectStationary = true;
		
	/*
	 * Type of projection construction
	 */	
	ProjConstructType pctype = UPPERTRIPROJ;
	
	/*
	 * Variable used to store the m of the RNG
	 */
	NTL::ZZ m; 
	
	/*
	 * Decides if the primal lattice shall be read out directly because 
	 * it is not necessary to apply a projection construction
	 */
	bool ReadOutPrimal = false;
	
	/*
	 * Variable containing the Weights for the FoM 
	 */
	WeightsUniform *weights; 
	        
    /*
     * Initialization of a FiguresOfMerit object
     */
    FiguresOfMerit(Int & i) {
	    weights = new WeightsUniform(1.0); // This just puts a weight of 1 to everything          
        //This is to initialize proj in order to avoid recreation
	}
	
	/*
	 * This function calculates the Figure of Merit M or Q based 
	 * on the chosen MeritType of a given lattice 'lat'. The vector 't'
	 * defines the set of dimensions for which the figure of merit
	 * is calculated needs to be passed as second input variable.  
	 * The normalizer and reducer object need to be passed to avoid
	 * creating these objects inside the FiguresOfMerit class.
	 * The IntLattice object 'proj' is needed for saving projections.
	 */
	double computeMerit(IntLatticeExt<Int, Real> & lat, Normalizer & norm, 
            Reducer<Int, Real> & red, IntLattice<Int, Real> & proj, const IntVec & t);
	
	/*
	 * This function calculates the Figure of Merit M and should be
	 * called by the user. The parameters have the same meaning as for
	 * computeMerit.
	 */
	double computeMeritM(IntLatticeExt<Int, Real> & lat, Normalizer & norm, 
            Reducer<Int, Real> & red, IntLattice<Int, Real> & proj, const IntVec & t);
	
	/*
	 * This function calculates the Figure of Merit Q and should be
	 * called by the user. The parameters have the same meaning as for
	 * computeMerit.
	 */
	double computeMeritQ(IntLatticeExt<Int, Real> & lat, Normalizer & norm, 
            Reducer<Int, Real> & red, IntLattice<Int, Real> & proj, const IntVec & t);	
	
	/*
	 * This function calculates the Figure of Merit for a single projection
	 * of a given lattice 'lat'. The variable 'Coord' sets the coordinates
	 * of the projection to use. The variable 'useLatBasis' indicates if the
	 * primal / dual basis stored in the Intlattice 'lat' can be directly used 
	 * because it is already upper triangular.
	 */
	double computeMeritProj(IntLatticeExt<Int, Real> & lat, Normalizer & norm, 
            Reducer<Int, Real> & red, IntLattice<Int, Real> & proj, const Coordinates & Coord, bool useLatBasis);

	/*
	 * This functions calculates the Figure of Merit for all projections 
	 * consisting of successive coordinates of the forms 
	 * {1, 2, ..., low} to {1, 2, ..., upp}
	 * Note that upp > low must hold 
	 */
	double computeMeritSucc(IntLatticeExt<Int, Real> & lat, Normalizer & norm, 
            Reducer<Int, Real> & red, IntLattice<Int, Real> & proj, const int64_t & low, const Int & upp);
	
	/*
	 * This functions calculates the Figure of Merit for all projections 
	 * consisting of non-successive coordinates.
	 */
	double computeMeritVec(IntLatticeExt<Int, Real> & lat, Normalizer & norm, 
            Reducer<Int, Real> & red, IntLattice<Int, Real> & proj, const IntVec & t);
	/*
	 * Stores the matrix of the projection basis
	 */
	IntMat projBasis; 
	
	/*
	 * CoordinateSets object is used to store sets of coordinates
	 */
    CoordinateSets::FromRanges *CoordRange;  	

};

//============================================================================
// Implementation

template<typename Int>
double FiguresOfMerit<Int>::computeMerit(IntLatticeExt<Int, Real> & lat, Normalizer & norm, 
        Reducer<Int, Real> & red, IntLattice<Int, Real> & proj, const IntVec & t) {
   double merit = 0;
   double minmerit = 1.0;
   Int maxDim;
   maxDim = lat.getDim();
   int64_t max_dim, low_dim;
   NTL::conv(max_dim, maxDim);
   Coordinates Coord;
   
   if (succCoordFirst == true) { 
	   low_dim = t.length();
	   minmerit = computeMeritSucc(lat, norm, red, proj, low_dim, t[0]);
	   if (minmerit < lowbound) return minmerit;
   }   
   
   merit = computeMeritVec(lat, norm, red, proj, t);
   if (merit < minmerit) minmerit = merit;
   if (minmerit < lowbound) return minmerit;

   if (succCoordFirst == false) {
	   low_dim = t.length();
	   merit = computeMeritSucc(lat, norm, red, proj, low_dim, t[0]);
	   if (merit < minmerit) minmerit = merit;
	   if (minmerit < lowbound) return minmerit;
   } 
   return minmerit; 
}

//=========================================================================

template<typename Int>
double FiguresOfMerit<Int>::computeMeritM(IntLatticeExt<Int, Real> & lat, Normalizer & norm, 
        Reducer<Int, Real> & red, IntLattice<Int, Real> & proj, const IntVec & t) {
   double out;
   fom = MERITM;
   out = computeMerit(lat, norm, red, proj, t);
   m = lat.getModulo();
   return out;
}

//=========================================================================

template<typename Int>
double FiguresOfMerit<Int>::computeMeritQ(IntLatticeExt<Int, Real> & lat, Normalizer & norm, 
        Reducer<Int, Real> & red, IntLattice<Int, Real> & proj, const IntVec & t) {
   double out;
   fom = MERITQ;
   out = computeMerit(lat, norm, red, proj, t);
   m = lat.getModulo();
   return out;
}

//=========================================================================

template<typename Int>
double FiguresOfMerit<Int>::computeMeritProj(IntLatticeExt<Int, Real> & lat, Normalizer & norm, 
        Reducer<Int, Real> & red, IntLattice<Int, Real> & proj, const Coordinates & Coord, bool useLatBasis) {
   double shortest, merit;
   merit = 0.0;
   m = lat.getModulo();

   if (dual == false) {
	   if (useLatBasis) {
           projBasis = lat.getBasis();
	   } else {
           BasisConstruction<Int>::projectionConstruction(lat.getBasis(), projBasis, Coord, m, pctype, delta); 
           IntMat projBasisDual;
           if (pctype == UPPERTRIPROJ) {
              BasisConstruction<Int>::mDualUpperTriangular(projBasis, projBasisDual, m);
           } else {                   
              BasisConstruction<Int>::mDualBasis(projBasis, projBasisDual, m);
           }      
	   }
   } else {
	   if (useLatBasis) {
           projBasis = lat.getDualBasis(); 
	   } else {
           IntMat projBasisDual;
		   if (ReadOutPrimal == true) {
			   projBasis = lat.getBasis();
		   } else {
			   BasisConstruction<Int>::projectionConstruction(lat.getBasis(), projBasis, Coord, m, pctype, delta); 
		   }
           if (pctype == UPPERTRIPROJ) {
              BasisConstruction<Int>::mDualUpperTriangular(projBasis, projBasisDual, m);
           } else {                   
              BasisConstruction<Int>::mDualBasis(projBasis, projBasisDual, m);
           }                         
           projBasis = projBasisDual;
		}
	 }		   
   proj = IntLattice<Int, Real> (projBasis, m, projBasis.NumCols());
   //double log_density=(double)(-log(abs(NTL::determinant(proj->getBasis())))); 
   //Normalizer* norma = new NormaBestLat(log_density, max_dim);
   proj.updateVecNorm();
   proj.sort(0); 
   if (reductionMethod == BKZ) {
      red.redBKZ(proj.getBasis(), delta, blocksize);  
   } else {
      red.redLLLNTL(proj.getBasis(), delta);  
   } 
   if (performBB == true) {
	  if (fom == MERITQ) {
         red.reductMinkowski(proj, 0);
         merit = NTL::conv<double>(red.getMinLength()) / NTL::conv<double>(red.getMaxLength());
      } else {
      red.shortestVector(proj);
      shortest = NTL::conv<double>(red.getMinLength());
      // std::cout << norma->getBound((Coord).size()) << "\n";
      merit = weights->getWeight(Coord) * shortest/norm.getBound((Coord).size());
      }
   } else { merit = 0.0;};
   return merit;

}

//=========================================================================

template<typename Int>
double FiguresOfMerit<Int>::computeMeritSucc(IntLatticeExt<Int, Real> & lat, Normalizer & norm, 
        Reducer<Int, Real> & red, IntLattice<Int, Real> & proj, const int64_t & low, const Int & upp) {
   Coordinates Coord;	
   double merit = 0;
   double minmerit = 1.0;
   double dim = lat.getBasis().NumCols();
   Coord.clear();
   for (int j = 0; j < low + 1; j++) Coord.insert(j+1);
   lat.buildBasis(low+1);
   merit = computeMeritProj(lat, norm, red, proj, Coord, true);
   for (int j = low +1; j < upp; j++)
   {
       Coord.insert(j+1);
	   lat.incDim();
     //lat.buildBasis(j+1);
	   merit = computeMeritProj(lat, norm, red, proj, Coord, true);
       if (merit < minmerit) minmerit = merit;
       if (merit < lowbound) {
          //std::cout << "Figure of merit is smaller than lower bound!!!";
          return minmerit;
       }
   }   
   if (upp < dim)
      lat.buildBasis(dim);
   return minmerit;
}

//=========================================================================

template<typename Int>
double FiguresOfMerit<Int>::computeMeritVec(IntLatticeExt<Int, Real> & lat, Normalizer & norm, 
        Reducer<Int, Real> & red, IntLattice<Int, Real> & proj, const IntVec & t) {
	Coordinates Coord;	
	double merit = 0;
	double minmerit = 1.0;
    Int maxDim;
	maxDim = lat.getDim();
	int64_t max_dim, min_dim;
	NTL::conv(max_dim, maxDim);   
	if (projectStationary) {
       min_dim = 2;
	} else min_dim = 1;

	for (int i = 1; i < t.length(); i++) {
       NTL::conv(max_dim, t[i]);
	   CoordinateSets::FromRanges CoordRange(i, i, min_dim, max_dim, projectStationary);  
       for (auto it = CoordRange.begin(); it != CoordRange.end(); it++){
          Coord = *it;
          merit = computeMeritProj(lat, norm, red, proj, Coord, false);
          if (merit < minmerit) minmerit = merit;
	      if (merit < lowbound) {
             //std::cout << "Figure of merit is smaller than lower bound!!!";
             return minmerit;
          }
       }
    }
	return minmerit;  
}

//=========================================================================
	
template class FiguresOfMerit<NTL::ZZ>;
//template class FiguresOfMerit<std::int64_t>;

} // end namespace LatticeTester

#endif
