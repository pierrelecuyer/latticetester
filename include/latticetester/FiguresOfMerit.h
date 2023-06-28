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
	 * Variable decides BB algorithm is performed
	 */
	bool performBB = true; 
	
	/*
	 * If set to true successive corrdinates a considered first when calculating FigureOfMerit
	 */
	bool succCoordFirst = false;
	
	/*
	 * Figure of merit is calculated for the dual lattice if forDual is true
	 */
	bool forDual = false; 
	
	/*
	 * Sets if first variable shall always be included for the non-successive coordinates
	 */
	bool projectStationary = true;
	
	/* 
	 * Defines the prereduction type
	 */
	PreReductionType prered = BKZ; 
	
	/* 
	 * Defines the type of the figure of merit
	 */	
	MeritType TypeMerit = MERITM; 
	
	/*
	 * delta-Parameter for LLL or BKZ reduction
	 */
	double delta = 0.9; 
	
	/*
	 * blocksize of BKZ algorithm
	 */
	int64_t blocksize = 10; 
	
	/*
	 * If FoM is below this bound, then calculation of FoM is stopped
	 */	
	double lowerbound = 0; // 
	
	/*
	 * Type of projection construction
	 */	
	ProjConstructType pctype = UPPERTRIPROJ;
	
	/*
	 * Variable used to store the m of the RNG
	 */
	NTL::ZZ m; 
	
	/*
	 * Variable containing the Weights for the FoM 
	 */
	WeightsUniform *weights; 
	        
    /*
     * Initialization of a FiguresOfMerit object
     */
    FiguresOfMerit(IntLatticeExt<Int, Real> & lat, Int & i, int64_t max_dim) {
		m = i;   
	    red = new Reducer<Int, Real>(max_dim);   
	    weights = new WeightsUniform(1.0); // This just puts a weight of 1 to everything          
        //This is to initialize proj in order to avoid recreation
        proj = new IntLattice<Int, Real> (lat.getBasis(), m, lat.getBasis().NumCols()); 
	}
	
    /*
     * This function calculates the normalizer based on whether
     * FiguresOfMerit object is used for primal or dual lattice
     */
	void calculNorma(IntLatticeExt<Int, Real> & lat, int64_t & dim) {
       if (forDual == true) {
    		IntMat BasisDual;
    		BasisConstruction<Int>::mDualBasis(lat.getBasis(), BasisDual, m);
    	    double log_density=(double)(-log(abs(NTL::determinant(BasisDual))));
    	    norma  = new NormaBestLat(log_density, dim);
       }
       else {
          double log_density=(double)(-log(abs(NTL::determinant(lat.getBasis()))));
          norma  = new NormaBestLat(log_density, dim);
       }
    }	
	
	/*
	 * This function calculates the Figure of Merit M or Q based 
	 * on the chosen MeritType of a given lattice 'lat'. The vector 't'
	 * which defes the set of dimensions for which the figure of merit
	 * is calculated needs to be passed as second input variable.  
	 */
	double computeMerit(IntLatticeExt<Int, Real> & lat, const IntVec & t);
		
	/*
	 * This function calculates the Figure of Merit for a single projection
	 * of a given lattice 'lat'. The variable 'Coord' sets the coordinates
	 * of the projection to use
	 */
	double computeMeritProj(IntLatticeExt<Int, Real> & lat, const Coordinates & Coord);

	/*
	 * This functions calculates the Figure of Merit for all projections 
	 * consisting of successive coordinates of the forms 
	 * {1, 2, ..., low} to {1, 2, ..., upp}
	 * Note that upp > low must hold 
	 */
	double computeMeritSucc(IntLatticeExt<Int, Real> & lat, const int64_t & low, const Int & upp);
	
	/*
	 * This functions calculates the Figure of Merit for all projections 
	 * consisting of non-successive coordinates.
	 */
	double computeMeritVec(IntLatticeExt<Int, Real> & lat, const IntVec & t);
	/*
	 * Stores the matrix of the projection basis
	 */
	IntMat projBasis; 
	
	/*
	 * IntLattice object to store the projection matrix
	 */
	IntLattice<Int, Real> *proj; 
	
	/*
	 * A normalizer object used to normalize figure of merit
	 */
    Normalizer *norma;
    
    /*
     * A reducer object used to perform reductions
     */
	Reducer<Int, Real> *red;  
	
	/*
	 * CoordinateSets object is used to store sets of coordinates
	 */
    CoordinateSets::FromRanges *CoordRange;  	

};

//============================================================================
// Implementation

template<typename Int>
double FiguresOfMerit<Int>::computeMerit(IntLatticeExt<Int, Real> & lat, const IntVec & t) {
   double merit = 0;
   double minmerit = 1.0;
   Int maxDim;
   maxDim = lat.getDim();
   int64_t max_dim, min_dim, low_dim;
   NTL::conv(max_dim, maxDim);
   Coordinates Coord;
   
   // Do the calculation for the successive coordinates first if succCoordFirst = true
   if (succCoordFirst == true) { 
	   low_dim = t.length();
	   minmerit = computeMeritSucc(lat, low_dim, t[0]);
	   if (minmerit < lowerbound) return minmerit;
   }   
   
   //Do the calculation for the other coordinate sets
   merit = computeMeritVec(lat, t);
   if (merit < minmerit) minmerit = merit;
   if (minmerit < lowerbound) return minmerit;

   // Do the calculation for the successive coordinates last if succCoordFirst = false
   if (succCoordFirst == false) {
	   low_dim = t.length();
	   minmerit = computeMeritSucc(lat, low_dim, t[0]);
	   if (minmerit < lowerbound) return minmerit;
   }
 
   return minmerit; 

}


template<typename Int>
double FiguresOfMerit<Int>::computeMeritProj(IntLatticeExt<Int, Real> & lat, const Coordinates & Coord) {
   double shortest, merit;
   merit = 0.0;
   
   BasisConstruction<Int>::projectionConstruction(lat.getBasis(), projBasis, Coord, m, pctype);
   if (forDual == true)
   { 
      IntMat projBasisDual;
	  if (pctype == UPPERTRIPROJ) {
	     BasisConstruction<Int>::mDualUpperTriangular(projBasis, projBasisDual, m);
	  } else
	     BasisConstruction<Int>::mDualBasis(projBasis, projBasisDual, m);
	  projBasis = projBasisDual;
   }
   *proj = IntLattice<Int, Real> (projBasis, m, projBasis.NumCols());
   //double log_density=(double)(-log(abs(NTL::determinant(proj->getBasis())))); 
   //Normalizer* norma = new NormaBestLat(log_density, max_dim);
   proj->updateVecNorm();
   proj->sort(0); 
   if (prered == BKZ) {
       red->redBKZ(proj->getBasis(), delta, blocksize);  
   } else {
      red->redLLLNTL(proj->getBasis(), delta);  
   } 
   if (performBB == true) {
	  if (TypeMerit == MERITQ) {
         red->reductMinkowski(*proj, 0);
         merit = NTL::conv<double>(red->getMinLength()) / NTL::conv<double>(red->getMaxLength());
      } else {
      red->shortestVector(*proj);
      shortest = NTL::conv<double>(red->getMinLength());
      // std::cout << norma->getBound((Coord).size()) << "\n";
      merit = weights->getWeight(Coord) * shortest/norma->getBound((Coord).size());
      }
   } else { merit = 0.0;};
   return merit;

}

template<typename Int>
double FiguresOfMerit<Int>::computeMeritSucc(IntLatticeExt<Int, Real> & lat, const int64_t & low, const Int & upp) {
   Coordinates Coord;	
   double merit = 0;
   double minmerit = 1.0;
   Coord.clear();
   for (int j = 0; j < low + 1; j++) Coord.insert(j+1);
   lat.buildBasis(low+1);
   merit = computeMeritProj(lat, Coord);
   for (int j = low +1; j < upp; j++)
   {
       Coord.insert(j+1);
	   lat.incDim();
	   merit = computeMeritProj(lat, Coord);
       if (merit < minmerit) minmerit = merit;
       if (merit < lowerbound) {
          //std::cout << "Figure of merit is smaller than lower bound!!!";
          return minmerit;
       }
   }   
   return minmerit;
}

template<typename Int>
double FiguresOfMerit<Int>::computeMeritVec(IntLatticeExt<Int, Real> & lat, const IntVec & t) {
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
          merit = computeMeritProj(lat, Coord);
          if (merit < minmerit) minmerit = merit;
	      if (merit < lowerbound) {
             //std::cout << "Figure of merit is smaller than lower bound!!!";
             return minmerit;
          }
       }
    }
	return minmerit;  
}
	
template class FiguresOfMerit<NTL::ZZ>;
//template class FiguresOfMerit<std::int64_t>;

} // end namespace LatticeTester

#endif
