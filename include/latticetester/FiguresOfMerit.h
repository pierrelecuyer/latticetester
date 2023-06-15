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

//#include <string>
//#include <sstream>
//#include <iostream>
//#include <iomanip>
//#include <vector>
//#include <array>
//#include <set>
//#include <map>
//#include <cmath>
//#include <cstdlib>
//#include <cstdint>
//#include <type_traits>
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
#include "latticetester/FiguresOfMerit.h"
#include "latticetester/Util.h"
#include "latticetester/ParamReader.h"
#include "latticetester/IntLattice.h"
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
	
	//Globally stored parameters of the FoM
	bool performBB = true; //Perform BB algorithm	- IMPLEMENTED
	bool succCoordFirst = false; //At first look at successive corrdinates when calculating FigureOfMerit 
	PreReductionType prered = BKZ; //Define the prereduction type - IMPLEMENTED
	double delta = 0.9; //delta-Parameter for LLL or BKZ - IMPLEMENTED
	int64_t blocksize = 10; //blocksize of BKZ algorithm - IMPLEMENTED
	double lowerbound = 0; //If FoM is below this bound, then calculation of FoM is stopped
	IntMat projBasis; //Matrix for basis of projection
	ProjConstructType pctype = UPPERTRIPROJ; //Type of projection construction
	NTL::ZZ m; //m of the RNG
	WeightsUniform *weights; //Weights for FoM 	
	
	//Set and get functions for global variables
	
    void set_m(Int i) { m = i; };    

	Int get_m() { return m;	};   
  
    void set_performBB(bool BB) { performBB = BB; };
    
    bool get_performBB() { return performBB; };
  
    void set_succCoordFirst(bool first) { succCoordFirst = first; };
    
    bool get_succCoordFirst() { return succCoordFirst; };
  
    void set_PreReductionType(PreReductionType pr) { prered = pr; };
    
    PreReductionType get_PreReductionType () { return prered; };
  
    void set_delta(double del) { delta = del;};
    
    double get_delta() { return delta; };
  
    void set_blocksize(int64_t bl) { blocksize = bl; };
    
    int64_t get_blocksize() { return blocksize; };
  
    void set_lowerbound(double lb) { lowerbound = lb; };
    
    double get_lowerbound() { return lowerbound; };
  
    void set_ProjConstructType (ProjConstructType pt) { pctype = pt; };     
    
    ProjConstructType get_ProjConstructType() { return pctype; };
    
    void set_Weights(double weight) { weights = new WeightsUniform(weight); };
    
    WeightsUniform get_Weights() { return *weights; };
    
    //Further Declarations
    
    //Initialization
    FiguresOfMerit(IntLattice<Int, Real> & lat, Int & i, int64_t max_dim) {
		m = i;   
	    red = new Reducer<Int, Real>(max_dim);   
	    weights = new WeightsUniform(1.0); // This just puts a weight of 1 to everything          
        //This is to initialize proj in order to avoid recreation
        proj = new IntLattice<Int, Real> (lat.getBasis(), m, lat.getBasis().NumCols()); 
	}
	
    //Function for calculating the normalizer
	  void calculNorma(IntLattice<Int, Real> & lat, int64_t & dim) {
     double log_density=(double)(-log(abs(NTL::determinant(lat.getBasis()))));
     norma  = new NormaBestLat(log_density, dim);
    }	

	//Function for calculating the Figure of Merit M
	double computeMeritM(IntLattice<Int, Real> & lat, const IntVec & t);
	
	//Function for calculating the Figure of Merit for a single projection
	double computeMeritProj(IntLattice<Int, Real> & lat, const Coordinates & Coord, const int64_t & max_dim);

	//Function for calculating the Figure of Merit Q
	static void computeMeritQ(const IntLattice<Int, Real> & Lat, const Int & t);	
	
	//Further global objects
	IntLattice<Int, Real> *proj; 
    Normalizer *norma;
	Reducer<Int, Real> *red;  
    CoordinateSets::FromRanges *CoordRange;  	

};

//============================================================================
// Implementation

template<typename Int>
double FiguresOfMerit<Int>::computeMeritM(IntLattice<Int, Real> & lat, const IntVec & t) {
   double merit;
   double minmerit = 1.0;
   Int maxDim;
   maxDim = lat.getDim();
   int64_t max_dim;
   NTL::conv(max_dim, maxDim);
   Coordinates Coord;
   //Do the calculation for the successive coordinates first if succCoordFirst = true
   if (succCoordFirst == true) {
      for (Int i = t[0]; i < maxDim; i++) {
         Coord.clear();
	     for (int j = 0; j < i+1; j++) Coord.insert(j+1);
		 merit = computeMeritProj(lat, Coord, max_dim);
		 if (merit < minmerit) minmerit = merit;
		 if (merit < lowerbound) {
			 //std::cout << "Figure of merit is smaller than lower bound!!!";
			 return minmerit;
		 }
      }
   }   
   //Do the calculation for the other coordinate sets
   for (int i = 1; i < t.length(); i++) {
      CoordinateSets::FromRanges CoordRange(i, i, 1, 5);  
	  for(auto it = CoordRange.begin(); it != CoordRange.end(); it++){
         Coord = *it;
         //According to the notation is the guide, the first coordinate is always included
		 Coord.insert(0);
		 merit = computeMeritProj(lat, Coord, max_dim);
		 if (merit < minmerit) minmerit = merit;
		 if (merit < lowerbound) {
			 //std::cout << "Figure of merit is smaller than lower bound!!!";
			 return minmerit;
		 }
	   }
    }
   //Do the calculation for the successive coordinates last if succCoordFirst = false
    if (succCoordFirst == false) {
      for (Int i = t[0]; i < maxDim; i++) {
         Coord.clear();
	     for (int j = 0; j < i+1; j++) Coord.insert(j+1);
		 merit = computeMeritProj(lat, Coord, max_dim);
		 if (merit < minmerit) minmerit = merit;
		 if (merit < lowerbound) {
			 //std::cout << "Figure of merit is smaller than lower bound!!!";
			 return minmerit;
		 }
      }
   }
 
   return minmerit; 

}

template<typename Int>
double FiguresOfMerit<Int>::computeMeritProj(IntLattice<Int, Real> & lat, const Coordinates & Coord, const int64_t & max_dim) {
   double shortest, merit;
   BasisConstruction<Int>::projectionConstruction(lat.getBasis(), projBasis, Coord, m, pctype);
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
      red->shortestVector(*proj);
      shortest = NTL::conv<double>(red->getMinLength());
      merit = weights->getWeight(Coord) * shortest/norma->getBound((Coord).size());
   } else { merit = 0.0;};
   return merit;

}

//template<typename Int>
//void FoMQ(const IntLattice<Int, Real> & Lat, const Int & t) {
//	std::cout << "Calculating Q";
//}
	
template class FiguresOfMerit<NTL::ZZ>;
//template class FiguresOfMerit<std::int64_t>;

} // end namespace LatticeTester

#endif
