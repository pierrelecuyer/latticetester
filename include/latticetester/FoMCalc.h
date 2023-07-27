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

#ifndef LATTICETESTER_FOMCALC_H
#define LATTICETESTER_FOMCALC_H

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
template<typename Int>
class FoMCalc: public FiguresOfMerit<Int>  {

private:
	typedef NTL::vector<Int> IntVec;
	typedef NTL::matrix<Int> IntMat;
public:	
	
    /*
     * Constructor of a FoMCalc object. For this purpose a reducer object 
     * needs to be passed. 
     */
    FoMCalc(Reducer<Int, Real> & red, const vector<int64_t> & t);
	
	
    /*
     * This function calculates the Figure of Merit M using method A: 
     * This method uses incDimBasis to increase the basis at each step, 
     * as explained in the 1997 paper.  
    */
    double computeMeritMSucc_MethodA (IntLatticeExt<Int, Real> & lat);
    
    /*
     * This function calculates the Figure of Merit M using method B:
     * if we do not maintain the dual, incDim will only give you a basis for 
     * the primal. Then you will need to do something else to get the dual basis, 
     * in case your FOM is in the dual. You may compute a triangular basis from 
     * the primal basis, then get the m-dual of it, which will of course be triangular, 
     * then reduce this triangular dual basis via LLL or BKZ, then apply BB.
    */
    double computeMeritMSucc_MethodB (IntLatticeExt<Int, Real> & lat);

    /*
     * This function calculates the Figure of Merit M using method C:
     * With this one, you take the primal basis, then compute the m-dual of 
     * it using the general mDualBasis method, then reduce this m-dual basis 
     * via LLL or BKZ, then apply BB.
    */
    double computeMeritMSucc_MethodC (IntLatticeExt<Int, Real> & lat);
        
    /*
    * This function calculates the Figure of Merit M using method D:
    * The simplest is to just reconstruct the basis for the m-dual lattice from 
    * scratch as in section 4.1.1 of the guide each time we increase the dimension t by 1.  
    * Then we apply LLL or BKZ + BB directly to this m-dual basis. 
    */
    double computeMeritMSucc_MethodD (IntLatticeExt<Int, Real> & lat);
            
    /*
    * This function calculates the Figure of Merit M using method E:
    * That is, we just update the m-dual basis, then apply LLL or BKZ to it, then BB
    */
    double computeMeritMSucc_MethodE (IntLatticeExt<Int, Real> & lat);
    
};

//============================================================================
// Implementation
template<typename Int>
FoMCalc<Int>::FoMCalc(Reducer<Int, Real> & red, const vector<int64_t> & t):
    FiguresOfMerit<Int> (red, t) {};


template<typename Int>
double FoMCalc<Int>::computeMeritMSucc_MethodA(IntLatticeExt<Int, Real> & lat) {
   double merit = 0;
   double minmerit = 1.0;
   int64_t lower_dim = static_cast<int64_t>(this->m_t.size());
   
   lat.buildBasis(lower_dim+1);
   merit = this->computeMeritMNoProj(lat);
   BasisConstruction<Int>::mDualBasis(lat.getDualBasis(), lat.getBasis(), lat.getModulo());  
   if (merit == 0) return merit;
   //uint64_t / size_t
   for (int64_t j = lower_dim+1; j < this->m_t[0]; j++)
   {
       lat.incDimBasis(); // Increase the dimension of the lattice and of its dual 
       merit = this->computeMeritMNoProj(lat);
       BasisConstruction<Int>::mDualBasis(lat.getDualBasis(), lat.getBasis(), lat.getModulo());  
       if (merit < minmerit) minmerit = merit;
       if (merit <= this->m_lowbound || merit > this->m_highbound) return 0;
   }   
   
   return minmerit; 
}

//=========================================================================

template<typename Int>
double FoMCalc<Int>::computeMeritMSucc_MethodB(IntLatticeExt<Int, Real> & lat) {
   double merit = 0;
   double minmerit = 1.0;
   int64_t lower_dim = static_cast<int64_t>(this->m_t.size());
   
   double dim = lat.getBasis().NumCols();
   lat.buildBasis(lower_dim+1);
   BasisConstruction<Int>::mDualBasis(lat.getBasis(), lat.getDualBasis(), lat.getModulo());
   this->m_pctype = UPPERTRIPROJ; // Use upper triangular algorithm to calculate the dual
   merit = this->computeMeritMNoProj(lat);
   if (merit == 0) return merit;
   for (int64_t j = lower_dim+1; j < this->m_t[0]; j++)
   {
       lat.incDimBasis(); // Increase the dimension of the lattice and of its dual 
       BasisConstruction<Int>::mDualBasis(lat.getBasis(), lat.getDualBasis(), lat.getModulo());
       merit = this->computeMeritMNoProj(lat);
       if (merit < minmerit) minmerit = merit;
       if (merit <= this->m_lowbound || merit > this->m_highbound) return 0;
   }   
   if (this->m_t[0] < dim)
      lat.buildBasis(dim);	     

   return minmerit; 
}

//=========================================================================

template<typename Int>
double FoMCalc<Int>::computeMeritMSucc_MethodC(IntLatticeExt<Int, Real> & lat) {
   double merit = 0;
   double minmerit = 1.0;
   int64_t lower_dim = static_cast<int64_t>(this->m_t.size());
   
   double dim = lat.getBasis().NumCols();
   lat.buildBasis(lower_dim+1);
   BasisConstruction<Int>::mDualBasis(lat.getBasis(), lat.getDualBasis(), lat.getModulo());
   this->m_pctype = LLLPROJ; // Use the general method to calculate the dual
   merit = this->computeMeritMNoProj(lat);
   if (merit == 0) return merit;
   for (int64_t j = lower_dim+1; j < this->m_t[0]; j++)
   {
       lat.incDimBasis(); // Increase the dimension of the lattice and of its dual
       BasisConstruction<Int>::mDualBasis(lat.getBasis(), lat.getDualBasis(), lat.getModulo());
       merit = this->computeMeritMNoProj(lat);
       if (merit < minmerit) minmerit = merit;
       if (merit <= this->m_lowbound || merit > this->m_highbound) return 0;
   }   
   if (this->m_t[0] < dim)
      lat.buildBasis(dim);	     

   return minmerit; 
}

//=========================================================================

template<typename Int>
double FoMCalc<Int>::computeMeritMSucc_MethodD(IntLatticeExt<Int, Real> & lat) {
   double merit = 0;
   double minmerit = 1.0;
   int64_t lower_dim = static_cast<int64_t>(this->m_t.size());
   
   double dim = lat.getBasis().NumCols();
   lat.buildDualBasis(lower_dim+1);
   merit = this->computeMeritMNoProj(lat);
   if (merit == 0) return merit;
   for (int64_t j = lower_dim+1; j < this->m_t[0]; j++)
   {
       lat.buildDualBasis(j+1); // Build (only) the dual basis from scratch 
       merit = this->computeMeritMNoProj(lat);
       if (merit < minmerit) minmerit = merit;
       if (merit <= this->m_lowbound || merit > this->m_highbound) return 0;
   }   
   if (this->m_t[0] < dim)
      lat.buildBasis(dim);	     

   return minmerit; 
}

//=========================================================================

template<typename Int>
double FoMCalc<Int>::computeMeritMSucc_MethodE(IntLatticeExt<Int, Real> & lat) {
   double merit = 0;
   double minmerit = 1.0;
   int64_t lower_dim = static_cast<int64_t>(this->m_t.size());
   
   double dim = lat.getBasis().NumCols();
   lat.buildDualBasis(lower_dim+1);
   merit = this->computeMeritMNoProj(lat);
   if (merit == 0) return merit;
   for (int64_t j = lower_dim+1; j < this->m_t[0]; j++)
   {
       lat.incDimDualBasis(); // Increase the dimension of (only) the dual basis 
       merit = this->computeMeritMNoProj(lat);
       if (merit < minmerit) minmerit = merit;
       if (merit <= this->m_lowbound || merit > this->m_highbound) return 0;
   }   
   if (this->m_t[0] < dim)
      lat.buildBasis(dim);    

   return minmerit; 
}

//=========================================================================
	
template class FoMCalc<NTL::ZZ>;
//template class FiguresOfMerit<std::int64_t>;

} // end namespace LatticeTester

#endif
