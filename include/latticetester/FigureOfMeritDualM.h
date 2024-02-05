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

#ifndef LATTICETESTER_FIGUREOFMERITDUALM_H
#define LATTICETESTER_FIGUREOFMERITDUALM_H

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
template<typename Int>
class FigureOfMeritDualM: public FigureOfMeritM<Int>  {

private:
	typedef NTL::vector<Int> IntVec;
	typedef NTL::matrix<Int> IntMat;
public:	
	
    /*
     * Constructor of a FoMCalc object. For this purpose a reducer object 
     * needs to be passed. 
     */
	FigureOfMeritDualM(const vector<int64_t> & t, ReductionType & meth, Reducer<Int, Real> & red);
	       
    /*
     * Same as computeMeritM for the dual lattice.
     */
    double computeMeritDualM (IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> *proj);
    
    /*
     * Same as computeMeritMSuccPrimal but for the dual lattice.
     */
    double computeMeritMSuccDual (IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> *proj);
    
    /*
     * Same as computeMeritMNonSuccPrimal but for the dual lattice.
     */
    double computeMeritMNonSuccDual (IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> *proj);
};

//============================================================================
// Implementation
template<typename Int>
FigureOfMeritDualM<Int>::FigureOfMeritDualM (const vector<int64_t> & t, ReductionType & meth, Reducer<Int, Real> & red):
    FigureOfMeritM<Int> (t, meth, red) {};

//=========================================================================
template<typename Int>
double FigureOfMeritDualM<Int>::computeMeritDualM (IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> *proj) {
   double merit = 0;
   double minmerit = 1.0;
   
//  if (m_succCoordFirst == true) {    
//	   m_b = new double[lat.getDim()];
//     minmerit = computeMeritMSuccDual_(lat, proj);
//     // In any of these cases the calcuation is stopped
//     if (minmerit == 0 || minmerit < m_lowbound || minmerit > m_highbound) return 0;
//   }   

   this->m_b = new double[this->m_t.size()];
   merit = computeMeritMNonSuccDual(lat, proj);
   if (merit < minmerit) minmerit = merit;
   // In any of these cases the calcuation is stopped
   if (minmerit == 0 || minmerit < this->m_lowbound || minmerit > this->m_highbound) return 0;
   
// if (m_succCoordFirst == false) {    
      this->m_b = new double[lat.getDim()];
      merit = computeMeritMSuccDual(lat, proj);
      if (merit < minmerit) minmerit = merit;
       // In any of these cases the calcuation is stopped
      if (minmerit == 0 || minmerit < this->m_lowbound || minmerit > this->m_highbound) return 0;
// }    
   return minmerit; 
}

//=========================================================================
template<typename Int>
double FigureOfMeritDualM<Int>::computeMeritMSuccDual (IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> *proj) {
   double merit = 0;
   double minmerit = 1.0;
   int64_t lower_dim = static_cast<int64_t>(this->m_t.size());  
   lat.buildDualBasis(lower_dim);
   for (int64_t j = lower_dim + 1; j < this->m_t[0] + 1; j++) {
       lat.incDimDualBasis();
       if (this->m_reductionMethod == BKZBB || this->m_reductionMethod == BKZ) {
          this->m_red->redBKZ(lat.getDualBasis(), this->m_delta, this->m_blocksize);
       } else if (this->m_reductionMethod == LLLBB || this->m_reductionMethod == LLL) {
         LLL_FPZZflex(lat.getDualBasis(), this->m_delta, j, j, this->m_b);
       } else if (this->m_reductionMethod == PAIRBB) {
         this->m_red->redDieter(0);
       }
       lat.dualize();
       if (!this->m_doingBB) {
           lat.updateSingleVecNorm(0,j);
           NTL::conv(merit, sqrt(lat.getVecNorm(0)) / this->m_norma->getBound(j));
       } else {
           if (!this->m_red->shortestVector(lat)) return 0;
           merit = NTL::conv<double>(this->m_red->getMinLength() / this->m_norma->getBound(j));
       }
       lat.dualize();
       if (merit < minmerit) minmerit = merit;
       if (minmerit <= this->m_lowbound || minmerit > this->m_highbound) return 0;
   }
   return minmerit;
   
//   int64_t lower_dim = static_cast<int64_t>(this->m_t.size());
//   lat.buildDualBasisFullMatrix(this->m_t[0]);   
   //lat.buildDualBasisFullMatrix(lower_dim+1);
//   if (this->m_reductionMethod == BKZBB || this->m_reductionMethod == BKZ) {
//      this->m_red->redBKZ(lat.getDualBasis(), this->m_delta, this->m_blocksize);  
//   } else if (this->m_reductionMethod == LLLBB || this->m_reductionMethod == LLL) {
//     LLL_FPZZflex(lat.getDualBasis(), this->m_delta, lower_dim+1, lower_dim+1, this->m_b);
//   } else if (this->m_reductionMethod == PAIRBB) {
//     this->m_red->redDieter(0);
//   }
//   if (!this->m_doingBB) {
//	      NTL::conv(merit, sqrt(this->m_b[0]) / this->m_norma->getBound(lower_dim+1));
//     } else {
//       this->m_projBasis.SetDims(lower_dim+1,lower_dim+1); //m_projBasis resized -> Can this be avoided?
//       for (int k = 0; k < lower_dim+1; k++) {
//          for (int l = 0; l < lower_dim+1; l++) this->m_projBasis[k][l] = lat.getDualBasis()[k][l];
//       }
//       proj.setBasis(this->m_projBasis, this->m_projBasis.NumCols());
//       if (!this->m_red->shortestVector(proj)) return 0;
//       merit = NTL::conv<double>(this->m_red->getMinLength() / this->m_norma->getBound(lower_dim+1));
//     }
//   if (merit == 0) return merit;
//   minmerit = merit;
//   for (int64_t j = lower_dim+2; j < this->m_t[0] + 1; j++)
//   {
	   //lat.incDimDualBasisFullMatrix(j);
//	   if (this->m_reductionMethod == BKZBB || this->m_reductionMethod == BKZ) {
//	       this->m_red->redBKZ(lat.getDualBasis(), this->m_delta, this->m_blocksize);
//	    } else if (this->m_reductionMethod == LLLBB || this->m_reductionMethod == LLL) {
//	    	LLL_FPZZflex(lat.getDualBasis(), this->m_delta, j, j, this->m_b);
//	    } else if (this->m_reductionMethod == PAIRBB) {
//	       this->m_red->redDieter(0);
//	    }
//     if (!this->m_doingBB) {
//	      NTL::conv(merit, sqrt(this->m_b[0]) / this->m_norma->getBound(j));
//     } else {
//       this->m_projBasis.SetDims(j,j); //m_projBasis resized -> Can this be avoided?
//       for (int k = 0; k < j; k++) {
//          for (int l = 0; l < j; l++) this->m_projBasis[k][l] = lat.getDualBasis()[k][l];
//       }
//       proj.setBasis(this->m_projBasis, this->m_projBasis.NumCols());
//       if (!this->m_red->shortestVector(proj)) return 0;
//       merit = NTL::conv<double>(this->m_red->getMinLength() / this->m_norma->getBound(j));
//     }
//       if (merit < minmerit) minmerit = merit;
//       if (minmerit <= this->m_lowbound || minmerit > this->m_highbound) return 0;
//   }  
//   return minmerit; 
}

//=========================================================================
template<typename Int>
double FigureOfMeritDualM<Int>::computeMeritMNonSuccDual (IntLatticeExt<Int, Real> & lat, IntLattice<Int, Real> *proj) {
    double merit = 0;
    double minmerit = 1.0;
    for (auto it = this->m_coordRange.begin(); it != this->m_coordRange.end(); it++){
        lat.buildProjection(proj, *it, this->m_delta);  
      	if (this->m_reductionMethod == BKZBB || this->m_reductionMethod == BKZ) {
      	   this->m_red->redBKZ(proj->getDualBasis(), this->m_delta, this->m_blocksize);
      	} else if (this->m_reductionMethod == LLLBB || this->m_reductionMethod == LLL) {
      	    LLL_FPZZflex(proj->getDualBasis(), this->m_delta, proj->getBasis().NumRows(), proj->getBasis().NumCols(), this->m_b);
      	} else if (this->m_reductionMethod == PAIRBB) {
      	   this->m_red->redDieter(0);
        }
        if (!this->m_doingBB) {
            proj->updateSingleDualVecNorm(0,proj->getDim());
            NTL::conv(merit, sqrt(proj->getDualVecNorm(0)) / this->m_norma->getBound(proj->getDim()));
        } else {
            proj->dualize();
            if (!this->m_red->shortestVector(*proj)) return 0;
            merit = NTL::conv<double>(this->m_red->getMinLength() / this->m_norma->getBound(proj->getDim()));
        }
        if (merit < minmerit) minmerit = merit;
        if (minmerit <= this->m_lowbound || minmerit > this->m_highbound) return 0;
    }
    return minmerit;
//   for (auto it = this->m_coordRange.begin(); it != this->m_coordRange.end(); it++){
//       BasisConstruction<Int>::projectMatrixDual(lat.getBasis(), this->m_projBasis, *it); //m_projBasis is resized by projectMatrixDual
       //QUESTION: Do we potentially need to rows here? Maybe this is missing?
       // Define IntLattice based on mdual basis
//       proj.setBasis(this->m_projBasis, this->m_projBasis.NumCols());  
       // Apply selected reduction technique
//       if (this->m_reductionMethod == BKZBB || this->m_reductionMethod == BKZ) {
//           this->m_red->redBKZ(proj.getBasis(), this->m_delta, this->m_blocksize);  
//       } else if (this->m_reductionMethod == LLLBB || this->m_reductionMethod == LLL) {
//           LLL_FPZZflex(proj.getBasis(), this->m_delta, proj.getBasis().NumRows(), proj.getBasis().NumCols(), this->m_b);    
//       } else if (this->m_reductionMethod == PAIRBB) {
//           this->m_red->redDieter(0);
//       }
//       // TODO: CHECK IF IT IS NECESSARY TO ADD DIMENSIONS TO THE REDUCED BASIS
//       if (!this->m_doingBB) {
//         // Use first basis vector as proxy for shortest length
//	       NTL::conv(merit, sqrt(this->m_b[0]) / this->m_norma->getBound(proj.getDim()));
//       } else {
//         if (!this->m_red->shortestVector(proj)) return 0;
//         merit = NTL::conv<double>(this->m_red->getMinLength() / this->m_norma->getBound(proj.getDim()));
//       }
//       if (merit < minmerit) minmerit = merit;
//       if (merit == 0 || merit < this->m_lowbound) return 0;  
//    }  
//    return minmerit;  
}

//=========================================================================


template class FigureOfMeritDualM<NTL::ZZ>;
//template class FiguresOfMeritM<std::int64_t>;

} // end namespace LatticeTester

#endif
