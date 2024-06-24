/**
 * This example shows ...
 */

#define TYPES_CODE  ZD  // ZZ + double

#include <iterator>
#include <iostream>
#include <cstdint>
#include <algorithm>
#include <vector>
#include <numeric>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>

#include "latticetester/FlexTypes.h"
#include "latticetester/EnumTypes.h"
#include "latticetester/Util.h"
#include "latticetester/IntLattice.h"
#include "latticetester/Rank1Lattice.h"
#include "latticetester/FigureOfMeritM.h"
#include "latticetester/FigureOfMeritDualM.h"
#include "latticetester/ReducerStatic.h"
#include "latticetester/LLL_lt.h"
#include "latticetester/Weights.h"
#include "latticetester/WeightsOrderDependent.h"

// Application specific headers
#include "latticetester/NormaBestLat.h"
#include "latticetester/NormaBestBound.h"
#include "latticetester/CoordinateSets.h"
#include "latticetester/WeightsUniform.h"
#include "latticetester/BasisConstruction.h"
#include "latticetester/MRGLattice.h"

using namespace LatticeTester;

  //These variables are global at the moment but some of them may be turned local at some point  
  Int m(1048573); // Modulus
  int64_t dim = 32; // Dimension of the lattice
  Int a(91); // Generating multiplier for the given modulus
  const int noBest = 5; // The number of best multipliers to keep
  double bestFoms[noBest]; // Array to store the best FoMs
  Int bestMultipliers[noBest]; // Array to store the multipliers corresponding to the best FoMs
  bool early_discard;
  int numRep = 1000; // Total number of multipliers to check 
  
  Rank1Lattice<Int, Real> lat(m, a, dim); // Rank1Lattice to store the current lattice for which the FoM is calculated 
  Rank1Lattice<Int, Real> proj(m, a, dim); // The IntLattice used to store projections  
  
  NTL::vector<int64_t> worstdims(dim); // vector to store the worst dims per lattice
   

template<typename Int, typename Real>
static void testLoopDim(FigureOfMeritM<Int, Real> fom) {
   Int g; // Variable to store the greatest common divisor
   double merit;
   double maxmerit = 0;
   Int currMultiplier;
   currMultiplier = a;      
  
   fom.setCollectLevel(2);
   // Empty vector for worst dimensions
   for (int64_t j = 0; j < dim; j++)
      worstdims[j] = 0;   
  
   // 1. Here we collect the worst dimensions for the primal/dual lattice
   for (int64_t j = 0; j < numRep; j++) {
      g = NTL::GCD(NTL::conv<Int>(j+1), m - 1);
      // Only look at the full period lattices
      if (g==1) {
         lat.seta(currMultiplier);       
         merit = fom.computeMerit (lat, proj);
         maxmerit = max(merit, maxmerit);
         fom.setLowBound(maxmerit);
         worstdims[fom.m_worstproj.size()-1] += 1;    
      }
      currMultiplier = currMultiplier * a % m;
   }
}

template<typename Int, typename Real>
static void testLoopBestFoms(FigureOfMeritM<Int, Real> fom) {
   Int g; // Variable to store the greatest common divisor
   double merit;
   Int currMultiplier;
   currMultiplier = a;
  
   
   // Remove information from prior calculations
   for (int i = 0; i < noBest; i++) {
      bestFoms[i] = 0;
      bestMultipliers[i] = 0;
   }
   // Initialize variables for minimum elements of the noBest largest FoMs, get the corresponding index
   // and also get the maximal value.
   auto mini = std::min_element(bestFoms, bestFoms + noBest);
   auto index = std::distance(std::begin(bestFoms), mini);
   index = 0;
  
   clock_t tmp; // Variables for measuring time elapsed

   tmp = clock();   
   for (int64_t j = 0; j < numRep; j++) {
     g = NTL::GCD(NTL::conv<Int>(j+1), m - 1);
     // Only look at the full period lattices
     if (g==1) {
         lat.seta(currMultiplier);
         merit = fom.computeMerit(lat, proj);
         
         if (merit > bestFoms[index]) {
             // Overwrite the currently smallest FoM in the stored list
             mini = std::min_element(bestFoms, bestFoms + noBest);
             index = std::distance(std::begin(bestFoms), mini);
             bestFoms[index] = merit;
             bestMultipliers[index] = currMultiplier;
             // Get the new smallest FoM in the list
             mini = std::min_element(bestFoms, bestFoms + noBest);
             index = std::distance(std::begin(bestFoms), mini);
             // Set the lower bound to the current smallest FoM of the stored ones  
            if (early_discard)
               fom.setLowBound(*mini);
         }
     }
     currMultiplier = currMultiplier * a % m;;
   }
   tmp = clock() - tmp;
   std::cout << "Best " << noBest << " multipliers found: ";
   for (int64_t j = 0; j < noBest; j++)
     std::cout << bestMultipliers[j] << ", ";
   std::cout << "\n";
   std::cout << "Corresponding FoMs: ";
   for (int64_t j = 0; j < noBest; j++)
     std::cout << bestFoms[j] << ", ";
   std::cout << "\n";
   std::cout << "Running Time: " << (double) tmp / (CLOCKS_PER_SEC) << "s \n\n";
 
}

template<typename Int, typename Real>
static void testLoopBestFromList(FigureOfMeritM<Int, Real> fom, Int List[]) {
   double merit;
   double maxmerit;
   Int winner;
   maxmerit = 0;

   for (int64_t j = 0; j < noBest; j++) {
      lat.seta(List[j]);
      lat.buildBasis(dim);
      merit = fom.computeMerit(lat, proj);   
      if (merit > maxmerit) {
         maxmerit = merit;
         winner = List[j];
      }
   }
   std::cout << "Best FoM: " << merit << "\n";
   std::cout << "Best mulitplier: " << winner << "\n\n";
 
}



int main() {

  
  WeightsUniform weights(1.0);
  NormaBestLat norma(log(m), 1, dim);  // Factors will be computed for primal.
  NormaBestLat normadual(-log(m), 1, dim);  // Factors will be computed for dual
  
  ReducerBB<Int, Real> red(dim);   // Reducer created for up to dim dimensions.

  NTL::vector<int64_t> t(4); // length of the t-vector
  t[0] = 16;
  t[1] = 32;
  t[2] = 16;
  t[3] = 12;  
  
  FigureOfMeritM<Int, Real> fom(t, weights, norma, &red); // FoM for the primal lattice with reducer
  FigureOfMeritM<Int, Real> fom_wo_red(t, weights, norma, 0); //FoM for the primal without reducer
  
  FigureOfMeritDualM<Int, Real> fomdual(t, weights, normadual, &red); // FoM for the dual lattice with reducer
  FigureOfMeritDualM<Int, Real> fomdual_wo_red(t, weights, normadual, 0); // FoM for the dual lattice without reducer
  
    
  numRep = 1048573;
  std::cout << "##################################" << "\n"; 
  std::cout << "1. WORST OR ELIMINATION DIMENSIONS" << "\n"; 
  std::cout << "##################################" << "\n\n";
  
  std::cout << "RESULTS FOR PRIMAL LATTICE" << "\n\n"; 
  testLoopDim<NTL::ZZ, double>(fom);
  std::cout << worstdims << "\n\n";
  
  std::cout << "RESULTS FOR DUAL LATTICE" << "\n\n";
  testLoopDim<NTL::ZZ, double>(fomdual);
  std::cout << worstdims << "\n\n";
  
 
  numRep = 1000;
  early_discard = true;
  std::cout << "##########################" << "\n"; 
  std::cout << "2. USE OF EARLY DISCARDING" << "\n";
  std::cout << "##########################" << "\n\n"; 
  
  std::cout << "RESULTS FOR PRIMAL LATTICE" << "\n\n";
  early_discard = true;
  std::cout << "WITH BB and WITH EARLY DISCARDING for " << numRep << " multipliers \n";
  testLoopBestFoms<NTL::ZZ, double>(fom);  
  std::cout << "WITHOUT BB and WITH EARLY DISCARDING for " << numRep << " multipliers \n";
  testLoopBestFoms<NTL::ZZ, double>(fom_wo_red);  
  early_discard = false;
  std::cout << "WITH BB and WITHOUT EARLY DISCARDING for " << numRep << " multipliers \n";
  testLoopBestFoms<NTL::ZZ, double>(fom);  
  std::cout << "WITHOUT BB and WITHOUT EARLY DISCARDING for " << numRep << " multipliers \n";
  testLoopBestFoms<NTL::ZZ, double>(fom_wo_red);
  
 
  std::cout << "RESULTS FOR DUAL LATTICE" << "\n\n";
  early_discard = true;
  std::cout << "WITH BB and WITH EARLY DISCARDING for " << numRep << " multipliers \n";
  testLoopBestFoms<NTL::ZZ, double>(fomdual);  
  std::cout << "WITHOUT BB and WITH EARLY DISCARDING for " << numRep << " multipliers \n";
  testLoopBestFoms<NTL::ZZ, double>(fomdual_wo_red);  
  early_discard = false;
  std::cout << "WITH BB and WITHOUT EARLY DISCARDING for " << numRep << " multipliers \n";
  testLoopBestFoms<NTL::ZZ, double>(fomdual);  
  std::cout << "WITHOUT BB and WITHOUT EARLY DISCARDING for " << numRep << " multipliers \n";
  testLoopBestFoms<NTL::ZZ, double>(fomdual_wo_red);
  
  
  std::cout << "###################################################" << "\n"; 
  std::cout << "3. EFFICIENT SEARCH - FIRST LLL THEN BEST FROM LIST" << "\n";
  std::cout << "###################################################" << "\n"; 
  early_discard = true;
  numRep = 1048573;
  
  std::cout << "RESULTS FOR PRIMAL LATTICE" << "\n\n";
  
  std::cout << "WITHOUT BB and WITH EARLY DISCARDING for all multipliers" << "\n";
  testLoopBestFoms<NTL::ZZ, double>(fom_wo_red);
  
  testLoopBestFromList<NTL::ZZ, double>(fom, bestMultipliers); 
  
  std::cout << "RESULTS FOR DUAL LATTICE" << "\n\n";  
  
  std::cout << "WITHOUT BB and WITH EARLY DISCARDING for all multipliers" << "\n";
  testLoopBestFoms<NTL::ZZ, double>(fomdual_wo_red);
  
  testLoopBestFromList<NTL::ZZ, double>(fomdual, bestMultipliers); 

  return 0;
}
