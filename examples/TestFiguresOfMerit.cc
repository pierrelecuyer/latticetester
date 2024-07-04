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
  int64_t dim = 16; // Maximal dimension of the lattice
  Int a(91); // Primitive element mod 1048573
  const int maxnobest = 1000; // The maximal number of best multipliers to keep
  bool early_discard;
  int numRep = 1000; // Total number of multipliers to check 
  
  Rank1Lattice<Int, Real> lat(m, dim); // Rank1Lattice to store the current lattice for which the FoM is calculated 
  Rank1Lattice<Int, Real> proj(m, dim); // The IntLattice used to store projections  
  
  NTL::vector<int64_t> worstdims(dim); // vector to store the worst dims per lattice
  
  NormaBestLat norma(log(m), 1, dim);  // Normalization factors computed for primal lattice
  NormaBestLat normadual(-log(m), 1, dim);  // Normalization Factors computed for dual lattice
  
  WeightsUniform weights(1.0);
  
  ReducerBB<Int, Real> red(dim);   // Reducer created for up to dim dimensions.

/*
* This function is a test loop over 'numMult' different multipliers which are powers 
* of the primitive element 'a' defined in the list above. It calculates the Figure of Merit
* according to the parameter FigureOfMeritM 'fom' and returns a list of the worst or
* elimination dimensions. The parameter 'fom' needs to be a pointer because it cannot ensured
* otherwise that the figure of merit for the dual is calculated for fom in the respective class.
*/
template<typename Int, typename Real>
static void testLoopDim(FigureOfMeritM<Int, Real> *fom, int numMult) {
   Int g; // Variable to store the greatest common divisor
   double merit = 0;
   double maxmerit = 0;
   Int currMultiplier;
   currMultiplier = a;  
   fom->setCollectLevel(2);
   
   // Empty vector for worst dimensions   
   for (int64_t j = 0; j < dim; j++)
      worstdims[j] = 0;   
  
   // Here we collect the worst dimensions for the primal/dual lattice
   for (int64_t j = 0; j < numMult; j++) {
      lat.seta(currMultiplier);       
      merit = fom->computeMerit (lat, proj);
      maxmerit = max(merit, maxmerit);
      fom->setLowBound(maxmerit);
      worstdims[fom->getMinMeritProj().size()-1]++;    
      currMultiplier = currMultiplier * a % m;
   }   
   
}


/*
* This function calculates the figure of merit 'fom' for all multipliers from the
* list of multipliers 'List'. Its output art the best FoM and the value of the best
* multiplier from the list of length 'noBest'.
*/
template<typename Int, typename Real>
static void testLoopBestFromList(FigureOfMeritM<Int, Real> *fom, Int List[], int noBest) {
   double merit;
   double maxmerit;
   Int winner;
   maxmerit = 0;

   for (int64_t j = 0; j < noBest; j++) {
      lat.seta(List[j]);
      //lat.buildBasis(dim);
      //lat.buildDualBasis(dim);
      merit = fom->computeMerit(lat, proj);   
      if (merit > maxmerit) {
         maxmerit = merit;
         winner = List[j];
      }
   }
   std::cout << "Best FoM: " << maxmerit << "\n";
   std::cout << "Best mulitplier: " << winner << "\n\n";
 
}

/*
* This function loops over 'numMult' many powers of the primitive element a and 
* calculates the figure of merit 'fom' for each of them. The 'noBest' number of
* best multipliers and their Foms are stored in the variables 'bestFoms' and 'bestMultipliers'.
* If 'fom2' is set, then the best of the multipliers is identified by applying
* fom2. A typcially example is that fom only applies the LLL algorithm and fom2 includes the BB algorithm.
*/

template<typename Int, typename Real>
static void testLoopBestFoms(FigureOfMeritM<Int, Real> *fom, int numMult, int noBest, FigureOfMeritM<Int, Real> *fom2 = 0) {
   double merit;
   Int currMultiplier;
   currMultiplier = a; 
    
   
   Int bestMultipliers[noBest]; // Array to store the multipliers corresponding to the best FoMs

   double bestFoms[maxnobest]; // Array to store the best FoMs
   
   // Remove information from prior calculations
   for (int i = 0; i < noBest; i++) {
      bestFoms[i] = 0;
      bestMultipliers[i] = 0;
   }
   // Initialize variables for minimum elements of the noBest largest FoMs, get the corresponding index
   // and also get the maximal value.
   
   auto mini = std::min_element(bestFoms, bestFoms + noBest);
   auto index = std::distance(bestFoms, mini);
   index = 0;
  
   clock_t tmp; // Variables for measuring time elapsed

   tmp = clock();   
   
   for (int64_t j = 0; j < numMult; j++) {
      lat.seta(currMultiplier);
      merit = fom->computeMerit(lat, proj);
         
      if (merit > bestFoms[index]) {
          // Overwrite the currently smallest FoM in the stored list
          mini = std::min_element(bestFoms, bestFoms + noBest);
          index = std::distance(std::begin(bestFoms), mini);
          bestFoms[index] = merit;
          bestMultipliers[index] = currMultiplier;
          // Get the new smallest FoM in the list
          mini = std::min_element(bestFoms, bestFoms + noBest);
          index = std::distance(bestFoms, mini);
          // Set the lower bound to the current smallest FoM of the stored ones  
         if (early_discard)
            fom->setLowBound(*mini);
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
   
   if (fom2)
     testLoopBestFromList<NTL::ZZ, double>(fom2, bestMultipliers, noBest); 
 
}

/*
* This function calls 'testLoopDim' for the primal and the dual lattice. The vector 't'
* used to define the figure of merit and the numbers of multipliers are parameters set by 
* the user (at a later point we might to decide to have more parameters). The BB algorithm
* is always apppied.
*/
template<typename Int, typename Real>
static void testDim(NTL::vector<int64_t> t, int numMult) {

  
  
  FigureOfMeritM<Int, Real> fomprimal(t, weights, norma, &red, true); // FoM for the primal lattice with reducer
  FigureOfMeritDualM<Int, Real> fomdualx(t, weights, normadual, &red, true); // FoM for the dual lattice with reducer
  
 
  std::cout << "##################################" << "\n"; 
  std::cout << "  WORST OR ELIMINATION DIMENSIONS " << "\n"; 
  std::cout << "##################################" << "\n\n";
   
   
   
  // Something goes wrong here - needs to be checked!
  
  std::cout << "RESULTS FOR PRIMAL LATTICE" << "\n\n"; 
  testLoopDim<Int, Real>(&fomprimal, numMult);
  std::cout << worstdims << "\n\n";
  
  std::cout << "RESULTS FOR DUAL LATTICE" << "\n\n";
  testLoopDim<Int, Real>(&fomdualx, numMult);
  std::cout << worstdims << "\n\n";
  
}

/*
* This functions tests the effect of early discaring. If we use early discarding then the calculation of the 
* figures merit is much faster because usually few dimensions suffice to see that the FoM is not the best
* one. The 'noBest' number of best multipliers are stored.
*/
template<typename Int, typename Real>
static void testBestFoms(NTL::vector<int64_t> t, int numMult, int noBest) {
  
  FigureOfMeritM<Int, Real> fomprimal(t, weights, norma, &red, true); // FoM for the primal lattice with reducer
  FigureOfMeritM<Int, Real> fomprimal_wo_red(t, weights, norma, 0, true); //FoM for the primal without reducer
  
  FigureOfMeritDualM<Int, Real> fomdual(t, weights, normadual, &red, true); // FoM for the dual lattice with reducer
  FigureOfMeritDualM<Int, Real> fomdual_wo_red(t, weights, normadual, 0, true); // FoM for the dual lattice without reducer
  
  early_discard = true;
  std::cout << "###########################" << "\n"; 
  std::cout << "  USE OF EARLY DISCARDING  " << "\n";
  std::cout << "###########################" << "\n\n"; 
  
  std::cout << "RESULTS FOR PRIMAL LATTICE" << "\n\n";
  early_discard = true;
  std::cout << "WITH BB and WITH EARLY DISCARDING for " << numMult << " multipliers \n";
  testLoopBestFoms<NTL::ZZ, double>(&fomprimal, numMult, noBest);  
  std::cout << "WITHOUT BB and WITH EARLY DISCARDING for " << numMult << " multipliers \n";
  testLoopBestFoms<NTL::ZZ, double>(&fomprimal_wo_red, numMult, noBest);  
  early_discard = false;
  std::cout << "WITH BB and WITHOUT EARLY DISCARDING for " << numMult << " multipliers \n";
  testLoopBestFoms<NTL::ZZ, double>(&fomprimal, numMult, noBest);  
  std::cout << "WITHOUT BB and WITHOUT EARLY DISCARDING for " << numMult << " multipliers \n";
  testLoopBestFoms<NTL::ZZ, double>(&fomprimal_wo_red, numMult, noBest);  
 
  std::cout << "RESULTS FOR DUAL LATTICE" << "\n\n";
  early_discard = true;
  std::cout << "WITH BB and WITH EARLY DISCARDING for " << numMult << " multipliers \n";
  testLoopBestFoms<NTL::ZZ, double>(&fomdual, numMult, noBest);  
  std::cout << "WITHOUT BB and WITH EARLY DISCARDING for " << numMult << " multipliers \n";
  testLoopBestFoms<NTL::ZZ, double>(&fomdual_wo_red, numMult, noBest);  
  early_discard = false;
  std::cout << "WITH BB and WITHOUT EARLY DISCARDING for " << numMult << " multipliers \n";
  testLoopBestFoms<NTL::ZZ, double>(&fomdual, numMult, noBest);  
  std::cout << "WITHOUT BB and WITHOUT EARLY DISCARDING for " << numMult << " multipliers \n";
  testLoopBestFoms<NTL::ZZ, double>(&fomdual_wo_red, numMult, noBest);

}

/*
* This functtions demonstrates a fast way to search for the mutliplier with the largest FoM.
* At first we calcualte the Figure Of Merit without applying the BB algorithm but only using
* LLL as prereduction. The best 'noBest' FoMs and multipliers are stored. Afterwards the BB algorithm 
* is applied to detect the best multiplier among them.
*/

template<typename Int, typename Real>
static void testSearch(NTL::vector<int64_t> t, int numMult, int noBest) {
  
  FigureOfMeritM<Int, Real> fomprimal(t, weights, norma, &red, true); // FoM for the primal lattice with reducer
  FigureOfMeritM<Int, Real> fomprimal_wo_red(t, weights, norma, 0, true); //FoM for the primal without reducer
  
  FigureOfMeritDualM<Int, Real> fomdual(t, weights, normadual, &red, true); // FoM for the dual lattice with reducer
  FigureOfMeritDualM<Int, Real> fomdual_wo_red(t, weights, normadual, 0, true); // FoM for the dual lattice without reducer  
    
  std::cout << "####################################################" << "\n"; 
  std::cout << "  EFFICIENT SEARCH - FIRST LLL THEN BEST FROM LIST  " << "\n";
  std::cout << "####################################################" << "\n"; 
  early_discard = true;
  std::cout << "RESULTS FOR PRIMAL LATTICE" << "\n\n";
  
  std::cout << "WITHOUT BB and WITH EARLY DISCARDING for " << numMult << " multipliers" << "\n";
  testLoopBestFoms<NTL::ZZ, double>(&fomprimal_wo_red, numMult, noBest, &fomprimal);
  
  std::cout << "RESULTS FOR DUAL LATTICE" << "\n\n";  
  
  std::cout << "WITHOUT BB and WITH EARLY DISCARDING for " << numMult << " multipliers" << "\n";
  testLoopBestFoms<NTL::ZZ, double>(&fomdual_wo_red, numMult, noBest, &fomdual);
  
}



int main() {
  
   NTL::vector<int64_t> t(4); // The t-vector for the FOM.
   t[0] = 16;    // We look at successive coordinates in up to t[0] dimensions.
   t[1] = 16;    // Then pairs and triples up to coord. 5.
   t[2] = 12;
   t[3] = 10;
   
  testDim<NTL::ZZ, double>(t, 500);
  
  testBestFoms<NTL::ZZ, double>(t, 100, 5);
  
  testSearch<NTL::ZZ, double>(t, 1000, 5);  

  return 0;
}
