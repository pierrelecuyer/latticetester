/**
 * This example shows how to use `LatticeTester` to calculate a figure of merit in practive
 * In a first loop an approximation of the FoM is calculated for many multipliers (numRep)
 * by using the chosen pre-reduction algorithm (meth). A chosen number (noBest) of the 
 * best multipliers according to this loop is stored. Afterwards the exact FoM for these 
 * stored multipliers is calculated by means of the BB algorihm.
 */

//#define NTL_TYPES_CODE 2
#define TYPES_CODE  ZD
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
#include "latticetester/FigureOfMeritM.h"
#include "latticetester/FigureOfMeritDualM.h"
#include "latticetester/Util.h"
#include "latticetester/ParamReader.h"
#include "latticetester/IntLattice.h"
#include "latticetester/Rank1Lattice.h"
#include "latticetester/Reducer.h"
#include "latticetester/LLL_FPInt.h"

// Application specific headers
#include "latticetester/NormaBestLat.h"
#include "latticetester/NormaBestBound.h"
#include "latticetester/CoordinateSets.h"
#include "latticetester/WeightsUniform.h"
#include "latticetester/BasisConstruction.h"

using namespace LatticeTester;


int main() {
  /*
   * The following settings may be changed by the user
   */
  Int m(1048573); // Modulus
  //Int m(1073741827); // Prime modulus near 2^{30}
  //Int m(1099511627791);  // Prime modulus near 2^{40}
  const int noBest = 50; // The number of best multipliers to keep for the BB algorithm
  const int numRep = 1048573; // Total number of multipliers to check 
  int64_t max_dim = 32; // Dimension of the lattice
  double delta = 0.9999; // Delta for the pre-reduction algorithm
  ReductionType meth = LLL; // Sets the reduction type
  //The t-vector of the FOM, here M_{16,32,16,12}
  vector<int64_t> t(4); // length of the t-vector
  t[0] = 16;
  t[1] = 32;
  t[2] = 16;
  t[3] = 12;

  /*
   * The following variables are technical and shall not be changed by the user
   */
  const long numMult = 10; // Number of stored primitive elements
  const long multipliers[numMult] = { 91, 93, 94, 96, 98, 102, 105, 115, 117, 118}; // These are all primitive elements mod 1048573
  double high = 1; // Higher bound when a multiplier is rejected (stays constant)
  double low = 0; // Lower bound when a multiplier is rejected (dynamically changes in the code)
  double f; // Variable for calculation current figure of merit
  double bestFoms[noBest]; // Array to store the best FoMs
  Int bestMultipliers[noBest]; // Array to store the multipliers corresponding to the best FoMs
  clock_t tmp, totalTime; // Variables for measuring time elapsed
  Int a; // Variable to store the current multiplier
  Int g; // Variable to store the greatest common divisor
  bool with_primal = true; // Shall the primal lattice be calculated?
  bool with_dual = true; // Shall the dual lattice be calculated?
  double bestFom; // Variable to store the best FoM after BB
  Int bestMult; // Varialbe to store the multiplier which yields the best FoM
  Rank1Lattice<Int, Real> *lat; // Rank1Lattice to store the current lattice for which the FoM is calculated 
  IntLattice<Int, Real> *proj; // The IntLattice used to store projections  
  Normalizer *norma; // Normalizer object (necessary to normalize FoMs)
  Reducer<Int, Real> *red; // Reducer object (necessary for BB)

  // Initialization of objects
  // Arrays to store information about the FoM must be set to zero
  for (int i = 0; i < noBest; i++) {
      bestFoms[i] = 0;
      bestMultipliers[i] = 0;
  }
  // Calculate the log-density and initilize the normalizer
  double log_density=(double)(-log(abs(m)));
  norma = new NormaBestLat(log_density, max_dim);
  // Initialize the Reducer
  red = new Reducer<Int, Real>(max_dim);
  // Initialize the IntMat objects
  lat = new Rank1Lattice<Int, Real>(m, a, max_dim, with_primal, with_dual);
  proj = new IntLattice<Int, Real> (m, lat->getBasis().NumCols(), with_primal, with_dual);
  // Initialize the FoM object
  FigureOfMeritDualM<Int> fom(t, meth, *red, *norma, true); 
  
  // Set the first multiplier equal to a primitive elements from the list
  a = multipliers[0];
  // Initialize variables for minimum elements of the noBest largest FoMs, get the corresponding index
  // and also get the maximal value.
  auto mini = std::min_element(bestFoms, bestFoms + noBest);
  auto index = std::distance(std::begin(bestFoms), mini);
  auto maxi = std::max_element(bestFoms, bestFoms + noBest);

  // At first find the noBest best multipliers when only applying LLL
  tmp = clock();
  fom.setReductionMethod(meth, delta);
  fom.setPrintDetails(false);
  index = 0;
  for (int64_t j = 0; j < numRep; j++) {
     g = NTL::GCD(NTL::conv<Int>(j+1), m - 1);
     // Only look at the full period lattices
     if (g==1) {
         lat->seta(a);
         lat->buildBasis(max_dim);
         fom.setBounds(low, high);
         f = fom.computeMeritDual(*lat, proj);
         if (f > bestFoms[index]) {
             // Overwrite the current smallest FoM in the stored list
             mini = std::min_element(bestFoms, bestFoms + noBest);
             index = std::distance(std::begin(bestFoms), mini);
             bestFoms[index] = f;
             bestMultipliers[index] = a;
             // Get the newst smallest number in the list
             mini = std::min_element(bestFoms, bestFoms + noBest);
             index = std::distance(std::begin(bestFoms), mini);
             // Set the lower bound to the current smallest FoM of the stored ones
            low = bestFoms[index];
         }
     }
     a = a * multipliers[0] % m;
    	 
  }
  totalTime = clock() - tmp;
  mini = std::min_element(bestFoms, bestFoms + noBest);
  maxi = std::max_element(bestFoms, bestFoms + noBest);
  index = std::distance(std::begin(bestFoms), maxi);
  std::cout << "Time elapsed for LLL: " << (double) totalTime / (CLOCKS_PER_SEC) << " seconds\n";
  std::cout << "Best multiplier found: " << bestMultipliers[index] << "\n";
  std::cout << "Best FoM: " << *maxi << "\n";
  std::cout << noBest <<"-th Best FoM: " << *mini << "\n";
  
  // Now switch the reduction type to LLL + BB, i.e. we perform the BB algorithm too
  lat = new Rank1Lattice<Int, Real>(m, a, max_dim, true, with_dual);
  proj = new IntLattice<Int, Real> (m, lat->getBasis().NumCols(), true, with_dual);
  meth = LLLBB;
  fom.setReductionMethod(meth, delta);
  low = 0;
  fom.setBounds(low, high);
  bestMult = 0;
  bestFom = 0;
  tmp = clock();
  for (int j = 0; j < noBest; j++)
  {
      a = bestMultipliers[j];
      lat->seta(a);
      f = fom.computeMeritDual(*lat, proj);
      if (f > bestFom) {
          bestFom = f;
          bestMult = a; 
      }
  }
  totalTime = clock() - tmp;
  std::cout << "Time elapsed for BB: " << (double) totalTime / (CLOCKS_PER_SEC) << " seconds\n";
  std::cout << "Best multiplier found: " << bestMult << "\n";
  std::cout << "Which has FoM: " << bestFom << "\n";
  
  return 0;
}
