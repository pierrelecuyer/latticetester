/**
 * This example shows how to use `LatticeTester` to calculate a figure of merit in practice
 * In a first loop an approximation of the FoM is calculated for many multipliers (numRep)
 * by using the chosen pre-reduction algorithm (meth). A chosen number (noBest) of the 
 * best multipliers according to this loop is stored. Afterwards the exact FoM for these 
 * stored multipliers is calculated by means of the BB algorihm.
 */

#define TYPES_CODE  ZD   // ZZ + double

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


int main() {
  /*
   * This defines the example from your paper with Raymond Couture
   */
  Int m(9223372036854773561); 
  Int a1, a2, a3;
  a1 = 1145902849652723;
  a2 = 0; 
  a3 = -1184153554609676;
  IntVec a;
  a.SetLength(3);
  a[0] = a1;
  a[1] = a2;
  a[2] = a3;
  int64_t dim = 12; // Dimension of the lattice
  
  int64_t max_dim = 32; // Maximal allowed dimension of the lattice
  // double delta = 0.9999; // Delta for the pre-reduction algorithm
  // ReductionType meth = LatticeTester::LLL; // Sets the reduction type
  //The t-vector of the FOM, here M_{16,32,16,12}
  NTL::vector<int64_t> t(3); // length of the t-vector
  t[0] = 12;
   t[1] = 2;
   t[2] = 3;

  /*
   * The following variables are technical and shall not be changed by the user
  */ 
  //double high = 1; // Higher bound when a multiplier is rejected (stays constant)
  //double low = 0; // Lower bound when a multiplier is rejected (dynamically changes in the code)
  double f; // Variable for calculation current figure of merit
  bool with_primal = true; // Shall the primal lattice be calculated?
  bool with_dual = true; // Shall the dual lattice be calculated?
  IntLattice<Int, Real> *proj; // The IntLattice used to store projections  
  MRGLattice<Int, Real> *lat; // Rank1Lattice to store the current lattice for which the FoM is calculated
  Normalizer *norma; // Normalizer object (necessary to normalize FoMs)
  ReducerBB<Int, Real> *red; // Reducer object (necessary for BB)
  WeightsOrderDependent weights; // Object for the weights applied to the FoM
  // WeightsUniform weights; // Object for the weights applied to the FoM

  // Calculate the log-density and initialize the normalizer
  double log_density=(double)(-log(abs(m)));
  norma = new NormaBestLat(log_density, a.length(), max_dim);
  // Initialize the Reducer
  red = new ReducerBB<Int, Real>(max_dim);
  // Set the default weight to 1
  weights.setDefaultWeight(1.0);
  // Initialize the FoM object
  // FigureOfMeritM<Int, Real> fom(t, weights, *norma, *red, true);
  FigureOfMeritDualM<Int, Real> fom(t, weights, *norma, *red, true);

  //Here the MRG Lattices for the lattice and its projections are defined
  lat = new MRGLattice<Int, Real>(m, a, dim, with_primal, with_dual);
  proj = new MRGLattice<Int, Real>(m, a, dim, with_primal, with_dual);
  
  // meth = LLLBB;
  // fom.setReductionMethod(meth, delta);
  fom.setPrintDetails(true);
  f = fom.computeMerit(lat*, proj);
  std::cout << "Figure of merit is: " << f << "\n";
  return 0;
  
}
