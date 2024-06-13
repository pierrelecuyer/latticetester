/**
 * This example shows how to use `LatticeTester` to calculate a figure of merit in practice
 * In a first loop an approximation of the FoM is calculated for many multipliers (numRep)
 * by using the chosen pre-reduction algorithm (meth). A chosen number (noBest) of the 
 * best multipliers according to this loop is stored. Afterwards the exact FoM for these 
 * stored multipliers is calculated by means of the BB algorithm.
 */

#define TYPES_CODE  ZQ  // ZZ + quad_float

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
   * This tests the MRG retained in Table VII of the LatMRG paper of L'Ecuyer and Couture (1997).
   */
  int64_t order(3);
  Int m(9223372036854773561); 
  NTL::vector<Int> a(order); // Vector a has size 3.
  a[0] = 1145902849652723;
  a[1] = 0;
  a[2] = -1184153554609676;
  int64_t maxdim = 8;  // Maximum dimension of the lattice
  NTL::vector<int64_t> t0(1); // The t-vector
  t0[0] = 8;  // We look at successive coordinates in up to t[0] dimensions.
  double merit;

  MRGLattice<Int, Real> lat(m, a, maxdim);
  MRGLattice<Int, Real> proj(m, a, 4);

  WeightsUniform weights(1.0);
  NormaBestLat norma(log(m), order, maxdim);  // Factors will be computed for primal.
  ReducerBB<Int, Real> red(maxdim);   // Reducer created for up to dim dimensions.
  FigureOfMeritM<Int, Real> fom(t0, weights, norma, &red);
  fom.setVerbosity(2);
  merit = fom.computeMeritSucc (lat);
  std::cout << "Figure of merit primal succ, with BB: " << merit << "\n\n";

  NTL::vector<int64_t> t(4); // Another t-vector
  // t.SetLength(4);
  t[0] = 8;  t[1] = 4;   t[2] = 4;  t[3] = 4;
  fom.setTVector(t);

  lat.buildBasis(maxdim);
  merit = fom.computeMeritNonSucc (lat, proj);
  std::cout << "Figure of merit primal succ + non-succ, with BB: " << merit << "\n\n";

  // Now the m-dual
  NormaBestLat normadual (-log(m), order, maxdim);  // Factors will be computed for m-dual.
  FigureOfMeritDualM<Int, Real> fomdual (t0, weights, normadual, &red);
  fomdual.setLLL(0.8);
  fomdual.setBKZ(0.9999999, 12);
  fomdual.setVerbosity(2);
  // lat.buildDualBasis(maxdim);
  merit = fomdual.computeMeritSucc (lat);
  std::cout << "Figure of merit for dual: " << merit << "\n\n";

  fomdual.setTVector(t);
  merit = fomdual.computeMeritNonSucc (lat, proj);
  std::cout << "Figure of merit dual non-succ, with BB: " << merit << "\n\n";

  return 0;
}
