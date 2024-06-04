/**
 * This example shows how to use `LatticeTester` to calculate a figure of merit in practice
 * In a first loop an approximation of the FoM is calculated for many multipliers (numRep)
 * by using the chosen pre-reduction algorithm (meth). A chosen number (noBest) of the 
 * best multipliers according to this loop is stored. Afterwards the exact FoM for these 
 * stored multipliers is calculated by means of the BB algorithm.
 */

#define TYPES_CODE  ZD  // ZZ + double

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
  int64_t dim = 12;  // Maximum dimension of the lattice
  NTL::vector<int64_t> t(1); // The t-vector
  t[0] = 12;  // We look at successive coordinates in up to 12 dimensions.

  bool with_primal = true;  // Needed?
  bool with_dual = true;
  MRGLattice<Int, Real> lat(m, a, dim, with_primal, with_dual);

  WeightsUniform weights(1.0);
  NormaBestLat norma(log(m), order, dim);  // Factors will be computed for primal.
  ReducerBB<Int, Real> red(dim);   // Reducer created for up to dim dimensions.
  // ReducerBB<Int, Real> red;   // To pass no red, it would need to be a pointer.
  FigureOfMeritM<Int, Real> fom(t, weights, norma, red, true);
  fom.setVerbosity(2);
  double merit = fom.computeMeritSucc (lat);
  std::cout << "Figure of merit is: " << merit << "\n\n";

  NormaBestLat normad (-log(m), order, dim);  // Factors will be computed for dual.
  FigureOfMeritDualM<Int, Real> fomd(t, weights, normad, red, true);
  fomd.setVerbosity(2);
  merit = fomd.computeMeritSucc (lat);
  std::cout << "Figure of merit is: " << merit << "\n\n";

  t.SetLength(4);
  t[1] = 0;   t[2] = 0;  t[3] = 8;
  fomd.setTVector(t);
  MRGLattice<Int, Real> proj(m, a, dim, with_primal, with_dual);
  merit = fomd.computeMeritNonSucc (lat, proj);
  std::cout << "Figure of merit is: " << merit << "\n\n";

  return 0;
}
