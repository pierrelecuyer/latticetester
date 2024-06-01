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
  NTL::vector<Int> a;
  a.SetLength(order);
  a[0] = 1145902849652723;
  a[1] = 0;
  a[2] = -1184153554609676;
  int64_t dim = 12; // Dimension of the lattice
  NTL::vector<int64_t> t(1); // The t-vector
  t[0] = 12;

  bool with_primal = true;
  bool with_dual = true;
  WeightsUniform weights(1.0);
  NormaBestLat norma(-log(m), order, dim);
  ReducerBB<Int, Real> red(dim);
  FigureOfMeritDualM<Int, Real> fom(t, weights, norma, red, true);
  fom.setPrintDetails(true);
  MRGLattice<Int, Real> lat(m, a, dim, with_primal, with_dual);

  double merit = fom.computeMeritSuccDual (lat);
  std::cout << "Figure of merit is: " << merit << "\n";
  return 0;
}
