/**
 * This example shows how to use `LatticeTester` to calculate a figure of merit in practice
 * In a first loop an approximation of the FoM is calculated for many multipliers (numRep)
 * by using the chosen pre-reduction algorithm (meth). A chosen number (noBest) of the 
 * best multipliers according to this loop is stored. Afterwards the exact FoM for these 
 * stored multipliers is calculated by means of the BB algorithm.
 */

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
   * This is to illustrate the problems with the primal projections described at the end of section 9 of the guide.
   * Some interesting values of (m, a_1, a_3):  (101, 12, 36),  (13, 4, 8), (13, 12, 2),
   */
  int64_t order(3);
  Int m(13);
  NTL::vector<Int> a(order); // Vector a has size 3.
  a[0] = 4;
  a[1] = 0;
  a[2] = 8;
  int64_t maxdim = 8;  // Maximum dimension of the lattice
  NTL::vector<int64_t> t(4); // The t-vector
  t[0] = 8;  // We look at successive coordinates in up to t[0] dimensions.
  double merit;

  MRGLattice<Int, Real> lat(m, a, maxdim);
  MRGLattice<Int, Real> proj(m, a, 4);

  WeightsUniform weights(1.0);
  NormaBestLat norma(log(m), order, maxdim);  // Normalization factors for primal.
  ReducerBB<Int, Real> red(maxdim);   // Reducer created for up to dim dimensions.
  FigureOfMeritM<Int, Real> fom(t, weights, norma, &red);
  fom.setVerbosity(2);
  // fom.setLLL(0.5);
  // fom.setBKZ(0.0, 0);
  std::cout << "\nFigure of merit primal succ with BB.\n";
  merit = fom.computeMeritSucc (lat);
  std::cout << "FOM value: " << merit << "\n\n";

  // NTL::vector<int64_t> t(4); // Another t-vector
  t.SetLength(4);
  t[0] = 8;  t[1] = 5;   t[2] = 4;  t[3] = 4;
  fom.setTVector(t);

  lat.buildBasis(maxdim);
  std::cout << " Basis B = \n" << lat.getBasis() << "\n";
  Coordinates coord;
  coord.insert(1);
  coord.insert(4);
  std::cout << " Building projection over coord {1,4}. \n";
  lat.buildProjection(proj, coord);
  std::cout << " proj basis (2 x 2) = \n" << proj.getBasis() << "\n";

  std::cout << "\nFigure of merit primal non-succ, with BB. \n";
  merit = fom.computeMeritNonSucc (lat, proj);
  std::cout << "FOM value: " << merit << "\n\n";

  NormaBestLat normadual (-log(m), order, maxdim);  // Factors will be computed for dual.
  FigureOfMeritDualM<Int, Real> fomdual (t, weights, normadual, &red);
  fomdual.setVerbosity(2);

  // lat.buildDualBasis(maxdim);
  std::cout << "\nFigure of merit for dual, succ coordinates.\n";
  merit = fomdual.computeMeritSucc (lat);
  std::cout << "FOM value: " << merit << "\n\n";

  // fomdual.setTVector(t);
  std::cout << "\nFigure of merit dual non-succ, with BB.\n";
  merit = fomdual.computeMeritNonSucc (lat, proj);
  std::cout << "FOM value: " << merit << "\n\n";

  return 0;
}
