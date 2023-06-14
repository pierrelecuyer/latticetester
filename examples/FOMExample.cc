/**
 * This example shows how to use `LatticeTester` to implement a figure of merit that
 * is a weighted sum, maximum, minimum or average of a measure on projections on
 * a lattice. It shows how to use the `Weights` classes, how to build
 * projections of a basis, and how to normalize a computation. Typically,
 * when building figures of merit, measures need to be rescaled to the same
 * interval to be compared with one another, this is what is called
 * normalization.
 * 
 * This program computes a simple spectral test a) on all projections of a lattice
 * in 10 dimensions and b) on all two- and three-dimensional projections, 
 * normalizes it between 0 and 1 and then takes the minimal value observed as a figure 
 * of merit for that lattice. 
 * */

//#define NTL_TYPES_CODE 2
#define TYPES_CODE  ZD
#include <iostream>
#include <cstdint>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>

#include "latticetester/FlexTypes.h"
#include "latticetester/EnumTypes.h"
#include "latticetester/FiguresOfMerit.h"
#include "latticetester/Util.h"
#include "latticetester/ParamReader.h"
#include "latticetester/IntLattice.h"
#include "latticetester/Rank1Lattice.h"
#include "latticetester/Reducer.h"

// Application specific headers
#include "latticetester/NormaBestLat.h"
#include "latticetester/NormaBestBound.h"
#include "latticetester/CoordinateSets.h"
#include "latticetester/WeightsUniform.h"
#include "latticetester/BasisConstruction.h"

using namespace LatticeTester;

int main() {
  // Reading a matrix to use the normalizer class
  int64_t max_dim = 10;
  Int m(1021);     // Modulus m = 1021
  Int a; 
  IntVec t;
  double f;
  a = 57;
  long dim;
  dim = max_dim;
  //Build a Korobov lattice
  Rank1Lattice<Int, Real> *lat;
  lat = new Rank1Lattice<Int, Real>(m, a, dim);
  lat->buildBasis(dim);   
  
  FiguresOfMerit<Int> fom(m, max_dim);
  fom.calculNorma(*lat, max_dim);
  //fom.set_ProjConstructType(LLLPROJ);
  fom.set_PreReductionType(LLL);
  //fom.set_lowerbound(0.5);

  t.SetLength(3);
  t[0] = 3;
  t[1] = 5;
  t[2] = 8;
  //Calculate figure for t = [3 5 8];
  f = fom.computeMeritM(*lat, t);
  std::cout << "CASE 1: Look at t = [3,5,8]:" << "\n";
  std::cout << "Figure of merit FoM M is: " << f << "\n";
  std::cout << "\n";
  t.SetLength(2);
  t[0] = 2;
  t[1] = 4;  
  std::cout << "CASE 2: Look at t = [2,4]:" << "\n";
  std::cout << "Figure of merit with BestLat: " << f << "\n";
  std::cout << "\n";
  //std::cout << "Figure of merit with BestBound: " << merit2 << std::endl;
  std::cout << "Figures of merit are different for different normalizers,"
               " weights and projections choices\n";
  std::cout << "\n";


  return 0;
}
