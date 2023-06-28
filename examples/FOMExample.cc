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
  int64_t max_dim = 28;
  Int m(2147483647);     // Modulus m = 1021
  Int a; 
  IntVec t;
  double f;
  a = 742938285;
  long dim;
  dim = max_dim;
  bool with_dual = true;
  Normalizer *norma;
  // Build a Korobov lattice
  Rank1Lattice<Int, Real> *lat;
  lat = new Rank1Lattice<Int, Real>(m, a, dim, with_dual);
  // The next four lines are just for test purposes
  lat->buildBasis(3);
  lat->incDimNew();
  //lat->incDim();
  lat->buildBasis(dim);

  FiguresOfMerit<Int> fom(*lat, m, max_dim);
  fom.succCoordFirst = true; // successive coordinates shall be calculated first 
  fom.prered = LLL; //Set pre-reduction to LLL
  fom.forDual = true; // Do calculations for dual
  fom.TypeMerit = MERITM; // Choose type of figure of merit
  
  if (fom.forDual == true) {
		IntMat BasisDual;
		BasisConstruction<Int>::mDualBasis(lat->getBasis(), BasisDual, m);
	    double log_density=(double)(-log(abs(NTL::determinant(BasisDual))));
	    norma  = new NormaBestLat(log_density, dim);
  }
  else {
     double log_density=(double)(-log(abs(NTL::determinant(lat->getBasis()))));
     norma  = new NormaBestLat(log_density, dim);
  }

  t.SetLength(1);
  t[0] = 28;
  //Calculate figure for t = [3 5 8];

  clock_t tmp;
  clock_t timer;
  tmp = clock();
  f = fom.computeMerit(*lat, *norma, t);
  timer = clock() - tmp;
  std::cout << "CASE 1: Look at t = " << t << ":\n";
  std::cout << "Figure of merit M is: " << f << "\n";
  std::cout << "\n";

  t.SetLength(2);
  t[0] = 10;
  t[1] = 8;
  f = fom.computeMerit(*lat, *norma, t);

  std::cout << "CASE 2: Look at t = " << t << ":\n";
  std::cout << "Figure of merit M is: " << f << "\n";
  std::cout << "\n";
  //std::cout << "Figure of merit with BestBound: " << merit2 << std::endl;
  std::cout << "Figures of merit are different for different normalizers,"
               " weights and projections choices\n";
  std::cout << "\n";


  return 0;
}
