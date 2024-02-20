/**
 * This example shows how to use `LatticeTester` to calculate a figure of merit that
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


const long numMult = 10;    // Number of multipliers
const long multipliers[numMult] = { 1597, 19021, 49109, 71904, 90941, 1090942, 809519, 371915, 1824915, 577841};
int numRep = 1; // Number of repetitions of the code
const long numDelta = 1; // Number of different deltas
const double deltas[numDelta] = { 0.8 };//, 0.9, 0.99999};
const long numMeth = 6;
clock_t timer[numMeth][2*numDelta];
clock_t tmp, totalTime; //variables for measuring time elapsed
Int m(1048573);
//Int m(1021);
//Int m(1073741827);
std::string names[numMeth] = { "0         ","A         ", "B         ", "C         ", "D         ", "D with lb "};
const long width = 10;
int64_t max_dim = 32;


using namespace LatticeTester;


int main() {
  // Set all necessary variables
  //Int m(1021);
  Int a; 
  double f; // Variable for figure of merit --- was double
  vector<int64_t> t(5); // t-Vector of the FOM
  long dim; 
  bool with_dual = true; // Decide whether calculation is used or not
  dim = max_dim; 
  Rank1Lattice<Int, Real> *lat; // Initialize variable for lattice 
  Normalizer *norma;
  Reducer<Int, Real> *red;
  red = new Reducer<Int, Real>(max_dim);
  ReductionType meth = LLL;
  
  // Set all necessary objects
  IntLattice<Int, Real> *proj; // The IntLattice used to store projections  
  IntMat tempMat; // Lattice used to store projections
  tempMat.SetDims(max_dim, max_dim);  
  
  //FOM M_{5,32,16,12,8}
  t.resize(5);
  t[0] = 5;
  t[1] = 32;
  t[2] = 16;
  t[3] = 12;
  t[4] = 8;
  

  double log_density=(double)(-log(abs(m)));
  norma = new NormaBestLat(log_density, dim);
  
  //FigureOfMeritM<Int> fom(*weights, *red); 
  FigureOfMeritDualM<Int> fom(t, meth, *red, *norma, true); 
  
  a = multipliers[0];
  lat = new Rank1Lattice<Int, Real>(m, a, dim, with_dual);
//  proj = new IntLattice<Int, Real> (lat->getBasis(), m, lat->getBasis().NumCols()); // initialize object for projection
  proj = new IntLattice<Int, Real> (m, lat->getBasis().NumCols(), true, true); // initialize object for projection

  
  for (int j = 0; j < numMult; j++) {
     a = multipliers[j]; 
     for (int k = 0; k < numDelta; k++) {
        lat = new Rank1Lattice<Int, Real>(m, a, dim, true, with_dual);
        lat->buildBasis(dim);         
        f = fom.computeMeritDual(*lat, proj);
        std::cout << "For dim = " << dim << " and a = " << a << " and t = " << t << "" << "\n";
        std::cout << "the figure of merit M is: " << f << "\n";
        std::cout << "\n";  
     }
  }
  return 0;
}
