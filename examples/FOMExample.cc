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

const long numMult = 10;    // Number of multipliers
const long multipliers[numMult] = { 1597, 19021, 49109, 71904, 90941, 1090942, 809519, 371915, 1824915, 577841};
int numRep = 10; // Number of repetitions of the code


using namespace LatticeTester;

int main() {
  // Set all necessary variables
  int64_t max_dim = 32; // Dimension max_im = 32
  Int m(1048573);     // Modulus m = 1048573
  Int a; // Multiplier
  double f; // Variable for figure of merit
  IntVec t; // t-Vector of the FOM
  long dim; // long variable for storing max_dim
  bool with_dual = true; // Decide whether bool variable is used or not
  Normalizer *norma;
  dim = max_dim; // Set dim = max_dim
  Rank1Lattice<Int, Real> *lat; // Initialize variable for lattice 
  clock_t tmp, timer; //variables for measuring time elapsed
  
  a = multipliers[0]; // Set a to initial value
  lat = new Rank1Lattice<Int, Real>(m, a, dim, with_dual); // build the Korobov lattice
  lat->buildBasis(dim); // initialize the basis

  FiguresOfMerit<Int> fom(m); // The FoM calculation currently needs these inputs. This might be changed at a later point in time.

  fom.succCoordFirst = true; // successive coordinates shall be calculated first 
  fom.reductionMethod = LLL; //Set pre-reduction to LLL
  fom.dual = true; // Do calculations for dual
  fom.pctype = UPPERTRIPROJ; // Define the projecton type
  fom.delta = 0.9; // Set delta-value for BKZ or LLL
  fom.ReadOutDual = true; // Set if the dual basis shall be read out directly from the IntLattice object for successive coordinates
  fom.fom = MERITM; // Choose type of figure of merit
  
  // Create all object which need to be passed to the FiguresOfMerit object 
  if (fom.dual == true) {
		IntMat BasisDual;
		BasisConstruction<Int>::mDualBasis(lat->getBasis(), BasisDual, m);
	    double log_density=(double)(-log(abs(NTL::determinant(BasisDual))));
	    norma  = new NormaBestLat(log_density, dim);
  }
  else {
     double log_density=(double)(-log(abs(NTL::determinant(lat->getBasis()))));
     norma  = new NormaBestLat(log_density, dim);
  }  
  Reducer<Int, Real> *red;            
  red = new Reducer<Int, Real>(max_dim);  
  IntLattice<Int, Real> *proj; //The IntLattice used to store projections  
  
  // Start clock
  tmp = clock();
  // Loop over all multipliers
  for (int j = 0; j < numRep; j++) {
     for (int i = 0; i < numMult; i++) {
        a = multipliers[i]; 
        lat = new Rank1Lattice<Int, Real>(m, a, dim, fom.dual); // build the Korobov lattice
        proj = new IntLattice<Int, Real> (lat->getBasis(), m, lat->getBasis().NumCols()); 
        lat->buildBasis(dim); // initialize the basis
  
        //FOM M_{32}
        t.SetLength(1); // Look at the first figure of merit
        t[0] = 32;
        f = fom.computeMeritM(*lat, *norma, *red, *proj, t);
//      std::cout << "CASE 1: Look at t = " << t << ":" << "\n";
//      std::cout << "Figure of merit M is: " << f << "\n";
//      std::cout << "\n";
        
        //FOM M_{5,32,16,12,8}
        t.SetLength(5);
        t[0] = 5;
        t[1] = 32;
        t[2] = 16;
        t[3] = 12;
        t[4] = 8;
        f = fom.computeMerit(*lat, *norma, *red, *proj, t);
//       std::cout << "CASE 2: Look at t = " << t << ":" << "\n";
//       std::cout << "Figure of merit M is: " << f << "\n";
//       std::cout << "\n";     
     }
  }
  timer = clock() - tmp;
  std::cout << "Time elapsed: " << (double) timer / (CLOCKS_PER_SEC) << " seconds\n";
  


  return 0;
}
