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
#include <ctime>  
#include <cstdint>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>

#include "latticetester/FlexTypes.h"
#include "latticetester/EnumTypes.h"
#include "latticetester/FiguresOfMerit.h"
#include "latticetester/Util.h"
// #include "latticetester/ParamReader.h"
#include "latticetester/IntLattice.h"
#include "latticetester/Rank1Lattice.h"
#include "latticetester/Reducer.h"

// Application specific headers
#include "latticetester/NormaBestLat.h"
#include "latticetester/NormaBestBound.h"
#include "latticetester/CoordinateSets.h"
#include "latticetester/WeightsUniform.h"
#include "latticetester/BasisConstruction.h"

const long numMult = 1;    // Number of multipliers
const long multipliers[numMult] = { 1597};//, 19021, 49109, 71904, 90941, 1090942, 809519, 371915, 1824915, 577841};
const long numDelta = 3; // Number of different deltas
const double deltas[numDelta] = { 0.8, 0.9, 0.99999};
const long numMeth = 5;
clock_t timer[numMeth][numDelta * 2];
int numRep = 1; // Number of repetitions of the code
std::string names[numMeth] = { "A         ", "B         ", "C         ", "D         ",
		"E         "};
const long width = 10;
Int m(1048573);     
clock_t tmp, totalTime; //variables for measuring time elapsed

using namespace LatticeTester;

void printOutput() {
	// Create output string
   std::cout << "Results of FOMExample.cc with m = " << m << "\n";
   std::cout << "Flexible types: " << strFlexTypes << "\n";
   std::cout << "Timings for finding shortest vector with BB, with various pre-reduction methods.\n";
   std::cout << "The times are in basic clock units. \n";
   std::cout << " dim:    ";
   std::cout << "LLL ";
   for (int d = 0; d < numDelta; d++) std::cout << std::setw(width) << deltas[d] << " ";
   std::cout << "      BKZ ";
   for (int d = 0; d < numDelta; d++) std::cout << std::setw(width) << deltas[d] << " ";
	std::cout << std::endl << std::endl;
	for (int meth = 0; meth < numMeth; meth++) {
		std::cout << names[meth] << " ";
		std::cout << "    ";
		for (int d = 0; d < numDelta; d++) {
			std::cout << std::setw(width) << timer[meth][d] << " ";
		}		
		std::cout << "          ";
		for (int d = 0; d < numDelta; d++) {
		std::cout << std::setw(width) << timer[meth][d+numDelta] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
  std::cout << "Total time: " << (double) (clock() - totalTime) / (CLOCKS_PER_SEC) << " seconds\n\n\n";
 
}


int main() {
  // Set all necessary variables
  int64_t max_dim = 32;
  Int a; 
  double f; // Variable for figure of merit
  IntVec t; // t-Vector of the FOM
  long dim; 
  bool with_dual = true; // Decide whether calculation is used or not
  dim = max_dim; 
  Rank1Lattice<Int, Real> *lat; // Initialize variable for lattice 
  Normalizer *norma;
  Weights *weights; 
  weights = new WeightsUniform(1.0); 
  Reducer<Int, Real> *red;
  red = new Reducer<Int, Real>(max_dim);
  
  // Set all necessary objects
  IntLattice<Int, Real> *proj; // The IntLattice used to store projections  
  
  FiguresOfMerit<Int> fom(*weights, *red); 
  
  fom.m_succCoordFirst = true; // successive coordinates shall be calculated first 
  fom.m_reductionMethod = BKZBB; //Set pre-reduction to BKZ
  fom.m_fomInDual = true; // Calculate FoM for the dual
  fom.m_incDualOnly = true; //Only the dual basis shall be calculated when increasing the dimension but not the primal
  fom.m_pctype = LLLPROJ; // Define the projecton type
  
  // Create all objects which need to be passed to the FiguresOfMerit object 
  if (fom.m_fomInDual == true) {
    double log_density=(double)(-log(abs(m)));
    norma = new NormaBestLat(log_density, dim);
    fom.setNormalizer(*norma);
  }
  else {
    Int det;
    det = 1;
    for (int i = 0; i < dim - 1; i ++) det = det * m;
    double log_density=(double)(-log(abs(det)));   // We should not compute a determinant!   *****
    norma = new NormaBestLat(log_density, dim);
    fom.setNormalizer(*norma);
  }  
  

  // Start clock
  totalTime = clock();
  for (int j = 0; j < numRep; j++) {
     for (int i = 0; i < numMult; i++) {
    	 for (int k = 0; k < numDelta; k++) {
           a = multipliers[i]; 
           lat = new Rank1Lattice<Int, Real>(m, a, dim, with_dual);
           lat->buildBasis(dim); 
           proj = new IntLattice<Int, Real> (lat->getBasis(), m, lat->getBasis().NumCols()); // initialize object for projection
           fom.m_delta = deltas[k];
           fom.m_reductionMethod = LLLBB; //Set pre-reduction to LLL

           //FOM M_{32}
           t.SetLength(1); 
           t[0] = 32;
           
           tmp = clock();
           f = fom.computeMeritMSucc_MethodA(*lat, *proj, t);
           timer[0][k] = clock() - tmp;
           
           tmp = clock();
           f = fom.computeMeritMSucc_MethodB(*lat, *proj, t);
           timer[1][k] = clock() - tmp;

           tmp = clock();
           f = fom.computeMeritMSucc_MethodC(*lat, *proj, t);
           timer[2][k] = clock() - tmp;
           
           tmp = clock();
           f = fom.computeMeritMSucc_MethodD(*lat, *proj, t);
           timer[3][k] = clock() - tmp;           

           tmp = clock();
           f = fom.computeMeritMSucc_MethodE(*lat, *proj, t);
           timer[4][k] = clock() - tmp;
           
           fom.m_reductionMethod = BKZBB; //Set pre-reduction to BKZ
           
           tmp = clock();
           f = fom.computeMeritMSucc_MethodA(*lat, *proj, t);
           timer[0][k+numDelta] = clock() - tmp;
           
           tmp = clock();
           f = fom.computeMeritMSucc_MethodB(*lat, *proj, t);
           timer[1][k+numDelta] = clock() - tmp;

           tmp = clock();
           f = fom.computeMeritMSucc_MethodC(*lat, *proj, t);
           timer[2][k+numDelta] = clock() - tmp;
           
           tmp = clock();
           f = fom.computeMeritMSucc_MethodD(*lat, *proj, t);
           timer[3][k+numDelta] = clock() - tmp;           

           tmp = clock();
           f = fom.computeMeritMSucc_MethodE(*lat, *proj, t);
           timer[4][k+numDelta] = clock() - tmp;

           
           //std::cout << "CASE 1: Look at t = " << t << ":" << "\n";
           //std::cout << "Figure of merit M is: " << f << "\n";
           //std::cout << "\n";
        
           //FOM M_{5,32,16,12,8}
           t.SetLength(5);
           t[0] = 5;
           t[1] = 32;
           t[2] = 16;
           t[3] = 12;
           t[4] = 8;
           //f = fom.computeMeritM(*lat, *proj, t);
           //std::cout << "CASE 2: Look at t = " << t << ":" << "\n";
           //std::cout << "Figure of merit M is: " << f << "\n";
           //std::cout << "\n";     
         }
     }
  }
  //std::cout << "Time elapsed: " << (double) (clock() - totalTime) / (CLOCKS_PER_SEC) << " seconds\n";
  printOutput();
  return 0;
}
