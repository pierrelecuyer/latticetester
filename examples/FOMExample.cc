/**
 * In this example, we compare different ways of computing the shortest vector
 * in a lattice or in its $m$-dual.  For this, we do the BB after some type of
 * pre-reduction of the basis. We compare the total times to do that
 * (pre-red. + BB) for various types of pre-reductions, with lattices that come
 * from Korobov lattice rules in 10 to 40 dimensions, with prime modulus `m`.
 * For the BB, we use the Cholesky decomposition.
 * The pre-reductions considered are:
 * none, pairwise, LLL with delta = 0.5, 0.8, and 0.99999, and BKZ with the
 * default parameters. We do this for two different prime values of `m`.
 * The timings are in terms of total number of clock ticks used by each method.
 *
 * Example of results:

Flexible types: Int = NTL::ZZ, Real = double

Timings for finding shortest vector with BB, with various pre-reduction methods.
The times are in basic clock units.

m = 1048573
dim \ pre-reduct:    None       pairwise        LLL5       LLL8    LLL99999     BKZ

10                    68550         2613         163        219         162     306
20                   532359       206594        1903       1163        1764    2259
30                     --            --       220817      48455       22932   22194
40                     --            --      1333125      83891       59178   43498
Total time: 2.66798 seconds

m = 1073741827
dim \ pre-reduct:    None       pairwise        LLL5     LLL8     LLL99999      BKZ

10                  2279295       103743         260        238         518     401
20                  5937701      1934427        2523       1331        1595    3100
30                        0                   584351      90625       16362   18358
40                        0                239600527   15685520     2797533 3390653
Total time: 272.465 seconds

**/

//#define TYPES_CODE  LD     // Int == int64_t
#define TYPES_CODE  ZD     // Int == ZZ

#include <iostream>
#include <ctime>  
#include <NTL/mat_GF2.h>
#include <cstdint>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include "latticetester/FlexTypes.h"    // This defines Int = int64_t
#include "latticetester/EnumTypes.h"
#include "latticetester/BasisConstruction.h"
#include "latticetester/Util.h"
#include "latticetester/ParamReader.h"
#include "latticetester/IntLattice.h"
#include "latticetester/Rank1Lattice.h"
#include "latticetester/Reducer.h"
#include "latticetester/FiguresOfMerit.h"

using namespace LatticeTester;

Int m, a;     // Modulus, LCG Multiplier
// Int m(1048573);  // Prime modulus near 2^{20}
// Int m(1073741827);  // Prime modulus near 2^{30}
// Int m(1099511627791);  // Prime modulus near 2^{40}
// Int m(1125899906842597);  // Prime modulus near 2^{50}

const long numMeth = 6;
const long numSizes = 4;    // Number of matrix sizes (choices of dimension).
const long numPrimes = 2;    // Number of matrix sizes (choices of dimension).
const long dimensions[numSizes] = { 10, 20, 30, 40 };
const long colwidth[numMeth] = { 25, 12, 11, 10, 11, 11 };
const long primes[numPrimes] = { 1048573, 1073741827 };
const DecompTypeBB decomp = CHOLESKY; // Decomposition type inside BB algorithm
const bool choiceDual = false; // Choose whether to use primal or dual basis

IntMat basis1;
Rank1Lattice<Int, double> *korlat;  // The original Korobov lattice.
IntLattice<Int, double> *lat;       // A copy of the lattice whose basis can be modified.
Reducer<Int, Real> *red;            // The Reducer object that we use.

long dim; //Saves the current dimension

clock_t tmp;
clock_t totalTime;  // Global timer for total time.
clock_t timer[numMeth][numSizes];
clock_t total_times[numSizes];

std::string names[numMeth] = { "    None ", "pairwise    ", "LLL5        ", "LLL8    ",
		"LLL99999  ", "BKZ   " };

void copyLattice (const bool & du) {
   // du decides if the primal (false) or dual (true) lattice is used
   if (!du)
      copy(korlat->getBasis(), basis1);
   else
      copy(korlat->getDualBasis(), basis1);
   lat = new IntLattice<Int, Real>(basis1, m, dim);
   red->setIntLattice(*lat);
   lat->updateVecNorm(); 	  
}

void shortestNone (IntLattice<Int, Real> & lattice) {
    // Finding shortest vector without pre-reduction takes too much in big dimensions
   if (dim < 30) {
      if (!red->shortestVector(*lat)) {
         std::cout << " shortestVector failed with no pre-reduction, dim  = " << dim << "\n";
      }
   }	  
}

void shortestDieter (IntLattice<Int, Real> & lattice) {
   red->redDieter(0);
   if (dim < 30) {
      if (!red->shortestVector(*lat)) {
  		   std::cout << " shortestVector failed for pairwise pre-reduction with dim  = " << dim << "\n";
      }
   }
}

void shortestBKZ (IntLattice<Int, Real> & lattice) {
   red->redBKZ(lattice.getBasis());
   if (!red->shortestVector(lattice))
      std::cout << " shortestVector failed for BKZ \n";
}

void shortestLLL (IntLattice<Int, Real> & lattice, double delta) {
   red->redLLLNTL(lattice.getBasis(), delta);
   if (!red->shortestVector(lattice))
      std::cout << " shortestVector failed for LLL " << delta << "\n";  
}

int main() {
  lat = new IntLattice<Int, Real>(m, dimensions[numSizes]);
  red = new Reducer<Int, Real>(*lat);
  red->setDecompTypeBB(decomp);
  for (int p = 0; p < numPrimes; p++) {
	  totalTime = clock();
	  // Go through the list of prime numbers and set a accordingly
	  m = primes[p];
  	  a = m/7;
  	  for (int d = 0; d < numSizes; d++) {
  
  		  dim = dimensions[d];
  
  		  // Build the Rank1Lattice in dimension dim 
  		  korlat = new Rank1Lattice<Int, Real>(m, a, dim, choiceDual);
  		  korlat->buildBasis(dim);
  		  basis1.SetDims(dim, dim);
  		  
        copyLattice(choiceDual);        
  	    tmp = clock();
  		  shortestNone(*lat);
  	    timer[0][d] = clock() - tmp;
     
        copyLattice(choiceDual);
  		  tmp = clock();
  		  shortestDieter(*lat);
		    timer[1][d] = clock() - tmp;

		    // Finding shortest vector with LLL reduction and different values of delta (0.5, 0.8, 0.99999)
        copyLattice(choiceDual);
        tmp = clock();
  		  shortestLLL(*lat, 0.5);
  		  timer[2][d] = clock() - tmp;
   
        copyLattice(choiceDual);
  		  tmp = clock();
  		  shortestLLL(*lat, 0.8);
        timer[3][d] = clock() - tmp;
	   
        copyLattice(choiceDual);
  		  tmp = clock();
  		  shortestLLL(*lat, 0.99999);
  		  timer[4][d] = clock() - tmp;
	  
  		  // Finding shortest vector using BKZ pre-reduction
        copyLattice(choiceDual);
  		  tmp = clock();
        shortestBKZ(*lat);
  		  timer[5][d] = clock() - tmp;
  	  }
  
      // Create output string
  	  std::cout << "Flexible types: " << strFlexTypes << "\n";
  	  std::cout << "m = " << m << "\n";
  	  std::cout << "Timings for finding shortest vector with BB, with various pre-reduction methods.\n";
  	  std::cout << "The times are in basic clock units. \n";
  	  std::cout << " dim \\ pre-reduct:  ";
  	  for (int meth = 0; meth < numMeth; meth++)
	  	  std::cout << names[meth] << " ";
  	  std::cout << std::endl << std::endl;
  	  for (int d = 0; d < numSizes; d++) {
  		  std::cout << dimensions[d] << std::setw(10); 
  		  for (int meth = 0; meth < numMeth; meth++) 
  			  std::cout << std::setw(colwidth[meth]) << timer[meth][d] << " ";	  
  		  std::cout << std::endl;
  	  }
  	  std::cout << "Total time: "
  		<< (double) (clock() - totalTime) / (CLOCKS_PER_SEC)
		<< " seconds\n\n\n";
  }
  return 0;
}
