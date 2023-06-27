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
 *
Types: Int = NTL::ZZ, Real = double
Prime: 1048573
Timings for finding shortest vector with different reduction methods, in basic clock units 
 dim \ meth:  w/o Reduction  LLL 0.5     LLL 0.8     LLL 0.9     BKZ         

10                    80367          168         173        190         363 
20                   564200         1469        1451       1349        2739 
30                  3264998       102590       54008      42885       23799 
40                        0       687940       95834      95443       49089 
Total time: 5.07194 seconds
Types: Int = NTL::ZZ, Real = double
Prime: 1073741827
Timings for finding shortest vector with different reduction methods, in basic clock units 
 dim \ meth:  w/o Reduction  LLL 0.5     LLL 0.8     LLL 0.9     BKZ         

10                  2615400          206         246        255         509 
20                  6593212         2631        1666       1713        3706 
30                 34397831       757446      100025      34368       20861 
40                        0    248498019    17401999   12248034     3865791 
Total time: 326.547 seconds
**
*/

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
const long colwidth[numMeth] = { 25, 12, 11, 10, 11 };
const long primes[numPrimes] = { 1048573, 1073741827 };

IntMat basis1;
Rank1Lattice<Int, double> *korlat;  // The original Korobov lattice.
IntLattice<Int, double> *lat;       // A copy of the lattice whose basis can be modified.
Reducer<Int, Real> *red;            // The Reducer object that we use.

clock_t tmp;
clock_t totalTime;  // Global timer for total time.
clock_t timer[numMeth][numSizes];
clock_t total_times[numSizes];

std::string names[numMeth] = { "  None  ", "pairwise   ", "LLL5    ", "LLL8    ",
		"LLL99999    ", "BKZ        " };


int main() {
  for (int p = 0; p < numPrimes; p++) {
	  totalTime = clock();
			  
	  m = primes[p];
  
  	  a = m/7;
  
  	
  	  for (int d = 0; d < numSizes; d++) {
  
  		  long dim = dimensions[d];
  
  		  korlat = new Rank1Lattice<Int, Real>(m, a, dim);
  		  korlat->buildBasis(dim);
  		  basis1.SetDims(dim, dim);
  		  copy(korlat->getBasis(), basis1);

  		  lat = new IntLattice<Int, Real>(basis1, m, dim);
  		  red = new Reducer<Int, Real>(*lat);
  		  lat->updateVecNorm(); 

  	      tmp = clock();
  	      if (dim < 30)
  			  if (!red->shortestVector(*lat))
  				  std::cout << " shortestVector failed with no pre-reduction, dim  = " << dim << "\n";
  		  timer[0][d] = clock() - tmp;
  		  //	  std::cout << " The shortest vector length wih Cholesky decomposition\n";
  		  //	  std::cout << red->getMinLength() << std::endl;
  		  //	  std::cout << "The time to compute shortest vector WITHOUT redcution= " <<(double)(tmp)/(CLOCKS_PER_SEC)<<"second"<<std::endl;
     
  		  tmp = clock();
  		  red->redDieter(0);
  	      if (dim < 30)
  		     if (!red->shortestVector(*lat))
  			    std::cout << " shortestVector failed for pairwise pre-red with dim  = " << dim << "\n";
		  timer[1][d] = clock() - tmp;
  		  //	  std::cout << " The shortest vector length wih Cholesky decomposition\n";
  		  //	  std::cout << red->getMinLength() << std::endl;
  		  //	  std::cout << "The time to compute shortest with LLL 0.5 reduction vector = " <<(double)(tmp)/(CLOCKS_PER_SEC)<<"second"<<std::endl;

  		  tmp = clock();
  		  red->redLLLNTL(lat->getBasis(), 0.5);
  		  if (!red->shortestVector(*lat))
  			  std::cout << " shortestVector failed for LLL 0.5 \n";
  		  timer[2][d] = clock() - tmp;
  		  //	  std::cout << " The shortest vector length wih Cholesky decomposition\n";
  		  //	  std::cout << red->getMinLength() << std::endl;
  		  //	  std::cout << "The time to compute shortest with LLL 0.5 reduction vector = " <<(double)(tmp)/(CLOCKS_PER_SEC)<<"second"<<std::endl;
   
  		  lat = new IntLattice<Int, Real>(basis1, m, dim);
  		  red = new Reducer<Int, Real>(*lat);
  		  lat->updateVecNorm();  	  
  		  tmp = clock();
  		  red->redLLLNTL(lat->getBasis(), 0.8);
  		  if (!red->shortestVector(*lat))
  			  std::cout << " shortestVector failed for LLL 0.8 \n";
  		  timer[3][d] = clock() - tmp;
  		  //	  std::cout << " The shortest vector length wih Cholesky decomposition\n";
  		  //	  std::cout << red->getMinLength() << std::endl;
  		  //	  std::cout << "The time to compute shortest with LLL 0.8 reduction vector = " <<(double)(tmp)/(CLOCKS_PER_SEC)<<"second"<<std::endl;
	   
  		  lat = new IntLattice<Int, Real>(basis1, m, dim);
  		  red = new Reducer<Int, Real>(*lat);
  		  lat->updateVecNorm();  	 	  
  		  tmp = clock();
  		  red->redLLLNTL(lat->getBasis(), 0.99999);
  		  if (!red->shortestVector(*lat))
  			  std::cout << " shortestVector failed for LLL 0.99999 \n";
  		  timer[4][d] = clock() - tmp;
  		  //	  std::cout << " The shortest vector length wih Cholesky decomposition\n";
  		  //	  std::cout << red->getMinLength() << std::endl;
  		  //	  std::cout << "The time to compute shortest with LLL 0.9 reduction vector = " <<(double)(tmp)/(CLOCKS_PER_SEC)<<"second"<<std::endl;
	  
  		  lat = new IntLattice<Int, Real>(basis1, m, dim);
  		  red = new Reducer<Int, Real>(*lat);
  		  lat->updateVecNorm();	    	  
  		  tmp = clock();
  		  red->redBKZ(lat->getBasis());
  		  if (!red->shortestVector(*lat))
  			  std::cout << " shortestVector failed for BKZ \n";
  		  timer[5][d] = clock() - tmp;
  		  //	  std::cout << " The shortest vector length wih Cholesky decomposition\n";
  		  //	  std::cout << red->getMinLength() << std::endl;
  		  //	  std::cout << "The time to compute shortest with BKZ reduction vector = " <<(double)(tmp)/(CLOCKS_PER_SEC)<<"second"<<std::endl;
  	  }
  
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
