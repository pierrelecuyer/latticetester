/**
 * This example illustrates the usage of the BasisConstruction module.
 *
 * TO REWRITE:   *****
 *
 * This reads matrices from files and builds a basis and a dual for an `IntLattice`
 * object. The files this is set to use are in the `bench.zip` archive. To
 * execute the program, the archive should be unziped and the `bench` folder
 * should be put in the same directory from which the executable is called.

 * The bases we used to showcase of BasisConstruction methods are in a folder named 
 * 'examples/bench/'. Each file in 'examples/bench/' folder contain a basis, and the 
 * file is nameed as follows: 'prime_dimBasis_exanpleNumber' where 'prime' is modulo 
 * value of the basis, 'dimBasis' is the dimension of the basis, and 'exampleNumber' 
 * is the number of the example for the bases of dimension 'dimBasis'.
 *
 * This example reads matrices from files and performs the different construction
 * algorithms in BasisConstruction on them. The program then prints the execution
 * time of the various algorithms. Note that the execution of the program is not
 * what you would expect in reality since bench contains random full matrices.
 *
 * We show a use of BasisContruction::upperTriangularBasis, 
 * BasisContruction::lowerTriangularBasis, BasisContruction::LLLConstruction with 
 * two differrent parametter of delta, BasisContruction::mDualUpperTriangular, 
 *  BasisContruction::mDualBasis
 *
 * In this example, we can compare the speed of BasisConstruction<Int>::calcDual method
 * which compute an m-dual basis using any basis in input,
 * and BasisConstruction::mDualUpperTriangular method which compute an m-dual basis
 * with an upper triangular basis.
 *
 * We can also compare the speed of 'BasisConstruction::upperTriangularBasis'
 * and the speed of 'BasisConstruction::LLLConstruction'
 *
 * Example of results with m = 1048573 (prime modulus near 2^{20}):
 *
Types: Int = NTL::ZZ, Real = double
Timings for different methods, in basic clock units
 dim:          10       20       30       40       50

LLL5          3053     6595    12300    20728    29432
LLL8          3562    13129    23420    38350    55714
LLL99999      4175    22895    47833    74681   108213
UppTri         829     2585     3682     4819     6492
Tri96         2016     7068    11316    15016    17853
mDualUT        258     1074     3368     7607    13729
mDualUT96      424     1342     3301     6734    11970
mDual         7081    33344    91112   180249   312506
Total time: 1.22782 seconds

With our LLL_FPZZtest:
LLL5          2543     6599    11889    19372    28548
LLL8          3370    12650    22586    36106    52300
LLL99999      3696    20513    45761    68970    99438

With our LLL_PFZZflex:
LLL5          2359     6908    11512    20036    29576
LLL8          3208    12281    22975    38997    53318
LLL99999      3789    21911    44474    69750    98466


Types: Int = int64_t, Real = double
Timings for different methods, in basic clock units
 dim:          10       20       30       40       50

LLL5          1424     5465    10151    18259    27243
LLL8          1678    11336    22571    41384    65910
LLL99999      2064    20915    55015    96772   156141
UppTri          92      319      520      746     1017
Tri96          133      591     1006     1522     2015
mDualUT         50      203      383      772     1303
mDualUT96       87      481     1142     2820     4460
Total time: 0.559324 seconds

**/

// #define TYPES_CODE  LD     // Int == int64_t
#define TYPES_CODE  ZD     // Int == ZZ

#include <iostream>
#include <cstdint>
#include <ctime>
#include <type_traits>
#include <typeinfo>

#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include "latticetester/FlexTypes.h"    // This defines Int and Real
#include "latticetester/EnumTypes.h"
#include "latticetester/BasisConstruction.h"
#include "latticetester/Util.h"
#include "latticetester/Chrono.h"
// #include "latticetester/ParamReader.h"
#include "latticetester/IntLattice.h"
#include "latticetester/Rank1Lattice.h"
//
#include "latticetester/LLL_FPZZflex.h"


using namespace LatticeTester;

//Int m(101);      // Modulus m = 101
//Int m(1021);     // Modulus m = 1021
Int m(1048573);  // Prime modulus near 2^{20}
//Int m(1073741827);  // Prime modulus near 2^{30}
//Int m(1099511627791);  // Prime modulus near 2^{40}
//Int m(1125899906842597);  // Prime modulus near 2^{50}
Int a;       // The LCG multiplier

const long numSizes = 5;    // Number of matrix sizes (choices of dimension).
const long dimensions[numSizes] = { 10, 20, 30, 40, 50 };
//const long dimensions[numSizes] = {5};

const long numRep = 50;  // Number of replications (multipliers) for each case.
const long numMeth = 8;    // Number of methods, and their names.
std::string names[numMeth] = { "LLL5     ",	"LLL8     ", "LLL99999 ",
		 "UppTri   ", "Tri96    ", "mDualUT  ", "mDualUT96", "mDual    " };


int main() {
	// Here we use ctime directly for the timings, to minimize overhead.
	clock_t totalTime = clock();  // Global timer for total time.
	clock_t timer[numMeth][numSizes];
	clock_t tmp;
	// Chrono totTime;
	// totTime.init();
	IntMat basis1, basis2, basis3, basisdual;
	// Int sqlength;
	Rank1Lattice<Int, Real> *korlat;    // Will be a Korobov lattice.

	long d;
	for (d = 0; d < numSizes; d++) {  // Each matrix size
		long dim = dimensions[d]; // The corresponding dimension.
		basis1.SetDims(dim, dim); // Will be initial triangular basis.
		basis2.SetDims(dim, dim); // Will be LLL-reduced basis.
		basis3.SetDims(dim, dim);
		basisdual.SetDims(dim, dim);  // m-dual basis.
		for (int64_t meth = 0; meth < numMeth; meth++)
			timer[meth][d] = 0;
		for (int64_t r = 0; r < numRep; r++) {
			a = (m / 5 + 10 * r) % m;   // The multiplier we use for this rep.
			korlat = new Rank1Lattice<Int, Real> (m, a, dim);
			korlat->buildBasis(dim);
			copy(korlat->getBasis(), basis1); // This initial basis is triangular.

			// We apply LLL to basis1 in basis2.
			copy(basis1, basis2);
			tmp = clock();
			BasisConstruction<Int>::LLLConstruction0(basis2, 0.5);
			timer[0][d] += clock() - tmp;

			copy(basis1, basis2);
			tmp = clock();
			BasisConstruction<Int>::LLLConstruction0(basis2, 0.8);
			timer[1][d] += clock() - tmp;

			copy(basis1, basis2);
			tmp = clock();
			BasisConstruction<Int>::LLLConstruction0(basis2, 0.99999);
			timer[2][d] += clock() - tmp;

            // We now construct an upper-triangular basis from basis2.
			// We copy basis2 into basis3, because it will be modified.
			copy(basis2, basis3);
			tmp = clock();
			BasisConstruction<Int>::upperTriangularBasis(basis3, basis1, m);
			timer[3][d] += clock() - tmp;

			// Again. This function is in Util.h, it is the old method from 1996.  ???
			copy(basis2, basis3);
			tmp = clock();
			Triangularization(basis3, basis1, dim, dim, m);
			timer[4][d] += clock() - tmp;
			// This basis1 is upper triangular.

            // Now we compute an m-dual basis with various methods.
			tmp = clock();
			BasisConstruction<Int>::mDualUpperTriangular(basis1, basisdual, m);
			timer[5][d] += clock() - tmp;

			tmp = clock();
			BasisConstruction<Int>::mDualUpperTriangular96(basis1, basisdual, m);
			timer[6][d] += clock() - tmp;

#if TYPES_CODE  ==  ZD
			// mDualBasis is currently implemented only for Int = ZZ.
			tmp = clock();
			BasisConstruction<Int>::mDualBasis(basis2, basisdual, m);
			timer[7][d] += clock() - tmp;
#endif
			delete korlat;
			}
	}
	std::cout << "Results of BasisManipulation.cc with m = " << m << "\n";
	std::cout << "Types: " << strFlexTypes << "\n";
	std::cout << "Timings for different methods, in basic clock units \n";
	std::cout << " dim:    ";
	for (d = 0; d < numSizes; d++)
		std::cout << std::setw(8) << dimensions[d] << " ";
	std::cout << std::endl << std::endl;
	for (int meth = 0; meth < numMeth; meth++) {
		std::cout << names[meth] << " ";
		for (d = 0; d < numSizes; d++)
			std::cout << std::setw(8) << timer[meth][d] << " ";
		std::cout << std::endl;
	}
	std::cout << std::endl;
	std::cout << "Total time: "
			<< (double) (clock() - totalTime) / (CLOCKS_PER_SEC)
			<< " seconds\n";
	// double time = totTime.val (Chrono::SEC);
	// std::cout << "Total time Chrono: " << time << " seconds\n";
	// totTime.write(Chrono::SEC);
	return 0;
}

