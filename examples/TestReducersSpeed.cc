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
 */

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include "latticetester/FlexTypes.h"
#include "latticetester/Util.h"
#include "latticetester/Rank1Lattice.h"
#include "latticetester/ReducerStatic.h"
#include "latticetester/ReducerBB.h"

using namespace LatticeTester;

// This should be made as a class, to reduce parameter passing in functions.
// template<typename Int, typename Real>
// class TestReducersSpeed {

// We cannot use Int or Real here, because they are not yet defined.
const long dimensions[] = { 5, 10, 20, 30, 40 };
const long numSizes = 5;   // Number of matrix sizes (choices of dimension).
long maxdim = dimensions[numSizes - 1];   // Maximum dimension

const long numMeth = 15;    // Number of methods, and their short names.
std::string names[] = {
      "LLL5          ", "LLL99999      ", "BKZ99999-10   ",
//      "Direct-BB     ", "Pairwise+BB   ",
      "LLL5+BB       ", "LLL8+BB       ", "LLL99999+BB   ",
      "BKZ99999-6+BB ", "BKZ99999-8+BB ", "BKZ99999-10+BB",
      "BKZ999-6+BB   ", "BKZ999-8+BB   ", "BKZ999-10+BB  ",
      "BKZ999-12+BB  ", "L8+BKZ-10+BB  ", "L99+BKZ-10+BB ",
};

// We use ctime directly for the timings, to minimize overhead.
clock_t tmp;
clock_t totalTime;  // Global timer for total time.
clock_t timer[numMeth][numSizes]; // Clock times in microseconds.
double sumSq[numMeth][numSizes]; // Sum of squares of vector lengths (for checking).
std::string stringTypes;  // To print the selected flexible types.

void printResults();  // Must be declared, because it has no parameters.

/* This function builds the dual basis of `korlat` in dim = dimensions[d]
 * and dualizes to put the dual as a primal basis. It then applies the
 * selected reductions and adds the CPU times and square length of shortest
 * vector to the appropriate sum. The reduction method has number `meth`,
 * its short name is `names[meth]`, and it is specified by the following
 * parameters: If `deltaLLL > 0` it applies LLL with this `delta`.
 * Then if `deltaBKZ > 0` it applies BKZ with this `delta` and the given `k`.
 * Then if `BB == true` it applies the BB algorithm to find a shortest vector.
 * The square length of the shortest basis vector is recovered in `len2`.
 */
template<typename Int, typename Real>
void performReduction(Rank1Lattice<Int, Real> *korlat, ReducerBB<Int, Real> *red,
      long d, long meth, double deltaLLL, double deltaBKZ, long k,
      bool BB, NTL::vector<Real> sqlen) {
   long dim = dimensions[d];
   double len2;
   korlat->setPrimalFlag(false);
   korlat->setDualFlag(true);
   korlat->buildDualBasis(dim);  // Rebuild the dual basis (only) anew.
   korlat->dualize();   // This also exchanges the primal / dual flags.
   tmp = clock();
   if (deltaLLL > 0.0)
      redLLL(korlat->getBasis(), deltaLLL, dim, &sqlen);
   if (deltaBKZ > 0.0)
      redBKZ(korlat->getBasis(), deltaBKZ, k, 0, dim, &sqlen);
   len2 = conv<double>(sqlen[0]);
   if (BB) {
      if (!red->shortestVector())
         std::cout << " shortestVector failed for " << names[meth] << "\n";
      len2 = conv<double>(korlat->getVecNorm(0));
   }
   timer[meth][d] += clock() - tmp;
   sumSq[meth][d] += len2;
}


// Speed test for dim = dimensions[d], with given matrices.
// We test several methods to approximate or find a shortest vector in the dual lattice.
// The same initial dual basis is rebuilt each time by `beforeReduct`.
template<typename Int, typename Real>
static void tryManyMethods(Rank1Lattice<Int, Real> *korlat,
      ReducerBB<Int, Real> *red, long d) {
   NTL::vector<Real> sqlen; // Cannot be global because it depends on Real.
   sqlen.SetLength(1);  // With store only the shortest vector square length.

   performReduction(korlat, red, d, 0, 0.5, 0.0, 1, false, sqlen);
   performReduction(korlat, red, d, 1, 0.99999, 0.0, 1, false, sqlen);
   performReduction(korlat, red, d, 2, 0.0, 0.99999, 10, false, sqlen);
   performReduction(korlat, red, d, 3, 0.5, 0.0, 1, true, sqlen);
   performReduction(korlat, red, d, 4, 0.8, 0.0, 1, true, sqlen);
   performReduction(korlat, red, d, 5, 0.99999, 0.0, 1, true, sqlen);
   performReduction(korlat, red, d, 6, 0.0, 0.99999, 6, true, sqlen);
   performReduction(korlat, red, d, 7, 0.0, 0.99999, 8, true, sqlen);
   performReduction(korlat, red, d, 8, 0.0, 0.99999, 10, true, sqlen);
   performReduction(korlat, red, d, 9, 0.0, 0.999, 6, true, sqlen);
   performReduction(korlat, red, d, 10, 0.0, 0.999, 8, true, sqlen);
   performReduction(korlat, red, d, 11, 0.0, 0.999, 10, true, sqlen);
   performReduction(korlat, red, d, 12, 0.0, 0.999, 12, true, sqlen);
   performReduction(korlat, red, d, 13, 0.8, 0.999, 10, true, sqlen);
   performReduction(korlat, red, d, 14, 0.99, 0.999, 10, true, sqlen);
}

// In this testing loop, we generate `numRep` multipliers `a` and for each one
template<typename Int, typename Real>
static void testLoop(Int m, long numRep) {
   strTypes<Int, Real>(stringTypes);  // Functions from FlexTypes
   std::cout << "****************************************************\n";
   std::cout << "Types: " << stringTypes << "\n\n";
   std::cout << "TestReducersSpeed with m = " << m << "\n";
   std::cout << "Results for `testLoop`\n";
   std::cout << "Timings (in microseconds) for different methods for " << numRep
         << " replications \n\n";
   long d;  // dim = dimensions[d].
   Rank1Lattice<Int, Real> *korlat; // We create a single lattice object.
   korlat = new Rank1Lattice<Int, Real>(m, maxdim, false, true);
   ReducerBB<Int, Real> *red;       // Also a single ReducerBB object.
   red = new ReducerBB<Int, Real>(*korlat);
   Int a;        // The LCG multiplier
   // IntMat basisdual1, basisdual2;
   // basisdual1.SetDims(maxdim, maxdim); // Will be initial triangular dual basis.
   // basisdual2.SetDims(maxdim, maxdim); // m-dual basis passed to LLL, etc.
   for (d = 0; d < numSizes; d++)   // Reset accumulators.
      for (int64_t meth = 0; meth < numMeth; meth++) {
         timer[meth][d] = 0;
         sumSq[meth][d] = 0.0;
      }
   totalTime = clock();
   for (int64_t r = 0; r < numRep; r++) {
      a = (m / 5 + 13 * r) % m;   // The multiplier we use for this rep.
      korlat->seta(a);
      for (d = 0; d < numSizes; d++) {  // Each matrix size
         tryManyMethods<Int, Real>(korlat, red, d);
      }
   }
   printResults();
}

void printResults() {
   long d;
   std::cout << " dim:        ";
   for (d = 0; d < numSizes; d++)
      std::cout << std::setw(8) << dimensions[d] << "  ";
   std::cout << "\n\n";
   for (int meth = 0; meth < numMeth; meth++) {
      if (timer[meth][0] > 0) { // Only if we have done this method.
         std::cout << names[meth] << " ";
         for (d = 0; d < numSizes; d++)
            std::cout << std::setw(9) << timer[meth][d] << " ";
         std::cout << "\n";
      }
   }
   std::cout << "\n";
   std::cout << "Sums of square lengths of shortest basis vector:\n";
   std::cout << " dim:     ";
   for (d = 0; d < numSizes; d++)
      std::cout << std::setw(10) << dimensions[d] << " ";
   std::cout << "\n\n";
   for (int meth = 0; meth < numMeth; meth++) {
      if (sumSq[meth][0] > 0) {
         std::cout << names[meth];
         for (d = 0; d < numSizes; d++)
            std::cout << std::setw(10) << sumSq[meth][d] << " ";
         std::cout << "\n";
      }
   }
   std::cout << "\n";
   std::cout << "Total time for everything: "
         << (double) (clock() - totalTime) / (CLOCKS_PER_SEC) << " seconds\n\n";
   std::cout
         << "We see that LLL or BKZ alone do not always find a shortest vector.\n\n\n";
}

int main() {
   // Here, Int and Real are not yet defined.
   NTL::ZZ m(1048573);  // Prime modulus near 2^{20}
   // NTL::ZZ m(1099511627791);  // Prime modulus near 2^{40}
   long numRep = 50;  // Number of replications (multipliers) for each case.

   testLoop<long, double>(conv<long>(m), numRep);
   testLoop<NTL::ZZ, double>(m, numRep);
   testLoop<NTL::ZZ, xdouble>(m, numRep);
   testLoop<NTL::ZZ, quad_float>(m, numRep);
   // testLoop<NTL::ZZ, NTL::RR>(m, numRep);
}
