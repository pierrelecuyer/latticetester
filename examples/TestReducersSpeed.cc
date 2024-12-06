/**
 * In this example, we compare different ways of computing the shortest vector
 * in a lattice or in its $m$-dual.  For this, we do the BB after some type of
 * pre-reduction of the basis. We compare the total times to do that
 * (pre-red. + BB) for various types of pre-reductions, with lattices that come
 * from Korobov lattice rules in 10 to 40 dimensions, with prime modulus `m`.
 * For the BB, we use the Cholesky decomposition.
 * See the Lattice Tester Guide for more explanations and results.
 *
 * The `main` must be changed and recompiled to change the value of `m` and
 * to switch between primal and m-dual.
 */

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// #define TYPES_CODE  LD     // Int = int64_t, Real = double

#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include "latticetester/FlexTypes.h"
#include "latticetester/Util.h"
#include "latticetester/Rank1Lattice.h"
#include "latticetester/ReducerStatic.h"
#include "latticetester/ReducerBB.h"

using namespace LatticeTester;

// This could also be made as a class, to reduce parameter passing in functions.
// template<typename Int, typename Real>
// class TestReducersSpeed {

// const long dimensions[] = { 5, 10, 20, 30, 40, 50, 60 };
const long dimensions[] = { 6, 8, 10, 12, 14, 16, 18 };
const long maxNumSizes = 7;   // Number of matrix sizes (choices of dimension).
long maxdim = dimensions[maxNumSizes - 1];   // Maximum dimension

const long numMeth = 18;   // Number of methods, and their short names.
std::string names[] = { "LLL5           ", "LLL99999       ", "BKZ99999-10    ", "L5+9+BKZ-10    ",
      "LLL5+BB        ", "LLL8+BB        ", "LLL99999+BB    ", "BKZ99999-6+BB  ", "BKZ99999-8+BB  ",
      "BKZ99999-10+BB ", "BKZ99999-12+BB ", "BKZ999-6+BB    ", "BKZ999-8+BB    ", "BKZ999-10+BB   ",
      "BKZ999-12+BB   ", "L8+BKZ-10+BB   ", "L5+9+BKZ-10+BB ", "L5+9+BKZ-12+BB ", };

// const long numMeth2 = 2;   // Number of methods, and their short names.
// std::string names2[] = { "LLL99999+BB    ", "BKZ999-10+BB   "};

// We use ctime directly for the timings.
clock_t tmp;
clock_t totalTime;  // Global timer for total time.
clock_t timer[numMeth][maxNumSizes]; // Clock times in microseconds.
double sumSq[numMeth][maxNumSizes]; // Sum of squares of vector lengths (for checking).

void printResults(long numMeth, long numSizes);  // Must be declared, because it has no template parameters.

/* This function builds the dual basis of `korlat` in dim = dimensions[d]
 * and dualizes to put the dual as a primal basis. It then applies the
 * selected reductions and adds the CPU times and square length of shortest
 * vector to the appropriate sum. The reduction method has number `meth`,
 * its short name is `names[meth]`, and it is specified by the following parameters:
 * If `deltaLLL1 > 0` it applies LLL with this first `delta`.
 * Then if `deltaLLL2 > 0` it applies LLL with this second `delta`.
 * Then if `deltaBKZ > 0` it applies BKZ with this `delta` and the given `k`.
 * Then if `BB == true` it applies the BB algorithm to find a shortest vector.
 * The square length of the shortest basis vector is recovered in `len2`.
 */
template<typename Int, typename Real>
void performReduction(Rank1Lattice<Int, Real> &korlat, ReducerBB<Int, Real> &red, bool inDual,
      long d, long meth, double deltaLLL1, double deltaLLL2, double deltaBKZ, long k, bool BB,
      NTL::Vec<Real> sqlen) {
   long dim = dimensions[d];
   if (inDual) {
      korlat.buildDualBasis(dim);  // Rebuild the dual basis (only) anew.
      korlat.dualize();
   } else {
      korlat.buildBasis(dim);  // Rebuild the primal basis anew.
   }
   tmp = clock();
   if (deltaLLL1 > 0.0) redLLL(korlat.getBasis(), deltaLLL1, dim, &sqlen);
   if (deltaLLL2 > 0.0) redLLL(korlat.getBasis(), deltaLLL2, dim, &sqlen);
   if (deltaBKZ > 0.0) redBKZ(korlat.getBasis(), deltaBKZ, k, 0, dim, &sqlen);
   double len2 = conv<double>(sqlen[0]);   // This is always the squared L2 norm.
   if (BB) {
      // Here we get the square norm of our choice, either L1 or l2.
      if (red.shortestVector(korlat)) len2 = conv<double>(red.getMinLength2());
      else std::cout << " shortestVector failed for " << names[meth] << "\n";
   }
   // std::cout << " korlat.getVecNorm(0): = " << korlat.getVecNorm(0) << "\n";
   // std::cout << " red.getShortVec(): = " << red.getShortVec() << "\n";
   timer[meth][d] += clock() - tmp;
   sumSq[meth][d] += len2;
}

// Speed test for dim = dimensions[d], with given matrices.
// We test several methods to approximate or find a shortest vector in the dual lattice.
// The same initial dual basis is rebuilt each time by `performReduction`.
template<typename Int, typename Real>
static void tryManyMethods1(Rank1Lattice<Int, Real> &korlat, ReducerBB<Int, Real> &red, bool inDual,
      long d) {
   NTL::Vec<Real> sqlen; // Cannot be global because it depends on Real.
   sqlen.SetLength(1);   // We retrieve only the shortest vector square length.

   // For the first 4 parameter choices, we take BB = false (no BB).
   // performReduction(korlat, red, inDual, d, 0, 0.0, 0.0, 0.0, 1, true, sqlen);
   performReduction(korlat, red, inDual, d, 0, 0.5, 0.0, 0.0, 1, false, sqlen);
   performReduction(korlat, red, inDual, d, 1, 0.99999, 0.0, 0.0, 1, false, sqlen);
   performReduction(korlat, red, inDual, d, 2, 0.0, 0.0, 0.99999, 10, false, sqlen);
   performReduction(korlat, red, inDual, d, 3, 0.5, 0.9, 0.99999, 10, false, sqlen);

   // For the other choices, we take BB = true.
   // We hide the next two cases for the large m, because they are much too slow
   // and they often fail.
   //performReduction(korlat, red, inDual, d, 4, 0.5, 0.0, 0.0, 1, true, sqlen);
   //performReduction(korlat, red, inDual, d, 5, 0.8, 0.0, 0.0, 1, true, sqlen);
   performReduction(korlat, red, inDual, d, 6, 0.99999, 0.0, 0.0, 1, true, sqlen);
   performReduction(korlat, red, inDual, d, 7, 0.0, 0.0, 0.99999, 6, true, sqlen);

   performReduction(korlat, red, inDual, d, 8, 0.0, 0.0, 0.99999, 8, true, sqlen);
   performReduction(korlat, red, inDual, d, 9, 0.0, 0.0, 0.99999, 10, true, sqlen);
   performReduction(korlat, red, inDual, d, 10, 0.0, 0.0, 0.99999, 12, true, sqlen);
   performReduction(korlat, red, inDual, d, 11, 0.0, 0.0, 0.999, 6, true, sqlen);
   performReduction(korlat, red, inDual, d, 12, 0.0, 0.0, 0.999, 8, true, sqlen);
   performReduction(korlat, red, inDual, d, 13, 0.0, 0.0, 0.999, 10, true, sqlen);
   performReduction(korlat, red, inDual, d, 14, 0.0, 0.0, 0.999, 12, true, sqlen);
   performReduction(korlat, red, inDual, d, 15, 0.8, 0.0, 0.99999, 10, true, sqlen);
   performReduction(korlat, red, inDual, d, 16, 0.5, 0.9, 0.99999, 10, true, sqlen);
   performReduction(korlat, red, inDual, d, 17, 0.5, 0.9, 0.99999, 12, true, sqlen);
}

// This version just tries two methods.  We use it to compare the speeds for L1 and L2 norms.
template<typename Int, typename Real>
static void tryManyMethods2(Rank1Lattice<Int, Real> &korlat, ReducerBB<Int, Real> &red, bool inDual,
      long d) {
   NTL::Vec<Real> sqlen; // Cannot be global because it depends on Real.
   sqlen.SetLength(1);   // We retrieve only the shortest vector square length.
   // Here we take BB = true.
   performReduction(korlat, red, inDual, d, 6, 0.99999, 0.0, 0.0, 1, true, sqlen);
   // performReduction(korlat, red, inDual, d, 7, 0.0, 0.0, 0.99999, 6, true, sqlen);
   performReduction(korlat, red, inDual, d, 13, 0.0, 0.0, 0.999, 10, true, sqlen);
}

// In this testing loop, we generate `numRep` multipliers `a` and for each one
// we call tryManyMathods.  We use the same sequence of multipliers `a` for all methods.
template<typename Int, typename Real>
static void testLoop(Int m, NormType norm, bool inDual, long numSet, long numSizes, long numRep) {
   std::string stringTypes;  // To print the selected flexible types.
   strTypes<Int, Real>(stringTypes);  // Functions from FlexTypes
   std::cout << "****************************************\n";
   std::cout << "Types: " << stringTypes << "\n";
   std::cout << "TestReducersSpeed with m = " << m;
   if (inDual) std::cout << ", in the dual lattice ";
   else std::cout << ", in the primal lattice ";
   std::cout << "with norm L" << norm << ".\n\n";
   std::cout << "Timings (in microseconds) for different methods for " << numRep
         << " replications. \n";
   long d;  // dim = dimensions[d].
   Rank1Lattice<Int, Real> korlat(m, maxdim, norm); // We use single lattice object.
   ReducerBB<Int, Real> red(korlat);   // Single ReducerBB with internal lattice `korlat`.
   Int a0(91);
   Int a(a0);   // For the LCG multiplier, we take successive powers of a0 mod m.
   for (d = 0; d < numSizes; d++)   // Reset the accumulators.
      for (int64_t meth = 0; meth < numMeth; meth++) {
         timer[meth][d] = 0;
         sumSq[meth][d] = 0.0;
      }
   totalTime = clock();
   for (int64_t r = 0; r < numRep; r++) {
      a = a * a0 % m;   // The multiplier we use for this rep.
      korlat.seta(a);
      for (d = 0; d < numSizes; d++) {  // Each matrix size.
         if (numSet == 1) tryManyMethods1<Int, Real>(korlat, red, inDual, d);
         else tryManyMethods2<Int, Real>(korlat, red, inDual, d);
      }
   }
   printResults(numMeth, numSizes);
}

void printResults(long numMeth, long numSizes) {
   long d;
   std::cout << "Num. dimensions:";
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
   std::cout << "Num. dimensions:";
   for (d = 0; d < numSizes; d++)
      std::cout << std::setw(8) << dimensions[d] << "   ";
   std::cout << "\n\n";
   for (int meth = 0; meth < numMeth; meth++) {
      if (sumSq[meth][0] > 0) {
         std::cout << names[meth];
         for (d = 0; d < numSizes; d++)
            std::cout << std::setw(10) << std::setprecision(10) << sumSq[meth][d] << " ";
         std::cout << "\n";
      }
   }
   std::cout << "\n";
   std::cout << "Total time for everything: " << (double) (clock() - totalTime) / (CLOCKS_PER_SEC)
         << " seconds\n\n";
}

// This function compares the speed with different `Real` types.
// For each type, it compares different pre-reduction strategies.
template<typename Int>
void compareRealTypes (Int m, NormType norm, bool inDual, long numSet, long numSizes, long numRep) {
   std::cout << "\n========================================================================\n";
   std::cout << "Compare different reduction strategies for different real types, with L2 norm, ";
   if (inDual) std::cout << ", in the dual lattice. \n\n";
   else std::cout << ", in the primal lattice. \n\n";
   testLoop<int64_t, double>(conv<int64_t>(m), norm, inDual, numSet, numSizes, numRep);
   testLoop<NTL::ZZ, double>(m, norm, inDual, numSet, numSizes, numRep);
   //testLoop<NTL::ZZ, xdouble>(m, norm, inDual, numSizes, numRep);
   //testLoop<NTL::ZZ, quad_float>(m, norm, inDual, numSizes, numRep);
   //testLoop<NTL::ZZ, NTL::RR>(m, norm, inDual, numSizes, numRep);
}

// Compares the speeds and results for the two norms, L1 and L2.
template<typename Int>
void compareL1L2 (Int m, bool inDual, long numSet, long numSizes, long numRep) {
    std::cout << "\n=======================================================================\n";
    std::cout << "Compare L2 vs L2 norms ";
    if (inDual) std::cout << ", in the dual lattice. \n\n";
    else std::cout << ", in the primal lattice. \n\n";
    testLoop<NTL::ZZ, double>(m, L2NORM, inDual, numSet, numSizes, numRep);
    testLoop<NTL::ZZ, double>(m, L1NORM, inDual, numSet, numSizes, numRep);
}

int main() {
   // Here, Int and Real are not yet defined.
   NTL::ZZ m(1048573);  // Prime modulus near 2^{20}
   // NTL::ZZ m(1099511627791);  // Prime modulus near 2^{40}
   // NormType norm = L2NORM;
   //long numSizes = 2;
   //long numRep = 2; // Number of replications (multipliers a) for each case.
   //bool inDual = false;  // Tests in dual lattice ?
   // long numSet = 1;  // The set of methods, either 1 or 2.

   // compareL1L2 (m, false, 2, 3, 10);
   compareL1L2 (m, true, 2, 5, 10);

   //compareRealTypes (m, L2NORM, false, 1, 3, 10);
   //compareRealTypes (m, L2NORM, true, 1, 3, 10);
}
