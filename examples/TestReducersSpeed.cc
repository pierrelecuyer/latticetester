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

#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include "latticetester/FlexTypes.h"
#include "latticetester/Util.h"
#include "latticetester/Rank1Lattice.h"
#include "latticetester/ReducerStatic.h"
#include "latticetester/ReducerBB.h"

using namespace LatticeTester;

const long dimensions[] = { 5, 10, 20, 30, 40, 50, 60, 70 };
const long maxNumSizes = 5; // Number of matrix sizes (choices of dimension), can be adjusted. ***

std::string methNames[] = { "LLL5           ", "LLL99999       ", "BKZ99999-10    ", "L5+9+BKZ-10    ",
      "LLL5+BB        ", "LLL8+BB        ", "LLL99999+BB    ", "BKZ99999-6+BB  ", "BKZ99999-8+BB  ",
      "BKZ99999-10+BB ", "BKZ99999-12+BB ", "BKZ999-6+BB    ", "BKZ999-8+BB    ", "BKZ999-10+BB   ",
      "BKZ999-12+BB   ", "L8+BKZ-10+BB   ", "L5+9+BKZ-10+BB ", "L5+9+BKZ-12+BB ", };
const long numMeth = 18;   // Number of methods, and their short names.

// We use ctime directly for the timings.
clock_t tmp;
clock_t totalTime;  // Global timer for total time.
clock_t timer[numMeth][maxNumSizes]; // Clock times in microseconds.
double sumSq[numMeth][maxNumSizes]; // Sum of squares of vector lengths (for checking).
long numNodes[numMeth][maxNumSizes]; // Total number of calls to tryZ.
long numLeaves[numMeth][maxNumSizes]; // Total number of leaves visited by the BB.

// void printResults(long numMeth, long numSizes);  // Must be declared if it has no template heading.

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
      long d, long dim, long meth, double deltaLLL1, double deltaLLL2, double deltaBKZ, long k, bool BB,
      NTL::Vec<Real> sqlen) {
   if (inDual) {
      korlat.buildDualBasis(dim);  // Rebuild the dual basis (only) anew.
      korlat.dualize();
   } else {
      korlat.buildBasis(dim);  // Rebuild the primal basis anew.
   }
   //std::cout << "Before pre-reduction, initial basis = \n" << korlat.getBasis() << "\n";
   tmp = clock();
   if (deltaLLL1 > 0.0) redLLL(korlat.getBasis(), deltaLLL1, dim, &sqlen);
   if (deltaLLL2 > 0.0) redLLL(korlat.getBasis(), deltaLLL2, dim, &sqlen);
   if (deltaBKZ > 0.0) redBKZ(korlat.getBasis(), deltaBKZ, k, 0, dim, &sqlen);
   double len2 = conv<double>(sqlen[0]);   // This is always the squared L2 norm.
   if (BB) {
      // Here we get the square norm of our choice, either L1 or l2.
      if (red.shortestVector(korlat)) {
         len2 = conv<double>(red.getMinLength2());
         numNodes[meth][d] += red.getCountNodes();
         numLeaves[meth][d] += red.getCountLeaves();
      }
      else std::cout << " shortestVector failed for " << methNames[meth] << "\n";
   }
   timer[meth][d] += clock() - tmp;
   sumSq[meth][d] += len2;
}

// Speed test for dim = dimensions[d], with given matrices.
// We test several methods to approximate or find a shortest vector in the dual lattice.
// The same initial dual basis is rebuilt each time by `performReduction`.
template<typename Int, typename Real>
static void compareManyReductions(Rank1Lattice<Int, Real> &korlat, ReducerBB<Int, Real> &red,
      bool inDual, long d, long dim, NTL::Vec<Real> sqlen) {
   // For the first 4 parameter choices, we take BB = false (no BB is done).
   // performReduction(korlat, red, inDual, d, 0, 0.0, 0.0, 0.0, 1, true, sqlen);
   performReduction(korlat, red, inDual, d, dim, 0, 0.5, 0.0, 0.0, 1, false, sqlen);
   performReduction(korlat, red, inDual, d, dim, 1, 0.99999, 0.0, 0.0, 1, false, sqlen);
   performReduction(korlat, red, inDual, d, dim, 2, 0.0, 0.0, 0.99999, 10, false, sqlen);
   performReduction(korlat, red, inDual, d, dim, 3, 0.5, 0.9, 0.99999, 10, false, sqlen);

   // For the other choices, we take BB = true.
   // For large m, the next two cases are much too slow and often fail.
   if (korlat.getModulus() <= 1024 * 1024) {
      performReduction(korlat, red, inDual, d, dim, 4, 0.5, 0.0, 0.0, 1, true, sqlen);
      performReduction(korlat, red, inDual, d, dim, 5, 0.8, 0.0, 0.0, 1, true, sqlen);
   }
   performReduction(korlat, red, inDual, d, dim, 6, 0.99999, 0.0, 0.0, 1, true, sqlen);
   performReduction(korlat, red, inDual, d, dim, 7, 0.0, 0.0, 0.99999, 6, true, sqlen);

   performReduction(korlat, red, inDual, d, dim, 8, 0.0, 0.0, 0.99999, 8, true, sqlen);
   performReduction(korlat, red, inDual, d, dim, 9, 0.0, 0.0, 0.99999, 10, true, sqlen);
   performReduction(korlat, red, inDual, d, dim, 10, 0.0, 0.0, 0.99999, 12, true, sqlen);
   performReduction(korlat, red, inDual, d, dim, 11, 0.0, 0.0, 0.999, 6, true, sqlen);
   performReduction(korlat, red, inDual, d, dim, 12, 0.0, 0.0, 0.999, 8, true, sqlen);
   performReduction(korlat, red, inDual, d, dim, 13, 0.0, 0.0, 0.999, 10, true, sqlen);
   performReduction(korlat, red, inDual, d, dim, 14, 0.0, 0.0, 0.999, 12, true, sqlen);
   performReduction(korlat, red, inDual, d, dim, 15, 0.8, 0.0, 0.99999, 10, true, sqlen);
   performReduction(korlat, red, inDual, d, dim, 16, 0.5, 0.9, 0.99999, 10, true, sqlen);
   performReduction(korlat, red, inDual, d, dim, 17, 0.5, 0.9, 0.99999, 12, true, sqlen);
}

// In this testing loop, we generate `numRep` multipliers `a` and for each one
// we call tryManyMathods.  We use the same sequence of multipliers `a` for all methods.
template<typename Int, typename Real>
static void testLoop(Int m, NormType norm, DecompTypeBB decomp, bool inDual,
      long numSizes, long numRep) {
   std::string stringTypes;  // To print the selected flexible types.
   strTypes<Int, Real>(stringTypes);  // Functions from FlexTypes
   std::cout << "****************************************************\n";
   std::cout << "TestReducersSpeed with m = " << m << "\n";
   std::cout << "Types: " << stringTypes << "\n";
   std::cout << "Number of replications (different multipliers a): " << numRep << "\n";
   if (inDual) std::cout << "DUAL lattice,  ";
   else std::cout << "PRIMAL lattice,  ";
   std::cout << "Norm: " << toStringNorm(norm) << ",  ";
   std::cout << "Decomposition: " << toStringDecomp(decomp) << ".\n\n";
   long d;  // dim = dimensions[d].
   long maxdim = dimensions[numSizes - 1];   // Maximum dimension
   Rank1Lattice<Int, Real> korlat(m, maxdim, norm); // We use single lattice object.
   ReducerBB<Int, Real> red(korlat);   // Single ReducerBB with internal lattice `korlat`.
   red.setDecompTypeBB(decomp);
   red.setVerbosity(4);  // **************

   NTL::Vec<Real> sqlen; // Cannot be global because it depends on Real.
   sqlen.SetLength(1);   // We retrieve only the shortest vector square length.
   Int a0(73);
   // Int a0(7);
   Int a(a0);   // For the LCG multiplier, we take successive powers of a0 mod m.
   for (d = 0; d < numSizes; d++)   // Reset the accumulators.
      for (int64_t meth = 0; meth < numMeth; meth++) {
         timer[meth][d] = 0;
         sumSq[meth][d] = 0.0;
         numNodes[meth][d] = 0;
         numLeaves[meth][d] = 0;
      }
   totalTime = clock();
   for (int64_t r = 0; r < numRep; r++) {
      korlat.seta(a);
      for (d = 0; d < numSizes; d++)   // Each matrix size.
         performReduction(korlat, red, inDual, d, dimensions[d], 6, 0.99999, 0.0, 0.0, 1, true, sqlen);
         // compareManyReductions<Int, Real>(korlat, red, inDual, d, dimensions[d], sqlen);
      a = a * a0 % m;   // The multiplier we use for this rep. First one is 113.
      }
   printTables<Int, Real>(numMeth, numSizes, numRep, dimensions);
}

template<typename Int, typename Real>
void printTables(long numMeth, long numSizes, long numRep, const long *dimensions) {
   long d;
   std::cout << "Computing times in microseconds:\n";
   std::cout << "Num. dimensions:";
   for (d = 0; d < numSizes; d++)
      std::cout << std::setw(8) << dimensions[d] << "  ";
   std::cout << "\n\n";
   for (int meth = 0; meth < numMeth; meth++) {
      if (timer[meth][0] > 0) { // Only if we have done this method.
         std::cout << methNames[meth] << " ";
         for (d = 0; d < numSizes; d++)
            std::cout << std::setw(9) << timer[meth][d] << " ";
         std::cout << "\n";
      }
   }
   std::cout << "\n";
   std::cout << "Average square length of shortest basis vector:\n";
   std::cout << "Num. dimensions: ";
   for (d = 0; d < numSizes; d++)
      std::cout << std::setw(8) << dimensions[d] << "   ";
   std::cout << "\n\n";
   for (int meth = 0; meth < numMeth; meth++) {
      if (sumSq[meth][0] > 0) {
         std::cout << methNames[meth];
         for (d = 0; d < numSizes; d++)
            std::cout << std::setw(10) << std::setprecision(10) << sumSq[meth][d] / numRep << " ";
         std::cout << "\n";
      }
   }
   std::cout << "\n";
   std::cout << "Average number of calls to the recursive BB procedure `tryZ`:\n";
   std::cout << "Num. dimensions: ";
   for (d = 0; d < numSizes; d++)
      std::cout << std::setw(8) << dimensions[d] << "   ";
   std::cout << "\n\n";
   for (int meth = 0; meth < numMeth; meth++) {
      if (numNodes[meth][0] > 0) {
         std::cout << methNames[meth];
         for (d = 0; d < numSizes; d++)
            std::cout << std::setw(10) << std::setprecision(10) << numNodes[meth][d] / numRep << " ";
         std::cout << "\n";
      }
   }
   std::cout << "\n";
   std::cout << "Average number of visited leaves in the BB procedure:\n";
   std::cout << "Num. dimensions: ";
   for (d = 0; d < numSizes; d++)
      std::cout << std::setw(8) << dimensions[d] << "   ";
   std::cout << "\n\n";
   for (int meth = 0; meth < numMeth; meth++) {
      if (numNodes[meth][0] > 0) {
         std::cout << methNames[meth];
         for (d = 0; d < numSizes; d++)
            std::cout << std::setw(10) << std::setprecision(10) << numLeaves[meth][d] / numRep << " ";
         std::cout << "\n";
      }
   }
   std::cout << "\n";
   std::cout << "Total time for everything: " << (double) (clock() - totalTime) / (CLOCKS_PER_SEC)
         << " seconds\n\n";
}

// This function compares the speed with different `Real` types.
// For each type, it compares different pre-reduction strategies.
template<typename Int, typename Real>
void comparePreRed (Int m, NormType norm, DecompTypeBB decomp, long numSizes, long numRep) {
   std::cout << "\n========================================================================\n";
   std::cout << "Compare different reduction strategies for different real types.\n";
   testLoop<Int, Real>(m, norm, decomp, false, numSizes, numRep);
   testLoop<Int, Real>(m, norm, decomp, true, numSizes, numRep);
}

int main() {
   // NTL::ZZ m(1021);  // Prime modulus near 2^{10}
   NTL::ZZ m(1048573);  // Prime modulus near 2^{20}
   // NTL::ZZ m(1099511627791);  // Prime modulus near 2^{40}
   DecompTypeBB decomp = CHOLESKY;
   NormType norm = L2NORM;
   long numSizes = 3;
   long numRep = 1;

   // comparePreRed<int64_t, double>(conv<int64_t>(m), norm, decomp, numSizes, numRep);
   comparePreRed<NTL::ZZ, double>(m, norm, decomp, numSizes, numRep);
   //comparePreRed<NTL::ZZ, xdouble>(m, norm, decomp, numSizes, numRep);
   //comparePreRed<NTL::ZZ, quad_float>(m, norm, decomp, numSizes, numRep);
   // comparePreRed<NTL::ZZ, NTL::RR>(m, norm, decomp, numSizes, numRep);
}
