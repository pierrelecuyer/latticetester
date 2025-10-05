/**
 * In this example, we compare different ways of computing the shortest vector
 * in a lattice or in its $m$-dual.  For this, we do the BB after some type of
 * pre-reduction of the basis. We compare the total times to do that
 * (pre-red. + BB) for various types of pre-reductions, with lattices that come
 * from Korobov lattice rules in 10 to 40 dimensions, with prime modulus `m`.
 * For the BB, we use the Cholesky decomposition.
 * See the Lattice Tester Guide for more explanations and results.
 *
 * The `main` must be changed and recompiled to change the value of `m`.
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

const long dimensions[] = { 4, 6, 8, 10, 12, 14, 16, 18 };
const long maxNumSizes = 8;   // Number of matrix sizes (choices of dimension).
long maxdim = dimensions[maxNumSizes - 1];   // Maximum dimension

// We use ctime directly for the timings.
clock_t tmp;
clock_t totalTime;  // Global timer for total time.
clock_t timer[maxNumSizes]; // Clock times in microseconds.
double sumSq[maxNumSizes]; // Sum of squares of vector lengths (for checking).
long numBranch[maxNumSizes]; // Total number of calls to tryZ.
long maxZ[maxNumSizes];      // Max absolute value of a z_j.

// Declaration required.
void printResultsNormsDecomp(long numSizes, const long *dimensions, long numRep);

/* This function builds the primal or dual basis of `korlat` in dim = dimensions[d].
 * If `deltaBKZ > 0` it applies BKZ with this `delta` and the given `k`.
 * Then it applies the BB algorithm to find a shortest vector.
 * The square length of the shortest basis vector is recovered in `len2`.
 */
template<typename Int, typename Real>
void performReduction(Rank1Lattice<Int, Real> &korlat, ReducerBB<Int, Real> &red, bool inDual,
      long d, long dim, double deltaBKZ, long k) {
   if (inDual) {
      korlat.buildDualBasis(dim);  // Rebuild the dual basis (only) anew.
      korlat.dualize();
   } else {
      korlat.buildBasis(dim);  // Rebuild the primal basis anew.
   }
   double len2;
   tmp = clock();
   if (deltaBKZ > 0.0) len2 = conv<double>(redBKZ<Int, Real>(korlat.getBasis(), deltaBKZ, k, 0, dim));
   // In what follows, we get the square norm of our choice, either L1 or l2.
   if (red.shortestVector(korlat)) {
         len2 = conv<double>(red.getMinLength2());
         numBranch[d] += red.getCountNodes();
         maxZ[d] = max (maxZ[d], red.getMaxZj());
         // std::cout << "  after shortestVector, square norm of shortest = " << len2 << "\n";
   }
   else std::cout << " shortestVector failed, nodesBB = " << red.getCountNodes() << "\n";
   timer[d] += clock() - tmp;
   sumSq[d] += len2;
}

// In this testing loop, we generate `numRep` multipliers `a` and for each one
// we call tryManyMathods.  We use the same sequence of multipliers `a` for all methods.
template<typename Int, typename Real>
static void testLoop(Int m, NormType norm, DecompTypeBB decomp, bool inDual,
      long numSizes, long numRep) {
   // std::cout << "*********************************************************\n";
   if (inDual) std::cout << "DUAL lattice,  ";
   else std::cout << "PRIMAL lattice,  ";
   std::cout << "Norm: " << toStringNorm(norm) << ",  ";
   std::cout << "Decomposition: " << toStringDecomp(decomp) << "\n";
   long d;  // dim = dimensions[d].
   Rank1Lattice<Int, Real> korlat(m, maxdim, norm); // We use single lattice object.
   ReducerBB<Int, Real> red(korlat);   // Single ReducerBB with internal lattice `korlat`.
   red.setDecompTypeBB(decomp);
   // red.maxNodesBB = 100000;   // 10^5 nodes max
   // red.setVerbosity(4);

   Int a0(73);
   // Int a0(426593);
   Int a(a0);   // For the LCG multiplier, we take successive powers of a0 mod m.
   for (d = 0; d < numSizes; d++) {  // Reset the accumulators.
      timer[d] = 0;
      sumSq[d] = 0.0;
      numBranch[d] = 0;
      maxZ[d] = 0;
   }
   totalTime = clock();
   for (int64_t r = 0; r < numRep; r++) {
      korlat.seta(a);
      for (d = 0; d < numSizes; d++)  // Each matrix size.
         if ((decomp != TRIANGULAR) | (d < 2) | (m < 20000))  // Skip triangular for large m and d.
            performReduction(korlat, red, inDual, d, dimensions[d], 0.999, 10);
      a = a * a0 % m;   // The multiplier we use for this rep. First one is 73.
      }
   printResultsNormsDecomp(numSizes, dimensions, numRep);
   }

// Compares the speeds and results for the two norms and two decompositions.
template<typename Int, typename Real>
void compareNormsDecomp (Int m, long numSizes, long numRep) {
    std::cout << "\n=======================================================================\n";
    std::cout << "Compare L2 vs L1 norms, and Cholesky vs Triangular decompositions.\n";
    std::string stringTypes;  // To print the selected flexible types.
    strTypes<Int, Real>(stringTypes);  // Functions from FlexTypes
    std::cout << "Types: " << stringTypes << "\n";
    std::cout << "TestNormDecomp with m = " << m << "\n";
    std::cout << "Number of replications (different multipliers a): " << numRep << "\n";
    std::cout << "Timings are in basic clock units (microseconds) \n";
    //  testLoop<Int, Real>(m, L2NORM, TRIANGULAR, true, numSizes, numRep);  return;
    for (bool inDual : {false, true}) {
       testLoop<Int, Real>(m, L2NORM, CHOLESKY, inDual, numSizes, numRep);
       testLoop<Int, Real>(m, L2NORM, TRIANGULAR, inDual, numSizes, numRep);
       testLoop<Int, Real>(m, L1NORM, CHOLESKY, inDual, numSizes, numRep);
       testLoop<Int, Real>(m, L1NORM, TRIANGULAR, inDual, numSizes, numRep);
    }
}

void printResultsNormsDecomp(long numSizes, const long *dimensions, long numRep) {
   long d;
   std::cout << "Num dimens:    ";
   for (d = 0; d < numSizes; d++)
      std::cout << std::setw(10) << dimensions[d] << " ";
   std::cout << "\n";
   std::cout << "Microseconds:  ";
   for (d = 0; d < numSizes; d++)
      std::cout << std::setw(10) << timer[d] << " ";
   std::cout << "\n";
   std::cout << "Aver. squares: ";
   for (d = 0; d < numSizes; d++)
      std::cout << std::setw(10) << std::setprecision(10) << sumSq[d]/numRep << " ";
   std::cout << "\n";
   std::cout << "Aver. calls BB:";
   for (d = 0; d < numSizes; d++)
      std::cout << std::setw(10) << std::setprecision(10) << numBranch[d]/numRep << " ";
   std::cout << "\n";
   std::cout << "Max |z_j|:     ";
   for (d = 0; d < numSizes; d++)
      std::cout << std::setw(10) << std::setprecision(10) << maxZ[d] << " ";
   std::cout << "\n";
   std::cout << "Total time for everything: " << (double) (clock() - totalTime) / (CLOCKS_PER_SEC)
         << " seconds\n\n";
}

int main() {
   // Here, Int and Real are not yet defined.
   NTL::ZZ m(1021);  // Prime modulus near 2^{10}
   //NTL::ZZ m(1048573);  // Prime modulus near 2^{20}
   // NTL::ZZ m(1099511627791);  // Prime modulus near 2^{40}

   compareNormsDecomp<NTL::ZZ, double> (m, 6, 1000);
   // compareNormsDecomp<NTL::ZZ, NTL::RR> (m, 4, 1000);
   // compareNormsDecomp<NTL::ZZ, double> (m, 2, 1);
}

