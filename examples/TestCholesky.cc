/**
 * This example tests the Cholesky decomposition with and without pre-reductions,
 * in the primal and/or m-dual, for Korobov lattices of different sizes.
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

const long dimensions[] = { 5, 30, 40 };
// const long dimensions[] = { 5, 10, 20, 30, 40, 50, 60, 70 };
const long maxNumSizes = 3; // Number of matrix sizes (choices of dimension), can be adjusted. ***

std::string methNames[] = { "No pre-reduction", "LLL5            ", "LLL999          ",
      "BKZ99999-12+BB ", "BKZ9999999-20+BB " };
const long numMeth = 5;   // Number of methods, and their short names.

// We use ctime directly for the timings.
clock_t tmp;
clock_t totalTime;  // Global timer for total time.
clock_t timer[numMeth][maxNumSizes]; // Clock times in microseconds.
double sumSq[numMeth][maxNumSizes]; // Sum of squares of vector lengths (for checking).
long numNodes[numMeth][maxNumSizes]; // Total number of calls to tryZ.
long numLeaves[numMeth][maxNumSizes]; // Total number of leaves visited by the BB.
long maxZ[numMeth][maxNumSizes];     // Max absolute value of z_j.

/* This function builds the dual basis of `korlat` in dim = dimensions[d]
 * and dualizes to put the dual as a primal basis. It then applies the
 * selected reductions. The reduction method has number `meth`,
 * its short name is `names[meth]`, and it is specified by the following parameters:
 * If `deltaLLL1 > 0` it applies LLL with this first `delta`.
 * Then if `deltaLLL2 > 0` it applies LLL with this second `delta`.
 * Then if `deltaBKZ > 0` it applies BKZ with this `delta` and the given `k`.
 * Then if `BB == true` it applies the BB algorithm to find a shortest vector.
 */
template<typename Int, typename Real>
void performReduction(Rank1Lattice<Int, Real> &korlat, ReducerBB<Int, Real> &red, bool inDual,
      long d, long dim, long meth, double deltaLLL1, double deltaLLL2, double deltaBKZ, long k,
      bool doBB, NTL::Vec<Real> sqlen) {
   if (inDual) {
      korlat.buildDualBasis(dim);  // Rebuild the dual basis (only) anew.
      korlat.dualize();
   } else {
      korlat.buildBasis(dim);  // Rebuild the primal basis anew.
   }
   std::cout << std::setprecision(10);
   //std::cout << "Before pre-reduction, initial basis = \n" << korlat.getBasis() << "\n";
   std::cout << "\n**************************************\n";
   std::cout << "Dimension: " << dim << "\n";
   std::cout << "Pre-reduction: " << methNames[meth] << "\n";
   if (deltaLLL1 > 0.0) redLLL(korlat.getBasis(), deltaLLL1, dim, &sqlen);
   if (deltaBKZ > 0.0) redBKZ(korlat.getBasis(), deltaBKZ, k, 0, dim, &sqlen);
   if (doBB && !red.shortestVector(korlat))
      std::cout << " shortestVector failed for " << methNames[meth] << "\n";
}

// This function compares the speed with different `Real` types.
// For each type, it compares different pre-reduction strategies, for the m-dual.
// It uses an LCG with modulus `m` and multiplier `a`.
// We use the same sequence of multipliers `a` for all methods.
template<typename Int, typename Real>
static void testLoop(Int m, Int a, NormType norm, DecompTypeBB decomp, bool inDual,
      long numSizes, bool doBB, double epsBounds, long verbose) {
   std::cout << "\n=============================================================================\n";
   std::cout << "TestCholesky: Tests the Cholesky decomposition and BB with different \n";
   std::cout << "  reduction strategies and different real types, with a trace.\n\n";
   std::string stringTypes;  // To print the selected flexible types.
   strTypes<Int, Real>(stringTypes);// Functions from FlexTypes
   std::cout << "Types: " << stringTypes << "\n\n";
   std::cout << "Modulus m = " << m << "\n";
   std::cout << "Multiplier a = " << a << ".\n";
   if (inDual) std::cout << "DUAL lattice,  ";
   else std::cout << "PRIMAL lattice,  ";
   std::cout << "Norm: " << toStringNorm(norm) << ",  ";
   std::cout << "Decomposition: " << toStringDecomp(decomp) << ".\n";
   long d;// dim = dimensions[d].
   long maxdim = dimensions[numSizes - 1];// Maximum dimension
   Rank1Lattice<Int, Real> korlat(m, maxdim, norm);// We use single lattice object.
   ReducerBB<Int, Real> red(korlat);// Single ReducerBB with internal lattice `korlat`.
   red.setDecompTypeBB(decomp);
   if (!doBB) red.maxNodesBB = 1;
   Real eps = Real(epsBounds);
   red.setEpsBounds(eps);// Safety margin on the bounds in the BB.
   std::cout << std::setprecision(10) << "Safety margin on the BB bounds: epsBounds = " << eps << ".\n";
   red.setVerbosity(verbose);
   NTL::Vec<Real> sqlen;// Cannot be global because it depends on Real.
   sqlen.SetLength(1);// We retrieve only the shortest vector square length.
   korlat.seta(a);
   for (d = 0; d < numSizes; d++) {  // Each matrix size.
      if (!doBB)
         performReduction(korlat, red, inDual, d, dimensions[d], 0, 0.0, 0.0, 0.0, 0, doBB, sqlen);
      // if ((korlat.getModulus() <= 1024 * 1024) || !doBB)
      performReduction(korlat, red, inDual, d, dimensions[d], 1, 0.5, 0.0, 0.0, 1, doBB, sqlen);
      performReduction(korlat, red, inDual, d, dimensions[d], 2, 0.999, 0.0, 0.0, 1, doBB, sqlen);
      performReduction(korlat, red, inDual, d, dimensions[d], 3, 0.99999, 0.0, 0.0, 12, doBB, sqlen);
      if (d > 0)
         performReduction(korlat, red, inDual, d, dimensions[d], 4, 0.9999999, 0.0, 0.0, 20, doBB, sqlen);
   }
}

int main() {
   //NTL::ZZ m(1021);  // Prime modulus near 2^{10}
   //NTL::ZZ m(1048573);  // Prime modulus near 2^{20}
   NTL::ZZ m(1099511627791);  // Prime modulus near 2^{40}
   NTL::ZZ a0(401173573);
   NTL::ZZ a(a0 % m);   // The LCG multiplier.
   DecompTypeBB decomp = CHOLESKY;
   NormType norm = L2NORM;
   long numSizes = 3;  // Number of values to test for the dimension.
   bool doBB = true;  // Perform the BB ?
   double epsBounds = 0.000001;  // Safety margin on the bounds in the BB.
   long verbose = 3;  // Level of detail in the output trace.

   testLoop<NTL::ZZ, double>(m, a, norm, decomp, true, numSizes, doBB, epsBounds, verbose);
   testLoop<NTL::ZZ, quad_float>(m, a, norm, decomp, true, numSizes, doBB, epsBounds, verbose);
   testLoop<NTL::ZZ, NTL::RR>(m, a, norm, decomp, true, numSizes, doBB, epsBounds, verbose);
   std::cout << "\n========================================================================\n";
   std::cout << "We now set the RR precision to 250.\n";
   NTL::RR::SetPrecision(250);
   testLoop<NTL::ZZ, NTL::RR>(m, a, norm, decomp, true, numSizes, doBB, epsBounds, verbose);
}
