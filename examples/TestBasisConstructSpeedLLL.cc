// File `testBasisConstructionLLL`

#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/xdouble.h>
#include <NTL/ZZ.h>
#include <NTL/mat_ZZ.h>
#include "latticetester/FlexTypes.h"
#include "latticetester/Util.h"
#include "latticetester/Rank1Lattice.h"
#include "latticetester/BasisConstruction.h"

/**
 * This example makes speed comparisons with the `BasisConstruction` functions,
 * with different combinations of types. See the Lattice Tester guide for more explanations.
 * This experiment concerns mostly LLL reduction in the primal and m-dual.
 * The `main` must be changed and recompiled to change the value of `m`.
 */
using namespace NTL;
using namespace LatticeTester;

// The types Int and Real are not yet defined here.
// They are passed as template parameters from the `main`.
const int64_t maxNumSizes = 8; // Number of matrix sizes (choices of dimensions).
const int64_t dimensions[maxNumSizes] = { 5, 10, 20, 30, 40, 50, 60, 70 };
const int64_t numMeth = 12;    // Number of methods to test for LLL, and their names.
std::string names[numMeth] = { "LLL5         ", "LLL8         ", "LLL99        ", "LLL99999     ",
      "LLL99999-pnew", "UppTri       ", "mDualUT      ", "LLL5-dual    ", "LLL8-dual    ",
      "LLL99-dual   ", "LLL99999-dual", "LLL99999-dnew" };

// We use `ctime` directly for the timings, to minimize overhead.
clock_t tmpTotal = 0;            // Global timer for total time.
clock_t timer[numMeth][maxNumSizes]; // Collects timings for each case.
double sumSq[numMeth][maxNumSizes];  // Sums of square lengths.

// Declaration required.
void printTables(int64_t numSizes);

// This function applies LLL to `basis` in `dim` dimensions.
// It also updates the cumulative times and sums of square lengths.
template<typename Int, typename Real>
void LLLTest(IntMat &basis, int64_t d, int64_t meth, double delta) {
   int64_t dim = dimensions[d];
   NTL::Vec<Real> sqlen; // Cannot be global, because it depends on Real.
   sqlen.SetLength(1);
   clock_t tmp = clock();
   // We apply LLL to the basis.
   // We could equivalently use `redLLL` from the file `ReducerStatic`.
   LLLConstruction0(basis, delta, dim, dim, &sqlen);
   timer[meth][d] += clock() - tmp;
   sumSq[meth][d] += conv<double>(sqlen[0]);
}

// Run a speed test for dim = dimensions[d], with given basis matrices.
// Only basis0 needs to be initialized; the other matrices are used only for copy.
template<typename Int, typename Real>
void transformBasesLLL(const Int &m, int64_t d, int64_t dim, IntMat &basis0,
      IntMat &basis1, IntMat &basis2, IntMat &basisdual) {
   clock_t tmp;
   CopyPartMat<IntMat>(basis1, basis0, dim, dim);  // Copy basis0 to basis1.

   // We apply LLL to basis2 with different values of `delta`, incrementally.
   // We start with delta=0.5, then continue with 0.8, then with 0.99, etc.
   LLLTest<Int, Real>(basis1, d, 0, 0.5);
   LLLTest<Int, Real>(basis1, d, 1, 0.8);
   LLLTest<Int, Real>(basis1, d, 2, 0.99);
   LLLTest<Int, Real>(basis1, d, 3, 0.99999);
   // Now we restart LLL from the initial triangular basis, with delta=0.99999.
   CopyPartMat(basis1, basis0, dim, dim);  // Copy basis0 to basis1.
   LLLTest<Int, Real>(basis1, d, 4, 0.99999);

   // We now construct an upper-triangular basis from basis1 into basis2.
   tmp = clock();
   upperTriangularBasis(basis2, basis1, m, dim, dim);
   timer[5][d] += clock() - tmp;
   // We compute basisdual, the m-dual basis to basis2.
   tmp = clock();
   mDualUpperTriangular(basisdual, basis2, m, dim);
   timer[6][d] += clock() - tmp;

   // We apply LLL to this m-dual basis, with delta = 0.5, 0.8, etc.
   LLLTest<Int, Real>(basisdual, d, 7, 0.5);
   LLLTest<Int, Real>(basisdual, d, 8, 0.8);
   LLLTest<Int, Real>(basisdual, d, 9, 0.99);
   LLLTest<Int, Real>(basisdual, d, 10, 0.99999);
   // We restart anew with delta = 0.99999.
   mDualUpperTriangular(basisdual, basis2, m, dim);
   LLLTest<Int, Real>(basisdual, d, 11, 0.99999);
}

// Testing loop. The `IntMat` and `Rank1Lattice` objects are created only once.
template<typename Int, typename Real>
void testLoop(const Int &mm, int64_t numSizes, int64_t numRep) {
   std::string stringTypes;  // To print the selected flexible types.
   strTypes<Int, Real>(stringTypes);  // Functions from FlexTypes
   std::cout << "****************************************************\n";
   std::cout << "TestBasisConstructSpeedLLL with m = " << mm << "\n";
   std::cout << "Types: " << stringTypes << "\n";
   std::cout << "Number of replications (different multipliers a): " << numRep << "\n";
   int64_t d, dim;  // Index of dimension, and dimension.
   Int m = conv<Int>(mm);
   Int a;
   IntMat basis0, basis1, basis2, basisdual;
   int64_t maxdim = dimensions[numSizes - 1];   // Maximum dimension
   basis0.SetDims(maxdim, maxdim); // Will be initial triangular basis.
   basis1.SetDims(maxdim, maxdim); // Another basis.
   basis2.SetDims(maxdim, maxdim); // Another basis.
   basisdual.SetDims(maxdim, maxdim);  // m-dual basis.
   // We create a single Korobov lattice object.
   Rank1Lattice<Int, Real> korlat(m, maxdim);

   for (d = 0; d < numSizes; d++)   // Reset accumulators.
      for (int64_t meth = 0; meth < numMeth; meth++) {
         timer[meth][d] = 0;
         sumSq[meth][d] = 0.0;
      }
   tmpTotal = clock();
   for (int64_t r = 0; r < numRep; r++) {
      a = (m / 5 + 17 * r) % m;   // The multiplier a used for this rep.
      korlat.seta(a);
      for (d = 0; d < numSizes; d++) {  // Each matrix size
         dim = dimensions[d]; // The corresponding dimension.
         korlat.buildBasis(dim);
         // Resizing the matrices as below makes no significant change in computing times.
         //basis0.SetDims(dim, dim);    basis1.SetDims(dim, dim);
         //basis2.SetDims(dim, dim);    basisdual.SetDims(dim, dim);
         CopyPartMat<IntMat>(basis0, korlat.getBasis(), dim, dim); // Copy korlat basis to basis0.
         transformBasesLLL<Int, Real>(m, d, dim, basis0, basis1, basis2, basisdual);
      }
   }
   printTables(numSizes);
}

void printTables(int64_t numSizes) {
   int64_t d;
   std::cout << "Total time: " << (double) (clock() - tmpTotal) / (CLOCKS_PER_SEC) << " seconds.\n\n";
   std::cout << "Timings for different methods, in basic clock units (microseconds) \n";
   std::cout << "Dimension:    ";
   for (d = 0; d < numSizes; d++)
      std::cout << std::setw(8) << dimensions[d] << "  ";
   std::cout << "\n\n";
   for (int meth = 0; meth < numMeth; meth++) {
      std::cout << names[meth] << " ";
      for (d = 0; d < numSizes; d++)
         std::cout << std::setw(9) << timer[meth][d] << " ";
      std::cout << "\n";
   }
   std::cout << "\n";
   std::cout << "Sums of square lengths of shortest basis vectors";
   std::cout << " (must be the same for all flexible types): \n";
   std::cout << "Dimension:";
   for (d = 0; d < numSizes; d++)
      std::cout << std::setw(13) << dimensions[d] << " ";
   std::cout << "\n\n";
   for (int meth = 0; meth < numMeth; meth++) {
      if (sumSq[meth][0] > 0.0) {
         std::cout << names[meth] << "  ";
         for (d = 0; d < numSizes; d++)
            std::cout << std::setw(13) << sumSq[meth][d] << " ";
         std::cout << "\n";
      }
   }
   std::cout << "\n";
}

int main() {
   // After changing the value of `mm` (the modulus) one must recompile.
   // Here, `Int` and `Real` are not yet defined, they will be passed as template parameters.
   //int64_t m(1048573);  // Prime modulus near 2^{20}
   //NTL::ZZ mm(1048573);  // Prime modulus near 2^{20}
   // The following values of `mm` work only with ZZ.
   // NTL::ZZ mm(1073741827);  // Prime modulus near 2^{30}
   NTL::ZZ mm(1099511627791);  // Prime modulus near 2^{40}
   // NTL::ZZ mm(1125899906842597);  // Prime modulus near 2^{50}
   int64_t numSizes = 8;
   int64_t numRep = 1000;   // Number of replications (multipliers) for each case.

   // Here we can test with any combination of types.
   // testLoop<int64_t, double>(conv<int64_t>(mm), numSizes, numRep);   // This one works only for the smaller m.
   testLoop<NTL::ZZ, double>(mm, numSizes, numRep);
   testLoop<NTL::ZZ, xdouble>(mm, numSizes, numRep);
   testLoop<NTL::ZZ, quad_float>(mm, numSizes, numRep);
   testLoop<NTL::ZZ, NTL::RR>(mm, numSizes, numRep);
   return 0;
}

