
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
 * with all five combinations of types. See the Lattice Tester guide for more explanations.
 */

using namespace NTL;
using namespace LatticeTester;

// The types Int and Real are not yet defined here.
// They are defined in the main via template parameters when we call the functions.
const int64_t numSizes = 5; // Number of matrix sizes (choices of dimensions).
const int64_t dimensions[numSizes] = { 4, 6, 10, 20, 30 };
const int64_t numMeth = 12;    // Number of methods to test, and their names.
std::string names[numMeth] = { "LLL5         ", "LLL8         ", "LLL99        ", "LLL99999     ",
      "LLL99999-pnew ", "UppTri       ", "mDualUT      ", "LLL5-dual    ", "LLL8-dual    ",
      "LLL99-dual   ", "LLL99999-dual", "LLL99999-dnew" };
// We use `ctime` directly for the timings, to minimize overhead.
clock_t totalTime = clock(); // Global timer for total time.
clock_t timer[numMeth][numSizes];
double sumSq[numMeth][numSizes];  // Sums of square lengths.

void printResults();  // Must be declared, because it has no parameters.

// This function applies LLL to `basis` in `dim` dimensions.
// It also updates the cumulative times and sums of square lengths.
template<typename Int, typename Real>
void LLLTest(IntMat &basis, int64_t d, int64_t meth, double delta) {
   int64_t dim = dimensions[d];
   NTL::Vec<Real> sqlen; // Cannot be global variable because it depends on Real.
   sqlen.SetLength(1);
   clock_t tmp = clock();
   // Here we apply LLL to the basis. Since we already have a basis,
   // we could equivalently use `redLLL` from the file `ReducerStatic`.
   LLLConstruction0(basis, delta, dim, dim, &sqlen);
   timer[meth][d] += clock() - tmp;
   sumSq[meth][d] += conv<double>(sqlen[0]);
}

// Runs a speed test for dim = dimensions[d], with given basis matrices.
// Only basis1 needs to be initialized; basis2 and basisdual are used only for copy.
template<typename Int, typename Real>
void transformBases(Int m, int64_t d, int64_t dim, IntMat &basis1, IntMat &basis2,
      IntMat &basisdual) {
   CopyPartMat<IntMat>(basis2, basis1, dim, dim);  // Copy basis1 to basis2.
   clock_t tmp;

   // We apply LLL to basis2 with different values of `delta`, incrementally.
   // We start with delta=0.5, then continue with 0.9, then with 0.99999.
   LLLTest<Int, Real>(basis2, d, 0, 0.5);
   // We continue the LLL process with larger values of `delta`.
   LLLTest<Int, Real>(basis2, d, 1, 0.8);
   LLLTest<Int, Real>(basis2, d, 2, 0.99);
   LLLTest<Int, Real>(basis2, d, 3, 0.99999);
   // Here we restart LLL from the initial triangular basis, with delta=0.99999.
   CopyPartMat(basis2, basis1, dim, dim);  // Copy basis1 to basis2.
   LLLTest<Int, Real>(basis2, d, 4, 0.99999);

   // We now construct an upper-triangular basis from basis2 into basis1.
   tmp = clock();
   upperTriangularBasis(basis2, basis1, m, dim, dim);
   timer[5][d] += clock() - tmp;
   // We compute an m-dual basis to basis1.
   tmp = clock();
   mDualUpperTriangular(basis1, basisdual, m, dim);
   timer[6][d] += clock() - tmp;

   // We apply LLL to this m-dual basis, with delta = 0.5, 0.8, etc.
   LLLTest<Int, Real>(basisdual, d, 7, 0.5);
   LLLTest<Int, Real>(basisdual, d, 8, 0.8);
   LLLTest<Int, Real>(basisdual, d, 9, 0.99);
   LLLTest<Int, Real>(basisdual, d, 10, 0.99999);
   // Restart anew with delta = 0.99999.
   mDualUpperTriangular(basis1, basisdual, m, dim);
   LLLTest<Int, Real>(basisdual, d, 11, 0.99999);
}

// In this testing loop, new `Rank1Lattice` objects are created
// and the  `IntMat` matrices are resized inside the loop.
// The timings turn out to be about the same.
/*
 template<typename Int, typename IntMat, typename Real>
 void testLoopResize(Int mm, int64_t numRep) {
 int64_t d, dim;
 Int m = conv<Int>(mm);
 Int a;       // The LCG multiplier
 IntMat basis1, basis2, basisdual;
 std::cout << "Results for `testLoopResize` (many objects are created or resized)\n";
 for (d = 0; d < numSizes; d++)      // Reset timers and sums.
 for (int64_t meth = 0; meth < numMeth; meth++) {
 timer[meth][d] = 0;
 sumSq[meth][d] = 0.0;  // NTL::conv<Real>(0.0);
 }
 totalTime = clock(); // Global timer for total time.
 for (int64_t r = 0; r < numRep; r++) {
 a = (m / 5 + 17 * r) % m;   // The multiplier we use for this rep.
 for (d = 0; d < numSizes; d++) {  // Each matrix size
 dim = dimensions[d]; // The corresponding dimension.
 basis1.SetDims(dim, dim); // Will be initial triangular basis.
 basis2.SetDims(dim, dim); // Will be LLL-reduced basis.
 basisdual.SetDims(dim, dim);  // m-dual basis.

 // *** The following does not work well with LLL_FPInt.h (when Int == int64_t).      *******
 Rank1Lattice<Int, Real> korlat(m, a, dim);  // Create a new one.
 // Rank1Lattice<Int, Real> *korlat = new Rank1Lattice<Int, Real> (m, a, dim);  // Create a new one.

 // std::cout << "Just created a new korlat \n";
 korlat.buildBasis(dim);
 basis1 = korlat.getBasis();
 // std::cout << " Basis B = \n" << basis1 << "\n";
 transformBases<Int, IntMat, Real>(m, d, dim, basis1, basis2, basisdual);
 // delete &korlat;
 }
 }
 printResults();
 basis1.kill();  // Since we create objects repeatedly,
 basis2.kill();  // it is a good idea to release the memory when we are done.
 basisdual.kill();
 }
 */

// Testing loop. The `IntMat` and `Rank1Lattice` objects are created only once.
template<typename Int, typename Real>
void testLoop(Int mm, int64_t numRep) {
   std::string stringTypes;  // To print the selected flexible types.
   strTypes<Int, Real>(stringTypes);  // Functions from FlexTypes
   std::cout << "****************************************************\n";
   std::cout << "Types: " << stringTypes << "\n\n";
   std::cout << "TestBasisConstructionSpeed with m = " << mm << "\n";
   std::cout << "Number of replications (different multipliers a): " << numRep << "\n\n";

   int64_t d, dim;  // Index of dimension.
   Int m = conv<Int>(mm);
   Int a;
   IntMat basis1, basis2, basisdual;
   int64_t maxdim = dimensions[numSizes - 1];   // Maximum dimension
   basis1.SetDims(maxdim, maxdim); // Will be initial triangular basis.
   basis2.SetDims(maxdim, maxdim); // Will be LLL-reduced basis.
   basisdual.SetDims(maxdim, maxdim);  // m-dual basis.
   // We create a single Korobov lattice object.
   Rank1Lattice<Int, Real> korlat(m, maxdim);

   for (d = 0; d < numSizes; d++)   // Reset accumulators.
      for (int64_t meth = 0; meth < numMeth; meth++) {
         timer[meth][d] = 0;
         sumSq[meth][d] = 0.0;
      }
   totalTime = clock(); // Global timer for total time.
   for (int64_t r = 0; r < numRep; r++) {
      a = (m / 5 + 17 * r) % m;   // The multiplier we use for this rep.
      korlat.seta(a);
      for (d = 0; d < numSizes; d++) {  // Each matrix size
         dim = dimensions[d]; // The corresponding dimension.
         korlat.buildBasis(dim);
         CopyPartMat<IntMat>(basis1, korlat.getBasis(), dim, dim); // Triangular basis.
         transformBases<Int, Real>(m, d, dim, basis1, basis2, basisdual);
      }
   }
   printResults();
   basis1.kill();
   basis2.kill();
   basisdual.kill();
}

void printResults() {
   int64_t d;
   std::cout << "Timings for different methods, in basic clock units (microseconds) \n\n";
   std::cout << " dim:    ";
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
   std::cout << "Sums of square lengths of shortest basis vector";
   std::cout << " (must be the same for all flexible types):\n";
   std::cout << " dim:    ";
   for (d = 0; d < numSizes; d++)
      std::cout << std::setw(13) << dimensions[d] << "  ";
   std::cout << "\n\n";
   for (int meth = 0; meth < numMeth; meth++) {
      if (sumSq[meth][0] > 0.0) {
         std::cout << names[meth] << "  ";
         for (d = 0; d < numSizes; d++)
            std::cout << std::setw(14) << sumSq[meth][d] << " ";
         std::cout << "\n";
      }
   }
   std::cout << "\n";
   std::cout << "Total time: " << (double) (clock() - totalTime) / (CLOCKS_PER_SEC)
         << " seconds\n\n\n";
}

int main() {

   // Here, `Int` and `Real` are not yet defined, they will be passed as template parameters.
   int64_t m(1048573);  // Prime modulus near 2^{20}
   NTL::ZZ mm(1048573);  // Prime modulus near 2^{20}
   // The following values of `mm` work only with ZZ.
   // NTL::ZZ mm(1073741827);  // Prime modulus near 2^{30}
   // NTL::ZZ mm(1099511627791);  // Prime modulus near 2^{40}
   // NTL::ZZ mm(1125899906842597);  // Prime modulus near 2^{50}
   int64_t numRep = 1000;   // Number of replications (multipliers) for each case.

   // Here we can test with any combination of types.
   testLoop<int64_t, double>(m, numRep);
   testLoop<NTL::ZZ, double>(mm, numRep);
   testLoop<NTL::ZZ, xdouble>(mm, numRep);
   testLoop<NTL::ZZ, quad_float>(mm, numRep);
   testLoop<NTL::ZZ, NTL::RR>(mm, numRep);
   return 0;
}

