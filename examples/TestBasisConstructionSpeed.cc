
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
 * One experiment concerns the construction of triangular bases (TRIE) and the other
 * concerns LLL reduction (LLLE).
 */

using namespace NTL;
using namespace LatticeTester;

enum ExperType { TRIE, LLLE };  // Two sets of experiments: triangular basis vs LLL.

// The types Int and Real are not yet defined here.
// They are passed as template parameters from the `main`.
const int64_t maxNumSizes = 5; // Number of matrix sizes (choices of dimensions).
const int64_t dimensions[maxNumSizes] = { 5, 10, 20, 30, 40 };
const int64_t maxNumMeth = 12;
int64_t numMeth;
std::string names[maxNumMeth];

const int64_t numMethTri = 9;    // Number of methods to test for TRIE, and their names.
std::string namesTri[numMethTri] = { "LLL 0.5   ", "LowerTri  ", "LowerTri  ", "UpperTri  ", "UpperTri2 ",
     "mDualUp   ", "LowerTriD ", "UpperTriD ", "UpperTriD2" };

const int64_t numMethLLL = 12;    // Number of methods to test for LLL, and their names.
std::string namesLLL[numMethLLL] = { "LLL5         ", "LLL8         ", "LLL99        ", "LLL99999     ",
      "LLL99999-pnew ", "UppTri       ", "mDualUT      ", "LLL5-dual    ", "LLL8-dual    ",
      "LLL99-dual   ", "LLL99999-dual", "LLL99999-dnew" };

// We use `ctime` directly for the timings, to minimize overhead.
clock_t totalTime = clock();            // Global timer for total time.
clock_t timer[maxNumMeth][maxNumSizes]; // Collects timings for each case.
double sumSq[maxNumMeth][maxNumSizes];  // Sums of square lengths.


//  void printResults(int64_t numSizes, ExperType exper);  // Must be declared first, because it has no template parameters.

// This function applies LLL to `basis` in `dim` dimensions.
// It also updates the cumulative times and sums of square lengths.
template<typename Int, typename Real>
void LLLTest(IntMat &basis, int64_t d, int64_t meth, double delta) {
   int64_t dim = dimensions[d];
   NTL::Vec<Real> sqlen; // Cannot be global, because it depends on Real.
   sqlen.SetLength(1);
   clock_t tmp = clock();
   // Here we apply LLL to the basis.
   // We could equivalently use `redLLL` from the file `ReducerStatic`.
   LLLConstruction0(basis, delta, dim, dim, &sqlen);
   timer[meth][d] += clock() - tmp;
   sumSq[meth][d] += conv<double>(sqlen[0]);
}

// Runs a speed test for dim = dimensions[d], with given basis matrices.
// Only basis1 needs to be initialized; basis2 and basisdual are used only for copy.
template<typename Int, typename Real>
void transformBasesLLL(const Int &m, int64_t d, int64_t dim, IntMat &basis1, IntMat &basis2,
      IntMat &basisdual) {
   clock_t tmp;
   CopyPartMat<IntMat>(basis2, basis1, dim, dim);  // Copy basis1 to basis2.

   // We apply LLL to basis2 with different values of `delta`, incrementally.
   // We start with delta=0.5, then continue with 0.8, then with 0.99, etc.
   LLLTest<Int, Real>(basis2, d, 0, 0.5);
   LLLTest<Int, Real>(basis2, d, 1, 0.8);
   LLLTest<Int, Real>(basis2, d, 2, 0.99);
   LLLTest<Int, Real>(basis2, d, 3, 0.99999);
   // Now we restart LLL from the initial triangular basis, with delta=0.99999.
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
   // We restart anew with delta = 0.99999.
   mDualUpperTriangular(basis1, basisdual, m, dim);
   LLLTest<Int, Real>(basisdual, d, 11, 0.99999);
}


// Runs a speed test for dim = dimensions[d], for triangular basis constructions,
// from given basis matrices
// Only basis1 needs to be initialized; basis2 and basisdual are used only for copy.
template<typename Int, typename Real>
void triangularBases(const Int &m, int64_t d, int64_t dim, IntMat &basis1, IntMat &basis2,
      IntMat &basisdual) {
   clock_t tmp;
   NTL::Vec<Real> sqlen; // Cannot be global variable because it depends on Real.
   sqlen.SetLength(1);
   CopyPartMat<IntMat>(basis2, basis1, dim, dim);  // Copy basis1 to basis2.

   // We first apply LLL to basis2.
   tmp = clock();
   LLLConstruction0<Int, Real>(basis2, 0.5, dim, dim, &sqlen);
   timer[0][d] += clock() - tmp;
   sumSq[0][d] += conv<double>(sqlen[0]);

   tmp = clock();
   lowerTriangularBasis<Int>(basis2, basis1, m, dim, dim);
   timer[1][d] += clock() - tmp;

   LLLConstruction0<Int, Real>(basis1, 0.5, dim, dim);
   tmp = clock();
   lowerTriangularBasis<Int>(basis1, basis2, m, dim, dim);
   timer[2][d] += clock() - tmp;

   LLLConstruction0<Int, Real>(basis2, 0.5, dim, dim, &sqlen);
   tmp = clock();
   upperTriangularBasis<Int>(basis2, basis1, m, dim, dim);
   timer[3][d] += clock() - tmp;
   sumSq[3][d] += conv<double>(sqlen[0]);

   tmp = clock();
   upperTriangularBasis<Int>(basis1, basis2, m, dim, dim);
   timer[4][d] += clock() - tmp;

   // We compute an m-dual basis to basis2.
   tmp = clock();
   mDualUpperTriangular(basis2, basisdual, m, dim);
   timer[5][d] += clock() - tmp;

   LLLConstruction0<Int, Real>(basisdual, 0.5, dim, dim, &sqlen);
   tmp = clock();
   lowerTriangularBasis<Int>(basisdual, basis1, m, dim, dim);
   timer[6][d] += clock() - tmp;
   sumSq[6][d] += conv<double>(sqlen[0]);

   LLLConstruction0<Int, Real>(basis1, 0.5, dim, dim, &sqlen);
   sumSq[7][d] += conv<double>(sqlen[0]);
   // std::cout << "Basis after LLL with delta=0.5: \n" << basis1 << "\n";
   // std::cout << "Square length of first basis vector: " << sqlen[0] << "\n\n";
   tmp = clock();
   upperTriangularBasis<Int>(basis1, basis2, m, dim, dim);
   timer[7][d] += clock() - tmp;

   tmp = clock();
   upperTriangularBasis<Int>(basis2, basis1, m, dim, dim);
   timer[8][d] += clock() - tmp;
}


// Testing loop. The `IntMat` and `Rank1Lattice` objects are created only once.
template<typename Int, typename Real>
void testLoop(Int &mm, int64_t numSizes, int64_t numRep, ExperType exper) {
   std::string stringTypes;  // To print the selected flexible types.
   strTypes<Int, Real>(stringTypes);  // Functions from FlexTypes
   std::cout << "****************************************************\n";
   std::cout << "Types: " << stringTypes << "\n\n";
   std::cout << "TestBasisConstructionSpeed with m = " << mm << "\n";
   std::cout << "Number of replications (different multipliers a): " << numRep << "\n";
   std::cout << "Type of Experiment: " << exper << "\n\n";

   if (exper == TRIE) {
      numMeth = numMethTri;
      for (int64_t meth = 0; meth < numMeth; meth++)
         names[meth] = namesTri[meth];
   }
   if (exper == LLLE) {
      numMeth = numMethLLL;
      for (int64_t meth = 0; meth < numMeth; meth++)
          names[meth] = namesLLL[meth];
   }
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
         // Resizing the matrices as below slows down the program by about 2-3 percent.
         // basis1.SetDims(dim, dim); // Will be initial triangular basis.
         // basis2.SetDims(dim, dim); // Will be LLL-reduced basis.
         // basisdual.SetDims(dim, dim); // Will be LLL-reduced basis.
         CopyPartMat<IntMat>(basis1, korlat.getBasis(), dim, dim); // Copy korlat basis to basis1.
         if (exper == TRIE)
            triangularBases<Int, Real>(m, d, dim, basis1, basis2, basisdual);
         if (exper == LLLE)
            transformBasesLLL<Int, Real>(m, d, dim, basis1, basis2, basisdual);
      }
   }
   printResults<Int, Real>(numSizes, exper);
}

template<typename Int, typename Real>
void printResults(int64_t numSizes, ExperType exper) {
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
   std::cout << " (must be the same for all flexible types) \n";
   if (exper == TRIE)
      std::cout << "For the `TRIE` case, the sums are for the shortest vectors after the last LLL \n";
   std::cout << " dim:    ";
   for (d = 0; d < numSizes; d++)
      std::cout << std::setw(13) << dimensions[d] << "  ";
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
   int64_t numSizes = 5;
   int64_t numRep = 1000;   // Number of replications (multipliers) for each case.

   ExperType exper(TRIE);
   // ExperType exper(LLLE);

   // Here we can test with any combination of types.
   testLoop<int64_t, double>(m, numSizes, numRep, exper);
   testLoop<NTL::ZZ, double>(mm, numSizes, numRep, exper);
   testLoop<NTL::ZZ, xdouble>(mm, numSizes, numRep, exper);
   //testLoop<NTL::ZZ, quad_float>(mm, numSizes, numRep, exper);
   //testLoop<NTL::ZZ, NTL::RR>(mm, numSizes, numRep, exper);
   return 0;
}

