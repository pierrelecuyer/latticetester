// File `testBasisConstructSpeedTri`

#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/xdouble.h>
#include <NTL/ZZ.h>
#include <NTL/mat_ZZ.h>
#include "latticetester/FlexTypes.h"
#include "latticetester/Util.h"
#include "latticetester/Rank1Lattice.h"
#include "latticetester/BasisConstruction.h"
// #include "latticetester/ReducerStatic.h"
// #include "latticetester/ReducerBB.h"

/**
 * This example concerns mostly the construction of triangular bases and their m-duals.
 * It makes speed comparisons with the `BasisConstruction` functions,
 * with different combinations of types. See the Lattice Tester guide for more explanations.
 * The `main` must be changed and recompiled to change the value of `m`.
 */
using namespace NTL;
using namespace LatticeTester;

// The types Int and Real are will be passed as template parameters from the `main`.
const int64_t maxNumSizes = 8; // Number of matrix sizes (choices of dimensions).
const int64_t dimensions[maxNumSizes] = { 5, 10, 20, 30, 40, 50, 60, 70 };
const int64_t numMeth = 14;    // Number of methods and their names.
std::string names[numMeth] = { "LLLPrimal    ", "mDualBasis   ", "LLLDualmDual " , "LowTriP      ", "mDualLow     ",
      "LLLDualUT    ", "UppTriP      ", "UppTriP2     ", "mDualUp      ", "LLLDualLT    ",
      "LowTriDual   ", "UppTriDual   ", "UppTriDualOld", "UppTriDual2  " };

// We use `ctime` directly for the timings, to minimize overhead.
clock_t tmpTotal = 0;            // Global timer for total time.
clock_t timer[numMeth][maxNumSizes]; // Collects timings for each case.
double sumSq[numMeth][maxNumSizes];  // Sums of square lengths, for a few cases.

// Declaration required.
void printTable(int64_t numSizes, double delta);

// Runs a speed test for dim = dimensions[d], for triangular basis constructions.
// Only basis0 needs to be initialized; the other matrices are used only for copy.
template<typename Int, typename Real>
void triangularBases(const Int &m, int64_t d, int64_t dim, double delta, IntMat &basis0,
      IntMat &basis1, IntMat &basis2, IntMat &basisdual) {  // , ReducerBB<Int, Real> &red, IntLattice<Int, Real> &lat2) {

   IntMat basis4, basis5;   // Used for the `mDualBasis` function.
   basis4.SetDims (dim, dim);
   basis5.SetDims (dim, dim);
   clock_t tmp;
   NTL::Vec<Real> sqlen; // Cannot be global variable because it depends on Real.
   sqlen.SetLength(1);

   // We first apply LLL to primal basis0.
   tmp = clock();
   LLLConstruction0<Int, Real>(basis0, delta, dim, dim, &sqlen);
   timer[0][d] += clock() - tmp;
   // std::cout << "basis primal after LLL5:\n" << basis0 << "\n";
   // std::cout << "shortest vector length: " <<  sqrt(conv<double>(sqlen[0])) << "\n\n";

   CopyPartMat<IntMat>(basis4, basis0, dim, dim);  // Copy basis0 to basis1.
   tmp = clock();
   mDualBasis (basis5, basis4, m);
   timer[1][d] += clock() - tmp;

   tmp = clock();
   LLLConstruction0<Int, Real>(basis5, delta, dim, dim, &sqlen);
   timer[2][d] += clock() - tmp;
   // lat2.setBasis(basisdual, dim);
   // red.setIntLattice(lat2);
   // red.shortestVector();
   // sumSq[2][d] += conv<double>(red.getMinLength2());
   sumSq[2][d] += conv<double>(sqlen[0]);

   CopyPartMat<IntMat>(basis1, basis0, dim, dim);  // Copy basis0 to basis1.
   tmp = clock();
   lowerTriangularBasis<Int>(basis2, basis1, m, dim, dim);
   timer[3][d] += clock() - tmp;
   // std::cout << "lower-triangular primal basis:\n" << basis2 << "\n";

   // Compute the m-dual of basis2, which will be upper-triangular.
   tmp = clock();
   mDualLowerTriangular(basisdual, basis2, m, dim);
   timer[4][d] += clock() - tmp;
   // std::cout << "upper-triangular m-dual basis before LLL5:\n" << basisdual << "\n";

   tmp = clock();
   LLLConstruction0<Int, Real>(basisdual, delta, dim, dim, &sqlen);
   timer[5][d] += clock() - tmp;
   sumSq[5][d] += conv<double>(sqlen[0]);

   CopyPartMat<IntMat>(basis1, basis0, dim, dim);  // Copy basis0 to basis1.
   tmp = clock();
   upperTriangularBasis<Int>(basis2, basis1, m, dim, dim);
   timer[6][d] += clock() - tmp;

   tmp = clock();
   upperTriangularBasis<Int>(basis1, basis2, m, dim, dim);
   timer[7][d] += clock() - tmp;

   // We compute an m-dual basis to basis1.
   tmp = clock();
   mDualUpperTriangular(basisdual, basis1, m, dim);
   timer[8][d] += clock() - tmp;

   // Here we start from a lower-triangular m-dual basis.
   tmp = clock();
   LLLConstruction0<Int, Real>(basisdual, delta, dim, dim, &sqlen);
   timer[9][d] += clock() - tmp;
   sumSq[9][d] += conv<double>(sqlen[0]);
   // std::cout << "basisdual after LLL:\n" << basisdual << "\n";
   // std::cout << "shortest vector length: " << sqrt(conv<double>(sqlen[0])) << "\n\n";

   CopyPartMat<IntMat>(basis1, basisdual, dim, dim);  // Copy basisdual to basis1.
   tmp = clock();
   lowerTriangularBasis<Int>(basis2, basis1, m, dim, dim);
   timer[10][d] += clock() - tmp;

   CopyPartMat<IntMat>(basis1, basisdual, dim, dim);  // Copy basisdual to basis1.
   tmp = clock();
   // This one is much slower than the previous (lower-triangular) one!
   upperTriangularBasis<Int>(basis2, basis1, m, dim, dim);
   timer[11][d] += clock() - tmp;

   CopyPartMat<IntMat>(basis1, basisdual, dim, dim);  // Copy basisdual to basis1.
   tmp = clock();
   upperTriangularBasisOld96<IntMat, Int>(basis2, basis1, m, dim, dim);
   timer[12][d] += clock() - tmp;

   tmp = clock();
   upperTriangularBasis<Int>(basis1, basis2, m, dim, dim);
   timer[13][d] += clock() - tmp;
}


// Testing loop. The `IntMat` and `Rank1Lattice` objects are created only once.
template<typename Int, typename Real>
void testLoop(const Int &m, int64_t numSizes, int64_t numRep, double delta) {
   std::string stringTypes;  // To print the selected flexible types.
   strTypes<Int, Real>(stringTypes);  // Functions from FlexTypes
   std::cout << "**************************************************************\n";
   std::cout << "TestBasisConstructSpeedTri with m = " << m << "\n";
   std::cout << "Types: " << stringTypes << "\n";
   std::cout << "LLL with delta = " << delta << "\n";
   std::cout << "Number of replications (different multipliers a): " << numRep << "\n";
   int64_t d, dim;  // Index of dimension, and dimension.
   Int a;
   IntMat basis0, basis1, basis2, basisdual;
   int64_t maxdim = dimensions[numSizes - 1];   // Maximum dimension
   basis0.SetDims(maxdim, maxdim); // Will be initial triangular basis.
   basis1.SetDims(maxdim, maxdim); // Another basis.
   basis2.SetDims(maxdim, maxdim); // Another basis.
   basisdual.SetDims(maxdim, maxdim);  // m-dual basis.
   // We create a single Korobov lattice object.
   Rank1Lattice<Int, Real> korlat(m, maxdim);
   // IntLattice<Int, Real> lat2(m, maxdim);  //
   // ReducerBB<Int, Real> red(lat2);  // In case we want to do the BB to get the true shortest length.

   for (d = 0; d < numSizes; d++)   // Reset accumulators.
      for (int64_t meth = 0; meth < numMeth; meth++) {
         timer[meth][d] = 0;
         sumSq[meth][d] = 0.0;
      }
   tmpTotal = clock();
   for (int64_t r = 0; r < numRep; r++) {
      a = (m / 5 + 17 * r) % m;   // The multiplier `a` used for this rep.
      korlat.seta(a);
      for (d = 0; d < numSizes; d++) {  // Each matrix size
         dim = dimensions[d]; // The corresponding dimension.
         korlat.buildBasis(dim);
         CopyPartMat<IntMat>(basis0, korlat.getBasis(), dim, dim); // Copy korlat basis to basis1.
         triangularBases<Int, Real>(m, d, dim, delta, basis0, basis1, basis2, basisdual);  // , red, lat2);
      }
   }
   printTable(numSizes, delta);
}

void printTable(int64_t numSizes, double delta) {
   int64_t d;
   std::cout << "Total time (includes time for other operations between triangularizations): "
             << (double) (clock() - tmpTotal) / (CLOCKS_PER_SEC) << " seconds.\n\n";
   std::cout << "Dimension:    ";
   for (d = 0; d < numSizes; d++)
      std::cout << std::setw(8) << dimensions[d] << "  ";
   std::cout << "\n\n";
   std::cout << "Timings for the different tasks, in basic clock units (microseconds): \n";
   for (int meth = 0; meth < numMeth; meth++) {
      std::cout << names[meth] << " ";
      for (d = 0; d < numSizes; d++)
         std::cout << std::setw(9) << timer[meth][d] << " ";
      std::cout << "\n";
   }
   std::cout << "\nSum of squares after each of the LLL in dual, with delta = " << delta << " :\n";
   for (int meth = 0; meth < numMeth; meth++) {
      if (sumSq[meth][0] > 0.0) {
         std::cout << names[meth] << " ";
         for (d = 0; d < numSizes; d++)
            std::cout << std::setw(9) << sumSq[meth][d] << " ";
         std::cout << "\n";
      }
   }
   std::cout << "\n";
}


int main() {

   // Here, `Int` and `Real` are not yet defined, they will be passed as template parameters.
   // NTL::ZZ mm(1021);  // Prime modulus near 2^{10}
   // int64_t m(1048573);  // Prime modulus near 2^{20}
   // NTL::ZZ mm(1048573);  // Prime modulus near 2^{20}
   // The following values of `mm` work only with ZZ.
   // NTL::ZZ mm(1073741827);  // Prime modulus near 2^{30}
   NTL::ZZ mm(1099511627791);  // Prime modulus near 2^{40}
   // NTL::ZZ mm(1125899906842597);  // Prime modulus near 2^{50}
   int64_t numSizes = 8;
   int64_t numRep = 1000;   // Number of replications (multipliers) for each case.
   double delta = 0.99999;

   // Here we can test with any combination of types.
   // testLoop<int64_t, double>(conv<int64_t>(mm), numSizes, numRep, delta);  // This one works only for the smaller m.
   testLoop<NTL::ZZ, double>(mm, numSizes, numRep, delta);
   //testLoop<NTL::ZZ, xdouble>(mm, numSizes, numRep, delta);
   //testLoop<NTL::ZZ, quad_float>(mm, numSizes, numRep, delta);
   //testLoop<NTL::ZZ, NTL::RR>(mm, numSizes, numRep, delta);
   return 0;
}

