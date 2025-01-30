// File `testBasisTriSmall`

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
 * This program examines what happens when we apply LLL and construct triangular bases
 * for the primal and m-dual lattices associated with an LCG.
 * We want to see for instance when `m e_i` becomes a shortest vector in the primal.
 */
using namespace NTL;
using namespace LatticeTester;

// The types Int and Real are will be passed as template parameters from the `main`.
const int64_t maxNumSizes = 5; // Number of matrix sizes (choices of dimensions).
const int64_t dimensions[maxNumSizes] = { 5, 10, 15, 20, 25 };

// Runs a test for dim = dimensions[d], for triangular basis constructions.
// Only basis0 needs to be initialized; the other matrices are used only for copy.
template<typename Int, typename Real>
void triangularBasesTrace(const Int &m, int64_t d, int64_t dim, IntMat &basis0,
      IntMat &basis1, IntMat &basis2, IntMat &basisdual) {
   NTL::Vec<Real> sqlen; // Cannot be global variable because it depends on Real.
   sqlen.SetLength(1);

   // We first apply LLL to basis0.
   std::cout << "*************************************************************\n";
   std::cout << "basis primal before LLL5:\n" << basis0 << "\n";
   LLLConstruction0<Int, Real>(basis0, 0.5, dim, dim, &sqlen);
   std::cout << "basis primal after LLL5:\n" << basis0 << "\n";
   std::cout << "shortest vector length: " <<  sqrt(conv<double>(sqlen[0])) << "\n\n";

   CopyPartMat<IntMat>(basis1, basis0, dim, dim);  // Copy basis0 to basis1.
   lowerTriangularBasis<Int>(basis2, basis1, m, dim, dim);
   mDualLowerTriangular(basisdual, basis2, m, dim);
   std::cout << "upper-triangular basisdual before LLL5:\n" << basisdual << "\n";
   LLLConstruction0<Int, Real>(basisdual, 0.5, dim, dim, &sqlen);
   std::cout << "basisdual after LLL5:\n" << basisdual << "\n";
   std::cout << "shortest vector length in m-dual: " << sqrt(conv<double>(sqlen[0])) << "\n\n";
}

// Testing loop. The `IntMat` and `Rank1Lattice` objects are created only once.
template<typename Int, typename Real>
void testLoop(const Int &mm, int64_t numSizes, int64_t numRep) {
   std::string stringTypes;  // To print the selected flexible types.
   strTypes<Int, Real>(stringTypes);  // Functions from FlexTypes
   std::cout << "**************************************************************\n";
   std::cout << "TestBasisTriSmall with m = " << mm << "\n";
   std::cout << "Types: " << stringTypes << "\n";
   std::cout << "Number of replications (different multipliers a): " << numRep << "\n";
   int64_t d, dim;  // Index of dimension, and dimension.
   Int m = conv<Int>(mm);
   Int a;
   IntMat basis0, basis1, basis2, basisdual;
   int64_t maxdim = dimensions[numSizes - 1];   // Maximum dimension
   Rank1Lattice<Int, Real> korlat(m, maxdim);
   for (int64_t r = 0; r < numRep; r++) {
      a = (m / 5 + 17 * r) % m;   // The multiplier `a` used for this rep.
      korlat.seta(a);
      for (d = 0; d < numSizes; d++) {  // Each matrix size
         dim = dimensions[d]; // The corresponding dimension.
         basisdual.SetDims(dim, dim);  // m-dual basis.
         basis0.SetDims(dim, dim);  // basis.
         basis1.SetDims(dim, dim);  // basis.
         basis2.SetDims(dim, dim);  // basis.
         korlat.buildBasis(dim);
         CopyPartMat<IntMat>(basis0, korlat.getBasis(), dim, dim); // Copy korlat basis to basis1.
         triangularBasesTrace<Int, Real>(m, d, dim, basis0, basis1, basis2, basisdual);
      }
   }
}

int main() {

   // Here, `Int` and `Real` are not yet defined, they will be passed as template parameters.
   NTL::ZZ mm(1021);  // Prime modulus near 2^{10}
   // int64_t m(1048573);  // Prime modulus near 2^{20}
   // NTL::ZZ mm(1048573);  // Prime modulus near 2^{20}
   // The following values of `mm` work only with ZZ.
   // NTL::ZZ mm(1073741827);  // Prime modulus near 2^{30}
   //int64_t numSizes = 8;
   //int64_t numRep = 1000;   // Number of replications (multipliers) for each case.
   int64_t numSizes = 5;
   int64_t numRep = 1;   // Number of replications (multipliers) for each case.

   // Here we can test with any combination of types.
   //testLoop<int64_t, double>(conv<int64_t>(mm), numSizes, numRep);  // This one works only for the smaller m.
   testLoop<NTL::ZZ, double>(mm, numSizes, numRep);
   //testLoop<NTL::ZZ, xdouble>(mm, numSizes, numRep);
   //testLoop<NTL::ZZ, quad_float>(mm, numSizes, numRep);
   // testLoop<NTL::ZZ, NTL::RR>(mm, numSizes, numRep);
   return 0;
}

