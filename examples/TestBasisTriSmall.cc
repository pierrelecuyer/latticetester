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
 * It shows for instance when `m e_i` becomes a shortest vector in the primal,
 * and also some particular property of the m-dual basis after LLL.
 */
using namespace NTL;
using namespace LatticeTester;

// The types Int and Real are will be passed as template parameters from the `main`.
const int64_t maxNumSizes = 7; // Number of matrix sizes (choices of dimensions).
const int64_t dimensions[maxNumSizes] = { 5, 10, 15, 20, 25, 30, 35 };

// Runs a test for dim = dimensions[d], for triangular basis constructions.
// Only basis0 needs to be initialized; the other matrices are used only for copy.
template<typename Int, typename Real>
void triangularBasesTrace(const Int &m, int64_t d, int64_t dim, double delta, IntMat &basis0,
IntMat &basis1, IntMat &basis2, IntMat &basisdual) {
   NTL::Vec<Real> sqlen; // Cannot be global variable because it depends on Real.
   sqlen.SetLength(1);

   // We first apply LLL to basis0.
   std::cout << "*************************************************************\n";
   std::cout << "Number of dimensions: " << dim << "\n\n";
   std::cout << "primal basis before LLL:\n" << basis0 << "\n";
   LLLConstruction0<Int, Real>(basis0, delta, dim, dim, &sqlen);
   std::cout << "primal basis after LLL:\n" << basis0 << "\n";
   std::cout << "shortest vector length: " << sqrt(conv<double>(sqlen[0])) << "\n\n";

   CopyPartMat<IntMat>(basis1, basis0, dim, dim);  // Copy basis0 to basis1.
   lowerTriangularBasis<Int>(basis2, basis1, m, dim, dim);
   std::cout << "lower-triangular primal basis before LLL:\n" << basis2 << "\n";
   mDualLowerTriangular(basisdual, basis2, m, dim);
   std::cout << "upper-triangular basisdual before LLL:\n" << basisdual << "\n";
   LLLConstruction0<Int, Real>(basisdual, delta, dim, dim, &sqlen);
   std::cout << "basisdual after LLL:\n" << basisdual << "\n";
   std::cout << "shortest vector length in m-dual: " << sqrt(conv<double>(sqlen[0])) << "\n\n";

   CopyPartMat<IntMat>(basis1, basis0, dim, dim);  // Copy basis0 to basis1.
   upperTriangularBasis<Int>(basis2, basis1, m, dim, dim);
   std::cout << "upper-triangular primal basis before LLL:\n" << basis2 << "\n";
   mDualUpperTriangular(basisdual, basis2, m, dim);
   std::cout << "lower-triangular basisdual before LLL:\n" << basisdual << "\n";
   LLLConstruction0<Int, Real>(basisdual, delta, dim, dim, &sqlen);
   std::cout << "basisdual after LLL:\n" << basisdual << "\n";
   std::cout << "shortest vector length in m-dual: " << sqrt(conv<double>(sqlen[0])) << "\n\n";
}

// Testing loop. The `IntMat` and `Rank1Lattice` objects are created only once.
template<typename Int, typename Real>
void testLoop(const Int &m, Int &a, int64_t numSizes, double delta) {
   std::string stringTypes;  // To print the selected flexible types.
   strTypes<Int, Real>(stringTypes);  // Functions from FlexTypes
   std::cout << "**************************************************************\n";
   std::cout << "TestBasisTriSmall with m = " << m << ",  a = " << a << "\n";
   std::cout << "Types: " << stringTypes << "\n";
   std::cout << "LLL wih delta = " << delta << "\n";
   int64_t d, dim;  // Index of dimension, and dimension.
   IntMat basis0, basis1, basis2, basisdual;
   int64_t maxdim = dimensions[numSizes - 1];   // Maximum dimension
   Rank1Lattice<Int, Real> korlat(m, maxdim);
   //  a = (m / 5 + 17 * r) % m;   // The multiplier `a` used for this rep.
   korlat.seta(a);
   for (d = 0; d < numSizes; d++) {  // Each matrix size
      dim = dimensions[d]; // The corresponding dimension.
      basisdual.SetDims(dim, dim);  // m-dual basis.
      basis0.SetDims(dim, dim);  // basis.
      basis1.SetDims(dim, dim);  // basis.
      basis2.SetDims(dim, dim);  // basis.
      korlat.buildBasis(dim);
      CopyPartMat<IntMat>(basis0, korlat.getBasis(), dim, dim); // Copy korlat basis to basis1.
      triangularBasesTrace<Int, Real>(m, d, dim, delta, basis0, basis1, basis2, basisdual);
   }
}

int main() {

// Here, `Int` and `Real` are not yet defined, they will be passed as template parameters.
NTL::ZZ mm(1021);  // Prime modulus near 2^{10}
//TL::ZZ mm(1048573);  // Prime modulus near 2^{20}
//NTL::ZZ mm(1073741827);  // Prime modulus near 2^{30}
// NTL::ZZ a(73245663);
NTL::ZZ a(73);
int64_t numSizes = 7;
double delta = 0.5;

testLoop<NTL::ZZ, double>(mm, a, numSizes, delta);
return 0;
}

