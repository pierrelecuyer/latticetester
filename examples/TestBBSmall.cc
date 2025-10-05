/**
 * Here we take a lattice obtained from a small LCG and we follow what happens when
 * we run `shortestVector` from `ReducerBB` to find a shortest vector.
 * We consider the L1 and L2 norms and compare he `CHOLESKY` and `TRIANGULAR`
 * decompositions, for both the primal and $m$-dual lattices.
 * For Cholesky, we first perform a weak LLL reduction before applying the BB.
 */

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#define TYPES_CODE  ZD     // Int = ZZ, Real = double

#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include "latticetester/FlexTypes.h"
#include "latticetester/Util.h"
#include "latticetester/Rank1Lattice.h"
#include "latticetester/ReducerStatic.h"
#include "latticetester/ReducerBB.h"

using namespace LatticeTester;

template<typename Int, typename Real>
void findShortest (Rank1Lattice<Int, Real> &korlat, ReducerBB<Int, Real> &red,
      NormType norm, DecompTypeBB decomp, bool inDual,
      long dim, double deltaLLL) {
   korlat.setNormType(norm);     // Select the norm.
   red.setDecompTypeBB(decomp);  // Select the decomposition type.
   if (inDual) {
      korlat.buildDualBasis(dim);  // Rebuild the dual basis (only) anew.
      korlat.dualize();
   } else {
      korlat.buildBasis(dim);  // Rebuild the primal basis anew.
      // std::cout << "After build basis = \n" << korlat.getBasis() << "\n";
   }
   std::cout << "\n**********************************\nPerforming Reduction: \n";
   if (inDual) std::cout << "DUAL lattice,  ";
   else std::cout << "PRIMAL lattice,  ";
   std::cout << "Norm: " << toStringNorm(norm) << ",  ";
   std::cout << "Decomposition: " << toStringDecomp(decomp) << ".\n\n";
   std::cout << "Before pre-reduction, initial basis = \n" << korlat.getBasis() << "\n";

   Real minSqlen = Real(0.0);
   if (deltaLLL > 0.0) minSqlen = redLLL<Int, Real>(korlat.getBasis(), deltaLLL, dim);
   std::cout << "After LLL reduction, squared L2 norm = " << conv<double>(minSqlen) << ",  basis = \n" << korlat.getBasis() << "\n";
   red.shortestVector(korlat);  // BB is applied here.  ***

   //std::cout << "After BB, square length of shortest vector: " << red.getMinLength2() << "\n";
   std::cout << "  length of shortest vector: " << sqrt(red.getMinLength2()) << "\n";
   //std::cout << "Number of calls to BB procedure `tryZ`: " << red.getCountNodes() << "\n\n";
}

int main() {
   // NTL::ZZ m(1021);  // Prime modulus near 2^{10}
   // Int a(73);
   NTL::ZZ m(1048573);  // Prime modulus near 2^{20}
   Int a(29873);
   // NTL::ZZ m(1073741827);  // Prime modulus near 2^{30}
   long dim = 6;

   std::cout << "=========================================================\n";
   std::string stringTypes;  // To print the selected flexible types.
   strTypes<Int, Real>(stringTypes);  // Functions from FlexTypes
   std::cout << "TestReducerBBSmall,  types: " << stringTypes << "\n";
   std::cout << "Modulo m = " << m;
   std::cout << ",   Multiplier a = " << a << "\n";

   Rank1Lattice<Int, Real> korlat(m, dim, L2NORM); // We use single lattice object.
   korlat.seta(a);
   ReducerBB<Int, Real> red(korlat);   // Single ReducerBB with internal lattice `korlat`.
   red.setVerbosity(2);

   bool inDual = false;  // Primal basis
   findShortest<Int, Real>(korlat, red, L2NORM, CHOLESKY, inDual, dim, 0.99);
   findShortest<Int, Real>(korlat, red, L2NORM, TRIANGULAR, inDual, dim, 0.99);
   findShortest<Int, Real>(korlat, red, L1NORM, CHOLESKY, inDual, dim, 0.99);
   findShortest<Int, Real>(korlat, red, L1NORM, TRIANGULAR, inDual, dim, 0.99);
   inDual = true;   // Dual basis
   findShortest(korlat, red, L2NORM, CHOLESKY, inDual, dim, 0.99);
   findShortest(korlat, red, L2NORM, TRIANGULAR, inDual, dim, 0.99);
   findShortest(korlat, red, L1NORM, CHOLESKY, inDual, dim, 0.99);
   findShortest(korlat, red, L1NORM, TRIANGULAR, inDual, dim, 0.99);
}
