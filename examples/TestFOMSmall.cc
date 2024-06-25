/**
 * In this example, we compare different ways of computing the shortest vector
 * in a lattice or in its $m$-dual.  For this, we do the BB after some type of
 * pre-reduction of the basis. We compare the total times to do that
 * (pre-red. + BB) for various types of pre-reductions, with lattices that come
 * from Korobov lattice rules in 10 to 40 dimensions, with prime modulus `m`.
 * For the BB, we use the Cholesky decomposition.
 * See the Lattice Tester Guide for more explanations and results.
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
#include "latticetester/FigureOfMeritM.h"
#include "latticetester/FigureOfMeritDualM.h"
#include "latticetester/CoordinateSets.h"
#include "latticetester/WeightsUniform.h"

using namespace LatticeTester;


// Applies some tests to an LCG with modulus m and multiplier a, for projections specified by t.
// A trace is printed on the terminal.
template<typename Int, typename Real>
//void testOneLCG (const Int m, const NTL::vector<Int> a, const NTL::vector<int64_t> t) {
void testOneLCG (Rank1Lattice<Int, Real> korlat, FigureOfMeritM<Int, Real> fomprimal,
       FigureOfMeritDualM<Int, Real> fomdual) {
}

// In this testing loop, we generate `numRep` multipliers `a` and for each one
// we call `testOneLCG`.
template<typename Int, typename Real>
void testLoop (const Int m, const NTL::vector<int64_t> t, long numRep) {
   std::string stringTypes;  // To print the selected flexible types.
   strTypes<Int, Real>(stringTypes);  // Functions from FlexTypes
   std::cout << "****************************************************\n";
   std::cout << "Types: " << stringTypes << "\n";
   std::cout << "We look at FOM properties in primal and m-dual for a small list of LCGs.\n";
   std::cout << "TestFOMSmall with m = " << m << ", vector t = " << t << "\n\n";
   // if (inDual) std::cout << ", in the dual lattice. \n\n";
   // else std::cout << ", in the primal lattice. \n\n";
   std::cout << "-----------------------------------------------------\n";
   int64_t maxdim = t[0];  // Maximum dimension of the lattice
   Rank1Lattice<Int, Real> korlat(m, maxdim); // We use single lattice object.
   ReducerBB<Int, Real> red(korlat);  // Single ReducerBB with internal lattice `korlat`.
   WeightsUniform weights(1.0);
   NormaBestLat normaprimal(log(m), 1, maxdim);  // Factors computed for primal.
   NormaBestLat normadual(-log(m), 1, maxdim);  // Factors computed for dual.
   normadual.computeBounds (-log(m), 1);
   FigureOfMeritM<Int, Real> fomprimal (t, weights, normaprimal, &red, true);
   FigureOfMeritDualM<Int, Real> fomdual (t, weights, normadual, &red, true);
   fomprimal.setVerbosity(2);
   fomdual.setVerbosity(2);
   IntLattice<Int, Real> proj(m, t.size());
   Int a;        // The LCG multiplier
   double merit;
   for (int64_t r = 0; r < numRep; r++) {
      // a = (m / 5 + 13 * r) % m;   // The multiplier we use for this rep.
      a = 33;
      korlat.seta(a);
      std::cout << "a = " << a << "\n";
      korlat.buildBasis(maxdim);
      korlat.buildDualBasis(maxdim);

      merit = fomprimal.computeMerit(korlat, proj);
      std::cout << "FOM value in primal: " << merit << ",  worst projection: " << fomprimal.getMinMeritProj() << "\n\n";
      // merit = fomdual.computeMeritSucc(korlat);
      std::cout << "Primal basis of korlat: \n" << korlat.getBasis() << "\n";
      std::cout << "Primal basis of last projection: \n" << proj.getBasis() << "\n";

      merit = fomdual.computeMerit(korlat, proj);
      std::cout << "FOM value in dual: " << merit << ",  worst projection: " << fomdual.getMinMeritProj() << "\n\n";
      // testOneLCG (korlat, fomprimal, fomdual);
      std::cout << "m-Dual basis of korlat: \n" << korlat.getDualBasis() << "\n";
      std::cout << "m-Dual basis of last projection (1,2): \n" << proj.getDualBasis() << "\n";
    }
}

int main() {

   NTL::vector<int64_t> t(2); // The t-vector for the FOM.
   t[0] = 5;    // We look at successive coordinates in up to t[0] dimensions.
   t[1] = 3;    // Then pairs and triples up to coord. 5.
   //t[2] = 3;

   // Here, Int and Real are not yet defined, they will be passed as template parameters.
   std::cout << "Types: NTL::ZZ, double \n";
   NTL::ZZ m(101);  // Prime modulus near 100
   // NTL::ZZ m(1048573);  // Prime modulus near 2^{20}
   // NTL::ZZ m(1099511627791);  // Prime modulus near 2^{40}
   long numRep = 1; // Number of replications (multipliers) for each case.
   // bool inDual = true;  // Tests in dual lattice ?

   // These functions apply the tests with the desired types.
   // testLoop<int64_t, double>(conv<int64_t>(m), numRep, inDual);
   testLoop<NTL::ZZ, double>(m, t, numRep);
   //testLoop<NTL::ZZ, xdouble>(m, t, numRep);
   //testLoop<NTL::ZZ, quad_float>(m, numRep, inDual);
   //testLoop<NTL::ZZ, NTL::RR>(m, numRep, inDual);
}
