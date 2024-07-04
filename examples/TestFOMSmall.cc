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

// Here we compute the FOM for an LCG with modulus m and multiplier a, for projections specified by t.
// We first generate `numRepVerb` multipliers `a` and for each one we compute and display
// (on the terminal) the FOM and the worst projection in both the primal and m-dual.
// Then we generate `numRepStat` different multipliers and we count how many times
// the worst-case proj. has s dimensions, for each s, for both the primal and m-dual.
// We also print in a file the pairs (primal FOM, dual FOM) for each a, to make a scatter plot.
template<typename Int, typename Real>
void testLoop (const Int m, const NTL::vector<int64_t> t, long numRepVerb, long numRepStat) {
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
   // fomprimal.setVerbosity(2);
   // fomdual.setVerbosity(2);
   IntLattice<Int, Real> proj(m, t.size());
   Int a;        // The LCG multiplier
   double merit;

   // This part is to visualize the FOM values for a few multipliers a.
   for (int64_t r = 0; r < numRepVerb; r++) {
      a = (m / 5 + 1021 * r) % m;   // The multiplier we use for this rep.
      korlat.seta(a);
      std::cout << "a = " << a << "\n";
      merit = fomprimal.computeMerit(korlat, proj);
      std::cout << "FOM in primal: " << merit << ",  worst projection: " << fomprimal.getMinMeritProj() << "\n";
      merit = fomdual.computeMerit(korlat, proj);
      std::cout << "FOM in dual: " << merit << ",  worst projection: " << fomdual.getMinMeritProj() << "\n\n";
   }

   // This part is to collect statistics for a larger sample of values of a.
   // int64_t countDimOfMinPrimal[maxdim+1];  // Vectors to store the counts.
   NTL::vector<int64_t> countDimOfMinPrimal(maxdim);
   NTL::vector<int64_t> countDimOfMinDual(maxdim);
   for (int64_t j = 0; j <= maxdim; j++) {
      countDimOfMinPrimal[j] = 0;
      countDimOfMinDual[j] = 0;
   }
   ofstream pairsFile;
   pairsFile.open ("pairsPrimalDualFOM.txt");

   for (int64_t r = numRepVerb; r < numRepVerb + numRepStat; r++) {
      a = (m / 5 + 1021 * r) % m;   // The multiplier we use for this rep.
      korlat.seta(a);
      merit = fomprimal.computeMerit(korlat, proj);
      countDimOfMinPrimal[(fomprimal.getMinMeritProj()).size()-2]++;
      pairsFile << merit << "  ";
      merit = fomdual.computeMerit(korlat, proj);
      countDimOfMinDual[(fomdual.getMinMeritProj()).size()-2]++;
      pairsFile << merit << "\n";
   }
   std::cout << "Counts for dimension of worst-case projection (starts at 2):" << "\n";
   std::cout << "  primal lattice: " << countDimOfMinPrimal << "\n";
   std::cout << "  dual lattice:   " << countDimOfMinDual << "\n";
   pairsFile.close();
}

int main() {

   NTL::vector<int64_t> t(4); // The t-vector for the FOM.
   t[0] = 16;    // We look at successive coordinates in up to t[0] dimensions.
   t[1] = 16;    // Then pairs and triples up to coord. 5.
   t[2] = 12;
   t[3] = 10;

   // Here, Int and Real are passed as template parameters.
   std::cout << "Types: NTL::ZZ, double \n";
   // NTL::ZZ m(101);  // Prime modulus near 100
   NTL::ZZ m(1048573);  // Prime modulus near 2^{20}
   // NTL::ZZ m(1099511627791);  // Prime modulus near 2^{40}

   // These functions apply the tests with the desired types.
   // testLoop<int64_t, double>(conv<int64_t>(m), numRep, inDual);
   testLoop<NTL::ZZ, double>(m, t, 20, 1000);
   //testLoop<NTL::ZZ, quad_float>(m, numRep, inDual);
}
