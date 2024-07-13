/**
 * In this example, we compute FOMs for LCGs with modulus m and multiplier a,
 * for projections specified by t, and `numRepStat` values of a.
 * For each a, we compute the FOM for the primal and for the m-dual, we find
 * the worst-case projection and its dimension in each case, and we print these values in
 * a data file that can be used to make a scatter plot (e.g., with PGFplots).
 * The file has one row for each a, and the row contains the FOM and the dimension of the
 * worst projection for both the primal and the m-dual, plus an integer that indicates
 * if the worst projection is the same for both the primal and dual.
 * The program also reports how many times the worst projection has s dimensions for each
 * integer s > 1, for both the primal and the m-dual.
 * And it displays these projections explicitly for the first few values of a.
 *
 * See the Lattice Tester Guide for more explanations and results.
 */

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include "latticetester/FlexTypes.h"
#include "latticetester/Util.h"
#include "latticetester/Chrono.h"
#include "latticetester/Rank1Lattice.h"
#include "latticetester/ReducerStatic.h"
#include "latticetester/ReducerBB.h"
#include "latticetester/FigureOfMeritM.h"
#include "latticetester/FigureOfMeritDualM.h"
#include "latticetester/CoordinateSets.h"
#include "latticetester/WeightsUniform.h"

using namespace LatticeTester;

/**
 * This is the main function that runs the experiment. The LCGs will have modulus m and
 * the multipliers a will be the successive powers of `a0` modulo `m`.
 * The FOM will be M_t with the given vector t, only for the projections that include coordinate 1.
 * We consider `numRepStat` values of a and information is displayed on the terminal for the
 * first `numRepVerb` of them. The results for the scatter plot will be in file `outFile`.
 */
template<typename Int, typename Real>
void testLoop(const Int &m, Int &a0, const NTL::vector<int64_t> t, long numRepVerb, long numRepStat,
        std::string &outFile) {
   std::string stringTypes;  // To print the selected flexible types.
   strTypes<Int, Real>(stringTypes);  // Function from FlexTypes
   std::cout << "\n************************************************************************\n";
   std::cout << "Types: " << stringTypes << "\n";
   std::cout << "We look at FOM properties in primal and m-dual for a small list of LCGs.\n";
   std::cout << "TestFOMPairs with m = " << m << ", vector t = " << t << "\n\n";
   std::cout << "-----------------------------------------------------\n";

   int64_t maxdim = t[0];  // Maximum dimension of the lattice
   Rank1Lattice<Int, Real> korlat(m, maxdim); // We use single lattice object.
   ReducerBB<Int, Real> red(korlat);  // Single ReducerBB with internal lattice `korlat`.
   WeightsUniform weights(1.0);
   NormaBestLat normaPrimal(log(m), 1, maxdim);  // Factors computed for primal.
   NormaBestLat normaDual(-log(m), 1, maxdim);  // Factors computed for dual.
   // normaDual.computeBounds(-log(m), 1);
   FigureOfMeritM<Int, Real> fomPrimal(t, weights, normaPrimal, &red, true);
   FigureOfMeritDualM<Int, Real> fomDual(t, weights, normaDual, &red, true);
   IntLattice<Int, Real> proj(m, t.size());
   Int a = a0;        // The LCG multiplier
   double meritPrimal, meritDual;         // The merit values, for primal and dual.
   int64_t dimMeritPrimal, dimMeritDual;  // The dimensions of the worst-case projections.
   int64_t sameProj;  // Dim of the worst projection when it is the same for primal and dual.
   Chrono timer;

   // This part is to collect statistics for a large sample of values of a.
   NTL::vector<int64_t> countDimOfMinPrimal(maxdim);
   NTL::vector<int64_t> countDimOfMinDual(maxdim);
   for (int64_t j = 0; j <= maxdim; j++) {
      countDimOfMinPrimal[j] = 0;
      countDimOfMinDual[j] = 0;
   }
   ofstream pairsFile;  // Output file.
   pairsFile.open(outFile);
   pairsFile << " primalDim  primalFOM  dualDim  dualFOM  sameProj \n";
   for (int64_t r = 0; r < numRepStat; r++) {
      a = a0 * a % m;   // The multiplier we use for this rep.
      korlat.seta(a);
      meritPrimal = fomPrimal.computeMerit(korlat, proj);
      dimMeritPrimal = (fomPrimal.getMinMeritProj()).size();
      countDimOfMinPrimal[dimMeritPrimal - 2]++;
      meritDual = fomDual.computeMerit(korlat, proj);
      dimMeritDual = (fomDual.getMinMeritProj()).size();
      countDimOfMinDual[dimMeritDual - 2]++;
      // If the worst projection is the same, we will print its dimension.
      if (fomPrimal.getMinMeritProj() == fomDual.getMinMeritProj()) sameProj = dimMeritPrimal;
      else sameProj = 0;
      pairsFile << "      " << dimMeritPrimal << "    " << meritPrimal << "      ";
      pairsFile << dimMeritDual << "    " << meritDual << "    " << sameProj << "\n";
      if (r < numRepVerb) {  // We display the worst projections.
         std::cout << "a = " << a << "\n";
         std::cout << "FOM in primal: " << meritPrimal << ",  worst projection: "
               << fomPrimal.getMinMeritProj() << "\n";
         std::cout << "FOM in dual: " << meritDual << ",  worst projection: "
               << fomDual.getMinMeritProj() << "\n\n";
      }
   }
   pairsFile.close();
   std::cout << "\nCounts for dimension of worst-case projection (starts at 2 dim.):" << "\n";
   std::cout << "  primal lattice: " << countDimOfMinPrimal << "\n";
   std::cout << "  dual lattice:   " << countDimOfMinDual << "\n\n";
   std::cout << "CPU time (seconds): " << timer.val(Chrono::SEC) << "\n";
}

int main() {

   // Here, Int and Real are passed as template parameters.
   std::cout << "Types: NTL::ZZ, double \n";
   NTL::ZZ m(1048573); // Prime modulus near 2^{20}
   NTL::ZZ a0(91);     // This a0 is a primitive element mod m=1048573.
   // NTL::ZZ m(1099511627791);  // Prime modulus near 2^{40}

   NTL::vector<int64_t> t(4); // The t-vector for the FOM.
   t[0] = 16;    // We look at successive coordinates in up to t[0] dimensions.
   t[1] = 16;    // Also pairs, triples, and quadruples.
   t[2] = 12;
   t[3] = 10;
   std::string outFile = "pairsPrimalDualFOM16.dat";
   testLoop<NTL::ZZ, double>(m, a0, t, 10, 1000, outFile);

   t.resize(5);
   t[0] = 32;    // We look at successive coordinates in up to t[0] dimensions.
   t[1] = 32;    // Also pairs, triples, etc.
   t[2] = 16;
   t[3] = 12;
   t[4] = 10;
   outFile = "pairsPrimalDualFOM32.dat";
   testLoop<NTL::ZZ, double>(m, a0, t, 10, 1000, outFile);
}
