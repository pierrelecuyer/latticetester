/**
 * This example illustrates and compares different ways of making a search for good RNG
 * parameters in terms of the lattice structure. We do this for simple LCGs with
 * modulus m and multiplier a, based on the FOM M with a given vector t = (t_1,...,t_d),
 * with coordinate 1 included in all the projections, for either the primal or the m-dual.
 * The FOM can be computed exactly using the BB, or only approximated via LLL.
 * When computing the FOM for a given a, we can either exit the procedure (early discard)
 * as soon as we know that the FOM will be too small, or always complete the computations.
 * The multipliers a that are considered can be read from a file, or they can be generated
 * as the successive powers of a given a0 modulo m.
 */

// #define TYPES_CODE  ZD  // ZZ + double
// Is everything below really needed?
#include <iterator>
#include <iostream>
#include <cstdint>
#include <algorithm>
#include <vector>
#include <numeric>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>

#include "latticetester/FlexTypes.h"
#include "latticetester/EnumTypes.h"
#include "latticetester/Util.h"
#include "latticetester/IntLattice.h"
#include "latticetester/Rank1Lattice.h"
#include "latticetester/FigureOfMeritM.h"
#include "latticetester/FigureOfMeritDualM.h"
#include "latticetester/ReducerStatic.h"
// #include "latticetester/LLL_lt.h"
#include "latticetester/Weights.h"

// Application specific headers
#include "latticetester/NormaBestLat.h"
//#include "latticetester/NormaBestBound.h"
#include "latticetester/CoordinateSets.h"
#include "latticetester/WeightsUniform.h"
#include "latticetester/BasisConstruction.h"
//#include "latticetester/MRGLattice.h"

using namespace LatticeTester;

/*
 * This static function loops over 'numMult' values of a for LCGs with modulus m,
 * and makes a list of the `numBest` ones based on the FOM 'fom'.
 * This FOM can be either for the primal or for the m-dual.
 * The variable `earlyDiscard` tells if we do early discard or not.
 * If `fromList` is true, the `numMult` values of a are taken from the array `inList`,
 * otherwise they are the successive powers of `a0`.
 * The `numBest` best candidates are returned in the list `outList`.
 * If `verbose` is true, they are also printed on the terminal, together with
 * the FOM values and the required CPU time.
 */
template<typename Int, typename Real>
static void findBestFOMs(const Int m, const Int a0, Rank1Lattice<Int, Real> &lat,
      IntLattice<Int, Real> &proj, FigureOfMeritM<Int, Real> *fom, bool earlyDiscard, bool fromList,
      Int inList[], int64_t numMult, Int outList[], int64_t numBest, bool verbose = false) {
   Int a = a0;
   double merit;
   fom->setLowBound(0.0);
   // Int outList[numBest]; // Array to store the numBest multipliers with the best FoMs.
   double bestFoms[numBest];     // Array to store the best FoMs.
   for (int64_t i = 0; i < numBest; i++) {
      outList[i] = 0;
      bestFoms[i] = 0.0;
   }
   // Initialize variables for minimum elements of the numBest largest FoMs, get the corresponding index
   // and also get the maximal value.  ????    They are not kept sorted?
   auto adrMin = std::min_element(bestFoms, bestFoms + numBest); // Iterator to smallest FOM value in the list.
   auto posMin = std::distance(bestFoms, bestFoms); // Position of smallest FOM value in the list, here 0.

   clock_t tmp; // Variables for measuring time elapsed
   tmp = clock();
   // We try numMult values of a.
   for (int64_t j = 0; j < numMult; j++) {
      if (fromList) a = inList[j];
      else a = a * a0 % m;
      lat.seta(a);
      merit = fom->computeMerit(lat, proj);  // Here we compute the FOM.
      if (merit > bestFoms[posMin]) {
         // Overwrite the currently smallest FoM and best a in the stored list.
         outList[posMin] = a;
         bestFoms[posMin] = merit;
         // Find the new smallest FoM in the list.
         adrMin = std::min_element(bestFoms, bestFoms + numBest); // Address of smallest FOM value
         posMin = std::distance(bestFoms, adrMin); // Its position in the list.
         // If earlyDiscard, set the new lower bound to the current smallest FoM value in the list.
         // This value represents the FOM for the worst candidate still in the list.
         if (earlyDiscard) fom->setLowBound(bestFoms[posMin]);
      }
   }
   std::cout << "\n--------------------------------------\n";
   std::cout << "Function findBestFOMs with m = " << m << "\n";
   std::cout << "Number of multipliers a examined: " << numMult << "\n";
   tmp = clock() - tmp;
   std::cout << "Running Time in seconds: " << (double) tmp / (CLOCKS_PER_SEC) << "\n\n";
   if (verbose) {
      std::cout << "Best " << numBest << " multipliers `a` found, and their FOMs:\n";
      for (int64_t j = 0; j < numBest; j++)
         std::cout << "  " << outList[j] << ", " << bestFoms[j] << "\n";
   }
}

/*
 * This function tests different ways of making a search for the `numBest2` multipliers a
 * among `numMultLong` candidates, for the FOM M with vector `t`, first for the primal
 * lattice, then for the m-dual.
 * We start with a naive (inefficient) method that computes everything (no early discarding)
 * and applies the BB for each multiplier a.  For that case, we examine only `numMultShort`
 * multipliers to avoid excessive computing times, but we can easily approximate what the
 * CPU time would be for `numMultLong` multipliers.
 * In the second method, we also apply BB, but we use early discarding.
 * In the third method, we use two stages, with early discarding at each stage,
 * In the first stage, we use only LLL and retain the `numBest1` multipliers in a list.
 * In the second stage, we test these retained multipliers using the BB and return the `numBest2`
 * best ones.
 */
template<typename Int, typename Real>
static void compareSearchMethods(const Int m, const Int a0, const NTL::vector<int64_t> t,
      int64_t numMultLong, int64_t numMultShort, int64_t numBest1, int64_t numBest2) {
   int64_t maxdim = t[0];  // Maximum dimension of the lattice
   WeightsUniform weights(1.0);
   Rank1Lattice<Int, Real> lat(m, maxdim);  // The current lattice for which the FoM is calculated.
   IntLattice<Int, Real> proj(m, t.size()); // Lattice used for projections.
   ReducerBB<Int, Real> red(lat);  // Single ReducerBB with internal lattice `lat`.
   NormaBestLat normaPrimal(log(m), 1, maxdim);  // Factors computed for primal.
   NormaBestLat normaDual(-log(m), 1, maxdim);  // Factors computed for dual.
   FigureOfMeritM<Int, Real> fomPrimal(t, weights, normaPrimal, &red, true); // FoM for primal lattice.
   FigureOfMeritDualM<Int, Real> fomDual(t, weights, normaDual, &red, true); // FoM for dual lattice.

   Int emptyList[0];
   Int inList[numBest1];
   Int outList[numBest1];

   std::cout << "\n************************************************\n";
   std::cout << "FOM experiments in primal lattice \n";
   // Full computation (no discard), with BB, for primal.
   findBestFOMs(m, a0, lat, proj, &fomPrimal, false, false, emptyList, numMultShort, outList,
         numBest2, true);

   // Early discard, with BB, or primal.
   findBestFOMs(m, a0, lat, proj, &fomPrimal, true, false, emptyList, numMultLong, outList,
         numBest2, true);

   // Early discard, two stages, with LLL only to retain numBest1 on first stage, for primal.
   fomPrimal.setBB(false);
   findBestFOMs(m, a0, lat, proj, &fomPrimal, true, false, emptyList, numMultLong, inList, numBest1,
         false);
   fomPrimal.setBB(true);
   findBestFOMs(m, a0, lat, proj, &fomPrimal, true, true, inList, numBest1, outList, numBest2,
         true);

   // Same experiment for the m-dual.
   std::cout << "\n************************************************\n";
   std::cout << "FOM experiments in m-dual lattice \n";
   findBestFOMs(m, a0, lat, proj, &fomDual, false, false, emptyList, numMultShort, outList,
         numBest2, true);
   findBestFOMs(m, a0, lat, proj, &fomDual, true, false, emptyList, numMultLong, outList, numBest2,
         true);
   fomDual.setBB(false);
   findBestFOMs(m, a0, lat, proj, &fomDual, true, false, emptyList, numMultLong, inList, numBest1,
         false);
   fomDual.setBB(true);
   findBestFOMs(m, a0, lat, proj, &fomDual, true, true, inList, numBest1, outList, numBest2, true);
}

int main() {

   // Here, Int and Real are passed as template parameters.
   std::cout << "Types: NTL::ZZ, double \n";
   NTL::ZZ m(1048573); // Prime modulus near 2^{20}
   NTL::ZZ a0(91);     // This a0 is a primitive element mod m=1048573.
   // NTL::ZZ m(1099511627791);  // Prime modulus near 2^{40}

   int64_t numMultLong = 10000;  // Total number of multipliers to examine.
   int64_t numMultShort = 10000;  // Number to examine in the `no early discard` case.
   int64_t numBest1 = 10;  // When doing two rounds, we retain `numBest1` from the first round.
   int64_t numBest2 = 3;   // We want the best three at the end.

   NTL::vector<int64_t> t(4); // The t-vector for the FOM.
   t[0] = 16;    // We look at successive coordinates in up to t[0] dimensions.
   t[1] = 16;    // Then pairs, triples, etc.
   t[2] = 12;
   t[3] = 10;

   compareSearchMethods<NTL::ZZ, double> (m, a0, t, numMultLong, numMultShort, numBest1, numBest2);
   return 0;
}
