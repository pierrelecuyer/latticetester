/**
 * This example illustrates and compares different ways of making a search for good RNG
 * parameters in terms of the lattice structure. We do this for simple LCGs with
 * modulus m and multiplier a, based on the FOM `M_t` with a given vector `t = (t_1,...,t_d)`,
 * with coordinate 1 included in all the projections, for either the primal or the m-dual.
 * The FOM can be computed exactly using the BB, or only approximated via LLL.
 * When computing the FOM for a given a, we can either exit the procedure (early discard)
 * as soon as we know that the FOM will be too small, or always complete the computations.
 * The multipliers a that are considered can be read from a file, or they can be generated
 * as the successive powers of a given a0 modulo m.
 */

#include <iostream>
#include <vector>
#include <NTL/vector.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>

#include "latticetester/IntLattice.h"
#include "latticetester/Rank1Lattice.h"
#include "latticetester/FigureOfMeritM.h"
#include "latticetester/FigureOfMeritDualM.h"
#include "latticetester/ReducerStatic.h"
#include "latticetester/Weights.h"
#include "latticetester/NormaBestLat.h"
#include "latticetester/CoordinateSets.h"
#include "latticetester/WeightsUniform.h"
#include "latticetester/BasisConstruction.h"

using namespace LatticeTester;

/*
 * This static function loops over 'numMult' values of a for LCGs with modulus m,
 * and makes a list of the `numBest` ones based on the FOM 'fom'.
 * This FOM can be either for the primal or for the m-dual.
 * The variable `earlyDiscard` tells if we do early discard or not.
 * If `fromList` is true, the `numMult` values of a are taken from the array `inList`
 * and `a0` is not used, otherwise they are the successive powers of `a0`.
 * The `numBest` best candidates are returned in the list `outList`.
 * If `verbose` is true, they are also printed on the terminal, together with the FOM values.
 */
template<typename Int, typename Real>
static void findBestFOMs(const Int m, const Int a0, Rank1Lattice<Int, Real> &lat,
      IntLattice<Int, Real> &proj, FigureOfMeritM<Int, Real> *fom, bool earlyDiscard, bool fromList,
      Int inList[], int64_t numMult, Int outList[], int64_t numBest, bool verbose = false,
      std::string strMethod = " ") {
   Int a = a0;
   double merit;
   fom->setLowBound(0.0);
   double bestFoms[numBest];  // Array to store the best FoM values, in decreasing order.
   for (int64_t i = 0; i < numBest; i++) {
      outList[i] = 0;
      bestFoms[i] = 0.0;
   }
   int64_t posComp;  // Position in the list of the current candidate to dislodge.
   clock_t tmp;      // To measure computing time.
   tmp = clock();
   // We try numMult values of a.
   for (int64_t j = 0; j < numMult; j++) {
      if (fromList) a = inList[j];
      else a = a * a0 % m;
      lat.seta(a);
      merit = fom->computeMerit(lat, proj);  // Here we compute the FOM.
      posComp = numBest - 1;  // The last index in the list.
      if (merit > bestFoms[posComp]) {
         // Overwrite the currently smallest FoM and corresponding a in the stored list.
         for (; (posComp > 0) && (merit > bestFoms[posComp - 1]); posComp--) {
            outList[posComp] = outList[posComp - 1];
            bestFoms[posComp] = bestFoms[posComp - 1];
         }
         outList[posComp] = a;
         bestFoms[posComp] = merit;
         // If earlyDiscard, set the new lower bound to the current smallest FoM value in the list.
         // This value represents the FOM for the worst candidate still in the list.
         if (earlyDiscard) fom->setLowBound(bestFoms[posComp]);
      }
   }
   tmp = clock() - tmp;
   std::cout << "\n--------------------------------------\n";
   std::cout << "Function findBestFOMs with m = " << m << "\n";
   std::cout << std::boolalpha << "Discard: " << earlyDiscard << ",  from List: " << fromList;
   std::cout << ",  method: " << strMethod << "\n";
   std::cout << "Number of multipliers a examined: " << numMult << "\n";
   std::cout << "Number retained: " << numBest << "\n";
   std::cout << "Running Time in seconds: " << (double) tmp / (CLOCKS_PER_SEC) << "\n";
   if (verbose) {
      std::cout << "\nBest " << numBest << " multipliers `a` found, and their FOMs:\n";
      for (int64_t j = 0; j < numBest; j++)
         std::cout << "  " << outList[j] << ", " << bestFoms[j] << "\n";
   }
}

/*
 * This function tests different ways of making a search for the `numBest2` multipliers a
 * among `numMultLong` candidates, for the FOM M with vector `t`, for either the primal
 * (if `dual` is false) or the m-dual (if `dual` is true).
 * We start by testing a naive (inefficient) method that computes everything (no early discarding)
 * and applies the BKZ + BB for each multiplier a.  For that case, we examine only `numMultShort`
 * multipliers to avoid excessive computing times, then we can approximate the CPU time
 * required for `numMultLong` multipliers by a simple rule of Three.
 * In the second method, we also apply BKZ + BB, but we use early discarding.
 * In the third method, we apply only LLL and use early discarding.
 * In the fourth method, we use two stages, with early discarding at each stage,
 * The first stage applies only LLL and retains the `numBest1` multipliers in a list.
 * The second stage tests the retained multipliers using the BB and return the `numBest2` best ones.
 * The fifth method is similar to the fourth, except that in the first stage we use a smaller
 * vector t, with a restricted number of projections, to try speed up the screening.
 */
template<typename Int, typename Real>
static void compareSearchMethods(FigureOfMeritM<Int, Real> *fom, const Int m, const Int a0,
      const NTL::Vec<int64_t> t, const NTL::Vec<int64_t> t0, int64_t numMultLong,
      int64_t numMultShort, int64_t numBest0, int64_t numBest, NormType norm = L2NORM) {
   int64_t maxdim = t[0];  // Maximum dimension of the lattice
   // WeightsUniform weights(1.0);
   Rank1Lattice<Int, Real> lat(m, maxdim, norm);  // The current lattice for which the FoM is calculated.   ****  L2NORM CW
   IntLattice<Int, Real> proj(m, t.length(), norm); // Lattice used for projections. CW
   Int emptyList[0];
   Int inList[numBest0];
   Int outList[numBest0];

   // 1. Full computation (no discard), with default BKZ + BB.
 //  findBestFOMs(m, a0, lat, proj, fom, false, false, emptyList, numMultShort, outList, numBest,
 //        true, "BKZ + BB");

   // 2. Early discard, with BKZ + BB.
   findBestFOMs(m, a0, lat, proj, fom, true, false, emptyList, numMultLong, outList, numBest, true,
         "BKZ + BB");

   fom->setLLL(0.99999);
   fom->setBKZ(0.0);
   fom->setBB(false);
   // 3. Early discard, approximation with LLL only.
   findBestFOMs(m, a0, lat, proj, fom, true, false, emptyList, numMultLong, outList, numBest, true,
         "LLL only");

   // 4. Early discard, two stages, retain numBest0 in inList, with LLL only in first stage.
   findBestFOMs(m, a0, lat, proj, fom, true, false, emptyList, numMultLong, inList, numBest0, false,
         "LLL only, stage 1, with vector t");
   // We use LLL + BB on second stage, only for the numBest0 retained.
   fom->setBB(true);
   findBestFOMs(m, a0, lat, proj, fom, true, true, inList, numBest0, outList, numBest, true,
         "LLL + BB, stage 2");

   // 5. Early discard, two stages, retain numBest0 in inList, with LLL only and vector `t0` in first stage.
   fom->setTVector(t0, true);
   findBestFOMs(m, a0, lat, proj, fom, true, false, emptyList, numMultLong, inList, numBest0, false,
         "LLL only, stage 1 with vector t0");
   // We use LLL + BB on second stage, only for the numBest0 retained.
   fom->setTVector(t, true);   fom->setBB(true);
   findBestFOMs(m, a0, lat, proj, fom, true, true, inList, numBest0, outList, numBest, true,
         "LLL + BB, stage 2");
}

int main() {
   typedef NTL::ZZ Int;
   typedef double Real;
   std::cout << "Types: Int = NTL::ZZ, Real = double \n";
   
   NormType norm = L1NORM;

   NTL::ZZ m(1048573); // Prime modulus near 2^{20}
   NTL::ZZ a0(91);     // This a0 is a primitive element mod m=1048573.
   //NTL::ZZ m(1099511627791);  // Prime modulus near 2^{40}

   NTL::Vec<int64_t> t; // The t-vector for the FOM.
   t.SetLength(5);
   t[0] = 32;    // We look at successive coordinates in up to t[0] dimensions.
   t[1] = 32;    // Then pairs, triples, etc.
   t[2] = 16;
   t[3] = 12;
   t[4] = 10;

   NTL::Vec<int64_t> t0; // A reduced t-vector for the FOM.
   t0.SetLength(4);
   t0[0] = 4;
   t0[1] = 32;
   t0[2] = 16;
   t0[3] = 12;

   int64_t maxdim = t[0];  // Maximum dimension of the lattice
   WeightsUniform weights(1.0);
   NormaBestLat normaPrimal(log(m), 1, maxdim, norm);  // Factors computed for primal. 
   NormaBestLat normaDual(-log(m), 1, maxdim, norm);  // Factors computed for dual. 
   ReducerBB<Int, Real> red(maxdim);   // Single ReducerBB with internal lattice `lat`.
   FigureOfMeritM<Int, Real> fomPrimal(t, weights, normaPrimal, &red, true);
   FigureOfMeritDualM<Int, Real> fomDual(t, weights, normaDual, &red, true); // FoM for dual lattice.

   int64_t numMultLong = 100;  // Total number of multipliers to examine.
   int64_t numMultShort = 10;  // Number to examine in the `no early discard` case.
   int64_t numBest0 = 10;  // When doing two rounds, we retain `numBest1` from the first round.
   int64_t numBest = 3;   // We want the best three at the end.

   // We do first the primal, then the dual.
   std::cout << "\n=========================================================\n";
   std::cout << "FOM experiments in primal lattices \n";
   compareSearchMethods<Int, Real>(&fomPrimal, m, a0, t, t0, numMultLong, numMultShort, numBest0,
         numBest, norm);

   std::cout << "\n=========================================================\n";
   std::cout << "FOM experiments in dual lattices \n";
   compareSearchMethods<Int, Real>(&fomDual, m, a0, t, t0, numMultLong, numMultShort, numBest0,
         numBest, norm);
   std::cout << "\n***     DONE     ***\n";
   return 0;
}
