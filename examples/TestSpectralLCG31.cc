/*
 * TestSpectralLCG.cc
 */

// This defines the Int type. We must recompile to change it.
#define TYPES_CODE  ZD     // Int = ZZ, Real = double

#include <NTL/ZZ.h>
#include "latticetester/FlexTypes.h"
#include "latticetester/Chrono.h"
#include "latticetester/BasisConstruction.h"
#include "latticetester/Rank1Lattice.h"
#include "latticetester/ReducerStatic.h"
#include "latticetester/NormaBestLat.h"
#include "latticetester/NormaRogers.h"
#include "latticetester/WeightsUniform.h"
#include "latticetester/FigureOfMeritDualM.h"

using namespace LatticeTester;

int main() {
   int64_t maxdim(30);  // Maximum dimension of the lattice.
   WeightsUniform weights(1.0);
   Int m(2147483647);   // Modulus.
   NormaBestLat norma(log(m), 1, maxdim);      // Normalization for primal.
   NormaBestLat normaDual(-log(m), 1, maxdim); // Normalization for m-dual.
   ReducerBB<Int, Real> red(maxdim);           // Single ReducerBB.
   NTL::Vec <int64_t> t; // The t-vector for the FOM, used only for pairs.
   t.SetLength(2);
   t[0] = 0;
   t[1] = maxdim;  // We look at pairs {1,s} for s up to maxdim.
   IntLattice<Int, Real> proj(m, t.length());
   FigureOfMeritM<Int, Real> fom(t, weights, norma, &red, true); // FoM for dual lattice.
   FigureOfMeritDualM<Int, Real> fomdual(t, weights, normaDual, &red, true); // FoM for dual lattice.
   fom.setTVector(t, true);
   fomdual.setTVector(t, true);
   fom.setVerbosity(3);
   fomdual.setVerbosity(3);
   std::cout << "\nResults from TestSpectralLCG.cc \n";

   Int a(16807);
   std::cout << "\n=============================================\n";
   std::cout << "LCG-16897, with BKZ+BB with default parameters, succ. coordinates \n";
   Rank1Lattice<Int, Real> lcg = Rank1Lattice<Int, Real>(m, a, maxdim, maxdim);
   std::cout << "\nPrimal lattice:\n";
   fom.computeMeritSucc(lcg, 2, maxdim);
   std::cout << "\nDual lattice:\n";
   fomdual.computeMeritSucc(lcg, 2, maxdim);
   std::cout << "\n=========================\n";
   std::cout << "LCG-16807, projections over pairs of coordinates \n";
   std::cout << "\nPrimal lattice:\n";
   fom.computeMerit(lcg, proj);
   std::cout << "\nDual lattice:\n";
   fomdual.computeMerit(lcg, proj);

   std::cout << "\n=============================================\n";
   std::cout << "LCG-742938285, with BKZ+BB with default parameters, succ. coordinates \n";
   a = 742938285;
   lcg.seta(a);
   std::cout << "\nPrimal lattice:\n";
   fom.computeMeritSucc(lcg, 2, maxdim);
   std::cout << "\nDual lattice:\n";
   fomdual.computeMeritSucc(lcg, 2, maxdim);
   std::cout << "\n===========================\n";
   std::cout << "LCG-742938285, projections over pairs of coordinates \n";
   std::cout << "\nPrimal lattice:\n";
   fom.computeMerit(lcg, proj);
   std::cout << "\nDual lattice:\n";
   fomdual.computeMerit(lcg, proj);
   std::cout << "\nDone";
   return 0;
}


