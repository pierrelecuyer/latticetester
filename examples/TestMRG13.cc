/**
 * This example shows how to use `LatticeTester` to calculate a figure of merit in practice
 * In a first loop an approximation of the FoM is calculated for many multipliers (numRep)
 * by using the chosen pre-reduction algorithm (meth). A chosen number (noBest) of the 
 * best multipliers according to this loop is stored. Afterwards the exact FoM for these 
 * stored multipliers is calculated by means of the BB algorithm.
 */

#define TYPES_CODE  ZD

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
#include "latticetester/LLL_lt.h"
#include "latticetester/Weights.h"
#include "latticetester/WeightsOrderDependent.h"

// Application specific headers
#include "latticetester/NormaBestLat.h"
#include "latticetester/NormaBestBound.h"
#include "latticetester/CoordinateSets.h"
#include "latticetester/WeightsUniform.h"
#include "latticetester/BasisConstruction.h"
#include "latticetester/MRGLattice.h"

using namespace LatticeTester;

int main() {
   /*
    * This is to illustrate the problems with the primal projections described at the end of section 9 of the guide.
    */
   std::cout << "Results of example TestMRG13 \n";
   std::cout << "Types: " << strFlexTypes << "\n\n";

   int64_t order(3);
   Int m(13);
   NTL::vector<Int> a(order); // Vector a has size 3.
   a[0] = 7;
   a[1] = 0;
   a[2] = 4;
   int64_t maxdim = 8;  // Maximum dimension of the lattice
   NTL::vector<int64_t> t(3); // The t-vector
   t[0] = 8;     // We look at successive coordinates in up to t[0] dimensions.
   t[1] = 5;
   t[2] = 4;  // Then pairs up to coord. 5 and triples up to coord. 4.
   double merit;

   // Defining the required lattices and tools.
   MRGLattice<Int, Real> lat(m, a, maxdim);
   IntLattice<Int, Real> proj(m, 3);
   WeightsUniform weights(1.0);
   NormaBestLat norma(log(m), order, maxdim);  // Factors will be computed for primal.
   ReducerBB<Int, Real> red(maxdim);   // Reducer created for up to dim dimensions.
   FigureOfMeritM<Int, Real> fom(t, weights, norma, &red);
   fom.setVerbosity(2);

   // Building full primal lattice and looking at projections over successive coordinates.
   std::cout << "===================================================\n";
   std::cout << "Looking at projections over successive coordinates.\n";
   std::cout << "\nFigure of merit primal succ, with BB.\n";
   merit = fom.computeMeritSucc(lat);
   std::cout << "FOM value: " << merit << "\n\n";

   // Now looking at other primal projections, placed in lattice proj.
   std::cout << "===================================================\n";
   std::cout << "Then we look at other primal projections, over pairs and triples.\n";
   std::cout << "\nFigure of merit primal non-succ, with BB.\n";
   merit = fom.computeMeritNonSucc(lat, proj);
   std::cout << "FOM value: " << merit << "\n\n";

   // We now examine the m-dual lattice, first for successive coordinates.
   std::cout << "===================================================\n";
   std::cout << "We now examine the m-dual lattice, first for successive coordinates.\n";
   NormaBestLat normadual(-log(m), order, maxdim);  // Factors will be computed for m-dual.
   FigureOfMeritDualM<Int, Real> fomdual(t, weights, normadual, &red);
   fomdual.setVerbosity(2);
   std::cout << "\nFigure of merit for dual.\n";
   merit = fomdual.computeMeritSucc(lat);
   std::cout << "FOM value: " << merit << "\n\n";

   // Then the m-duals of the projections.
   std::cout << "===================================================\n";
   std::cout << "Then we look at the m-duals of the other projections.\n";
   std::cout << "\nFigure of merit dual non-succ, with BB.\n";
   merit = fomdual.computeMeritNonSucc(lat, proj);
   std::cout << "FOM value: " << merit << "\n\n";

   // A closer look at projection {1,3,4}
   std::cout << "===================================================\n";
   std::cout << "Let us have a closer look at the projection over coordinates {1,3,4}. \n";
   Coordinates coord;
   coord.insert(1);
   coord.insert(3);
   coord.insert(4);
   lat.buildBasis(maxdim);
   std::cout << "Full basis B before taking projection {1,3,4}: \n" << lat.getBasis() << "\n";
   lat.buildProjection(proj, coord);
   std::cout << "Basis for projection {1,3,4}: \n" << proj.getBasis() << "\n";
   fom.setVerbosity(0);
   fom.computeMeritOneProj(proj, coord);

   lat.buildProjectionDual(proj, coord);
   proj.dualize();
   std::cout << "Basis for m-dual projection {1,3,4}: \n" << proj.getBasis() << "\n";
   fomdual.setVerbosity(0);
   fomdual.computeMeritOneProj(proj, coord);
   proj.dualize();

   /*
    for (int a1 = 1; a1 < 13; a1++) {
    a[0] = a1;
    for (int a3 = 1; a3 < 13; a3++) {
    a[2] = a3;
    lat.setaa(a);
    lat.buildBasis(4);
    lat.buildProjection(proj, coord);
    merit = fom.computeMeritOneProj (proj, coord);
    if (merit > 1.99) {
    std::cout << " merit = " << merit << ",  a1 = " << a1 << ",  a3 = " << a3 << "\n";
    lat.buildProjection(proj, coord);
    std::cout << " proj basis = \n" << proj.getBasis() << "\n";
    }
    }
    }
    */

   return 0;
}
