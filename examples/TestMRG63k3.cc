/**
 * This example shows how to use `LatticeTester` to calculate a figure of merit in practice
 * In a first loop an approximation of the FoM is calculated for many multipliers (numRep)
 * by using the chosen pre-reduction algorithm (meth). A chosen number (noBest) of the 
 * best multipliers according to this loop is stored. Afterwards the exact FoM for these 
 * stored multipliers is calculated by means of the BB algorithm.
 */

#define TYPES_CODE  ZD  // ZZ + quad_float

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
    * This tests the MRG retained in Table VII of the LatMRG paper of L'Ecuyer and Couture (1997).
    */
   std::cout << "Results of example TestMRG63k3 \n";
   std::cout << "Types: " << strFlexTypes << "\n\n";

   int64_t order(3);
   Int m(9223372036854773561);
   NTL::vector<Int> a(order); // Vector a has size 3.
   a[0] = 1145902849652723;
   a[1] = 0;
   a[2] = -1184153554609676;
   int64_t maxdim = 8;  // Maximum dimension of the lattice
   NTL::vector<int64_t> t(3); // The t-vector
   // t.SetLength(4);
   t[0] = 8;
   t[1] = 4;
   t[2] = 4;
   double merit;

   // Building full primal lattice and looking at projections over successive coordinates.
   std::cout << "===================================================\n";
   std::cout << "We build the lattice and look at projections over successive coordinates.\n";
   MRGLattice<Int, Real> lat(m, a, maxdim);
   WeightsUniform weights(1.0);
   NormaBestLat norma(log(m), order, maxdim);  // Factors will be computed for primal.
   ReducerBB<Int, Real> red(maxdim);   // Reducer created for up to dim dimensions.
   FigureOfMeritM<Int, Real> fom(t, weights, norma, &red);
   fom.setVerbosity(2);
   std::cout << "\nFigure of merit primal succ, with BB.\n";
   merit = fom.computeMeritSucc(lat);
   std::cout << "FOM value: " << merit << "\n\n";

   // Now looking at other primal projections, placed in lattice proj.
   std::cout << "===================================================\n";
   std::cout << "Then we look at other primal projections, over pairs and triples.\n";
   IntLattice<Int, Real> proj(m, 3);
   lat.buildBasis(maxdim);
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

   return 0;
}
