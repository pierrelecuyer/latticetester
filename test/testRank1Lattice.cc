#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include <iostream>
#include <cstdint>
#include <cmath>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include "latticetester/FlexTypes.h"
#include "latticetester/EnumTypes.h"
#include "latticetester/Rank1Lattice.h"
#include "latticetester/BasisConstruction.h"

namespace LatticeTester {

/**
 * Takes a modulus 'm' and a multiplier 'a' and builds and checks if the
 * basis and the dual basis from dimensions 'dim' to 'maxdim' are built
 * correctly. Also checks if the projection to all existing odd dimensions
 * works fine.
 */
template<typename Int, typename Real>
void BaseTest(Int m, Int a, int64_t maxdim, int dim) {
     Rank1Lattice<Int, Real> lat(m, maxdim);
     IntLattice<Int, Real> proj(m, maxdim);
     // First we check the creation of the basis and dual basis in all given dimensions
     lat.seta(a);
     lat.buildBasis(dim);
     lat.buildDualBasis(dim);
     NTL::Mat<NTL::ZZ> basis, dualbasis;
     basis.SetDims(maxdim, maxdim);
     dualbasis.SetDims(maxdim, maxdim);
     basis[0][0] = 1;
     for (int i = 1; i < dim; i++)
          basis[i][i] = m;
     basis[0][1] = a;
     for (int i = 2; i < dim; i++) 
          basis[0][i] = basis[0][i-1] * a % m;
     CHECK(lat.getBasis() == basis); 
     mDualUpperTriangular(dualbasis, lat.getBasis(), m, dim);
     CHECK(lat.getDualBasis() == dualbasis); 
     for (int i = dim; i < maxdim; i++) {
          lat.incDimBasis();
          basis[0][i] = basis[0][i-1] * a % m;
          basis[i][i] = m;
          CHECK(lat.getBasis() == basis); 
          mDualUpperTriangular(dualbasis, lat.getBasis(), m, i+1); 
          lat.incDimDualBasis();         
          CHECK(lat.getDualBasis() == dualbasis); 
     }
     // Now we check the projection 
     Coordinates coord;
     for (int i = 0; i < maxdim/2; i++)
          coord.insert(2*i+1);
     // Reset the matrix for basis and the dualbasis 
     for (int i = 0; i < maxdim; i++) {
          for (int j = 0; j < maxdim; j++) {
               basis[i][j] = 0;
               dualbasis[i][j] = 0;
          }
     }
     basis[0][0] = 1;
     basis[0][1] = a * a % m;     
     for (int i = 0; i < maxdim/2-1; i++) {
          basis[i+1][i+1] = m;
          basis[0][i+1] = basis[0][i] * a * a % m;
     }
     lat.buildProjection(proj, coord);
     CHECK(proj.getBasis() == basis);
     lat.buildProjectionDual(proj, coord);
     mDualUpperTriangular(dualbasis, lat.getBasis(), m, coord.size()); 
     CHECK(proj.getDualBasis() == dualbasis);
}


TEST_CASE("Test basis and dual basis") {
     int64_t maxdim = 8;
     int dim = 4;
     NTL::ZZ m; m = 1021;
     NTL::ZZ a; a = 822;
     BaseTest<NTL::ZZ, double> (m, a, maxdim, dim);
     BaseTest<NTL::ZZ, NTL::RR> (m, a, maxdim, dim);

}

};