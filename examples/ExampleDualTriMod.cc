// File `ExampleDualTriMod.cc`

#define TYPES_CODE  LD     // Int = int64_t, Real = double
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include "latticetester/FlexTypes.h"
#include "latticetester/Util.h"
#include "latticetester/IntLattice.h"
#include "latticetester/BasisConstruction.h"

/**
 * This gives a small example in which taking the m-dual of a triangular basis
 * can give non-diagonal entries larger than m.
 * The initial basis is for a small MRG with k=3, m=13, a_1=7, a_2=5, a_3=11.
 * We set up a basis, compute a triangular basis, then its m-dual, 
 * and then reduce it mod m towards 0.
 */
using namespace NTL;
using namespace LatticeTester;

const int64_t dim = 4;
const int64_t k(3);
const int64_t m(13);
const int64_t aa[k] = {7, 5, 11};
int64_t yy[dim+k-1] = {0, 0, 1, 0, 0, 0};

int main() {
   int64_t i, j;
   std::cout << "*************************************************\n";
   std::cout << "ExampleDualTriMod with:  k = " << k << ", m = " << m << "\n";
   for (i = 0; i < k; i++) 
      std::cout << " a_" << i+1 << " = " << aa[i] << ",  ";
   std::cout << "\n";
   IntMat basis0, basis1, basisdual;
   basis0.SetDims(dim, dim);
   basis1.SetDims(dim, dim);
   basisdual.SetDims(dim, dim);
   for (i = k; i < dim+k-1; i++) {
      yy[i] = (aa[0] * yy[i-1] + aa[1] * yy[i-2] + aa[2] * yy[i-3]) % m;
      std::cout << " y_" << i << " = " << yy[i] << ",  ";
      }
   std::cout << "\n\n";
   for (i = 0; i < dim; i++) {
      for (j = 0; j < dim; j++) {
         if (i < k)
            basis0[i][j] = yy[k-1-i+j];
         else
            if (i == j) basis0[i][j] = m;
            else basis0[i][j] = 0;
      }
   }
   std::cout << "Initial basis0:\n" << basis0 << "\n";
   upperTriangularBasis<Int>(basis1, basis0, m, dim, dim);
   std::cout << "upper-triangular primal basis V:\n" << basis1 << "\n";

   mDualUpperTriangular(basisdual, basis1, m, dim);
   std::cout << "lower-triangular basisdual W:\n" << basisdual << "\n";
   // Here we verify that basisdual is indeed the m-dual of basis1. 
   bool inverse = checkInverseModm<Int> (basisdual, basis1, m);
   std::cout << "Is basisdual the m-dual of basis1? " << inverse << "\n";
   
   mDualUpperTriangularMod0(basisdual, basis1, m, dim);
   std::cout << "lower-triangular basisdual W' after mod_0:\n" << basisdual << "\n";
   // Here basisdual is NOT the m-dual of basis1, but it is another basis of the m-dual lattice. 
   inverse = checkInverseModm<Int> (basisdual, basis1, m);
   std::cout << "Is basisdual the m-dual of basis1? " << inverse << "\n";
   return 0;
}

