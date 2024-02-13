/**
 * This example illustrates the use of functions from the `BasisConstruction` class,
 * with a small five-dimensional lattice obtained from an LCG with a small modulus.
 * The basis is printed between function calls, to show what is going on.
 * An initial upper-triangular basis is constructed by `Rank1Lattice` in a standard way.
 * We then call LLL with `delta = 0.5` to obtain a basis with shorter vectors,
 * which is not triangular. Calling `upperTriangulatBasis` transforms this basis
 * to an upper triangular basis, which happens to be the same as the initial one.
 * Different triangularization methods are compared.  Then LLL with different
 * values of `delta`. Then a dual basis is computed in different ways.
 *
 * After that, we look at projections of this lattice over subsets of coordinates.
 * We show how to construct a basis for such a projection, compute the corresponding
 * dual basis, and compute a shortest vector in this dual basis.
 *
 * This program can be run with the two different types of `Int` (just change the first
 * line below), and with various choices of the modulus `m` and multiplier `a`.
 **/

// #define TYPES_CODE  LD     // int64_t
#define TYPES_CODE  ZD     // ZZ

#include <iostream>
#include <cstdint>
// #include <ctime>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include "latticetester/FlexTypes.h"    // This defines Int = int64_t
#include "latticetester/EnumTypes.h"
#include "latticetester/BasisConstruction.h"
#include "latticetester/Util.h"
#include "latticetester/IntLattice.h"
#include "latticetester/Rank1Lattice.h"
#include "latticetester/Reducer.h"
#include "latticetester/Coordinates.h"

using namespace LatticeTester;

Int m(101);      // Modulus m = 101
//Int m(1021);     // Modulus m = 1021
//Int m(1048573);  // Modulus m = 1048573 (prime number near 2^{20})
Int a(33);       // The LCG multiplier
const long dim = 5;  // Dimension
const long dimProj = 3;  // Dimension

int main() {
    std::cout << "Types: " << strFlexTypes << "\n";

    // All the IntMat objects are created in 5 dimensions, but we may use less.
    IntMat basis1, basis2, basisProj, basisDual;
    basis1.SetDims(dim, dim);
    basis2.SetDims(dim, dim);
    basisProj.SetDims(dim, dim);
    basisDual.SetDims(dim, dim);
    Int sqlength;

    // We construct a Korobov lattice in dim dimensions.
    Rank1Lattice<Int, double> *korlat;
    korlat = new Rank1Lattice<Int, Real>(m, a, dim);
    korlat->buildBasis(dim);
    copy(korlat->getBasis(), basis1);  // This initial basis is triangular.
    std::cout << "Initial Korobov lattice basis = \n" << basis1 << "\n";
    ProdScal<Int>(basis1[0], basis1[0], dim, sqlength);
    std::cout << "Square length of first basis vector: " << sqlength << "\n\n";

    // We apply LLL to reduce basis1.
    BasisConstruction<Int>::LLLConstruction0(basis1, 0.5);
    std::cout << "Basis after LLL with delta=0.5: \n" << basis1 << "\n";
    ProdScal<Int>(basis1[0], basis1[0], dim, sqlength);
    std::cout << "Square length of first basis vector: " << sqlength << "\n\n";

    BasisConstruction<Int>::LLLConstruction0(basis1, 0.8);
    std::cout << "Basis after LLL with delta=0.8: \n" << basis1 << "\n";
    ProdScal<Int>(basis1[0], basis1[0], dim, sqlength);
    std::cout << "Square length of first basis vector: " << sqlength << "\n\n";

    BasisConstruction<Int>::LLLConstruction0(basis1, 0.99999);
    std::cout << "Basis after LLL with delta=0.99999: \n" << basis1 << "\n";
    ProdScal<Int>(basis1[0], basis1[0], dim, sqlength);
    std::cout << "Square length of first basis vector: " << sqlength << "\n\n";

    // We now transform basis1 to the upper-triangular basis2.
    // Note that after this, basis1 contains only garbage.
    BasisConstruction<Int>::upperTriangularBasis(basis1, basis2, m);
    std::cout << "After `upperTriangularBasis`: \n" << basis2 << "\n\n";

    // Then we compute the m-dual of basis2 and put it in basisDual.
    BasisConstruction<Int>::mDualUpperTriangular(basis2, basisDual, m);
    std::cout << "m-dual upperTriangular: \n" << basisDual << "\n\n";

    // Here we compute the guaranteed shortest vector, with BB.
    Reducer<Int, Real> *red = new Reducer<Int, Real>(*korlat);
    red->shortestVector(); // For this, we need to create a Reducer object.
    std::cout << "Shortest vector: " << korlat->getBasis()[0] << "\n";
    std::cout << "Its square length: " << korlat->getVecNorm[0] << "\n\n";

    // We now investigate the projection over coordinates {1, 3, 5}.
    // We first insert those three coordinates one by one in `proj`.
    // We then compute a basis for this projection in two ways.
    Coordinates proj;
    proj.insert(1);
    proj.insert(3);
    proj.insert(5);
    std::cout << "Lattice projection over coordinates " << proj << ".\n";
    BasisConstruction<Int>::projectionConstructionLLL(basis2, basisProj, proj,         m);
    std::cout << "Basis for this projection, with LLL: \n" << basisProj << "\n";
    BasisConstruction<Int>::projectionConstructionUpperTri(basis2, basisProj,
            proj);
    std::cout << "Upper-triangular basis for this proj.: \n" << basisProj
            << "\n";

    // We use only three coordinates of these matrices for the projection.
    BasisConstruction<Int>::mDualUpperTriangular(basisProj, basisDual, m, 3);
    std::cout << "Triangular basis for m-dual of this projection: \n"
            << basisDual << "\n";
    BasisConstruction<Int>::LLLConstruction0(basisdual, 0.99999, 3, 3);
    std::cout << "m-dual basis after LLL with delta=0.99999: \n" << basisDual
            << "\n";
    return 0;
}

