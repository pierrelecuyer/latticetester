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
 * We show how to construct a basis for such a projection and compute the corresponding
 * dual basis.
 *
 * This program can be run with the three different types of `Int` (just change the first
 * line below), and with various choices of the modulus `m` and multiplier `a`.
 **/

#define TYPES_CODE  LD     // int64_t + double
//#define TYPES_CODE  ZD     // ZZ + double
//#define TYPES_CODE  ZR     // ZZ + RR

#include <iostream>
#include <cstdint>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include "latticetester/FlexTypes.h"    // This defines Int = int64_t
#include "latticetester/EnumTypes.h"
#include "latticetester/Util.h"
#include "latticetester/Rank1Lattice.h"
#include "latticetester/BasisConstruction.h"
#include "latticetester/Coordinates.h"

typedef NTL::vector<Real> RealVec;

using namespace LatticeTester;

//Int m(101);      // Modulus m = 101
Int m(1021);     // Modulus m = 1021
//Int m(1048573);  // Modulus m = 1048573 (prime number near 2^{20})
Int a(33);       // The LCG multiplier
const long dim = 5;  // Dimension of lattice.
const long dimProj = 3;  // Dimension of projection.

int main() {
    std::cout << "Types: " << strFlexTypes << "\n";
    std::cout << "TestBasisConstructionSmall \n\n";

    // All the IntMat objects are created in 5 dimensions, but we may use less.
    IntMat basis1, basis2, basisProj, basisDual, basisDualProj;
    basis1.SetDims(dim, dim);
    basis2.SetDims(dim, dim);
    basisDual.SetDims(dim, dim);
    basisProj.SetDims(dim, dimProj);
    basisDualProj.SetDims(dimProj, dimProj);
    Int sqlength;
    RealVec* sqlen = 0;

    // We construct a Korobov lattice in dim dimensions.
    Rank1Lattice<Int, Real> *korlat;
    korlat = new Rank1Lattice<Int, Real>(m, a, dim, true, true);
    korlat->buildBasis(dim);
    basis1 = korlat->getBasis();
    // copy(korlat->getBasis(), basis1);  // This initial basis is triangular.
    std::cout << "Initial Korobov lattice basis (triangular) = \n" << basis1 << "\n";
    ProdScal<Int>(basis1[0], basis1[0], dim, sqlength);
    std::cout << "Square length of first basis vector: " << sqlength << "\n\n";

    // We apply LLL to reduce basis1.
    LLLConstruction0<IntMat, RealVec>(basis1, sqlen, 0.5, 0, 0);
    std::cout << "Basis after LLL with delta=0.5: \n" << basis1 << "\n";
    ProdScal<Int>(basis1[0], basis1[0], dim, sqlength);
    std::cout << "Square length of first basis vector: " << sqlength << "\n\n";

    LLLConstruction0<IntMat, RealVec>(basis1, sqlen, 0.99999, 0, 0);
    std::cout << "Basis after LLL with delta=0.99999: \n" << basis1 << "\n";
    ProdScal<Int>(basis1[0], basis1[0], dim, sqlength);
    std::cout << "Square length of first basis vector: " << sqlength << "\n\n";

    // We now transform basis1 to the upper-triangular basis2.
    // Note that after this, basis1 contains only garbage.
    upperTriangularBasis(basis1, basis2, m);
    std::cout << "After `upperTriangularBasis`: \n" << basis2 << "\n\n";
    // Then we compute the m-dual of basis2 and put it in basisDual.
    mDualUpperTriangular(basis2, basisDual, m);
    std::cout << "m-dual of upper-triangular basis: \n" << basisDual << "\n\n";
    // We reduce this basisDual with LLL.
    LLLConstruction0(basisDual, sqlen, 0.99999);
    std::cout << "m-dual basis after LLL with delta=0.99999: \n" << basisDual << "\n";
    ProdScal<Int>(basisDual[0], basisDual[0], dim, sqlength);
    std::cout << "Square length of first dual basis vector: " << sqlength << "\n\n";

    // We now investigate the projection over coordinates {1, 3, 5}.
    // We first insert those three coordinates one by one in `proj`.
    // We then compute a basis for this projection in two ways.
    Coordinates proj;
    proj.insert(1);
    proj.insert(3);
    proj.insert(5);
    std::cout << "Lattice projection over coordinates " << proj << ".\n";
    std::cout << "In the following basisProj matrices, we need 5 rows and 3 columns\n";
    std::cout << " to make the projection, then 3 rows and 3 columns for the basis.\n";
    std::cout << " When part of matrix is not used, it must be ignored.\n\n";

    // Basis construction with LLL.
    projectMatrix(basis2, basisProj, proj, dim);
    std::cout << "basisProj after projectMatrix (the generating vectors): \n" << basisProj << "\n";
    LLLBasisConstruction(basisProj, m, sqlen, 0.5, dim);
    std::cout << "Basis for this projection, with LLL: \n" << basisProj << "\n";

    // Basis construction with upper-triangular method, using 5 rows.
    projectionConstructionUpperTri(basis2, basisProj, proj, m, 5);
    std::cout << "Upper-triangular basis for this proj.: \n" << basisProj
            << "\n";

    // Use first three rows of `basisProj` basis matrix to construct an m-dual basis.
    mDualUpperTriangular(basisProj, basisDualProj, m, 3);
    std::cout << "Triangular basis for the m-dual of this projection: \n"
            << basisDualProj << "\n";
    LLLConstruction0(basisDualProj, sqlen, 0.99999, 3, 3);
    std::cout << "m-dual basis of proj after LLL with delta=0.99999: \n" << basisDualProj
            << "\n";
    ProdScal<Int>(basisDualProj[0], basisDualProj[0], dim, sqlength);
    std::cout << "Square length of first m-dual basis vector: " << sqlength << "\n\n";
    return 0;
}

