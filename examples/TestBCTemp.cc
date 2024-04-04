/**
 * This example illustrates the usage of the BasisConstruction module.
 *
 * TO REWRITE WHEN DONE AND STABLE:   *****
 *
 * This one is to test using templates to avoid recompiling when
 * changing the types codes.
 *
 **/

// Select the flexible types Int and Real here.
#define TYPES_CODE  LD     // Int == int64_t
//#define TYPES_CODE  ZD     // Int == ZZ, Real = double
//#define TYPES_CODE  ZQ     // Int == ZZ, Real = quad_float
//#define TYPES_CODE  ZX     // Int == ZZ, Real = xdouble
//#define TYPES_CODE  ZR     // Int == ZZ, Real = RR

#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>

#include "latticetester/FlexTypes.h"    // This defines Int and Real
#include "latticetester/EnumTypes.h"
#include "latticetester/Util.h"
#include "latticetester/Chrono.h"
#include "latticetester/IntLattice.h"
#include "latticetester/Rank1Lattice.h"
#include "latticetester/BasisConstruction.h"

using namespace LatticeTester;


template<typename Int, typename IntMat, typename Real>
void testBC


//Int m(1021);     // Modulus m = 1021
Int m(1048573);  // Prime modulus near 2^{20}
//Int m(1073741827);  // Prime modulus near 2^{30}
//Int m(1099511627791);  // Prime modulus near 2^{40}
//Int m(1125899906842597);  // Prime modulus near 2^{50}
Int a;       // The LCG multiplier

const long numSizes = 5;    // Number of matrix sizes (choices of dimension).
//const long dimensions[numSizes] = { 10, 15, 20, 25, 30 };
const long dimensions[numSizes] = { 4, 6, 10, 20, 30 };
long maxdim = dimensions[numSizes - 1];   // Maximum dimension
const long numMeth = 6;    // Number of methods to test, and their names.
std::string names[numMeth] = { "LLL5      ", "LLL9      ", "LLL99999  ",
                               "LLL99999-R", "UppTri    ", "mDualUT   "};

// We use `ctime` directly for the timings, to minimize overhead.
clock_t tmp;
clock_t totalTime;  // Global timer for total time.
clock_t timer[numMeth][numSizes];
Real sumSq[numMeth][numSizes];  // Sums of square lengths, in Real type.
NTL::vector<Real> sqlen;        // This vector is created in the main.

// Run a speed test for dim = dimensions[d], with given basis matrices.
template<typename Int, typename IntMat, typename Real>
void transformBases (long d, long dim, IntMat &basis1, IntMat &basis2,
        IntMat &basisdual) {
    CopyPartMat (basis2, basis1, dim, dim);  // Copy basis1 to basis2.

    // We apply LLL to basis1 with different values of `delta`, incrementally.
    // We start with delta=0.5, then continue with 0.9, then with 0.99999.
    tmp = clock();
    LLLConstruction0(basis2, 0.5, dim, dim, &sqlen);
    timer[0][d] += clock() - tmp;
    sumSq[0][d] += sqlen[0];

    // We continue the LLL process with a larger `delta`.
    tmp = clock();
    LLLConstruction0(basis2,  0.9, dim, dim, &sqlen);
    timer[1][d] += clock() - tmp;
    sumSq[1][d] += (sqlen)[0];

    tmp = clock();
    LLLConstruction0(basis2, 0.99999, dim, dim, &sqlen);
    timer[2][d] += clock() - tmp;
    sumSq[2][d] += (sqlen)[0];

    // Here we restart LLL from the initial triangular basis, with delta=0.99999.
    CopyPartMat (basis2, basis1, dim, dim);  // Copy basis1 to basis2.
    tmp = clock();
    LLLConstruction0(basis2, 0.99999, dim, dim, &sqlen);
    timer[3][d] += clock() - tmp;
    sumSq[3][d] += (sqlen)[0];

    // We now construct an upper-triangular basis from basis2 into basis1.
    tmp = clock();
    upperTriangularBasis(basis2, basis1, m, dim, dim);
    timer[4][d] += clock() - tmp;

    // We compute an m-dual basis to basis1.
    tmp = clock();
    mDualUpperTriangular(basis1, basisdual, m, dim);
    timer[5][d] += clock() - tmp;
}

// In this testing loop, new `Rank1Lattice` objects are created
// and the  `IntMat` matrices are resized inside the loop.
template<typename Int, typename IntMat, typename Real>
void testLoopResize(long numRep) {
    long d, dim;
    IntMat basis1, basis2, basisdual;
    Rank1Lattice<Int, Real> *korlat;    // Will be a Korobov lattice.
    for (d = 0; d < numSizes; d++)      // Reset timers and sums.
        for (int64_t meth = 0; meth < numMeth; meth++) {
            timer[meth][d] = 0;   sumSq[meth][d] = 0.0;  // NTL::conv<Real>(0.0);
        }
    totalTime = clock();
    for (int64_t r = 0; r < numRep; r++) {
        a = (m / 5 + 17 * r) % m;   // The multiplier we use for this rep.
        //for (d = 0; d < 2; d++) {  // Each matrix size
        for (d = 0; d < numSizes; d++) {  // Each matrix size
            dim = dimensions[d]; // The corresponding dimension.
            basis1.SetDims(dim, dim); // Will be initial triangular basis.
            basis2.SetDims(dim, dim); // Will be LLL-reduced basis.
            basisdual.SetDims(dim, dim);  // m-dual basis.
            korlat = new Rank1Lattice<Int, Real>(m, a, dim, true, false);
            korlat->buildBasis(dim);
            basis1 = korlat->getBasis();
            transformBases(d, dim, basis1, basis2, basisdual);
            delete korlat;
        }
    }
    basis1.kill();  // Since we create objects repeatedly,
    basis2.kill();  // it is a good idea to release the memory when we are done.
    basisdual.kill();
}

// In this testing loop, we try to minimize the creation of objects.
// The `IntMat` and `Rank1Lattice` objects are created only once.
template<typename Int, typename IntMat, typename Real>
static void testLoopNoResize(long numRep) {
    long d, dim;  // Index of dimension.
    IntMat basis1, basis2, basisdual;
    basis1.SetDims(maxdim, maxdim); // Will be initial triangular basis.
    basis2.SetDims(maxdim, maxdim); // Will be LLL-reduced basis.
    basisdual.SetDims(maxdim, maxdim);  // m-dual basis.
    Rank1Lattice<Int, Real> *korlat;  // We create a single Korobov lattice object.
    korlat = new Rank1Lattice<Int, Real>(m, maxdim, true, false);

    for (d = 0; d < numSizes; d++)   // Reset accumulators.
        for (int64_t meth = 0; meth < numMeth; meth++) {
            timer[meth][d] = 0;     sumSq[meth][d] = 0.0;  // NTL::conv<Real>(0.0);
        }
    totalTime = clock();
    for (int64_t r = 0; r < numRep; r++) {
        a = (m / 5 + 17 * r) % m;   // The multiplier we use for this rep.
        korlat->seta(a);
        for (d = 0; d < numSizes; d++) {  // Each matrix size
            dim = dimensions[d]; // The corresponding dimension.
            korlat->buildBasis(dim);
            // std::cout << "a = " << a << ",  dim = " << dimensions[d] << "\n";
            copy(korlat->getBasis(), basis1, dim, dim); // Triangular basis.
            transformBases(d, dim, basis1, basis2, basisdual);
        }
    }
}

static void printResults() {
    long d;
    std::cout << " dim:    ";
    for (d = 0; d < numSizes; d++)
        std::cout << std::setw(8) << dimensions[d] << "  ";
    std::cout << "\n\n";
    for (int meth = 0; meth < numMeth; meth++) {
        std::cout << names[meth] << " ";
        for (d = 0; d < numSizes; d++)
            std::cout << std::setw(9) << timer[meth][d] << " ";
        std::cout << "\n";
    }
    std::cout << "\n";
    std::cout << "Sums of square lengths of shortest basis vector";
    std::cout << " (must be the same across all implementations):\n";
    std::cout << " dim:    ";
    for (d = 0; d < numSizes; d++)
        std::cout << std::setw(13) << dimensions[d] << "  ";
    std::cout << "\n\n";
    for (int meth = 0; meth < numMeth-2; meth++) {
        std::cout << names[meth] << "  ";
        for (d = 0; d < numSizes; d++)
            std::cout << std::setw(14) << sumSq[meth][d] << " ";
        std::cout << "\n";
    }
    std::cout << "\n";
    std::cout << "Total time: "
            << (double) (clock() - totalTime) / (CLOCKS_PER_SEC)
            << " seconds\n\n\n";
}

int main() {
    long numRep = 1000;   // Number of replications (multipliers) for each case.

    sqlen.SetLength(1);   // Done here because cannot be done in preamble.
    std::cout << "Types: " << strFlexTypes << "\n";
    // std::cout << "PrecisionType: " << prec << "\n\n";
    std::cout << "Results of TestBasisConstructionSpeed.cc with m = " << m << "\n";
    std::cout << "Timings for different methods, in basic clock units \n\n";
    testLoopResize(numRep);
    std::cout << "Timings for `testLoop Resize` (many objects are created or resized)\n";
    printResults();
    testLoopNoResize(numRep);
    std::cout << "Timings for `testLoop No Resize`\n";
    printResults();
}

