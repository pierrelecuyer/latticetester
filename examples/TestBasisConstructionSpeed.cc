/**
 * This example makes speed comparisons with the BasisConstruction functions,
 * with all five combinations of types.
 * See the Lattice Tester guide for more explanations.
 *
 * ****  Add more explanations later.  *****
 *
 *
 **/

#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/xdouble.h>
#include <NTL/ZZ.h>
#include <NTL/mat_ZZ.h>
#include "latticetester/FlexTypes.h"
#include "latticetester/Util.h"
#include "latticetester/Rank1Lattice.h"
#include "latticetester/BasisConstruction.h"

using namespace NTL;
using namespace LatticeTester;

// We cannot use Int or Real here, because they are not yet defined.
const long numSizes = 5; // Number of matrix sizes (choices of dimension).
const long dimensions[numSizes] = { 4, 6, 10, 20, 30 };
const long numMeth = 6;    // Number of methods to test, and their names.
std::string names[numMeth] = { "LLL5      ", "LLL9      ", "LLL99999  ",
                               "LLL99999-R", "UppTri    ", "mDualUT   "};
// We use `ctime` directly for the timings, to minimize overhead.
clock_t totalTime = clock(); // Global timer for total time.
clock_t timer[numMeth][numSizes];
double sumSq[numMeth][numSizes];  // Sums of square lengths, in Real type.
std::string stringTypes;  // To print the selected flexible types.

void printResults();  // Must be declared, because it has no parameters.

// Run a speed test for dim = dimensions[d], with given basis matrices.
template<typename Int, typename IntMat, typename Real>
void transformBases (Int m, long d, long dim, IntMat &basis1, IntMat &basis2,
        IntMat &basisdual) {
    CopyPartMat<IntMat> (basis2, basis1, dim, dim);  // Copy basis1 to basis2.
    clock_t tmp;
    NTL::vector<Real> sqlen; // Cannot be global variable because it depends on Real.
    sqlen.SetLength(1);

    // We apply LLL to basis1 with different values of `delta`, incrementally.
    // We start with delta=0.5, then continue with 0.9, then with 0.99999.
    tmp = clock();
    LLLConstruction0(basis2, 0.5, dim, dim, &sqlen);
    timer[0][d] += clock() - tmp;
    sumSq[0][d] += conv<double>(sqlen[0]);

    // We continue the LLL process with a larger `delta`.
    tmp = clock();
    LLLConstruction0(basis2,  0.9, dim, dim, &sqlen);
    timer[1][d] += clock() - tmp;
    sumSq[1][d] += conv<double>(sqlen[0]);

    tmp = clock();
    LLLConstruction0(basis2, 0.99999, dim, dim, &sqlen);
    timer[2][d] += clock() - tmp;
    sumSq[2][d] += conv<double>(sqlen[0]);

    // Here we restart LLL from the initial triangular basis, with delta=0.99999.
    CopyPartMat (basis2, basis1, dim, dim);  // Copy basis1 to basis2.
    tmp = clock();
    LLLConstruction0(basis2, 0.99999, dim, dim, &sqlen);
    timer[3][d] += clock() - tmp;
    sumSq[3][d] += conv<double>(sqlen[0]);

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
void testLoopResize(NTL::ZZ mm, long numRep) {
    long d, dim;
    Int m = conv<Int>(mm);
    Int a;       // The LCG multiplier
    IntMat basis1, basis2, basisdual;
    Rank1Lattice<Int, Real> *korlat;    // Will be a Korobov lattice.
    std::cout << "TestBasisConstructionSpeed with m = " << m << "\n";
    std::cout << "Results for `testLoopResize` (many objects are created or resized)\n";
    for (d = 0; d < numSizes; d++)      // Reset timers and sums.
        for (int64_t meth = 0; meth < numMeth; meth++) {
           timer[meth][d] = 0;   sumSq[meth][d] = 0.0;  // NTL::conv<Real>(0.0);
        }
    totalTime = clock(); // Global timer for total time.
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
            transformBases<Int, IntMat, Real>(m, d, dim, basis1, basis2, basisdual);
            delete korlat;
        }
    }
    printResults();
    basis1.kill();  // Since we create objects repeatedly,
    basis2.kill();  // it is a good idea to release the memory when we are done.
    basisdual.kill();
}

// In this testing loop, we try to minimize the creation of objects.
// The `IntMat` and `Rank1Lattice` objects are created only once.
template<typename Int, typename IntMat, typename Real>
void testLoopNoResize(NTL::ZZ mm, long numRep) {
    long d, dim;  // Index of dimension.
    Int m = conv<Int>(mm);
    Int a;
    IntMat basis1, basis2, basisdual;
    long maxdim = dimensions[numSizes - 1];   // Maximum dimension
    basis1.SetDims(maxdim, maxdim); // Will be initial triangular basis.
    basis2.SetDims(maxdim, maxdim); // Will be LLL-reduced basis.
    basisdual.SetDims(maxdim, maxdim);  // m-dual basis.
    Rank1Lattice<Int, Real> *korlat;  // We create a single Korobov lattice object.
    korlat = new Rank1Lattice<Int, Real>(m, maxdim, true, false);
    std::cout << "TestBasisConstructionSpeed with m = " << m << "\n";
    std::cout << "Results for `testLoop No Resize`\n";

    for (d = 0; d < numSizes; d++)   // Reset accumulators.
        for (int64_t meth = 0; meth < numMeth; meth++) {
            timer[meth][d] = 0;     sumSq[meth][d] = 0.0;  // NTL::conv<Real>(0.0);
        }
    totalTime = clock(); // Global timer for total time.
    for (int64_t r = 0; r < numRep; r++) {
        a = (m / 5 + 17 * r) % m;   // The multiplier we use for this rep.
        korlat->seta(a);
        for (d = 0; d < numSizes; d++) {  // Each matrix size
            dim = dimensions[d]; // The corresponding dimension.
            korlat->buildBasis(dim);
            // std::cout << "a = " << a << ",  dim = " << dimensions[d] << "\n";
            CopyPartMat<IntMat> (basis1, korlat->getBasis(), dim, dim); // Triangular basis.
            transformBases<Int, IntMat, Real>(m, d, dim, basis1, basis2, basisdual);
        }
    }
    printResults();
    basis1.kill();
    basis2.kill();
    basisdual.kill();
}


// This function runs the two types of test loops.
template<typename Int, typename IntMat, typename Real>
void testTwoLoops(NTL::ZZ mm, long numRep) {
   strTypes<Int, Real>(stringTypes);  // Functions from FlexTypes
   std::cout << "****************************************************\n";
   std::cout << "Types: " << stringTypes << "\n\n";
   testLoopResize<Int, IntMat, Real>(mm, numRep);
   testLoopNoResize<Int, IntMat, Real>(mm, numRep);
}

void printResults() {
    long d;
    std::cout << "Timings for different methods, in basic clock units (microseconds) \n\n";
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
   // Here, Int and Real are not yet defined.
   //Int m = 1021;  // Prime modulus near 2^{10}
    NTL::ZZ mm(1048573);  // Prime modulus near 2^{20}
   //Int m = 1073741827;  // Prime modulus near 2^{30}
   //Int m = 1099511627791;  // Prime modulus near 2^{40}
   //Int m = 1125899906842597;  // Prime modulus near 2^{50}
    //long mlong = mm;
    long numRep = 1000;   // Number of replications (multipliers) for each case.

    testTwoLoops<long, NTL::matrix<long>, double>(mm, numRep);
    testTwoLoops<NTL::ZZ, NTL::matrix<NTL::ZZ>, double>(mm, numRep);
    testTwoLoops<NTL::ZZ, NTL::matrix<NTL::ZZ>, xdouble>(mm, numRep);
    testTwoLoops<NTL::ZZ, NTL::matrix<NTL::ZZ>, quad_float>(mm, numRep);
    testTwoLoops<NTL::ZZ, NTL::matrix<NTL::ZZ>, NTL::RR>(mm, numRep);
    return 0;
};

