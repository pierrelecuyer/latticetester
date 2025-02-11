// File `tesMatrixCreationSpped`

#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/xdouble.h>
#include <NTL/ZZ.h>
#include <NTL/mat_ZZ.h>
#include "latticetester/FlexTypes.h"
#include "latticetester/Util.h"
#include "latticetester/Rank1Lattice.h"
#include "latticetester/BasisConstruction.h"

/**
 * This program examines what is the fastest way to use IntMat objects
 * when many many matrices of different dimensions are used througout the program.
 * We can either just use one big matrix object of the biggest possible dimension
 * and only fill the upper left part of the big matrix to store the values of 
 * smaller matrices. Alternatively we can just use one matrix object and resize it.
 * As LatticeTester releis on NTL this means that a new object is created whenever
 * there is an actual resizing, i.e. the dimension of the object really changes.
 * The last alternative would be to create an IntMat object whenever we want to 
 * use one.
 */
using namespace NTL;
using namespace LatticeTester;

// The types Int and Real are will be passed as template parameters from the `main`.
const int64_t maxNumSizes = 7; // Number of matrix sizes (choices of dimensions).
const int64_t dimensions[maxNumSizes] = { 5, 10, 15, 20, 25, 30, 35 };
const int64_t numMeth = 5; // Number of differnt methods to be tested.
const int64_t numChanges = 5; // Stores the sqrt of the numbers of manipulations applied to each matrix
clock_t timer[numMeth][maxNumSizes];
std::string names[numMeth] = { "One Matrix, no Resize       ", "One Matrix, Resize if nec.  ", "One Matrix, numRep Resize   ", "One Matrix, num Rows const  ", "numRep Matrices             "};

void printTable(int64_t numSizes) {
  
   int64_t d;
   std::cout << "Dimension:                 ";
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
}


// Testing loop. The `numRepetitions' number of repetitions, 'numSizes' of different matrix sizes and 'numFlex'^2 number of matrix changes
// In the seconde loop a fixed number of matrix changes is applied according to 'numFlex'^2
template<typename Int, typename Real>
void testLoop(const int64_t & numRepetitions, const int64_t & numSizes, const int64_t & numFlex) {
   clock_t tmp;
   int64_t i, j, k, d;
      
   std::string stringTypes;  // To print the selected flexible types.
   strTypes<Int, Real>(stringTypes);  // Functions from FlexTypes   
   std::cout << "**************************************************************\n";
   std::cout << "TestMatrixCreationSpeed \n";
   std::cout << "Types: " << stringTypes << "\n";
   std::cout << "Number of repetitions: " << numRepetitions << "\n";
   std::cout << "With flexible number of changes" << "\n";
   std::cout << "Timings for the different tasks, in basic clock units (microseconds): \n";    
   
   // The first part of the test examines the case when the upp
   
   // Method 1: Only create matrix object once, set the dimension once and fill all entries of the matrix
   IntMat A;
   A.SetDims(dimensions[numSizes-1], dimensions[numSizes-1]);
   for (d = 0; d < numSizes; d++) {
      tmp = clock();
      for (i = 0; i < numRepetitions; i++) {
         for (j = 0; j < numFlex; j++) {
            for (k = 0; k < numFlex; k++) {
               A[j % dimensions[d]][k % dimensions[d]] = i+j+k;
            }
         }
      }
      timer[0][d] = clock() - tmp;
   }
   
      
   // Method 2: Only create matrix object once, resize it if necessary and fill the entries
   IntMat B;   
   for (d = 0; d < numSizes; d++) {
      tmp = clock(); 
      for (i = 0; i < numRepetitions; i++) {  
         B.SetDims(dimensions[d], dimensions[d]);   
         for (j = 0; j < numFlex; j++) {
            for (k = 0; k < numFlex; k++) {
               B[j % dimensions[d]][k % dimensions[d]] = i+j+k;
            }
         }
      }
      timer[1][d] = clock() - tmp;
   }
   
   // Method 3: Only create matrix object once, always resize it and fill the entries  
   for (d = 0; d < numSizes; d++) {
      tmp = clock(); 
      for (i = 0; i < numRepetitions; i++) {  
         B.SetDims(0, 0);
         B.SetDims(dimensions[d], dimensions[d]);   
         for (j = 0; j < numFlex; j++) {
            for (k = 0; k < numFlex; k++) {
               B[j % dimensions[d]][k % dimensions[d]] = i+j+k;
            }
         }
      }
      timer[2][d] = clock() - tmp;
   }
   
   // Method 4: Only create matrix object once, always resize the number of rows but not the number of columns  
   for (d = 0; d < numSizes; d++) {
      tmp = clock(); 
      for (i = 0; i < numRepetitions; i++) {  
         B.SetDims(0, dimensions[numSizes - 1]);
         B.SetDims(dimensions[d], dimensions[numSizes - 1]);   
         for (j = 0; j < numFlex; j++) {
            for (k = 0; k < numFlex; k++) {
               B[j % dimensions[d]][k % dimensions[d]] = i+j+k;
            }
         }
      }
      timer[3][d] = clock() - tmp;
   }
   
   // Method 5: Create the matrix again and again and fill the entries
   
   for (d = 0; d < numSizes; d++) {
      tmp = clock();
      for (i = 0; i < numRepetitions; i++) {    
      IntMat C;
      C.SetDims(dimensions[d], dimensions[d]);     
         for (j = 0; j < numChanges; j++) {
            for (k = 0; k < numChanges; k++) {
               C[j % dimensions[d]][k % dimensions[d]] = i+j+k;
            }
         }
      }
      timer[4][d] = clock() - tmp;
   }
   
   
   printTable(numSizes);   
   
   // To better understand what is the influence of the for loop to refill the entries,
   // we also look at the case where the number of changes is not set by the user via
   // a local parameter but by a globar variable. This has significant effect.


   std::cout << "**************************************************************\n";
   std::cout << "TestMatrixCreationSpeed \n";
   std::cout << "Types: " << stringTypes << "\n";
   std::cout << "Number of repetitions: " << numRepetitions << "\n";
   std::cout << "With number of changes fixed in advance" << "\n";
   std::cout << "Timings for the different tasks, in basic clock units (microseconds): \n";    
   
   IntMat D;
   D.SetDims(dimensions[numSizes-1], dimensions[numSizes-1]);
   for (d = 0; d < numSizes; d++) {
      tmp = clock();
      for (i = 0; i < numRepetitions; i++) {
         for (j = 0; j < numChanges; j++) {
            for (k = 0; k < numChanges; k++) {
               D[j % dimensions[d]][k % dimensions[d]] = i+j+k;
            }
         }
      }
      timer[0][d] = clock() - tmp;
   }
   
      
   // Method 2: Only create matrix object once, resize it if necessary and fill the entries
   IntMat E;   
   for (d = 0; d < numSizes; d++) {
      tmp = clock(); 
      for (i = 0; i < numRepetitions; i++) {
         E.SetDims(dimensions[d], dimensions[d]);   
         for (j = 0; j < numChanges; j++) {
            for (k = 0; k < numChanges; k++) {
               E[j % dimensions[d]][k % dimensions[d]] = i+j+k;
            }
         }
      }
      timer[1][d] = clock() - tmp;
   }
   
   
   // Method 3: Only create matrix object once, always resize it and fill the entries
   for (d = 0; d < numSizes; d++) {
      tmp = clock(); 
      for (i = 0; i < numRepetitions; i++) {
         E.SetDims(0, 0);
         E.SetDims(dimensions[d], dimensions[d]);   
         for (j = 0; j < numChanges; j++) {
            for (k = 0; k < numChanges; k++) {
               E[j % dimensions[d]][k % dimensions[d]] = i+j+k;
            }
         }
      }
      timer[2][d] = clock() - tmp;
   }
   
   // Method 4: Only create matrix object once, always resize it and fill the entries
   for (d = 0; d < numSizes; d++) {
      tmp = clock(); 
      for (i = 0; i < numRepetitions; i++) {
         E.SetDims(0, dimensions[numSizes - 1]);
         E.SetDims(dimensions[d], dimensions[numSizes - 1]);   
         for (j = 0; j < numChanges; j++) {
            for (k = 0; k < numChanges; k++) {
               E[j % dimensions[d]][k % dimensions[d]] = i+j+k;
            }
         }
      }
      timer[3][d] = clock() - tmp;
   }
   
   // Method 5: Create the matrix again and again and fill the entries
   for (d = 0; d < numSizes; d++) {
      tmp = clock();
      for (i = 0; i < numRepetitions; i++) {    
      IntMat F;
      F.SetDims(dimensions[d], dimensions[d]);     
         for (j = 0; j < numChanges; j++) {
            for (k = 0; k < numChanges; k++) {
               F[j % dimensions[d]][k % dimensions[d]] = i+j+k;
            }
         }
      }
      timer[4][d] = clock() - tmp;
   }
   
   printTable(numSizes);   
}

int main() {

// Here, `Int` and `Real` are not yet defined, they will be passed as template parameters.
const int64_t numSizes = 7;
const int64_t numRepetitions = 100000;
const int64_t numberOfChanges = 5;

testLoop<int64_t, double>(numRepetitions, numSizes, numberOfChanges);
return 0;
}

