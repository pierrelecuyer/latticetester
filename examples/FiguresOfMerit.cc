/**
 * This example shows how to use `LatticeTester` to implement a figure of merit that
 * is a weighted sum, maximum, minimum or average of a measure on projections on
 * a lattice. It shows how to use the `Weights` classes, how to build
 * projections of a basis, and how to normalize a computation. Typically,
 * when building figures of merit, measures need to be rescaled to the same
 * interval to be compared with one another, this is what is called
 * normalization.
 * 
 * This program computes a simple spectral test on all projections of a lattice
 * in 10 dimensions, normalizes it between 0 and 1 and then takes the minimal
 * value observed as a figure of merit for that lattice.
 * It outputs only the figure of merit for the lattice for two
 * different normalizers. This is not really interesting in itself, hence it is
 * not included here (???). To get interesting informations on a figure of merit like
 * this one, it would be possible to store an few of the worst projections and
 * print them after the test. 
 * */

//#define NTL_TYPES_CODE 2
#define TYPES_CODE  ZD
#include <iostream>
#include <cstdint>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>

#include "latticetester/FlexTypes.h"
#include "latticetester/EnumTypes.h"
#include "latticetester/BasisConstruction.h"
#include "latticetester/Util.h"
#include "latticetester/ParamReader.h"
#include "latticetester/IntLattice.h"
#include "latticetester/Rank1Lattice.h"
#include "latticetester/Reducer.h"

// Application specific headers
#include "latticetester/NormaBestLat.h"
#include "latticetester/NormaBestBound.h"
#include "latticetester/CoordinateSets.h"
#include "latticetester/WeightsUniform.h"
#include "latticetester/BasisConstruction.h"

using namespace LatticeTester;

int main() {
  //! Reading a matrix to use the normalizer class
  int64_t min_dim = 0, max_dim = 10;
  ParamReader<Int, Real> reader("./44matrixEx.dat");
  reader.getLines();
  IntMat matrix(max_dim,max_dim);
  uint64_t ln = 0;
  reader.readBMat(matrix, ln, 0, max_dim);
  Int m;
  determinant(m, matrix);
  int64_t n;
  NTL::conv(n,m);
  IntLattice<Int, Real> lat(matrix, m, max_dim);   // New IntLattice.
  double merit1 = 1.0, merit2 = 1.0;

  // The variables specific to the construction of a figure of merit
   WeightsUniform weights(1.0); // This just puts a weight of 1 to everything
   BasisConstruction<Int> constructor; // Computes projections basis
   IntLattice<Int, Real> proj(m, max_dim);    // Another IntLattice to store projections

  // Creates an iterator over a set of projections.
   CoordinateSets::FromRanges coord(min_dim+1, max_dim, min_dim, max_dim-1);
   
   Reducer<Int, Real> *red = new Reducer<Int,Real>(max_dim);
   //Reducer<Int, Real> red(max_dim);
   
   
  // Loop over the selected set of projections.
  for(auto it = coord.begin(); it != coord.end(); it++){
    // Computing a basis for the projection, using LLL.     
    constructor.projectionConstructionLLL(lat, proj, *it, n);
    
    //! Computing the shortest vector in the lattice spanned by matrix
    proj.updateVecNorm();
    proj.sort(0);
    //red.redBKZ(proj->getBasis());
  
//    std::string ch("cholesky");
//    red.shortestVector(L2NORM,ch);
    red->shortestVector(proj);
    double shortest = NTL::conv<double>(red->getMinLength());

    // Instantiating the normalizers
    // The prefered way of doing this is described in Normalizer documentation
    double log_density=(double)(-log(abs(NTL::determinant(proj.getBasis()))));
    Normalizer* norma = new NormaBestLat(log_density, max_dim);

    // Computing the figure of merit for this projection
    double merit = weights.getWeight(*it) * shortest/norma->getBound((*it).size());
    // Testing if it is the minimum as of now
    if (merit < merit1) merit1 = merit;
    delete norma;
    norma = new NormaBestBound(log_density, max_dim);
    merit = weights.getWeight(*it) * shortest/norma->getBound((*it).size());
    if (merit < merit2) merit2 = merit;
    delete norma;
  }

  //! Printing the results in three simple lines
  std::cout << "Figure of merit with BestLat: " << merit1 << "\n";
  std::cout << "Figure of merit with BestBound: " << merit2 << std::endl;
  std::cout << "Figures of merit are different for different normalizers,"
               " weights and projections choices\n";
  

  return 0;
}
