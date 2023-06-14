/**
 * This example shows how to use `LatticeTester` to implement a figure of merit that
 * is a weighted sum, maximum, minimum or average of a measure on projections on
 * a lattice. It shows how to use the `Weights` classes, how to build
 * projections of a basis, and how to normalize a computation. Typically,
 * when building figures of merit, measures need to be rescaled to the same
 * interval to be compared with one another, this is what is called
 * normalization.
 * 
 * This program computes a simple spectral test a) on all projections of a lattice
 * in 10 dimensions and b) on all two- and three-dimensional projections, 
 * normalizes it between 0 and 1 and then takes the minimal value observed as a figure 
 * of merit for that lattice. 
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
  Int m(1021);     // Modulus m = 1021
  Int a; 
  a = m /5;
  long dim;
  dim = max_dim;
  //Build a Korobov lattice
  Rank1Lattice<Int, Real> *lat;
  lat = new Rank1Lattice<Int, Real>(m, a, dim);
  lat->buildBasis(dim);   
  ProjConstructType pctype;
  //Choose which projection construction is used LLLPROJ or UPPERTRIPROJ
  pctype = UPPERTRIPROJ;
  
  IntLattice<Int, Real> *proj; // Another IntLattice to store projections
  
  double merit1 = 1.0, merit2 = 1.0;

  IntMat projBasis; //Matrix for basis of projection
  
  Reducer<Int, Real> *red = new Reducer<Int, Real>(max_dim);   

  // The variables specific to the construction of a figure of merit
   WeightsUniform weights(1.0); // This just puts a weight of 1 to everything
   
  // EXAMPLE 1: LOOK AT PROJECTIONS TO ALL COORDINATES //
   
  // Creates an iterator over a set of projections.
  // It contains all sets of size at least min_dim + 2 and at most max_dim. 
  // within the coordinates {min_dim, min_dim + 1, ..., max_dim - 1}
  // Example (2, 4, 0, 3) would be all subsets of {0, 1, 2, 3} of cardinaly between 2 and 4.
   //CoordinateSets::FromRanges coord(min_dim+2, max_dim, min_dim, max_dim - 1);     
   CoordinateSets::FromRanges allcoord(min_dim+2, max_dim, min_dim, max_dim - 1);      
   
      
  // Loop over the selected set of projections.
  for(auto it = allcoord.begin(); it != allcoord.end(); it++){
    // Computing a basis for the projection, using chosen projection type 	  
	BasisConstruction<Int>::projectionConstruction(lat->getBasis(), projBasis, *it, m, pctype);	
	proj = new IntLattice<Int, Real> (projBasis, m, projBasis.NumCols());
    
    //! Computing the shortest vector in the lattice spanned by matrix
    proj->updateVecNorm();
    proj->sort(0);
    
    red->redBKZ(proj->getBasis());
  
    red->shortestVector(*proj);
    double shortest = NTL::conv<double>(red->getMinLength());

    // Instantiating the normalizers
    // The prefered way of doing this is described in Normalizer documentation
    double log_density=(double)(-log(abs(NTL::determinant(proj->getBasis()))));
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
  std::cout << "CASE 1: Look at all possible subsets with cardinality > 1:" << "\n";
  std::cout << "Figure of merit with BestLat: " << merit1 << "\n";
  std::cout << "Figure of merit with BestBound: " << merit2 << std::endl;
  std::cout << "Figures of merit are different for different normalizers,"
               " weights and projections choices\n";
  std::cout << "\n";

  // ============================================================================

  // EXAMPLE 2: LOOK AT TWO- AND THREE-DIMENSIONAL PROJECTIONS ONLY
  
  // Reset figures of merit and basis of projection
  merit1 = 1.0, merit2 = 1.0;
  projBasis.resize(0,0);
  
  // Pick projections
  CoordinateSets::FromRanges ldcoord(2, 3, min_dim, max_dim - 1);  
  
  // Loop over the selected set of projections.
  for(auto it = ldcoord.begin(); it != ldcoord.end(); it++){
    // Computing a basis for the projection, using chosen projection type.
    BasisConstruction<Int>::projectionConstruction(lat->getBasis(), projBasis, *it, m, pctype);
	proj = new IntLattice<Int, Real> (projBasis, m, projBasis.NumCols());
    //! Computing the shortest vector in the lattice spanned by matrix
    proj->updateVecNorm();
    proj->sort(0);
    
    red->redBKZ(proj->getBasis());
  
    red->shortestVector(*proj);
    double shortest = NTL::conv<double>(red->getMinLength());

    // Instantiating the normalizers
    // The prefered way of doing this is described in Normalizer documentation
    double log_density=(double)(-log(abs(NTL::determinant(proj->getBasis()))));
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
  std::cout << "CASE 2: Only look at all two- and three-dimensional projections:" << "\n";
  std::cout << "Figure of merit with BestLat: " << merit1 << "\n";
  std::cout << "Figure of merit with BestBound: " << merit2 << std::endl;
  std::cout << "Figures of merit are different for different normalizers,"
               " weights and projections choices\n";
  


  return 0;
}
