/**
 * This program computes the various bounds in the `Normalizer` subclasses and display
 * them in two tables, for the L1 and L2 norms, for comparison.  It also prints a table that
 * gives the center densities \f$\delta_t\f$ that correspond to the given estimates of \f$gamma_t\f$.
 * See the guide for more details.
 */

#include <iostream>
#include <vector>
#include <NTL/vector.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>

#include "latticetester/Normalizer.h"
#include "latticetester/NormaBestLat.h"
#include "latticetester/NormaLaminated.h"
#include "latticetester/NormaMinkHlaw.h"
#include "latticetester/NormaBestUpBound.h"
#include "latticetester/NormaRogers.h"
#include "latticetester/NormaMinkL1.h"


using namespace LatticeTester;

int main() {

   std::string typeBound[] = {
         "BestLat   ",
         "Laminated ",
         "MinkHlaw  ",
         "BestUp    ",
         "Rogers    ",
         "MinkL1    "
   };
   const long numTypes = 6;   // Number of types of bounds.
   const long maxDim = 48;
   // double valGam[numTypes][maxDim];  // The values of gamma_t

   std::cout << "TestNormBounds: Approximations of the Hermite constants \n\n";
   NormaBestLat     norBest(1.0, 48);
   NormaLaminated   norLam(1.0, 48);
   NormaMinkHlaw    norMH(1.0, 48);
   NormaBestUpBound norBestUp(1.0, 48);
   NormaRogers      norRog(1.0, 48);
   NormaMinkL1      norML1(1.0, 48);

   std::cout << "For the L2 norm (the \\gamma_t): \n\n";
   std::cout << " Dim";
   for (long d = 0; d < numTypes-1; d++)
       std::cout << " &    " << std::setw(10) << typeBound[d];
   std::cout << " \\\\ \n";
   for (int64_t t = 1; t <= maxDim; t++) {
      std::cout << std::setw(4) << std::setprecision(2) << t << " & ";
      std::cout << std::setw(13) << std::setprecision(10) << norBest.getGamma(t) << " & ";
      std::cout << std::setw(13) << std::setprecision(10) << norLam.getGamma(t) << " & ";
      std::cout << std::setw(13) << std::setprecision(10) << norMH.getGamma(t) << " & ";
      std::cout << std::setw(13) << std::setprecision(10) << norBestUp.getGamma(t) << " & ";
      std::cout << std::setw(13) << std::setprecision(10) << norRog.getGamma(t) << " \\\\ ";
      std::cout << "\n";
   }
   std::cout << "\n\n";

   std::cout << "For the L1 norm (the \\gamma_t^{(1)} = t * \\gamma_t): \n\n";
   std::cout << " Dim";
   for (long d = 0; d < numTypes-1; d++)
       std::cout << " &    " << std::setw(10) << typeBound[d];
   std::cout << " & " << std::setw(10) << typeBound[5] << "\\\\ \n";
   for (int64_t t = 1; t <= maxDim; t++) {
      std::cout << std::setw(4) << std::setprecision(2) << t << " & ";
      std::cout << std::setw(13) << std::setprecision(10) << norBest.getGamma(t) * t << " & ";
      std::cout << std::setw(13) << std::setprecision(10) << norLam.getGamma(t)*t << " & ";
      std::cout << std::setw(13) << std::setprecision(10) << norMH.getGamma(t)*t << " & ";
      std::cout << std::setw(13) << std::setprecision(10) << norBestUp.getGamma(t)*t << " & ";
      std::cout << std::setw(13) << std::setprecision(10) << norRog.getGamma(t)*t << " & ";
      std::cout << std::setw(13) << std::setprecision(10) << norML1.getGamma(t) << " \\\\";
      std::cout << "\n";
   }
   std::cout << "\n\n";

   std::cout << "Center density \\delta_t that corresponds to the given estimate of \\gamma_t: \n\n";
   std::cout << " Dim";
   for (long d = 0; d < numTypes-1; d++)
       std::cout << " &    " << std::setw(10) << typeBound[d];
   std::cout << " \\\\ \n";
   for (int64_t t = 1; t <= maxDim; t++) {
//      for (auto type : { norBest, norLam, norMH, norBestUp, norRog, norML1 })
//         std::cout << std::setw(10) << std::setprecision(10) << type.getGamma(t) << " ";
      std::cout << std::setw(4) << std::setprecision(2) << t << " & ";
      std::cout << std::setw(13) << std::setprecision(10) << pow(norBest.getGamma(t)/4.0, t/2.0) << " & ";
      std::cout << std::setw(13) << std::setprecision(10) << pow(norLam.getGamma(t)/4.0, t/2.0) << " & ";
      std::cout << std::setw(13) << std::setprecision(10) << pow(norMH.getGamma(t)/4.0, t/2.0) << " & ";
      std::cout << std::setw(13) << std::setprecision(10) << pow(norBestUp.getGamma(t)/4.0, t/2.0) << " & ";
      std::cout << std::setw(13) << std::setprecision(10) << pow(norRog.getGamma(t)/4.0, t/2.0) << " \\\\ ";
      std::cout << "\n";
   }

   return 0;
}
