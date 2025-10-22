//
// Copyright (C) 2012-2023  The LatticeTester authors, under the supervision
// of Pierre L'Ecuyer at Universit� de Montr�al.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef LATTICETESTER_FIGUREOFMERITM_H
#define LATTICETESTER_FIGUREOFMERITM_H

#include <iostream>
#include <cstdint>
#include <float.h>

#include "NTL/tools.h"
#include "NTL/ZZ.h"
#include "NTL/vec_ZZ.h"
#include "NTL/mat_ZZ.h"
#include "NTL/LLL.h"
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/RR.h>

#include "latticetester/FlexTypes.h"
#include "latticetester/EnumTypes.h"
#include "latticetester/Util.h"
#include "latticetester/Chrono.h"
#include "latticetester/BasisConstruction.h"
#include "latticetester/IntLattice.h"
#include "latticetester/IntLatticeExt.h"
#include "latticetester/LLL_lt.h"
#include "latticetester/Rank1Lattice.h"
#include "latticetester/ReducerStatic.h"
#include "latticetester/ReducerBB.h"
#include "latticetester/Normalizer.h"
#include "latticetester/CoordinateSets.h"
#include "latticetester/Weights.h"
#include "latticetester/WeightsOrderDependent.h"

namespace LatticeTester {

/**
 * \class FigureOfMeritM
 *
 * This class provides tools to calculate the *figure of merit* (FOM)
 * \f[
 *    M_{t_1,\dots,t_d} = \min\left[ \min_{I\in S(t_1)} \frac{\ell_I}{ \omega_I\, \ell_I^*(\eta_I)},\;
 *    \min_{2\le s\le d}\, \min_{I\in S(s,t_s)} \frac{\ell_I}{\omega_I \,\ell_I^*(\eta_I)} \right]
 * \f]
 * for any given `IntLatticeExt` object.
 * This FOM is computed for the (rescaled) primal lattice and the m-dual is never used.
 * To compute the FOM for the m-dual, see the class `FigureOfMeritMDual`.
 * Just using the present class for the dualized primal lattice object will not give correct results in general,
 * because the projections of the m-dual lattice are not the same as the m-duals of the projections.
 *
 * In the formula, the projections in \f$S(t_1)\f$ are those over successive coordinates in up to \f$t_1\f$
 * dimensions, while those is @f$S(s,t_s)@f$ are projections over \f$s\f$ distinct coordinates
 * that are non necessarily successive and are all in the set \f$\{1,\dots,t_s\}\f$,
 * for each order @f$s > 1@f$.
 * There are two variants for the latter: the first (default) variant takes
 * @f$S(s,t_s)@f$ as just defined, and the other considers
 * only the set @f$S^{(1)}(s,t_s)@f$ of projections that contain coordinate 1.
 * The parameter `includeFirst` in the constructor determines which variant is taken:
 * the latter option is taken when `includeFirst` is set to `true`.
 * See the guide for more details and examples.
 *
 * The lengths of the shortest vectors in the projections can be calculated exactly by using the
 * BB algorithm after applying pre-reductions such as LLL or BKZ, or they can be just approximated
 * by the lengths of the shortest basis vector after applying the pre-reductions only.
 * The latter is faster but may not give a shortest vector.
 * The norm used to compute the vector lengths is always the one inside the `IntLattice` object for which
 * we compute the FOM.
 *
 * The constructor has two template parameters to specify which `Int` and `Real` types are used.
 * It also requires the vector @f$(t_1,\dots,t_d)@f$,
 * a `Weights` object to specify the weights \f$\omega_I\f$,
 * a `Normalizer` object used to normalize the merit values,
 * a `ReducerBB` object used for the reduction in case we want to apply BB, and the
 * `includeFirst` parameter in case we want to change it to `true`.
 * The last two parameters are optional.
 * The BB is applied if and only if a (nonzero) `ReducerBB` is given.
 * Otherwise, we just use static methods for the reduction and need no `Reducer`.
 * By default, the pre-reduction method is BKZ with `delta = 0.99999` and `blocksize = 10`.
 * To change these values and/or apply LLL, one can use `setBKZ` and/or `setLLL`.
 * The reductions are always applied in the order: LLL, BKZ, BB.
 * To remove LLL or BKZ, it suffices to set its `delta` parameter to 0.0.
 *
 * After a `FigureOfMeritM` has been created by the constructor, one can use `setTVector`
 * to set (or reset) the vector  @f$(t_1,\dots,t_d)@f$, `setWeights` to change the weights,
 * `setLLL` or `setBKZ` to change the LLL or BKZ parameters, etc.
 *
 * The function `computeMerit` computes the FOM for a given lattice.
 * The computation is stopped (early exit) as soon as we know that the value of the FOM
 * will be below the lower bound, whose value can be changed via `setLowBound`.
 * The default value is 0.
 * The functions `computeMeritSucc` compute the FOM only for the projections in the set \f$S(t_1)\f$.
 * One version does not use the vector `t`.
 */
template<typename Int, typename Real>
class FigureOfMeritM {

public:

   /**
    * The first constructor will call `setTVector (t, includeFirst)`,
    * then set the 'Weights', `Normalizer`, and `ReducerBB` to the given values.
    * See the text above for other default values.
    * The second constructor does not set the vector `t`, which is not used
    * for example when we only call `computeMeritSucc`.
    */
   FigureOfMeritM(const NTL::Vec<int64_t> &t, Weights &w, Normalizer &norma,
         ReducerBB<Int, Real> *red = 0, bool includeFirst = false);

   FigureOfMeritM(Weights &w, Normalizer &norma, ReducerBB<Int, Real> *red = 0);


   //===========================================================================

   virtual ~FigureOfMeritM();

   /**
    * Sets the vector @f$(t_1,..., t_d)@f$ in the FOM definition to the vector `t`.
    * The values of @f$t_1,..., t_d@f$ are taken from `t[0],...,t[d-1]`,
    * respectively.  When `includeFirst` is `true`, we consider only the non-successive
    * projections that contain coordinate 1.
    * See the doc of the class `FromRanges` in `CoordinateSets` for more details.
    */
   void setTVector(const NTL::Vec<int64_t> &t, bool includeFirst);

   /**
    * Counts and returns the total number of projections that are considered for the given `t`.
    */
   int64_t countProjections();

   /**
    * Sets the weights used for calculating the FoM
    */
   void setWeights(Weights &w) {
      m_weights = &w;
   }

   /**
    * Sets the normalizer to `norma`.
    */
   void setNormalizer(Normalizer &norma) {
      m_norma = &norma;
   }

   /**
    * Sets the parameters for the LLL reduction. If `delta = 0`, LLL is not applied.
    */
   void setLLL(double delta = 0.99999) {
      m_deltaLLL = delta;
   }

   /**
    * Sets the parameters for the BKZ reduction. If `delta = 0`, BKZ is not applied.
    */
   void setBKZ(double delta = 0.99999, int64_t blocksize = 10) {
      m_deltaBKZ = delta;
      m_blocksizeBKZ = blocksize;
   }

   /**
    * The BB method will be applied iff this flag is set to true (the default value).
    */
   void setBB(bool redBB) {
      m_redBB = redBB;
   }

   /**
    * Sets the `ReducerBB` object that will be used for BB.
    */
   void setReducerBB(ReducerBB<Int, Real> *red) {
      m_red = red;
   }

   /**
    * Sets the low bound for the FOM.
    * The FOM computation is stopped as soon as we know it is below this bound.
    * The default value is 0.
    */
   void setLowBound(double low) {
      m_lowbound = low;
   }

   /**
    * The level of verbosity in the terminal output.
    * The default value is 0 (minimal output).
    * Values from 1 to 5 give increasingly more details.
    * When the merit is smaller than `meritLevel` for one projection,
    * then the shortest vector is also printed for that projection.
    */
   void setVerbosity(int64_t verbose, double meritLevel = 0.0) {
      m_verbose = verbose;
      m_meritLevel = meritLevel;
   }

   /**
    * The level of details in the informations that are collected, such as the worst-case projection,
    * the corresponding shortest vector, its square length, etc.
    * The default value is 0 (minimal information, only the FOM).
    * With values from 1 to 4, we collect increasingly more details.
    */
   void setCollectLevel(int64_t collect) {
      m_collectLevel = collect;
   }

   /**
    * Returns the square length of the shortest vector for the worst-case projection.
    * The returned value is valid only if the collection level was at least 1.
    */
   double getMinMeritSqlen() {
      return m_minMeritSqlen;
   }

   /**
    * Returns the projection that gives the worst merit value.
    * The returned value is valid only if the collection level was at least 1.
    */
   Coordinates getMinMeritProj() {
      return m_minMeritProj;
   }

   /**
    * This function computes and returns the merit value for a single projection
    * represented in lattice `proj`, in `dim` dimensions.
    * It returns 0 if the computation is not completed for any reason.
    */
   double computeMeritOneProj(IntLattice<Int, Real> &proj, const Coordinates &coord,
         double minmerit = DBL_MAX);

   /**
    * This function computes and returns the FOM or the projections
    * over sets of successive coordinates of the form `{lowDim, lowDim+1, ..., j}`
    * for `j <= highDim`, and returns the minimum.
    * It uses `incDimBasis` to increase the dimension by 1 at each step.
    * It returns 0 if the computation was not completed for some reason.
    */
   virtual double computeMeritSucc(IntLatticeExt<Int, Real> &lat, int64_t lowDim, int64_t highDim,
         double minmerit = DBL_MAX);

   /**
    * This version takes `lowDim = m_t.length()` and `highDim = m_t[0]}.
    */
   virtual double computeMeritSucc(IntLatticeExt<Int, Real> &lat, double minmerit = DBL_MAX);

   /**
    * Same as `computeMeritSucc`, except that it rebuilds the basis anew each time
    * the dimension is increased, instead of using `incDimBasis`.
    * This approach is inefficient and this function is for experimentation only.
    */
   virtual double computeMeritSuccRebuild(IntLatticeExt<Int, Real> &lat, int64_t lowDim,
         int64_t highDim, double minmerit = DBL_MAX);
   virtual double computeMeritSuccRebuild(IntLatticeExt<Int, Real> &lat, double minmerit = DBL_MAX);

   /**
    * This function computes and returns the FOM only for the projections
    * over sets on non-successive coordinates determined by
    *  @f$S(s,t_s)@f$ or  @f$S^{(1)}(s,t_s)@f$.
    * It returns 0 if the computation was not completed for some reason.
    * The parameter `proj` is like for `computeMerit`.
    */
   virtual double computeMeritNonSucc(IntLatticeExt<Int, Real> &lat, IntLattice<Int, Real> &proj,
         double minmerit = DBL_MAX);

   /**
    * This function computes and returns the value of the FOM for the given lattice 'lat',
    * using the norm associated with that lattice.
    * The function returns 0 if the computation was not completed for some reason
    * (early exit, error, etc.).
    * The parameter `proj` points to a secondary `IntLattice` object used to store the
    * projections when computing the FOM.  The `maxDim` in this object must be large
    * enough so it can store any of the projections: `maxDim`\f$\ge \max(s,t_1)\f$.
    * Re-using this object permits one to avoid creating new objects internally.
    * We need a full `IntLattice` object (not only a basis) for when we apply BB to the projection.
    */
   double computeMerit(IntLatticeExt<Int, Real> &lat, IntLattice<Int, Real> &proj, double minmerit =
         DBL_MAX);

// The following variables should not be accessed directly.
protected:

   /*
    * This 'm_t' specifies the set of projections for which the FOM is computed.
    * `m_t[s-1]` represents t_s, for s=1,...,d.  `m_tsize is the value of d.
    */
   NTL::Vec<int64_t> m_t;
   int64_t m_tsize = 0;

   /*
    * Specifies the weights assigned to the projections.
    */
   Weights *m_weights;

   /*
    * Internal normalizer object for storing normalizing values in FiguresOfMeritM class
    */
   Normalizer *m_norma;

   /*
    * The parameters `delta` and `blocksize` for the LLL and BKZ reductions.
    * When `delta = 0`, the method is not applied.
    */
   double m_deltaLLL = 0.0;
   double m_deltaBKZ = 0.99999;
   int64_t m_blocksizeBKZ = 10;

   /*
    * The BB is applied iff this flag is set to true.
    */
   bool m_redBB = true;

   /*
    * Internal `CoordinateSets` object used to store the set of projections
    * that are considered for the non-successive coordinates.
    * It is constructed and populated by the `setTVector` function.
    * This object specifies a set of projections of each order.
    */
   CoordinateSets::FromRanges *m_coordRange;

   /*
    * Internal `reducerBB` object used when we apply BB to find the shortest vector in a projection.
    * We define it as a pointer because we want to pass a null pointer if BB is not applied.
    */
   ReducerBB<Int, Real> *m_red;

   /*
    * Pointer to a vector used to store the square Euclidean lengths of the basis vectors
    * after an LLL or BKZ pre-reduction via the static methods of `LLL_FP_lt`.
    */
   Real m_minSqlen;

   /*
    * Indicates if the first coordinate will always be included in all the projections
    * over the non-successive coordinates. This is used only when building  `*m_coordRange`,
    * but the variable could be useful for printouts of results.
    */
   bool m_includeFirst = false;

   /*
    * As soon as we know the FOM is above this bound, its calculation is stopped.
    * The default value is `infinite` to make sure that it is never exceeded.
    * Values larger than 1 might be possible because the normalizers
    * provide only approximations in some cases.
    */
   //  double m_highbound = DBL_MAX;
   /*
    * As soon as we know the FOM is below this bound, its calculation is stopped.
    */
   double m_lowbound = 0.0;

   /*
    * Variables to save the values at which the min merit is reached.
    */
   double m_minMerit = DBL_MAX; // The worst merit value.
   double m_minMeritSqlen = 0;  // Square length of worst-case vector.
   Coordinates m_minMeritProj;  // Stores the worst-case projection..

   /*
    * Variable to store the projection with the smallest FoM
    */
   //Coordinates m_worstproj;

   /*
    * Global timer for this class.
    */
   clock_t m_clock = 0;

   /*
    * Indicates how much details of FoM calculations are printed on the screen.
    */
   int64_t m_verbose = 0;

   /*
    * The function `computeMeritOneProj` prints the shortest vector when the merit
    * is smaller than this value and `m_verbose > 3`.
    */
   double m_meritLevel = 0.0;

   /*
    * Indicates how much detailed information is collected.
    */
   int64_t m_collectLevel = 0;

};

//============================================================================
// Implementation

template<typename Int, typename Real>
FigureOfMeritM<Int, Real>::FigureOfMeritM(const NTL::Vec<int64_t> &t, Weights &w, Normalizer &norma,
      ReducerBB<Int, Real> *red, bool includeFirst) {
   setTVector(t, includeFirst);
   m_weights = &w;
   setNormalizer(norma);
   m_red = red;
}

template<typename Int, typename Real>
FigureOfMeritM<Int, Real>::FigureOfMeritM(Weights &w, Normalizer &norma,
      ReducerBB<Int, Real> *red) {
   m_weights = &w;
   setNormalizer(norma);
   m_red = red;
}

//===========================================================================

template<typename Int, typename Real>
FigureOfMeritM<Int, Real>::~FigureOfMeritM() {
}

//=========================================================================

template<typename Int, typename Real>
void FigureOfMeritM<Int, Real>::setTVector(const NTL::Vec<int64_t> &t, bool includeFirst) {
   m_t = t;
   m_tsize = static_cast<int64_t>(t.length()); // Number of orders for the projections.
   m_includeFirst = includeFirst;
   // Reconstructs the CoordinateSets object.
   m_coordRange = new CoordinateSets::FromRanges;
   /* Defines the lower bound for the range of coordinates that belong to a projection.
    * It is 2 if the first coordinate belongs to all the projections, because we do
    * not have to consider coordinate 1. Otherwise it is 1.
    * In the first case, coordinate 1 is added to all the projections as an extra coord.
    * See the doc of the class `FromRanges` in `CoordinateSets`.
    */
   int64_t min_dim = 1;
   if (includeFirst) min_dim = 2;
   for (int64_t i = 1; i < m_tsize; i++) {
      // Adds the set of projections of order i, if non-empty.
      if (t[i] >= min_dim + i - includeFirst) m_coordRange->includeOrder(i + 1 - includeFirst,
            min_dim, t[i], includeFirst);
      else m_coordRange->excludeOrder(i + 1 - includeFirst);
   }
}

//=========================================================================

template<typename Int, typename Real>
int64_t FigureOfMeritM<Int, Real>::countProjections() {
   int64_t lower_dim = static_cast<int64_t>(this->m_t.length()) + 1;  // We start in d+1 dimensions.
   int64_t numProj = 1 + m_t[0] - lower_dim;
   for (auto it = m_coordRange->begin(); it != m_coordRange->end(); it++)
      numProj++;
   return numProj;
}

//=========================================================================
// Computes the merit value for one projection in `dim` dimensions.
// The dimension of `proj` must equal the size of coord.
template<typename Int, typename Real>
double FigureOfMeritM<Int, Real>::computeMeritOneProj(IntLattice<Int, Real> &proj,
      const Coordinates &coord, double minmerit) {
   int64_t dim = coord.size();  // Dimension of the projection.
   if (dim != proj.getDim())
      std::cout << "Error computeMeritOneProj:  size of coord = " << dim
            << ",  dim of projection = " << proj.getDim() << "\n";
   // assert (dim == proj.getDim());  // The dimensions must agree.  ******
   if (m_deltaLLL > 0.0) m_minSqlen = redLLL<Int, Real>(proj.getBasis(), m_deltaLLL, dim);
   if (m_deltaBKZ > 0.0)
      m_minSqlen = redBKZ<Int, Real>(proj.getBasis(), m_deltaBKZ, m_blocksizeBKZ, 0, dim);
   if (m_redBB) {       // We do the BB.
      // If `proj` is already the internal lattice for `m_red`, the following does nothing.
      // Otherwise it just sets a pointer to `proj`, and enlarges the arrays in m_red if needed.
      // m_red->setIntLattice(proj);
      if (!m_red->shortestVector(proj)) m_minSqlen = 0;
      else m_minSqlen = m_red->getMinLength2();
   }
   double merit;
   if (proj.getNormType() == L2NORM) NTL::conv(merit, sqrt(m_minSqlen) / m_norma->getBound(dim));
   else if (m_redBB) NTL::conv(merit, m_red->getMinLength() / m_norma->getBound(dim)); // For L1 norm.
   else {   // L1 norm and no BB
      proj.updateSingleVecNorm(0, dim);
      NTL::conv(merit, proj.getVecNorm(0) / m_norma->getBound(dim));
   }
   merit *= m_weights->getWeight(coord);
   if (merit < minmerit) {
      m_minMerit = merit;
      m_minMeritProj = coord;
      NTL::conv(m_minMeritSqlen, m_minSqlen);
   }
   if (m_verbose > 2) {
      if (dim < 8) std::cout << coord << std::setw(16 - 2 * dim) << " ";
      else std::cout << std::left << std::setw(2) << "{1,...," << dim << "}       ";
      if (7 < dim && dim < 10) std::cout << " ";
      std::cout << std::setw(12) << conv<double>(m_minSqlen) << "  "
            << std::setw(12) << conv<double>(1.0 / sqrt(m_minSqlen))
            << "  " << std::setw(10) << merit << "  " << m_minMerit << "    "
            << (double) (clock() - m_clock) / (CLOCKS_PER_SEC) << "\n";
   }
   if ((m_verbose > 4) || ((m_verbose > 3) && (merit < m_meritLevel))) {
      IntVec sv;
      sv.SetLength(dim);
      for (int64_t j = 0; j < dim; j++) sv[j] = proj.getBasis()[0][j];
      std::cout << "  short vec: " << sv  << "\n";
   }
   return merit;
}

//=========================================================================
template<typename Int, typename Real>
double FigureOfMeritM<Int, Real>::computeMeritSucc(IntLatticeExt<Int, Real> &lat, int64_t lowDim,
      int64_t highDim, double minmerit) {
   m_minMerit = minmerit;
   if (lowDim > highDim) return m_minMerit;  // No succ projection to look at.
   m_clock = clock();
   if (m_verbose > 2) {
      std::cout << "coordinates      sqlen         1/len      merit       minmerit    cumul sec \n";
   }
   Coordinates coord;
   for (int64_t j = 1; j <= lowDim; j++)
      coord.insert(j);
   lat.buildBasis(lowDim);
   // The dimension of `lat` here is equal to the size of coord.
   computeMeritOneProj(lat, coord, m_minMerit);   // This is the merit in lowDim dimensions.
   if (m_minMerit < this->m_lowbound) return 0;
   for (int64_t j = lowDim + 1; j <= highDim; j++) {
      coord.insert(j);
      lat.incDimBasis();
      computeMeritOneProj(lat, coord, m_minMerit);
      if (m_minMerit < this->m_lowbound) return 0;
   }
   return m_minMerit;
}

//=========================================================================
template<typename Int, typename Real>
double FigureOfMeritM<Int, Real>::computeMeritSucc(IntLatticeExt<Int, Real> &lat, double minmerit) {
   int64_t lowDim = static_cast<int64_t>(this->m_t.length()) + 1;  // We start in d+1 dimensions.
   return computeMeritSucc(lat, lowDim, this->m_t[0], minmerit);
}

//=========================================================================
template<typename Int, typename Real>
double FigureOfMeritM<Int, Real>::computeMeritSuccRebuild(IntLatticeExt<Int, Real> &lat,
      int64_t lowDim, int64_t highDim, double minmerit) {
   m_minMerit = minmerit;
   if (lowDim > highDim) return m_minMerit;  // No succ projection to look at.
   m_clock = clock();
   if (m_verbose > 2) {
      std::cout << "coordinates      sqlen         1/len      merit       minmerit    cumul sec \n";
   }
   Coordinates coord;
   for (int64_t j = 1; j <= lowDim; j++)
      coord.insert(j);
   lat.buildBasis(lowDim);
   // The dimension of `lat` here is equal to the size of coord.
   computeMeritOneProj(lat, coord, m_minMerit);
   if (m_minMerit < this->m_lowbound) return 0;
   for (int64_t j = lowDim + 1; j <= highDim; j++) {
      coord.insert(j);
      lat.buildBasis(j);
      computeMeritOneProj(lat, coord, m_minMerit);
      if (m_minMerit < this->m_lowbound) return 0;
   }
   return m_minMerit;
}

//=========================================================================
template<typename Int, typename Real>
double FigureOfMeritM<Int, Real>::computeMeritSuccRebuild(IntLatticeExt<Int, Real> &lat,
      double minmerit) {
   int64_t lowDim = static_cast<int64_t>(this->m_t.length()) + 1;  // We start in d+1 dimensions.
   return computeMeritSuccRebuild(lat, lowDim, this->m_t[0], minmerit);
}

//=========================================================================
template<typename Int, typename Real>
double FigureOfMeritM<Int, Real>::computeMeritNonSucc(IntLatticeExt<Int, Real> &lat,
      IntLattice<Int, Real> &proj, double minmerit) {
   // std::cout << "Start of computeMeritNonSucc in primal \n";
   m_minMerit = minmerit;
   m_clock = clock();
   if (m_verbose > 2) {
      std::cout << "coordinates      sqlen         1/len      merit       minmerit    cumul sec \n";
   }
   // assert (proj.getMaxDim() > m_tsize);
   Coordinates coord;
   for (auto it = m_coordRange->begin(); it != m_coordRange->end(); it++) {
      coord = *it;
      // We must make sure that we have all the required coordinates in the basis lat.
      while (*coord.end() > uint64_t(lat.getDim()))
         lat.incDimBasis();
      lat.buildProjection(proj, coord);
      // After this `buildProjection`, the dimensions of proj and coord must agree.
      computeMeritOneProj(proj, coord, m_minMerit);
      // std::cout << " Basis B after computeMerit = \n" << proj.getBasis() << "\n";
      if (m_minMerit <= this->m_lowbound) return 0;
   }
   return m_minMerit;
}

//=========================================================================

template<typename Int, typename Real>
double FigureOfMeritM<Int, Real>::computeMerit(IntLatticeExt<Int, Real> &lat,
      IntLattice<Int, Real> &proj, double minmerit) {
   // std::cout << "Start of computeMerit  \n";
   this->computeMeritNonSucc(lat, proj, minmerit);
   if (m_minMerit == 0) return 0;
   this->computeMeritSucc(lat, m_minMerit);
   // if (m_minMerit > this->m_highbound) return 0; // Removed. We want to see the large values!
   return m_minMerit;
}

template class FigureOfMeritM<std::int64_t, double> ;
template class FigureOfMeritM<NTL::ZZ, double> ;
template class FigureOfMeritM<NTL::ZZ, xdouble> ;
template class FigureOfMeritM<NTL::ZZ, quad_float> ;
template class FigureOfMeritM<NTL::ZZ, NTL::RR> ;

} // end namespace LatticeTester

#endif
