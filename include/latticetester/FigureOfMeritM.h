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
#include "latticetester/BasisConstruction.h"
#include "latticetester/Util.h"
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
 * This class provides a tool to calculate the figure of merit (FOM)
 * \f{equation}{
 *    M_{t_1,\dots,t_d} = \min\left[ \min_{I\in S(t_1)} \omega_I \ell_I/\ell_I^*(\eta_I),\;
 *    \min_{2\le s\le d}\, \min_{I\in S(s,t_s)} \omega_I \ell_I/\ell_I^*(\eta_I) \right],
 * \f}
 * defined in Section 10 of the guide, for any given `IntLatticeExt` object.
 * The FOM is computed only for the (rescaled) primal lattice, the m-dual is never used.
 * To compute the FOM for the dual, one should first dualize the lattice.
 * The projections in \f$S(t_1)$\f are those over successive coordinates in up to \f$t_1$\f
 * dimensions, while those is @f$S(s,t_s)@f$ are projections over \f$s$\f distinct coordinates
 * that are non necessarily successive and are all in the set \f$\{1,\dots,t_s\}$\f,
 * for each order @f$s > 1@f$.
 * There are two variants for the latter: the first (default) variant takes
 * @f$S(s,t_s)@f$ as just defined (also defined in the guide), and the other considers
 * only the set @f$S^{(1)}(s,t_s)@f$ of projections that contain coordinate 1.
 * The parameter `includeFirst` in the constructor determines which variant is taken:
 * the latter option is taken when `includeFirst` is set to `true`.
 *
 * The lengths of the shortest vectors in the projections can be calculated exactly by using the
 * BB algorithm after applying some pre-reduction, or they can be just approximated
 * by the lengths of the shortest basis vector obtained after applying some pre-reduction
 * such as LLL or BKZ.  The latter is much faster but not exact.
 *
 * The constructor has two template parameters to specify which `Int` and `Real` types are used.
 * It also requires the vector @f$(t_1,\dots,t_d)@f$,
 * a `Weights` object that gives a weight to each projection,
 * a `Normalizer` object used to normalize the merit values,
 * a `ReducerBB` object used for the reduction in case we want to apply BB, and the
 * `includeFirst` parameter in case we want to change it to `false`.
 * The last two parameters are optional.
 * The BB is applied if and only if a (nonzero) `ReducerBB` is given.
 * Otherwise, we just use static methods for the reduction and need no `Reducer`.
 * By default, the pre-reduction method is only BKZ with `delta = 0.99999` and `blocksize = 10`.
 * To change these values and/or apply LLL, one can use `setBKZ` and/or `setLLL`.
 * The reductions are always applied in the order: LLL, BKZ, BB.
 * To remove LLL or BKZ, it suffices to set its `delta` parameter to 0.0.
 *
 * After a `FigureOfMeritM` has been created by the constructor, one can use `setTVector`
 * to set (or reset) the vector  @f$(t_1,\dots,t_d)@f$, `setWeights` to change the weights,
 * `setLLL` or `setBKZ` to change the LLL or BKZ parameters, etc.
 *
 * The method `computeMerit` computes the FOM for a given lattice.
 * The computation is stopped (early exit) as soon as we know that the value of the FOM
 * will be outside the interval `[low, high]`.  By default, these two bounds are 0 and 2
 * so there is no early stopping, but they can be changed via the function `setBounds`.
 *
 * *****  WARNING: For now, this class works only for the Euclidean norm, right?    *****
 */

template<typename Int, typename Real>
class FigureOfMeritM {

private:
   typedef NTL::vector<Int> IntVec;
   typedef NTL::matrix<Int> IntMat;
   // typedef NTL::vector<Real> RealVec;
public:

   /*
    * This constructor will call `setTVector` with the given vector `t`
    * and `includeFirst` variable,
    * then set the 'Weights', `Normalizer`, and `ReducerBB` to the given values.
    */
   FigureOfMeritM(const NTL::vector<int64_t> &t, Weights &w, Normalizer &norma,
         ReducerBB<Int, Real> *red = 0, bool includeFirst = true);

   //===========================================================================

   virtual ~FigureOfMeritM();

   /*
    * Sets the vector @f$(t_1,..., t_d)@f$ in the FOM definition to the vector `t`.
    * Note that the values of @f$t_1,..., t_d@f$ are taken from `t[0],...,t[d-1]`,
    * respectively.  When `includeFirst` is `true`, we consider only the non-successive
    * projections that contain coordinate 1.
    * See the doc of the class `FromRanges` in `CoordinateSets` for more details.
    */
   void setTVector(const NTL::vector<int64_t> &t, bool includeFirst = true);

   /*
    * Sets the weights used for calculating the FoM
    */
   void setWeights(Weights &w) {
      m_weights = &w;
   }

   /*
    * Sets the normalizer to `norma`.
    */
   void setNormalizer(Normalizer &norma) {
      m_norma = &norma;
   }

   /*
    * Sets the parameters for the LLL reduction. If `delta = 0`, LLL is not applied.
    */
   void setLLL(double delta = 0.99999) {
      m_deltaLLL = delta;
   }

   /*
    * Sets the parameters for the BKZ reduction. If `delta = 0`, BKZ is not applied.
    */
   void setBKZ(double delta = 0.99999, int64_t blocksize = 10) {
      m_deltaBKZ = delta;
      m_blocksizeBKZ = blocksize;
   }

   /*
    * The BB method will be applied iff this flag is set to true (the default value).
    */
   void setBB(bool redBB) {
      m_redBB = redBB;
   }

   /*
    * Sets the `ReducerBB` object that will be used for BB.
    */
   void setReducerBB(ReducerBB<Int, Real> *red) {
      *m_red = *red;
   }

   /*
    * Sets the low and the high bound for the FOM.
    * The FOM computation is stopped as soon as we know it is outside these bounds.
    * The default values are 0 and 2.
    */
   void setBounds(double &low, double &high) {
      m_highbound = high;
      m_lowbound = low;
   }

   /*
    * The level of verbosity in the terminal output.
    * The default value is 0 (minimal output).
    * Values from 1 to 4 give increasingly more details.
    */
   void setVerbosity(int64_t verbose) {
      m_verbose = verbose;
   }

   /*
    * The level of details in the informations that are collected, such as the worst-case projection,
    * the corresponding shortest vector, its square length, etc.
    * The default value is 0 (minimal information, only the FOM).
    * With values from 1 to 4, we collect increasingly more details.
    */
   void setCollectLevel(int64_t collect) {
      m_collectLevel = collect;
   }

   /*
    * This function computes and returns the value of the FOM for the given lattice 'lat'.
    * The function returns 0 if the computation was not completed for some reason
    * (early exit, error, etc.).
    * The parameter `proj` points to a secondary `IntLattice` object used to store the
    * projections when computing the FOM.  The `maxDim` in this object must be large
    * enough so it can store any of the projections: `maxDim`\f$\ge \max(s,t_1)$\f.
    * Re-using this object permits one to avoid creating new objects internally.
    * We need a full `IntLattice` object (not only a basis) for when we apply BB to the projection.
    */
   double computeMerit(IntLatticeExt<Int, Real> &lat, IntLattice<Int, Real> &proj, double minmerit = 1000.0);

   /*
    * This function computes and returns the FOM only for the projections
    * over sets of successive coordinates of the form
    * {1, 2, ..., m_t.size()} to {1, 2, ..., m_t[0]}, and returns the minimum.
    * It returns 0 if the computation was not completed for some reason.
    */
   virtual double computeMeritSucc(IntLatticeExt<Int, Real> &lat, double minmerit = 1000.0);

   /*
    * This function computes and returns the FOM only for the projections
    * over sets on non-successive coordinates determined by
    *  @f$S(s,t_s)@f$ or  @f$S^{(1)}(s,t_s)@f$.
    * It returns 0 if the computation was not completed for some reason.
    * The parameter `proj` is like for `computeMerit`.
    */
   virtual double computeMeritNonSucc(IntLatticeExt<Int, Real> &lat, IntLattice<Int, Real> &proj, double minmerit = 1000.0);

   /*
    * This function computes and returns the merit value for a single projection
    * represented in lattice `proj`, in `dim` dimensions.
    * It returns 0 if the computation is not completed for any reason.
    */
   double computeMeritOneProj(IntLattice<Int, Real> &proj, const Coordinates &coord, double minmerit = 1000.0);

   /*
    * This 'm_t' specifies the set of projections for which the FOM is computed.
    * `m_t[s-1]` represents t_s, for s=1,...,d.  `m_tsize` is the value of d.
    */
   NTL::vector<int64_t> m_t;
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
   NTL::vector<Real> m_sqlen;

   /*
    * Indicates if the first coordinate will always be included in all the projections
    * over the non-successive coordinates. This is used only when building  `*m_coordRange`,
    * but the variable could be useful for printouts of results.
    */
   bool m_includeFirst = true;

   /*
    * As soon as we know the FOM is above this bound, its calculation is stopped.
    * The default value is 1000 to make sure that it is never exceeded.
    * Values slightly larger than 1 might be possible because the normalizers
    * provide only approximations in some cases.
    */
   double m_highbound = 1000.0;

   /*
    * As soon as we know the FOM is below this bound, its calculation is stopped.
    */
   double m_lowbound = 0.0;


   /*
    * Variables to save the values at which the min merit is reached.
    */
   double m_minMerit = 1000.0;
   double m_minMeritSqlen = 0;  // Square length of worst-case vector.
   Coordinates m_minMeritProj;  // Initialized to empty set of coordinates.

   /*
    * Indicates how much details of FoM calculations are printed on the screen.
    */
   int64_t m_verbose = 0;

   /*
    * Indicates how much detailed information is collected.
    */
   int64_t m_collectLevel = 0;

};

//============================================================================
// Implementation

template<typename Int, typename Real>
FigureOfMeritM<Int, Real>::FigureOfMeritM(const NTL::vector<int64_t> &t, Weights &w,
      Normalizer &norma, ReducerBB<Int, Real> *red, bool includeFirst) {
   setTVector(t, includeFirst);
   m_weights = &w;
   setNormalizer(norma);
   m_red = red;
   m_sqlen.SetLength(1); // We will retrieve only the square length of the shortest.
}

//===========================================================================

template<typename Int, typename Real>
FigureOfMeritM<Int, Real>::~FigureOfMeritM() {
}

//=========================================================================

template<typename Int, typename Real>
void FigureOfMeritM<Int, Real>::setTVector(const NTL::vector<int64_t> &t, bool includeFirst) {
   m_t = t;
   m_tsize = static_cast<int64_t>(t.size()); // Number of orders for the projections.
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
// Computes the merit value for one projection in dim dimensions.
template<typename Int, typename Real>
double FigureOfMeritM<Int, Real>::computeMeritOneProj(IntLattice<Int, Real> &proj,
      const Coordinates &coord, double minmerit) {
   int64_t dim = coord.size();  // Dimension of the projection.
   if (m_deltaLLL > 0.0)
      redLLL<IntMat, NTL::vector<Real>>(proj.getBasis(), m_deltaLLL, dim, &m_sqlen);
   if (m_deltaBKZ > 0.0)
      redBKZ<IntMat, NTL::vector<Real>>(proj.getBasis(), m_deltaBKZ, m_blocksizeBKZ, 0, dim,
            &m_sqlen);
   if (m_red) {       // We do the BB.
      // If `proj` is already the internal lattice for `m_red`, the following does nothing.
      // Otherwise it just sets a pointer to `proj`, and enlarges the arrays in m_red if needed.
      // m_red->setIntLattice(proj);
      if (!m_red->shortestVector(proj)) return 0;
      m_sqlen[0] = m_red->getMinLength2();
   }
   double merit = 0.0;
   if (proj.getNormType() == L2NORM)
      NTL::conv(merit, sqrt(m_sqlen[0]) / m_norma->getBound(dim));
   else
      NTL::conv(merit, m_red->getMinLength() / m_norma->getBound(dim)); // For L1 norm.
   merit *= m_weights->getWeight(coord);
   if (merit < minmerit) {
      m_minMerit = merit;
      m_minMeritProj = coord;
      NTL::conv(m_minMeritSqlen, m_sqlen[0]);
      if (m_collectLevel > 0) {};   // Maybe store the shortest vector.
   }
   if (m_verbose > 1) {
      if (dim < 8) std::cout << coord << std::setw(15 - 2*dim) << " ";
      else std::cout << std::left << std::setw(2) << "{1,...," << dim << "}       ";
      std::cout << merit << "  " << m_sqlen[0] << "  " << minmerit << "\n";
   }
   return merit;
}

//=========================================================================

template<typename Int, typename Real>
double FigureOfMeritM<Int, Real>::computeMerit(IntLatticeExt<Int, Real> &lat,
      IntLattice<Int, Real> &proj, double minmerit) {
   // m_minMerit = minmerit;
   this->computeMeritNonSucc(lat, proj, minmerit);
   if (m_minMerit == 0) return 0;
   this->computeMeritSucc(lat, m_minMerit);
   if (m_minMerit > this->m_highbound) return 0;   // Maybe this could be removed.
   return m_minMerit;
}

//=========================================================================
template<typename Int, typename Real>
double FigureOfMeritM<Int, Real>::computeMeritSucc(IntLatticeExt<Int, Real> &lat, double minmerit) {
   m_minMerit = minmerit;
   if (m_verbose > 1) std::cout << "\n coordinates    merit           sqlen   minmerit \n";
   Coordinates coord;
   int64_t lower_dim = static_cast<int64_t>(this->m_t.size());  // Could be made more flexible.
   lat.buildBasis(lower_dim);
   for (int64_t j = 1; j < lower_dim + 1; j++)
      coord.insert(j);
   for (int64_t j = lower_dim + 1; j < this->m_t[0] + 1; j++) {
      coord.insert(j);
      lat.incDimBasis();
      // minmerit = min(minmerit, computeMeritOneProj(lat, coord, minmerit));
      computeMeritOneProj(lat, coord, m_minMerit);
      if (m_minMerit <= this->m_lowbound) return 0;
   }
   return m_minMerit;
}

//=========================================================================
template<typename Int, typename Real>
double FigureOfMeritM<Int, Real>::computeMeritNonSucc(IntLatticeExt<Int, Real> &lat,
      IntLattice<Int, Real> &proj, double minmerit) {
   m_minMerit = minmerit;
   if (m_verbose > 1) std::cout << "\n coordinates    merit           sqlen   minmerit  \n";
   // std::cout << "proj.getDim = " << proj.getDim() << ",  m_tsize = " << m_tsize << "\n";
   // assert (proj.getMaxDim() > m_tsize);
   Coordinates coord;
   for (auto it = m_coordRange->begin(); it != m_coordRange->end(); it++) {
      coord = *it;
      // We must make sure that we have all the required coordinates in the basis.
      while (*coord.end() > uint64_t(lat.getDim())) lat.incDimBasis();
      lat.buildProjection(proj, coord);
      // std::cout << " Basis B = \n" << proj.getBasis() << "\n";
      // minmerit = min(minmerit, computeMeritOneProj(proj, coord));
      computeMeritOneProj(proj, coord, m_minMerit);
      // std::cout << " Basis B after computeMerit = \n" << proj.getBasis() << "\n";
      if (m_minMerit <= this->m_lowbound) return 0;
   }
   return m_minMerit;
}

template class FigureOfMeritM<std::int64_t, double> ;
template class FigureOfMeritM<NTL::ZZ, double> ;
template class FigureOfMeritM<NTL::ZZ, xdouble> ;
template class FigureOfMeritM<NTL::ZZ, quad_float> ;
template class FigureOfMeritM<NTL::ZZ, NTL::RR> ;

} // end namespace LatticeTester

#endif
