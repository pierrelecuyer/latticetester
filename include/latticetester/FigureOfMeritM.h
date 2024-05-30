// This file is part of LatticeTester.
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
 * `includeFirst` parameter in case we want to put it to `true`.
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
 * will be outside the interval `[low, high]`.  By default,
 * these two bounds are 0 and 1, but they can be changed via the function `setBounds`.
 *
 * *****  WARNING: For now, this class works only for the Euclidean norm!    *****
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
	// template<typename Int, typename Real>
	FigureOfMeritM(const NTL::vector<int64_t> &t, Weights &w, Normalizer &norma,
			ReducerBB<Int, Real> &red = 0, bool includeFirst = false);

	/*
	 * Sets the vector @f$(t_1,..., t_d)@f$ in the FOM definition to the vector `t`.
	 * Note that the values of @f$t_1,..., t_d@f$ are taken from `t[0],...,t[d-1]`,
	 * respectively.  When `includeFirst` is `true`, we consider only the non-successive
	 * projections that contain coordinate 1.
	 * See the doc of the class `FromRanges` in `CoordinateSets` for more details.
	 */
	void setTVector(const NTL::vector<int64_t> &t, bool includeFirst = false);

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
	 * Sets the low and the high bound for the FOM.
	 * The FOM computation is stopped as soon as we know it is outside these bounds.
	 */
	void setBounds(double &low, double &high) {
		m_highbound = high;
		m_lowbound = low;
	}

	/*
	 * Sets if details of FoM shall be printed on the screen
	 */
	void setPrintDetails(bool print) {
		m_printDetails = print;
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
	 *
	 * *** In this function and the two that follow, one lattice is passed by reference (&lat)
	 *     and the other via a pointer (*proj).  Why is it different?  Is there a reason?
	 *     It makes the code more complicated.                                            *****
	 */
	double computeMerit(IntLatticeExt<Int, Real> &lat,
			IntLattice<Int, Real> &proj);

	/*
	 * This function computes and returns the FOM only for the projections
	 * over sets of successive coordinates of the form
	 * {1, 2, ..., m_t.size()} to {1, 2, ..., m_t[0]}.
	 * It returns 0 if the computation was not completed for some reason.
	 * The parameter `proj` is like for `computeMerit`.
	 */
	double computeMeritSucc(IntLatticeExt<Int, Real> &lat);

	/*
	 * This function computes and returns the FOM only for the projections
	 * over sets on non-successive coordinates determined by
	 *  @f$S(s,t_s)@f$ or  @f$S^{(1)}(s,t_s)@f$.
	 * It returns 0 if the computation was not completed for some reason.
	 * The parameter `proj` is like for `computeMerit`.
	 */
	double computeMeritNonSucc(IntLatticeExt<Int, Real> &lat,
			IntLattice<Int, Real> &proj);

	/*
	 * This function computes and returns the merit value for a single projection
	 * represented in lattice `proj`, in `dim` dimensions.
	 * It returns 0 if the computation is not completed for any reason.
	 */
	double computeMeritOneProj(IntLattice<Int, Real> &proj, const Coordinates &coord);

	/*
	 * This 'm_t' specifies the set of projections for which the FOM is computed.
	 * `m_t[s-1]` represents t_s, for s=1,...,d.
	 */
	NTL::vector<int64_t> m_t;

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
	 * The parameters `delta` used to build projections via LLL.
	 */
	double m_deltaProj = 0.5;

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
	bool m_includeFirst = false;

	/*
	 * As soon as we know the FOM is above this bound, its calculation is stopped.
	 */
	double m_highbound = 1.0;

	/*
	 * As soon as we know the FOM is below this bound, its calculation is stopped.
	 */
	double m_lowbound = 0.0;

	/*
	 * Indicates if details of FoM calculation are printed on the screen.
	 */
	bool m_printDetails = false;

};

//============================================================================
// Implementation

template<typename Int, typename Real>
FigureOfMeritM<Int, Real>::FigureOfMeritM(const NTL::vector<int64_t> &t,
		Weights &w, Normalizer &norma, ReducerBB<Int, Real> &red,
		bool includeFirst) {
	setTVector(t, includeFirst);
	m_weights = &w;
	setNormalizer(norma);
	m_red = &red;
	m_sqlen.SetLength(1);  // We will retrieve only the square length of the shortest.
}

//=========================================================================

template<typename Int, typename Real>
void FigureOfMeritM<Int, Real>::setTVector(const NTL::vector<int64_t> &t,
		bool includeFirst) {
	m_t = t;
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
	if (includeFirst)
		min_dim = 2;
	int64_t d = static_cast<int64_t>(t.size()); // Number of orders for the projections.
	for (int64_t i = 1; i < d; i++) {
		// Adds the set of projections of order i, if non-empty.
		if (t[i] >= min_dim + i - includeFirst) {
			m_coordRange->includeOrder(i + 1 - includeFirst, min_dim, t[i],
					includeFirst);
		}
	}
}

//=========================================================================

template<typename Int, typename Real>
double FigureOfMeritM<Int, Real>::computeMerit(IntLatticeExt<Int, Real> &lat,
		IntLattice<Int, Real> &proj) {
	double merit = 0;
	double minmerit = 1.0;
	// this->m_sqlen.SetLength(1);

	merit = computeMeritNonSucc(lat, proj);
	if (merit < minmerit)
		minmerit = merit;
	// If the calculation was stopped, we return 0.
	if (minmerit == 0)
		return 0;

	merit = computeMeritSucc(lat);
	if (merit < minmerit)
		minmerit = merit;
	// In any of these cases the calculation is stopped.
	if (minmerit == 0 || merit > m_highbound)
		return 0;
	return minmerit;
}

//=========================================================================
// Computes the merit value for one projection in dim dimensions.
template<typename Int, typename Real>
double FigureOfMeritM<Int, Real>::computeMeritOneProj(IntLattice<Int, Real> &proj, const Coordinates &coord) {
	int64_t dim = proj.getDim();
	if (m_deltaLLL > 0.0)
		redLLL<IntMat, NTL::vector<Real>>(proj.getBasis(), m_deltaLLL, dim,
				&m_sqlen);
	if (m_deltaBKZ > 0.0)
		redBKZ<IntMat, NTL::vector<Real>>(proj.getBasis(), m_deltaBKZ,
				m_blocksizeBKZ, 0, dim, &m_sqlen);
	if (this->m_red) {  // If m_red is not a null pointer, we do the BB.
		this->m_red->setIntLattice(proj); // Do we need to redo this each time?   Is this safe?   *****
		if (!this->m_red->shortestVector())
			return 0;
		this->m_sqlen[0] = this->m_red->getMinLength2();
	}
	double merit = 0.0;
	if (proj.getNormType() == L2NORM)
		NTL::conv(merit, sqrt(this->m_sqlen[0]) / this->m_norma->getBound(dim));
	else
		// This is in case we use the L1 norm.
		NTL::conv(merit, this->m_red->getMinLength() / this->m_norma->getBound(dim));
	merit *= m_weights->getWeight(coord);
	if (m_printDetails)
		std::cout << "Coordinates: {1,...," << dim << "}, FoM: " << merit << "\n";
	return merit;
}

//=========================================================================
template<typename Int, typename Real>
double FigureOfMeritM<Int, Real>::computeMeritSucc(
		IntLatticeExt<Int, Real> &lat) {
	double minmerit = 1.0;
	double merit;
	Coordinates coord;
	int64_t lower_dim = static_cast<int64_t>(this->m_t.size());
	lat.buildBasis(lower_dim);
	for (int64_t j = 1; j < lower_dim + 1; j++)
		coord.insert(j);
	for (int64_t j = lower_dim + 1; j < this->m_t[0] + 1; j++) {
		coord.insert(j);
		lat.incDimBasis();
		merit = computeMeritOneProj (lat, coord);
		if (merit < minmerit)
			minmerit = merit;
		if (minmerit <= this->m_lowbound)
			return 0;
	}
	return minmerit;
}

/*
//=========================================================================
template<typename Int, typename Real>
// double FigureOfMeritM<Int>::computeMeritSucc(IntLatticeExt<Int, Real> &lat,
double FigureOfMeritM<Int, Real>::computeMeritSuccOld(
		IntLatticeExt<Int, Real> &lat) {
	double merit = 0;
	double minmerit = 1.0;
	Coordinates coord;
	int64_t lower_dim = static_cast<int64_t>(this->m_t.size());
	lat.buildBasis(lower_dim);

	for (int64_t j = 1; j < lower_dim + 1; j++)
		coord.insert(j);
	for (int64_t j = lower_dim + 1; j < this->m_t[0] + 1; j++) {
		coord.insert(j);
		lat.incDimBasis();
		if (this->m_deltaLLL > 0.0)
			redLLL<IntMat, RealVec>(lat.getBasis(), this->m_deltaLLL, j,
					&this->m_sqlen);
		if (this->m_deltaBKZ > 0.0)
			redBKZ<IntMat, RealVec>(lat.getBasis(), this->m_deltaBKZ,
					this->m_blocksizeBKZ, 0, j, &this->m_sqlen);
		if (this->m_red) {
			// We do the BB.
			this->m_red->setIntLattice(lat);
			if (!this->m_red->shortestVector())
				return 0;
			NTL::conv(merit,
					this->m_red->getMinLength() / this->m_norma->getBound(j));
		} else {
			if (lat.getNormType() == L2NORM)
				NTL::conv(merit,
						sqrt(this->m_sqlen[0]) / this->m_norma->getBound(j));
			else
				NTL::conv(merit, this->m_sqlen[0] / this->m_norma->getBound(j));
		}
		merit = merit * this->m_weights->getWeight(coord);
		if (m_printDetails)
			std::cout << "Coordinates: {1,...," << j << "}, FoM: " << merit
					<< "\n";
		if (merit < minmerit)
			minmerit = merit;
		if (minmerit <= this->m_lowbound)
			return 0;
	}
	return minmerit;
}
*/

//=========================================================================
template<typename Int, typename Real>
double FigureOfMeritM<Int, Real>::computeMeritNonSucc(
		IntLatticeExt<Int, Real> &lat, IntLattice<Int, Real> &proj) {
	double merit = 0;
	double minmerit = 1.0;
	Coordinates coord;
	for (auto it = m_coordRange->begin(); it != m_coordRange->end(); it++) {
		coord = *it;
		lat.buildProjection(proj, coord, this->m_deltaProj); // Done with LLL. Must have withDual = false   ***
		merit = computeMeritOneProj (proj, coord);
		if (merit < minmerit)
			minmerit = merit;
		if (minmerit <= this->m_lowbound)
			return 0;
	}
	return minmerit;
}

template class FigureOfMeritM<std::int64_t, double> ;
template class FigureOfMeritM<NTL::ZZ, double> ;
template class FigureOfMeritM<NTL::ZZ, xdouble> ;
template class FigureOfMeritM<NTL::ZZ, quad_float> ;
template class FigureOfMeritM<NTL::ZZ, NTL::RR> ;

} // end namespace LatticeTester

#endif
