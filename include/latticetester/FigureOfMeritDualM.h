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

#ifndef LATTICETESTER_FIGUREOFMERITDUALM_H
#define LATTICETESTER_FIGUREOFMERITDUALM_H

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
#include "latticetester/Weights.h"
// #include "latticetester/WeightsOrderDependent.h"

#include "latticetester/NormaBestLat.h"
#include "latticetester/NormaBestBound.h"
#include "latticetester/WeightsUniform.h"
#include "latticetester/CoordinateSets.h"

namespace LatticeTester {

/**
 * This class offers funcitons to calculate the figure of merit (FOM)
 * for the m-dual of any given 'IntLatticeExt' object. It is implemented 
 * as a subclass of 'FigureOfMeritM.h' and thus requires the same input for the 
 * constructor. We refer the reader to 'FigureOfMeritM.h' for further details.
 * The main method of this class is 'computeMeritM' which calculates the actual 
 * FOM for the m-dual lattice of a given lattice. The computation is stopped 
 * (early exit) as soon as we know that the value of the FOM will be outside 
 * the interval `[low, high]`.
 * 
 * Note: The class works only for the case where  "PrecisionType == DOUBLE".   
 * This is a limitation.
 */
template<typename Int, typename Real>
class FigureOfMeritDualM: public FigureOfMeritM<Int, Real> {

private:
	typedef NTL::vector<Int> IntVec;
	typedef NTL::matrix<Int> IntMat;
	// typedef NTL::vector<Real> RealVec;
public:

	/*
	 * This constructor will call `setTVector` with the given vector `t`
	 * and `includeFirst` variable,
	 * then set the `Weights`, `ReductionType`, `Reducer`, and `Normalizer` to the given values.
	 */
	FigureOfMeritDualM(const NTL::vector<int64_t> &t, Weights &w,
			Normalizer &norma, ReducerBB<Int, Real> *red = 0,
			bool includeFirst = false);

	/*
	 * This function computes and returns the value of the FOM for the dual of the
	 * given lattice 'lat'. The function returns 0 if the computation was not completed
	 * for some reason (early exit, error, etc.).
	 * The parameter `proj` points to an `IntLattice` object that is used to store the
	 * projections when computing the FOM. The `maxDim` in this object must be large
	 * enough so it can store any of the projections: it must be at least \f$d$\f and
	 * at least \f$t_1$\f.
	 * Re-using this object permits one to avoid creating new objects internally.
	 */
	//double computeMeritDual (IntLatticeExt<Int, Real> & lat,
	//        IntLattice<Int, Real> &proj);
	/*
	 * Same as computeMeritSucc but for the m-dual lattice.
	 */
	double computeMeritSucc(IntLatticeExt<Int, Real> &lat) override;

	/*
	 * Same as computeMeritNonSucc but for the m-dual lattice.
	 */
	double computeMeritNonSucc(IntLatticeExt<Int, Real> &lat,
			IntLattice<Int, Real> &proj) override;
};

//============================================================================
// Implementation

template<typename Int, typename Real>
FigureOfMeritDualM<Int, Real>::FigureOfMeritDualM(const NTL::vector<int64_t> &t,
		Weights &w, Normalizer &norma, ReducerBB<Int, Real> *red,
		bool includeFirst) :
		FigureOfMeritM<Int, Real>(t, w, norma, red, includeFirst) {
}
;

/*
 //=========================================================================
 template<typename Int, typename Real>
 double FigureOfMeritDualM<Int, Real>::computeMeritDual (IntLatticeExt<Int, Real> &lat,
 IntLattice<Int, Real> &proj) {
 double minmerit = 1.0;
 // this->m_sqlen.SetLength(this->m_t.size() + 1);
 minmerit = min (minmerit, computeMeritNonSuccDual(lat, proj));
 // if (merit < minmerit) minmerit = merit;
 if (minmerit == 0) return 0;
 //this->m_sqlen.SetLength(this->m_t[0]);
 minmerit = min (minmerit, computeMeritSuccDual(lat));
 // if (merit < minmerit) minmerit = merit;
 // In any of these cases the calcuation is stopped
 if (minmerit == 0 || minmerit > this->m_highbound) return 0;
 return minmerit;
 }
 */

//=========================================================================
template<typename Int, typename Real>
double FigureOfMeritDualM<Int, Real>::computeMeritSucc(
		IntLatticeExt<Int, Real> &lat) {
	double minmerit = 1000.0;
	Coordinates coord;
	int64_t lower_dim = static_cast<int64_t>(this->m_t.size()); // We start in d dimensions.
	lat.buildDualBasis(lower_dim);
	for (int64_t j = 1; j < lower_dim + 1; j++)
		coord.insert(j);
	for (int64_t j = lower_dim + 1; j < this->m_t[0] + 1; j++) {
		coord.insert(j);
		lat.incDimDualBasis();
		lat.dualize();
		minmerit = min(minmerit, this->computeMeritOneProj(lat, coord));
		lat.dualize();
		if (minmerit <= this->m_lowbound)
			return 0;
	}
	return minmerit;
}

//=========================================================================
template<typename Int, typename Real>
double FigureOfMeritDualM<Int, Real>::computeMeritNonSucc(
		IntLatticeExt<Int, Real> &lat, IntLattice<Int, Real> &proj) {
	double minmerit = 1000.0;
	Coordinates coord;
	for (auto it = this->m_coordRange->begin(); it != this->m_coordRange->end();
			it++) {
		coord = *it;
		// The following builds a triangular basis for proj, takes its dual, and dualize.
		lat.buildProjectionDual(proj, coord);
		proj.dualize();
		minmerit = min(minmerit, this->computeMeritOneProj(proj, coord));
		proj.dualize();
		if (minmerit <= this->m_lowbound)
			return 0;
	}
	return minmerit;
}

//=========================================================================

template class FigureOfMeritDualM<std::int64_t, double> ;
template class FigureOfMeritDualM<NTL::ZZ, double> ;
template class FigureOfMeritDualM<NTL::ZZ, xdouble> ;
template class FigureOfMeritDualM<NTL::ZZ, quad_float> ;
template class FigureOfMeritDualM<NTL::ZZ, NTL::RR> ;

} // end namespace LatticeTester

#endif
