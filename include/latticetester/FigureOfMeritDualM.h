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
#include "latticetester/BasisConstruction.h"
#include "latticetester/IntLattice.h"
#include "latticetester/IntLatticeExt.h"
#include "latticetester/Rank1Lattice.h"
#include "latticetester/ReducerStatic.h"
#include "latticetester/ReducerBB.h"
#include "latticetester/Weights.h"
#include "latticetester/NormaBestLat.h"
#include "latticetester/NormaBestUpBound.h"
#include "latticetester/WeightsUniform.h"
#include "latticetester/CoordinateSets.h"

namespace LatticeTester {

/**
 * \class FigureOfMeritDualM
 *
 * This class offers tools to calculate the same figure of merit (FOM) as `FigureOfMerit`,
 * but for the m-duals of the projections.
 *
 * That is, for each projection, a shortest vector
 * is computed for the m-dual of the projection, as explained in the guide,
 * and not for the projection of the $m$-dual lattice.
 * The only change from the parent class `FigureOfMerit` is in the
 * two functions `computeMeritSucc` and `computeMeritNonSucc`.
 * The rest is essentially the same.  The computation is stopped
 * (early exit) as soon as we know that the value of the FOM will be too small.
 */
template<typename Int, typename Real>
class FigureOfMeritDualM: public FigureOfMeritM<Int, Real> {

private:
	// typedef NTL::Vec<Int> IntVec;
	// typedef NTL::Mat<Int> IntMat;
	// typedef NTL::Vec<Real> RealVec;
public:

	/**
	 * This constructor will call `setTVector` with the given vector `t`
	 * and `includeFirst` variable, then set the `Weights`, `Normalizer` and `ReducerBB` to the given values.
	 */
	FigureOfMeritDualM(const NTL::Vec<int64_t> &t, Weights &w, Normalizer &norma,
	      ReducerBB<Int, Real> *red = 0, bool includeFirst = false);

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

	/**
	 * Same as `computeMeritSucc` in parent class, but for the m-dual lattice.
	 */
	double computeMeritSucc(IntLatticeExt<Int, Real> &lat, double minmerit = DBL_MAX) override;

	/**
	 * Same as `computeMeritNonSucc` in parent class, but for the m-dual lattice.
	 */
	double computeMeritNonSucc(IntLatticeExt<Int, Real> &lat,
			IntLattice<Int, Real> &proj, double minmerit = DBL_MAX) override;
};

//============================================================================
// Implementation

template<typename Int, typename Real>
FigureOfMeritDualM<Int, Real>::FigureOfMeritDualM(const NTL::Vec<int64_t> &t,
		Weights &w, Normalizer &norma, ReducerBB<Int, Real> *red, bool includeFirst) :
		FigureOfMeritM<Int, Real>(t, w, norma, red, includeFirst) {
};

//=========================================================================
template<typename Int, typename Real>
double FigureOfMeritDualM<Int, Real>::computeMeritSucc(
		IntLatticeExt<Int, Real> &lat, double minmerit) {
   this->m_minMerit = minmerit;
	Coordinates coord;
   if (this->m_verbose > 1) std::cout << "coordinates      merit        sqlen   minmerit  \n";
	int64_t lower_dim = static_cast<int64_t>(this->m_t.length()) + 1; // We start in d+1 dimensions.
   for (int64_t j = 1; j <= lower_dim; j++)
      coord.insert(j);
   lat.buildDualBasis(lower_dim);
   lat.dualize();
   this->computeMeritOneProj(lat, coord, this->m_minMerit);
   lat.dualize();
   if (this->m_minMerit < this->m_lowbound) return 0;
   for (int64_t j = lower_dim + 1; j < this->m_t[0] + 1; j++) {
		coord.insert(j);
		// std::cout << "inDimDualBasis with j = " << j << ", dim = " << lat.getDim() << ", maxdim = " << lat.getMaxDim() << "\n";
		lat.incDimDualBasis();
		lat.dualize();
      this->computeMeritOneProj(lat, coord, this->m_minMerit);
		lat.dualize();
		if (this->m_minMerit <= this->m_lowbound)  return 0;
	}
	return this->m_minMerit;
}

//=========================================================================
template<typename Int, typename Real>
double FigureOfMeritDualM<Int, Real>::computeMeritNonSucc(
		IntLatticeExt<Int, Real> &lat, IntLattice<Int, Real> &proj, double minmerit) {
   this->m_minMerit = minmerit;
	Coordinates coord;
   if (this->m_verbose > 1) std::cout << "coordinates      merit        sqlen   minmerit  \n";
	for (auto it = this->m_coordRange->begin(); it != this->m_coordRange->end();
			it++) {
		coord = *it;
		// The following builds a triangular basis for proj, takes its dual, and dualize.
      lat.buildProjectionDual(proj, coord);
      assert(proj.getDim() == proj.getDimDual());
      assert(proj.getDim() == (unsigned) coord.size());
      proj.dualize();
      this->computeMeritOneProj(proj, coord, this->m_minMerit);
      proj.dualize();
		if (this->m_minMerit <= this->m_lowbound)
			return 0;
	}
	return this->m_minMerit;
}

//=========================================================================

template class FigureOfMeritDualM<std::int64_t, double> ;
template class FigureOfMeritDualM<NTL::ZZ, double> ;
template class FigureOfMeritDualM<NTL::ZZ, xdouble> ;
template class FigureOfMeritDualM<NTL::ZZ, quad_float> ;
template class FigureOfMeritDualM<NTL::ZZ, NTL::RR> ;

} // end namespace LatticeTester

#endif
