// This file is part of LatticeTester.
//
// Copyright (C) 2012-2022  The LatticeTester authors, under the supervision
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

#ifndef LATTICETESTER_NORMABESTUPBOUND_H
#define LATTICETESTER_NORMABESTUPBOUND_H

#include "latticetester/Normalizer.h"

#include <stdexcept>

namespace LatticeTester {

/**
 * \class NormaBestUpBound
 *
 * In this normalizer, the Hermite constants \f$\gamma_s\f$ are approximated using the
 * best upper bounds that are available.
 * For dimensions 0 through 8, and 24, the Hermite constant are known exactly.
 * For dimensions 9 through 36, except 24, they are computed using the article \cite mCOH03a.
 * For dimensions 37 through 48 we use Rogers's bounds.
 * The values given here are for the L2NORM.
 * For the L1NORM, the bounds must be multiplied by `s` in `s` dimensions.
 */
class NormaBestUpBound: public Normalizer {
public:

    /**
     * This constructor assumes that the rescaled primal lattice has scaling factor \f$m\f$
     * and order \f$k\f$, so its density is \f$m^{k-t}\f$ for \f$t\geq k\f$, and cannot
     * exceed 1 for projections in \f$s < k\f$ dimensions.
     */
    NormaBestUpBound(double logm, int k, int maxDim, NormType norm = L2NORM);

    /**
     * Constructs a `NormaBestUpBound` for up to `maxDim` dimensions, by assuming that the
     * log density is `logDensity` in all dimensions and the lattice was not rescaled.
     * Restriction: `maxDim`\f$ \le 48\f$.
     */
    NormaBestUpBound(double logDensity, int maxDim, NormType norm = L2NORM);

    /**
     * Destructor.
     */
    ~NormaBestUpBound();

    /**
     * Returns the value of the bound on the Hermite's constant \f$\gamma_j\f$
     * in dimension \f$j\f$.
     */
    double getGamma(int j) const;

private:

    /**
     * Precomputed lattice constants \f$\gamma_j\f$ for the best bound
     * in each dimension \f$j \le 48\f$.
     */
    static const double m_gamma[1 + Normalizer::MAX_DIM];

};
// End class NormaBestUpBound

//===========================================================================

/*
 * Values 1 through 36 are calculated from the article "A Conceptual
 * Breakthrough in sphere packing" from Henry Cohn (2016). They are computed
 * by taking gamma[n] = 4 * ( bound / V_n )^(n/2) where V_n is the volume of
 * an n dimensional sphere of radius 1.
 */
const double NormaBestUpBound::m_gamma[] = {
/* Gamma[0] = */0.0,
/* Gamma[1] = */1.0,
/* Gamma[2] = */1.1547005383793,   // exact
/* Gamma[3] = */1.2599210498949,   // exact
/* Gamma[4] = */1.4142135623731,   // exact
/* Gamma[5] = */1.5157165665104,   // exact
/* Gamma[6] = */1.6653663553112,   // exact
/* Gamma[7] = */1.8114473285278,   // exact
/* Gamma[8] = */2.0,               // exact
/* Gamma[9] = */2.1324942979385826,
/* Gamma[10] = */2.263452600069123,
/* Gamma[11] = */2.3930576417841327,
/* Gamma[12] = */2.5214594724429005,
/* Gamma[13] = */2.6487829527745377,
/* Gamma[14] = */2.7751332142995158,
/* Gamma[15] = */2.9005997190604442,
/* Gamma[16] = */3.025259312863684,
/* Gamma[17] = */3.149178557041412,
/* Gamma[18] = */3.2724155399984154,
/* Gamma[19] = */3.395021298677218,
/* Gamma[20] = */3.517040950462465,
/* Gamma[21] = */3.6385145917472883,
/* Gamma[22] = */3.7594780693416627,
/* Gamma[23] = */3.8799635383208293,
/* Gamma[24] = */4.0,          // exact
/* Gamma[25] = */4.119613703661868,
/* Gamma[26] = */4.238828491527929,
/* Gamma[27] = */4.357666099112316,
/* Gamma[28] = */4.476146404322391,
/* Gamma[29] = */4.594287637369152,
/* Gamma[30] = */4.7121065670603075,
/* Gamma[31] = */4.8296186570429365,
/* Gamma[32] = */4.94683820125262,
/* Gamma[33] = */5.063778445850659,
/* Gamma[34] = */5.180451691100012,
/* Gamma[35] = */5.296869382898868,
/* Gamma[36] = */5.413042186941963,
/* Gamma[37] = */5.7020718581143,
/* Gamma[38] = */5.8255656070255,
/* Gamma[39] = */5.9489276473284,
/* Gamma[40] = */6.0721635670068,
/* Gamma[41] = */6.1952785955803,
/* Gamma[42] = */6.3182776348,
/* Gamma[43] = */6.4411652860615,
/* Gamma[44] = */6.5639458749555,
/* Gamma[45] = */6.6866234733141,
/* Gamma[46] = */6.8092019190592,
/* Gamma[47] = */6.9316848341156,
/* Gamma[48] = */7.0540756406128 };

// If we use the upper bounds instead of exact values:
// Gamma[2] = 1.1547005395033834,
// Gamma[3] = 1.304077209522991,
// Gamma[4] = 1.4491512970689033,
// Gamma[5] = 1.5907345457787878,
// Gamma[6] = 1.7294425778558589,
// Gamma[7] = 1.8657438977611798,
// Gamma[8] = 2.0000000001950413,

/*=======================================================================*/

NormaBestUpBound::NormaBestUpBound(double logDensity, int maxDim, NormType norm) :
        Normalizer(maxDim, norm) {
   m_name = "NormaBestUpBound";
   Normalizer::computeBounds(logDensity);
}

/*=========================================================================*/

NormaBestUpBound::NormaBestUpBound(double logm, int k, int maxDim, NormType norm) :
        Normalizer(maxDim, norm) {
    if (maxDim > this->MAX_DIM)
        throw std::invalid_argument("NormaBestLat:   dimension > MAXDIM");
    m_name = "NormaBestUpBound";
    Normalizer::computeBounds(logm, k);
}

/*=========================================================================*/

NormaBestUpBound::~NormaBestUpBound() {
}

/*=========================================================================*/

inline double NormaBestUpBound::getGamma(int j) const {
    if (j < 1 || j > this->m_maxDim)
        throw std::out_of_range("NormaBestUpBound::getGamma");
    if (m_norm == L2NORM)
       return m_gamma[j];
    else if (m_norm == L1NORM)
       return m_gamma[j] * j;
    else
       throw std::domain_error("NormaBestUpBound::getGamma with wrong norm");
}

}

#endif
