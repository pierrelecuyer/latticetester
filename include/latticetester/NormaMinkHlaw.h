// This file is part of LatticeTester.
//
// Copyright (C) 2012-2022  The LatticeTester authors, under the occasional supervision
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

#ifndef LATTICETESTER_NORMAMINKHLAW_H
#define LATTICETESTER_NORMAMINKHLAW_H

#include "latticetester/Normalizer.h"

#include <stdexcept>

namespace LatticeTester {

/**
 * \class NormaMinkHlaw
 *
 * This class implements *lower bounds* on the Hermite constants based on the
 * Minkowski-Hlawka theorem \cite mHLA43a. The Hermite constants are approximated
 * by the \f$\gamma_t^{(Z)}\f$ given in the guide.
 * These lower bounds are all smaller than the lower bounds from `NormaBestLat`.
 */
class NormaMinkHlaw: public Normalizer {
public:

    /**
     * This constructor assumes that the rescaled primal lattice has scaling factor \f$m\f$
     * and order \f$k\f$, so its density is \f$m^{k-t}\f$ for \f$t\geq k\f$, and cannot
     * exceed 1 for projections in \f$s < k\f$ dimensions.
     */
    NormaMinkHlaw(double logm, int64_t k, int64_t maxDim, NormType norm = L2NORM);

    /**
     * Constructs a `NormaMinkHlaw` for up to `maxDim` dimensions, by assuming that the
     * log density is `logDensity` in all dimensions and the lattice was not rescaled.
     * Restriction: `maxDim`\f$ \le 48\f$.
     */
    NormaMinkHlaw(double logDensity, int64_t maxDim, NormType norm = L2NORM);


    /**
     * Creates the Normalizer object without computing any bounds.
     */
    NormaMinkHlaw(int64_t maxDim, NormType norm = L2NORM);

    /**
     * Returns the value of the lattice constant \f$\gamma_j\f$ in
     * dimension \f$j\f$.
     */
    double getGamma(int64_t j) const;

private:

    /**
     * Constants \f$\gamma_j\f$ for the Minkowski-Hlawka bounds in each
     * dimension \f$j\f$.
     */
    static const double m_gamma[1 + Normalizer::MAX_DIM];
};
// End class NormaMinkHlaw

//===========================================================================

/*
 * This is (2/V_n)^(2/n) which seems wrong ?
 */
const double NormaMinkHlaw::m_gamma[] = {
/* GamMinkowski[0] = */0.00000000000000,
/* GamMinkowski[1] = */0.00000000000000,
/* GamMinkowski[2] = */1.04719756392815,
/* GamMinkowski[3] = */0.69062814882789,
/* GamMinkowski[4] = */0.66230588438641,
/* GamMinkowski[5] = */0.68895682802115,
/* GamMinkowski[6] = */0.73293650323622,
/* GamMinkowski[7] = */0.78407967437680,
/* GamMinkowski[8] = */0.83869147754813,
/* GamMinkowski[9] = */0.89517405737802,
/* GamMinkowski[10] = */0.95274948000853,
/* GamMinkowski[11] = */1.01100354118401,
/* GamMinkowski[12] = */1.06969923726369,
/* GamMinkowski[13] = */1.12869260847180,
/* GamMinkowski[14] = */1.18789173406828,
/* GamMinkowski[15] = */1.24723543612412,
/* GamMinkowski[16] = */1.30668158848970,
/* GamMinkowski[17] = */1.36620037299243,
/* GamMinkowski[18] = */1.42577021185918,
/* GamMinkowski[19] = */1.48537521407207,
/* GamMinkowski[20] = */1.54500351420012,
/* GamMinkowski[21] = */1.60464615768560,
/* GamMinkowski[22] = */1.66429633246555,
/* GamMinkowski[23] = */1.72394882699329,
/* GamMinkowski[24] = */1.78359964036430,
/* GamMinkowski[25] = */1.84324569710679,
/* GamMinkowski[26] = */1.90288463550686,
/* GamMinkowski[27] = */1.96251464853672,
/* GamMinkowski[28] = */2.02213436300702,
/* GamMinkowski[29] = */2.08174274687719,
/* GamMinkowski[30] = */2.14133903755970,
/* GamMinkowski[31] = */2.20092268604477,
/* GamMinkowski[32] = */2.26049331306133,
/* GamMinkowski[33] = */2.32005067447401,
/* GamMinkowski[34] = */2.37959463382297,
/* GamMinkowski[35] = */2.43912514042718,
/* GamMinkowski[36] = */2.49864221184901,
/* GamMinkowski[37] = */2.55814591979856,
/* GamMinkowski[38] = */2.61763637876553,
/* GamMinkowski[39] = */2.67711373682499,
/* GamMinkowski[40] = */2.73657816818394,
/* GamMinkowski[41] = */2.79602986712778,
/* GamMinkowski[42] = */2.85546904309690,
/* GamMinkowski[43] = */2.91489591667914,
/* GamMinkowski[44] = */2.97431071634673,
/* GamMinkowski[45] = */3.03371367580026,
/* GamMinkowski[46] = */3.09310503180911,
/* GamMinkowski[47] = */3.15248502245848,
/* GamMinkowski[48] = */3.21185388573071 };

/*=========================================================================*/

NormaMinkHlaw::NormaMinkHlaw(double logDensity, int64_t maxDim, NormType norm) :
        NormaMinkHlaw(maxDim, norm) {
    Normalizer::computeBounds(logDensity);
}

/*=========================================================================*/

NormaMinkHlaw::NormaMinkHlaw(double logm, int64_t k, int64_t maxDim, NormType norm) :
        NormaMinkHlaw(maxDim, norm) {
    Normalizer::computeBounds(logm, k);
}

/*=========================================================================*/

NormaMinkHlaw::NormaMinkHlaw(int64_t maxDim, NormType norm) :
        Normalizer(maxDim, norm) {
    if (maxDim > this->MAX_DIM)
        throw std::invalid_argument("NormaMinkHlaw:   dimension > MAXDIM");
    m_name = "NormaMinkHlaw";
}

/*=========================================================================*/

inline double NormaMinkHlaw::getGamma(int64_t j) const {
    if (j < 1 || j > this->MAX_DIM)
        throw std::out_of_range("NormaMinkHlaw::getGamma");
    if (m_norm == L2NORM)
       return m_gamma[j];
    else if (m_norm == L1NORM)
       return m_gamma[j] * j;
    else
       throw std::domain_error("NormaMinkHlaw::getGamma with wrong norm");
}

}

#endif

