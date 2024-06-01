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

#ifndef LATTICETESTER_NORMALIZER_H
#define LATTICETESTER_NORMALIZER_H

#include "latticetester/Util.h"
#include "latticetester/EnumTypes.h"

#include <cassert>
#include <string>
#include <typeinfo>

namespace LatticeTester {

/**
 * This is a base class for implementing normalization constants used in figures of merit,
 * to normalize the length of the shortest nonzero vector in either the primal or dual lattice.
 * These constants are based on upper bounds (or approximations) on the best possible length,
 * for given lattice density and dimension.
 * The constructors in the various subclasses compute the specific bounds,
 * as explained in Section 9 of the user's guide.
 *
 * The preferred constructor assumes that the rescaled primal lattice has scaling factor \f$m\f$
 * and order \f$k\f$, so its density is \f$\min(1, m^{k-t})\f$ in \f$t\geq 1\f$ dimensions,
 * and cannot exceed 1 for projections in \f$s < k\f$ dimensions.
 * This constructor take \f$ \log m \f$ and  \f$k\f$ as inputs.
 * It invokes `computeBounds` internally to compute the bounds.
 * We work with the log density instead of the density itself
 * because the latter is sometimes extremely large or extremely small.
 *
 * Given a (non-rescaled) lattice of density \f$ \eta_t \f$ in \f$t\f$ dimensions,
 * the Euclidean length \f$ d_t \f$ of a shortest nonzero lattice vector is upper-bounded as follows:
 * \f[
 *    d_t \le  d_t^*(\eta_t) \eqdef \gamma_t^{1/2} \eta_t^{-1/t},
 * \f]
 * where the constants  \f$ \gamma_t \f$ are known exactly only for  \f$ t\leq 8\f$.
 * For larger values of  \f$ t \f$, we use bounds or approximations of these constants,
 * which are defined in subclasses of `Normalizer`.
 * These bounds also hold for the  \f$ L^1\f$ lengths of the shortest vectors,
 * but with different constants  \f$ \gamma_t \f$.
 * In the dual lattice, the density is \f$1/\eta_t\f$ instead, and the Euclidean length
 * \f$ \ell_t \f$ of a shortest nonzero lattice vector then obeys:
 * \f[
 *    \ell_t \le  \ell_t^*(\eta_t) \eqdef \gamma_t^{1/2} \eta_t^{1/t}.
 * \f]
 * The bounds \f$ d_t^*(\eta_t)\f$ and \f$ \ell_t^*(\eta_t)\f$ are used to normalize
 * the values of \f$ d_t \f$ or \f$ \ell_t \f$ computed by the software for given lattices
 * to obtain standardized measures that lie between 0 and 1.
 * For the best lattices, these measures should be close to 1.
 *
 * In subclasses, it suffices to implement the `getGamma(int) const` function
 * and the specialized constructors.
 * The constructors in the present abstract class do not compute any bounds,
 * they basically just reserve the space (an array) for the bounds.
 *
 * The preferred usage for this class is to declare a pointer to a Normalizer
 * and to instantiate a subclass with dynamically allocated memory:
 * \code{.cpp}
 *   NormaBestLat norma(logm, k, maxDim);
 * \endcode
 *
 * When making a search and examining millions of lattices, it is important to re-use
 * the same `Normalizer` object, and *not* construct a new one (and recompute the constants)
 * for each lattice that is examined.
 */

//template<typename Int>
class Normalizer {

public:
    /**
     * The maximum dimension of the lattices for which this class can give
     * an upper bound.
     */
    static const int64_t MAX_DIM = 48;

    /**
     * This is the preferred constructor in all the subclasses.
     * It assumes that the rescaled primal lattice has scaling factor \f$m\f$
     * and order \f$k\f$, so its density is \f$m^{k-t}\f$ for \f$t\geq k\f$, and cannot
     * exceed 1 for projections in \f$s < k\f$ dimensions.
     * The bounds \f$ d_t^*(\eta)\f$ will then be computed by assuming those densities.
     * To compute bounds for the m-dual, pass `-logm` instead of `logm`.
     * The values of \f$\log m\f$ (in natural basis) and \f$k\f$ must be given as inputs.
     */
    Normalizer(double logm, int64_t k, int64_t maxDim, std::string name = "",
            NormType norm = L2NORM);

    /**
     * *****  DEPRECATED  *****
     *
     * This old (legacy) constructor assumes that the lattice is not rescaled.
     * It creates a `Normalizer` by assuming that the density is
     * \f$\eta=\exp(\text{logDensity})\f$ in all dimensions.
     * This is (probably) never true when there is rescaling.
     * The bounds will be computed in up to `maxDim` dimensions.
     * To compute bounds for the dual, use `-logDensity` instead of `logDensity`.
     * Only subclasses can actually compute the bounds.
     * The `name` parameter gives name that will be printed by the ToString() method.
     * The `norm` parameter is the `NormType` used by this object.
     */
    Normalizer(double logDensity, int64_t maxDim, std::string name = "",
            NormType norm = L2NORM);

    /**
     * This constructor creates a Normalizer object without computing any bounds.
     * This is used in the case of rank 1 lattices in the `NormaPalpha`
     * class with a prime density.
     */
    Normalizer(int64_t maxDim, std::string name, NormType norm = L2NORM);

    /**
     * Destructor.
     */
    virtual ~Normalizer() {
        delete[] m_bounds;
    }

    /**
     * This method is rarely applicable!
     * It computes bounds by assuming that the log density is `logDensity`
     * for all dimensions up to the maximal dimension `maxDim`.
     * It ignores all rescaling.  This will become the new `logDensity` in this object.
     * To compute bounds for the dual, use `-logDensity` instead of `logDensity`.
     * This method is called by the constructors of subclasses, in which the constants
     * `gamma_t` are known.
     */
    virtual void computeBounds(double logDensity);

    /**
     *  This method is called by the constructors of subclasses, in which the constants
     * `gamma_t` are known. It computes the bounds that this normalizer will return, by assuming that
     * the primal lattice was rescaled by the factor \f$m\f$ and has order \f$k\f$,
     * so its density is \f$m^{k-t}\f$ for \f$t\geq k\f$, and cannot
     * exceed \f$m^s\f$ for projections in \f$s < k\f$ dimensions.
     * The values of \f$\log m\f$ (in natural basis) and \f$k\f$ must be given as inputs.
     * To compute bounds for the m-dual, just pass `-logm` instead of `logm`.
     * The bounds can be retrieved via `getBounds()` or `getBound(j)`.
     */
    void computeBounds(double logm, int64_t k);

    /**
     * Returns a string that describes this object.
     */
    std::string toString() const;

    /**
     * Returns the norm associated with this object.
     */
    NormType getNorm() const {
        return m_norm;
    }

    /**
     * Returns the maximal dimension for this object. This is the `maxDim`
     * parameter of the constructors.
     */
    int64_t getMaxDim() const {
        return m_maxDim;
    }

    /**
     * Returns the array that contains the bounds on the lengths of the shortest nonzero vectors.
     * Element `j` of the array is for `j` dimensions (element 0 is not used).
     * These bounds must have been computed before.
     */
    virtual double* getBounds() const;

    /**
     * Returns the bound on the length of the shortest nonzero vector in `j` dimensions.
     */
    virtual double getBound(int64_t j) const;

    /**
     * Returns the value of a lattice constant \f$\gamma\f$ in
     * dimension \f$j\f$. These constants can be used by subclasses to
     * implement `computeBounds()`. For this base class, always returns 1.
     */
    virtual double getGamma(int64_t j) const;

protected:
    /**
     * Name of the normalizer.
     */
    std::string m_name;

    /**
     * Norm type associated with this object.
     */
    NormType m_norm;

    /**
     * Maximum dimension. Only elements 1 to <tt>m_maxDim</tt> (inclusive) of m_bounds below
     * will be pre-computed.
     */
    int64_t m_maxDim;

    /**
     * Contains the bounds on the length of the shortest nonzero vector in
     * the lattice in each dimension up to `maxDim`. This array is initialized by the
     * `computeBounds` methods and its values can be accessed via `getBound`.
     */
    double *m_bounds;

private:
    /**
     * Use of the copy-constructor is forbidden.
     */
    Normalizer(const Normalizer&);

    /**
     * Use of assignment is forbidden.
     */
    Normalizer& operator=(const Normalizer&);

    /**
     * Return the min between k et j
     */
    // int64_t min(int64_t k, int64_t j);

}
;
// End class Normalizer

//===========================================================================

Normalizer::Normalizer(int64_t maxDim, std::string name, NormType norm) :
        m_name(name), m_norm(norm), m_maxDim(maxDim) {
    m_bounds = new double[maxDim + 1];
}

/*-------------------------------------------------------------------------*/

// This one makes sense only when the density is the same in all dimensions.
void Normalizer::computeBounds(double logDensity) {
    double x;
    for (int64_t j = 1; j <= m_maxDim; j++) {
        x = 0.5 * log(getGamma(j)) - logDensity / j;
        m_bounds[j] = exp(x);
    }
}

/*-------------------------------------------------------------------------*/

void Normalizer::computeBounds(double logm, int64_t k) {
    // The bounds \f$ d_t^*(\eta_t/m^t)\f$ (for primal) or \f$ \ell_t^*(\eta_t)\f$ (for dual).
    // These are the bounds when the primal lattice is *rescaled* by a factor m.
    double x;
    for (int64_t j = 1; j <= m_maxDim; j++) {
        x = 0.5 * log(getGamma(j)) - logm * ((k < j) ? (double) k / j : 1);
        // x = 0.5 * log(getGamma(j)) - logm * min(k/j, 1);
        if (logm > 0)   // Primal lattice: we need to add logm.
            x += logm;
        m_bounds[j] = exp(x);
    }
}

/*-------------------------------------------------------------------------*/
//int64_t Normalizer::min(int64_t k, int64_t j) {
//    return (k < j) ? k : j;
//}

/*-------------------------------------------------------------------------*/

std::string Normalizer::toString() const {
    std::ostringstream os;
    os << "-----------------------------\n"
            << "Content of Normalizer object:\n\n Normalizer = " << m_name;
    //  os << "\n n = " << exp(m_logDensity);
    //  os << "(log(n) = " << m_logDensity << ")";
    os << "\n\n";

    //   os.setf(std::ios::left);
    os << std::setprecision(13);
    for (int64_t t = 1; t <= m_maxDim; t++) {
        os << " Bound[" << std::setw(2) << std::right << t << "] = "
                << std::setw(14) << std::left << m_bounds[t] << "\n";
    }
    os << "\n";
    return os.str();
}

/*-------------------------------------------------------------------------*/

double Normalizer::getGamma(int64_t) const {
    // In this abstract class, the gamma_t's are undefined.
    return -1.0;
}

/*-------------------------------------------------------------------------*/

double Normalizer::getBound(int64_t j) const {
    return m_bounds[j];
}

/*-------------------------------------------------------------------------*/

double* Normalizer::getBounds() const {
    return m_bounds;
}
}
// end namespace LatticeTester

#endif
