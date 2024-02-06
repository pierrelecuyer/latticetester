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

#ifndef LATTICETESTER_RANK1LATTICE_H
#define LATTICETESTER_RANK1LATTICE_H

#include "latticetester/Util.h"
#include "latticetester/EnumTypes.h"
#include "latticetester/IntLatticeExt.h"

namespace LatticeTester {

    /**
     * This subclass of `IntLatticeExt` defines a general rank 1 lattice rule in \f$d\f$ dimensions,
     * whose points \f$\mathbb{u}_i\f$ are defined by
     * \f{equation}{
     *    \mathbf{u}_i = (i \mathbf{a} \mod m)/m,
     * \f}
     * where \f$\mathbf{a} \in \mathbb{Z}_m^d\f$ is called the generating vector.
     * The lattice is rescaled simply by removing the division by \f$m\f$.
     * Given \f$\mathbf{a}\f$, a basis of the rescaled (integer) lattice is given by
     * \f{align*}{
     *    \mathbf{v}_1 & = a \\
   *    \mathbf{v}_2 & = m \mathbf{e}_2 \\
   *    \vdots & \\
   *    \mathbf{v}_d & = m \mathbf{e}_v \\
   * \f}
     * where \f$\mathbf{e}_i\f$ is the \f$i^\text{th}\f$ unit vector.
     *
     * A condition that is often required when building a rank 1 lattice is that
     * \f$\gcd(a_i, m) = 1,\ \forall 1\leq i \leq d\f$. When this condition is
     * verified, each lower-dimensional projection of the lattice contains the same number
     * of points (has the same density \f$m\f$) as the full lattice.
     * When searching for lattices that satisfy this condition, one may assume
     * without loss of generality generality that \f$a_1 = 1\f$.
     * Under that condition, it is straightforward to construct a basis for a projection that
     * contains the first coordinate, and also its m-dual basis.  We exploit this.
     */

    template<typename Int, typename Real>
    class Rank1Lattice: public IntLatticeExt<Int, Real> {

    private:
        typedef NTL::vector<Int> IntVec;
        typedef NTL::matrix<Int> IntMat;
        typedef NTL::vector<Real> RealVec;

    public:

        /**
         * This constructor takes as input the modulus `m`, the generating vector `aa`,
         * and the norm used to measure the vector lengths.
         * The maximal dimension `maxDim` will be the length of the vector `aa`.
         * This constructor does not build the basis, to leave
         * more flexibility in the dimension when doing so.
         */
        Rank1Lattice (const Int & m, const IntVec & aa,
                bool withDual=false, NormType norm = L2NORM);

        /**
         * Constructor for the special case of a Korobov lattice.
         * Here the generating vector has the form aa = (1, a, a^2 mod m, a^3 mod m, ...)
         * where a is an integer such that 1 < a < m.
         */
        Rank1Lattice (const Int & m, const Int & a, int64_t maxDim,
                bool withDual=false, NormType norm = L2NORM);

        /**
         * This constructor does not specify the generating vector.
         */
        Rank1Lattice (const Int & m, int64_t maxDim,
                bool withDual=false, NormType norm = L2NORM);

        /**
         * Copy constructor.
         */
        Rank1Lattice (const Rank1Lattice<Int, Real> & Lat);

        /**
         * Assigns `Lat` to this object.
         */
        Rank1Lattice & operator= (const Rank1Lattice<Int, Real> & Lat);

        /**
         * Destructor.
         */
        ~Rank1Lattice();

        /**
         * Sets the generating vector to `aa`.  The maximum dimension will be set to its dimension.
         */
        void setaa (const IntVec & aa);

        /**
         * Sets this lattice to a Korobov lattice with multiplier `a`.
         * The maximum dimension is unchanged.
         */
        void seta (const Int & a);

        /**
         * Returns the first components of the generating vector \f$\ba\f$ as a string.
         * The number of components in the string will be the current dimension of the lattice.
         */
        std::string toStringCoef() const;

        /**
         * Builds a basis in `dim` dimensions. This `dim` must not exceed `this->maxDim()`.
         * This initial basis will be upper triangular.
         */
        void buildBasis (int64_t dim);

        /**
         * Builds only an m-dual lower triangular basis (and not the primal) directly
         * in `dim` dimensions.  This `dim` must not exceed `maxDim`.
         */
        void buildDualBasis (int64_t dim);

        /**
         * Increases the current dimension of the primal basis by 1 and updates the basis.
         * If `withDual`, it also increases the m-dual basis.
         * The new increased dimension must not exceed `maxDim`.
         */
        void incDimBasis ();

        /**
         * Increases the current dimension of only the m-dual basis by 1.
         * The primal basis is left unchanged (not updated).
         * The new increased dimension must not exceed `maxDim`.
         * This method uses the simplified method given in the lattice tester guide:
         * the new m-dual basis vector is simply  w_d = (-a_d, 0, ..., 0, 1).
         */
        void incDimDualBasis ();

        /**
         * This method overrides its namesake in `IntLattice`. The projection of this
         * `Rank1Lattice` over the coordinates in `proj` is returned in `projLattice`.
         * The implementation used here exploits the rank-1 lattice structure and it
         * is simpler and faster than the general one. See Section 5.5 of the guide.
         * When the first coordinate is 1 and belongs to the projection, both the primal
         * and m-dual constructions are direct, just by selecting the rows and columns
         * whose indices are in `proj`. Otherwise, those rows and columns plus the
         * first row form a set of `proj.size()+1` generating vectors for the primal.
         * This number must not exceed the `maxDim` of `projLattice`.
         * We reduce them to a primal basis, using LLL if we only want the primal.
         * If we want an m-dual basis, then we find an upper triangular basis for the primal,
         * and the corresponding lower-triangular m-dual basis.
         * In contrast to the same method in the parent class, this specialized version
         * does not require that a basis for the whole lattice has been constructed before.
         */
        void buildProjection (IntLattice<Int, Real> *projLattice,
                const Coordinates &proj, double delta = 0.99) override;

        /**
         * Returns the first `dim` components of the generating vector \f$\ba\f$ as a string,
         * where `dim` is the current lattice dimension.
         */
        //std::string toStringCoef() const;

    protected:

        /**
         * Vector of multipliers (generating vector) of the rank 1 lattice rule.
         * They are stored for up to `maxDim()` dimensions.
         * The first coordinate has index 0.
         */
        IntVec m_a;
    };

//============================================================================

    template<typename Int, typename Real>
    Rank1Lattice<Int, Real>::Rank1Lattice (
            const Int & m, int64_t maxDim, bool withDual, NormType norm):
    IntLatticeExt<Int, Real> (m, maxDim, withDual, norm) {
        this->m_a.SetLength(maxDim);
    }

//============================================================================

    template<typename Int, typename Real>
    Rank1Lattice<Int, Real>::Rank1Lattice (
            const Int & m, const IntVec & aa, bool withDual, NormType norm):
    IntLatticeExt<Int, Real> (m, aa.length(), withDual, norm) {
        this->m_maxDim = aa.length();
        this->m_a = aa;
    }

//============================================================================

    template<typename Int, typename Real>
    Rank1Lattice<Int, Real>::Rank1Lattice (
            const Int & m, const Int & a, int64_t maxDim, bool withDual, NormType norm):
    IntLatticeExt<Int, Real> (m, maxDim, withDual, norm) {
        this->m_a.SetLength(maxDim);
        this->seta(a);    // Set to a Korobov lattice.
    }

//============================================================================

    template<typename Int, typename Real>
    Rank1Lattice<Int, Real>::~Rank1Lattice() {
        this->m_a.kill ();
    }

    //============================================================================

    /*
     template<typename Int, typename Real>
     void Rank1Lattice<Int, Real>::initProj() {
     IntLatticeExt<Int, Real>::initProj();
     // for (int64_t r = 1; r < this->getDim(); r++)
     //   this->m_lgVolDual2[r] = this->m_lgVolDual2[r - 1];
     }
     */

    //============================================================================
    template<typename Int, typename Real>
    Rank1Lattice<Int, Real> & Rank1Lattice<Int, Real>::operator= (
            const Rank1Lattice<Int, Real> & lat) {
        if (this == &lat) return *this;
        this->copy (lat);
        this->m_a = lat.m_a;
        return *this;
    }

    //============================================================================

    template<typename Int, typename Real>
    Rank1Lattice<Int, Real>::Rank1Lattice (const Rank1Lattice<Int, Real> & lat):
    IntLatticeExt<Int, Real> (
            lat.m_modulo, lat.getDim (), lat.m_withDual, lat.getNormType ()) {
        this->m_a = lat.m_a;
        // Should also copy the basis and all other variables!
    }

    //============================================================================

    /**
     * Sets the generating vector to `aa`.
     */
    template<typename Int, typename Real>
    void Rank1Lattice<Int, Real>::setaa (const IntVec & aa) {
        this->m_a = aa;
    }

    //============================================================================

    /**
     * Sets this lattice to a Korobov lattice with multiplier `a`, with dim = maxDim.
     */
    template<typename Int, typename Real>
    void Rank1Lattice<Int, Real>::seta (const Int & a) {
        this->m_a[0] = 1;
        for (int64_t i=1; i < this->m_maxDim; i++) {
            this->m_a[i] = (a * this->m_a[i-1]) % this->m_modulo;
        }
    }

    //============================================================================

    // An upper-triangular basis is built directly, as explained in the guide of Lattice Tester.
    // The dimension `maxDim` of the `IntMat` array is unchanged.
    // In case `withDual` is true, the m-dual basis is also constructed directly.

    template<typename Int, typename Real>
    void Rank1Lattice<Int, Real>::buildBasis (int64_t d) {
        assert(d <= this->m_maxDim);
        this->m_dim = d;     // The dimension of the lattice is set to d.
        int64_t i, j;

        // This builds an upper-triangular basis in a standard way.
        for (j = 0; j < d; j++) {
            this->m_basis[0][j] = this->m_a[j];
        }
        for (i = 1; i < d; i++) {
            for (j = 0; j < d; j++) {
                if (i == j) this->m_basis[i][i] = this->m_modulo;
                else this->m_basis[i][j] = 0;
            }
        }
        this->setNegativeNorm ();
        if (!this->m_withDual) return;  // We are done.

        // If `withDual`, we construct the m-dual basis also in a direct way.
        this->m_dualbasis[0][0] = this->m_modulo;
        for (j = 1; j < d; j++) this->m_dualbasis[0][j] = 0;
        for (i = 1; i < d; i++) {
            this->m_dualbasis[i][0] = -this->m_a[i];
            for (j = 1; j < d; j++) {
                if (i == j) this->m_dualbasis[i][i] = 1;
                else this->m_dualbasis[i][j] = 0;
            }
        }
        this->setDualNegativeNorm ();
    }

    //============================================================================

    // This one builds only the m-dual basis, also in a direct way.
    template<typename Int, typename Real>
    void Rank1Lattice<Int, Real>::buildDualBasis (int64_t d) {
        assert(d <= this->m_maxDim);
        this->setDim (d);
        int64_t i, j;
        this->m_dualbasis[0][0] = this->m_modulo;
        for (j = 1; j < d; j++) this->m_dualbasis[0][j] = 0;
        for (i = 1; i < d; i++) {
            this->m_dualbasis[i][0] = -this->m_a[i];
            for (j = 1; j < d; j++) {
                if (i == j) this->m_dualbasis[i][i] = 1;
                else this->m_dualbasis[i][j] = 0;
            }
        }
        this->setDualNegativeNorm ();
    }

//============================================================================

    template<typename Int, typename Real>
    void Rank1Lattice<Int, Real>::incDimBasis () {
        int64_t d = 1 + this->getDim();  // New current dimension.
        assert(d <= this->m_maxDim);
        this->setDim (d);
        int64_t i, j;
        Int m_add;

        // Update new row and new column.
        for (j = 0; j < d-1; j++) this->m_basis[d-1][j] = 0;
        this->m_basis[d-1][d-1] = this->m_modulo;
        for (i = 0; i < d-1; i++) {
            this->m_basis[i][d-1] = (this->m_a[d-1] * this->m_basis[i][0]) % this->m_modulo;
        }
        this->setNegativeNorm ();

        // If m-dual basis is maintained, we also increase its dimension.
        if (this->m_withDual) {
            // Add extra coordinate to each vector
            for (i = 0; i < d; i++) {
                this->m_dualbasis[i][d-1] = 0;
                this->m_dualbasis[d-1][i] = 0;
            }
            this->m_dualbasis[d-1][d-1] = 1;
            for (j = 0; j < d-1; j++) {
                m_add = 0;
                for (int i = 0; i < d-1; i++) {
                    m_add = m_add - this->m_basis[i][d-1] * this->m_dualbasis[i][j];
                }
                m_add = m_add / this->m_modulo;
                this->m_dualbasis[d-1][j] = this->m_dualbasis[d-1][j] + m_add;
            }
            this->setDualNegativeNorm ();
        }
    }

    //============================================================================

    template<typename Int, typename Real>
    void Rank1Lattice<Int, Real>::incDimDualBasis () {
        int64_t d = 1 + this->getDim();
        assert(d <= this->m_maxDim);
        this->m_dim = d;
        int64_t i;
        // Add one extra coordinate to each vector.
        for (i = 0; i < d; i++) {
            this->m_dualbasis[i][d-1] = 0;
            this->m_dualbasis[d-1][i] = 0;
        }
        this->m_dualbasis[d-1][0] = -m_a[d-1];
        this->m_dualbasis[d-1][d-1] = 1;
        this->setDualNegativeNorm ();
    }

    //============================================================================
    template<typename Int, typename Real>
    void Rank1Lattice<Int, Real>::buildProjection(IntLattice<Int, Real> *projLattice,
            const Coordinates &proj, double delta) {
    	
        // We use the method described in the Lattice Tester guide, section 5.5.
        // Does not assume that a basis for `this` has been computed before.
        // bool case1 = proj.contains(1) && m_a[0] == 1;// First coord. selected and a_1 = 1.
        bool case1 = m_a[*proj.begin()-1] == 1;
        long i, j;
        long d = proj.size();// Number of coordinates in the projection.
        projLattice->setDim (d);
        Int c1, b1, b2;      
        IntMat &basis = projLattice->getBasis();
        IntMat &dualBasis = projLattice->getDualBasis();

        if (projLattice->withPrimal()) { // Build a primal basis.
            if (case1) {
                j = 0; // CW
                for (auto it = proj.begin(); it != proj.end(); it++, j++) {
                    basis[0][j] = m_a[*it - 1];  // First row.
                }
            }
            else {
                // XGCD (g, c, d, a, b) does g = gcd(a, b) = a*c + b*d.
                NTL::XGCD(c1, b1, b2, m_a[*proj.begin()-1], this->m_modulo);
                j = 0; //CW
                for (auto it = proj.begin(); it != proj.end(); it++, j++) {                	
                    // NTL::MulMod (projLattice->m_basis[0][j], m_a[*it - 1], b1, this->m_modulo); // First row - does not work here
                    basis[0][j] = m_a[*it - 1] * b1 % this->m_modulo;
                }
                basis[0][0] = c1;
            }
            for (i = 1; i < d; i++) {
                for (j = 0; j < d; j++) {
                    if (i == j) basis[i][i] = this->m_modulo;
                    else basis[i][j] = 0;
                }
            }
        }
        if (projLattice->withDual()) { // Compute m-dual basis directly.
            if (case1) {
                i = 0; // CW 
                for (auto it = proj.begin(); it != proj.end(); ++it, ++i) {
                    dualBasis[i][0] = m_a[*it - 1]; // First column.
                }
                dualBasis[0][0] = this->m_modulo;
            }
            else {
            	if (!projLattice->withPrimal())
                     NTL::XGCD(c1, b1, b2, m_a[*proj.begin()-1], this->m_modulo);   //  c1 and b1 were not computed.
            	j = 0; // CW
            	for (auto it = proj.begin(); it != proj.end(); it++, j++) {
            	   dualBasis[0][j] = - m_a[*it - 1] * b1 / c1;  // Warning: If c1 > 1, this may not be an integer.
            	}
            	dualBasis[0][0] = this->m_modulo / c1;
            }
            for (i = 0; i < d; i++) {
                for (j = 1; j < d; j++) {
                    if (i == j) dualBasis[i][i] = 1;
                    else dualBasis[i][j] = 0;
                }
            }
        }
        
    }

//============================================================================
    template<typename Int, typename Real>
    std::string Rank1Lattice<Int, Real>::toStringCoef() const {
        return toString(this->m_a, 0, this->getDim());
    }

//============================================================================

    template class Rank1Lattice<std::int64_t, double>;
    template class Rank1Lattice<NTL::ZZ, double>;
    template class Rank1Lattice<NTL::ZZ, NTL::RR>;

} // End namespace LatticeTester

#endif
