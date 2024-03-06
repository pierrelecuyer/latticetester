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

#ifndef LATTICETESTER_REDUCER_H
#define LATTICETESTER_REDUCER_H

#include "NTL/LLL.h"
#include "NTL/tools.h"
#include "NTL/ZZ.h"
#include "NTL/RR.h"

#include "latticetester/EnumTypes.h"
#include "latticetester/Util.h"
#include "latticetester/IntLattice.h"
#include "latticetester/BasisConstruction.h"
#include "latticetester/NTLWrap.h"
#include "latticetester/LLL_FPInt.h"

#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <cstdint>
#include <iostream>
#include <ctime>
#include <string>
#include <iomanip>
#include <set>
#include <map>
#include <cmath>
#include <cstdlib>
#include <type_traits>

namespace LatticeTester {

/**
 * This file provides only static functions to reduce a lattice basis in some ways
 * (pairwise, LLL, BKZ \cite rDIE75a, \cite mLEN82a, \cite mSCH91a).
 * The LLL and BKZ ones are wrappers to slightly modified versions of the NTL
 * functions described at https://libntl.org/doc/LLL.cpp.html .
 * These functions do not require the creation of an object like in `ReducerBB`.
 * They take an `IntMat` object that contains a basis and return the reduced basis
 * in the same object. In our versions, the basis can occupy only part of the `IntMat`
 * object so we can use the same object for several bases of various sizes,
 * the basis entries can be of either ZZ or `int64_t` type,
 * the shortest basis vector is always placed in the first row,
 * and the squared vector lengths are also returned in a vector.
 * The norm to measure the vector lengths is always taken as the Euclidean one.
 * The computations can be done either in `double` or in `RR` for the real numbers.
 *
 * All these functions take as input (and return as output) an `IntMat` object
 * whose first `dim` rows and columns (a square matrix) are assumed to form
 * be a set of independent vectors that form a basis for the lattice.
 * The same `IntMat` object can then be used for several lattices of different sizes.
 */

template<typename Int, typename Real>
class ReducerStatic {

// using namespace LatticeTester;

private:
    // Local typedefs for matrix and vector types needed in the class.
    typedef NTL::vector<Int> IntVec;
    typedef NTL::matrix<Int> IntMat;
    typedef NTL::vector<Real> RealVec;
    typedef NTL::matrix<Real> RealMat;

public:

     /**
     * This method performs pairwise reduction sequentially on all vectors
     * of the basis whose indices are greater of equal to `dim >=0`,
     * as proposed in \cite rDIE75a.
     * The boolean vector `taboo[]` is used internally in Minkowski reduction only:
     * when this vector is not `NULL`, `taboo[j]=true' means that the j-th vector
     * should not be modified.
     */
    static void redDieter(int64_t dim, bool taboo[] = NULL);

    /**
     * This is the NTL implementation of the floating point version of the
     * LLL reduction algorithm presented in \cite mSCH91a.
     * The factor `delta` has the same meaning as in `redLLLOld`.
     * The parameter `precision` specifies the precision of the floating point numbers
     * that the algorithm will use. `EnumTypes.h` provides a list of the possible values,
     * and their description is done in the module `LLL` of NTL.
     * The reduction is always applied to the entire basis (all vectors and coordinates).
     */

    /**
     * A static version of the same function for which the basis is passed as a parameter.
     * In this version, when `dim > 0`, the reduction is applied only to the square submatrix
     * comprised of the first `dim` rows and `dim` columns of the `basis` matrix, which is actually
     * considered as the basis. The other rows and columns are ignored.
     * This is useful if one wants to reuse the same `basis` object several times
     * to pass lattice bases having different dimensions.
     * When `dim = 0`, the number of columns of `gen` is used for the dimension.
     * Here, the precision type is always `DOUBLE`.
     */
    static void redLLLNTL(IntMat &basis, double delta = 0.999999,
            long dim = 0, Vec<double> *sqlen = 0, PrecisionType precision = DOUBLE);

    /**
     * This static function implements an exact algorithm from NTL to perform the original LLL reduction.
     * This is slower than `redLLLNTL`, but more accurate.
     * There is also a static version which takes the basis as input.
     */
    // void redLLLNTLExact(double delta = 0.999999);

    static void redLLLNTLExact(IntMat &basis, double delta = 0.999999
            long dim = 0, Vec<double> *sqlen = 0);

    /**
     * This calls the NTL implementation of the floating point version of the
     * BKZ reduction algorithm presented in \cite mSCH91a,
     * with reduction factor `delta` and block size `blocksize`; see \cite iLEC22l.
     * The factor `delta` has a similar meaning as in `redLLLOld`.
     * The `precision` and `dim` parameters have the same meaning as in `redLLLNTL`.
     * The parameter `blocksize` gives the size of the blocks in the BKZ
     * reduction. Roughly, larger blocks means a stronger condition.
     * A `blocksize` of 2 is equivalent to LLL reduction.
     * There is also a static version that takes the lattice as input.
     */
//    void redBKZ(double delta = 0.999999, int64_t blocksize = 10,
//            long dim=0, double *sqlen=0, PrecisionType prec = DOUBLE);

    /**
     * A static version of the previous method.  The lattice basis is passed as a parameter.
     */
    static void redBKZ(IntMat &basis, double delta = 0.999999, int64_t blocksize = 10,
            long prune = 0, long dim = 0, Vec<double> *sqlen = 0, PrecisionType prec = DOUBLE);

};


//============================================================================
// Implementation

//===========================================================================


/*
// This is the general implementation.
template<typename Int, typename Real>
void Reducer<Int, Real>::redLLLNTL(double delta, PrecisionType precision) {
    redLLLNTL(m_lat->getBasis(), delta, m_lat->getDim(), 0, precision);
    //        m_lat->getVecNorm(), precision);
    // MyExit(1, "redLLLNTL only works with integers.");
}

// A specialization for the case where Int = ZZ.
template<typename Real>
void redLLLNTL(Reducer<NTL::ZZ, Real> &red, double delta,
        PrecisionType precision) {
    IntLattice<NTL::ZZ>* lat = red.getIntLattice();
    redLLLNTL(lat->getBasis(), delta, lat->getDim(), 0, precision);
}
*/

// Static version: general.
template<typename Int, typename Real>
void Reducer<Int, Real>::redLLLNTL(IntMat &basis, Real delta,
        long dim, RealVec *sqlen, PrecisionType precision) {
    if (precision == DOUBLE) {
        NTL::LLL_FPInt(basis, delta, dim, dim, sqlen);
        return;
    }
    MyExit(1, "redLLLNTL: Int = int64_t with precision != DOUBLE.\n");
}
MyExit(1, "redLLLNTL: General version not implemented.\n");
}

// Static version: a specialization for the case where Int = int64_t.
template<typename Real>
void Reducer<Int, Real>::redLLLNTL(NTL::matrix<int64_t> &basis, double delta,
        long dim, Vec<double> *sqlen, PrecisionType precision) {
    if (precision == DOUBLE) {
        // NTL::LLL_FPZZflex(basis, delta, dim, dim, sqlen);
        NTL::LLL_FPInt(basis, delta, dim, dim, sqlen);
        return;
    }
    MyExit(1, "redLLLNTL: Int = int64_t with precision != DOUBLE.\n");
}

// Static version: a specialization for the case where Int = ZZ.
// See `https://github.com/u-u-h/NTL/blob/master/doc/LLL.txt` for details
// about the `PrecisionType` choices.
template<typename Real>
void Reducer<Int, Real>::redLLLNTL(NTL::matrix<NTL::ZZ> &basis, double delta,
        long dim, Vec<double> *sqlen, PrecisionType precision) {
    if (precision == DOUBLE) {
        // NTL::LLL_FPZZflex(basis, delta, dim, dim, sqlen);
        NTL::LLL_FPInt(basis, delta, dim, dim, sqlen);
        return;
    }
    // If precision is not DOUBLE.
    NTL::matrix<NTL::ZZ> cpbasis;
    case QUADRUPLE:
        cpbasis.SetDims (dim, dim);
        copy (basis, cpbasis, dim, dim);  // From Util
        NTL::LLL_QP(cpbasis, delta);
        copy (cpbasis, basis, dim, dim);
        break;
    case XDOUBLE:
        cpbasis.SetDims (dim, dim);
        copy (basis, cpbasis, dim, dim);  // From Util
        NTL::LLL_XD(cpbasis, delta);
        copy (cpbasis, basis, dim, dim);
        break;
    case RR:
        NTL::LLL_RRext(basis, delta, dim, dim, sqlen);
        break;
    default:
        MyExit(1, "redLLLNTL: undefined PrecisionType.");
    }
}



//=========================================================================

// Static version for the general case.
template<typename Int, typename Real>
void Reducer<Int, Real>::redLLLNTLExact(IntMat &basis, double delta) {
    MyExit(1, "redLLLNTLExact works only for Int = ZZ");
}

// Static version, for Int = ZZ.
template<>
void Reducer<Int, Real>::redLLLNTLExact(NTL::matrix<NTL::ZZ> &basis,
        double delta) {
    NTL::ZZ det(0);
    int64_t denum;
    denum = round(1.0 / (1.0 - delta)); // We want (denum-1)/denum \approx delta.
    NTL::LLL(det, basis, denum - 1, denum);
}

//=========================================================================


// Static version, general.
void redBKZ(IntMat &basis, double delta, int64_t blocksize,
            long prune, long r, long c, Vec<double> *sqlen, PrecisionType prec);



// Static version: specialization for Int = int64_t.
template<typename Real>
void Reducer<int64_t, Real>::redBKZ(NTL::matrix<int64_t> &basis, double delta,
        std::int64_t blocksize, long dim, double *sqlen, PrecisionType precision) {
    if (precision == DOUBLE) {
        NTL::BKZ_FPInt(basis, delta, blocksize, dim, dim, sqlen);
        return;
    }
    MyExit(1, "redBKZ: Int = int64_t with precision != DOUBLE.\n");
}

// Static version, for Int = ZZ.
template<typename Real>
void Reducer<NTL::ZZ, Real>::redBKZ(NTL::matrix<NTL::ZZ> &basis, double delta,
        std::int64_t blocksize, long dim, double *sqlen, PrecisionType precision) {
    if (precision == DOUBLE) {
        // NTL::LLL_FPZZflex(basis, delta, dim, dim, sqlen);
        NTL::BKZ_FPInt(basis, delta, blocksize, dim, dim, sqlen);
        return;
    }
    // If precision is not DOUBLE, dim is taken as the dimension of basis.
    NTL::matrix<NTL::ZZ> cpbasis;
    // if ((dim > 0) & (dim != basis.NumRows())) {
        cpbasis.SetDims (dim, dim);
        LatticeTester::copy (basis, cpbasis, dim, dim);  // From Util
    //} else
    //    cpbasis = &basis;
    switch (precision) {
    case DOUBLE:
        NTL::BKZ_FPInt(basis, delta, blocksize, dim, dim, sqlen);
        return;
    case QUADRUPLE:
        NTL::BKZ_QP(cpbasis, delta, blocksize);
        break;
    case XDOUBLE:
        NTL::BKZ_XD(cpbasis, delta, blocksize);
        break;
    case RR:
        NTL::BKZ_RR(cpbasis, delta, blocksize);
        break;
    default:
        MyExit(1, "Undefined precision type for redBKZ");
    }
    LatticeTester::copy (cpbasis, basis);
}

//============================================================================

template class Reducer<std::int64_t, double> ;
template class Reducer<NTL::ZZ, double> ;
template class Reducer<std::int64_t, NTL::RR> ;
template class Reducer<NTL::ZZ, NTL::RR> ;

}   // namespace LatticeTester

#endif 
