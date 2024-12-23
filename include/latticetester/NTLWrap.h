// This file is part of LatticeTester.
//
// LatticeTester
// Copyright (C) 2012-2018  Pierre L'Ecuyer and Universite de Montreal
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

#ifndef LATTICETESTER__NTLWRAP_H
#define LATTICETESTER__NTLWRAP_H

#include <cstdint>
#include <cmath>

#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>

using namespace NTL;

/**
 * \file This module extends the `Vec` and `Mat` classes of NTL. It was previously
 * necessary because NTL and boost (an old dependency) did not use the same 
 * function names and indices.
 * Most of the methods below only define alias names for some NTL methods.
 * This was done at the time to have the same names for methods that are in both
 * boost and NTL, allowing LatticeTester to work with either the boost or NTL library,
 * depending on pre-processing statements.  These alias could probably be removed
 * because we no longer use boost in LatticeTester.
 * However, they may still be used in LatNet Builder (?)              **************
 *
 * New functions have also been implemented in this module as a way to overload a
 * few operators and methods of NTL (especially on matrix and vector types) to
 * the usage of `NTL::Mat<std::int64_t>` because some basic utilities that we need
 * for those integers are not offered in NTL.
 */

namespace NTL {


/**
 * \name Pointers to NTL matrix rows.
 * @{
 * An extension of `NTL::Vec<T>` implemented in this module to be used as
 * a reference to a matrix row. This is essentially the same as defining a
 * variable which is a pointer to the row vector.
 */
template <class M>
  class Mat_row : public Vec<typename M::value_type> {
    public:
      /**
       * No new vector is created, only a reference to the row.
       */
      inline Mat_row(M& data, long i) {
        this->_vec__rep = (typename M::value_type*&) data[i]._vec__rep;
      }
      /**
       * Empty constructor.
       * */
      inline ~Mat_row() {
        this->_vec__rep = 0; /* avoid destruction in parent class */
      }
  };


//============================================================================

/**
 * \name Conversion functions for compatibility with NTL.
 * @{
 * These functions perform conversions between different types. Most of them
 * do not really need explanations, but sometimes a specific logic is used
 * when doing the conversion.
 */

/**
 * Converts the array of characters (string) `c` into an `std::int64_t` `l`
 * using the strtol() function of cstdlib.h.
 */
inline void conv(std::int64_t &l, const char *c) {
    l = strtol(c, (char**) NULL, 10);
}

/**
 * Converts the array of characters (string) `c` into a `double` `r` using the
 * strtod() function of cstdlib.h.
 */
inline void conv(double &r, const char *c) {
    r = strtod(c, (char**) NULL);
}

/**
 * Converts a `int64_t` to a `double`.
 *
 inline void conv(double &x, int64_t a) {
 x = static_cast<double>(a);
 }

 * Converts a `double` to a `int64_t`. This truncates the decimals of a.
 *
 inline void conv(int64_t &x, double a) {
 x = static_cast<int64_t>(a);
 }

 * Converts a `int64_t` to a `NTL::ZZ`.
 *
 inline void conv(ZZ &x, int64_t a) {
 x = NTL::conv<ZZ>(a);
 }

 * Converts a `NTL::ZZ` to a `int64_t`. This will truncate a if it has to
 * many digits.
 *
 inline void conv(int64_t &x, ZZ a) {
 x = NTL::conv<int64_t>(a);
 }

 * Converts a `int64_t` to a `int64_t`. This will truncate a if it has to
 * many digits.
 *
 inline void conv(int64_t &x, int64_t a) {
 x = static_cast<int64_t>(a);
 }

 * Converts a `int64_t` to a `int64_t`.
 *
 inline void conv(int64_t &x, int64_t a) {
 x = static_cast<int64_t>(a);
 }

 * Since both are of the same type, this assigns a to x.
 *
 inline void conv(int64_t &x, int64_t a) {
 x = a;
 }

 inline void conv(int64_t &x, NTL::RR a) {
 x = static_cast<int64_t>(NTL::conv<double>(a));
 }
 */

/**
 * @}
 * \name Function overloads
 * @{
 * These functions are already implemented in NTL for NTL::ZZ and NTL::RR
 * types, but not for the other standard types used below. These overloads allow
 * us to make a simple call to the function in the `NTL` namespace without
 * worrying about types and still have working algorithms.
 */

/**
 * Returns the `bool` resulting of the statement `x == 0`. `IsZero` is already
 * defined for the type `NTL::ZZ` in NTL, but not for `std::int64_t`.
 */
inline bool IsZero(const std::int64_t &x) {
    return x == 0;
}

/**
 * Sets `x` to 0.
 */
inline void clear(double &x) {
    x = 0;
}

/**
 * Sets `x` to 0.
 */
inline void clear(std::int64_t &x) {
    x = 0;
}

/**
 * Tests if `x` is odd. Returns 1 if it is odd, and 0 if it is even.
 */
inline std::int64_t IsOdd(const std::int64_t &x) {
    return x & 1;
}

/**
 * Sets `x` to 1.
 */
inline void set(std::int64_t &x) {
    x = 1;
}

/**
 * \name Mathematical functions
 * @{ These are complementary overloads to NTL power functions.
 */

/**
 * Returns \f$p^i\f$.
 */
inline std::int64_t power(std::int64_t p, std::int64_t i) {
    return NTL::power_long(p, i);
}

/**
 * Sets \f$z = 2^i\f$.
 */
inline void power2(std::int64_t &z, std::int64_t i) {
    z = NTL::power_long(2, i);
}
/**
 * Sets \f$z = 2^i\f$.
 */
inline void power2(NTL::ZZ &z, std::int64_t i) {
    z = NTL::power_ZZ(2, i);
}

//inline double sqrt(const double &a) { return std::sqrt(a); }

//inline double log(const double x) {
//	return std::log(x);
//}

inline double inv(const double x) {
    return 1.0 / x;
}

/**
 * Transposes `A` into `X`.
 * This is a template overload of the transpose function of NTL.
 * It does basically the same thing as in NTL. It might
 * be necessary to implement the swap function for the type `T` for this to work.
 *
 * This seems a bad idea, because it has the same signature as the NTL method,
 *  and it is probably slower.  Remove?                     !!!!!      *****
 *
 */
template<typename T>
static void transpose(NTL::Mat<T> &X, const NTL::Mat<T> &A) {
    int64_t n = A.NumRows();
    int64_t m = A.NumCols();
    int64_t i, j;

    // If both matrices have the same address, we need to transpose in place
    if (&X == &A) {
        if (n == m)
            for (i = 0; i < n; i++)
                for (j = 0; j < m; j++)
                    std::swap(X[i][j], X[j][i]);
        else {
            NTL::Mat<T> tmp;
            tmp.SetDims(m, n);
            for (i = 0; i < n; i++)
                for (j = 0; j < m; j++)
                    tmp[j][i] = A[i][j];
            X.kill();
            X = tmp;
        }
    } else {
        X.SetDims(m, n);
        for (i = 0; i < n; i++)
            for (j = 0; j < m; j++)
                X[j][i] = A[i][j];
    }
}

/**
 * Another implementation of the `transpose` function.  Returns the transpose of `a`.
 * No need for this:  Just call NTL::transpose(x, a) directly!  *********
 * */
template<typename T>
static inline NTL::Mat<T> transpose(const NTL::Mat<T> &A) {
    NTL::Mat<T> X;
    transpose(X, A);
    return X;
}


//============================================================================

// Defining some aliases.

// typedef std::int64_t int64_t;

//typedef NTL::Mat<int64_t> matrix64;
//typedef NTL::Vec<int64_t> vector64;

//typedef Mat<std::int64_t> Mat_64;
//typedef Vec<std::int64_t> Vec_64;

//==========================================

/** \name Inline functions for certain operators, again for compatibility with NTL.
* @{
*/

inline static void add(int64_t &x, const int64_t a, const int64_t b) {
    x = a + b;
}

inline static void sub(int64_t &x, const int64_t a, const int64_t b) {
    x = a - b;
}

// x = a - b;  assumes a >= b >= 0.
inline static void SubPos(int64_t &x, const int64_t a, const int64_t b) {
    x = a - b;
}

inline static void negate(int64_t &x, const int64_t a) {
    x = -a;
}

//inline static void abs(int64_t& x, const int64_t& a)
//   { x = abs(a); }

inline static void mul(int64_t &x, const int64_t a, const int64_t b) {
    x = a * b;
}

// Integer division.
inline static void div(int64_t &x, const int64_t a, const int64_t b) {
    x = a / b;
}

// Modulo. Not defined for `int64_t` in NTL.
// This one always returns a non-negative value, like for ZZ in NTL.
inline static void rem(int64_t &x, const int64_t a, const int64_t b) {
    x = a % b;
    if (x < 0) x += b;
}

inline static void sqr(int64_t &x, const int64_t a) {
    x = a * a;
}

inline static void MulAddTo(int64_t &x, const int64_t a, const int64_t b) {
    x += a * b;
}

inline static void MulSubFrom(int64_t &x, int64_t a, int64_t b) {
    x -= a * b;
}

inline static void LeftShift(int64_t &x, const int64_t a, int64_t k) {
    x = (a << k);
}

inline static void RightShift(int64_t &x, const int64_t a, int64_t k) {
    x = (a >> k);
}

// inline static bool IsZero(int64_t x)
//   {  return (x == 0); }

/*
 static double InnerProduct(double *a, double *b, int64_t n) {
 register double s=0;
 for (int64_t i = 0; i < n; i++)
 s += a[i]*b[i];
 return s;
 }


 static void InnerProduct(int64_t& xx, const vector64& a, const vector64& b) {
 register int64_t x = 0;
 int64_t n = min(a.length(), b.length());
 int64_t i;
 for (i = 0; i < n; i++) {
 x += a[i] * b[i];
 }
 xx = x;
 }
 */


/*
void static mul(Vec_64 &x, const Vec_64 &a, const int64_t b) {
    int64_t n = a.length();
    x.SetLength(n);
    int64_t i;
    for (i = 0; i < n; i++)
        mul(x[i], a[i], b);
}

void static add(Vec_64 &x, const Vec_64 &a, const Vec_64 &b) {
    int64_t n = a.length();
    if (b.length() != n)
        LogicError("vector add: dimension mismatch");
    x.SetLength(n);
    int64_t i;
    for (i = 0; i < n; i++)
        add(x[i], a[i], b[i]);
}

void static sub(Vec_64 &x, const Vec_64 &a, const Vec_64 &b) {
    int64_t n = a.length();
    if (b.length() != n)
        LogicError("vector sub: dimension mismatch");
    x.SetLength(n);
    int64_t i;
    for (i = 0; i < n; i++)
        sub(x[i], a[i], b[i]);
}
*/

/**
 * These are operator overloads for Mat_64 and Vec_64 types. Only the
 * overloads we currently use are defined.
 */

/*
Vec_64 operator*(const Vec_64 &vec, std::int64_t a);
Vec_64 operator*(std::int64_t a, const Vec_64 &vec);
std::int64_t operator*(const Vec_64 &vec1, const Vec_64 &vec2);
Vec_64& operator+=(Vec_64 &vec1, const Vec_64 &vec2);
Vec_64& operator-=(Vec_64 &vec1, const Vec_64 &vec2);
Vec_64& operator*=(Vec_64 &vec, std::int64_t a);
Mat_64& operator*=(Mat_64 &mat, std::int64_t a);
Mat_64 operator*(const Mat_64 &mat1, const Mat_64 &mat2);
*/

/**
 * Transforms `mat` into the identity matrix of dimensions
 * \f$\text{dim}\times\text{dim}\f$.
 */
//void ident(Mat_64 &mat, int64_t dim);

/**
 * Computes and returns the determinant of `mat'.
 */
//double determinant(const Mat_64 &mat);

} // End namespace NTL

#endif // LATTICETESTER__NTLWRAP_H
