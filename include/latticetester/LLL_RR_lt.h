//
// This file is part of LatticeTester, although most of it is just a modified
// version of the `LLL_FP` module of NTL available at https://libntl.org/.
// It was modified because we wanted extra flexibility in the functions to improve
// the performance of our tools that use these functions.

#ifndef NTL_LLL_RR_lt__H
#define NTL_LLL_RR_lt__H

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <set>
#include <map>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <type_traits>

#include <NTL/tools.h>
#include <NTL/fileio.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/vec_double.h>
#include <NTL/ZZ.h>
#include <NTL/LLL.h>

#include <NTL/LLL_RR.cpp>

#include <latticetester/Util.h>


/**
 * This module is a slight modification of `LLL_RR.cpp` from NTL.
 * The modifications are similar to those in `LLL_FPInt`.
 * With the modified functions, we can apply LLL or BKZ to a submatrix
 * (first `r` rows and `c` columns) of the matrix `B` that is passed in and returned.
 * The returned basis will have `c` columns and at most `max(r,c)` rows
 * (the rank of the basis matrix), so it may not occupy the entire space in `B`.
 * We can also recover an array `sqlen` that gives the square
 * Euclidean lengths of the basis vectors, either in `double` or `RR`.
 * Each function returns the dimension of the computed basis (number of independent rows).
 * Important: This basis is always returned in the upper-left corner of the matrix `B`.
 * This differs from the `LLL_RR` functions, which returns the zero vectors at the top.
 */

typedef NTL::vector<Int> IntVec;
typedef NTL::matrix<Int> IntMat;

NTL_OPEN_NNS

/**
 * These two functions are wrappers of `LLL_RR` in NTL, with a modified interface
 * similar to that of `LLL_FPInt`.  When `r` is positive, only the first `r`
 * rows of the matrix `B` are considered, and when `c` is positive, only the
 * first `c` columns of `B` are considered.
 * The other elements of `B` are ignored.
 * The LLL-reduced basis is returned in the upper left corner of `B`,
 * with the shortest basis vector always in the first row.
 * If `r=0`, then all the rows of the `IntMat` object are taken.
 * If `c=0`, then all the columns are taken.
 * The square lengths of the returned basis vectors are returned in the
 * vector `sqlen` if this vector is given (nonzero).
 * The functions return the dimension of the computed basis (the number of independent rows).
 */
static long LLL_RR_lt(mat_ZZ& B, const RR& delta = 0.99, long r = 0, long c = 0,
        vec_RR* sqlen = 0);

static long LLL_RR_lt(mat_ZZ &B, const double delta = 0.99, long r = 0, long c = 0,
        Vec<double}* sqlen = 0);

/**
 * These two functions are wrappers of `BKZ_RR` in NTL, with the same modifications
 * as in `LLL_RR_lt` above.
 */
static long BKZ_RR_lt(mat_ZZ &BB, const RR& delta = 0.99, long blocksize = 10,
        long prune = 0, long r = 0, long c = 0, vec_RR* sqlen = 0);

static long BKZ_RR_lt(mat_ZZ &BB, const double delta = 0.99, long blocksize = 10,
        long prune = 0, long r = 0, long c = 0, Vec<double>* sqlen = 0);


NTL_CLOSE_NNS

/* ============================================================== */

NTL_START_IMPL

// Here, `delta` and `sqlen` are in `RR`, as in NTL.
long LLL_RR_lt(mat_ZZ& B, const RR& delta, long m, long n, vec_RR* sqlen) {
       if (m == 0) m = B.NumRows();
       if (n == 0) n = B.NumCols();
       long i, j;
       RR s;
       ZZ MU, T1
       RR mu1, t1;

       NumSwaps = 0;
       if (delta < 0.50 || delta >= 1) LogicError("LLL_RR: bad delta");
       init_red_fudge();
       mat_RR B1;  // approximates B
       B1.SetDims(m, n);
       mat_RR mu;
       mu.SetDims(m, m);
       vec_RR c;  // squared lengths of Gramm-Schmidt basis vectors
       c.SetLength(m);

       vec_RR sqlen2; // squared lengths of basis vectors
       sqlen2.SetLength(m);
       // if (!sqlen) sqlen.SetLength(m);
       // if (sqlen.length() < m) sqlen.SetLength(m);

       for (i = 0; i < m; i++)
          for (j = 0; j < n; j++)
             conv(B1[i][j], B[i][j]);
       for (i = 0; i < m; i++) {
          InnerProduct(sqlen2[i], B1[i], B1[i]);
       }
       long new_m = ll_LLL_RR(B, 0, delta, 0, 0, B1, mu, sqlen2, c, m, 1, 0);

       // In this version, we leave the zero rows at the bottom.
       // The new_m independent basis vectors will be at the top of `B`.
       // We also put the shortest nonzero vector in first place.
       long imin = 0;
       RR minSqlen = sqlen2[0];
       for (i = 1; i < new_m; i++)
           if (sqlen2[i] < minSqlen) {
               minSqlen = sqlen2[i];
               imin = i;
           };
       if (imin > 0) {
           NTL::swap(B[0], B[imin]);
           NTL::swap(sqlen2[0], sqlen2[imin]);
       }
       if (sqlen) {
           if (sqlen->length() < new_m) sqlen->SetLength(new_m);
           for (i = 0; i < new_m; i++)  (sqlen*)[i] = sqlen2[i];
       }
       return new_m;
    }

// Here, `delta` and `sqlen` are in `double`.  We need to create a `VecRR` each time
// to call the other function.
long LLL_RR_lt(mat_ZZ& B, const double delta, long m, long n, Vec<double>* sqlen) {
    vec_RR sqlenRR;
    sqlenRR.SetDims(m);
    RR Delta;
    conv(Delta, delta);
    long new_m = LLL_RR_lt(B, Delta, m, n, sqlenRR);
    if (sqlen) sqlen = conv<double> (sqlenRR);
    return new_m;
}


// This is for BKZ with RR.
long BKZ_RR_lt(mat_ZZ& BB, const RR& delta, long beta, long prune,
        long m, long n, Vec<double>* sqlen) {
   NTL_TLS_GLOBAL_ACCESS(red_fudge);
   NTL_TLS_GLOBAL_ACCESS(BKZThresh);

   if (m == 0) m = B.NumRows();
   if (n == 0) n = B.NumCols();
   long m_orig = m;
   long i, j;
   ZZ MU;
   RR t1, t2;
   ZZ T1;

   init_red_fudge();
   mat_ZZ B;
   B = BB;
   B.SetDims(m+1, n);

   mat_RR B1;
   B1.SetDims(m+1, n);
   mat_RR mu;
   mu.SetDims(m+1, m);

   vec_RR c;
   c.SetLength(m+1);

   vec_RR b;
   b.SetLength(m+1);

   RR cbar;

   vec_RR ctilda;
   ctilda.SetLength(m+1);

   vec_RR vvec;
   vvec.SetLength(m+1);

   vec_RR yvec;
   yvec.SetLength(m+1);

   vec_RR uvec;
   uvec.SetLength(m+1);

   vec_RR utildavec;
   utildavec.SetLength(m+1);

   vec_long Deltavec;
   Deltavec.SetLength(m+1);

   vec_long deltavec;
   deltavec.SetLength(m+1);

   long quit;
   long new_m;
   long z, jj, kk;
   long s, t;
   long h;

   for (i = 0; i < m; i++)
      for (j = 0; j < n; j++)
         conv(B1[i][j], B[i][j]);
   for (i = 0; i < m; i++) {
      InnerProduct(b[i], B1[i], B1[i]);
   }
   m = ll_LLL_RR(B, 0, delta, 0, 0, B1, mu, b, c, m, 1, quit);

   double tt;
   double enum_time = 0;
   unsigned long NumIterations = 0;
   unsigned long NumTrivial = 0;
   unsigned long NumNonTrivial = 0;
   unsigned long NumNoOps = 0;

   if (m < m_orig) {
      for (i = m_orig; i >= m+1; i--) {
         swap(B[i], B[i-1]);
      }
   }
   long clean = 1;
   if (m > 1) {
      if (beta > m) beta = m;
      if (prune > 0)
         ComputeBKZConstant(beta, prune);
      z = 0;
      jj = 0;
      while (z < m-1) {
         jj++;
         kk = min(jj+beta-1, m);
         if (jj == m) {
            jj = 1;
            kk = beta;
            clean = 1;
         }
         if (prune > 0)
            ComputeBKZThresh(&c(jj), kk-jj+1);
         cbar = c(jj);
         conv(utildavec(jj), 1);
         conv(uvec(jj), 1);
         conv(yvec(jj), 0);
         conv(vvec(jj), 0);
         Deltavec(jj) = 0;
         s = t = jj;
         deltavec(jj) = 1;

         for (i = jj+1; i <= kk+1; i++) {
            conv(ctilda(i), 0);
            conv(uvec(i), 0);
            conv(utildavec(i), 0);
            conv(yvec(i), 0);
            Deltavec(i) = 0;
            conv(vvec(i), 0);
            deltavec(i) = 1;
         }
         long enum_cnt = 0;
         while (t <= kk) {
            add(t1, yvec(t), utildavec(t));
            sqr(t1, t1);
            mul(t1, t1, c(t));
            add(ctilda(t), ctilda(t+1), t1);

            if (prune > 0 && t > jj)
               sub(t1, cbar, BKZThresh(t-jj));
            else
               t1 = cbar;
            if (ctilda(t) <t1) {
               if (t > jj) {
                  t--;
                  clear(t1);
                  for (i = t+1; i <= s; i++) {
                     mul(t2, utildavec(i), mu(i,t));
                     add(t1, t1, t2);
                  }
                  yvec(t) = t1;
                  negate(t1, t1);
                  if (sign(t1) >= 0) {
                     sub(t1, t1, 0.5);
                     ceil(t1, t1);
                  }
                  else {
                     add(t1, t1, 0.5);
                     floor(t1, t1);
                  }
                  utildavec(t) = t1;
                  vvec(t) = t1;
                  Deltavec(t) = 0;
                  negate(t1, t1);

                  if (t1 < yvec(t))
                     deltavec(t) = -1;
                  else
                     deltavec(t) = 1;
               }
               else {
                  cbar = ctilda(jj);
                  for (i = jj; i <= kk; i++) {
                     uvec(i) = utildavec(i);
                  }
               }
            }
            else {
               t++;
               s = max(s, t);
               if (t < s) Deltavec(t) = -Deltavec(t);
               if (Deltavec(t)*deltavec(t) >= 0) Deltavec(t) += deltavec(t);
               add(utildavec(t), vvec(t), Deltavec(t));
            }
         }
         NumIterations++;
         h = min(kk+1, m);
         mul(t1, red_fudge, -8);
         add(t1, t1, delta);
         mul(t1, t1, c(jj));
         if (t1 > cbar) {
            clean = 0;
            // we treat the case that the new vector is b_s (jj < s <= kk)
            // as a special case that appears to occur most of the time.
            s = 0;
            for (i = jj+1; i <= kk; i++) {
               if (uvec(i) != 0) {
                  if (s == 0)
                     s = i;
                  else
                     s = -1;
               }
            }
            if (s == 0) LogicError("BKZ_RR: internal error, s==0");
            if (s > 0) {
               // special case
               // cerr << "special case\n";
               NumTrivial++;
               for (i = s; i > jj; i--) {
                  // swap i, i-1
                  swap(B(i-1), B(i));
                  swap(B1(i-1), B1(i));
                  swap(b(i-1), b(i));
               }
               new_m = ll_LLL_RR(B, 0, delta, 0, 0,
                                 B1, mu, b, c, h, jj, quit);
               if (new_m != h) LogicError("BKZ_RR: internal error");
               if (quit) break;
            }
            else {
               // the general case
               NumNonTrivial++;
               for (i = 1; i <= n; i++) conv(B(m+1, i), 0);
               for (i = jj; i <= kk; i++) {
                  if (uvec(i) == 0) continue;
                  conv(MU, uvec(i));
                  RowTransform2(B(m+1), B(i), MU);
               }
               for (i = m+1; i >= jj+1; i--) {
                  // swap i, i-1
                  swap(B(i-1), B(i));
                  swap(B1(i-1), B1(i));
                  swap(b(i-1), b(i));
               }
               for (i = 1; i <= n; i++)
                  conv(B1(jj, i), B(jj, i));
               InnerProduct(b(jj), B1(jj), B1(jj));
               if (b(jj) == 0) LogicError("BKZ_RR: internal error, b(jj) == 0");

               // remove linear dependencies
               // cerr << "general case\n";
               new_m = ll_LLL_RR(B, 0, delta, 0, 0, B1, mu, b, c, kk+1, jj, quit);
               if (new_m != kk) LogicError("BKZ_RR: internal error");

               // remove zero vector
               for (i = kk+2; i <= m+1; i++) {
                  // swap i, i-1
                  swap(B(i-1), B(i));
                  swap(B1(i-1), B1(i));
                  swap(b(i-1), b(i));
               }
               quit = 0;
               if (check) {
                  for (i = 1; i <= kk; i++)
                     if ((*check)(B(i))) {
                        quit = 1;
                        break;
                     }
               }
               if (quit) break;
               if (h > kk) {
                  // extend reduced basis
                  new_m = ll_LLL_RR(B, 0, delta, 0, 0,
                                   B1, mu, b, c, h, h, quit);
                  if (new_m != h) LogicError("BKZ_RR: internal error");
                  if (quit) break;
               }
            }
            z = 0;
         }
         else {
            NumNoOps++;
            if (!clean) {
               new_m =
                  ll_LLL_RR(B, 0, delta, 0, 0, B1, mu, b, c, h, h, quit);
               if (new_m != h) LogicError("BKZ_RR: internal error");
               if (quit) break;
            }
            z++;
         }
      }
   }
// In this version, we do not move the zero vectors to the top.
// We also do not change the dimensions of BB.
    for (i = 0; i < m_orig; i++) {
        for (j = 0; j < n; j++) {
            BB[i][j] = B[i][j];
        }
    }
// Put the shortest nonzero vector in first place.
    long imin = 0;
    RR minlen = b[0];
    for (i = 1; i < m; i++)
        if (b[i] < minlen) {
            minlen = b[i];
            imin = i;
        };
    if (imin > 0) {
        swap(BB[0], BB[imin]);
        swap(b[1], b[imin]);
    }
    if (sqlen)
        for (i = 0; i < m; i++)
            sqlen[i] = b[i];
    return m;    // Number of rows in basis.
}

// Here, `delta` and `sqlen` are in `double`.  We need to create a `VecRR` each time
// to call the other function.
long BKZ_RR_lt(mat_ZZ& BB, const double delta, long beta, long prune,
         long m, long n, Vec<double>* sqlen) {
    vec_RR sqlenRR;
    sqlenRR.SetDims(m);
    RR Delta;
    conv(Delta, delta);
    long new_m = BKZ_RR_lt(BB, Delta, beta, prune, m, n, sqlenRR);
    if (sqlen) sqlen = conv<double> (sqlenRR);
    return new_m;
}


NTL_END_IMPL

#endif
