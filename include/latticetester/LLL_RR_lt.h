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
#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>

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

static long LLL_RR_lt(mat_ZZ& B, const RR& delta = conv<RR>(0.99999), long r = 0, long c = 0,
        vec_RR* sqlen = 0);

static long LLL_RR_lt(mat_ZZ &B, const double delta = 0.99999, long r = 0, long c = 0,
        NTL::Vec<double>* sqlen = NULL);

/**
 * These two functions are wrappers of `BKZ_RR` in NTL, with the same modifications
 * as in `LLL_RR_lt` above.
 */
static long BKZ_RR_lt(mat_ZZ &BB, const RR& delta = conv<RR>(0.99999), long blocksize = 10,
        long prune = 0, long r = 0, long c = 0, NTL::vec_RR* sqlen = 0);

static long BKZ_RR_lt(mat_ZZ &BB, const double delta = 0.99999, long blocksize = 10,
        long prune = 0, long r = 0, long c = 0, NTL::Vec<double>* sqlen = 0);


NTL_CLOSE_NNS

/* ============================================================== */

NTL_START_IMPL

NTL_TLS_GLOBAL_DECL(RR, red_fudge_RR)

static NTL_CHEAP_THREAD_LOCAL long log_red_RR = 0;

static void init_red_fudge_RR()
{
   NTL_TLS_GLOBAL_ACCESS(red_fudge_RR);
   log_red_RR = long(0.50*RR::precision());
   power2(red_fudge_RR, -log_red_RR);
}

static void inc_red_fudge_RR()
{
   NTL_TLS_GLOBAL_ACCESS(red_fudge_RR);
   mul(red_fudge_RR, red_fudge_RR, 2);
   log_red_RR--;
   cerr << "LLL_RR: warning--relaxing reduction (" << log_red_RR << ")\n";
   if (log_red_RR < 4)
      ResourceError("LLL_RR: can not continue, log_red_RR < 4, ...sorry");
}

long ll_LLL_RR_lt(mat_ZZ& B, mat_ZZ* U, const RR& delta, long deep,
           LLLCheckFct check, mat_RR& B1, mat_RR& mu,
           vec_RR& b, vec_RR& c, long m, long init_k, long &quit);

// Here, `delta` and `sqlen` are in `RR`, as in NTL.
long LLL_RR_lt(NTL::mat_ZZ& B, const RR& delta, long m, long n, vec_RR* sqlen) {
       if (m == 0) m = B.NumRows();
       if (n == 0) n = B.NumCols();
       long i, j, new_m, quit;
       RR s;
       ZZ MU, T1;
       RR mu1, t1;

       NumSwaps = 0;
       if (delta < 0.50 || delta >= 1) LogicError("LLL_RR: bad delta");
       init_red_fudge_RR();
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
       new_m = ll_LLL_RR_lt(B, 0, delta, 0, 0, B1, mu, sqlen2, c, m, 1, quit);
       // new_m = ll_LLL_RR(B, U, delta, deep, check, B1, mu, b, c, m, 1, quit);

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
           for (i = 0; i < new_m; i++)  (*sqlen)[i] = sqlen2[i];
       }
       return new_m;
    }

// Here, `delta` and `sqlen` are in `double`.  We need to create a `VecRR` each time
// to call the other function.
long LLL_RR_lt(mat_ZZ& B, const double delta, long m, long n, NTL::Vec<double>* sqlen) {
    vec_RR* sqlenRR;
    if (sqlen) (*sqlenRR).SetLength(m);
    // RR Delta;
    // conv(Delta, delta);
    long new_m = LLL_RR_lt(B, conv<RR>(delta), m, n, sqlenRR);
    if (sqlen) {
        if (sqlen->length() < new_m) sqlen->SetLength(new_m);
        for (int i = 0; i < new_m; i++)  (*sqlen)[i] = conv<double> ((*sqlenRR)[i]);
    }
    return new_m;
}



NTL_TLS_GLOBAL_DECL(vec_RR, BKZConstant)
NTL_TLS_GLOBAL_DECL(vec_RR, BKZThresh)

// This is for BKZ with RR.
long BKZ_RR_lt(mat_ZZ& BB, const RR& delta, long beta, long prune,
        long m, long n, NTL::Vec<double>* sqlen) {

   NTL_TLS_GLOBAL_ACCESS(red_fudge_RR);
   NTL_TLS_GLOBAL_ACCESS(BKZThresh);
   if (m == 0) m = BB.NumRows();
   if (n == 0) n = BB.NumCols();
   long m_orig = m;
   long i, j;
   ZZ MU;
   RR t1, t2;
   ZZ T1;

   init_red_fudge_RR();
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
   m = ll_LLL_RR_lt(B, 0, delta, 0, 0, B1, mu, b, c, m, 1, quit);

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
         mul(t1, red_fudge_RR, -8);
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
               new_m = ll_LLL_RR_lt(B, 0, delta, 0, 0,
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
               new_m = ll_LLL_RR_lt(B, 0, delta, 0, 0, B1, mu, b, c, kk+1, jj, quit);
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
                  new_m = ll_LLL_RR_lt(B, 0, delta, 0, 0,
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
                  ll_LLL_RR_lt(B, 0, delta, 0, 0, B1, mu, b, c, h, h, quit);
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
         long m, long n, NTL::Vec<double>& sqlen) {
    vec_RR sqlenRR;
    if (sqlen) sqlenRR.SetDims(m);
    // RR Delta;
    // conv(Delta, delta);
    long new_m = BKZ_RR_lt(BB, conv<RR>(delta), beta, prune, m, n, sqlenRR);
    if (sqlen) {
        if (sqlen.length() < new_m) sqlen.SetLength(new_m);
        for (int i = 0; i < new_m; i++)  sqlen[i] = conv<double> (sqlenRR[i]);
    }
    return new_m;
}

// ===========================================================================

// Copied from NTL::LLL_RR

// NTL_START_IMPL


static void RowTransform(vec_ZZ& A, vec_ZZ& B, const ZZ& MU1)
// x = x - y*MU
{
   NTL_ZZRegister(T);
   NTL_ZZRegister(MU);
   long k;

   long n = A.length();
   long i;

   MU = MU1;

   if (MU == 1) {
      for (i = 1; i <= n; i++)
         sub(A(i), A(i), B(i));

      return;
   }

   if (MU == -1) {
      for (i = 1; i <= n; i++)
         add(A(i), A(i), B(i));

      return;
   }

   if (MU == 0) return;

   if (NumTwos(MU) >= NTL_ZZ_NBITS)
      k = MakeOdd(MU);
   else
      k = 0;


   if (MU.WideSinglePrecision()) {
      long mu1;
      conv(mu1, MU);

      for (i = 1; i <= n; i++) {
         mul(T, B(i), mu1);
         if (k > 0) LeftShift(T, T, k);
         sub(A(i), A(i), T);
      }
   }
   else {
      for (i = 1; i <= n; i++) {
         mul(T, B(i), MU);
         if (k > 0) LeftShift(T, T, k);
         sub(A(i), A(i), T);
      }
   }
}

static void RowTransform2(vec_ZZ& A, vec_ZZ& B, const ZZ& MU1)
// x = x + y*MU
{
   NTL_ZZRegister(T);
   NTL_ZZRegister(MU);
   long k;

   long n = A.length();
   long i;

   MU = MU1;

   if (MU == 1) {
      for (i = 1; i <= n; i++)
         add(A(i), A(i), B(i));

      return;
   }

   if (MU == -1) {
      for (i = 1; i <= n; i++)
         sub(A(i), A(i), B(i));

      return;
   }

   if (MU == 0) return;

   if (NumTwos(MU) >= NTL_ZZ_NBITS)
      k = MakeOdd(MU);
   else
      k = 0;

   if (MU.WideSinglePrecision()) {
      long mu1;
      conv(mu1, MU);

      for (i = 1; i <= n; i++) {
         mul(T, B(i), mu1);
         if (k > 0) LeftShift(T, T, k);
         add(A(i), A(i), T);
      }
   }
   else {
      for (i = 1; i <= n; i++) {
         mul(T, B(i), MU);
         if (k > 0) LeftShift(T, T, k);
         add(A(i), A(i), T);
      }
   }
}

void ComputeGS(const mat_ZZ& B, mat_RR& B1,
               mat_RR& mu, vec_RR& b,
               vec_RR& c, long k, const RR& bound, long st,
               vec_RR& buf, const RR& bound2)
{
   long i, j;
   RR s, t, t1;
   ZZ T1;

   if (st < k) {
      for (i = 1; i < st; i++)
         mul(buf(i), mu(k,i), c(i));
   }

   for (j = st; j <= k-1; j++) {
      InnerProduct(s, B1(k), B1(j));

      sqr(t1, s);
      mul(t1, t1, bound);
      mul(t, b(k), b(j));

      if (t >= bound2 && t >= t1) {
         InnerProduct(T1, B(k), B(j));
         conv(s, T1);
      }

      clear(t1);
      for (i = 1; i <= j-1; i++) {
         mul(t, mu(j, i), buf(i));
         add(t1, t1, t);
      }

      sub(t, s, t1);
      buf(j) = t;
      div(mu(k, j), t, c(j));
   }


   clear(s);
   for (j = 1; j <= k-1; j++) {
      mul(t, mu(k, j), buf(j));
      add(s, s, t);
   }

   sub(c(k), b(k), s);
}


static NTL_CHEAP_THREAD_LOCAL long verbose = 0;
static NTL_CHEAP_THREAD_LOCAL unsigned long NumSwaps = 0;
static NTL_CHEAP_THREAD_LOCAL double StartTime = 0;
static NTL_CHEAP_THREAD_LOCAL double LastTime = 0;



static void LLLStatus(long max_k, double t, long m, const mat_ZZ& B)
{
   cerr << "---- LLL_RR status ----\n";
   cerr << "elapsed time: ";
   PrintTime(cerr, t-StartTime);
   cerr << ", stage: " << max_k;
   cerr << ", rank: " << m;
   cerr << ", swaps: " << NumSwaps << "\n";

   ZZ t1;
   long i;
   double prodlen = 0;

   for (i = 1; i <= m; i++) {
      InnerProduct(t1, B(i), B(i));
      if (!IsZero(t1))
         prodlen += log(t1);
   }

   cerr << "log of prod of lengths: " << prodlen/(2.0*log(2.0)) << "\n";

   if (LLLDumpFile64) {
      cerr << "dumping to " << LLLDumpFile64 << "...";

      ofstream f;
      OpenWrite(f, LLLDumpFile64);

      f << "[";
      for (i = 1; i <= m; i++) {
         f << B(i) << "\n";
      }
      f << "]\n";

      f.close();

      cerr << "\n";
   }
   LastTime = t;
}


static
long ll_LLL_RR_lt(mat_ZZ& B, mat_ZZ* U, const RR& delta, long deep,
           LLLCheckFct check, mat_RR& B1, mat_RR& mu,
           vec_RR& b, vec_RR& c, long m, long init_k, long &quit) {
   NTL_TLS_GLOBAL_ACCESS(red_fudge_RR);

   long n = B.NumCols();

   long i, j, k, Fc1;
   ZZ MU;
   RR mu1, t1, t2, cc;
   ZZ T1;

   RR bound;

      // we tolerate a 15% loss of precision in computing
      // inner products in ComputeGS.

   power2(bound, 2*long(0.15*RR::precision()));


   RR bound2;

   power2(bound2, 2*RR::precision());


   quit = 0;
   k = init_k;

   vec_long st_mem;
   st_mem.SetLength(m+2);
   long *st = st_mem.elts();

   for (i = 1; i < k; i++)
      st[i] = i;

   for (i = k; i <= m+1; i++)
      st[i] = 1;

   vec_RR buf;
   buf.SetLength(m);

   long rst;
   long counter;

   long trigger_index;
   long small_trigger;
   long cnt;

   RR half;
   conv(half,  0.5);
   RR half_plus_fudge;
   add(half_plus_fudge, half, red_fudge_RR);

   long max_k = 0;
   double tt;

   while (k <= m) {

      if (k > max_k) {
         max_k = k;
      }

      if (verbose) {
         tt = GetTime();

         if (tt > LastTime + LLLStatusInterval64)
            LLLStatus(max_k, tt, m, B);
      }


      if (st[k] == k)
         rst = 1;
      else
         rst = k;

      if (st[k] < st[k+1]) st[k+1] = st[k];
      ComputeGS(B, B1, mu, b, c, k, bound, st[k], buf, bound2);
      st[k] = k;

      counter = 0;
      trigger_index = k;
      small_trigger = 0;
      cnt = 0;

      do {
         // size reduction

         counter++;
         if (counter > 10000) {
            cerr << "LLL_XD: warning--possible infinite loop\n";
            counter = 0;
         }


         Fc1 = 0;

         for (j = rst-1; j >= 1; j--) {
            abs(t1, mu(k,j));
            if (t1 > half_plus_fudge) {

               if (!Fc1) {
                  if (j > trigger_index ||
                      (j == trigger_index && small_trigger)) {

                     cnt++;

                     if (cnt > 10) {
                        inc_red_fudge_RR();
                        add(half_plus_fudge, half, red_fudge_RR);
                        cnt = 0;
                     }
                  }

                  trigger_index = j;
                  small_trigger = (t1 < 4);
               }

               Fc1 = 1;

               mu1 = mu(k,j);
               if (sign(mu1) >= 0) {
                  sub(mu1, mu1, half);
                  ceil(mu1, mu1);
               }
               else {
                  add(mu1, mu1, half);
                  floor(mu1, mu1);
               }

               if (mu1 == 1) {
                  for (i = 1; i <= j-1; i++)
                     sub(mu(k,i), mu(k,i), mu(j,i));
               }
               else if (mu1 == -1) {
                  for (i = 1; i <= j-1; i++)
                     add(mu(k,i), mu(k,i), mu(j,i));
               }
               else {
                  for (i = 1; i <= j-1; i++) {
                     mul(t2, mu1, mu(j,i));
                     sub(mu(k,i), mu(k,i), t2);
                  }
               }


               conv(MU, mu1);

               sub(mu(k,j), mu(k,j), mu1);

               RowTransform(B(k), B(j), MU);
               if (U) RowTransform((*U)(k), (*U)(j), MU);
            }
         }

         if (Fc1) {
            for (i = 1; i <= n; i++)
               conv(B1(k, i), B(k, i));

            InnerProduct(b(k), B1(k), B1(k));
            ComputeGS(B, B1, mu, b, c, k, bound, 1, buf, bound2);
         }
      } while (Fc1);

      if (check && (*check)(B(k)))
         quit = 1;

      if (IsZero(b(k))) {
         for (i = k; i < m; i++) {
            // swap i, i+1
            swap(B(i), B(i+1));
            swap(B1(i), B1(i+1));
            swap(b(i), b(i+1));
            if (U) swap((*U)(i), (*U)(i+1));
         }

         for (i = k; i <= m+1; i++) st[i] = 1;

         m--;
         if (quit) break;
         continue;
      }

      if (quit) break;

      if (deep > 0) {
         // deep insertions

         cc = b(k);
         long l = 1;
         while (l <= k-1) {
            mul(t1, delta, c(l));
            if (t1 > cc) break;
            sqr(t1, mu(k,l));
            mul(t1, t1, c(l));
            sub(cc, cc, t1);
            l++;
         }

         if (l <= k-1 && (l <= deep || k-l <= deep)) {
            // deep insertion at position l

            for (i = k; i > l; i--) {
               // swap rows i, i-1
               swap(B(i), B(i-1));
               swap(B1(i), B1(i-1));
               swap(mu(i), mu(i-1));
               swap(b(i), b(i-1));
               if (U) swap((*U)(i), (*U)(i-1));
            }

            k = l;
            continue;
         }
      } // end deep insertions

      // test LLL reduction condition

      if (k <= 1) {
         k++;
      }
      else {
         sqr(t1, mu(k,k-1));
         mul(t1, t1, c(k-1));
         add(t1, t1, c(k));
         mul(t2, delta, c(k-1));
         if (t2 > t1) {
            // swap rows k, k-1
            swap(B(k), B(k-1));
            swap(B1(k), B1(k-1));
            swap(mu(k), mu(k-1));
            swap(b(k), b(k-1));
            if (U) swap((*U)(k), (*U)(k-1));

            k--;
            NumSwaps++;
         }
         else {
            k++;
         }
      }
   }

   if (verbose) {
      LLLStatus(m+1, GetTime(), m, B);
   }


   return m;
}


static
void ComputeBKZConstant(long beta, long p)
{
   NTL_TLS_GLOBAL_ACCESS(BKZConstant);

   RR c_PI;
   ComputePi(c_PI);

   RR LogPI = log(c_PI);

   BKZConstant.SetLength(beta-1);

   vec_RR Log;
   Log.SetLength(beta);


   long i, j, k;
   RR x, y;

   for (j = 1; j <= beta; j++)
      Log(j) = log(to_RR(j));

   for (i = 1; i <= beta-1; i++) {
      // First, we compute x = gamma(i/2)^{2/i}

      k = i/2;

      if ((i & 1) == 0) { // i even
         x = 0;
         for (j = 1; j <= k; j++)
            x += Log(j);

         x = exp(x/k);

      }
      else { // i odd
         x = 0;
         for (j = k + 2; j <= 2*k + 2; j++)
            x += Log(j);

         x += 0.5*LogPI - 2*(k+1)*Log(2);

         x = exp(2*x/i);
      }

      // Second, we compute y = 2^{2*p/i}

      y = exp(-(2*p/to_RR(i))*Log(2));

      BKZConstant(i) = x*y/c_PI;
   }

}

static
void ComputeBKZThresh(RR *c, long beta)
{
   NTL_TLS_GLOBAL_ACCESS(BKZConstant);
   NTL_TLS_GLOBAL_ACCESS(BKZThresh);

   BKZThresh.SetLength(beta-1);

   long i;
   RR x;
   RR t1;

   x = 0;

   for (i = 1; i <= beta-1; i++) {
      log(t1, c[i-1]);
      add(x, x, t1);
      div(t1, x, i);
      exp(t1, t1);
      mul(BKZThresh(i), t1, BKZConstant(i));
   }
}


static
void BKZStatus(double tt, double enum_time, unsigned long NumIterations,
               unsigned long NumTrivial, unsigned long NumNonTrivial,
               unsigned long NumNoOps, long m,
               const mat_ZZ& B)
{
   cerr << "---- BKZ_RR status ----\n";
   cerr << "elapsed time: ";
   PrintTime(cerr, tt-StartTime);
   cerr << ", enum time: ";
   PrintTime(cerr, enum_time);
   cerr << ", iter: " << NumIterations << "\n";
   cerr << "triv: " << NumTrivial;
   cerr << ", nontriv: " << NumNonTrivial;
   cerr << ", no ops: " << NumNoOps;
   cerr << ", rank: " << m;
   cerr << ", swaps: " << NumSwaps << "\n";
   ZZ t1;
   long i;
   double prodlen = 0;

   for (i = 1; i <= m; i++) {
      InnerProduct(t1, B(i), B(i));
      if (!IsZero(t1))
         prodlen += log(t1);
   }
   cerr << "log of prod of lengths: " << prodlen/(2.0*log(2.0)) << "\n";

   if (LLLDumpFile64) {
      cerr << "dumping to " << LLLDumpFile64 << "...";

      ofstream f;
      OpenWrite(f, LLLDumpFile64);

      f << "[";
      for (i = 1; i <= m; i++) {
         f << B(i) << "\n";
      }
      f << "]\n";

      f.close();
      cerr << "\n";
   }
   LastTime = tt;
}


NTL_END_IMPL

#endif
