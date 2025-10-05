
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
#include <NTL/LLL.h>

#include <latticetester/Util.h>
#include <latticetester/LLL_FP64.h>

/* ============================================================== */

// This implementation is modified from NTL.
// Some array indices start at 1 in NTL and at 0 here, but not all of them.
// Both here and in NTL, some indices start at 0 and others start at 1.
// This makes the code complicated and not so easy to modify.

// The macro below is defined in NTL/tools.h
NTL_START_IMPL

// Just to be safe!!
#define TR_BND (NTL_FDOUBLE_PRECISION/2.0)

static inline void CheckFinite(double *p) {
   if (!IsFinite(p))
      ResourceError("LLL_FP64: numbers too big...use LLL_XD \n");
}

// Returns the largest absolute value of the coordinates of v.
static double max_abs(double *v, long n) {
   double res, t;
   res = 0.0;
   for (long i = 0; i < n; i++) {
      t = fabs(v[i]);
      if (t > res)
         res = t;
   }
   return res;
}

// Returns the inner product of two arrays of double of size n.
static double InnerProductD(double *a, double *b, long n) {
   double s = 0.0;
   for (long i = 0; i < n; i++)
      s += a[i] * b[i];
   return s;
}

void InnerProductV(int64_t &prod, const NTL::Vec<int64_t> &a,
      const NTL::Vec<int64_t> &b, long n) {
   long x = 0;
   for (long i = 0; i < n; i++) {
      x += a[i] * b[i];
   }
   prod = x;
}

// A version for RR.
static void InnerProductR(RR &xx, const vec_RR &a, const vec_RR &b, long n) {
   RR t1, x;
   clear(x);
   for (long i = 0; i < n; i++) {
      mul(t1, a[i], b[i]);
      add(x, x, t1);
   }
   xx = x;
}

// ----------------------------------------------------------
// Returns 1 in `in_float` iff all the coefficients `a[i]` are within the bounds.
static void RowTransformStart64(double *a, long *in_a, long &in_float, long n) {
   long inf = 1;
   for (long i = 0; i < n; i++) {
      in_a[i] = (a[i] < TR_BND && a[i] > -TR_BND);
      inf = inf & in_a[i];
   }
   in_float = inf;
}

static void RowTransformFinish64(NTL::Vec<int64_t> &A, double *a, long *in_a, long n) {
   for (long i = 0; i < n; i++) {
      if (in_a[i]) {
         conv(A[i], a[i]);
      } else {
         conv(a[i], A[i]);
         CheckFinite(&a[i]);
      }
   }
}

// ---------------------------------------------------------
// A = A - B*MU1  for the first n vector entries only.
// The change is on the vector A.
void RowTransform64(NTL::Vec<long> &A, NTL::Vec<long> &B, const int64_t &MU1,
      long n, double *a, double *b, long *in_a, double &max_a, double max_b,
      int64_t &in_float) {
   int64_t T, MU;
   // int64_t k;
   double mu;
   conv(mu, MU1);
   CheckFinite(&mu);
   int64_t i;
   if (in_float) {
      double mu_abs = fabs(mu);
      if (mu_abs > 0 && max_b > 0 && (mu_abs >= TR_BND || max_b >= TR_BND)) {
         in_float = 0;
      } else {
         max_a += mu_abs * max_b;
         if (max_a >= TR_BND)
            in_float = 0;
      }
   }
   if (in_float) {
      if (mu == 1) {
         for (i = 0; i < n; i++)
            a[i] -= b[i];
         return;
      }
      if (mu == -1) {
         for (i = 0; i < n; i++)
            a[i] += b[i];
         return;
      }
      if (mu == 0)
         return;
      for (i = 0; i < n; i++)
         a[i] -= mu * b[i];
      return;
   }
   // std::cout << "RowTransform64, not in_float! \n";
   MU = MU1;
   if (MU == 1) {
      for (i = 0; i < n; i++) {
         if (in_a[i] && a[i] < TR_BND && a[i] > -TR_BND && b[i] < TR_BND
               && b[i] > -TR_BND) {
            a[i] -= b[i];
         } else {
            if (in_a[i]) {
               conv(A[i], a[i]);
               in_a[i] = 0;
            }
            sub(A[i], A[i], B[i]);
         }
      }
      return;
   }
   if (MU == -1) {
      for (i = 0; i < n; i++) {
         if (in_a[i] && a[i] < TR_BND && a[i] > -TR_BND && b[i] < TR_BND
               && b[i] > -TR_BND) {
            a[i] += b[i];
         } else {
            if (in_a[i]) {
               conv(A[i], a[i]);
               in_a[i] = 0;
            }
            add(A[i], A[i], B[i]);
         }
      }
      return;
   }
   if (MU == 0)
      return;

   double b_bnd = fabs(TR_BND / mu) - 1;
   if (b_bnd < 0)
      b_bnd = 0;
   for (i = 0; i < n; i++) {
      if (in_a[i]) {
         conv(A[i], a[i]);
         in_a[i] = 0;
      }
      // A[i] = A[i] - B[i] * MU;
      mul(T, B[i], MU);
      sub(A[i], A[i], T);
      //if ((A[i] > modulus64) || (A[i] < -modulus64))
      //   std::cout << "RowTransform64: i = " << i << ",  A[i] = " << A[i] << "\n";
   }
}

// -------------------------------------------------------
// A = A + B*MU  for the first n vector entries only.
// The change is on the vector A.  This is used once in BKZ.
void RowTransformAdd(NTL::Vec<int64_t> &A, NTL::Vec<int64_t> &B,
      const int64_t &MU1, long n) {
   int64_t T, MU = MU1;
   int64_t i;
   if (MU == 1) {
      for (i = 0; i < n; i++)
         add(A[i], A[i], B[i]);
      return;
   }
   if (MU == -1) {
      for (i = 0; i < n; i++)
         sub(A[i], A[i], B[i]);
      return;
   }
   if (MU == 0)
      return;
   for (i = 0; i < n; i++) {
      T = MU * B[i];
      A[i] += T;
   }
}

// ----------------------------------------------------
// This function computes mu[k] and c[k].
// This version works with arrays of `double`.
static void ComputeGS64(NTL::Mat<long> &B, double **B1, double **mu, double *b,
      double *c, long k, long n, double bound, long st, double *buf) {
   // The indices in B, B1, mu, b, c, buf, all start at 0.
   // The integers k and st are both reduced by 1 compared with NTL version.
   long i, j;
   double s, t1, y, t;
   long T1;
   long test;
   double *mu_k = mu[k];

   // std::cout << "ComputeGS64, st = " << st << " k = " << k << "\n";
   if (st < k) {
      for (i = 0; i < st; i++)
         buf[i] = mu_k[i] * c[i];
   }
   for (j = st; j < k; j++) {
      s = InnerProductD(B1[k], B1[j], n);  // Returns a double.
      // std::cout << "ComputeGS64, j = " << j << " Inner product s = " << s << "\n";

      // test = b[k]*b[j] >= NTL_FDOUBLE_PRECISION^2
      test = (b[k] / NTL_FDOUBLE_PRECISION >= NTL_FDOUBLE_PRECISION / b[j]);
      // test = test && s^2 <= b[k]*b[j]/bound,
      // but we compute it in a strange way to avoid overflow
      if (test && (y = fabs(s)) != 0) {
         t = y / b[j];
         t1 = b[k] / y;
         // std::cout << "ComputeGS, t1 = " << t1 << "\n";
         if (t <= 1)
            test = (t * bound <= t1);
         else if (t1 >= 1)
            test = (t <= t1 / bound);
         else
            test = 0;
      }
      if (test) {
         InnerProductV(T1, B[k], B[j], n);
         conv(s, T1);
         // std::cout << "ComputeGS, T1 = s = " << s << "\n";
      }
      double *mu_j = mu[j];
      t1 = 0;
      for (i = 0; i < j; i++) {
         t1 += mu_j[i] * buf[i];
      }
      mu_k[j] = (buf[j] = (s - t1)) / c[j];
      // std::cout << "ComputeGS, mu_k[j] = " << mu_k[j] << "\n";
   }

#if (!NTL_EXT_DOUBLE)
   // Kahan summation
   double c1;
   s = c1 = 0;
   for (j = 0; j < k; j++) {
      y = mu_k[j] * buf[j] - c1;
      t = s + y;
      c1 = t - s;
      c1 = c1 - y;
      s = t;
   }
#else
   s = 0;
   for (j = 0; j < k; j++)
      s += mu_k[j] * buf[j];
#endif
   c[k] = b[k] - s;
}

static NTL_CHEAP_THREAD_LOCAL double red_fudge = 0;
static NTL_CHEAP_THREAD_LOCAL long log_red = 0;
static NTL_CHEAP_THREAD_LOCAL unsigned long NumSwaps = 0;
//static NTL_CHEAP_THREAD_LOCAL double StartTime = 0;
//static NTL_CHEAP_THREAD_LOCAL double LastTime = 0;

static void init_red_fudge() {
   long i;
   log_red = long(0.50 * NTL_DOUBLE_PRECISION);
   red_fudge = 1;
   for (i = log_red; i > 0; i--)
      red_fudge = red_fudge * 0.5;
}

static void inc_red_fudge() {
   red_fudge = red_fudge * 2;
   log_red--;
   cerr << "LLL_FP: inc_red_fudge warning--relaxing reduction (" << log_red
         << ")\n";
   if (log_red < 4)
      ResourceError("LLL_FP: too much loss of precision...stop!");
}

void printBB1(NTL::Mat<long> &B, double **B1, long m) {
   std::cout << " Matrix B = \n" << B << "\n";
   for (int i = 0; i < m; i++) {
      std::cout << " B1[" << i << "] = [" << B1[i][0] << ", " << B1[i][1]
            << ", " << B1[i][2];
      std::cout << ", " << B1[i][3] << ", " << B1[i][4] << ", " << B1[i][5]
            << " ...\n";
   }
}


// The main LLL procedure.
// The int64_t version.
// In NTL, the indices of B start at 0, but those of B1, b, c, st, start at 1.
// Here they all start at 0.
// Also init_k and k start at 1 in NTL and at 0 (they are one less) here.
long ll_LLL_FP64 (Mat<int64_t> &B, double delta, double **B1, double **mu,
      double *b, double *c, int64_t m, int64_t n, int64_t init_k, long quit) {
   int64_t i, j, k, Fc1;
   int64_t MU;
   double mu1;
   double t1;
   double *tp;
   // we tolerate a 15% loss of precision in computing
   // inner products in ComputeGS64.
   static double bound = 1;
   for (i = 2 * int64_t(0.15 * NTL_DOUBLE_PRECISION); i > 0; i--)
      bound = bound * 2;
   double half_plus_fudge = 0.5 + red_fudge;
   quit = 0;
   k = init_k;

   NTL::Vec<int64_t> st_mem;
   st_mem.SetLength(m + 2);
   int64_t *st = st_mem.elts();    // An array of integers.

   for (i = 0; i < k; i++)
      st[i] = i;
   for (i = k; i <= m; i++)
      st[i] = 0;    // With k=0, we do only this.

   UniqueArray<double> buf_store;
   buf_store.SetLength(m + 1);
   double *buf = buf_store.get();

   NTL::Vec<int64_t> in_vec_mem;
   in_vec_mem.SetLength(n + 1);
   int64_t *in_vec = in_vec_mem.elts();

   UniqueArray<double> max_b_store;
   max_b_store.SetLength(m + 1);
   double *max_b = max_b_store.get();
   // std::cout << "LLL: after creating UniqueArray's \n";

   for (i = 0; i < m; i++)
      max_b[i] = max_abs(B1[i], n);

   int64_t in_float = 0;
   int64_t rst;
   int64_t counter;
   int64_t start_over;
   int64_t trigger_index;
   int64_t small_trigger;
   int64_t cnt;
   // int64_t m_orig = m;
   int64_t rr_st = 0;   // One less than in NTL.
   int64_t max_k = -1;
   int64_t swap_cnt = 0;
   // long prec = RR::precision();

   while (k < m) {
      // std::cout << "ll_LLL FP64 enter while k < m with k = " << k << "\n";
      if (k > max_k) {
         max_k = k;
         swap_cnt = 0;
      }
      if (k < rr_st)
         rr_st = k;  // Both are 1 less than in NTL.
      if (st[k] == k)
         rst = 0;   // 1 less than in NTL.
      else
         rst = k;
      if (st[k] < st[k + 1])
         st[k + 1] = st[k];

      //std::cout << "LLL64: before ComputeGS, B1[k][1] = " << B1[k][1] << "\n";
      ComputeGS64(B, B1, mu, b, c, k, n, bound, st[k], buf);
      CheckFinite(&c[k]);
      st[k] = k;
      //std::cout << "After ComputeGS64, mu[k] = " << mu[k][0] << "  " << mu[k][1] << "  " << mu[k][2] << "  " << mu[k][3] << "\n";

      if (swap_cnt > 200000) {
         cerr << "LLL_FP64: swap loop? \n";
         abort();
         // In NTL, there is more stuff here. Removed for the `int64_t` version.
      }
      counter = 0;
      trigger_index = k;
      small_trigger = 0;
      cnt = 0;
      int64_t thresh = 10;
      int64_t sz = 0, new_sz;
      do { // size reduction
           //std::cout << "do loop: k = " <<  k << ",  counter = " << counter << "\n";
         counter++;
         if ((counter & 127) == 0) {   // Should be 127
            new_sz = 0;
            for (j = 0; j < n; j++)
               new_sz += NumBits(B[k][j]);
            if ((counter >> 7) == 1 || new_sz < sz) {  // Should be 7
               sz = new_sz;
            } else {
               cerr << "LLL_FP64 sz = " << sz
                     << " not smaller; infinite loop? \n";
               cerr << "new_sz = " << new_sz << ",  counter = " << counter
                     << ",  k = " << k << "\n";
               std::cout << " Matrix B = \n" << B << "\n";
               abort();
            }
         }
         Fc1 = 0;
         start_over = 0;

         for (j = rst - 1; j >= 0; j--) { // both j and rst are 1 less than in NTL.
            t1 = fabs(mu[k][j]);
            //std::cout << "entered for loop: j =  " <<  j << "  \n";
            //std::cout << "mu[k,j] =  " <<  mu[k][j] << "  t1 =  " <<  t1 << "  \n";
            if (t1 > half_plus_fudge) {
               // std::cout << "we have t1 > half_plus_fudge, j =  " <<  j << ",  Fc1 = " << Fc1 << "\n";
               if (!Fc1) {
                  // std::cout << "Before j > trigger_index, we have j = " << j << ",  trigger_index = " << trigger_index << "\n";
                  if (j > trigger_index
                        || (j == trigger_index && small_trigger)) {
                     cnt++;
                     if (cnt > thresh) {
                        //std::cout << "inc_red_fudge():  cnt= " << cnt << ", thresh = " << thresh
                        //      << ",  dim= " << n << " \n";
                        if (log_red <= 15) {
                           while (log_red > 10)
                              inc_red_fudge();
                           half_plus_fudge = 0.5 + red_fudge;
                           std::cout
                                 << "Skipping the part removed for the `long` version \n";
                           // Some part removed here for the `int64_t` version
                        } else {
                           std::cout << "... log_red > 15, in the else, put cnt = 0. \n";
                           inc_red_fudge();
                           half_plus_fudge = 0.5 + red_fudge;
                           cnt = 0;
                        }
                     }
                  }
                  trigger_index = j;
                  small_trigger = (t1 < 4);
                  Fc1 = 1;
                  if (k < rr_st)
                     rr_st = k;
                  //std::cout << "ll_LLL FP64 calling RowTransformStart, k = "
                  //      << k << ",  rr_st = " << rr_st << "\n";
                  RowTransformStart64(B1[k], in_vec, in_float, n);
               }
               mu1 = mu[k][j];
               //std::cout << "Before row transform, mu1 = " << mu1 << " \n";
               if (mu1 >= 0)
                  mu1 = ceil(mu1 - 0.5);
               else
                  mu1 = floor(mu1 + 0.5);

               double *mu_k = mu[k];
               double *mu_j = mu[j];

               if (mu1 == 1) {
                  for (i = 0; i < j - 1; i++)
                     mu_k[i] -= mu_j[i];
               } else if (mu1 == -1) {
                  for (i = 0; i < j - 1; i++)
                     mu_k[i] += mu_j[i];
               } else {
                  for (i = 0; i < j - 1; i++)
                     mu_k[i] -= mu1 * mu_j[i];
               }
               mu_k[j] -= mu1;
               conv(MU, mu1);
               //std::cout << "Before row transform, mu1 = " << mu1 << " \n";

               int64_t T, MU2 = MU;
               for (i = 0; i < n; i++) {
                  T = MU2 * B1[j][i];
                  B1[k][i] -= T;
               }
               // std::cout << "Before row transform, mu1 = " << mu1 << " \n";
               // We have `in_float=1` if all entries of B1[k] are in [-TR_BND, TR_BND].
               // The change must be on vector B[k].
               // RowTransform64(B[k], B[j], MU, n);
               //  RowTransform64(B[k], B[j], MU, n, B1[k], B1[j], in_vec,
               //         max_b[k], max_b[j], in_float);
               // std::cout << "After row transform, MU = " << MU << " \n";
               // std::cout << "Basis after row transform: \n" << B << "\n";
            }
         }
         //std::cout << "ll_LLL FP64 before if Fc1 \n";
         if (Fc1) {
            // std::cout << "ll_LLL FP64 inside `if(Fc1)` \n";
            NTL::Vec<int64_t> temp = B[k];
            RowTransformFinish64(temp, B1[k], in_vec, n);
            B[k] = temp;
            max_b[k] = max_abs(B1[k], n);
            b[k] = InnerProductD(B1[k], B1[k], n);
            CheckFinite(&b[k]);
            ComputeGS64(B, B1, mu, b, c, k, n, bound, 0, buf);
            CheckFinite(&c[k]);
            rst = k;
            // std::cout << "After ComputeGS in (Fc1), rst = " << rst << ",  max_b[k]= "
            //     << max_b[k] << ", did_rr_gs=  " << ", b[k]=  " << b[k] << "\n";
            // std::cout << "After ComputeGS in (Fc1), mu[k] = " << mu[k][0] << "  " << mu[k][1] << "  " << mu[k][2] << "  " << mu[k][3] << "\n";
         }
         // std::cout << "End of loop, B = " <<  B << "  \n";
      } while (Fc1 || start_over);  // End of `do` loop.
      // std::cout << "ll_LLL FP64 after while Fc1, k = " << k << "  b[k] = " << b[k] << "\n";
      // std::cout << "Basis after while Fc1 \n" << B << "\n";

      if (b[k] == 0) {
         for (i = k; i < m - 1; i++) {
            // swap(B[i], B[i+1]);
            B[i].swap(B[i + 1]);
            tp = B1[i];
            B1[i] = B1[i + 1];
            B1[i + 1] = tp;
            t1 = b[i];
            b[i] = b[i + 1];
            b[i + 1] = t1;
            t1 = max_b[i];
            max_b[i] = max_b[i + 1];
            max_b[i + 1] = t1;
         }
         for (i = k; i <= m; i++)
            st[i] = 0;
         if (k < rr_st)
            rr_st = k;
         m--;
         continue;
      }
// test LLL reduction condition
// std::cout << "Test reduction condition  \n";
      if (k > 0
            && delta * c[k - 1]
                  > c[k] + mu[k][k - 1] * mu[k][k - 1] * c[k - 1]) {
         // swap rows k, k-1
         swap(B[k], B[k - 1]);
         tp = B1[k];
         B1[k] = B1[k - 1];
         B1[k - 1] = tp;
         tp = mu[k];
         mu[k] = mu[k - 1];
         mu[k - 1] = tp;
         t1 = b[k];
         b[k] = b[k - 1];
         b[k - 1] = t1;
         t1 = max_b[k];
         max_b[k] = max_b[k - 1];
         max_b[k - 1] = t1;
         k--;
         NumSwaps++;
         swap_cnt++;
      } else {
         k++;
      }
   }
   return m;
}

double LLL_FP64(NTL::Mat<long> &B, double delta, long m, long n, NTL::Vec<double> *sqlen) {
   if (m == 0)
      m = B.NumRows();
   if (n == 0)
      n = B.NumCols();
   NumSwaps = 0;
   if (delta < 0.50 || delta >= 1)
      LogicError("LLL_FP: bad delta");
   long i, j;
   long new_m;
   long quit = 0;
   init_red_fudge();

   Unique2DArray<double> B1_store;
   B1_store.SetDims(m, n);
   double **B1 = B1_store.get();  // approximates B by a Mat of `double`

   Unique2DArray<double> mu_store;
   mu_store.SetDims(m, m);
   double **mu = mu_store.get();

   UniqueArray<double> c_store;
   c_store.SetLength(m);
   double *c = c_store.get();  // squared lengths of Gramm-Schmidt basis vectors

   UniqueArray<double> b_store;
   b_store.SetLength(m + 1);
   double *sqlen2 = b_store.get();  // squared lengths of basis vectors
   // This sqlen2 is usually the same as sqlen, but we have to do this
   // because we are not sure how much space has been reserved for sqlen.
   // and we do not want to change the pointer sqlen that is passed!

   for (i = 0; i < m; i++)
      for (j = 0; j < n; j++) {
         conv(B1[i][j], B[i][j]);  // Converts from int_64 to double
         CheckFinite(&B1[i][j]);
      }
   for (i = 0; i < m; i++) {
      sqlen2[i] = InnerProductD(B1[i], B1[i], n);  // Square norms in double.
      CheckFinite(&sqlen2[i]);
   }
   // std::cout << "LLL FP64 before ll_LLL  \n";
   // Indices in B1, mu, sqlen2 start at 0, which is 1 less than in NTL.
   // Note that the Mat B passed here may be larger than m x n.
   new_m = ll_LLL_FP64(B, delta, B1, mu, sqlen2, c, m, n, 0, quit);
   // std::cout << "LLL FP64 after ll_LLL  \n";

   // In this version, we leave the zero rows at the bottom.
   // The new_m independent basis vectors will be at the top of `B`.
   // Put shortest nonzero vector in first place.
   long imin = 0;
   double minSqlen = sqlen2[0];
   for (i = 1; i < new_m; i++)
      if (sqlen2[i] < minSqlen) {
         minSqlen = sqlen2[i];
         imin = i;
      };
   if (imin > 0) {
      NTL::swap(B[0], B[imin]);
      std::swap(sqlen2[0], sqlen2[imin]);
   }
   if (sqlen) {
      //if (sqlen->length() < new_m)
      //   sqlen->SetLength(new_m);
      for (i = 0; i < min(new_m, sqlen->length()); i++)
         (*sqlen)[i] = sqlen2[i];
   }
   // std::cout << "In LLL FP64 after the swaps:  ";
   // std::cout << "sqlen2[0] = " << sqlen2[0] << "\n";
   // std::cout << "Inside LLL, after swaps: sqlen[0] = " << sqlen[0] << "\n";
   return sqlen2[0];
}

//  BKZ   =====================================================================

static vec_double BKZConstant;

// Used only if `prune > 0`.
static void ComputeBKZConstant(long beta, long p) {
   const double c_PI = 3.14159265358979323846264338328;
   const double LogPI = 1.14472988584940017414342735135;

   BKZConstant.SetLength(beta - 1);   // Index starts at 0.
   vec_double Log;
   Log.SetLength(beta);
   long i, j, k;
   double x, y;

   for (i = 0; i < beta; i++)
      Log[i] = log(double(i + 1));
   for (i = 1; i <= beta - 1; i++) {
// First, we compute x = gamma(i/2)^{2/i}
      k = i / 2;
      if ((i & 1) == 0) { // i even
         x = 0;
         for (j = 1; j <= k; j++)
            x = x + Log[j - 1];
         x = x * (1 / double(k));
         x = exp(x);
      } else { // i odd
         x = 0;
         for (j = k + 2; j <= 2 * k + 2; j++)
            x = x + Log[j - 1];
         x = 0.5 * LogPI + x - 2 * (k + 1) * Log[1]; // ln 2
         x = x * (2.0 / double(i));
         x = exp(x);
      }
// Second, we compute y = 2^{2*p/i}
      y = -(2 * p / double(i)) * Log[1];
      y = exp(y);
      BKZConstant[i - 1] = x * y / c_PI;
   }
}

static vec_double BKZThresh;

// Used only if `prune > 0`.  Same beta as in NTL.
static void ComputeBKZThresh(double *c, long beta) {
   BKZThresh.SetLength(beta - 1);
   long i;
   double x;
   x = 0;
   for (i = 0; i < beta - 1; i++) {
      x += log(c[i]);
      BKZThresh[i] = exp(x / double(i + 1)) * BKZConstant[i];
      if (!IsFinite(&BKZThresh[i]))
         BKZThresh[i] = 0;
   }
}


double BKZ_FP64(NTL::Mat<long> &BB, double delta, long beta, long prune,
        long m, long n, Vec<double> *sqlen) {
   if (m == 0)
      m = BB.NumRows();
   if (n == 0)
      n = BB.NumCols();
   NumSwaps = 0;
   if (delta < 0.50 || delta >= 1)
      LogicError("BKZ_FPZZ: bad delta");
   if (beta < 2)
      LogicError("BKZ_FPZZ: bad block size");
   long m_orig = m;
   long i, j;
   long MU;
   double t1;
   double *tp;
   init_red_fudge();

   NTL::Mat<long> B;  // A copy of the used part of BB, plus one extra row.
   B.SetDims(m + 1, n);
   for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
         B[i][j] = BB[i][j];
      }
   }
//std::cout << " In BKZ, Matrix B = \n" << B << "\n";

// Here, the entries used in B1, mu, c, b, all start at 0.
// In NTL, they start at 1.
   Unique2DArray<double> B1_store;
   B1_store.SetDims(m + 2, n);
   double **B1 = B1_store.get();  // approximates B

   Unique2DArray<double> mu_store;
   mu_store.SetDims(m + 2, m);
   double **mu = mu_store.get();

   UniqueArray<double> c_store;
   c_store.SetLength(m + 2);
   double *c = c_store.get();  // squared lengths of Gramm-Schmidt basis vectors

   UniqueArray<double> b_store;
   b_store.SetLength(m + 2);
   double *b = b_store.get();  // squared lengths of basis vectors

   double cbar;

// The used entries of these arrays start at 1.
   UniqueArray<double> ctilda_store;
   ctilda_store.SetLength(m + 2);
   double *ctilda = ctilda_store.get();

   UniqueArray<double> vvec_store;
   vvec_store.SetLength(m + 2);
   double *vvec = vvec_store.get();

   UniqueArray<double> yvec_store;
   yvec_store.SetLength(m + 2);
   double *yvec = yvec_store.get();

   UniqueArray<double> uvec_store;
   uvec_store.SetLength(m + 2);
   double *uvec = uvec_store.get();

   UniqueArray<double> utildavec_store;
   utildavec_store.SetLength(m + 2);
   double *utildavec = utildavec_store.get();

   UniqueArray<long> Deltavec_store;
   Deltavec_store.SetLength(m + 2);
   long *Deltavec = Deltavec_store.get();

   UniqueArray<long> deltavec_store;
   deltavec_store.SetLength(m + 2);
   long *deltavec = deltavec_store.get();

   long quit = 0;
   double eta;  // Always 0 here.
   long new_m;
// Here, i, j, jj will be 1 less than in NTL.
// z, jj, kk, h, s, t, beta  are the same.
   long z, jj, kk;
   long s, t;
   long h;

// Indices for B1 and b start at 0 here.
   for (i = 0; i <= m; i++)
      for (j = 0; j < n; j++) {
         conv(B1[i][j], B[i][j]);
         CheckFinite(&B1[i][j]);
      }
   for (i = 0; i < m; i++) {
      b[i] = InnerProductD(B1[i], B1[i], n);
      CheckFinite(&b[i]);
   }
// Here we first perform LLL on the initial set of vectors.
   m = ll_LLL_FP64(B, delta, B1, mu, b, c, m, n, 0, quit);

   unsigned long NumIterations = 0;
   unsigned long NumTrivial = 0;
   unsigned long NumNonTrivial = 0;
   unsigned long NumNoOps = 0;
   long clean = 1;

   if (m < m_orig) {
      for (i = m_orig; i >= m + 1; i--) {
         // swap i, i-1
         swap(B[i], B[i - 1]);
      }
   }
   if (!quit && m > 1) {
      if (beta > m)
         beta = m;  // Block size must not exceed dimension.
      if (prune > 0)
         ComputeBKZConstant(beta, prune);
      z = 0;
      jj = -1;  // Same values as in NTL for z, kk, one less for jj.
      while (z < m - 1) {
         jj++;       // Will start with jj = 0.
         kk = min(jj + beta, m);
         if (jj == m - 1) {
            jj = 0;
            kk = beta;
            clean = 1;
         }
         if (prune > 0)
            ComputeBKZThresh(&c[jj], kk - jj);
         cbar = c[jj];
         utildavec[jj + 1] = uvec[jj + 1] = 1;
         yvec[jj + 1] = vvec[jj + 1] = 0;
         Deltavec[jj + 1] = 0;
         s = t = jj + 1;  // Same values as in NTL.
         deltavec[jj + 1] = 1;
         for (i = jj + 2; i <= kk + 1; i++) {
            // Same index i as in NTL.
            ctilda[i] = uvec[i] = utildavec[i] = yvec[i] = 0;
            Deltavec[i] = 0;
            vvec[i] = 0;
            deltavec[i] = 1;
         }
         // std::cout << " In BKZ, cbar = " << cbar <<  "\n";
         // long enum_cnt = 0;
         while (t <= kk) {
            ctilda[t] = ctilda[t + 1]
                  + (yvec[t] + utildavec[t]) * (yvec[t] + utildavec[t])
                        * c[t - 1];
            ForceToMem(&ctilda[t]);  // prevents an infinite loop
            if (prune > 0 && t > jj + 1)
               eta = BKZThresh(t - jj - 1);
            else
               eta = 0.0;
            if (ctilda[t] < cbar - eta) {
               if (t > jj + 1) {
                  //std::cout << " Decrease t = " << t << "\n";
                  t--;  // t decreases here
                  t1 = 0;
                  for (i = t + 1; i <= s; i++)
                     t1 += utildavec[i] * mu[i - 1][t - 1]; // Indices in mu are 1 less.
                  yvec[t] = t1;
                  t1 = -t1;
                  if (t1 >= 0)
                     t1 = ceil(t1 - 0.5);
                  else
                     t1 = floor(t1 + 0.5);
                  utildavec[t] = vvec[t] = t1;
                  Deltavec[t] = 0;
                  if (utildavec[t] > -yvec[t])
                     deltavec[t] = -1;
                  else
                     deltavec[t] = 1;
               } else {
                  cbar = ctilda[jj + 1];
                  for (i = jj + 1; i <= kk; i++) {
                     uvec[i] = utildavec[i];
                  }
               }
            } else {
               t++;   // t increases here
               s = max(s, t);
               if (t < s)
                  Deltavec[t] = -Deltavec[t];
               if (Deltavec[t] * deltavec[t] >= 0)
                  Deltavec[t] += deltavec[t];
               utildavec[t] = vvec[t] + Deltavec[t];
            }
         }
         NumIterations++;
         //std::cout << " NumIterations = " << NumIterations << "  ********** \n";
         h = min(kk + 1, m);   // A number of rows, same value as in NTL.
         if ((delta - 8 * red_fudge) * c[jj] > cbar) {
            clean = 0;
            // we treat the case that the new vector is b_s (jj+1 < s <= kk)
            // as a special case that appears to occur most of the time.
            s = 0;
            // std::cout << " BKZ_FPZZ find s, jj = " << jj << ", kk = " << kk << ", i = " << i << "\n";
            for (i = jj + 2; i <= kk; i++) {
               if (uvec[i] != 0) {
                  if (s == 0)
                     s = i;
                  else
                     s = -1;
               }
            }
            if (s == 0)
               LogicError("BKZ_FPZZ: internal error, s==0.");
            if (s > 0) {
               // special case
               NumTrivial++;
               for (i = s - 1; i > jj; i--) {
                  // This i is one less than in NTL.
                  // swap i, i-1
                  swap(B[i - 1], B[i]);
                  // B1 and b start at 0 instead of 1.
                  tp = B1[i - 1];
                  B1[i - 1] = B1[i];
                  B1[i] = tp;
                  t1 = b[i - 1];
                  b[i - 1] = b[i];
                  b[i] = t1;
               }
               // cerr << "special case\n";
               // h is the same as in NTL.
               new_m = ll_LLL_FP64(B, delta, B1, mu, b, c, h, n, jj, quit);
               if (new_m != h)
                  LogicError("BKZ_FPZZ: internal error, new_m != h");
               if (quit)
                  break;
            } else {
               // the general case, here s < 0
               NumNonTrivial++;
               for (i = 0; i < n; i++)
                  conv(B[m][i], 0);            // Put 0 in last row.
               for (i = jj + 1; i <= kk; i++) {
                  if (uvec[i] == 0)
                     continue;
                  conv(MU, uvec[i]);
                  // This changes row B[m].
                  RowTransformAdd(B[m], B[i - 1], MU, n);
               }
               //printBB1(B, B1, m);
               for (i = m; i > jj; i--) {
                  // swap i, i-1
                  swap(B[i - 1], B[i]);
                  tp = B1[i - 1];
                  B1[i - 1] = B1[i];
                  B1[i] = tp;
                  t1 = b[i - 1];
                  b[i - 1] = b[i];
                  b[i] = t1;
               }
               // std::cout << " After swap, jj = " << jj << "\n";
               //printBB1(B, B1, m);
               for (i = 0; i < n; i++) {
                  conv(B1[jj][i], B[jj][i]);
                  CheckFinite(&B1[jj][i]);
               }
               b[jj] = InnerProductD(B1[jj], B1[jj], n);
               CheckFinite(&b[jj]);
               //std::cout << " After B1[jj] <-- B[jj]  \n";
               //printBB1(B, B1, m);

               if (b[jj] == 0)
                  LogicError("BKZ_FP64: internal error, b[jj]==0");
               // remove linear dependencies
               // cerr << "general case\n";
               // Here, jj starts at 0.
               // The matrices B and B1 must have one more row!
               new_m = ll_LLL_FP64(B, delta, B1, mu, b, c, kk + 1, n, jj,
                     quit);
               if (new_m != kk) {
                  //std::cout << " new_m = " << new_m << ", kk = " << kk << "\n";
                  LogicError("BKZ_FP64: internal error, new_m != kk");
               }
               // remove zero vectors
               for (i = kk + 1; i <= m; i++) {
                  // swap i, i-1
                  swap(B[i - 1], B[i]);
                  tp = B1[i - 1];
                  B1[i - 1] = B1[i];
                  B1[i] = tp;
                  t1 = b[i - 1];
                  b[i - 1] = b[i];
                  b[i] = t1;
               }
               if (h > kk) {
                  // extend reduced basis
                  new_m = ll_LLL_FP64(B, delta, B1, mu, b, c, h, n, h - 1,
                        quit);
                  if (new_m != h)
                     LogicError("BKZ_FP64: internal error, new_m != h");
               }
            }
            z = 0;
         } else {
            NumNoOps++;
            if (!clean) {
               new_m = ll_LLL_FP64(B, delta, B1, mu, b, c, h, n, h - 1, quit);
               if (new_m != h)
                  LogicError("BKZ_FP64: internal error, new_m != h");
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
   double minlen = b[0];
   for (i = 1; i < m; i++)
      if (b[i] < minlen) {
         minlen = b[i];
         imin = i;
      };
   if (imin > 0) {
      swap(BB[0], BB[imin]);
      std::swap(b[0], b[imin]);
   }
   if (sqlen) {
//      if (sqlen->length() < m)
//         sqlen->SetLength(m);
      for (i = 0; i < min(m, sqlen->length()); i++)
         (*sqlen)[i] = b[i];
   }
// std::cout << " End of BKZ in LLL_FP64, Matrix B = \n" << B << "\n";
   return b[0];            // Square length of shortest vector.
}

NTL_END_IMPL

