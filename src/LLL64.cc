#include <cmath>
//#include <cstdlib>
//#include <cstdint>

#include <NTL/fileio.h>
#include <NTL/vec_double.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <latticetester/NTLWrap.h>
#include <latticetester/LLL64.h>

NTL_START_IMPL

//typedef NTL::Mat<int64_t> mat_long;
//typedef NTL::Vec<int64_t> vec_long;

typedef NTL::matrix<int64_t> mat_long;
typedef NTL::vector<int64_t> vec_long;

// static inline bool IsZero(long x)
//   {  return (x == 0); }

//static inline void MulSubFrom(long & x, long a, long b)
//   { x -= a * b; }

static inline void CheckFinite(double *p) {
   if (!IsFinite(p)) ResourceError("LLL64_FP: numbers too big...use LLL_XD");
}

static double InnerProduct(double *a, double *b, long n) {
   double s = 0;
   long i;
   for (i = 1; i <= n; i++)
      s += a[i]*b[i];
   return s;
}

static void InnerProduct(long &xx, const vec_long& a, const vec_long& b) {
   long x = 0;
   long n = min(a.length(), b.length());
   long i;
   for (i = 0; i < n; i++) {
      x += a[i] * b[i];
   }
   xx = x;
}

static void RowTransformSimple (vec_long& A, vec_long& B, long& MU1)
// x = x - y*MU
{
   register long MU = MU1;
   long n = A.length();
   long i;

   if (MU == 1) {
      for (i = 0; i < n; i++)
         sub(A[i], A[i], B[i]);
      return;
   }
   if (MU == -1) {
      for (i = 0; i < n; i++)
         add(A[i], A[i], B[i]);
      return;
   }
   if (MU == 0) return;

   for (i = 0; i < n; i++) {
	  A[i] -= MU * B[i];
      // MulSubFrom(A[i], B[i], MU);
   }
}


#define TR_BND (NTL_FDOUBLE_PRECISION/2.0)
// Just to be safe!!

static double max_abs(double *v, long n)
{
   long i;
   double res, t;
   res = 0;
   for (i = 0; i < n; i++) {
      t = fabs(v[i]);
      if (t > res) res = t;
   }
   return res;
}


static void RowTransformStart(double *a, long *in_a, long& in_float, long n) {
   long i;
   long inf = 1;

   for (i = 0; i < n; i++) {
      in_a[i] = (a[i] < TR_BND && a[i] > -TR_BND);
      inf = inf & in_a[i];
   }
   in_float = inf;    // Returns 1 if no a[i] is too large.
}


static void RowTransformFinish(vec_long& A, double *a, long *in_a) {
   long n = A.length();
   long i;

   for (i = 0; i < n; i++) {
      if (in_a[i])  {
         conv(A[i], a[i]);
      }
      else {
         conv(a[i], A[i]);
         CheckFinite(&a[i]);
      }
   }
}


static void RowTransform2(vec_long& A, vec_long& B, const long& MU1) {
// x = x + y*MU

   register long T, MU;
   // long k;
   long n = A.length();
   long i;
   MU = MU1;

   if (MU == 1) {
      for (i = 1; i <= n; i++)
         add(A[i], A[i], B[i]);
      return;
   }
   if (MU == -1) {
      for (i = 1; i <= n; i++)
         sub(A[i], A[i], B[i]);
      return;
   }
   if (MU == 0) return;

   for (i = 0; i < n; i++) {
       mul(T, B[i], MU);
       add(A[i], A[i], T);
   }
}

static void ComputeGS(mat_long& B, double **B1, double **mu, double *b,
               double *c, long k, double bound, long st, double *buf) {
   long n = B.NumCols();
   long i, j;
   double s, t1, y, t;

   long T1;
   long test;

   double *mu_k = mu[k];

   if (st < k) {
      for (i = 0; i < st; i++)
         buf[i] = mu_k[i]*c[i];
   }

   for (j = st; j <= k-1; j++) {
      s = InnerProduct(B1[k], B1[j], n);
      s=0;
      for (i = 0; i < n; i++)
         s += B1[k][i] * B1[j][i];

      // test = b[k]*b[j] >= NTL_FDOUBLE_PRECISION^2
      test = (b[k]/NTL_FDOUBLE_PRECISION >= NTL_FDOUBLE_PRECISION/b[j]);

      // test = test && s^2 <= b[k]*b[j]/bound,
      // but we compute it in a strange way to avoid overflow
      if (test && (y = fabs(s)) != 0) {
         t = y/b[j];
         t1 = b[k]/y;
         if (t <= 1)
            test = (t*bound <= t1);
         else if (t1 >= 1)
            test = (t <= t1/bound);
         else
            test = 0;
      }

      if (test) {
         InnerProduct(T1, (vec_long)B(k), (vec_long)B(j));
         conv(s, T1);
      }

      double *mu_j = mu[j];

      t1 = 0;
      for (i = 1; i <= j-1; i++) {
         t1 += mu_j[i]*buf[i];
      }
  
      mu_k[j] = (buf[j] = (s - t1))/c[j];
   }

#if (!NTL_EXT_DOUBLE)

   // Kahan summation 

   double c1;

   s = c1 = 0;
   for (j = 1; j <= k-1; j++) {
      y = mu_k[j]*buf[j] - c1;
      t = s+y;
      c1 = t-s;
      c1 = c1-y;
      s = t;
   }


#else

   s = 0;
   for (j = 1; j <= k-1; j++)
      s += mu_k[j]*buf[j];

#endif

   c[k] = b[k] - s;
}

NTL_CHEAP_THREAD_LOCAL double LLLStatusInterval = 900.0;
NTL_CHEAP_THREAD_LOCAL char *LLLDumpFile = 0;

static NTL_CHEAP_THREAD_LOCAL double red_fudge = 0;
static NTL_CHEAP_THREAD_LOCAL long log_red = 0;
// static NTL_CHEAP_THREAD_LOCAL long verbose = 0;

static NTL_CHEAP_THREAD_LOCAL unsigned long NumSwaps = 0;
static NTL_CHEAP_THREAD_LOCAL double RR_GS_time = 0;
static NTL_CHEAP_THREAD_LOCAL double StartTime = 0;
static NTL_CHEAP_THREAD_LOCAL double LastTime = 0;

void LLLStatus(long max_k, double t, long m, const mat_long& B) {
   cerr << "---- LLL_FP status ----\n";
   cerr << "elapsed time: ";
   PrintTime(cerr, t-StartTime);
   cerr << ", stage: " << max_k;
   cerr << ", rank: " << m;
   cerr << ", swaps: " << NumSwaps << "\n";

   long t1;
   long i;
   double prodlen = 0;

   for (i = 0; i < m; i++) {
      InnerProduct(t1, B[i], B[i]);
      if (t1 != 0)
         prodlen += log(t1);
   }
   const double two=2.0;
   cerr << "log of prod of lengths: " << prodlen/(two*log(two)) << "\n";

   if (LLLDumpFile) {
      cerr << "dumping to " << LLLDumpFile << "...";
      ofstream f;
      OpenWrite(f, LLLDumpFile);
      f << "[";
      for (i = 0; i < m; i++) {
         f << B[i] << "\n";
      }
      f << "]\n";
      f.close();
      cerr << "\n";
   }
   LastTime = t;
   
}

static void init_red_fudge()
{
   long i;
   log_red = long(0.50*NTL_DOUBLE_PRECISION);
   red_fudge = 1;
   for (i = log_red; i > 0; i--)
      red_fudge = red_fudge*0.5;
}

static void inc_red_fudge()
{
   red_fudge = red_fudge * 2;
   log_red--;
   cerr << "LLL_FP: warning--relaxing reduction (" << log_red << ")\n";
   if (log_red < 4)
      ResourceError("LLL_FP: too much loss of precision...stop!");
}


static long ll_LLL_FP(mat_long& B, double delta,
           double **B1, double **mu,
           double *b, double *c,
           long m, long init_k, long &quit) {
   long n = B.NumCols();
   long i, j, k, Fc1;
   long MU;
   double mu1;

   double t1;
   // long T1;
   double *tp;
   static double bound = 0;

   if (bound == 0) {
      // we tolerate a 15% loss of precision in computing
      // inner products in ComputeGS.

      bound = 1;
      for (i = 2*long(0.15*NTL_DOUBLE_PRECISION); i > 0; i--)
         bound = bound * 2;
   }

   double half_plus_fudge = 0.5 + red_fudge;

   quit = 0;
   k = init_k;


   vec_long st_mem;
   st_mem.SetLength(m+2);
   long *st = st_mem.elts();

   for (i = 1; i < k; i++)
      st[i] = i;

   for (i = k; i <= m+1; i++)
      st[i] = 1;

   UniqueArray<double> buf_store;
   buf_store.SetLength(m+1);
   double *buf = buf_store.get();

   vec_long in_vec_mem;
   in_vec_mem.SetLength(n+1);
   long *in_vec = in_vec_mem.elts();

   UniqueArray<double> max_b_store;
   max_b_store.SetLength(m+1);
   double *max_b = max_b_store.get();


   for (i = 1; i <= m; i++)
      max_b[i] = max_abs(B1[i], n);

   long in_float;

   long rst;
   long counter;
   long start_over;

   long trigger_index;
   long small_trigger;
   long cnt;

   //mat_RR rr_B1;
   //mat_RR rr_mu;
   //vec_RR rr_c;
   //vec_RR rr_b;

   // long m_orig = m;

   long rr_st = 1;

   long max_k = 0;

   // double tt;

   long swap_cnt = 0;


   while (k <= m) {

      if (k > max_k) {
         max_k = k;
         swap_cnt = 0;
      }

      if (k < rr_st) rr_st = k;

      if (st[k] == k)
         rst = 1;
      else
         rst = k;

      if (st[k] < st[k+1]) st[k+1] = st[k];
      ComputeGS(B, B1, mu, b, c, k, bound, st[k], buf);
      CheckFinite(&c[k]);
      st[k] = k;

      if (swap_cnt > 200000) {
         cerr << "LLL_FP: swap loop?\n";
      }

      counter = 0;
      trigger_index = k;
      small_trigger = 0;
      cnt = 0;

      long thresh = 10;
      long sz=0, new_sz;

      // long did_rr_gs = 0;


      do {
         // size reduction

         counter++;
         if ((counter & 127) == 0) {

            new_sz = 0;
            for (j = 1; j <= n; j++)
               new_sz += NumBits(B(k,j));

            if ((counter >> 7) == 1 || new_sz < sz) {
               sz = new_sz;
            }
            else {
               cerr << "LLL_FP: warning--infinite loop?\n";
            }
         }

         Fc1 = 0;
         start_over = 0;
   
         for (j = rst-1; j >= 1; j--) {
            t1 = fabs(mu[k][j]);
            if (t1 > half_plus_fudge) { 


               if (!Fc1) {
                  if (j > trigger_index || 
                      (j == trigger_index && small_trigger)) {

                     cnt++;

                     if (cnt > thresh) {
                        if (log_red <= 15) { 

                           while (log_red > 10)
                              inc_red_fudge();

                           half_plus_fudge = 0.5 + red_fudge;
                        }
                        else {
                           inc_red_fudge();
                           half_plus_fudge = 0.5 + red_fudge;
                           cnt = 0;
                        }
                     }
                  }

                  trigger_index = j;
                  small_trigger = (t1 < 4);

                  Fc1 = 1;
                  if (k < rr_st) rr_st = k;
                  RowTransformStart(B1[k], in_vec, in_float, n);
               }
                  

               mu1 = mu[k][j];
               if (mu1 >= 0)
                  mu1 = ceil(mu1-0.5);
               else
                  mu1 = floor(mu1+0.5);
   
               double *mu_k = mu[k];
               double *mu_j = mu[j];
   
               if (mu1 == 1) {
                  for (i = 1; i <= j-1; i++)
                     mu_k[i] -= mu_j[i];
               }
               else if (mu1 == -1) {
                  for (i = 1; i <= j-1; i++)
                     mu_k[i] += mu_j[i];
               }
               else {
                  for (i = 1; i <= j-1; i++)
                     mu_k[i] -= mu1*mu_j[i];
               }
   
               mu_k[j] -= mu1;
   
               conv(MU, mu1);
               RowTransformSimple ((vec_long)B[k], (vec_long)B[j], MU);
               // RowTransform(B(k), B(j), MU, B1[k], B1[j], in_vec,
               //             max_b[k], max_b[j], in_float);
            }
         }


         if (Fc1) {
            RowTransformFinish(B[k], B1[k], in_vec);
            max_b[k] = max_abs(B1[k], n);

               b[k] = InnerProduct(B1[k], B1[k], n);
               CheckFinite(&b[k]);

               ComputeGS(B, B1, mu, b, c, k, bound, 1, buf);
               CheckFinite(&c[k]);
            rst = k;
         }
      } while (Fc1 || start_over);

      if (b[k] == 0) {
         for (i = k; i < m; i++) {
            // swap i, i+1
            NTL::swap(B[i], B[i+1]);
            tp = B1[i]; B1[i] = B1[i+1]; B1[i+1] = tp;
            t1 = b[i]; b[i] = b[i+1]; b[i+1] = t1;
            t1 = max_b[i]; max_b[i] = max_b[i+1]; max_b[i+1] = t1;
         }

         for (i = k; i <= m+1; i++) st[i] = 1;
         if (k < rr_st) rr_st = k;

         m--;
         if (quit) break;
         continue;
      }
      if (quit) break;

      // test LLL reduction condition

      if (k > 1 && delta*c[k-1] > c[k] + mu[k][k-1]*mu[k][k-1]*c[k-1]) {
         // swap rows k, k-1
         swap(B[k], B[k-1]);
         tp = B1[k]; B1[k] = B1[k-1]; B1[k-1] = tp;
         tp = mu[k]; mu[k] = mu[k-1]; mu[k-1] = tp;
         t1 = b[k]; b[k] = b[k-1]; b[k-1] = t1;
         t1 = max_b[k]; max_b[k] = max_b[k-1]; max_b[k-1] = t1;

         k--;
         NumSwaps++;
         swap_cnt++;
         // cout << "-\n";
      }
      else {
         k++;
         // cout << "+\n";
      }
   }
   return m;
}


long LLL64_FP(NTL::matrix<int64_t>& B, double delta) {
   long m = B.NumRows();
   long n = B.NumCols();

   long i, j;
   long new_m, dep, quit;
   // long MU;
   // long T1;

   RR_GS_time = 0;
   NumSwaps = 0;
   if (delta < 0.50 || delta >= 1) LogicError("LLL_FP: bad delta");

   init_red_fudge();

   Unique2DArray<double> B1_store;
   B1_store.SetDimsFrom1(m+1, n+1);
   double **B1 = B1_store.get();  // approximates B

   Unique2DArray<double> mu_store;
   mu_store.SetDimsFrom1(m+1, m+1);
   double **mu = mu_store.get();

   UniqueArray<double> c_store;
   c_store.SetLength(m+1);
   double *c = c_store.get(); // squared lengths of Gramm-Schmidt basis vectors

   UniqueArray<double> b_store;
   b_store.SetLength(m+1);
   double *b = b_store.get(); // squared lengths of basis vectors

   for (i = 1; i <=m; i++)
      for (j = 1; j <= n; j++) {
         conv(B1[i][j], B(i, j));
         CheckFinite(&B1[i][j]);
      }
         
   for (i = 1; i <= m; i++) {
      b[i] = InnerProduct(B1[i], B1[i], n);
      CheckFinite(&b[i]);
   }

   new_m = ll_LLL_FP(B, delta, B1, mu, b, c, m, 1, quit);
   dep = m - new_m;
   m = new_m;

   if (dep > 0) {
      // for consistency, we move all of the zero rows to the front

      for (i = 0; i < m; i++) {
         swap(B[m+dep-i], B[m-i]);
      }
   }
   return m;
}


static vec_double BKZConstant;

static void ComputeBKZConstant(long beta, long p)
{
   const double c_PI = 3.14159265358979323846264338328;
   const double LogPI = 1.14472988584940017414342735135;

   BKZConstant.SetLength(beta-1);

   vec_double Log;
   Log.SetLength(beta);


   long i, j, k;
   double x, y;

   for (j = 1; j <= beta; j++) {
      x = log(double(j));  Log(j) = x;
   }
   for (i = 1; i <= beta-1; i++) {
      // First, we compute x = gamma(i/2)^{2/i}

      k = i/2;

      if ((i & 1) == 0) { // i even
         x = 0;
         for (j = 1; j <= k; j++)
            x = x + Log(j);
          
         x = x * (1/double(k));

         x = exp(x);
      }
      else { // i odd
         x = 0;
         for (j = k + 2; j <= 2*k + 2; j++)
            x = x + Log(j);

         x = 0.5*LogPI + x - 2*(k+1)*Log(2);
         x = x * (2.0/double(i));
         x = exp(x);
      }

      // Second, we compute y = 2^{2*p/i}

      y = -(2*p/double(i))*Log(2);
      y = exp(y);

      BKZConstant[i] = x*y/c_PI;
   }
}

/******************************************************/

long BKZ64_FP(mat_long& BB, double delta, long beta) {
   long m = BB.NumRows();
   long n = BB.NumCols();
   long m_orig = m;
   long i, j;
   long MU;

   double t1;
   // long T1;
   double *tp;

   RR_GS_time = 0;
   NumSwaps = 0;
   if (delta < 0.50 || delta >= 1) LogicError("BKZ_FP: bad delta");
   if (beta < 2) LogicError("BKZ_FP: bad block size");

   init_red_fudge();

   mat_long B;
   B = BB;

   B.SetDims(m+1, n);

   Unique2DArray<double> B1_store;
   B1_store.SetDimsFrom1(m+2, n+1);
   double **B1 = B1_store.get();  // approximates B


   Unique2DArray<double> mu_store;
   mu_store.SetDimsFrom1(m+2, m+1);
   double **mu = mu_store.get();

   UniqueArray<double> c_store;
   c_store.SetLength(m+2);
   double *c = c_store.get(); // squared lengths of Gramm-Schmidt basis vectors

   UniqueArray<double> b_store;
   b_store.SetLength(m+2);
   double *b = b_store.get(); // squared lengths of basis vectors

   double cbar;


   UniqueArray<double> ctilda_store;
   ctilda_store.SetLength(m+2);
   double *ctilda = ctilda_store.get();


   UniqueArray<double> vvec_store;
   vvec_store.SetLength(m+2);
   double *vvec = vvec_store.get();

   UniqueArray<double> yvec_store;
   yvec_store.SetLength(m+2);
   double *yvec = yvec_store.get();

   UniqueArray<double> uvec_store;
   uvec_store.SetLength(m+2);
   double *uvec = uvec_store.get();

   UniqueArray<double> utildavec_store;
   utildavec_store.SetLength(m+2);
   double *utildavec = utildavec_store.get();

   UniqueArray<long> Deltavec_store;
   Deltavec_store.SetLength(m+2);
   long *Deltavec = Deltavec_store.get();

   UniqueArray<long> deltavec_store;
   deltavec_store.SetLength(m+2);
   long *deltavec = deltavec_store.get();;

   long quit;
   long new_m;
   long z, jj, kk;
   long s, t;
   long h;
   double eta;


   for (i = 0; i <m; i++)
      for (j = 0; j < n; j++) {
         conv(B1[i][j], B[i][j]);
         CheckFinite(&B1[i][j]);
      }

         
   for (i = 0; i < m; i++) {
      b[i] = InnerProduct(B1[i], B1[i], n);
      CheckFinite(&b[i]);
   }
   m = ll_LLL_FP(B, delta, B1, mu, b, c, m, 1, quit);

   // double tt;
   // double enum_time = 0;
   unsigned long NumIterations = 0;
   unsigned long NumTrivial = 0;
   unsigned long NumNonTrivial = 0;
   unsigned long NumNoOps = 0;
   long clean = 1;

   if (m < m_orig) {
      for (i = m_orig; i >= m+1; i--) {
         // swap i, i-1
         swap(B[i], B[i-1]);
      }
   }

   if (!quit && m > 1) {
      if (beta > m) beta = m;
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
   
         cbar = c[jj];
         utildavec[jj] = uvec[jj] = 1;
   
         yvec[jj] = vvec[jj] = 0;
         Deltavec[jj] = 0;
   
         s = t = jj;
         deltavec[jj] = 1;
   
         for (i = jj; i <= kk; i++) {
            ctilda[i] = uvec[i] = utildavec[i] = yvec[i] = 0;
            Deltavec[i] = 0;
            vvec[i] = 0;
            deltavec[i] = 1;
         }

         while (t <= kk) {
            ctilda[t] = ctilda[t+1] + 
               (yvec[t]+utildavec[t])*(yvec[t]+utildavec[t])*c[t];

            ForceToMem(&ctilda[t]);  // prevents an infinite loop
   
            eta = 0;
            if (ctilda[t] < cbar - eta) {
               if (t > jj) {
                  t--;
                  t1 = 0;
                  for (i = t+1; i <= s; i++)
                     t1 += utildavec[i]*mu[i][t];
                  yvec[t] = t1;
                  t1 = -t1;
                  if (t1 >= 0)
                     t1 = ceil(t1-0.5);
                  else
                     t1 = floor(t1+0.5);
                  utildavec[t] = vvec[t] = t1;
                  Deltavec[t] = 0;
                  if (utildavec[t] > -yvec[t]) 
                     deltavec[t] = -1;
                  else
                     deltavec[t] = 1;
               }
               else {
                  cbar = ctilda[jj];
                  for (i = jj; i <= kk; i++) {
                     uvec[i] = utildavec[i];
                  }
               }
            }
            else {
               t++;
               s = max(s, t);
               if (t < s) Deltavec[t] = -Deltavec[t];
               if (Deltavec[t]*deltavec[t] >= 0) Deltavec[t] += deltavec[t];
               utildavec[t] = vvec[t] + Deltavec[t];
            }
         }
         NumIterations++;
   
         h = min(kk+1, m);
   
         if ((delta - 8*red_fudge)*c[jj] > cbar) {

            clean = 0;
            // we treat the case that the new vector is b_s (jj < s <= kk)
            // as a special case that appears to occur most of the time.
            s = 0;
            for (i = jj+1; i <= kk; i++) {
               if (uvec[i] != 0) {
                  if (s == 0)
                     s = i;
                  else
                     s = -1;
               }
            }
            if (s == 0) LogicError("BKZ_FP: internal error");
   
            if (s > 0) {
               // special case
               NumTrivial++;
               for (i = s-1; i >= jj; i--) {
                  // swap i, i-1
                  swap(B[i-1], B[i]);
                  tp = B1[i-1]; B1[i-1] = B1[i]; B1[i] = tp;
                  t1 = b[i-1]; b[i-1] = b[i]; b[i] = t1;
               }
               // cerr << "special case\n";
               new_m = ll_LLL_FP(B, delta, B1, mu, b, c, h, jj, quit);
               if (new_m != h) LogicError("BKZ_FP: internal error");
               if (quit) break;
            }
            else {
               // the general case

               NumNonTrivial++;
               for (i = 0; i < n; i++) conv(B[m][i], 0);
               for (i = jj-1; i < kk; i++) {
                  if (uvec[i] == 0) continue;
                  conv(MU, uvec[i]);
                  RowTransform2(B[m], B[i], MU);
               }
      
               for (i = m; i >= jj; i--) {
                  // swap i, i-1
                  swap(B[i-1], B[i]);
                  tp = B1[i-1]; B1[i-1] = B1[i]; B1[i] = tp;
                  t1 = b[i-1]; b[i-1] = b[i]; b[i] = t1;
               }
      
               for (i = 0; i < n; i++) {
                  conv(B1[jj][i], B[jj-1][i]);
                  CheckFinite(&B1[jj][i]);
               }
      
               b[jj] = InnerProduct(B1[jj], B1[jj], n);
               CheckFinite(&b[jj]);
      
               if (b[jj] == 0) LogicError("BKZ_FP: internal error, b[jj]==0");
      
               // remove linear dependencies
   
               // cerr << "general case\n";
               new_m = ll_LLL_FP(B, delta, B1, mu, b, c, kk+1, jj, quit);
              
               if (new_m != kk) LogicError("BKZ_FP: internal error, new_m != kk");

               // remove zero vector
      
               for (i = kk+1; i <= m; i++) {
                  // swap i, i-1
                  swap(B[i-1], B[i]);
                  tp = B1[i-1]; B1[i-1] = B1[i]; B1[i] = tp;
                  t1 = b[i-1]; b[i-1] = b[i]; b[i] = t1;
               }
               if (h > kk) {
                  // extend reduced basis
                  new_m = ll_LLL_FP(B, delta, B1, mu, b, c, h, h, quit);
   
                  if (new_m != h) LogicError("BKZ_FP: internal error, new_m != h");
                  if (quit) break;
               }
            }
            z = 0;
         }
         else {
            // LLL_FP
            // cerr << "progress\n";

            NumNoOps++;
            if (!clean) {
               new_m = 
                  ll_LLL_FP(B, delta, B1, mu, b, c, h, h, quit);
               if (new_m != h) LogicError("BKZ_FP: internal error, new_m != h");
               if (quit) break;
            }
            z++;
         }
      }
   }

   // clean up
   if (m_orig > m) {
      // for consistency, we move zero vectors to the front
      for (i = m+1; i <= m_orig; i++) {
         swap(B[i], B[i+1]);
      }
      for (i = 0; i < m; i++) {
         swap(B[m_orig-i], B[m-i]);
      }
   }
   B.SetDims(m_orig, n);
   BB = B;
   return m;
}


NTL_END_IMPL
