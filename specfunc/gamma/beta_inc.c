/* specfunc/gamma/beta_inc.c
 *
 * Copyright (C) 2007 Brian Gough
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 * Copyright (C) 2026 Patrick Alken
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/*
 * Regularized incomplete beta function I_x(a,b) and its complement
 * 1 - I_x(a,b) = I_{1-x}(b,a).
 *
 * References:
 *
 * The algorithms are implemented independently (no code was copied)
 * from the published descriptions in
 *
 *   V. Egorova, A. Gil, J. Segura, N. M. Temme, "Computation of the
 *   regularized incomplete Beta function", Dolomites Res. Notes Approx.
 *   16 (2023) 10-16 -- the continued fraction, the power series and the
 *   form of the power term x^a (1-x)^b / B(a,b) in terms of the scaled
 *   gamma function (their eq. 7);
 *
 *   A. R. DiDonato and A. H. Morris, Jr., "Algorithm 708: Significant
 *   digit computation of the incomplete beta function ratios", ACM
 *   Trans. Math. Software 18 (1992) 360-373 -- the choice of computing
 *   the smaller of I and 1 - I, and the region of the power series;
 *
 *   NIST Digital Library of Mathematical Functions, sections 8.17(v)
 *   (series 8.17.7, continued fraction 8.17.22).
 *
 * For a <= 0 or b <= 0 (not a negative integer) the generalized
 * function x^a 2F1(a, 1-b; a+1; x) / (a B(a,b)) is computed as before.
 */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_dd.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_log.h>

#include "error.h"
#include "constants.h"

#define BETA_INC_MAXITER 100000

static int
isnegint (const double x)
{
  return (x < 0) && (x == floor (x));
}

/* power term x^a y^b / B(a,b), with x and y = 1 - x given exactly as
   double-double numbers; returns the value, and its relative error in
   *relerr.  With Gamma(z) = sqrt(2 pi/z) z^z e^-z Gamma*(z),

     x^a y^b / B(a,b) = sqrt(ab / (2 pi (a+b))) Gamma*(a+b) /
                        (Gamma*(a) Gamma*(b)) exp(E),
     E = a (log(1+s) - s) + b (log(1+t) - t),
     s = x(a+b)/a - 1, t = y(a+b)/b - 1,

   (Egorova et al eq. 7; a s + b t = 0).  E is computed in double-double
   arithmetic, since for large |E| its rounding error in double precision
   would be |E| eps relative in the result. */
/* log(1+s) - s for s = x (p+q)/q - 1 (s = d/a or d/-b above): for
   s < -1/2, 1 + s may be tiny or subnormal, and log(1+s) is computed as
   log(x) + log1p(p/q) without forming 1 + s */
static gsl_dd
beta_inc_log1pmx (const gsl_dd s, const gsl_dd xd, const double p,
                  const double q)
{
  if (s.hi > -0.5)
    {
      return gsl_dd_log1pmx (s);
    }
  else
    {
      gsl_dd r;

      r.hi = p;
      r.lo = 0.0;
      r = gsl_dd_log1p (gsl_dd_div_d (r, q));

      return gsl_dd_sub (gsl_dd_add (gsl_dd_log (xd), r), s);
    }
}

static double
beta_inc_power_term (const double a, const double b, const gsl_dd xd,
                     const gsl_dd yd, double *relerr)
{
  const gsl_dd d = gsl_dd_sub (gsl_dd_mul_d (xd, b), gsl_dd_mul_d (yd, a));
  const gsl_dd s = gsl_dd_div_d (d, a);
  const gsl_dd t = gsl_dd_div_d (d, -b);
  const gsl_dd E = gsl_dd_add (gsl_dd_mul_d (beta_inc_log1pmx (s, xd, b, a), a),
                               gsl_dd_mul_d (beta_inc_log1pmx (t, yd, a, b), b));
  const double g = gsl_sf_gammastar (a + b)
    / (gsl_sf_gammastar (a) * gsl_sf_gammastar (b));
  const double pre = sqrt (a / (2.0 * M_PI)) * sqrt (b / (a + b)) * g;
  const gsl_dd ex = gsl_dd_exp (E);

  *relerr = 8.0 * GSL_DBL_EPSILON;

  return pre * (ex.hi + ex.lo);
}

/* continued fraction for I_x(a,b) x^-a y^-b a B(a,b), x < (a+1)/(a+b+2).
   DLMF 8.17.22 is 1/(1 + d1/(1 + d2/(1 + ...))); its even part is

     1/(B_1 + A_1/(B_2 + A_2/(B_3 + ...))),
     B_1 = 1 + d_1, B_{k+1} = 1 + d_{2k} + d_{2k+1}, A_k = -d_{2k-1} d_{2k}.

   Near the transition point 1 + d_1 = 1 - (a+b)x/(a+1) suffers
   cancellation (and so do the leading 1 + d_{2k+1}).  Written with
   lambda = a - (a+b) x = a y - b x, computed from the exact x and y,

     B_1 = (lambda + 1)/(a + 1),
     B_{m+1} = d_{2m} + [(a+m) lambda + m (2a + 3m) + a + 2m + m (a+m) y]
                        / ((a+2m)(a+2m+1)),
     A_k = (a+k-1)(a+b+k-1) k (b-k) x^2 / ((a+2k-2)(a+2k-1)^2 (a+2k)),

   with d_{2m} = m (b-m) x / ((a+2m-1)(a+2m)), and for lambda >= 0 the
   bracket has no cancellation.  Evaluated by the modified Lentz method. */
static int
beta_inc_cf (const double a, const double b, const gsl_dd xd,
             const gsl_dd yd, double *cf, int *nterms)
{
  const double tiny = 1.0e-300;
  const double x = xd.hi + xd.lo;
  const double y = yd.hi + yd.lo;
  const gsl_dd lam_dd = gsl_dd_sub (gsl_dd_mul_d (yd, a), gsl_dd_mul_d (xd, b));
  const double lam = lam_dd.hi + lam_dd.lo;
  double f, C, D;
  int k;

  f = (lam + 1.0) / (a + 1.0);
  if (fabs (f) < tiny)
    f = tiny;
  C = f;
  D = 0.0;

  for (k = 1; k <= BETA_INC_MAXITER; k++)
    {
      const double a2k = a + 2.0 * k;
      const double Ak = (a + k - 1.0) * (a + b + k - 1.0) * k * (b - k) * x * x
        / ((a2k - 2.0) * (a2k - 1.0) * (a2k - 1.0) * a2k);
      const double d2k = k * (b - k) * x / ((a2k - 1.0) * a2k);
      const double Bk = d2k + ((a + k) * lam + k * (2.0 * a + 3.0 * k) + a2k
                               + k * (a + k) * y) / (a2k * (a2k + 1.0));
      double delta;

      D = Bk + Ak * D;
      if (fabs (D) < tiny)
        D = tiny;
      C = Bk + Ak / C;
      if (fabs (C) < tiny)
        C = tiny;
      D = 1.0 / D;
      delta = C * D;
      f *= delta;

      if (fabs (delta - 1.0) < GSL_DBL_EPSILON)
        {
          *cf = 1.0 / f;
          *nterms = k;
          return GSL_SUCCESS;
        }
    }

  *cf = 1.0 / f;
  *nterms = k;
  return GSL_EMAXITER;
}

/* power series for I_x(a,b) x^-a a B(a,b) = sum_n (1-b)_n x^n a / ((a+n) n!),
   DLMF 8.17.7; for b x <= 0.7 or b <= 1, with x <= 0.7 */
static int
beta_inc_series (const double a, const double b, const double x,
                 double *sum, int *nterms)
{
  double t = 1.0, s = 1.0;
  int n;

  for (n = 1; n <= BETA_INC_MAXITER; n++)
    {
      double term;

      t *= (n - b) * x / n;
      term = t * a / (a + n);
      s += term;

      if (fabs (term) <= 0.5 * GSL_DBL_EPSILON * fabs (s))
        {
          *sum = s;
          *nterms = n;
          return GSL_SUCCESS;
        }
    }

  *sum = s;
  *nterms = n;
  return GSL_EMAXITER;
}

/* asymptotic expansion for large a and b <= 1, DLMF 8.18.3-8.18.7
   (uniform for x in (0,1]):

     I_x(a,b) = Gamma(a+b) / (Gamma(a) a^b) sum_k d_k Fs_k,

   with Fs_k = a^b F_k, xi = -ln x, z = a xi,

     Fs_0 = Q(b, z),  Fs_1 = (b - z)/a Fs_0 + z^b e^-z / (a Gamma(b)),
     a Fs_{k+1} = (k + b - z) Fs_k + k xi Fs_{k-1},

   and d_k the Taylor coefficients of phi(t)^(b-1) at t = xi, where
   phi(t) = (1 - e^-t)/t = int_0^1 e^-tu du, so that phi has the Taylor
   coefficients c_k = (-1)^k J_k / k!, J_k = int_0^1 u^k e^-xi u du
   (summed as a power series in xi; used for x > 0.7, xi < 0.36). */
#define BETA_INC_BGRAT_NMAX 40
#define BETA_INC_BGRAT_AMIN 15.0

static int
beta_inc_bgrat (const double a, const double b, const gsl_dd xd,
                double *val, double *relerr, int *nterms)
{
  const gsl_dd lnx = gsl_dd_log (xd);             /* -xi */
  const gsl_dd zd = gsl_dd_mul_d (lnx, -a);       /* a xi */
  const double xi = -(lnx.hi + lnx.lo);
  const double z = zd.hi + zd.lo;
  const double beta = b - 1.0;
  double c[BETA_INC_BGRAT_NMAX + 1], d[BETA_INC_BGRAT_NMAX + 1];
  double G, r, Fm1, F0, F1, S, prev;
  gsl_sf_result Qr;
  gsl_dd E, rd;
  int k, n, status;

  /* Taylor coefficients of phi at xi */
  for (k = 0; k <= BETA_INC_BGRAT_NMAX; k++)
    {
      double t = 1.0, J = 1.0 / (k + 1.0), fact = 1.0;
      int m;

      for (n = 1; n < 200; n++)
        {
          t *= -xi / n;
          J += t / (k + n + 1.0);
          if (fabs (t) < 1.0e-18 * J)
            break;
        }

      for (m = 2; m <= k; m++)
        fact *= m;

      c[k] = ((k & 1) ? -J : J) / fact;
    }

  /* d = phi^beta (J. C. P. Miller's recurrence for powers of a series) */
  d[0] = pow (c[0], beta);
  for (k = 1; k <= BETA_INC_BGRAT_NMAX; k++)
    {
      double sum = 0.0;
      int j;

      for (j = 1; j <= k; j++)
        sum += ((beta + 1.0) * j - k) * c[j] * d[k - j];

      d[k] = sum / (k * c[0]);
    }

  /* Gamma(a+b) / (Gamma(a) a^b) */
  G = sqrt (a / (a + b)) * gsl_sf_gammastar (a + b) / gsl_sf_gammastar (a)
    * exp ((a + b) * gsl_sf_log1p (b / a) - b);

  /* z^b e^-z / Gamma(b), with z in double-double */
  E = gsl_dd_sub (gsl_dd_mul_d (gsl_dd_log (zd), b), zd);
  E = gsl_dd_add_d (E, -gsl_sf_lngamma (b));
  rd = gsl_dd_exp (E);
  r = rd.hi + rd.lo;

  status = gsl_sf_gamma_inc_Q_e (b, z, &Qr);

  F0 = Qr.val;
  F1 = (b - z) / a * F0 + r / a;
  S = d[0] * F0 + d[1] * F1;
  prev = d[1] * F1;
  Fm1 = F0;

  for (k = 1; k < BETA_INC_BGRAT_NMAX; k++)
    {
      const double F2 = ((k + b - z) * F1 + k * xi * Fm1) / a;
      const double term = d[k + 1] * F2;

      S += term;
      Fm1 = F1;
      F1 = F2;

      /* for small b the odd d_k are nearly zero, so a single small term
         does not mean convergence: require two in succession */
      if (fabs (term) <= 0.5 * GSL_DBL_EPSILON * fabs (S)
          && fabs (prev) <= 0.5 * GSL_DBL_EPSILON * fabs (S))
        break;

      prev = term;
    }

  *val = G * S;
  *relerr = (10.0 + k) * GSL_DBL_EPSILON + Qr.err / fabs (Qr.val);
  *nterms = k;

  return (k >= BETA_INC_BGRAT_NMAX) ? GSL_EMAXITER : status;
}

/* I_x(a,b) for a, b > 0; x and y = 1 - x exact as double-double
   numbers, lny = ln y */
static int
beta_inc_lower (const double a, const double b, const gsl_dd xd,
                const gsl_dd yd, const double lny,
                double *val, double *err)
{
  const double x = xd.hi + xd.lo;
  double v, relerr;
  int nterms, status;

  if ((b <= 1.0 || b * x <= 0.7) && x <= 0.7)
    {
      double s, pre;

      status = beta_inc_series (a, b, x, &s, &nterms);

      /* x^a / (a B(a,b)) = power term / (a y^b) */
      pre = beta_inc_power_term (a, b, xd, yd, &relerr)
        * exp (-b * lny) / a;
      v = pre * s;
      relerr += (2.0 + b * fabs (lny) + nterms) * GSL_DBL_EPSILON;
    }
  else if (x * (a + b + 2.0) >= a + 1.0 && b <= 1.0)
    {
      if (a >= BETA_INC_BGRAT_AMIN)
        {
          status = beta_inc_bgrat (a, b, xd, &v, &relerr, &nterms);
        }
      else
        {
          /* shift a to a + n >= BETA_INC_BGRAT_AMIN with DLMF 8.17.20:
             I_x(a,b) = I_x(a+n,b) + sum_{i<n} t_i, t_0 = x^a y^b / (a B(a,b)),
             t_{i+1} = t_i x (a+b+i) / (a+i+1); all terms are positive */
          const int n = (int) ceil (BETA_INC_BGRAT_AMIN - a);
          double t, sum, relerr2;
          int i;

          status = beta_inc_bgrat (a + n, b, xd, &v, &relerr, &nterms);

          t = beta_inc_power_term (a, b, xd, yd, &relerr2) / a;
          sum = t;
          for (i = 0; i < n - 1; i++)
            {
              t *= x * (a + b + i) / (a + i + 1.0);
              sum += t;
            }

          relerr = (relerr * v + (relerr2 + n * GSL_DBL_EPSILON) * sum) / (v + sum);
          v += sum;
        }
    }
  else
    {
      double cf;

      status = beta_inc_cf (a, b, xd, yd, &cf, &nterms);
      v = beta_inc_power_term (a, b, xd, yd, &relerr) * cf / a;
      relerr += (2.0 + 2.0 * nterms) * GSL_DBL_EPSILON;
    }

  *val = v;
  *err = relerr * fabs (v);

  return status;
}

/* whether beta_inc_lower can evaluate I_x(a,b) in its region of rapid
   convergence: the power series, or the continued fraction for
   x < (a+1)/(a+b+2) (DLMF 8.17.22) */
static int
beta_inc_lower_ok (const double a, const double b, const double x)
{
  return ((b <= 1.0 || b * x <= 0.7) && x <= 0.7)
    || x * (a + b + 2.0) < a + 1.0
    || b <= 1.0;
}

/* both I = I_x(a,b) and J = 1 - I for a, b > 0, 0 < x < 1 */
static int
beta_inc_IJ (const double a, const double b, const double x,
             gsl_sf_result * I, gsl_sf_result * J)
{
  const gsl_dd yd = gsl_dd_two_sum (1.0, -x);   /* 1 - x exactly */
  const double y = yd.hi;
  const double lnx = log (x);
  const double lny = gsl_sf_log1p (-x);
  gsl_dd xd;
  int direct_I;
  double v, e;
  int status;

  xd.hi = x;
  xd.lo = 0.0;

  /* compute the smaller of I and J directly (I for x <= a/(a+b)), unless
     only the other one can be computed in its region of rapid
     convergence */
  direct_I = (x * (a + b) <= a);

  if (direct_I && !beta_inc_lower_ok (a, b, x) && beta_inc_lower_ok (b, a, y))
    direct_I = 0;
  else if (!direct_I && !beta_inc_lower_ok (b, a, y) && beta_inc_lower_ok (a, b, x))
    direct_I = 1;

  if (direct_I)
    status = beta_inc_lower (a, b, xd, yd, lny, &v, &e);
  else
    status = beta_inc_lower (b, a, yd, xd, lnx, &v, &e);

  /* a/(a+b) does not always separate the tails (e.g. for small b); if
     the value computed is close to 1, compute the other directly when
     possible, since 1 - v would amplify its error by v/(1-v) > 9 */
  if (v > 0.9 && (direct_I ? beta_inc_lower_ok (b, a, y) : beta_inc_lower_ok (a, b, x)))
    {
      direct_I = !direct_I;

      if (direct_I)
        status = beta_inc_lower (a, b, xd, yd, lny, &v, &e);
      else
        status = beta_inc_lower (b, a, yd, xd, lnx, &v, &e);
    }

  if (direct_I)
    {
      I->val = v;
      I->err = e;
      J->val = 1.0 - v;
      J->err = e + GSL_DBL_EPSILON * fabs (J->val);
    }
  else
    {
      J->val = v;
      J->err = e;
      I->val = 1.0 - v;
      I->err = e + GSL_DBL_EPSILON * fabs (I->val);
    }

  return status;
}

/* generalized function for a <= 0 or b <= 0 (not negative integers) */
static int
beta_inc_nonpositive (const double a, const double b, const double x,
                      gsl_sf_result * result)
{
  gsl_sf_result f, beta;
  int stat;
  const int stat_f = gsl_sf_hyperg_2F1_e (a, 1 - b, a + 1, x, &f);
  const int stat_beta = gsl_sf_beta_e (a, b, &beta);
  double prefactor = (pow (x, a) / a);

  result->val = prefactor * f.val / beta.val;
  result->err = fabs (prefactor) * f.err / fabs (beta.val)
    + fabs (result->val / beta.val) * beta.err;

  stat = GSL_ERROR_SELECT_2 (stat_f, stat_beta);

  return stat;
}

/* common driver: want_Q selects I (0) or 1 - I (1) */
static int
beta_inc_impl (const double a, const double b, const double x,
               const int want_Q, gsl_sf_result * result)
{
  if (x < 0.0 || x > 1.0)
    {
      DOMAIN_ERROR (result);
    }
  else if (isnegint (a) || isnegint (b) || isnegint (a + b))
    {
      DOMAIN_ERROR (result);
    }
  else if (x == 0.0)
    {
      result->val = want_Q ? 1.0 : 0.0;
      result->err = 0.0;
      return GSL_SUCCESS;
    }
  else if (x == 1.0)
    {
      result->val = want_Q ? 0.0 : 1.0;
      result->err = 0.0;
      return GSL_SUCCESS;
    }
  else if (a <= 0 || b <= 0)
    {
      int stat = beta_inc_nonpositive (a, b, x, result);

      if (want_Q)
        {
          result->val = 1.0 - result->val;
          result->err += GSL_DBL_EPSILON * fabs (result->val);
        }

      return stat;
    }
  else if (a == b && x == 0.5)
    {
      /* I_{1/2}(a,a) = 1/2 by symmetry (DLMF 8.17.4) */
      result->val = 0.5;
      result->err = 0.0;
      return GSL_SUCCESS;
    }
  else
    {
      gsl_sf_result I, J;
      int stat = beta_inc_IJ (a, b, x, &I, &J);

      *result = want_Q ? J : I;

      if (stat == GSL_EMAXITER)
        {
          GSL_ERROR ("too many iterations", GSL_EMAXITER);
        }

      /* as for gsl_sf_gamma_inc_P/Q, a result below GSL_DBL_MIN is
         returned without an underflow error; it is accurate to about one
         unit of the smallest subnormal, which the relative estimate does
         not capture */
      result->err += 2.0 * GSL_SF_TRUE_MIN;

      return GSL_SUCCESS;
    }
}

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_beta_inc_e (const double a, const double b, const double x,
                   gsl_sf_result * result)
{
  return beta_inc_impl (a, b, x, 0, result);
}

int
gsl_sf_beta_inc_Q_e (const double a, const double b, const double x,
                     gsl_sf_result * result)
{
  return beta_inc_impl (a, b, x, 1, result);
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

double
gsl_sf_beta_inc (const double a, const double b, const double x)
{
  EVAL_RESULT (gsl_sf_beta_inc_e (a, b, x, &result));
}

double
gsl_sf_beta_inc_Q (const double a, const double b, const double x)
{
  EVAL_RESULT (gsl_sf_beta_inc_Q_e (a, b, x, &result));
}
