/* cdf/betainv.c
 *
 * Copyright (C) 2004 Free Software Foundation, Inc.
 * Copyright (C) 2006, 2007 Brian Gough
 * Copyright (C) 2026 Patrick Alken
 * Written by Jason H. Stover.
 * Modified for GSL by Brian Gough
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
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Inverse of the beta distribution: solve I_x(a,b) = p or
   1 - I_x(a,b) = q for x.  The starting values follow

     A. Gil, J. Segura, N. M. Temme, "Efficient algorithms for the
     inversion of the cumulative central beta distribution", Numer.
     Algorithms 74(1), 77-91 (2017), sections 2.2.1 and 2.3,

   which use the uniform asymptotic expansion of N. M. Temme,
   "Asymptotic inversion of the incomplete beta function", J. Comput.
   Appl. Math. 41, 145-157 (1992), and are refined by the safeguarded
   Halley iteration in invert.h, applied to whichever of P and Q is
   smaller. */

#include <config.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_cdf.h>

#include "error.h"
#include "invert.h"

/* below this, P or Q is treated as 0 without evaluating it */
#define BETAINV_LOG_ZERO (-760.0)

/* ln of half the smallest subnormal number, and of 2^-54: a solution x
   (or 1 - x) below these rounds to 0 (or x rounds to 1) */
#define BETAINV_LOG_HALF_TRUE_MIN (-745.13)
#define BETAINV_LOG_HALF_ULP_ONE (-37.43)

typedef struct
{
  double a;
  double b;
  double lnbeta;                /* ln B(a,b) */
  double t;                     /* target value of P (or Q if upper) */
  double logt;
  int upper;                    /* solve Q(x) = t instead of P(x) = t */
} betainv_params;

/* residual for cdf_invert: g = P/p - 1 (or 1 - Q/q), u = f/f', d = f''/f' */
static int
betainv_fn (const double x, void *vparams, double *g, double *u, double *d)
{
  const betainv_params *p = (const betainv_params *) vparams;
  const double a = p->a, b = p->b;
  const double logx = log (x);
  const double log1mx = gsl_sf_log1p (-x);
  const double logpdf = (a - 1.0) * logx + (b - 1.0) * log1mx - p->lnbeta;
  double T = 0.0, logbound, r;

  if (p->upper)
    {
      /* Q(x) <= (1-x)^b / (b B(a,b)) max(1, x^(a-1)) */
      logbound = b * log1mx - log (b) - p->lnbeta + ((a < 1.0) ? (a - 1.0) * logx : 0.0);

      if (logbound > BETAINV_LOG_ZERO)
        T = gsl_cdf_beta_Q (x, a, b);
    }
  else
    {
      /* P(x) <= x^a / (a B(a,b)) max(1, (1-x)^(b-1)) */
      logbound = a * logx - log (a) - p->lnbeta + ((b < 1.0) ? (b - 1.0) * log1mx : 0.0);

      if (logbound > BETAINV_LOG_ZERO)
        T = gsl_cdf_beta_P (x, a, b);
    }

  r = T / p->t - 1.0;

  /* f = T - t, f' = +pdf (P) or -pdf (Q); (T - t)/pdf = r * t/pdf,
     computed in logarithms since r and t/pdf may be very large or small */
  *u = (r == 0.0) ? 0.0 : exp (log (fabs (r)) + p->logt - logpdf);
  if ((r < 0.0) != (p->upper != 0))
    *u = -*u;

  *g = p->upper ? -r : r;
  *d = (a - 1.0) / x - (b - 1.0) / (1.0 - x);

  return GSL_SUCCESS;
}

/* Solve s2 ln(x/s2) + c2 ln((1-x)/c2) = -eta^2/2 (GST 2017 eq. 2.10) for
   the branch x < s2 (eta < 0); returns ln(x).  For eta > 0 the caller
   swaps s2 and c2 and obtains ln(1-x). */
static double
betainv_logx_of_eta (const double eta, const double s2, const double c2)
{
  const double e2 = 0.5 * eta * eta;
  const double ls2 = log (s2), lc2 = log (c2);
  double hi = ls2;
  double lo = ls2 - (e2 - c2 * lc2) / s2 - 1.0;   /* F(lo) < 0 */
  double v;
  int i;

  /* leading term of the small x expansion (GST 2017 eq. 2.17-2.18) */
  v = (-e2 + s2 * ls2 + c2 * lc2) / s2;
  if (!(v > lo && v < hi))
    v = 0.5 * (lo + hi);

  for (i = 0; i < 100; i++)
    {
      const double ev = exp (v);
      const double F = s2 * (v - ls2) + c2 * (gsl_sf_log1p (-ev) - lc2) + e2;
      const double dF = s2 - c2 * ev / (1.0 - ev);
      double vn;

      if (F < 0.0)
        lo = v;
      else
        hi = v;

      vn = v - F / dF;
      if (!(vn > lo && vn < hi))
        vn = 0.5 * (lo + hi);

      /* Newton converges quadratically, so after a step below sqrt(eps)
         the error is at the eps level */
      if (fabs (vn - v) <= GSL_SQRT_DBL_EPSILON * (1.0 + fabs (v)))
        return vn;

      v = vn;
    }

  return v;
}

/* x from eta, for given s2 = a/(a+b), c2 = b/(a+b) */
static double
betainv_x_of_eta (const double eta, const double s2, const double c2)
{
  if (fabs (eta) < 1.0e-8)
    {
      return s2 + sqrt (s2 * c2) * eta;
    }
  else if (eta < 0.0)
    {
      return exp (betainv_logx_of_eta (eta, s2, c2));
    }
  else
    {
      return -gsl_sf_expm1 (betainv_logx_of_eta (-eta, c2, s2));
    }
}

/* Starting value for I_y(a,b) = t, t <= 1/2 (lower tail in y).  Returns
   y, or a value <= 0 if no estimate could be made. */
static double
betainv_start (const double t, const double a, const double b,
               const double lnbeta)
{
  int use_bounds;

  /* choice of method, after GST 2017 section 4 (scheme 1) */
  if (t <= 0.01)
    {
      if (a < 1.0)
        use_bounds = 1;
      else if (a < 30.0)
        use_bounds = (b < 1.0);
      else if (b < 0.5)
        use_bounds = 1;
      else if (b < 5.0)
        use_bounds = (t < 1.0e-4);
      else
        use_bounds = 0;
    }
  else
    {
      use_bounds = (a < 1.0 || b < 1.0);
    }

  if (use_bounds)
    {
      /* upper bound of GST 2017 eq. (2.30) by fixed point iteration from
         x = 0, computed in logarithms */
      const double c1 = (a + b) / (a + 1.0);
      const double c2 = (a + b) * (a + b + 1.0) / ((a + 1.0) * (a + 2.0));
      const double lnum = log (t) + log (a) + lnbeta;
      double x = 0.0;
      int i;

      for (i = 0; i < 3; i++)
        {
          double lx = (lnum - log (1.0 + x * (c1 + c2 * x))
                       - b * gsl_sf_log1p (-x)) / a;

          if (!(lx < 0.0))
            return -1.0;

          x = exp (lx);
        }

      return x;
    }
  else
    {
      /* uniform asymptotic inversion with the error function, GST 2017
         eqs. (2.9)-(2.14) */
      const double r = a + b;
      const double s2 = a / r, c2 = b / r;
      const double sc = sqrt (s2 * c2);
      const double eta0 = gsl_cdf_ugaussian_Pinv (t) / sqrt (r);
      const double x0 = betainv_x_of_eta (eta0, s2, c2);
      double eta1;

      if (fabs (eta0) < 1.0e-4 || x0 == s2)
        eta1 = -(1.0 - 2.0 * s2) / (3.0 * sc);
      else
        eta1 = log (eta0 * sc / (x0 - s2)) / eta0;

      return betainv_x_of_eta (eta0 + eta1 / r, s2, c2);
    }
}

/* solve I_x(a,b) = p (upper = 0) or 1 - I_x(a,b) = q (upper = 1), p + q = 1 */
static int
betainv (const double p, const double q, const double a, const double b,
         const int upper, double *x)
{
  betainv_params params;
  double x0;

  params.a = a;
  params.b = b;
  params.lnbeta = gsl_sf_lnbeta (a, b);
  params.t = upper ? q : p;
  params.logt = log (params.t);
  params.upper = upper;

  /* For small x, I_x(a,b) >= x^a / (a B(a,b)), so x <= (p a B(a,b))^(1/a);
     if this bound underflows, so does the solution.  Likewise
     1 - x <= (q b B(a,b))^(1/b), and if this is below 2^-54 the solution
     rounds to 1. */
  if ((log (p) + log (a) + params.lnbeta) / a < BETAINV_LOG_HALF_TRUE_MIN)
    {
      *x = 0.0;
      return GSL_SUCCESS;
    }
  else if ((log (q) + log (b) + params.lnbeta) / b < BETAINV_LOG_HALF_ULP_ONE)
    {
      *x = 1.0;
      return GSL_SUCCESS;
    }

  if (upper)
    {
      /* lower tail of 1 - x, with a and b exchanged */
      double y0 = betainv_start (q, b, a, params.lnbeta);
      x0 = (y0 > 0.0) ? 1.0 - y0 : 0.5;
      if (x0 >= 1.0)
        x0 = 1.0 - 0.5 * GSL_DBL_EPSILON;
    }
  else
    {
      x0 = betainv_start (p, a, b, params.lnbeta);
      if (!(x0 > 0.0))
        x0 = 0.5;
    }

  return cdf_invert (&betainv_fn, &params, x0, 0.0, 1.0, x);
}

double
gsl_cdf_beta_Pinv (const double P, const double a, const double b)
{
  double x;
  int status;

  if (!(P >= 0.0 && P <= 1.0))
    {
      CDF_ERROR ("P must be in range 0 <= P <= 1", GSL_EDOM);
    }
  else if (!(a > 0.0))
    {
      CDF_ERROR ("a must be positive", GSL_EDOM);
    }
  else if (!(b > 0.0))
    {
      CDF_ERROR ("b must be positive", GSL_EDOM);
    }
  else if (P == 0.0)
    {
      return 0.0;
    }
  else if (P == 1.0)
    {
      return 1.0;
    }

  /* for P > 1/2, 1 - P is exact and the upper tail is inverted */
  status = betainv (P, 1.0 - P, a, b, P > 0.5, &x);

  if (status)
    {
      CDF_ERROR ("inverse failed to converge", GSL_EFAILED);
    }

  return x;
}

double
gsl_cdf_beta_Qinv (const double Q, const double a, const double b)
{
  double x;
  int status;

  if (!(Q >= 0.0 && Q <= 1.0))
    {
      CDF_ERROR ("Q must be in range 0 <= Q <= 1", GSL_EDOM);
    }
  else if (!(a > 0.0))
    {
      CDF_ERROR ("a must be positive", GSL_EDOM);
    }
  else if (!(b > 0.0))
    {
      CDF_ERROR ("b must be positive", GSL_EDOM);
    }
  else if (Q == 0.0)
    {
      return 1.0;
    }
  else if (Q == 1.0)
    {
      return 0.0;
    }

  /* for Q > 1/2, 1 - Q is exact and the lower tail is inverted */
  status = betainv (1.0 - Q, Q, a, b, Q <= 0.5, &x);

  if (status)
    {
      CDF_ERROR ("inverse failed to converge", GSL_EFAILED);
    }

  return x;
}
