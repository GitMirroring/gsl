/* cdf/gammainv.c
 *
 * Copyright (C) 2003, 2007 Brian Gough
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
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Inverse of the gamma distribution: solve P(a, x) = p or Q(a, x) = q
   for x.  The starting values follow

     A. Gil, J. Segura, N. M. Temme, "Efficient and accurate algorithms
     for the computation and inversion of the incomplete gamma function
     ratios", SIAM J. Sci. Comput. 34(6), A2965-A2981 (2012), section 3,

   and are refined by the safeguarded Halley iteration in invert.h,
   applied to whichever of P and Q is smaller. */

#include <config.h>
#include <math.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_log.h>

#include "error.h"
#include "invert.h"

/* below this, P or Q is treated as 0 without calling gsl_sf_gamma_inc
   (which would report underflow); the smallest subnormal is e^-744.4 */
#define GAMMAINV_LOG_ZERO (-750.0)

typedef struct
{
  double a;
  double lngam;                 /* ln Gamma(a) */
  double lngam1;                /* ln Gamma(a+1) */
  double t;                     /* target value of P (or Q if upper) */
  double logt;
  int upper;                    /* solve Q(a,x) = t instead of P(a,x) = t */
} gammainv_params;

/* residual for cdf_invert: g = P/p - 1 (or 1 - Q/q), u = f/f', d = f''/f' */
static int
gammainv_fn (const double x, void *vparams, double *g, double *u, double *d)
{
  const gammainv_params *p = (const gammainv_params *) vparams;
  const double a = p->a;
  const double logx = log (x);
  const double logpdf = -x + (a - 1.0) * logx - p->lngam;
  double T = 0.0, r;
  gsl_sf_result result;

  if (p->upper)
    {
      /* Q(a,x) <= pdf * x / (x - a + 1) for x > a - 1, and <= pdf for a < 1 */
      double logbound = (x > a) ? logpdf + log (x / (x - a + 1.0)) : 0.0;

      if (logbound > GAMMAINV_LOG_ZERO)
        {
          gsl_sf_gamma_inc_Q_e (a, x, &result);
          T = result.val;
        }
    }
  else
    {
      /* P(a,x) < x^a / Gamma(a+1) */
      if (a * logx - p->lngam1 > GAMMAINV_LOG_ZERO)
        {
          gsl_sf_gamma_inc_P_e (a, x, &result);
          T = result.val;
        }
    }

  r = T / p->t - 1.0;

  /* f = T - t, f' = +pdf (P) or -pdf (Q); (T - t)/pdf = r * t/pdf,
     computed in logarithms since r and t/pdf may be very large or small */
  *u = (r == 0.0) ? 0.0 : exp (log (fabs (r)) + p->logt - logpdf);
  if ((r < 0.0) != (p->upper != 0))
    *u = -*u;

  *g = p->upper ? -r : r;
  *d = (a - 1.0) / x - 1.0;

  return GSL_SUCCESS;
}

/* solve lambda - 1 - ln(lambda) = eta^2 / 2 with sign(lambda - 1) = sign(eta) */
static double
gammainv_lambda (const double eta)
{
  const double s = 0.5 * eta * eta;
  double lambda;
  int i;

  if (fabs (eta) < 1.0e-8)
    return 1.0 + eta;

  if (eta < 0.0)
    {
      lambda = exp (-1.0 - s);
      if (eta > -1.0)
        lambda = GSL_MAX_DBL (lambda, 1.0 + eta + eta * eta / 3.0);
    }
  else
    {
      lambda = (eta < 1.0) ? 1.0 + eta + eta * eta / 3.0 : 1.0 + s + log (1.0 + s);
    }

  for (i = 0; i < 100; i++)
    {
      double F = lambda - 1.0 - log (lambda) - s;
      double dF = 1.0 - 1.0 / lambda;
      double ln = lambda - F / dF;

      if (ln <= 0.0)
        ln = 0.5 * lambda;

      /* Newton converges quadratically, so after a step below sqrt(eps)
         the error is at the eps level */
      if (fabs (ln - lambda) <= GSL_SQRT_DBL_EPSILON * lambda)
        return ln;

      lambda = ln;
    }

  return lambda;
}

/* starting value for P(a,x) = p, Q(a,x) = q, p + q = 1; sets *zero if the
   solution underflows */
static double
gammainv_start (const double p, const double q, const double a,
                const double lngam, const double lngam1, const int upper,
                int *zero)
{
  const double logp = upper ? gsl_sf_log1p (-q) : log (p);
  const double logq = upper ? log (q) : gsl_sf_log1p (-p);
  const double logr = (logp + lngam1) / a;  /* r = (p Gamma(1+a))^(1/a) */

  *zero = 0;

  if (logr < log (0.2 * (1.0 + a)))
    {
      /* small p: x = r + c2 r^2 + ... + c5 r^5, GST (3.2)-(3.3) */
      const double a1 = a + 1.0, a2 = a + 2.0, a3 = a + 3.0, a4 = a + 4.0;
      const double c2 = 1.0 / a1;
      const double c3 = (3.0 * a + 5.0) / (2.0 * a1 * a1 * a2);
      const double c4 = (8.0 * a * a + 33.0 * a + 31.0) / (3.0 * a1 * a1 * a1 * a2 * a3);
      const double c5 = ((((125.0 * a + 1179.0) * a + 3971.0) * a + 5661.0) * a + 2888.0)
        / (24.0 * a1 * a1 * a1 * a1 * a2 * a2 * a3 * a4);
      double r;

      if (logr < GAMMAINV_LOG_ZERO)
        {
          *zero = 1;
          return 0.0;
        }

      r = exp (logr);
      return r * (1.0 + r * (c2 + r * (c3 + r * (c4 + r * c5))));
    }
  else if (a < 10.0 && logq < -0.5 * a - lngam1)
    {
      /* small q: Q(a,x) ~ x^(a-1) e^-x / Gamma(a), so
         x - (a-1) ln x = -ln(q Gamma(a)), solved by fixed point iteration */
      const double eta0 = -logq - lngam;
      double x = eta0;
      int i;

      for (i = 0; i < 5 && x > 1.0; i++)
        x = eta0 + (a - 1.0) * log (x);

      return (x > 0.0) ? x : 1.0;
    }
  else if (a < 1.0)
    {
      /* lower bound x_l = (p Gamma(a+1))^(1/a), GST (3.8) */
      return exp (logr);
    }
  else
    {
      /* large a: uniform asymptotic inversion, GST (3.10)-(3.13);
         eta0 = z / sqrt(a) with z the normal quantile of p (or q) */
      const double z = upper ? gsl_cdf_ugaussian_Qinv (q) : gsl_cdf_ugaussian_Pinv (p);
      const double eta0 = z / sqrt (a);
      const double lambda0 = gammainv_lambda (eta0);
      double eps1, eta;

      if (fabs (eta0) < 1.0e-4)
        eps1 = -1.0 / 3.0;
      else
        eps1 = log (eta0 / (lambda0 - 1.0)) / eta0;

      eta = eta0 + eps1 / a;

      return a * gammainv_lambda (eta);
    }
}

/* solve P(a,x) = p (upper = 0) or Q(a,x) = q (upper = 1) with p + q = 1 */
static int
gammainv (const double p, const double q, const double a, const int upper,
          double *x)
{
  gammainv_params params;
  int zero;
  double x0;

  params.a = a;
  params.lngam = gsl_sf_lngamma (a);
  params.lngam1 = params.lngam + log (a);
  params.t = upper ? q : p;
  params.logt = log (params.t);
  params.upper = upper;

  x0 = gammainv_start (p, q, a, params.lngam, params.lngam1, upper, &zero);

  if (zero)
    {
      /* solution below the smallest representable number */
      *x = 0.0;
      return GSL_SUCCESS;
    }

  return cdf_invert (&gammainv_fn, &params, x0, 0.0, GSL_POSINF, x);
}

double
gsl_cdf_gamma_Pinv (const double P, const double a, const double b)
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
      return GSL_POSINF;
    }

  /* for P > 1/2, 1 - P is exact and the upper tail is inverted */
  status = gammainv (P, 1.0 - P, a, P > 0.5, &x);

  if (status)
    {
      CDF_ERROR ("inverse failed to converge", GSL_EFAILED);
    }

  return b * x;
}

double
gsl_cdf_gamma_Qinv (const double Q, const double a, const double b)
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
      return GSL_POSINF;
    }
  else if (Q == 1.0)
    {
      return 0.0;
    }

  /* for Q > 1/2, 1 - Q is exact and the lower tail is inverted */
  status = gammainv (1.0 - Q, Q, a, Q <= 0.5, &x);

  if (status)
    {
      CDF_ERROR ("inverse failed to converge", GSL_EFAILED);
    }

  return b * x;
}
