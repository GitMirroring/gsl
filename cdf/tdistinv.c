/* cdf/tdistinv.c
 *
 * Copyright (C) 2007, 2010 Brian Gough
 * Copyright (C) 2002 Jason H. Stover.
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

/* Inverse of the Student's t distribution.  By symmetry the tail
   equation Q(y) = p, y = |x| >= 0, p = min(P, 1 - P) is solved.  For
   nu = 1 and nu = 2 closed forms are used.  Otherwise the starting value
   is the inverse Cornish-Fisher expansion in the body of the
   distribution, or the leading term of the tail,

     Q(y) = 1/2 I_e(nu/2, 1/2) ~ e^(nu/2) / (nu B(nu/2, 1/2)),
     e = nu / (nu + y^2),

   and is refined by the safeguarded Halley iteration in invert.h
   applied to gsl_cdf_tdist_Q. */

#include <config.h>
#include <math.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_log.h>

#include "error.h"
#include "invert.h"

/* nu above which the inverse Cornish-Fisher expansion gives the
   starting value (gsl_cdf_tdist_P uses the expansion there, see
   tdist.c) */
#define TDISTINV_NU_CF 5e3

/* below this ln(e), the leading term of the series for I_e(nu/2, 1/2)
   is exact to rounding (as in tdist.c) */
#define TDISTINV_LOG_E_SMALL (-666.0)

static double
inv_cornish_fisher (double z, double nu)
{
  double a = 1 / (nu - 0.5);
  double b = 48.0 / (a * a);

  double cf1 = z * (3 + z * z);
  double cf2 = z * (945 + z * z * (360 + z * z * (63 + z * z * 4)));

  double y = z - cf1 / b + cf2 / (10 * b * b);

  double t = GSL_SIGN (z) * sqrt (nu * gsl_sf_expm1 (a * y * y));

  return t;
}

typedef struct
{
  double nu;
  double lognorm;               /* ln of the pdf normalisation */
  double t;                     /* target value of Q */
  double logt;
} tdistinv_params;

/* residual for cdf_invert: solve Q(y) = t for y >= 0 */
static int
tdistinv_fn (const double y, void *vparams, double *g, double *u, double *d)
{
  const tdistinv_params *p = (const tdistinv_params *) vparams;
  const double nu = p->nu;
  double log1pyy, logpdf, r;

  /* ln(1 + y^2/nu) without forming y^2 */
  if (y > sqrt (nu))
    {
      const double s = sqrt (nu) / y;
      log1pyy = 2.0 * log (y / sqrt (nu)) + gsl_sf_log1p (s * s);
    }
  else
    {
      log1pyy = gsl_sf_log1p (y * y / nu);
    }

  logpdf = p->lognorm - 0.5 * (nu + 1.0) * log1pyy;
  r = gsl_cdf_tdist_Q (y, nu) / p->t - 1.0;

  /* f = Q - t, f' = -pdf: u = (Q - t)/(-pdf) = -r t/pdf */
  *u = (r == 0.0) ? 0.0 : exp (log (fabs (r)) + p->logt - logpdf);
  if (r > 0.0)
    *u = -*u;

  *g = -r;
  *d = -(nu + 1.0) / (y + nu / y);     /* pdf'/pdf */

  return GSL_SUCCESS;
}

/* tail inverse: y >= 0 with Q(y) = p, 0 < p <= 1/2 */
static int
tdistinv_tail (const double p, const double nu, double *y)
{
  const double a = 0.5 * nu;
  const double lnbeta = gsl_sf_lnbeta (a, 0.5);
  const double loge0 = (log (2.0 * p) + log (a) + lnbeta) / a;
  tdistinv_params params;
  double y0;

  if (p == 0.5)
    {
      *y = 0.0;
      return GSL_SUCCESS;
    }
  else if (nu == 1.0)
    {
      /* Cauchy: y = cot(pi p), written to keep the argument <= pi/4 */
      *y = (p < 0.25) ? 1.0 / tan (M_PI * p) : tan (M_PI * (0.5 - p));
      return GSL_SUCCESS;
    }
  else if (nu == 2.0)
    {
      *y = (1.0 - 2.0 * p) / sqrt (2.0 * p * (1.0 - p));
      return GSL_SUCCESS;
    }
  else if (loge0 < TDISTINV_LOG_E_SMALL)
    {
      /* 2p = e^a / (a B(a, 1/2)) to rounding, and y = sqrt(nu/e); may
         overflow to +Inf */
      *y = exp (0.5 * (log (nu) - loge0));
      return GSL_SUCCESS;
    }

  if (nu > 1.0 && sqrt (M_PI * nu / 2.0) * p > pow (0.05, nu / 2.0))
    {
      /* body of the distribution: inverse Cornish-Fisher expansion */
      y0 = inv_cornish_fisher (gsl_cdf_ugaussian_Qinv (p), nu);

      /* for large nu gsl_cdf_tdist_Q uses gsl_cdf_ugaussian_Q, which is 0
         below GSL_DBL_MIN, so the expansion is the result */
      if (nu >= TDISTINV_NU_CF && p < GSL_DBL_MIN)
        {
          *y = y0;
          return GSL_SUCCESS;
        }
    }
  else
    {
      /* tail: y = sqrt(nu/e) from the leading term, with nu replaced by
         nu/(1 + nu/y^2) to account for the next terms */
      y0 = exp (0.5 * (log (nu) - loge0));
      y0 /= sqrt (1.0 + nu / (y0 * y0));
    }

  params.nu = nu;
  params.lognorm = gsl_sf_lngamma (0.5 * (nu + 1.0)) - gsl_sf_lngamma (a)
    - 0.5 * log (nu * M_PI);
  params.t = p;
  params.logt = log (p);

  return cdf_invert (&tdistinv_fn, &params, y0, 0.0, GSL_POSINF, y);
}

double
gsl_cdf_tdist_Pinv (const double P, const double nu)
{
  double y;
  int status;

  if (!(P >= 0.0 && P <= 1.0))
    {
      CDF_ERROR ("P must be in range 0 <= P <= 1", GSL_EDOM);
    }
  else if (!(nu > 0.0))
    {
      CDF_ERROR ("nu must be positive", GSL_EDOM);
    }
  else if (P == 0.0)
    {
      return GSL_NEGINF;
    }
  else if (P == 1.0)
    {
      return GSL_POSINF;
    }

  /* for P > 1/2, 1 - P is exact */
  status = tdistinv_tail ((P <= 0.5) ? P : 1.0 - P, nu, &y);

  if (status)
    {
      CDF_ERROR ("inverse failed to converge", GSL_EFAILED);
    }

  return (P < 0.5) ? -y : y;
}

double
gsl_cdf_tdist_Qinv (const double Q, const double nu)
{
  double y;
  int status;

  if (!(Q >= 0.0 && Q <= 1.0))
    {
      CDF_ERROR ("Q must be in range 0 <= Q <= 1", GSL_EDOM);
    }
  else if (!(nu > 0.0))
    {
      CDF_ERROR ("nu must be positive", GSL_EDOM);
    }
  else if (Q == 0.0)
    {
      return GSL_POSINF;
    }
  else if (Q == 1.0)
    {
      return GSL_NEGINF;
    }

  /* for Q > 1/2, 1 - Q is exact */
  status = tdistinv_tail ((Q <= 0.5) ? Q : 1.0 - Q, nu, &y);

  if (status)
    {
      CDF_ERROR ("inverse failed to converge", GSL_EFAILED);
    }

  return (Q <= 0.5) ? y : -y;
}
