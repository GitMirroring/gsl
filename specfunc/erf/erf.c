/* specfunc/erf/erf.c
 *
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

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_erf.h>

#include "constants.h"
#include "error.h"
#include "eval.h"

#ifndef HAVE_ERF
#include "erf_fallback.c"
#endif

static void square_twoprod(const double x, double * hi, double * lo);
static double exp_half_residual(const double r);
static double scaled_exp_half_square(const double h, const double r, const double erfcx);
static int erf_Q_positive(const double x, gsl_sf_result * result);

/*
 * Compute the error function erf(x).
 *
 * If a system erf() implementation is available, use it directly and
 * provide an estimate of its rounding error.  Otherwise, use the internal
 * fallback implementation.
 *
 * Returns GSL_SUCCESS on successful evaluation.
 */
int
gsl_sf_erf_e(const double x, gsl_sf_result * result)
{
#ifdef HAVE_ERF

  result->val = erf(x);
  result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
  return GSL_SUCCESS;

#else

  return erf_fallback(x, result);

#endif
}

/*
 * Compute the standard normal probability density
 *
 *   Z(x) = exp(-x^2/2) / sqrt(2*pi).
 *
 * The square x^2 is evaluated as a two-product decomposition to preserve
 * low-order bits in the exponential argument.  In the extreme tail, the
 * exponent is shifted by an integer multiple of log(2) and the result is
 * scaled back with ldexp() to preserve subnormal values.
 */

int
gsl_sf_erf_Z_e(const double x, gsl_sf_result * result)
{
  const double zero_x = 3.8580157609028419e1;
  const int scale_K = 48;
  const double ax = fabs(x);
  double h, r;

  if (ax >= zero_x)
    {
      result->val = 0.0;
      result->err = GSL_SF_TRUE_MIN;
      return GSL_SUCCESS;
    }

  square_twoprod(ax, &h, &r);

  if (ax < 35.0)
    {
      result->val = GSL_SF_INV_SQRT_2PI * exp(-0.5 * h) * exp(-0.5 * r);
      result->err = 8.0 * GSL_DBL_EPSILON * fabs(result->val);
      if (result->err == 0.0 && result->val != 0.0)
        result->err = GSL_SF_TRUE_MIN;

      return GSL_SUCCESS;
    }
  else
    {
      const double e1 = exp((double) scale_K * GSL_SF_LN2_HI - 0.5 * h);
      const double e2 = exp((double) scale_K * GSL_SF_LN2_LO - 0.5 * r);

      result->val = ldexp(GSL_SF_INV_SQRT_2PI * e1 * e2, -scale_K);
      result->err = 12.0 * GSL_DBL_EPSILON * fabs(result->val);
      if (result->err == 0.0 && result->val != 0.0)
        result->err = GSL_SF_TRUE_MIN;

      return GSL_SUCCESS;
    }
}

/*
 * Compute the standard normal upper-tail probability
 *
 *   Q(x) = 1/2 * erfc(x / sqrt(2)).
 *
 * Positive arguments are evaluated by erf_Q_positive().  For negative
 * arguments, the symmetry Q(-x) = 1 - Q(x) is used.  This avoids evaluating
 * erfc at a negative argument and gives the correct saturation to unity in
 * the negative tail.
 */

int
gsl_sf_erf_Q_e(const double x, gsl_sf_result * result)
{
  if (x == 0.0)
    {
      result->val = 0.5;
      result->err = 0.0;
      return GSL_SUCCESS;
    }

  if (x < 0.0)
    {
      int status;
      gsl_sf_result r;

      status = erf_Q_positive(-x, &r);
      if (status)
        return status;

      result->val = 1.0 - r.val;
      result->err = r.err + 2.0 * GSL_DBL_EPSILON * fabs(result->val);

      return GSL_SUCCESS;
    }
  else
    {
      return erf_Q_positive(x, result);
    }
}

/*
 * Compute the standard normal hazard function
 *
 *   h(x) = Z(x) / Q(x),
 *
 * where Z(x) is the standard normal probability density and Q(x) is the
 * upper-tail probability.
 *
 * For x < 0, evaluate Z(x)/Q(x) directly.  This is well conditioned because
 * Q(x) approaches unity in the negative tail.
 *
 * For x > 0, use the equivalent identity
 *
 *   h(x) = sqrt(2/pi) / erfcx(x/sqrt(2)),
 *
 * which avoids forming exponentially small quantities in the positive tail.
 * For sufficiently large positive x, the asymptotic behavior
 *
 *   h(x) = x + 1/x + O(x^-3)
 *
 * implies that h(x) rounds to x in binary64, so x is returned directly.
 */

int
gsl_sf_hazard_e(const double x, gsl_sf_result * result)
{
  int status;
  gsl_sf_result r, s;

  if (x == 0.0)
    {
      result->val = GSL_SF_SQRT2_PI;
      result->err = 0.0;
      return GSL_SUCCESS;
    }

  if (x < 0.0)
    {
      status = gsl_sf_erf_Z_e(x, &r);
      if (status)
        return status;

      status = gsl_sf_erf_Q_e(x, &s);
      if (status)
        return status;

      result->val = r.val / s.val;

      if (r.val != 0.0)
        {
          result->err = fabs(result->val) * (r.err / fabs(r.val) + s.err / fabs(s.val)) +
                        4.0 * GSL_DBL_EPSILON * fabs(result->val);
        }
      else
        {
          result->err = r.err / fabs(s.val) + 4.0 * GSL_DBL_EPSILON * fabs(result->val);
        }

      if (result->err == 0.0 && result->val != 0.0)
        result->err = GSL_SF_TRUE_MIN;

      return GSL_SUCCESS;
    }
  else
    {
      /* h(x) = x + 1/x + O(x^-3).  For x >= 2e8 the correction 1/x is
       * below one-half of a binary64 ulp of x, so the correctly rounded
       * result is x.  This also avoids overflow in the reciprocal of the
       * asymptotic erfcx value near DBL_MAX. */
      if (x >= 2.0e8)
        {
          result->val = x;
          result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
          return GSL_SUCCESS;
        }

      {
        const double z = x * M_SQRT1_2;

        status = gsl_sf_erfcx_e(z, &r);
        if (status)
          return status;

        result->val = GSL_SF_SQRT2_PI / r.val;
        result->err = fabs(result->val) * (r.err / fabs(r.val)) +
                      4.0 * GSL_DBL_EPSILON * fabs(result->val);

        return GSL_SUCCESS;
      }
    }
}

double gsl_sf_erf(const double x) { EVAL_RESULT(gsl_sf_erf_e(x, &result)); }
double gsl_sf_erf_Z(const double x) { EVAL_RESULT(gsl_sf_erf_Z_e(x, &result)); }
double gsl_sf_erf_Q(const double x) { EVAL_RESULT(gsl_sf_erf_Q_e(x, &result)); }
double gsl_sf_hazard(const double x) { EVAL_RESULT(gsl_sf_hazard_e(x, &result)); }

/*************************************************
 * INTERNAL ROUTINES                             *
 *************************************************/

/* product decomposition x*x = hi + lo. */
static void
square_twoprod(const double x, double * hi, double * lo)
{
  const double p = x * x;
  const double c = 134217729.0 * x; /* 2^27 + 1 */
  const double xh = c - (c - x);
  const double xl = x - xh;
  const double e = ((xh * xh - p) + 2.0 * xh * xl) + xl * xl;

  *hi = p;
  *lo = e;
}

/*
 * Compute exp(-r/2) for the small residual r in the decomposition
 * x*x = h + r, using a second-order Taylor expansion.
 */
static double
exp_half_residual(const double r)
{
  const double u = -0.5 * r;
  return 1.0 + u + 0.5 * u * u;
}

/*
 * Compute
 *
 *   1/2 * exp(-x^2/2) * erfcx
 *
 * from x*x = h + r.  The exponential is shifted by K binary powers before
 * multiplication by erfcx and restored with ldexp() to avoid premature
 * underflow in the extreme tail.
 */
static double
scaled_exp_half_square(const double h, const double r, const double erfcx)
{
  const int scale_K = 80;
  const double u = (double) scale_K * GSL_SF_LN2_LO - 0.5 * r;
  const double e1 = exp((double) scale_K * GSL_SF_LN2_HI - 0.5 * h);
  const double e2 = 1.0 + u + 0.5 * u * u;

  return ldexp(0.5 * e1 * e2 * erfcx, -scale_K);
}

/*
 * Compute Q(x) for x >= 0.
 *
 * In the central range, evaluate erfc(x/sqrt(2)) directly.  In the positive
 * tail, use
 *
 *   Q(x) = 1/2 * exp(-x^2/2) * erfcx(x/sqrt(2)),
 *
 * with a two-product decomposition of x*x to preserve the low-order bits of
 * the exponential argument.  The small residual exponential is evaluated by
 * a short polynomial.  In the extreme tail, the exponential is shifted by
 * an integer multiple of log(2) and restored with ldexp() to preserve
 * representable subnormal values.  Once Q(x) is below the binary64 rounding
 * threshold, zero is returned directly.
 */

static int
erf_Q_positive(const double x, gsl_sf_result * result)
{
  int status;
  const double z = x * M_SQRT1_2;
  const double zero_x = 3.848540833556734e1;
  gsl_sf_result s;
  double h, r;

  if (x >= zero_x)
    {
      result->val = 0.0;
      result->err = GSL_SF_TRUE_MIN;
      return GSL_SUCCESS;
    }

  if (x < 1.0)
    {
      status = gsl_sf_erfc_e(z, &s);
      if (status)
        return status;

      result->val = 0.5 * s.val;
      result->err = 0.5 * s.err + 4.0 * GSL_DBL_EPSILON * fabs(result->val);

      return GSL_SUCCESS;
    }

  square_twoprod(x, &h, &r);

  if (x < 37.0)
    {
      status = gsl_sf_erfcx_e(z, &s);
      if (status)
        return status;

      result->val = 0.5 * exp(-0.5 * h) * exp_half_residual(r) * s.val;
      result->err = fabs(result->val) * (s.err / fabs(s.val)) + 12.0 * GSL_DBL_EPSILON * fabs(result->val);

      return GSL_SUCCESS;
    }
  else
    {
      status = gsl_sf_erfcx_e(z, &s);
      if (status)
        return status;

      result->val = scaled_exp_half_square(h, r, s.val);
      result->err = fabs(result->val) * (s.err / fabs(s.val)) + 16.0 * GSL_DBL_EPSILON * fabs(result->val);

      if (result->err == 0.0 && result->val != 0.0)
        result->err = GSL_SF_TRUE_MIN;

      return GSL_SUCCESS;
    }
}
