/* specfunc/exp/exp_mult.c
 *
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
 * y exp(x), with and without input errors, as a double or in the
 * extended-range form val * 10^e10.
 *
 * The product is formed without log|y|: with y = m 2^ey (frexp, exact)
 * and x = k log(2) + f, |f| <= log(2)/2, reduced in double-double,
 *
 *   y exp(x) = (m exp(f)) 2^(ey + k),
 *
 * where the scaling by 2^(ey+k) is exact (ldexp).  The only roundings
 * are those of exp(f) and of the product m exp(f), so the result is
 * accurate to a few ulp however large |x| or |log|y|| is.  (Going
 * through exp(x + log|y|) loses eps |log|y|| to the rounding of the
 * logarithm.)  For the e10 form the same reduction is done with
 * log(10) in place of log(2).
 */

#include <config.h>
#include <limits.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_dd.h>
#include <gsl/gsl_sf_exp.h>

#include "constants.h"
#include "error.h"
#include "eval.h"

/* bound of the rounding error of m exp(f), in units of eps |result| */
#define EXP_MULT_REL      4.0

/* t - n (c + clo) in double-double, for a double-double t and an
   integer-valued n with |t - n c| small */
static gsl_dd
exp_mult_reduce(const gsl_dd t, const double n, const double c, const double clo)
{
  const gsl_dd p = gsl_dd_two_prod(n, c);
  const gsl_dd s = gsl_dd_two_sum(t.hi, -p.hi);
  return gsl_dd_two_sum(s.hi, s.lo + (t.lo - p.lo - n * clo));
}

/*
 * y exp(x) for finite y != 0: GSL_SUCCESS, or GSL_EOVRFLW /
 * GSL_EUNDRFLW (result beyond the normal range), without calling the
 * error handler.
 */
static int
exp_mult_impl(const double x, const double y, double * val)
{
  int ey;
  const double m = frexp(y, &ey);   /* 1/2 <= |m| < 1 */
  double kd;
  gsl_dd t, f;

  /* |y exp(x)| >= 2^(ey-1) exp(x), < 2^ey exp(x) */
  if (x > 1500.0)
    return GSL_EOVRFLW;
  if (x < -1500.0)
    return GSL_EUNDRFLW;

  t.hi = x;
  t.lo = 0.0;
  kd = floor(x / M_LN2 + 0.5);
  f = exp_mult_reduce(t, kd, M_LN2, GSL_SF_LN2_DD_LO);

  *val = ldexp(m * (exp(f.hi) * (1.0 + f.lo)), ey + (int) kd);

  if (!gsl_finite(*val))
    return GSL_EOVRFLW;
  if (fabs(*val) < GSL_DBL_MIN)
    return GSL_EUNDRFLW;
  return GSL_SUCCESS;
}

/*
 * y exp(x) = val 10^e10 for finite y != 0, with 1 <= |val| < 10 (up to
 * rounding at the boundaries); GSL_EOVRFLW / GSL_EUNDRFLW when e10 is
 * outside the int range.
 */
static int
exp_mult_e10_impl(const double x, const double y, double * val, int * e10)
{
  int ey;
  const double m = frexp(y, &ey);
  double nd;
  gsl_dd t, r;

  /* log|y exp(x)| = t + log|m|, t = x + ey log(2) */
  {
    const gsl_dd p = gsl_dd_two_prod((double) ey, M_LN2);
    t = gsl_dd_two_sum(x, p.hi);
    t.lo += p.lo + ey * GSL_SF_LN2_DD_LO;
  }

  nd = floor((t.hi + log(fabs(m))) / M_LN10);
  if (nd > INT_MAX - 1)
    return GSL_EOVRFLW;
  if (nd < INT_MIN + 1)
    return GSL_EUNDRFLW;

  r = exp_mult_reduce(t, nd, M_LN10, GSL_SF_LN10_DD_LO);
  *val = m * (exp(r.hi) * (1.0 + r.lo));
  *e10 = (int) nd;
  return GSL_SUCCESS;
}

/* no extended range needed (as in the original implementation, so that
   e10 = 0 for moderate arguments) */
static int
exp_mult_moderate(const double x, const double ay)
{
  return (x < 0.5 * GSL_LOG_DBL_MAX && x > 0.5 * GSL_LOG_DBL_MIN)
         && (ay < 0.8 * GSL_SQRT_DBL_MAX && ay > 1.2 * GSL_SQRT_DBL_MIN);
}

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_exp_mult_e(const double x, const double y, gsl_sf_result * result)
{
  int stat;

  if (gsl_isnan(x) || gsl_isnan(y))
    {
      DOMAIN_ERROR(result);
    }
  else if (y == 0.0)
    {
      result->val = 0.0;
      result->err = 0.0;
      return GSL_SUCCESS;
    }

  stat = exp_mult_impl(x, y, &result->val);
  if (stat == GSL_EOVRFLW)
    {
      OVERFLOW_ERROR(result);
    }
  else if (stat == GSL_EUNDRFLW)
    {
      UNDERFLOW_ERROR(result);
    }

  result->err = EXP_MULT_REL * GSL_DBL_EPSILON * fabs(result->val);
  return GSL_SUCCESS;
}

int
gsl_sf_exp_mult_err_e(const double x, const double dx,
                      const double y, const double dy,
                      gsl_sf_result * result)
{
  int stat;

  if (gsl_isnan(x) || gsl_isnan(y))
    {
      DOMAIN_ERROR(result);
    }
  else if (y == 0.0)
    {
      double e = 0.0;
      result->val = 0.0;
      if (dy != 0.0)
        {
          stat = exp_mult_impl(x, fabs(dy), &e);
          if (stat == GSL_EOVRFLW)
            e = GSL_POSINF;
          else if (stat == GSL_EUNDRFLW)
            e = GSL_DBL_MIN;
        }
      result->err = e;
      return GSL_SUCCESS;
    }

  stat = exp_mult_impl(x, y, &result->val);
  if (stat == GSL_EOVRFLW)
    {
      OVERFLOW_ERROR(result);
    }
  else if (stat == GSL_EUNDRFLW)
    {
      UNDERFLOW_ERROR(result);
    }

  result->err = fabs(result->val) * (EXP_MULT_REL * GSL_DBL_EPSILON + fabs(dx) + fabs(dy / y));
  return GSL_SUCCESS;
}

int
gsl_sf_exp_mult_e10_e(const double x, const double y, gsl_sf_result_e10 * result)
{
  int stat;

  if (gsl_isnan(x) || gsl_isnan(y))
    {
      DOMAIN_ERROR_E10(result);
    }
  else if (y == 0.0)
    {
      result->val = 0.0;
      result->err = 0.0;
      result->e10 = 0;
      return GSL_SUCCESS;
    }

  if (exp_mult_moderate(x, fabs(y)))
    {
      exp_mult_impl(x, y, &result->val);
      result->e10 = 0;
    }
  else
    {
      stat = exp_mult_e10_impl(x, y, &result->val, &result->e10);
      if (stat == GSL_EOVRFLW)
        {
          OVERFLOW_ERROR_E10(result);
        }
      else if (stat == GSL_EUNDRFLW)
        {
          UNDERFLOW_ERROR_E10(result);
        }
    }

  result->err = EXP_MULT_REL * GSL_DBL_EPSILON * fabs(result->val);
  return GSL_SUCCESS;
}

int
gsl_sf_exp_mult_err_e10_e(const double x, const double dx,
                          const double y, const double dy,
                          gsl_sf_result_e10 * result)
{
  int stat;

  if (gsl_isnan(x) || gsl_isnan(y))
    {
      DOMAIN_ERROR_E10(result);
    }
  else if (y == 0.0)
    {
      /* the error dy exp(x), in the scale of its own e10 */
      result->val = 0.0;
      result->err = 0.0;
      result->e10 = 0;
      if (dy != 0.0
          && exp_mult_e10_impl(x, fabs(dy), &result->err, &result->e10) != GSL_SUCCESS)
        {
          result->err = 0.0;
          result->e10 = 0;
        }
      return GSL_SUCCESS;
    }

  if (exp_mult_moderate(x, fabs(y)))
    {
      exp_mult_impl(x, y, &result->val);
      result->e10 = 0;
    }
  else
    {
      stat = exp_mult_e10_impl(x, y, &result->val, &result->e10);
      if (stat == GSL_EOVRFLW)
        {
          OVERFLOW_ERROR_E10(result);
        }
      else if (stat == GSL_EUNDRFLW)
        {
          UNDERFLOW_ERROR_E10(result);
        }
    }

  result->err = fabs(result->val) * (EXP_MULT_REL * GSL_DBL_EPSILON + fabs(dx) + fabs(dy / y));
  return GSL_SUCCESS;
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

double
gsl_sf_exp_mult(const double x, const double y)
{
  EVAL_RESULT(gsl_sf_exp_mult_e(x, y, &result));
}
