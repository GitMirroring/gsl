/* specfunc/erf/erfc.c
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
#include <gsl/gsl_sf_log.h>

#include "constants.h"
#include "error.h"
#include "eval.h"

#ifndef HAVE_ERFC
#include "erfc_fallback.c"
#endif

/*
 * Compute the complementary error function erfc(x).
 *
 * If a system erfc() implementation is available, use it directly and
 * provide an estimate of its rounding error.  Otherwise, use the internal
 * fallback implementation.
 *
 * Returns GSL_SUCCESS on successful evaluation.
 */
int
gsl_sf_erfc_e(const double x, gsl_sf_result * result)
{
#ifdef HAVE_ERFC

  result->val = erfc(x);
  result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
  return GSL_SUCCESS;

#else

  return erfc_fallback(x, result);

#endif
}

/*
 * First positive binary64 x for which x*x overflows.  At and above this
 * point log(erfc(x)) is below the finite binary64 range, so the result is
 * -infinity and x*x must not be evaluated.
 */
#define LOG_ERFC_X_OVERFLOW 1.3407807929942597e154

/*
 * Compute log(erfc(x)).
 *
 * Near the origin, use log1p(-erf(x)) to avoid loss of precision when
 * erfc(x) is close to 1.  For negative x, use the symmetry relation
 * erfc(x) = 2 - erfc(-x) in a log1p formulation, with saturation to
 * log(2) when the correction is below double precision.  In the
 * ordinary positive range, evaluate log(erfc(x)) directly.  In the
 * positive tail, use
 *
 *   log(erfc(x)) = log(erfcx(x)) - x^2
 *
 * to avoid underflow in erfc(x).
 *
 * Returns GSL_SUCCESS on successful evaluation.
 */
int
gsl_sf_log_erfc_e(const double x, gsl_sf_result * result)
{
  int status;
  const double ax = fabs(x);
  gsl_sf_result r, s;

  /* preserve both +0 and -0 as the exact mathematical value 0 */
  if (x == 0.0)
    {
      result->val = 0.0;
      result->err = 0.0;
      return GSL_SUCCESS;
    }

  /* erfc(x) is near unity here; do not form log(erfc(x)). */
  if (ax < 0.5)
    {
      status = gsl_sf_erf_e(x, &r);
      if (status)
        return status;

      status = gsl_sf_log1p_e(-r.val, &s);
      if (status)
        return status;

      result->val = s.val;
      result->err = s.err + r.err / fabs(1.0 - r.val);

      if (result->err == 0.0 && result->val != 0.0)
        result->err = GSL_SF_TRUE_MIN;

      return GSL_SUCCESS;
    }

  /* On the negative side erfc(x) approaches 2. The log1p form avoids
     first rounding 2 - erfc(-x) when the correction becomes small. */
  if (x < 0.0)
    {
      const double axv = -x;

      if (axv >= 6.0)
        {
          result->val = M_LN2;
          result->err = GSL_DBL_EPSILON * fabs(result->val);
          return GSL_SUCCESS;
        }

      status = gsl_sf_erfc_e(axv, &r);
      if (status)
        return status;

      status = gsl_sf_log1p_e(-0.5 * r.val, &s);
      if (status)
        return status;

      result->val = M_LN2 + s.val;
      result->err = s.err + r.err / fabs(2.0 - r.val) + GSL_DBL_EPSILON * fabs(result->val);

      return GSL_SUCCESS;
    }

  /* In the ordinary positive range erfc itself is well conditioned. */
  if (x <= 5.0)
    {
      status = gsl_sf_erfc_e(x, &r);
      if (status)
        return status;

      result->val = log(r.val);
      result->err = r.err / fabs(r.val) + GSL_DBL_EPSILON * fabs(result->val);

      return GSL_SUCCESS;
    }

  /* In the positive tail, erfc underflows but erfcx remains representable.
     Once x^2 exceeds GSL_DBL_MAX, the true log(erfc(x)) is itself below the
     finite binary64 range, so return -infinity without evaluating x*x. */
  if (x >= LOG_ERFC_X_OVERFLOW)
    {
      result->val = GSL_NEGINF;
      result->err = 0.0;
      return GSL_SUCCESS;
    }

  status = gsl_sf_erfcx_e(x, &r);
  if (status)
    return status;

  {
    const double logr = log(r.val);
    const double x2 = x * x;

    result->val = logr - x2;
    result->err = r.err / fabs(r.val) +
                  GSL_DBL_EPSILON * fabs(logr) +
                  GSL_DBL_EPSILON * x2 +
                  GSL_DBL_EPSILON * fabs(result->val);
  }

  return GSL_SUCCESS;
}

double gsl_sf_erfc(const double x) { EVAL_RESULT(gsl_sf_erfc_e(x, &result)); }
double gsl_sf_log_erfc(const double x) { EVAL_RESULT(gsl_sf_log_erfc_e(x, &result)); }
