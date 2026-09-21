/* specfunc/gamma/gamma_fallback.c
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
#include <gsl/gsl_sf_gamma.h>

#include "error.h"
#include "eval.h"

/*
 * Portable implementation of Gamma(x).
 *
 * Computes Gamma(x) over the real line using the kernels in
 * gamma_kernels.c.  This routine is used when the system tgamma()
 * function is unavailable.
 *
 * The argument is dispatched as follows:
 *
 *   x < -10
 *     Use the reflection formula with a scaled positive Gamma
 *     evaluation.
 *
 *   -10 <= x < -1
 *     Use recurrence from the central Gamma representation.
 *
 *   -1 <= x < 0
 *     Use the small-negative Gamma kernel.
 *
 *   0 < x < 1
 *     Evaluate
 *
 *       Gamma(x) = 1/x + (x-1) H(x),
 *
 *     avoiding division of Gamma(1+x) by x.
 *
 *   1 <= x <= 2
 *     Use the central representation
 *
 *       Gamma(x) = 1 + (x-1)(x-2) H(x-1).
 *
 *   2 < x <= 10
 *     Use recurrence from the central representation.
 *
 *   10 < x <= GAMMA_SCALED_X_MAX
 *     Use the scaled large-argument Gamma kernel and reconstruct the
 *     final value with ldexp().
 *
 *   x > GAMMA_SCALED_X_MAX
 *     Return overflow.
 *
 * Zero and the negative integers are poles and return GSL_EDOM.
 * Each successful branch returns an absolute error estimate in
 * result->err.
 */

static int
gamma_fallback_e(const double x, gsl_sf_result * result)
{
  int status;
  gsl_sf_result h, g;

  if (x == 0.0 || (x < 0.0 && x == floor(x)))
    DOMAIN_ERROR(result);
  else if (x < -10.0)
    return gamma_neg_reflection(x, result);
  else if (x < -1.0)
    return gamma_neg_recurrence(x, result);
  else if (x < 0.0)
    return gamma_neg_small(x, result);
  else if (x < 1.0)
    {
      const double one_over_x = 1.0 / x;
      const double eone_over_x = GSL_DBL_EPSILON * fabs(one_over_x);
      double term, eterm;

      status = gamma1p_H(x, &h);
      if (status)
        return status;

      term = (x - 1.0) * h.val;
      eterm = fabs(x - 1.0) * h.err
            + 2.0 * GSL_DBL_EPSILON * fabs(term)
            + GSL_DBL_EPSILON * fabs(x - 1.0) * fabs(h.val);

      result->val = one_over_x + term;
      result->err = eone_over_x + eterm
                  + GSL_DBL_EPSILON * fabs(result->val);

      if (gsl_isinf(result->val))
        {
          result->err = GSL_POSINF;
          return GSL_EOVRFLW;
        }

      return GSL_SUCCESS;
    }
  else if (x <= 2.0)
    return gamma1p_central(x - 1.0, result);
  else if (x <= 10.0)
    return gamma_low(x, result);
  else if (x <= GAMMA_SCALED_X_MAX)
    {
      int ge;

      status = gamma_2_inf_scaled(x, &g, &ge);
      if (status)
        return status;

      result->val = ldexp(g.val, ge);
      result->err = ldexp(g.err, ge);

      if (gsl_isinf(result->val))
        {
          result->err = GSL_POSINF;
          return GSL_EOVRFLW;
        }

      return GSL_SUCCESS;
    }
  else
    OVERFLOW_ERROR(result);
}
