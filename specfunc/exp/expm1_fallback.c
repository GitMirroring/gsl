/* specfunc/exp/expm1_fallback.c
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
#include <gsl/gsl_sf_exp.h>

#include "constants.h"
#include "error.h"

static int expm1_fallback(const double x, gsl_sf_result * result);

/*
 * Portable fallback for exp(x) - 1.
 *
 * This implementation is based on the SunPro/fdlibm expm1 algorithm.
 * It avoids the severe cancellation in exp(x) - 1 for small x by
 * combining argument reduction, a rational approximation on a small
 * interval, and reconstruction formulas chosen to preserve accuracy.
 *
 * For arguments outside the primary interval, reduce
 *
 *   x = k log(2) + r + c,
 *
 * where |r| <= log(2)/2 and c is the correction for rounding in the
 * reduction.  The logarithm is split as
 *
 *   log(2) = LN2_HI + LN2_LO
 *
 * so that r and c can be formed accurately.
 *
 * On the reduced interval, the algorithm approximates exp(r) - 1
 * through the identity
 *
 *   r (exp(r) + 1) / (exp(r) - 1) = 2 + r^2/6 - r^4 R(r^2)/360,
 *
 * where R is represented by the polynomial with coefficients Q[].
 * The resulting correction is evaluated in a form designed to limit
 * rounding error.
 *
 * Finally, exp(x) - 1 is reconstructed from the reduced result and
 * 2^k.  Separate formulas are used for k = 0, k = +/-1, moderate k,
 * and large |k| to avoid cancellation.  ldexp() is used for exact
 * power-of-two scaling rather than direct manipulation of the binary64
 * exponent field used by the original fdlibm implementation.
 *
 * For |x| < 2^-54, expm1(x) rounds sufficiently closely to x that the
 * function returns x directly.  For sufficiently large negative x,
 * the correctly rounded value is -1.
 */

static int
expm1_fallback(const double x, gsl_sf_result * result)
{
  /*
   * Polynomial coefficients for the reduced-argument correction R.
   * These are the original SunPro/fdlibm coefficients.
   */
  static const double Q[5] = { -3.33333333333331316428e-02,  1.58730158725481460165e-03, -7.93650757867487942473e-05,
                                4.00821782732936239552e-06, -2.01099218183624371326e-07 };
  double z = x;
  double y, hi, lo, c, t, e, hxs, hfx, r1;
  int k;

  /* for x < -56 log(2), exp(x) < 2^-56 and expm1(x) rounds to -1 */
  if (x < -56.0 * GSL_SF_LN2_HI)
    {
      result->val = -1.0;
      result->err = 0.5 * GSL_DBL_EPSILON;
      return GSL_SUCCESS;
    }

  /*
   * Argument reduction.
   *
   * Write x = k log(2) + z + c, with |z| <= log(2)/2.
   * Splitting log(2) into high and low parts limits cancellation;
   * c retains the rounding residual from forming the reduced argument.
   */
  if (fabs(z) > 0.5 * GSL_SF_LN2_HI)
    {
      if (fabs(z) < 1.5 * GSL_SF_LN2_HI)
        {
          if (z > 0.0)
            {
              hi = z - GSL_SF_LN2_HI;
              lo = GSL_SF_LN2_LO;
              k = 1;
            }
          else
            {
              hi = z + GSL_SF_LN2_HI;
              lo = -GSL_SF_LN2_LO;
              k = -1;
            }
        }
      else
        {
          k = (int) (M_LOG2E * z + (z > 0.0 ? 0.5 : -0.5));
          t = (double) k;
          hi = z - t * GSL_SF_LN2_HI;
          lo = t * GSL_SF_LN2_LO;
        }

      z = hi - lo;
      c = (hi - z) - lo;
    }
  else if (fabs(z) < GSL_SF_DBL_EPSILON_4)
    {
      /*
       * For |x| < 2^-54, the correction expm1(x) - x begins with
       * x^2/2 and is too small to require the full kernel.
       */
      result->val = z;
      result->err = 0.5 * z * z / (1.0 - fabs(z));
      return GSL_SUCCESS;
    }
  else
    {
      k = 0;
      c = 0.0;
    }

  /*
   * Evaluate the reduced-argument kernel.  Here hxs = z^2/2 and r1
   * evaluates the polynomial correction.  The quantity e is arranged
   * so that, for k = 0,
   *
   *   expm1(z) = z - (z*e - z^2/2),
   *
   * with reduced rounding error.
   */
  hfx = 0.5 * z;
  hxs = z * hfx;
  r1 = 1.0 + hxs * (Q[0] + hxs * (Q[1] + hxs * (Q[2] + hxs * (Q[3] + hxs * Q[4]))));
  t = 3.0 - r1 * hfx;
  e = hxs * ((r1 - t) / (6.0 - z * t));

  if (k == 0)
    {
      result->val = z - (z * e - hxs);
      result->err = 2.0 * GSL_DBL_EPSILON * (fabs(z) + fabs(z*e) + fabs(hxs) + fabs(result->val));
      return GSL_SUCCESS;
    }

  /*
   * Fold the argument-reduction correction c into the kernel error
   * before reconstructing the result at the original scale.
   */
  e = (z * (e - c) - c);
  e -= hxs;

  if (k == -1)
    {
      /* Specialized reconstruction for scaling by 2^-1 */
      result->val = 0.5 * (z - e) - 0.5;
      result->err = 2.0 * GSL_DBL_EPSILON * (0.5 * (fabs(z) + fabs(e)) + 0.5 + fabs(result->val));
      return GSL_SUCCESS;
    }

  if (k == 1)
    {
      /*
       * Specialized reconstruction for scaling by 2.  Two algebraically
       * equivalent forms are used to avoid cancellation in different
       * parts of the reduced interval.
       */
      if (z < -0.25)
        {
          result->val = -2.0 * (e - (z + 0.5));
          result->err = 2.0 * GSL_DBL_EPSILON * (2.0 * (fabs(e) + fabs(z + 0.5)) + fabs(result->val));
        }
      else
        {
          result->val = 1.0 + 2.0 * (z - e);
          result->err = 2.0 * GSL_DBL_EPSILON * (1.0 + 2.0 * (fabs(z) + fabs(e)) + fabs(result->val));
        }

      return GSL_SUCCESS;
    }

  /*
   * Reconstruction for k <= -2 or k > 56.  In these ranges the direct
   * scaled form 2^k (1 + expm1(z)) - 1 is numerically suitable; the
   * intermediate ranges below require rearrangements that retain the
   * subtraction of 1 more accurately.
   */
  if (k <= -2 || k > 56)
    {
      y = 1.0 - (e - z);
      y = ldexp(y, k);
      result->val = y - 1.0;
      result->err = 3.0 * (GSL_DBL_EPSILON * fabs(y) +
                           GSL_DBL_EPSILON +
                           GSL_DBL_EPSILON * fabs(result->val));
      return GSL_SUCCESS;
    }

  if (k < 20)
    {
      /*
       * For moderate k, incorporate 1 - 2^-k before scaling.  This
       * preserves the subtraction of 1 in exp(x) - 1 without forming
       * exp(x) and subtracting 1 afterward.
       */
      t = 1.0 - ldexp(1.0, -k);
      y = t - (e - z);
      y = ldexp(y, k);
      result->val = y;
      result->err = 3.0 * GSL_DBL_EPSILON * (fabs(t) + fabs(e) + fabs(z) + fabs(y) + fabs(result->val));
    }
  else
    {
      /*
       * For 20 <= k <= 56, 2^-k is small enough that the alternate
       * reconstruction retains the low-order correction more accurately.
       */
      t = ldexp(1.0, -k);
      y = z - (e + t);
      y += 1.0;
      y = ldexp(y, k);
      result->val = y;
      result->err = 3.0 * (GSL_DBL_EPSILON * fabs(z) +
                           GSL_DBL_EPSILON * fabs(e) +
                           GSL_DBL_EPSILON * fabs(t) +
                           GSL_DBL_EPSILON * fabs(y) +
                           GSL_DBL_EPSILON * fabs(result->val));
    }

  return GSL_SUCCESS;
}
