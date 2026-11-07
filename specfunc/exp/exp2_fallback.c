/* specfunc/exp/exp2_fallback.c
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

static int exp2_fallback(const double x, gsl_sf_result * result);

/*
 * C89 fallback for the base-2 exponential
 *
 *   exp2(x) = 2^x.
 *
 * The argument is reduced as
 *
 *   x = n + r,
 *
 * where n is the nearest integer and |r| <= 1/2.  The reduced
 * exponential is evaluated as
 *
 *   2^r = exp(r log(2)),
 *
 * using a high/low split of log(2) to reduce rounding error.  The
 * integer power 2^n is then applied with ldexp(), so the scaling is
 * exact except when the final result is subnormal.
 *
 * Near the upper and lower limits of the binary64 range, the reduction
 * is adjusted so that the intermediate scaling remains representable.
 * Exact integer powers of two are returned directly with ldexp().
 */

static int
exp2_fallback(const double x, gsl_sf_result * result)
{
  const double ln2_hi = 6.93147180559945286227e-01;
  const double ln2_lo = 2.31904681384629955842e-17;
  const double n_d = floor(x + 0.5);
  double r, y;
  int n;

  if (n_d >= 1024.0)
    {
      n = 1024;
      r = x - 1024.0;
    }
  else if (n_d <= -1075.0)
    {
      /* Use exponent -1074 and leave the remaining factor in r. */
      n = -1074;
      r = x + 1074.0;
    }
  else
    {
      n = (int) n_d;
      r = x - n_d;
    }

  if (r == 0.0)
    {
      result->val = ldexp(1.0, n);
      result->err = 0.0;
      return GSL_SUCCESS;
    }

  y = exp(r * ln2_hi + r * ln2_lo);
  result->val = ldexp(y, n);

  if (fabs(result->val) < GSL_DBL_MIN)
    result->err = GSL_SF_TRUE_MIN;
  else
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);

  return GSL_SUCCESS;
}
