/* specfunc/log/log1p_fallback.c
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
#include <gsl/gsl_sf_log.h>

#include "constants.h"

#define LOG1P_SMALL  1.86264514923095703125e-9   /* 2^-29 */

/* C89 fallback for log(1+x). The caller must ensure x > -1 */
static int
log1p_fallback(const double x, gsl_sf_result * result)
{
  const double ax = fabs(x);

  /*
   * For |x| < 2^-54, log(1+x) rounds to x.  Up to 2^-29,
   * x - x^2/2 avoids loss from forming 1+x and is sufficiently
   * accurate for binary64.
   */
  if (ax < GSL_SF_DBL_EPSILON_4)
    {
      result->val = x;
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    }
  else if (ax < LOG1P_SMALL)
    {
      result->val = x - 0.5 * x * x;
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    }
  else
    {
      const double u = 1.0 + x;
      const double c = (u - 1.0) - x;

      result->val = log(u) - c / u;
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    }

  return GSL_SUCCESS;
}
