/* specfunc/exp/exp2.c
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
#include "eval.h"

#ifndef HAVE_EXP2
#include "exp2_fallback.c"
#endif

/*
 * Compute the base-2 exponential
 *
 *   exp2(x) = 2^x.
 *
 * The system exp2() function is used when available.  On systems without
 * exp2(), a portable C89 fallback is used
 */

int
gsl_sf_exp2_e(const double x, gsl_sf_result * result)
{
  /* check for overflow */
  if (x >= 1024.0)
    {
      OVERFLOW_ERROR(result);
    }
  else if (x <= -1075.0)
    {
      result->val = 0.0;
      result->err = GSL_SF_TRUE_MIN;
      return GSL_SUCCESS;
    }
  else
    {
#ifdef HAVE_EXP2

      result->val = exp2(x);

      if (fabs(result->val) < GSL_DBL_MIN)
        result->err = GSL_SF_TRUE_MIN;
      else
        result->err = 3.0 * GSL_DBL_EPSILON * fabs(result->val);

      return GSL_SUCCESS;

#else

      return exp2_fallback(x, result);

#endif
    }
}

double gsl_sf_exp2(const double x) { EVAL_RESULT(gsl_sf_exp2_e(x, &result)); }
