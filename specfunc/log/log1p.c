/* specfunc/log/log1p.c
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
#include <gsl/gsl_sf_log.h>

#include "error.h"
#include "eval.h"

#ifndef HAVE_LOG1P
#include "log1p_fallback.c"
#endif

/*
 * Compute
 *
 *   log1p(x) = log(1+x),    x >= -1.
 *
 * Direct evaluation as log(1+x) loses accuracy for small x because
 * the addition 1+x may discard significant low-order bits.  When the
 * system provides log1p(), use it directly.  Otherwise use the
 * portable C89 fallback, which compensates for the rounding error in
 * forming 1+x.
 *
 * The function returns -Inf at x = -1 and a domain error for x < -1.
 */
int
gsl_sf_log1p_e(const double x, gsl_sf_result * result)
{
  if (x < -1.0)
    {
      DOMAIN_ERROR(result);
    }
  else if (x == -1.0)
    {
      result->val = GSL_NEGINF;
      result->err = 0.0;
      return GSL_SUCCESS;
    }
  else
    {
#ifdef HAVE_LOG1P
      result->val = log1p(x);
      result->err = 3.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
#else
      return log1p_fallback(x, result);
#endif
    }
}

double gsl_sf_log1p(const double x) { EVAL_RESULT(gsl_sf_log1p_e(x, &result)); }
