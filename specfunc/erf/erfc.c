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
gsl_new_sf_erfc_e(const double x, gsl_sf_result * result)
{
#ifdef HAVE_ERFC

  result->val = erfc(x);
  result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
  return GSL_SUCCESS;

#else

  return erfc_fallback(x, result);

#endif
}

double gsl_new_sf_erfc(const double x) { EVAL_RESULT(gsl_new_sf_erfc_e(x, &result)); }
