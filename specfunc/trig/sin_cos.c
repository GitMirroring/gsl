/* specfunc/trig/sin_cos.c
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

#include <config.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_trig.h>

#include "error.h"

/*
 * sin(x) and cos(x) are taken from the C library, which reduces the
 * argument exactly for every finite double and is accurate to within
 * about 1 ulp on current platforms (glibc, musl, the BSDs, macOS, the
 * Microsoft UCRT).
 *
 * The error estimate is for the exact argument x, 2 eps |result| (an
 * error of up to 1 ulp with a factor 2 of margin).  The sensitivity to
 * an uncertainty dx in x is added by gsl_sf_sin_err_e() and
 * gsl_sf_cos_err_e().
 */

/*
gsl_sf_sin_e()
  Compute sin(x) with an estimate of the absolute numerical error.

Inputs: x      - real argument
        result - (output) function value and absolute error estimate

Return: GSL_SUCCESS on success
        GSL_EDOM if x is infinite or NaN (result NaN)
*/

int
gsl_sf_sin_e(const double x, gsl_sf_result * result)
{
  /* sin(+-Inf) is undefined, and NaN propagates (as in C99, which
     returns NaN for both) */
  if (!gsl_finite(x))
    {
      DOMAIN_ERROR(result);
    }

  result->val = sin(x);
  result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
  return GSL_SUCCESS;
}

/*
gsl_sf_cos_e()
  Compute cos(x) with an estimate of the absolute numerical error.

Inputs: x      - real argument
        result - (output) function value and absolute error estimate

Return: GSL_SUCCESS on success
        GSL_EDOM if x is infinite or NaN (result NaN)
*/

int
gsl_sf_cos_e(const double x, gsl_sf_result * result)
{
  if (!gsl_finite(x))
    {
      DOMAIN_ERROR(result);
    }

  result->val = cos(x);
  result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
  return GSL_SUCCESS;
}

/*
gsl_sf_sin_err_e()
  Compute sin(x) for an argument x with absolute uncertainty dx.

Inputs: x      - real argument
        dx     - absolute error of x
        result - (output) function value and absolute error estimate,
                 including the propagated error |cos(x)| dx

Return: as gsl_sf_sin_e()
*/

int
gsl_sf_sin_err_e(const double x, const double dx, gsl_sf_result * result)
{
  int stat_s = gsl_sf_sin_e(x, result);
  result->err += fabs(cos(x) * dx);
  result->err += GSL_DBL_EPSILON * fabs(result->val);
  return stat_s;
}

/*
gsl_sf_cos_err_e()
  Compute cos(x) for an argument x with absolute uncertainty dx.

Inputs: x      - real argument
        dx     - absolute error of x
        result - (output) function value and absolute error estimate,
                 including the propagated error |sin(x)| dx

Return: as gsl_sf_cos_e()
*/

int
gsl_sf_cos_err_e(const double x, const double dx, gsl_sf_result * result)
{
  int stat_c = gsl_sf_cos_e(x, result);
  result->err += fabs(sin(x) * dx);
  result->err += GSL_DBL_EPSILON * fabs(result->val);
  return stat_c;
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

double
gsl_sf_sin(const double x)
{
  EVAL_RESULT(gsl_sf_sin_e(x, &result));
}

double
gsl_sf_cos(const double x)
{
  EVAL_RESULT(gsl_sf_cos_e(x, &result));
}
