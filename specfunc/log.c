/* specfunc/log.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
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

/* Author:  G. Jungman */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_log.h>

#include "error.h"

#include "chebyshev.h"
#include "cheb_eval.c"

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_log_e(const double x, gsl_sf_result * result)
{
  if(x <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else {
    result->val = log(x);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}

int
gsl_sf_log_abs_e(const double x, gsl_sf_result * result)
{
  if(x == 0.0) {
    DOMAIN_ERROR(result);
  }
  else {
    result->val = log(fabs(x));
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}

int
gsl_sf_complex_log_e(const double zr, const double zi, gsl_sf_result * lnr, gsl_sf_result * theta)
{
  if(zr != 0.0 || zi != 0.0) {
    const double ax = fabs(zr);
    const double ay = fabs(zi);
    const double min = GSL_MIN(ax, ay);
    const double max = GSL_MAX(ax, ay);
    lnr->val = log(max) + 0.5 * log(1.0 + (min/max)*(min/max));
    lnr->err = 2.0 * GSL_DBL_EPSILON * fabs(lnr->val);
    theta->val = atan2(zi, zr);
    theta->err = GSL_DBL_EPSILON * fabs(lnr->val);
    return GSL_SUCCESS;
  }
  else {
    DOMAIN_ERROR_2(lnr, theta);
  }
}

#ifndef GSL_DISABLE_DEPRECATED

int
gsl_sf_log_1plusx_e(const double x, gsl_sf_result * result)
{
  return gsl_sf_log1p_e(x, result);
}

int
gsl_sf_log_1plusx_mx_e(const double x, gsl_sf_result * result)
{
  return gsl_sf_log1pmx_e(x, result);
}

#endif /* !GSL_DISABLE_DEPRECATED */

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

double gsl_sf_log(const double x)
{
  EVAL_RESULT(gsl_sf_log_e(x, &result));
}

double gsl_sf_log_abs(const double x)
{
  EVAL_RESULT(gsl_sf_log_abs_e(x, &result));
}

#ifndef GSL_DISABLE_DEPRECATED

double gsl_sf_log_1plusx(const double x)
{
  EVAL_RESULT(gsl_sf_log_1plusx_e(x, &result));
}

double gsl_sf_log_1plusx_mx(const double x)
{
  EVAL_RESULT(gsl_sf_log_1plusx_mx_e(x, &result));
}

#endif /* !GSL_DISABLE_DEPRECATED */
