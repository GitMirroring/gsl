/* specfunc/beta.c
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
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>

#include "error.h"

static double
isnegint (const double x) 
{
  return (x < 0) && (x == floor(x));
}

int
gsl_sf_beta_e(const double x, const double y, gsl_sf_result * result)
{
  if((x > 0 && y > 0) && x < 50.0 && y < 50.0) {
    /* Handle the easy case */
    gsl_sf_result gx, gy, gxy;
    gsl_sf_gamma_e(x, &gx);
    gsl_sf_gamma_e(y, &gy);
    gsl_sf_gamma_e(x+y, &gxy);
    result->val  = (gx.val*gy.val)/gxy.val;
    result->err  = gx.err * fabs(gy.val/gxy.val);
    result->err += gy.err * fabs(gx.val/gxy.val);
    result->err += fabs((gx.val*gy.val)/(gxy.val*gxy.val)) * gxy.err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if (isnegint(x) || isnegint(y)) {
    DOMAIN_ERROR(result);
  } else if (isnegint(x+y)) {  /* infinity in the denominator */
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  } else {
    gsl_sf_result lb;
    double sgn;
    int stat_lb = gsl_sf_lnbeta_sgn_e(x, y, &lb, &sgn);
    if(stat_lb == GSL_SUCCESS) {
      int status = gsl_sf_exp_err_e(lb.val, lb.err, result);
      result->val *= sgn;
      return status;
    }
    else {
      result->val = 0.0;
      result->err = 0.0;
      return stat_lb;
    }
  }
}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

double gsl_sf_beta(const double x, const double y)
{
  EVAL_RESULT(gsl_sf_beta_e(x, y, &result));
}
