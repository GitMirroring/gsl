/* specfunc/trig.c
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
#include <gsl/gsl_sf_trig.h>

#include "error.h"


/* sinh(x) series
 * double-precision for |x| < 1.0
 */
inline
static
int
sinh_series(const double x, double * result)
{
  const double y = x*x;
  const double c0 = 1.0/6.0;
  const double c1 = 1.0/120.0;
  const double c2 = 1.0/5040.0;
  const double c3 = 1.0/362880.0;
  const double c4 = 1.0/39916800.0;
  const double c5 = 1.0/6227020800.0;
  const double c6 = 1.0/1307674368000.0;
  const double c7 = 1.0/355687428096000.0;
  *result = x*(1.0 + y*(c0+y*(c1+y*(c2+y*(c3+y*(c4+y*(c5+y*(c6+y*c7))))))));
  return GSL_SUCCESS;
}


/* cosh(x)-1 series
 * double-precision for |x| < 1.0
 */
inline
static
int
cosh_m1_series(const double x, double * result)
{
  const double y = x*x;
  const double c0 = 0.5;
  const double c1 = 1.0/24.0;
  const double c2 = 1.0/720.0;
  const double c3 = 1.0/40320.0;
  const double c4 = 1.0/3628800.0;
  const double c5 = 1.0/479001600.0;
  const double c6 = 1.0/87178291200.0;
  const double c7 = 1.0/20922789888000.0;
  const double c8 = 1.0/6402373705728000.0;
  *result = y*(c0+y*(c1+y*(c2+y*(c3+y*(c4+y*(c5+y*(c6+y*(c7+y*c8))))))));
  return GSL_SUCCESS;
}


int
gsl_sf_hypot_e(const double x, const double y, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(gsl_isinf(x) || gsl_isinf(y)) {
    /* +Inf even if the other argument is NaN (C99 hypot) */
    result->val = GSL_POSINF;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(gsl_isnan(x) || gsl_isnan(y)) {
    /* the comparisons below would discard a NaN argument */
    DOMAIN_ERROR(result);
  }
  else if(x == 0.0 && y == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else {
    const double a = fabs(x);
    const double b = fabs(y);
    const double min = GSL_MIN_DBL(a,b);
    const double max = GSL_MAX_DBL(a,b);
    const double rat = min/max;
    const double root_term = sqrt(1.0 + rat*rat);

    if(max < GSL_DBL_MAX/root_term) {
      result->val = max * root_term;
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else {
      OVERFLOW_ERROR(result);
    }
  }
}


int
gsl_sf_complex_sin_e(const double zr, const double zi,
                        gsl_sf_result * szr, gsl_sf_result * szi)
{
  /* CHECK_POINTER(szr) */
  /* CHECK_POINTER(szi) */

  if(fabs(zi) < 1.0) {
    double ch_m1, sh;
    sinh_series(zi, &sh);
    cosh_m1_series(zi, &ch_m1);
    szr->val = sin(zr)*(ch_m1 + 1.0);
    szi->val = cos(zr)*sh;
    szr->err = 2.0 * GSL_DBL_EPSILON * fabs(szr->val);
    szi->err = 2.0 * GSL_DBL_EPSILON * fabs(szi->val);
    return GSL_SUCCESS;
  }
  else if(fabs(zi) < GSL_LOG_DBL_MAX) {
    double ex = exp(zi);
    double ch = 0.5*(ex+1.0/ex);
    double sh = 0.5*(ex-1.0/ex);
    szr->val = sin(zr)*ch;
    szi->val = cos(zr)*sh;
    szr->err = 2.0 * GSL_DBL_EPSILON * fabs(szr->val);
    szi->err = 2.0 * GSL_DBL_EPSILON * fabs(szi->val);
    return GSL_SUCCESS;
  }
  else {
    OVERFLOW_ERROR_2(szr, szi);
  }
}


int
gsl_sf_complex_cos_e(const double zr, const double zi,
                        gsl_sf_result * czr, gsl_sf_result * czi)
{
  /* CHECK_POINTER(czr) */
  /* CHECK_POINTER(czi) */

  if(fabs(zi) < 1.0) {
    double ch_m1, sh;
    sinh_series(zi, &sh);
    cosh_m1_series(zi, &ch_m1);
    czr->val =  cos(zr)*(ch_m1 + 1.0);
    czi->val = -sin(zr)*sh;
    czr->err = 2.0 * GSL_DBL_EPSILON * fabs(czr->val);
    czi->err = 2.0 * GSL_DBL_EPSILON * fabs(czi->val);
    return GSL_SUCCESS;
  }
  else if(fabs(zi) < GSL_LOG_DBL_MAX) {
    double ex = exp(zi);
    double ch = 0.5*(ex+1.0/ex);
    double sh = 0.5*(ex-1.0/ex);
    czr->val =  cos(zr)*ch;
    czi->val = -sin(zr)*sh;
    czr->err = 2.0 * GSL_DBL_EPSILON * fabs(czr->val);
    czi->err = 2.0 * GSL_DBL_EPSILON * fabs(czi->val);
    return GSL_SUCCESS;
  }
  else {
    OVERFLOW_ERROR_2(czr,czi);
  }
}


int
gsl_sf_complex_logsin_e(const double zr, const double zi,
                           gsl_sf_result * lszr, gsl_sf_result * lszi)
{
  /* CHECK_POINTER(lszr) */
  /* CHECK_POINTER(lszi) */

  if(zi > 60.0) {
    lszr->val = -M_LN2 + zi;
    lszi->val =  0.5*M_PI - zr;
    lszr->err = 2.0 * GSL_DBL_EPSILON * fabs(lszr->val);
    lszi->err = 2.0 * GSL_DBL_EPSILON * fabs(lszi->val);
  }
  else if(zi < -60.0) {
    lszr->val = -M_LN2 - zi;
    lszi->val = -0.5*M_PI + zr;
    lszr->err = 2.0 * GSL_DBL_EPSILON * fabs(lszr->val);
    lszi->err = 2.0 * GSL_DBL_EPSILON * fabs(lszi->val);
  }
  else {
    gsl_sf_result sin_r, sin_i;
    int status;
    gsl_sf_complex_sin_e(zr, zi, &sin_r, &sin_i); /* ok by construction */
    status = gsl_sf_complex_log_e(sin_r.val, sin_i.val, lszr, lszi);
    if(status == GSL_EDOM) {
      DOMAIN_ERROR_2(lszr, lszi);
    }
  }
  return gsl_sf_angle_restrict_symm_e(&(lszi->val));
}


int
gsl_sf_lnsinh_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(fabs(x) < 1.0) {
    double eps;
    sinh_series(x, &eps);
    result->val = log(eps);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < -0.5*GSL_LOG_DBL_EPSILON) {
    result->val = x + log(0.5*(1.0 - exp(-2.0*x)));
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    result->val = -M_LN2 + x;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int gsl_sf_lncosh_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(fabs(x) < 1.0) {
    double eps;
    cosh_m1_series(x, &eps);
    return gsl_sf_log1p_e(eps, result);
  }
  else if(fabs(x) < -0.5*GSL_LOG_DBL_EPSILON) {
    result->val = fabs(x) + log(0.5*(1.0 + exp(-2.0*fabs(x))));
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    result->val = -M_LN2 + fabs(x);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}

int
gsl_sf_polar_to_rect(const double r, const double theta,
                          gsl_sf_result * x, gsl_sf_result * y)
{
  double t   = theta;
  int status = gsl_sf_angle_restrict_symm_e(&t);
  double c = cos(t);
  double s = sin(t);
  x->val = r * cos(t);
  y->val = r * sin(t);
  x->err  = r * fabs(s * GSL_DBL_EPSILON * t);
  x->err += 2.0 * GSL_DBL_EPSILON * fabs(x->val);
  y->err  = r * fabs(c * GSL_DBL_EPSILON * t);
  y->err += 2.0 * GSL_DBL_EPSILON * fabs(y->val);
  return status;
}


int
gsl_sf_rect_to_polar(const double x, const double y,
                          gsl_sf_result * r, gsl_sf_result * theta)
{
  int stat_h = gsl_sf_hypot_e(x, y, r);
  if(r->val > 0.0) {
    theta->val = atan2(y, x);
    theta->err = 2.0 * GSL_DBL_EPSILON * fabs(theta->val);
    return stat_h;
  }
  else {
    DOMAIN_ERROR(theta);
  }
}


int gsl_sf_angle_restrict_symm_err_e(const double theta, gsl_sf_result * result)
{
  /* synthetic extended precision constants */
  const double P1 = 4 * 7.8539812564849853515625e-01;
  const double P2 = 4 * 3.7748947079307981766760e-08;
  const double P3 = 4 * 2.6951514290790594840552e-15;
  const double TwoPi = 2*(P1 + P2 + P3);

  const double y = GSL_SIGN(theta) * 2 * floor(fabs(theta)/TwoPi);
  double r = ((theta - y*P1) - y*P2) - y*P3;

  if(r >  M_PI) { r = (((r-2*P1)-2*P2)-2*P3); }  /* r-TwoPi */
  else if (r < -M_PI) r = (((r+2*P1)+2*P2)+2*P3); /* r+TwoPi */

  result->val = r;

  if(fabs(theta) > 0.0625/GSL_DBL_EPSILON) {
    result->val = GSL_NAN;
    result->err = GSL_NAN;
    GSL_ERROR ("error", GSL_ELOSS);
  }
  else if(fabs(theta) > 0.0625/GSL_SQRT_DBL_EPSILON) {
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val - theta);
    return GSL_SUCCESS;
  }
  else {
    double delta = fabs(result->val - theta);
    result->err = 2.0 * GSL_DBL_EPSILON * ((delta < M_PI) ? delta : M_PI);
    return GSL_SUCCESS;
  }
}


int gsl_sf_angle_restrict_pos_err_e(const double theta, gsl_sf_result * result)
{
  /* synthetic extended precision constants */
  const double P1 = 4 * 7.85398125648498535156e-01;
  const double P2 = 4 * 3.77489470793079817668e-08;
  const double P3 = 4 * 2.69515142907905952645e-15;
  const double TwoPi = 2*(P1 + P2 + P3);

  const double y = 2*floor(theta/TwoPi);

  double r = ((theta - y*P1) - y*P2) - y*P3;

  if(r > TwoPi) {r = (((r-2*P1)-2*P2)-2*P3); }  /* r-TwoPi */
  else if (r < 0) { /* may happen due to FP rounding */
    r = (((r+2*P1)+2*P2)+2*P3); /* r+TwoPi */
  }

  result->val = r;

  if(fabs(theta) > 0.0625/GSL_DBL_EPSILON) {
    result->val = GSL_NAN;
    result->err = fabs(result->val);
    GSL_ERROR ("error", GSL_ELOSS);
  }
  else if(fabs(theta) > 0.0625/GSL_SQRT_DBL_EPSILON) {
    result->err = GSL_DBL_EPSILON * fabs(result->val - theta);
    return GSL_SUCCESS;
  }
  else {
    double delta = fabs(result->val - theta);
    result->err = 2.0 * GSL_DBL_EPSILON * ((delta < M_PI) ? delta : M_PI);
    return GSL_SUCCESS;
  }
}


int gsl_sf_angle_restrict_symm_e(double * theta)
{
  gsl_sf_result r;
  int stat = gsl_sf_angle_restrict_symm_err_e(*theta, &r);
  *theta = r.val;
  return stat;
}


int gsl_sf_angle_restrict_pos_e(double * theta)
{
  gsl_sf_result r;
  int stat = gsl_sf_angle_restrict_pos_err_e(*theta, &r);
  *theta = r.val;
  return stat;
}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

double gsl_sf_hypot(const double x, const double y)
{
  EVAL_RESULT(gsl_sf_hypot_e(x, y, &result));
}

double gsl_sf_lnsinh(const double x)
{
  EVAL_RESULT(gsl_sf_lnsinh_e(x, &result));
}

double gsl_sf_lncosh(const double x)
{
  EVAL_RESULT(gsl_sf_lncosh_e(x, &result));
}

double gsl_sf_angle_restrict_symm(const double theta)
{
  double result = theta;
  EVAL_DOUBLE(gsl_sf_angle_restrict_symm_e(&result));
}

double gsl_sf_angle_restrict_pos(const double theta)
{
  double result = theta;
  EVAL_DOUBLE(gsl_sf_angle_restrict_pos_e(&result));
}
