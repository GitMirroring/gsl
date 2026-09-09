/* specfunc/log/log1pmx.c
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
#include "chebyshev.h"
#include "cheb_eval.c"

static int kernel_cheb(const double x, gsl_sf_result * result);
static int kernel_log1pmx(const double x, gsl_sf_result * result);
static int kernel_logmxp1(const double x, gsl_sf_result * result);
static int kernel_poly(const double x, gsl_sf_result * result);

#define LOGMXP1_C1    3.55829253011726237e-1 /* log(8/3) - 5/8 */
#define LOGMXP1_C2    1.93147180559945309e-1 /* log(2) - 1/2 */

/*
 * Compute
 *
 *   log1pmx(x) = log(1+x) - x,    x >= -1.
 *
 * Direct evaluation loses accuracy near x = 0 because log(1+x)
 * and x nearly cancel.  Different kernels are therefore used in
 * regions where they give the best accuracy:
 *
 *   -1 < x < -0.5       logmxp1(1+x) transformation
 *   -0.5 <= x <= -0.425 Chebyshev approximation
 *   -0.425 < x < 0.7    transformed log1pmx polynomial
 *    0.7 <= x <= 1      dedicated degree-16 polynomial
 *    x > 1              log1p(x) - x
 */
int
gsl_sf_log1pmx_e(const double x, gsl_sf_result * result)
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
  else if (-0.425 < x && x < 0.7)
    {
      return kernel_log1pmx(x, result);
    }
  else if (-0.5 <= x && x <= -0.425)
    {
      return kernel_cheb(x, result);
    }
  else if (x < -0.5)
    {
      return kernel_logmxp1(1.0 + x, result);
    }
  else if (x <= 1.0)
    {
      return kernel_poly(x, result);
    }
  else
    {
      const double lterm = log1p(x);
      result->val = lterm - x;
      result->err = 2.0 * GSL_DBL_EPSILON * (fabs(lterm) + fabs(x));
      return GSL_SUCCESS;
    }
}

double gsl_sf_log1pmx(const double x) { EVAL_RESULT(gsl_sf_log1pmx_e(x, &result)); }

/********************************************
 * INTERNAL ROUTINES                        *
 ********************************************/

/*
 * Evaluate log(1+x) - x using the original GSL Chebyshev
 * approximation to
 *
 *   (log(1+x) - x) / x^2.
 *
 * The rational transformation
 *
 *   t = (8x + 1) / (2(x + 2))
 *
 * maps the approximation interval onto [-1,1].  This kernel is
 * used on -0.5 <= x <= -0.425
 */
static int
kernel_cheb(const double x, gsl_sf_result * result)
{
  /* Existing GSL Chebyshev approximation for
   *
   *   (log(1+x) - x) / x^2,  -1/2 < x < 1/2.
   *
   * with t = (8x+1)/(2(x+2)).
   */
  static const double lopxmx_data[20] = { -1.12100231323744103373737274541,      0.19553462773379386241549597019,    -0.01467470453808083971825344956,
                                           0.00166678250474365477643629067,     -0.00018543356147700369785746902,     0.00002280154021771635036301071,
                                          -2.8031253116633521699214134172e-06,   3.5936568872522162983669541401e-07, -4.6241857041062060284381167925e-08,
                                           6.08226374594039910117102904131e-09, -8.0339824424815790302621320732e-10,  1.0751718277499375044851551587e-10,
                                          -1.4445310914224613448759230882e-11,   1.9573912180610336168921438426e-12, -2.6614436796793061741564104510e-13,
                                           3.6402634315269586532158344584e-14,  -4.9937495922755006545809120531e-15,  6.8802890218846809524646902703e-16,
                                          -9.5034129794804273611403251480e-17,   1.3170135013050997157326965813e-17 };
  const cheb_series lopxmx_cs = { (double *) lopxmx_data, 19, -1, 1, 9 };
  const double t = 0.5 * (8.0 * x + 1.0) / (x + 2.0);
  const double x2 = x * x;
  gsl_sf_result r;

  cheb_eval_e(&lopxmx_cs, t, &r);

  result->val = x2 * r.val;
  result->err = x2 * r.err + 2.0 * GSL_DBL_EPSILON * fabs(result->val);

  return GSL_SUCCESS;
}

/*
 * Evaluate log(1+x) - x using the Float64 polynomial kernel from
 * JuliaStats/LogExpFunctions.
 *
 * With
 *
 *   r = x / (x + 2),  t = r^2,
 *
 * a degree-8 polynomial in t approximates the transformed remainder.
 * The polynomial is evaluated by Horner's method for portability.
 *
 * This kernel is used on -0.425 < x < 0.7.
 */
static int
kernel_log1pmx(const double x, gsl_sf_result * result)
{
  static const double p[9] = { 0.6666666666666669, 0.3999999999997768,  0.2857142857784595,
                               0.2222222142048249, 0.18181870670924566, 0.15382646727504887,
                               0.1337701340211177, 0.11201972567415432, 0.143418239946679 };
  const double r = x / (x + 2.0);
  const double t = r * r;
  const double hxsq = x * x / 2.0;
  double w = p[0] + t * (p[1] + t * (p[2] + t * (p[3] + t * (p[4] + t * (p[5] + t * (p[6] + t * (p[7] + t * p[8])))))));

  result->val = r * (w * t + hxsq) - hxsq;
  result->err = 4.0 * GSL_DBL_EPSILON * fabs(result->val);

  return GSL_SUCCESS;
}

/*
 * Evaluate
 *
 *   logmxp1(x) = log(x) - x + 1
 *
 * using the Float64 formulation from JuliaStats/LogExpFunctions.
 * For x near 3/8 and 1/2, the argument is transformed so that the
 * accurate log1pmx kernel can be reused.
 *
 * In gsl_sf_log1pmx_e this kernel is called with x = 1 + z for
 * -1 < z < -0.5, so only 0 < x < 0.5 is required.
 */
static int
kernel_logmxp1(const double x, gsl_sf_result * result)
{
  int status;
  gsl_sf_result r;

  if (x <= 0.3)
    {
      const double lterm = log(x);

      result->val = (lterm + 1.0) - x;
      result->err = 2.0 * GSL_DBL_EPSILON * (fabs(lterm) + 1.0 + fabs(x));
    }
  else if (x <= 0.4)
    {
      const double u = (x - 0.375) / 0.375;

      status = kernel_log1pmx(u, &r);
      if (status)
        return status;

      result->val = r.val - LOGMXP1_C1 + 0.625 * u;
      result->err = r.err + 2.0 * GSL_DBL_EPSILON * (fabs(r.val) + LOGMXP1_C1 + fabs(0.625 * u));
    }
  else
    {
      const double u = 2.0 * (x - 0.5);

      status = kernel_log1pmx(u, &r);
      if (status)
        return status;

      result->val = r.val - LOGMXP1_C2 + 0.5 * u;
      result->err = r.err + 2.0 * GSL_DBL_EPSILON * (fabs(r.val) + LOGMXP1_C2 + fabs(0.5 * u));
    }

  return GSL_SUCCESS;
}

/*
 * Evaluate log(1+x) - x on 0.7 <= x <= 1 using a dedicated
 * degree-16 polynomial.
 *
 * The centered variable
 *
 *   z = (x - 0.85) / 0.15
 *
 * maps the interval [0.7,1] onto [-1,1].  The polynomial is
 * evaluated by Horner's method.
 */
static int
kernel_poly(const double x, gsl_sf_result * result)
{
  static const double c[17] = { -2.2356519305974621e-19, 2.9405863435785961e-18, -3.7899601098817681e-17,
                                 5.0338661774506289e-16, -6.7274999406176179e-15, 9.0515450170917596e-14,
                                -1.227991342100528e-12, 1.6828029505940208e-11, -2.3348891025953613e-10,
                                 3.2910817827043572e-09, -4.735501009531037e-08, 7.0085414941059387e-07,
                                -1.0804834803413368e-05, 1.7767950565613093e-04, -3.2870708546384221e-03,
                                -6.8918918918918923e-02, -2.3481436090976654e-01 };
  const double z = (x - 0.85) / 0.15;
  double p = c[0];
  int i;

  for (i = 1; i < 17; ++i)
    p = p * z + c[i];

  result->val = p;
  result->err = 3.0 * GSL_DBL_EPSILON * fabs(result->val);

  return GSL_SUCCESS;
}
