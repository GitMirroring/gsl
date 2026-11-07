/* specfunc/trig/sinc.c
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
#include <gsl/gsl_sf_trig.h>

#include "constants.h"
#include "eval.h"

static int sinc_small(const double x, gsl_sf_result * result);
static int h_small(const double r, gsl_sf_result * result);
static void reduce_nearest(const double a, double *r, int *sign);

/*
gsl_sf_sinc_e()
  Compute the normalized sinc function

    sinc(x) = sin(pi*x)/(pi*x),

  with sinc(0) = 1, together with an estimate of the absolute
  numerical error.

  For |x| <= 0.5, a polynomial approximation centered at the
  origin is used.  For larger finite arguments, x is reduced to

    |x| = n + r,  |r| <= 0.5,

  and the result is evaluated from

    sinc(x) = (-1)^n h(r)/|x|,

  where h(r) = sin(pi*r)/pi.  Integer arguments are therefore
  represented by exact zeros.

Inputs: x      - real argument
        result - (output) function value and absolute error estimate

Return: GSL_SUCCESS on success
*/

int
gsl_sf_sinc_e(const double x, gsl_sf_result * result)
{
  int status;
  const double a = fabs(x);
  gsl_sf_result s;
  double r;
  double y;
  int sign;

  if (a <= 0.5)
    return sinc_small(x, result);

  if (a >= GSL_SF_1_DBL_EPSILON)
    {
      result->val = 0.0;
      result->err = 0.0;
      return GSL_SUCCESS;
    }

  reduce_nearest(a, &r, &sign);

  if (r == 0.0)
    {
      result->val = 0.0;
      result->err = 0.0;
      return GSL_SUCCESS;
    }

  status = h_small(r, &s);
  if (status)
    return status;

  y = s.val / a;

  result->val = (sign < 0) ? -y : y;
  result->err = s.err / a + GSL_DBL_EPSILON * fabs(result->val);

  return GSL_SUCCESS;
}

double gsl_sf_sinc(const double x) { EVAL_RESULT(gsl_sf_sinc_e(x, &result)); }

/**********************************
 * INTERNAL ROUTINES              *
 **********************************/

/*
sinc_small()
  Compute sinc(x) for |x| <= 0.5 using

    sinc(x) = 1 + z P(z),  z = x^2,

  where P is a degree-7 polynomial optimized for this interval.

  The form preserves sinc(0) = 1 exactly and avoids loss of
  accuracy from evaluating sin(pi*x)/(pi*x) near the origin.

Inputs: x      - reduced argument, |x| <= 0.5
        result - (output) function value and absolute error estimate

Return: GSL_SUCCESS
*/

static int
sinc_small(const double x, gsl_sf_result * result)
{
  static const double P[8] = { -1.64493406684822643314990746814573509233,
                                0.8117424252833519424346725555778806101073,
                               -0.1907518241219412791488859925125754376453,
                                0.02614784781307880722803542791312872550967,
                               -0.002346080963492399610245885136837096952275,
                                0.0001484281781669361754022771515825064989143,
                               -0.000006972960555564220376450033245902340235954,
                                0.0000002459180384815939379991674853390642302522 };
  const double z = x*x;
  const double p =  P[0] + z*(P[1] + z*(P[2] + z*(P[3] + z*(P[4] + z*(P[5] + z*(P[6] + z*P[7]))))));

  result->val = 1.0 + z*p;
  result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);

  return GSL_SUCCESS;
}

/*
h_small()
  Compute

    h(r) = sin(pi*r)/pi

  for |r| <= 0.5 using

    h(r) = r + r*z P(z),  z = r^2,

  where P is a degree-7 polynomial optimized independently for
  this kernel.

  The evaluation is intentionally left in the unfactored form
  r + (r*z)*P(z), which gives lower binary64 error than factoring
  out r.

Inputs: r      - reduced argument, |r| <= 0.5
        result - (output) function value and absolute error estimate

Return: GSL_SUCCESS
*/

static int
h_small(const double r, gsl_sf_result * result)
{
  static const double P[8] = { -1.644934066848226433159663043874384515665,
                                0.8117424252833519462583650184684866440375,
                               -0.1907518241219415238415438757794380326119,
                                0.02614784781308467941037526869529110099555,
                               -0.002346080963559506230814768229248132056388,
                                0.0001484281785606089877304174225951365211313,
                               -0.000006972961700747292354852591125151905966227,
                                0.0000002459193472168749781163852097615637766256 };
  const double z = r*r;
  const double p =  P[0] + z*(P[1] + z*(P[2] + z*(P[3] + z*(P[4] + z*(P[5] + z*(P[6] + z*P[7]))))));

  result->val = r + (r*z)*p;
  result->err = 3.0 * GSL_DBL_EPSILON * fabs(result->val);

  return GSL_SUCCESS;
}

/*
reduce_nearest()
  Reduce a nonnegative argument to its nearest integer:

    a = n + r,  |r| <= 0.5,

  and determine the sign (-1)^n required by

    sin(pi*a) = (-1)^n sin(pi*r).

  The caller guarantees a < 2^52, so the integer part n is exactly
  representable as a double.  Only its parity is required.

Input:  a    - nonnegative argument, a < 2^52
        r    - (output) reduced argument in [-0.5, 0.5]
        sign - (output) +1 for even n, -1 for odd n
*/

static void
reduce_nearest(const double a, double *r, int *sign)
{
  double n;
  double f;

  n = floor(a);
  f = a - n;

  if (f > 0.5) {
    n += 1.0;
    f -= 1.0;
  }

  *r = f;
  *sign = (fmod(n, 2.0) != 0.0) ? -1 : 1;
}
