/* specfunc/gamma/beta.c
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

/*
 * References: the algorithms are original to GSL, using the reflection
 * formula, Stirling's series (NIST DLMF 5.11.1) and the exact
 * argument-reduction techniques of pow_exp.c; the reflection formulas
 * follow lnbeta.c.  No external code was used.
 */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_dd.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_sincos_pi.h>

#include "constants.h"
#include "error.h"
#include "eval.h"
#include "stirling.c"
#include "pow_exp.c"

/*
 * relative error bounds (in units of GSL_DBL_EPSILON) of the
 * positive-argument kernels, from their rounding error counts
 */
#define BETA_POS_ERR   6.0
#define BETA_NEG_ERR  12.0

static double beta_pos(const double a0, const double b0);
static double beta_pos_small(const double a, const double b);
static double beta_pos_mixed(const double a, const double b);
static double beta_pos_large(const double a, const double b);
static double beta_reduce(const double x, const gsl_dd * xy, double * x0);
static double beta_pos_dd(const gsl_dd x, const gsl_dd y);
static double beta_sin_pi(const gsl_dd x);
static int beta_isnonposint(const double x);

/*
 * Beta function
 *
 *   B(a,b) = Gamma(a) Gamma(b) / Gamma(a+b).
 *
 * B has poles where a or b is a non-positive integer (GSL_EDOM) and is
 * zero where only a+b is one.  B is formed directly as a product, not
 * as exp(log B), whose relative error grows with |log B|.
 *
 * Positive arguments are handled by beta_pos().  Otherwise, with a <= b,
 * the reflection formula Gamma(x) = pi / (sin(pi x) Gamma(1-x)) reduces
 * B to a positive-argument Beta function and elementary factors,
 * s(x) = sin(pi x):
 *
 *   a < 0 < b, a+b > 0:  B(a,b) = pi / (b s(a) B(1-a, a+b))
 *   a < 0 < b, a+b < 0:  B(a,b) = s(a+b) / s(a) B(b, 1-a-b)
 *   a, b < 0:            B(a,b) = pi s(a+b) / (s(a) s(b) (1-a-b) B(1-a, 1-b))
 *
 * The sums a+b, 1-a, 1-a-b, 1-b are formed exactly, and their low parts
 * are carried to first order.
 */
int
gsl_sf_beta_e(const double a0, const double b0, gsl_sf_result * result)
{
  const double a = GSL_MIN(a0, b0);
  const double b = GSL_MAX(a0, b0);
  double v, rel;

  if (beta_isnonposint(a) || beta_isnonposint(b))
    {
      DOMAIN_ERROR(result);
    }
  else if (a > 0.0)
    {
      v = beta_pos(a, b);
      rel = BETA_POS_ERR;
    }
  else
    {
      const gsl_dd one = { 1.0, 0.0 };
      const gsl_dd s = gsl_dd_two_sum(a, b);
      gsl_dd ad, bd;

      if (s.lo == 0.0 && beta_isnonposint(s.hi))
        {
          /* 1/Gamma(a+b) = 0 */
          result->val = 0.0;
          result->err = 0.0;
          return GSL_SUCCESS;
        }

      ad.hi = a;
      ad.lo = 0.0;
      bd.hi = b;
      bd.lo = 0.0;

      if (b > 0.0 && (s.hi > 0.0 || (s.hi == 0.0 && s.lo > 0.0)))
        v = M_PI / ((b * beta_sin_pi(ad)) * beta_pos_dd(gsl_dd_sub(one, ad), s));
      else if (b > 0.0)
        v = (beta_sin_pi(s) / beta_sin_pi(ad)) * beta_pos_dd(bd, gsl_dd_sub(one, s));
      else
        {
          const gsl_dd y = gsl_dd_sub(one, s);
          const double d = (beta_sin_pi(ad) * beta_sin_pi(bd)) * (y.hi + y.lo);

          v = (M_PI * beta_sin_pi(s)) / (d * beta_pos_dd(gsl_dd_sub(one, ad), gsl_dd_sub(one, bd)));
        }

      rel = BETA_NEG_ERR;
    }

  if (!gsl_finite(v))
    {
      OVERFLOW_ERROR(result);
    }
  else if (fabs(v) < GSL_DBL_MIN)
    {
      UNDERFLOW_ERROR(result);
    }

  result->val = v;
  result->err = rel * GSL_DBL_EPSILON * fabs(v);

  return GSL_SUCCESS;
}

double gsl_sf_beta(const double a, const double b) { EVAL_RESULT(gsl_sf_beta_e(a, b, &result)); }

/*********************************************************
 * INTERNAL ROUTINES: POSITIVE ARGUMENTS                 *
 *********************************************************/

/*
 * B(a,b) for a, b > 0.  The result may overflow to inf or underflow
 * to zero; the caller signals that.
 */
static double
beta_pos(const double a0, const double b0)
{
  const double a = GSL_MIN(a0, b0);
  const double b = GSL_MAX(a0, b0);

  if (a >= 10.0)
    return beta_pos_large(a, b);
  else if (b >= 10.0)
    return beta_pos_mixed(a, b);
  else
    return beta_pos_small(a, b);
}

/*
 * Reduce x > 0 to x0 = x - n in (0,1] (exactly) with the recurrence
 *
 *   B(x, y) = (x-1)/(x+y-1) B(x-1, y),
 *
 * where xy = x + y is given exactly.  Returns the product of the
 * factors, formed in double-double.
 */
static double
beta_reduce(const double x, const gsl_dd * xy, double * x0)
{
  const int n = (x > 1.0) ? (int) ceil(x) - 1 : 0;
  gsl_dd f = { 1.0, 0.0 };
  int k;

  for (k = 1; k <= n; ++k)
    {
      gsl_dd num;

      num.hi = x - (double) k;      /* exact */
      num.lo = 0.0;
      f = gsl_dd_mul(f, gsl_dd_div(num, gsl_dd_add_d(*xy, -(double) k)));
    }

  *x0 = x - (double) n;            /* exact (Sterbenz) */

  return f.hi + f.lo;
}

/*
 * B(a,b) for 0 < a <= b < 10.  Both arguments are reduced to (0,1]
 * (beta_reduce()), and
 *
 *   B(a0,b0) = (1/a0 + 1/b0) Gamma(1+a0) Gamma(1+b0) / Gamma(1+a0+b0)
 *            = (1/a0 + 1/b0) exp(lngamma1p(a0) + lngamma1p(b0) - lngamma1p(a0+b0)),
 *
 * where the exponent is small (|.| < 1) and a0 + b0 is exact; the
 * rounding of its leading part is corrected with
 * d/dx lngamma1p(x) = psi(1+x).
 */
static double
beta_pos_small(const double a, const double b)
{
  const gsl_dd ab = gsl_dd_two_sum(a, b);
  double a0, b0, fa, fb, e;
  gsl_dd ab0, q;

  fa = beta_reduce(a, &ab, &a0);
  ab0 = gsl_dd_two_sum(a0, b);
  fb = beta_reduce(b, &ab0, &b0);

  ab0 = gsl_dd_two_sum(a0, b0);
  e = gsl_sf_lngamma1p(a0) + gsl_sf_lngamma1p(b0) - gsl_sf_lngamma1p(ab0.hi);
  if (ab0.lo != 0.0)
    e -= gsl_sf_psi(1.0 + ab0.hi) * ab0.lo;

  {
    /* (a0+b0)/(a0 b0) = 1/a0 + 1/b0: no underflow of a0 b0 for tiny arguments */
    const gsl_dd one = { 1.0, 0.0 };
    gsl_dd ad, bd;

    ad.hi = a0;
    ad.lo = 0.0;
    bd.hi = b0;
    bd.lo = 0.0;
    q = gsl_dd_add(gsl_dd_div(one, ad), gsl_dd_div(one, bd));
  }

  return ((fa * fb) * (q.hi + q.lo)) * exp(e);
}

/*
 * B(a,b) for 0 < a < 10 <= b.  a is reduced to (0,1] (beta_reduce()),
 * and with s = a0 + b, t = a0/s and Stirling's formula,
 *
 *   B(a0,b) = Gamma(a0) Gamma(b)/Gamma(s)
 *           = exp(lngamma1p(a0)) / a0 * s^(-a0) * exp(E),
 *   E = (b - 1/2) log1pmx(-t) + a0 (a0 + 1/2)/s + Delta(b) - Delta(s),
 *
 * where E is small (the large terms (b-1/2) log(1-t) and a0 cancel
 * analytically).  s = s.hi + s.lo exactly; the low part enters
 * s^(-a0) as the factor exp(-a0 s.lo/s.hi).
 */
static double
beta_pos_mixed(const double a, const double b)
{
  const gsl_dd ab = gsl_dd_two_sum(a, b);
  double a0, f, t, e;
  gsl_dd s;

  f = beta_reduce(a, &ab, &a0);
  s = gsl_dd_two_sum(a0, b);
  t = a0 / s.hi;

  e = (b - 0.5) * gsl_sf_log1pmx(-t) + a0 * (a0 + 0.5) / s.hi +
      (gamma_stirling_delta(b) - gamma_stirling_delta(s.hi)) - a0 * (s.lo / s.hi);

  return (f / a0) * (pow(s.hi, -a0) * exp(gsl_sf_lngamma1p(a0) + e));
}

/*
 * B(a,b) for 10 <= a <= b, from Stirling's formula:
 *
 *   B(a,b) = sqrt(2 pi) sqrt(1/a + 1/b) q1^a q2^b exp(E),
 *   q1 = a/s, q2 = b/s, s = a + b,
 *   E = Delta(a) + Delta(b) - Delta(s) + (corrections),
 *
 * where q1^a and q2^b are formed by pow() from q = fl(a/s.hi) with the
 * exact remainder rho = a - q s.hi (gamma_div_rem()), so that
 * (a/s)^a = q^a exp(a (rho/(q s.hi) - s.lo/s.hi)); the large exponents
 * are never rounded.  Near underflow of the powers, fall back to
 * exp(log B) (the result then underflows or is barely normal).
 */
static double
beta_pos_large(const double a, const double b)
{
  const gsl_dd s = gsl_dd_two_sum(a, b);
  const double q1 = a / s.hi;
  const double q2 = b / s.hi;
  const double m1 = pow(q1, a);
  const double m2 = pow(q2, b);
  double e;

  if (m1 < 1.0e-290 || m2 < 1.0e-290)
    return exp(gsl_sf_lnbeta(a, b));

  e = a * (gamma_div_rem(a, s.hi, q1) / (q1 * s.hi)) +
      b * (gamma_div_rem(b, s.hi, q2) / (q2 * s.hi)) - s.lo +
      gamma_stirling_delta(a) + (gamma_stirling_delta(b) - gamma_stirling_delta(s.hi));

  return ((GSL_SF_SQRT_2PI * sqrt(1.0 / a + 1.0 / b)) * (m1 * m2)) * exp(e);
}

/*********************************************************
 * INTERNAL ROUTINES: NEGATIVE ARGUMENTS                 *
 *********************************************************/

/* return 1 if x is zero or a negative integer */
static int
beta_isnonposint(const double x)
{
  return (x <= 0.0 && x == floor(x));
}

/*
 * sin(pi x) for x = x.hi + x.lo, using
 * sin(pi (x.hi + x.lo)) = sin(pi x.hi) + pi x.lo cos(pi x.hi) to first
 * order in the low part, which matters near the zeros of sin.
 */
static double
beta_sin_pi(const gsl_dd x)
{
  double v = gsl_sf_sin_pi(x.hi);

  if (x.lo != 0.0)
    v += M_PI * x.lo * gsl_sf_cos_pi(x.hi);

  return v;
}

/*
 * B(x,y) for positive x = x.hi + x.lo, y = y.hi + y.lo: beta_pos() at
 * the leading parts, corrected to first order for the low parts,
 *
 *   d log B / dx = psi(x) - psi(x+y),   d log B / dy = psi(y) - psi(x+y).
 */
static double
beta_pos_dd(const gsl_dd x, const gsl_dd y)
{
  const double v = beta_pos(x.hi, y.hi);

  if (x.lo != 0.0 || y.lo != 0.0)
    {
      const double pxy = gsl_sf_psi(x.hi + y.hi);
      const double c = x.lo * (gsl_sf_psi(x.hi) - pxy) + y.lo * (gsl_sf_psi(y.hi) - pxy);

      return v * (1.0 + c);
    }

  return v;
}
