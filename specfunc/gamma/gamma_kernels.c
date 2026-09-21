/* specfunc/gamma/gamma_kernels.c
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
 * The large-argument Gamma reconstruction and recurrence-product
 * strategy are based on the glibc tgamma implementation
 * (sysdeps/ieee754/dbl-64/e_gamma_r.c).
 * The implementation here has been adapted to return a scaled
 * Gamma value, use GSL double-double arithmetic, and propagate
 * GSL error estimates.
 */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_dd.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_sincos_pi.h>

#include "constants.h"
#include "error.h"
#include "chebyshev.h"
#include "cheb_eval.c"

static int gamma1p_H(const double x, gsl_sf_result * result);
static int gamma1p_central(const double x, gsl_sf_result * result);

static int gamma_low(const double x, gsl_sf_result * result);
static void gamma_product(const double x, const double x_eps, const int n, double *prod, double *eps);
static double gamma_stirling_exp_adj(const double x_adj, const double x_eps);
static int gamma_high_scaled(const double x, gsl_sf_result * mant, int * exp2);
static int gamma_high_scaled_eps(const double x, const double x_eps_in,
                                 gsl_sf_result * mant, int * exp2);
static int gamma_2_inf_scaled(const double x, gsl_sf_result * mant, int * exp2);
static int gamma_neg_small(const double x, gsl_sf_result * result);
static int gamma_neg_recurrence(const double x, gsl_sf_result * result);
static int gamma_neg_reflection(const double x, gsl_sf_result * result);

#define GAMMA_SCALED_X_MAX 2.0000000000000000e+02

/******************************************************
 * GAMMA1P KERNELS                                    *
 ******************************************************/

/*
 * Evaluate the auxiliary function H(x) defined by
 *
 *   Gamma(1+x) = 1 + x (x-1) H(x).
 *
 * H(x) is represented as
 *
 *   H(x) = 1/[2(1+x)] + P(t)/Q(t),
 *   t = (4x-1)/3,
 *
 * where P and Q are Chebyshev series.  The special values H(0) and
 * H(1) are returned directly.
 *
 * On success, result contains H(x) and an absolute error estimate.
 */
static int
gamma1p_H(const double x, gsl_sf_result * result)
{
  int status;
  static const double num[11] = { 1.9997086976398870e-01,  6.8934535348472520e-02,
                                  4.4780262376513234e-04,  1.1880463370853996e-03,
                                 -1.0382518859652589e-04,  2.1913973239773272e-05,
                                 -8.7934426582967198e-06, -7.1890233176185573e-07,
                                  7.4625416052559938e-08, -4.0879206065796960e-09,
                                 -3.0616099707490186e-11 };
  static const double den[11] = { 2.0000000000000000e+00, -1.0475694402387846e-02,
                                 -2.6019107670532705e-03, -2.3141923498623592e-04,
                                 -7.5954756490615848e-05,  6.7796212148716871e-06,
                                 -5.3416572705901811e-05,  4.8314994520769372e-06,
                                  4.7079768303290990e-07, -7.3451459845921029e-08,
                                  2.5236078802005865e-09 };
  const cheb_series num_cs = { (double *) num, 10, -1.0, 1.0, 10 };
  const cheb_series den_cs = { (double *) den, 10, -1.0, 1.0, 10 };
  double t;
  gsl_sf_result p, q;

  if (x == 0.0)
    {
      result->val = M_EULER;
      result->err = 0.0;
      return GSL_SUCCESS;
    }

  if (x == 1.0)
    {
      result->val = GSL_SF_1_MINUS_EULER;
      result->err = 0.0;
      return GSL_SUCCESS;
    }

  t = (4.0 * x - 1.0) / 3.0;

  status = cheb_eval_e(&num_cs, t, &p);
  if (status)
    return status;

  status = cheb_eval_e(&den_cs, t, &q);
  if (status)
    return status;

  result->val = 0.5 / (1.0 + x) + p.val / q.val;
  result->err = p.err / fabs(q.val) + fabs(p.val) * q.err / (q.val * q.val) +
                GSL_DBL_EPSILON * fabs(0.5 / (1.0 + x));

  return GSL_SUCCESS;
}

/*
 * Evaluate Gamma(1+x) on the central interval using
 *
 *   Gamma(1+x) = 1 + x (x-1) H(x).
 *
 * The error estimate from gamma1p_H() is propagated through the
 * reconstruction.
 */
static int
gamma1p_central(const double x, gsl_sf_result * result)
{
  int status;
  gsl_sf_result r;

  status = gamma1p_H(x, &r);
  if (status)
    return status;

  result->val = 1.0 + x * (x - 1.0) * r.val;
  result->err = fabs(x * (x - 1.0)) * r.err +
                2.0 * GSL_DBL_EPSILON * fabs(x * (x - 1.0) * r.val) +
                GSL_DBL_EPSILON * fabs(result->val);

  return GSL_SUCCESS;
}

/******************************************************
 * GAMMA KERNELS                                      *
 ******************************************************/

/*
 * Evaluate Gamma(x) for 2 <= x <= 10.
 *
 * Write
 *
 *   x = n + r,  n = floor(x),  0 <= r < 1,
 *
 * and use the recurrence
 *
 *   Gamma(x) = Gamma(1+r) prod_{k=1}^{n-1} (r+k).
 *
 * Gamma(1+r) is evaluated by gamma1p_central().  The returned error
 * estimate includes the central-kernel error and rounding errors in
 * forming and applying the recurrence product.
 */
static int
gamma_low(const double x, gsl_sf_result * result)
{
  int status;
  gsl_sf_result tmp;
  int n;
  int k;
  double r;
  double prod;

  n = (int) floor(x);
  r = x - (double) n;

  status = gamma1p_central(r, &tmp);
  if (status)
    return status;

  prod = 1.0;
  for (k = 0; k < n - 1; ++k)
    prod *= r + (double) (k + 1);

  result->val = tmp.val * prod;
  result->err = fabs(prod) * tmp.err +
                fabs(tmp.val) * (2.0 * (double) (n - 1) * GSL_DBL_EPSILON * fabs(prod)) +
                GSL_DBL_EPSILON * fabs(result->val);

  return GSL_SUCCESS;
}

/*
 * Form
 *
 *   prod = prod_{i=0}^{n-1} (x+i)
 *
 * using double-double arithmetic.
 *
 * The returned eps is a first-order relative correction containing
 * both the low component of the double-double product and the effect
 * of the argument correction x_eps:
 *
 *   eps = prod_lo/prod_hi + x_eps sum_{i=0}^{n-1} 1/(x+i).
 *
 * Thus prod*(1+eps) represents the corrected recurrence product to
 * first order.
 */
static void
gamma_product(const double x, const double x_eps, const int n,
              double *prod, double *eps)
{
  int i;
  gsl_dd p = { 1.0, 0.0 };
  double rel = 0.0;

  for (i = 0; i < n; ++i)
    {
      const double f = x + (double) i;

      p = gsl_dd_mul_d(p, f);
      rel += 1.0 / f;
    }

  *prod = p.hi;
  *eps = p.lo / p.hi + x_eps * rel;
}

/*
 * Evaluate the exponential correction in the Stirling expansion used
 * by gamma_high_scaled().
 *
 * For the adjusted argument x_adj and correction x_eps, return
 *
 *   x_eps log(x_adj) + S(x_adj),
 *
 * where S is the truncated odd inverse-power Stirling series
 *
 *   S(x) = c0/x + c1/x^3 + c2/x^5 + ...
 *
 * evaluated by Horner's method in 1/x^2.
 */
static double
gamma_stirling_exp_adj(const double x_adj, const double x_eps)
{
  static const double c[6] = { 8.333333333333332870740406e-02, -2.777777777777777883788657e-03,
                               7.936507936507936501052685e-04, -5.952380952380952917890600e-04,
                               8.417508417508417139715760e-04, -1.917526917526917633674555e-03 };
  double x2, b;
  int i;

  x2 = x_adj * x_adj;
  b = c[5];
  for (i = 1; i <= 5; ++i)
    b = b / x2 + c[5 - i];

  return x_eps * log(x_adj) + b / x_adj;
}

/*
 * Evaluate Gamma(x) in scaled form for x > 10,
 *
 *   Gamma(x) = mant * 2^exp2,
 *
 * with 0.5 <= |mant| < 1.
 *
 * The computation uses a Stirling reconstruction based on
 *
 *   Gamma(x) ~ sqrt(2*pi/x) (x/e)^x exp(S(x)).
 *
 * The power x^x is split into a normalized mantissa and an exact
 * power-of-two scale.  The remaining factors are combined using
 * double-double arithmetic to reduce intermediate rounding error.
 *
 * For 10 < x < 12, the argument is shifted upward and the recurrence
 * product required to recover Gamma(x) is formed by gamma_product().
 *
 * The exp2 and expm1 error estimates are propagated explicitly.
 * Additional allowance is made for the remaining floating-point
 * operations and, when applicable, the recurrence-product correction.
 *
 * On success, mant contains the normalized mantissa and its absolute
 * error estimate, while *exp2 contains the corresponding binary
 * exponent.
 */
static int
gamma_high_scaled(const double x, gsl_sf_result * mant, int * exp2)
{
  return gamma_high_scaled_eps(x, 0.0, mant, exp2);
}

/*
 * Evaluate Gamma(x+x_eps_in) in scaled form for x > 10,
 *
 *   Gamma(x+x_eps_in) = mant * 2^exp2,
 *
 * where x is the principal double-precision argument and x_eps_in is a
 * small correction to that argument.
 *
 * This is the corrected-argument form of gamma_high_scaled().  For
 * 10 < x < 12, the argument is shifted upward and the correction is
 * propagated through the recurrence product.  For larger x it enters
 * directly into the Stirling reconstruction.
 *
 * The computation uses double-double arithmetic for the principal
 * products and divisions.  Error estimates from exp2() and expm1(),
 * together with the recurrence-product correction and ordinary
 * floating-point rounding errors, are propagated to the returned
 * mantissa.
 *
 * On success, mant contains the normalized mantissa and its absolute
 * error estimate, while *exp2 contains the corresponding binary
 * exponent.
 */

static int
gamma_high_scaled_eps(const double x, const double x_eps_in,
                      gsl_sf_result * mant, int * exp2)
{
  int status;
  double x_adj = x;
  double prod = 1.0;
  double eps = 0.0;
  double x_eps, xi, xf;
  double mant0, ea, relerr;
  gsl_sf_result rexp2;
  gsl_sf_result rexpm1;
  gsl_dd a, b, c, d;
  int n = 0;
  int e2;

  if (x < 12.0)
    {
      n = (int) ceil(12.0 - x);
      x_adj = x + (double) n;
      x_eps = x - (x_adj - (double) n) + x_eps_in;
      gamma_product(x_adj - (double) n, x_eps, n, &prod, &eps);
    }
  else
    {
      x_eps = x_eps_in;
    }

  xi = floor(x_adj + 0.5);
  xf = x_adj - xi;

  mant0 = frexp(x_adj, &e2);
  if (mant0 < M_SQRT1_2)
    {
      --e2;
      mant0 *= 2.0;
    }

  /* a = mant0^x_adj * 2^(e2*xf). */
  a.hi = pow(mant0, x_adj);
  a.lo = 0.0;

  status = gsl_sf_exp2_e((double) e2 * xf, &rexp2);
  if (status)
    return status;

  a = gsl_dd_mul_d(a, rexp2.val);

  /* b = exp(-x_adj) * sqrt(2*pi/x_adj). */
  b.hi = exp(-x_adj);
  b.lo = 0.0;
  b = gsl_dd_mul_d(b, sqrt(2.0 * M_PI / x_adj));

  /* c = a*b. */
  c = gsl_dd_mul(a, b);

  /* Divide by the recurrence product. */
  d.hi = prod;
  d.lo = prod * eps;
  c = gsl_dd_div(c, d);

  /* Apply Stirling correction, including the x_eps correction. */
  ea = gamma_stirling_exp_adj(x_adj, x_eps);

  status = gsl_sf_expm1_e(ea, &rexpm1);
  if (status)
    return status;

  c.lo += c.hi * rexpm1.val;

  /* Return normalized mantissa. */
  mant->val = frexp(c.hi + c.lo, exp2);

  /* Relative error estimate. */
  relerr = 3.0 * GSL_DBL_EPSILON
         + rexp2.err / fabs(rexp2.val)
         + rexpm1.err / fabs(1.0 + rexpm1.val);

  if (n > 0)
    relerr += fabs(eps);

  mant->err = fabs(mant->val) * relerr;

  *exp2 += e2 * (int) xi;

  return GSL_SUCCESS;
}

/*
 * Evaluate Gamma(x) in scaled form for 2 <= x <= GAMMA_SCALED_X_MAX,
 *
 *   Gamma(x) = mant * 2^exp2,
 *
 * with 0.5 <= |mant| < 1.
 *
 * For x <= 10, Gamma(x) is evaluated by gamma_low() and then
 * normalized with frexp().  For larger arguments,
 * gamma_high_scaled() is used directly.
 *
 * The absolute error estimate is expressed in the same scale as the
 * returned mantissa.  Arguments beyond GAMMA_SCALED_X_MAX return
 * GSL_EOVRFLW.
 */
static int
gamma_2_inf_scaled(const double x, gsl_sf_result * mant, int * exp2)
{
  int status;
  gsl_sf_result r;

  if (x <= 10.0)
    {
      status = gamma_low(x, &r);
      if (status)
        return status;

      mant->val = frexp(r.val, exp2);
      mant->err = ldexp(r.err, -(*exp2));

      return GSL_SUCCESS;
    }

  if (x > GAMMA_SCALED_X_MAX)
    {
      return GSL_EOVRFLW;
    }

  return gamma_high_scaled(x, mant, exp2);
}

/*
 * Evaluate Gamma(x) for -1 < x < 0.
 *
 * Write r = x+1 and use
 *
 *   Gamma(x) = 1/[x(x+1)] + H(x+1),
 *
 * where H is defined by
 *
 *   Gamma(1+r) = 1 + r(r-1) H(r).
 *
 * This form avoids forming Gamma(x+1) followed by division by x,
 * which would amplify rounding error near x = 0.  The error estimate
 * includes the reciprocal-product error, the error in H, and the
 * final addition.
 */
static int
gamma_neg_small(const double x, gsl_sf_result * result)
{
  int status;
  const double r = x + 1.0;
  const double p = x * r;
  gsl_sf_result h;
  double a, ea, b, eb;

  status = gamma1p_H(r, &h);
  if (status)
    return status;

  a = 1.0 / p;
  ea = 3.0 * GSL_DBL_EPSILON * fabs(a);
  b = h.val;
  eb = h.err;

  result->val = a + b;
  result->err = ea + eb + GSL_DBL_EPSILON * fabs(result->val);

  if (gsl_isinf(result->val))
    {
      result->err = GSL_POSINF;
      return GSL_EOVRFLW;
    }

  return GSL_SUCCESS;
}

/*
 * Evaluate Gamma(x) for -10 <= x < -1 using recurrence.
 *
 * Shift the argument by
 *
 *   r = x+n,  0 < r < 1,
 *
 * and evaluate
 *
 *   Gamma(r) = 1/r + (r-1) H(r).
 *
 * The recurrence relation
 *
 *   Gamma(x) = Gamma(r) / [x (x+1) ... (x+n-1)]
 *
 * is then applied.  The denominator product and final division are
 * performed in double-double arithmetic to reduce accumulated
 * rounding error.
 *
 * The returned error estimate allows for rounding accumulated over
 * the recurrence.
 */
static int
gamma_neg_recurrence(const double x, gsl_sf_result * result)
{
  int status;
  const int n = (int) ceil(-x);
  const double r = x + (double) n;
  const double one_over_r = 1.0 / r;
  gsl_sf_result h;
  gsl_dd p = { 1.0, 0.0 };
  gsl_dd a, q;
  double term;
  int k;

  /* Gamma(r) = 1/r + (r-1) H(r), 0 < r < 1. */
  status = gamma1p_H(r, &h);
  if (status)
    return status;

  term = (r - 1.0) * h.val;

  a.hi = one_over_r + term; a.lo = 0.0;

  for (k = 0; k < n; ++k)
    p = gsl_dd_mul_d(p, x + (double) k);

  q = gsl_dd_div(a, p);
  result->val = q.hi + q.lo;
  result->err = 4.0 * (double) n * GSL_DBL_EPSILON * fabs(result->val);

  if (gsl_isinf(result->val))
    {
      result->err = GSL_POSINF;
      return GSL_EOVRFLW;
    }

  return GSL_SUCCESS;
}

/*
 * Evaluate Gamma(x) for x < -10 using scaled reflection.
 *
 * The reflection formula
 *
 *   Gamma(x) = pi / [sin(pi*x) Gamma(1-x)]
 *
 * is evaluated using gsl_sf_sin_pi_e() and the scaled positive-Gamma
 * kernel gamma_high_scaled_eps().  The correction y_eps preserves the
 * rounding residual introduced when forming
 *
 *   y = 1-x.
 *
 * Keeping Gamma(1-x) in scaled form avoids intermediate overflow in
 * the positive Gamma factor.  The errors in sin(pi*x), the scaled
 * Gamma value, the division, and the final scaling are propagated to
 * result->err.
 *
 * Negative integer arguments are poles and return GSL_EDOM.  Results
 * which round to zero return GSL_EUNDRFLW.  Nonzero subnormal results
 * use a conservative absolute error bound of 2*GSL_SF_TRUE_MIN,
 * since scaling the propagated error may itself underflow.
 */
static int
gamma_neg_reflection(const double x, gsl_sf_result * result)
{
  int status;
  const double y = 1.0 - x;
  const double y_eps = -x + (1.0 - y);
  gsl_sf_result s, gm;
  double num, epsnum, qmant, eqmant;
  int ge;

  status = gsl_sf_sin_pi_e(x, &s);
  if (status)
    return status;

  if (s.val == 0.0)
    DOMAIN_ERROR(result);

  if (y > GAMMA_SCALED_X_MAX)
    {
      result->val = 0.0;
      result->err = GSL_SF_TRUE_MIN;
      return GSL_EUNDRFLW;
    }

  status = gamma_high_scaled_eps(y, y_eps, &gm, &ge);
  if (status)
    return status;

  /* Gamma(x) = pi / (sin(pi*x) Gamma(1-x)). */
  num = M_PI / s.val;
  epsnum = M_PI * s.err / (s.val * s.val)
         + GSL_DBL_EPSILON * fabs(num);

  qmant = num / gm.val;
  eqmant = epsnum / fabs(gm.val)
         + fabs(num) * gm.err / (gm.val * gm.val)
         + 2.0 * GSL_DBL_EPSILON * fabs(qmant);

  result->val = ldexp(qmant, -ge);

  if (gsl_isinf(result->val))
    {
      result->err = GSL_POSINF;
      return GSL_EOVRFLW;
    }
  else if (result->val == 0.0)
    {
      result->err = GSL_SF_TRUE_MIN;
      return GSL_EUNDRFLW;
    }
  else if (fabs(result->val) < GSL_DBL_MIN)
    {
      /*
       * Scaling the propagated error can itself underflow.  Use a
       * conservative absolute error bound in the subnormal range.
       */
      result->err = 2.0 * GSL_SF_TRUE_MIN;
      return GSL_SUCCESS;
    }
  else
    {
      result->err = ldexp(eqmant, -ge);
      return GSL_SUCCESS;
    }
}
