/* specfunc/gamma/gamma_inc.c
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
#include <gsl/gsl_dd.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_gamma.h>

#include "constants.h"
#include "error.h"
#include "eval.h"

#define GAMMA_INC_MAXITER         10000

typedef enum { TARGET_P = 0, TARGET_Q = 1 } gamma_inc_target_t;

/* Taylor kernel functions */
static int gamma_inc_P_taylor_x_sum(const double a, const double x, double * temp, double * err);
static int gamma_inc_P_taylor_x(const double a, const double x, gsl_sf_result * result);

/* P series kernel functions */
static int gamma_inc_p_series_sum(const double a, const double x, double * sum, int * nterms);
static int gamma_inc_P_series(const double a, const double x, gsl_sf_result * result);

/* Q kernel functions */
static int gamma_inc_Q_small_a(const double a, const double x, gsl_sf_result * result);
static int gamma_inc_Q_cf(const double a, const double x, gsl_sf_result * result);
static int gamma_inc_Q_asymp(const double a, const double x, gsl_sf_result * result);

/* Temme kernel functions */
static double gamma_inc_half_erfc_dd(const gsl_dd z);
static int gamma_inc_use_temme(const double a, const double x);
static int gamma_inc_temme_eval(const double a, const double x, const int target, gsl_sf_result * result);

/* rgammax primitive functions */
static gsl_dd gamma_inc_log_ratio_m_ratio_p1(const double a, const double x);
static gsl_dd gamma_inc_log_ratio_near1(const double a, const double x);
static double gamma_inc_rgammax_large(const double a, const double x);
static double gamma_inc_rgammax(const double a, const double x);

/* double-double routines */
static gsl_dd dd_div_diff(const double a, const double x);

/*
 * Compute the regularized lower incomplete gamma function
 *
 *   P(a,x) = gamma(a,x) / Gamma(a),
 *
 * for a > 0 and x >= 0.
 *
 * For small a and x, use the dedicated Taylor kernel.  Near x = a
 * for sufficiently large a, use the Temme uniform asymptotic expansion.
 * On the lower side of the transition, evaluate P directly with the
 * positive series.  On the upper side, evaluate the complementary
 * quantity Q with the validated Q dispatcher and form
 *
 *   P(a,x) = 1 - Q(a,x).
 *
 * Using the complementary upper-tail representation for x > a avoids
 * the severe loss of accuracy that occurs when the positive P series
 * is used far into the upper tail.
 *
 * The limiting value at x = 0 is P(a,0) = 0.  Arguments with a <= 0
 * or x < 0 are outside the real-valued domain of this interface and
 * produce a domain error.
 */

int
gsl_new_sf_gamma_inc_P_e(const double a, const double x,
                         gsl_sf_result * result)
{
  int status;

  if (a <= 0.0 || x < 0.0)
    {
      DOMAIN_ERROR(result);
    }
  else if (x == 0.0)
    {
      result->val = 0.0;
      result->err = 0.0;
      return GSL_SUCCESS;
    }
  else if (a < 1.0 && x < 1.1)
    {
      return gamma_inc_P_taylor_x(a, x, result);
    }
  else if (gamma_inc_use_temme(a, x))
    {
      return gamma_inc_temme_eval(a, x, TARGET_P, result);
    }
  else if (x <= a)
    {
      return gamma_inc_P_series(a, x, result);
    }
  else
    {
      gsl_sf_result q;

      status = gsl_new_sf_gamma_inc_Q_e(a, x, &q);
      if (status)
        return status;

      result->val = 1.0 - q.val;

      /*
       * Protect against a tiny over/undershoot introduced by the
       * subtraction or by the Q kernel's final rounding.
       */
      result->val = GSL_MAX(result->val, 0.0);
      result->val = GSL_MIN(result->val, 1.0);

      result->err = q.err + 2.0 * GSL_DBL_EPSILON * fabs(result->val) + GSL_SF_TRUE_MIN;

      return GSL_SUCCESS;
    }
}

/*
 * Compute the regularized upper incomplete gamma function
 *
 *   Q(a,x) = Gamma(a,x) / Gamma(a),
 *
 * for a > 0 and x >= 0.
 *
 * The parameter plane is divided among several numerical methods.  For
 * 0 < a <= 2, use the direct small-a Taylor/GRAT1 kernel.  Near the
 * transition x = a for sufficiently large a, use the Temme uniform
 * asymptotic expansion.  On the lower side of the transition, form
 * Q = 1-P from the positive P series.  On the upper side, use either
 * P complementation, the normalized continued fraction, or the large-x
 * asymptotic expansion according to the validated dispatcher boundaries.
 *
 * The limiting value at x = 0 is Q(a,0) = 1.  Arguments with a <= 0
 * or x < 0 are outside the real-valued domain of this interface and
 * produce a domain error.
 */

int
gsl_new_sf_gamma_inc_Q_e(const double a, const double x,
                         gsl_sf_result * result)
{
  int status;

  if (a <= 0.0 || x < 0.0)
    {
      DOMAIN_ERROR(result);
    }
  else if (x == 0.0)
    {
      result->val = 1.0;
      result->err = 0.0;
      return GSL_SUCCESS;
    }
  else if (a <= 2.0)
    {
      return gamma_inc_Q_small_a(a, x, result);
    }
  else if (gamma_inc_use_temme(a, x))
    {
      return gamma_inc_temme_eval(a, x, TARGET_Q, result);
    }
  else if (x <= a)
    {
      gsl_sf_result p;

      status = gamma_inc_P_series(a, x, &p);
      if (status)
        return status;

      result->val = 1.0 - p.val;
      result->val = GSL_MAX(result->val, 0.0);
      result->val = GSL_MIN(result->val, 1.0);

      result->err = p.err + 2.0 * GSL_DBL_EPSILON * fabs(result->val) + GSL_SF_TRUE_MIN;

      return GSL_SUCCESS;
    }
  else if (a < 12.0)
    {
      if (x >= 40.0)
        return gamma_inc_Q_asymp(a, x, result);
      else
        {
          gsl_sf_result p;

          status = gamma_inc_P_series(a, x, &p);
          if (status)
            return status;

          if (1.0 - p.val >= 0.34)
            {
              result->val = 1.0 - p.val;
              result->val = GSL_MAX(result->val, 0.0);
              result->val = GSL_MIN(result->val, 1.0);

              result->err = p.err + 2.0 * GSL_DBL_EPSILON * fabs(result->val) + GSL_SF_TRUE_MIN;

              return GSL_SUCCESS;
            }

          return gamma_inc_Q_cf(a, x, result);
        }
    }
  else
    {
      const double s = (x - a) / a;

      if (x >= 200.0 || (a >= 20.0 && s >= 1.5))
        return gamma_inc_Q_asymp(a, x, result);

      return gamma_inc_Q_cf(a, x, result);
    }
}

double gsl_new_sf_gamma_inc_P(const double a, const double x) { EVAL_RESULT(gsl_new_sf_gamma_inc_P_e(a, x, &result)); }
double gsl_new_sf_gamma_inc_Q(const double a, const double x) { EVAL_RESULT(gsl_new_sf_gamma_inc_Q_e(a, x, &result)); }

/*****************************************
 * INTERNAL ROUTINES (Taylor kernel)     *
 *****************************************/

/*
 * Compute the Taylor correction used in gamma_inc_P_taylor_x().
 *
 * For a > 0, write
 *
 *   P(a,x) = x^a / Gamma(a+1) * (1 - J(a,x)),
 *
 * where
 *
 *   J(a,x) = -a sum_{n=1}^infinity (-x)^n / ((a+n) n!).
 *
 * The first two terms are evaluated explicitly and the remaining
 * convergent alternating series is accumulated recursively.  The
 * returned value temp is J(a,x), and err is an estimate of its
 * absolute error.
 *
 * This expansion is based on the small-x algorithm of
 * DiDonato and Morris (1986).
 *
 * Intended domain: 0 < a < 1, 0 < x < 1.1.
 */

static int
gamma_inc_P_taylor_x_sum(const double a, const double x,
                         double * temp, double * err)
{
  const double tol = 3.0 * GSL_DBL_EPSILON / (a + 1.0);
  double c;
  double sm;
  double t;
  double l;

  /* Accumulate the remainder beginning with the n=3 term. */
  l = 3.0;
  c = x;
  sm = x / (a + 3.0);

  while (1)
    {
      l += 1.0;
      c *= -x / l;
      t = c / (a + l);
      sm += t;
      if (fabs(t) <= tol)
        break;

      if (l >= GAMMA_INC_MAXITER)
        return GSL_EMAXITER;
    }

  *temp = a * x * ((sm / 6.0 - 0.5 / (a + 2.0)) * x + 1.0 / (a + 1.0));

  /*
   * First-pass series estimate.  Julia's original implementation does
   * not return a GSL-style error estimate; this estimate is for comparison
   * diagnostics only.
   */
  *err = fabs(t) + 8.0 * GSL_DBL_EPSILON * fabs(*temp);

  return GSL_SUCCESS;
}

/*
 * Compute P(a,x) for small a and x using a Taylor expansion in x.
 *
 * The representation is
 *
 *   P(a,x) = x^a / Gamma(a+1) * (1 - J(a,x)),
 *
 * where J(a,x) is evaluated by gamma_inc_P_taylor_x_sum().
 *
 * The reciprocal Gamma factor is formed from gsl_sf_lngamma1p_e(),
 *
 *   1 / Gamma(a+1) = exp(-lngamma1p(a)),
 *
 * and x^a is evaluated with pow() in the direct-P branch.  When P
 * is close to unity, an algebraically equivalent formulation using
 * expm1() is used to avoid loss of accuracy from cancellation.
 *
 * Intended domain: 0 < a < 1, 0 < x < 1.1.
 */

static int
gamma_inc_P_taylor_x(const double a, const double x, gsl_sf_result * result)
{
  int status;
  gsl_sf_result lg1p;
  double temp, temp_err;
  double z;
  double g;
  double p;

  status = gamma_inc_P_taylor_x_sum(a, x, &temp, &temp_err);
  if (status)
    return status;

  z = a * log(x);
  status = gsl_sf_lngamma1p_e(a, &lg1p);
  if (status)
    return status;

  g = exp(-lg1p.val);

  /*
   * Use an equivalent Q-based reconstruction when P is close to one
   * to avoid cancellation in the direct evaluation.
   */
  if ((x < 0.25 && z > -0.13394) || a < x / 2.59)
    {
      const double em1 = gsl_sf_expm1(z);
      const double w = 1.0 + em1;
      const double q = GSL_MAX((w * temp - em1) * g - (g - 1.0), 0.0);
      p = 1.0 - q;
    }
  else
    {
      const double w = pow(x, a);
      p = w * g * (1.0 - temp);
    }

  result->val = p;
  result->err = fabs(p) * (fabs(z) * 2.0 * GSL_DBL_EPSILON +
                fabs(temp_err) + 16.0 * GSL_DBL_EPSILON) +
                fabs(p) * lg1p.err + GSL_SF_TRUE_MIN;

  return GSL_SUCCESS;
}

/*****************************************
 * INTERNAL ROUTINES (P series kernel)   *
 *****************************************/

/*
 * Compute the positive series
 *
 *   S(a,x) = 1 + x/(a+1) + x^2/[(a+1)(a+2)] + ...
 *
 * appearing in the representation
 *
 *   P(a,x) = rgammax(a,x) S(a,x) / a,
 *
 * where
 *
 *   rgammax(a,x) = x^a exp(-x) / Gamma(a).
 *
 * The terms are generated recursively and accumulated using compensated
 * summation.  The number of accumulated terms is returned in nterms.
 *
 * Returns GSL_EMAXITER if the convergence criterion is not satisfied
 * within GAMMA_INC_MAXITER iterations.
 */

static int
gamma_inc_p_series_sum(const double a, const double x, double * sum, int * nterms)
{
  double term = 1.0;
  double s = 1.0;
  double c = 0.0;
  int n;

  for (n = 1; n < GAMMA_INC_MAXITER; ++n)
    {
      double y, t;

      term *= x / (a + (double)n);
      y = term - c;
      t = s + y;
      c = (t - s) - y;
      s = t;

      if (fabs(term) <= GSL_DBL_EPSILON * fabs(s))
        {
          *sum = s;
          *nterms = n + 1;
          return GSL_SUCCESS;
        }
    }

  *sum = s;
  *nterms = n;

  return GSL_EMAXITER;
}

/*
 * Compute the regularized lower incomplete gamma function from the
 * positive P series,
 *
 *   P(a,x) = rgammax(a,x) S(a,x) / a,
 *
 * with S(a,x) evaluated by gamma_inc_p_series_sum() and the scaled
 * gamma factor evaluated by gamma_inc_rgammax().
 *
 * This method is intended for the lower-tail region outside the
 * small-a Taylor and Temme regions selected by the dispatcher.
 */

static int
gamma_inc_P_series(const double a, const double x, gsl_sf_result * result)
{
  int status;
  double s;
  int n;

  status = gamma_inc_p_series_sum(a, x, &s, &n);
  if (status)
    return status;

  result->val = (gamma_inc_rgammax(a,x) / a) * s;
  result->err = (n + 16.0) * GSL_DBL_EPSILON * fabs(result->val) + GSL_SF_TRUE_MIN;

  return GSL_SUCCESS;
}

/*****************************************
 * INTERNAL ROUTINES (Q kernels)         *
 *****************************************/

/*
 * Compute Q(a,x) directly for
 *
 *   0 < a <= 2.
 *
 * For x < 1.1, use the small-x representation derived from
 *
 *   P(a,x) = x^a / Gamma(a+1) * (1 - J(a,x)),
 *
 * and reconstruct Q = 1-P using double-double arithmetic in the
 * cancellation-sensitive terms.  The logarithm of x, the exponential
 * correction, and the Taylor correction J(a,x) are retained in
 * double-double form where needed.
 *
 * For x >= 1.1, use the GRAT1 continued-fraction representation for Q.
 * The continued fraction and its exponential prefactor are evaluated
 * in double-double arithmetic.
 *
 * Intended domain: 0 < a <= 2, x > 0.
 */

static int
gamma_inc_Q_small_a(const double a, const double x, gsl_sf_result * result)
{
  int status;
  const double tol = 0.1 * GSL_DBL_EPSILON / (a + 1.0);
  double L, h;
  gsl_sf_result lg1p;
  int n;

  status = gsl_sf_lngamma1p_e(a, &lg1p);
  if (status)
    return status;

  L = lg1p.val;
  h = gsl_sf_expm1(-L);

  if (x < 1.1)
    {
      gsl_dd an, c, sum, t, j, z, logx;
      gsl_dd q1, q2, q3, oneg, hd, jd, ld, pd, omj, gdd;

      an.hi = 3.0;
      an.lo = 0.0;
      c.hi = x;
      c.lo = 0.0;
      sum = gsl_dd_div_d(c, a + 3.0);

      do
        {
          an.hi += 1.0;
          c = gsl_dd_mul_d(c, -x);
          c = gsl_dd_div_d(c, an.hi);
          t = gsl_dd_div_d(c, a + an.hi);
          sum = gsl_dd_add(sum, t);
          if (an.hi >= GAMMA_INC_MAXITER)
            return GSL_EMAXITER;
        }
      while (fabs(t.hi + t.lo) > tol);

      q1 = gsl_dd_mul_d(sum, 1.0/6.0);
      q2 = gsl_dd_add(q1, gsl_dd_div_d((gsl_dd){-0.5, 0.0}, a + 2.0));
      q3 = gsl_dd_mul_d(q2, x);
      q3 = gsl_dd_add(q3, gsl_dd_div_d((gsl_dd){1.0, 0.0}, a + 1.0));
      j = gsl_dd_mul_d(q3, a);
      j = gsl_dd_mul_d(j, x);

      if (x == 1.0)
        {
          logx.hi = 0.0;
          logx.lo = 0.0;
        }
      else
        logx = gsl_dd_log((gsl_dd){x, 0.0});

      z = gsl_dd_mul_d(logx, a);

      hd.hi = h;
      hd.lo = 0.0;

      gdd = gsl_dd_add((gsl_dd){1.0, 0.0}, hd);
      oneg = gdd;
      jd = j;

      if (x < 0.25 && z.hi + z.lo <= -0.13394)
        {
          gsl_dd ez;

          ez = gsl_dd_exp(z);
          omj.hi = 1.0 - jd.hi;
          omj.lo = -jd.lo;
          pd = gsl_dd_mul(ez, omj);
          pd = gsl_dd_mul(pd, oneg);

          {
            gsl_dd one_minus_p;
            one_minus_p = gsl_dd_add((gsl_dd){1.0, 0.0}, (gsl_dd){-pd.hi, -pd.lo});
            result->val = one_minus_p.hi + one_minus_p.lo;
          }

          result->err = fabs(pd.hi + pd.lo) * (fabs(z.hi + z.lo) + 16.0) * GSL_DBL_EPSILON +
                        fabs(result->val) * 4.0 * GSL_DBL_EPSILON + GSL_SF_TRUE_MIN;
        }
      else
        {
          ld = gsl_dd_expm1(z);

          {
            gsl_dd one_plus_l, term1, term2;
            one_plus_l = gsl_dd_add((gsl_dd){1.0, 0.0}, ld);
            term1 = gsl_dd_mul(one_plus_l, jd);
            term2 = gsl_dd_mul_d(ld, -1.0);
            term1 = gsl_dd_add(term1, term2);
            term1 = gsl_dd_mul(term1, oneg);
            hd.hi = -hd.hi; hd.lo = -hd.lo;
            term1 = gsl_dd_add(term1, hd);

            result->val = term1.hi + term1.lo;
            result->err = 64.0 * GSL_DBL_EPSILON * fabs(result->val) +
                          fabs(ld.hi) * GSL_DBL_EPSILON + fabs(jd.hi) * GSL_DBL_EPSILON + fabs(L) * GSL_DBL_EPSILON +
                          GSL_SF_TRUE_MIN;
          }
        }

      return 0;
    }

  /* GRAT1 continued fraction, entirely in double-double arithmetic. */
  {
    gsl_dd A2nm1, A2n, B2nm1, B2n, am0, an0, tmp, logx, logcore, rrdd;
    double c = 1.0;

    A2nm1.hi = 1.0; A2nm1.lo = 0.0;
    A2n.hi = 1.0; A2n.lo = 0.0;
    B2nm1.hi = x; B2nm1.lo = 0.0;
    B2n.hi = x + 1.0 - a; B2n.lo = 0.0;

    logx = gsl_dd_log((gsl_dd){x, 0.0});
    logcore = gsl_dd_mul_d(logx,a);
    logcore = gsl_dd_add(logcore, (gsl_dd){-x, 0.0});
    logcore = gsl_dd_add(logcore, (gsl_dd){-L, 0.0});
    rrdd.hi = a; rrdd.lo = 0.0;
    tmp = gsl_dd_exp(logcore);
    rrdd = gsl_dd_mul(rrdd, tmp);

    for (n = 0; n < GAMMA_INC_MAXITER; n++)
      {
        A2nm1 = gsl_dd_add(gsl_dd_mul_d(A2n, x), gsl_dd_mul_d(A2nm1, c));
        B2nm1 = gsl_dd_add(gsl_dd_mul_d(B2n, x), gsl_dd_mul_d(B2nm1, c));
        am0 = gsl_dd_div(A2nm1, B2nm1);
        c += 1.0;
        tmp = gsl_dd_mul_d(A2n, c - a);
        A2n = gsl_dd_add(A2nm1, tmp);
        tmp = gsl_dd_mul_d(B2n, c - a);
        B2n = gsl_dd_add(B2nm1, tmp);
        an0 = gsl_dd_div(A2n, B2n);

        if (fabs((an0.hi - am0.hi) + (an0.lo - am0.lo)) < tol * fabs(an0.hi + an0.lo))
          break;
      }

    if (n == GAMMA_INC_MAXITER)
      return GSL_EMAXITER;

    rrdd = gsl_dd_mul(rrdd, an0);
    result->val = rrdd.hi + rrdd.lo;
    result->err = 32.0 * GSL_DBL_EPSILON * fabs(result->val) + GSL_DBL_EPSILON * fabs(rrdd.hi) +
                  GSL_SF_TRUE_MIN;
  }

  return GSL_SUCCESS;
}

/*
 * Compute Q(a,x) on the upper side of the transition using the
 * normalized continued fraction
 *
 *   Q(a,x) = rgammax(a,x) C(a,x),
 *
 * where
 *
 *   rgammax(a,x) = x^a exp(-x) / Gamma(a)
 *
 * and C(a,x) is evaluated by the normalized Julia continued-fraction
 * recurrence.  The recurrence state and normalization divisions are
 * retained in double-double arithmetic to reduce rounding error.
 *
 * Iteration stops when successive normalized convergents agree to
 * approximately two machine epsilons.  Returns GSL_EMAXITER if the
 * convergence criterion is not satisfied within GAMMA_INC_MAXITER
 * iterations.
 *
 * The caller is responsible for selecting the validated continued-
 * fraction region.
 */

static int
gamma_inc_Q_cf(const double a, const double x, gsl_sf_result * result)
{
  gsl_dd a2nm1 = {1.0, 0.0}, a2n = {1.0, 0.0};
  gsl_dd b2nm1 = {x, 0.0}, b2n = {x + 1.0 - a, 0.0};
  gsl_dd prev = {0.0, 0.0};
  double c = 1.0;
  int n;

  for (n = 0; n < GAMMA_INC_MAXITER; n++)
    {
      gsl_dd curr, diff;

      a2nm1 = gsl_dd_add(gsl_dd_mul_d(a2n, x), gsl_dd_mul_d(a2nm1, c));
      b2nm1 = gsl_dd_add(gsl_dd_mul_d(b2n, x), gsl_dd_mul_d(b2nm1, c));

      c += 1.0;

      a2n = gsl_dd_add(a2nm1, gsl_dd_mul_d(a2n, c - a));
      b2n = gsl_dd_add(b2nm1, gsl_dd_mul_d(b2n, c - a));

      a2nm1 = gsl_dd_div(a2nm1, b2n);
      b2nm1 = gsl_dd_div(b2nm1, b2n);
      a2n = gsl_dd_div(a2n, b2n);

      b2n.hi = 1.0; b2n.lo = 0.0;
      curr = a2n;

      ++n;
      if (n > 1)
        {
          diff = gsl_dd_sub(curr, prev);
          if (fabs(diff.hi + diff.lo) <= 2.0 * GSL_DBL_EPSILON * fabs(curr.hi + curr.lo))
            break;
        }

      prev = curr;
    }

  if (n >= GAMMA_INC_MAXITER)
    return GSL_EMAXITER;

  result->val = gamma_inc_rgammax(a, x) * (a2n.hi + a2n.lo);
  result->err = 32.0 * GSL_DBL_EPSILON * fabs(result->val) + GSL_SF_TRUE_MIN;

  return GSL_SUCCESS;
}

/*
 * Compute Q(a,x) in the large-x upper tail using the asymptotic
 * expansion
 *
 *   Q(a,x) ~ rgammax(a,x) / x * [1 + (a-1)/x + (a-1)(a-2)/x^2 + ...],
 *
 * where
 *
 *   rgammax(a,x) = x^a exp(-x) / Gamma(a).
 *
 * Terms are accumulated until they reach machine precision or begin
 * increasing in magnitude, giving optimal truncation of the asymptotic
 * series.
 *
 * For moderate a and very large x, the ordinary rgammax evaluation can
 * underflow internally before Q itself becomes unrepresentable.  In
 * this region the logarithmic prefactor a*log(x)-x is reconstructed in
 * double-double arithmetic before exponentiation.
 *
 * The caller is responsible for selecting the validated asymptotic
 * region.
 */

static int
gamma_inc_Q_asymp(const double a, const double x, gsl_sf_result * result)
{
  double term = 1.0, sum = 1.0, prev = 1.0;
  int n;
  double r;

  for (n = 1; n < GAMMA_INC_MAXITER; n++)
    {
      term *= (a - (double)n) / x;
      if (fabs(term) > prev)
        break;

      sum += term;
      prev = fabs(term);
      if (term == 0.0 || fabs(term) <= GSL_DBL_EPSILON * fabs(sum))
        break;
    }

  if (n >= GAMMA_INC_MAXITER)
    return GSL_EMAXITER;

  r = gamma_inc_rgammax(a, x);
  if(a <= 20.0 && x > 700.0)
    {
      /* The moderate-a rgammax branch forms exp(-x) directly.  In the
       * far tail that can enter the subnormal range before Q does.
       * Reconstruct log(R) in binary64/DD arithmetic instead. */
      gsl_dd lx, lr, ta, negx, er;
      lx = gsl_dd_log((gsl_dd){x, 0.0});
      ta = gsl_dd_mul_d(lx, a);
      negx.hi = -x; negx.lo = 0.0;
      lr = gsl_dd_add(ta, negx);
      if(lr.hi > GSL_SF_LOG_TRUE_MIN)
        {
          er = gsl_dd_exp(lr);
          r = (er.hi + er.lo) / tgamma(a);
        }
      else
        r = 0.0;
    }

  result->val = (r / x) * sum;
  result->err = (n + 16.0) * GSL_DBL_EPSILON * fabs(result->val) + GSL_SF_TRUE_MIN;

  return GSL_SUCCESS;
}

/*****************************************
 * INTERNAL ROUTINES (Temme kernel)      *
 *****************************************/

/*
 * Evaluate one half of erfc(z) for a double-double argument
 *
 *   z = z.hi + z.lo.
 *
 * The leading value is evaluated at z.hi and corrected to first order
 * for z.lo,
 *
 *   erfc(z.hi + z.lo)/2 = erfc(z.hi)/2 - exp(-z.hi^2) z.lo / sqrt(pi) + O(z.lo^2).
 *
 * Retaining this correction is important in the small tail, where
 * rounding the erfc argument to a single double can noticeably amplify
 * the final ULP error.
 */
static double
gamma_inc_half_erfc_dd(const gsl_dd z)
{
  return 0.5 * gsl_sf_erfc(z.hi) - exp(-z.hi * z.hi) / M_SQRTPI * z.lo;
}

/*
 * Return nonzero if the Temme uniform asymptotic expansion should be
 * used for the point (a,x).
 *
 * Define
 *
 *   s = (x-a)/a,
 *   eta^2 / 2 = s - log(1+s),
 *
 * with sign(eta) = sign(x-a).  The present dispatcher uses Temme for
 *
 *   a >= 12,  |eta| <= 0.5.
 *
 * Only |eta| is needed here, so its sign is not formed explicitly.
 *
 * This routine is a Boolean predicate; its return value is not a GSL
 * error status.
 */

static int
gamma_inc_use_temme(const double a, const double x)
{
  gsl_dd s;
  gsl_dd q;
  double eta;

  if (a < 12.0 || x <= 0.0)
    return 0;

  s = dd_div_diff(a, x);
  q = gsl_dd_log1pmx(s);

  eta = sqrt(GSL_MAX(-2.0 * (q.hi + q.lo), 0.0));

  return eta <= 0.5;
}

/*
 * Compute P(a,x) or Q(a,x) using the Temme uniform asymptotic expansion
 * for the incomplete gamma ratios.
 *
 * With
 *
 *   s = (x-a)/a,
 *   phi = s - log(1+s) = eta^2/2,
 *
 * the leading transition term is an erfc function of
 *
 *   eta sqrt(a/2),
 *
 * with an asymptotic correction evaluated using the order-24 backward
 * recurrence and the coefficients gamma_inc_temme_d[].
 *
 * TARGET_P selects the lower ratio P(a,x), while TARGET_Q selects the
 * upper ratio Q(a,x).  The reduced variables and exponential scaling
 * are retained in double-double arithmetic where necessary to avoid
 * loss of accuracy for large a.
 *
 * The caller is responsible for ensuring that (a,x) lies in the
 * validated Temme region.
 */

static int
gamma_inc_temme_eval(const double a, const double x, const int target, gsl_sf_result * result)
{
  /*
   * Universal coefficients for the order-24 Temme asymptotic correction.
   * These coefficients enter the backward recurrence used by
   * gamma_inc_temme_eval().
   */
  static const double gamma_inc_temme_d[26] = { 1.0,                    -0.33333333333333331,     0.083333333333333329,   -0.014814814814814815,
                                                0.0011574074074074073,   0.00035273368606701942, -0.0001787551440329218,   3.9192631785224377e-05,
                                               -2.185448510679992e-06,  -1.85406221071516e-06,    8.2967113409530865e-07, -1.7665952736826078e-07,
                                                6.7078535434014984e-09,  1.0261809784240309e-08, -4.3820360184533529e-09,  9.1476995822367902e-10,
                                               -2.5514193994946248e-11, -5.8307721325504256e-11,  2.4361948020667415e-11, -5.0276692801141755e-12,
                                                1.1004392031956135e-13,  3.3717632624009851e-13, -1.3923887224181621e-13,  2.8534893807047445e-14,
                                               -5.1391118342425723e-16, -1.9752288294349442e-15 };

  const int N = 24;
  gsl_dd s, q, phi, aphi, zdd;
  double beta[27], sum, scale, R, value, last;
  int n;

  s = dd_div_diff(a, x);
  q = gsl_dd_log1pmx(s);
  phi.hi = -q.hi;
  phi.lo = -q.lo;

  if (phi.hi + phi.lo < 0)
    {
      phi.hi = 0;
      phi.lo = 0;
    }

  zdd = gsl_dd_sqrt(gsl_dd_mul_d(phi, 2.0));

  if (x < a)
    {
      zdd.hi = -zdd.hi;
      zdd.lo = -zdd.lo;
    }

  for (n = N + 1; n <= N + 2; ++n)
    beta[n] = 0;

  for (n = N; n >= 0; --n)
    beta[n] = ((n + 2.0) * beta[n + 2]) / a + gamma_inc_temme_d[n + 1];

  sum = beta[N];
  for (n = N - 1; n >= 0; --n)
    sum = sum * zdd.hi + beta[n];

  sum *= a / (a + beta[1]);
  aphi = gsl_dd_mul_d(phi, a);
  scale = exp(-aphi.hi) * exp(-aphi.lo) / sqrt(2*M_PI*a);
  R = scale * sum;

  if (aphi.hi == 0.0 && aphi.lo == 0.0)
    {
      zdd.hi = 0.0;
      zdd.lo = 0.0;
    }
  else
    {
      zdd = gsl_dd_sqrt(aphi);
    }

  if (x < a)
    {
      zdd.hi = -zdd.hi;
      zdd.lo = -zdd.lo;
    }

  if (target == TARGET_Q)
    value = gamma_inc_half_erfc_dd(zdd) + R;
  else
    {
      zdd.hi = -zdd.hi;
      zdd.lo = -zdd.lo;
      value = gamma_inc_half_erfc_dd(zdd) - R;
    }

  last = fabs(beta[N] * pow(zdd.hi < 0 ? -zdd.hi : zdd.hi, (double)N));

  result->val = value;
  result->err = 4.0 * GSL_DBL_EPSILON * fabs(value) + fabs(scale) * last +
                64.0 * GSL_DBL_EPSILON * fabs(R) + GSL_SF_TRUE_MIN;

  return GSL_SUCCESS;
}

/*****************************************
 * INTERNAL ROUTINES (rgammax primitive) *
 *****************************************/

/* q = log(x/a)-x/a+1 using a compensated quotient directly. */
static gsl_dd
gamma_inc_log_ratio_m_ratio_p1(const double a, const double x)
{
  gsl_dd t, lg, nt, one;

  t.hi = x;
  t.lo = 0.0;
  t = gsl_dd_div_d(t, a);
  lg = gsl_dd_log(t);
  nt.hi = -t.hi;
  nt.lo = -t.lo;
  one.hi = 1.0;
  one.lo = 0.0;

  return gsl_dd_add(gsl_dd_add(lg, nt), one);
}

/* q = log(1+s)-s, s=(x-a)/a, using the compensated near-one form. */
static gsl_dd
gamma_inc_log_ratio_near1(const double a, const double x)
{
  const gsl_dd s = dd_div_diff(a, x);
  return gsl_dd_log1pmx(s);
}

/*
 * Compute the Stirling correction to log Gamma(a),
 *
 *   E(a) = log Gamma(a) - [(a - 1/2) log(a) - a + log(sqrt(2*pi))].
 *
 * For a < 1000, use the rational approximation
 *
 *   E(a) = p(z) / [(c + z) a],
 *
 * where z = 1/a^2 and p(z) is evaluated by Horner's method.
 *
 * For a >= 1000, use the asymptotic expansion
 *
 *   E(a) = [1/12 - z/360 + z^2/1260 - z^3/1680] / a.
 */
static double
gamma_inc_stirling_error(const double a)
{
  const double denom = 0.30865217988013567769;
  static const double c[6] = {
    0.25721014990011306473e-1,
    0.82475966166999631057e-1,
   -0.25328157302663562668e-2,
    0.60992926669463371e-3,
   -0.33543297638406e-3,
    0.250505279903e-3
  };

  const double z = 1.0 / (a * a);

  if (a < 1000.0)
    {
      double p = c[0] + z * (c[1] + z * (c[2] + z * (c[3] + z * (c[4] + z * c[5]))));
      return p / ((denom + z) * a);
    }
  else
    {
      double p;

      p = -1.0 / 1680.0;
      p =  1.0 / 1260.0 + z * p;
      p = -1.0 / 360.0  + z * p;
      p =  1.0 / 12.0   + z * p;

      return p / a;
    }
}

/*
 * Compute
 *
 *   rgammax(a,x) = x^a exp(-x) / Gamma(a)
 *
 * for large a using the scaled representation
 *
 *   rgammax(a,x) = sqrt(a/(2*pi)) exp[a q - E(a)],
 *
 * where
 *
 *   q = log(x/a) - x/a + 1
 *
 * and E(a) is the Stirling correction returned by
 * gamma_inc_stirling_error().
 *
 * The reduced logarithmic quantity q and the product a*q are retained
 * in double-double arithmetic.  Near x=a, q is formed from
 * s=(x-a)/a to preserve the accurately computed difference x-a.
 *
 * When exp(a*q-E) approaches underflow, the sqrt(a/(2*pi)) prefactor
 * is absorbed into the exponent before evaluation.
 */
static double
gamma_inc_rgammax_large(const double a, const double x)
{
  const double E = gamma_inc_stirling_error(a);
  gsl_dd q, aq, w;

  /* Near x=a use the exact-difference representation.  Away from it,
   * form the quotient x/a directly as a compensated pair. */
  if (x >= 0.5 * a && x <= 2.0 * a)
    q = gamma_inc_log_ratio_near1(a, x);
  else
    q = gamma_inc_log_ratio_m_ratio_p1(a, x);

  aq = gsl_dd_mul_d(q, a);
  w = gsl_dd_sub_d(aq, E);

  /* Keep exp(w) normal whenever possible.  If it would approach
   * underflow, absorb sqrt(a/(2*pi)) into the exponent first. */
  if (w.hi > log(GSL_DBL_MIN) + 8.0)
    {
      const double ew = exp(w.hi) * exp(w.lo);
      return sqrt(a / (2.0 * M_PI)) * ew;
    }
  else
    {
      const double base = 0.5 * log(a) - GSL_SF_LOG_SQRT_2PI;
      gsl_dd z = gsl_dd_add_d(w, base);
      return exp(z.hi) * exp(z.lo);
    }
}

/*
 * Compute the scaled incomplete-gamma factor
 *
 *   rgammax(a,x) = x^a exp(-x) / Gamma(a).
 *
 * For a > 20, use the compensated large-a representation implemented
 * by gamma_inc_rgammax_large().  For moderate a, use the direct
 * representation.
 *
 * This factor is used to scale the positive P series and is also
 * suitable for reuse by upper-tail incomplete-gamma kernels.
 */
static double
gamma_inc_rgammax(const double a, const double x)
{
  if (a < 1.0 || x < 0.0)
    return GSL_NAN;

  if (x == 0.0)
    return 0.0;

  if (a > 20.0)
    return gamma_inc_rgammax_large(a, x);

  return pow(x, a) * exp(-x) / tgamma(a);
}

/**********************************************
 * INTERNAL ROUTINES (double-double routines) *
 **********************************************/

/*
 * Compute
 *
 *   (x-a)/a
 *
 * as a double-double number.  Near x=a, the subtraction x-a is exact
 * by Sterbenz's lemma, and the low component retains the residual from
 * the subsequent division by a.
 */
static gsl_dd
dd_div_diff(const double a, const double x)
{
  gsl_dd d;
  d.hi = x - a;
  d.lo = 0;
  return gsl_dd_div_d(d, a);
}
