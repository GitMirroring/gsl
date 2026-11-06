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
#include <gsl/gsl_sf_log.h>

#include "constants.h"
#include "error.h"
#include "eval.h"
#include "stirling.c"

#define GAMMA_INC_MAXITER         10000

/* Temme: largest a*phi for which phi in double precision suffices */
#define GAMMA_INC_TEMME_APHI_DD   0.25

typedef enum { TARGET_P = 0, TARGET_Q = 1 } gamma_inc_target_t;

/* Taylor kernel functions */
static int gamma_inc_P_taylor_x_sum(const double a, const double x, double * temp, double * err);
static int gamma_inc_P_taylor_x(const double a, const double x, gsl_sf_result * result);

/* P series kernel functions */
static int gamma_inc_p_series_sum(const double a, const double x, double * sum, int * nterms);
static int gamma_inc_P_series(const double a, const double x, gsl_sf_result * result);

/* Q kernel functions */
static int gamma_inc_Q_small_a(const double a, const double x, gsl_sf_result * result);
static int gamma_inc_complement(const gsl_sf_result * r, gsl_sf_result * result);
static double gamma_inc_rgamma1pm1(const double a);
static int gamma_inc_Q_cf(const double a, const double x, gsl_sf_result * result);
static int gamma_inc_Q_asymp(const double a, const double x, gsl_sf_result * result);

/* Temme kernel functions */
static int gamma_inc_use_temme(const double a, const double x);
static int gamma_inc_temme_eval(const double a, const double x, const int target, gsl_sf_result * result);

/* rgammax primitive functions */
static gsl_dd gamma_inc_log_ratio_m_ratio_p1(const double a, const double x);
static gsl_dd gamma_inc_phi_dd(const double a, const double x);
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
 * positive series.  On the upper side, evaluate Q with the Q dispatcher
 * and form P = 1 - Q when Q <= 1/2; otherwise (just above x = a) the
 * P series is still efficient and avoids the loss of an ulp of Q.
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

      /* complement only when P >= 1/2, where it loses at most half an ulp */
      if (q.val > 0.5)
        return gamma_inc_P_series(a, x, result);

      return gamma_inc_complement(&q, result);
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
 * asymptotic expansion.  For x <= a, and for a < 12, x < 40, form
 * Q = 1-P from the positive P series when P <= 1/2, and use the
 * continued fraction otherwise (1-P would lose up to an ulp of P,
 * twice an ulp of Q).  Elsewhere use the continued fraction or the
 * large-x asymptotic expansion.
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
  else if (x <= a || (a < 12.0 && x < 40.0))
    {
      gsl_sf_result p;

      status = gamma_inc_P_series(a, x, &p);
      if (status)
        return status;

      /* complement only when Q >= 1/2, where it loses at most half an ulp */
      if (p.val > 0.5)
        return gamma_inc_Q_cf(a, x, result);

      return gamma_inc_complement(&p, result);
    }
  else if (a < 12.0)
    {
      return gamma_inc_Q_asymp(a, x, result);
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

/*
 * result = 1 - r for r <= 1/2: the subtraction loses at most half an ulp
 * (both 1 - r and its ulp are >= 1/2 and 2^-53).
 */
static int
gamma_inc_complement(const gsl_sf_result * r, gsl_sf_result * result)
{
  result->val = 1.0 - r->val;
  result->err = r->err + GSL_DBL_EPSILON * result->val;
  return GSL_SUCCESS;
}

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
        GSL_ERROR ("gamma_inc: series or continued fraction failed to converge", GSL_EMAXITER);
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

  GSL_ERROR ("gamma_inc: series or continued fraction failed to converge", GSL_EMAXITER);
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
 * Compute 1/Gamma(1+a) - 1 for |a| < 0.25 from its Taylor series
 *
 *   1/Gamma(1+a) - 1 = sum_{k>=1} c_k a^k,   c_1 = Euler's constant, ...
 *
 * truncated after 19 terms (remainder below 1e-22 relative) and summed by
 * Horner's method in double precision (cf. "gam1" in DiDonato and Morris,
 * ACM TOMS 12, 1986).  The result has a relative error of about one ulp;
 * forming it from lngamma1p and expm1 instead roughly triples the error,
 * which matters because gamma_inc_Q_small_a() has Q(a,x) ~ a here.
 */
static double
gamma_inc_rgamma1pm1(const double a)
{
  static const double c[19] = {
    +5.77215664901532866e-01, -6.55878071520253902e-01, -4.20026350340952370e-02,
    +1.66538611382291479e-01, -4.21977345555443334e-02, -9.62197152787697303e-03,
    +7.21894324666309990e-03, -1.16516759185906517e-03, -2.15241674114950975e-04,
    +1.28050282388116196e-04, -2.01348547807882387e-05, -1.25049348214267063e-06,
    +1.13302723198169593e-06, -2.05633841697760707e-07, +6.11609510448141609e-09,
    +5.00200764446922295e-09, -1.18127457048702004e-09, +1.04342671169110054e-10,
    +7.78226343990507081e-12
  };
  double p = c[18];
  int k;

  for (k = 17; k >= 0; --k)
    p = p * a + c[k];

  return p * a;
}

static int
gamma_inc_Q_small_a(const double a, const double x, gsl_sf_result * result)
{
  /*
   * Q(a,x) for 0 < a <= 2, x > 0.  Use the continued fraction of
   * gamma_inc_Q_cf() for x >= 1.1, and for x >= 0.5 when a < 0.25 (see below).
   * Otherwise use the series of DiDonato and Morris (TOMS 12, 1986, "grat1"):
   *
   *   Q(a,x) = 1 - x^a (1 - J(a,x)) / Gamma(1+a)
   *          = J (1+l)(1+h) - (l + h + l h),
   *
   * with l = x^a - 1 = expm1(a log x), h = 1/Gamma(1+a) - 1 and J from
   * gamma_inc_P_taylor_x_sum().  The second form avoids the cancellation in
   * 1 - x^a/Gamma(1+a), but for small a its terms still nearly cancel as x
   * grows: Q -> a E1(x) while the terms are ~ a log(x), a*Euler and
   * a Ein(x), so the rounding errors are amplified by up to ~6 at x = 1
   * (~3 at x = 0.5).  The continued fraction is accurate there, so it is
   * used for x >= 1.1, and for x >= 0.5 when a < 0.25; it needs ~250 terms
   * at x = 0.5.  For a >= 0.25 the amplification is small.
   */
  double j, jerr, z, h;
  int status;

  if (x >= 1.1 || (a < 0.25 && x >= 0.5))
    return gamma_inc_Q_cf(a, x, result);

  status = gamma_inc_P_taylor_x_sum(a, x, &j, &jerr);
  if (status)
    return status;

  if (fabs(a) < 0.25)
    h = gamma_inc_rgamma1pm1(a);
  else
    h = gsl_sf_expm1(-gsl_sf_lngamma1p(a));   /* not 1/Gamma(1+a) - 1, which cancels */

  z = a * log(x);

  {
    const double l = gsl_sf_expm1(z);

    result->val = j * (1.0 + l) * (1.0 + h) - (l + h + l * h);
    result->err = 8.0 * GSL_DBL_EPSILON * (fabs(result->val) + fabs(j) + fabs(l) + fabs(h)) +
                  jerr + GSL_SF_TRUE_MIN;
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
  /*
   * Legendre continued fraction for x > a,
   *
   *   Q(a,x) = x^a e^-x / Gamma(a) * 1/g,
   *   g = b_0 + a_1/(b_1 + a_2/(b_2 + ...)),
   *   b_k = x + 2k + 1 - a,  a_k = -k (k - a),
   *
   * in double precision.  The modified Lentz algorithm only estimates the
   * number of terms n; g is then evaluated backward from b_N, which damps
   * the rounding errors of each step (Lentz's running product of correction
   * factors accumulates them instead: up to ~5 eps after a few terms).
   * Lentz's test |delta - 1| <= eps stops too early when convergence is
   * slow (x near 1: truncation errors up to ~7 eps), so the backward pass
   * uses N = 3n/2 + 2 terms, which makes the truncation error negligible.
   */
  const double tiny = 1.0e-300;
  double g = x + 1.0 - a;
  double C, D, delta;
  int n, k;

  if (g == 0.0)
    g = tiny;

  C = g;
  D = 0.0;

  for (n = 1; n < GAMMA_INC_MAXITER; n++)
    {
      const double an = -(double) n * ((double) n - a);
      const double bn = x + 2.0 * n + 1.0 - a;

      D = bn + an * D;
      if (D == 0.0)
        D = tiny;
      C = bn + an / C;
      if (C == 0.0)
        C = tiny;
      D = 1.0 / D;
      delta = C * D;
      g *= delta;

      if (fabs(delta - 1.0) <= GSL_DBL_EPSILON)
        break;
    }

  if (n >= GAMMA_INC_MAXITER)
    GSL_ERROR ("gamma_inc: series or continued fraction failed to converge", GSL_EMAXITER);

  /* backward evaluation */
  n += n / 2 + 2;
  g = x + 2.0 * n + 1.0 - a;
  for (k = n; k >= 1; k--)
    g = (x + 2.0 * (k - 1) + 1.0 - a) - (double) k * ((double) k - a) / g;

  result->val = gamma_inc_rgammax(a, x) / g;
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
  /* compensated (Kahan) summation: up to ~40 terms for x ~ 2.5a, whose
     rounding errors would otherwise add up to a few ulp */
  double term = 1.0, prev = 1.0;
  double sum = 1.0, comp = 0.0;
  int n;
  double r;

  for (n = 1; n < GAMMA_INC_MAXITER; n++)
    {
      double y, t;

      term *= (a - (double)n) / x;
      if (fabs(term) > prev)
        break;

      y = term - comp;
      t = sum + y;
      comp = (t - sum) - y;
      sum = t;

      prev = fabs(term);
      if (term == 0.0 || fabs(term) <= GSL_DBL_EPSILON * fabs(sum))
        break;
    }

  if (n >= GAMMA_INC_MAXITER)
    GSL_ERROR ("gamma_inc: series or continued fraction failed to converge", GSL_EMAXITER);

  r = gamma_inc_rgammax(a, x);

  result->val = (r / x) * sum;
  result->err = (n + 16.0) * GSL_DBL_EPSILON * fabs(result->val) + GSL_SF_TRUE_MIN;

  return GSL_SUCCESS;
}

/*****************************************
 * INTERNAL ROUTINES (Temme kernel)      *
 *****************************************/

/*
 * Return nonzero if the Temme uniform asymptotic expansion should be
 * used for the point (a,x).
 *
 * Define
 *
 *   s = (x-a)/a,
 *   eta^2 / 2 = phi = s - log(1+s),
 *
 * with sign(eta) = sign(x-a).  The present dispatcher uses Temme for
 *
 *   a >= 12,  |eta| <= 0.5.
 *
 * This routine is a Boolean predicate; its return value is not a GSL
 * error status.
 */

static int
gamma_inc_use_temme(const double a, const double x)
{
  if (a < 12.0 || x <= 0.0)
    return 0;

  /* eta <= 0.5 <=> phi <= 1/8; phi in double is accurate enough here */
  return -gsl_sf_log1pmx((x - a) / a) <= 0.125;
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
 * upper ratio Q(a,x).  a phi is formed in double-double only when it is
 * larger than GAMMA_INC_TEMME_APHI_DD (see below).
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
  const double phi = GSL_MAX(-gsl_sf_log1pmx((x - a) / a), 0.0);
  const double eta = (x < a) ? -sqrt(2.0 * phi) : sqrt(2.0 * phi);
  gsl_dd aphi, z;
  double beta[27], sum, scale, e, R, value, last;
  int n;

  /*
   * The leading term erfc(z)/2, z^2 = a phi, has relative error ~ |d(a phi)|,
   * so a phi is needed to absolute accuracy.  The double phi (relative
   * error ~ 2 eps, from rounding s) suffices while a phi is small;
   * otherwise form phi = s - log1p(s) in double-double from the exact
   * difference x - a.
   */
  aphi.hi = a * phi;
  aphi.lo = 0.0;

  if (aphi.hi > GAMMA_INC_TEMME_APHI_DD)
    aphi = gsl_dd_mul_d(gamma_inc_phi_dd(a, x), a);

  z = gsl_dd_sqrt(aphi);
  if (x < a)
    {
      z.hi = -z.hi;
      z.lo = -z.lo;
    }

  for (n = N + 1; n <= N + 2; ++n)
    beta[n] = 0;

  for (n = N; n >= 0; --n)
    beta[n] = ((n + 2.0) * beta[n + 2]) / a + gamma_inc_temme_d[n + 1];

  sum = beta[N];
  for (n = N - 1; n >= 0; --n)
    sum = sum * eta + beta[n];

  sum *= a / (a + beta[1]);

  /* e = exp(-a phi) = exp(-z^2), shared by R and the erfc correction */
  e = exp(-aphi.hi) * (1.0 - aphi.lo);
  scale = e / sqrt(2.0 * M_PI * a);
  R = scale * sum;

  /*
   * erfc(z.hi + z.lo)/2 = erfc(z.hi)/2 - exp(-z^2) z.lo / sqrt(pi) + O(z.lo^2);
   * the correction matters in the tail, where erfc amplifies the rounding
   * of its argument.
   */
  if (target == TARGET_Q)
    value = 0.5 * gsl_sf_erfc(z.hi) - e * z.lo / M_SQRTPI + R;
  else
    value = 0.5 * gsl_sf_erfc(-z.hi) + e * z.lo / M_SQRTPI - R;

  last = fabs(beta[N] * pow(fabs(eta), (double) N));

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

/*
 * phi = s - log(1+s), s = (x-a)/a, in double-double.  For x in [a/2, 2a]
 * the difference x - a is exact (Sterbenz), and dd_div_diff() keeps the
 * rounding of the division.
 */
static gsl_dd
gamma_inc_phi_dd(const double a, const double x)
{
  const gsl_dd s = dd_div_diff(a, x);
  return gsl_dd_sub(s, gsl_dd_log1p(s));
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
 * gamma_stirling_delta() (stirling.c).
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
  const double E = gamma_stirling_delta(a);
  gsl_dd q, aq, w;

  /* Near x=a use the exact-difference representation.  Away from it,
   * form the quotient x/a directly as a compensated pair. */
  if (x >= 0.5 * a && x <= 2.0 * a)
    q = gsl_dd_mul_d(gamma_inc_phi_dd(a, x), -1.0);
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
  if (a <= 0.0 || x < 0.0)
    return GSL_NAN;

  if (x == 0.0)
    return 0.0;

  if (a > 20.0)
    return gamma_inc_rgammax_large(a, x);

  if (x > 700.0)
    {
      /* exp(-x) alone would lose precision or underflow before the product
         does: form log(x^a e^-x) = a log(x) - x in double-double */
      gsl_dd lr = gsl_dd_sub_d(gsl_dd_mul_d(gsl_dd_log((gsl_dd){x, 0.0}), a), x);
      gsl_dd er;

      if (lr.hi <= GSL_SF_LOG_TRUE_MIN)
        return 0.0;

      er = gsl_dd_exp(lr);
      return (er.hi + er.lo) / gsl_new_sf_gamma(a);
    }

  if (a < 0.25)
    {
      /* 1/Gamma(a) = a (1 + h) is more accurate than dividing by Gamma(a) */
      return pow(x, a) * exp(-x) * (a + a * gamma_inc_rgamma1pm1(a));
    }

  return pow(x, a) * exp(-x) / gsl_new_sf_gamma(a);
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
