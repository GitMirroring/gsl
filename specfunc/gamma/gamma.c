/* specfunc/gamma/gamma.c
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
#include <gsl/gsl_sf_gamma.h>

#include "error.h"
#include "eval.h"

#include "chebyshev.h"
#include "cheb_eval.c"

#include "gamma_kernels.c"
#include "gamma_fallback.c"
#include "gammastar_kernels.c"

#ifdef HAVE_TGAMMA
static int gamma_system_e(const double x, gsl_sf_result * result);
#endif

/*
 * Gamma function
 *
 * Computes
 *
 *   Gamma(x)
 *
 * over the real line.  Negative integers and zero are poles and return
 * GSL_EDOM.  Results which overflow or underflow the double-precision
 * range return GSL_EOVRFLW or GSL_EUNDRFLW, respectively.
 *
 * If the system tgamma() function is available, it is used directly.
 * Otherwise gamma_fallback_e() provides a portable implementation using
 * the Gamma kernels in gamma_kernels.c.
 */

int
gsl_new_sf_gamma_e(const double x, gsl_sf_result * result)
{
#ifdef HAVE_TGAMMA

  /* use system tgamma() */
  return gamma_system_e(x, result);

#else

  return gamma_fallback_e(x, result);

#endif
}

/*
 * Reciprocal gamma function
 *
 * Computes
 *
 *   1 / Gamma(x)
 *
 * over the real line.  Several representations are used according to
 * the argument:
 *
 *   x < -10
 *     The reflection formula
 *
 *       1/Gamma(x) = (-x) sin(pi*x) Gamma(-x) / pi
 *
 *     is evaluated using sin(pi*x) and a scaled positive Gamma kernel.
 *     Negative integer poles of Gamma are returned as exact zeros.
 *
 *   -10 <= x < -1/2
 *     The argument is shifted into the central interval and the
 *     recurrence relation is used.
 *
 *   -1/2 <= x <= 2
 *     A rational approximation to the auxiliary function H(x) is used
 *     to evaluate reciprocal Gamma near unity without loss of accuracy.
 *
 *   x > 2
 *     A scaled positive Gamma representation is inverted directly,
 *     avoiding intermediate overflow of Gamma(x).
 */

int
gsl_new_sf_gammainv_e(const double x, gsl_sf_result * result)
{
  int status;
  gsl_sf_result tmp;
  gsl_sf_result gm;
  double r, base, p;
  double d, ed, term;
  double ebase, eprod;
  double qmant, eqmant;
  int n, k;
  int ge, qexp;

  if (x < -10.0)
    {
      /* Reflection: 1/Gamma(x) = (-x) sin(pi*x) Gamma(-x) / pi */
      status = gsl_sf_sin_pi_e(x, &tmp);
      if (status)
        return status;

      /* Gamma poles: 1/Gamma(-n) = 0 exactly */
      if (tmp.val == 0.0)
        {
          result->val = 0.0;
          result->err = 0.0;
          return GSL_SUCCESS;
        }

      /* The remaining nonintegral negative tail eventually overflows */
      if (-x > GAMMA_SCALED_X_MAX)
        {
          result->val = tmp.val < 0.0 ? GSL_NEGINF : GSL_POSINF;
          result->err = GSL_POSINF;
          return GSL_EOVRFLW;
        }

      status = gamma_2_inf_scaled(-x, &gm, &ge);
      if (status)
        return status;

      /*
       * qmant = (-x) sin(pi*x) Gamma(-x) / pi, with the Gamma
       * factor represented in scaled form.
       */
      qmant = ((-x) * tmp.val / M_PI) * gm.val;
      qmant = frexp(qmant, &qexp);
      qexp += ge;

      result->val = ldexp(qmant, qexp);

      if (gsl_isinf(result->val))
        {
          result->err = GSL_POSINF;
          return GSL_EOVRFLW;
        }

      /* Propagate the sin(pi*x) and scaled Gamma errors */
      eqmant = fabs((-x) / M_PI) * tmp.err * fabs(gm.val) +
               fabs(((-x) * tmp.val / M_PI)) * gm.err +
               3.0 * GSL_DBL_EPSILON * fabs(qmant);

      result->err = ldexp(eqmant, qexp);

      return GSL_SUCCESS;
    }
  else if (x < -0.5)
    {
      /* Moderate negative x: 1/Gamma(x) = x...(x+n-1) / Gamma(x+n) */
      n = (int) ceil(-x - 0.5);
      if (n < 1)
        n = 1;

      r = x + (double) n;

      status = gamma1p_H(r, &tmp);
      if (status)
        return status;

      d = 1.0 + r * (r - 1.0) * tmp.val;
      ed = fabs(r * (r - 1.0)) * tmp.err +
           2.0 * GSL_DBL_EPSILON *
           fabs(r * (r - 1.0) * tmp.val) +
           GSL_DBL_EPSILON * fabs(d);

      base = r / d;
      ebase = fabs(base) * (GSL_DBL_EPSILON + ed / fabs(d));

      p = 1.0;
      for (k = 0; k < n; ++k)
        p *= x + (double) k;

      eprod = (double) n * GSL_DBL_EPSILON * fabs(p);

      result->val = p * base;
      result->err = fabs(base) * eprod +
                    fabs(p) * ebase +
                    GSL_DBL_EPSILON * fabs(result->val);

      return GSL_SUCCESS;
    }
  else if (x <= 1.0)
    {
      /* Central H representation: 1/Gamma(x) = x / [1 + x(x-1) H(x)] */
      status = gamma1p_H(x, &tmp);
      if (status)
        return status;

      d = 1.0 + x * (x - 1.0) * tmp.val;
      ed = fabs(x * (x - 1.0)) * tmp.err +
           2.0 * GSL_DBL_EPSILON *
           fabs(x * (x - 1.0) * tmp.val) +
           GSL_DBL_EPSILON * fabs(d);

      result->val = x / d;
      result->err = fabs(result->val) * (GSL_DBL_EPSILON + ed / fabs(d));

      return GSL_SUCCESS;
    }
  else if (x <= 2.0)
    {
      /* 1/Gamma(x) = 1 / [1 + r(r-1)H(r)], r = x-1 */
      r = x - 1.0;

      status = gamma1p_H(r, &tmp);
      if (status)
        return status;

      term = r * (r - 1.0) * tmp.val;
      d = 1.0 + term;
      ed = fabs(r * (r - 1.0)) * tmp.err +
           2.0 * GSL_DBL_EPSILON * fabs(term) +
           GSL_DBL_EPSILON * fabs(d);

      result->val = 1.0 / d;
      result->err = fabs(result->val) * (GSL_DBL_EPSILON + ed / fabs(d));

      return GSL_SUCCESS;
    }
  else if (x <= GAMMA_SCALED_X_MAX)
    {
      status = gamma_2_inf_scaled(x, &gm, &ge);
      if (status)
        return status;

      result->val = ldexp(1.0 / gm.val, -ge);

      if (fabs(result->val) < GSL_DBL_MIN)
        result->err = GSL_SF_TRUE_MIN;
      else
        result->err = fabs(result->val) * (GSL_DBL_EPSILON + gm.err / fabs(gm.val));

      return GSL_SUCCESS;
    }
  else
    {
      /* 1/Gamma(x) has underflowed to zero. */
      result->val = 0.0;
      result->err = GSL_SF_TRUE_MIN;
      return GSL_SUCCESS;
    }
}

/*
 * Scaled Gamma function
 *
 * Computes
 *
 *   Gamma*(x) = Gamma(x) / [sqrt(2*pi/x) (x/e)^x]
 *
 * for x > 0.  The implementation uses several representations according
 * to the argument:
 *
 *   0 < x < 0.4
 *     Remove the square-root singularity at the origin and evaluate
 *
 *       Gamma*(x) = exp[log(Gamma(1+x)) + x(1-log(x))] / [sqrt(2*pi) sqrt(x)].
 *
 *     log(Gamma(1+x)) is evaluated by the small-argument rational kernel,
 *     and expm1() is used to reconstruct the exponential accurately.
 *
 *   0.4 <= x < 2
 *     Use a Chebyshev approximation to Gamma*(x).
 *
 *   2 <= x < 16
 *     Use a second Chebyshev approximation to Gamma*(x).
 *
 *   x >= 16
 *     Use the direct Stirling asymptotic expansion for Gamma*(x),
 *     truncated at its least term and evaluated with compensated
 *     summation.
 *
 * Gamma*(x) is defined here only for x > 0; nonpositive arguments
 * return GSL_EDOM.
 */

int
gsl_new_sf_gammastar_e(const double x, gsl_sf_result * result)
{
  if (x <= 0.0)
    DOMAIN_ERROR(result);

  if (x < GAMMASTAR_X_SMALL)
    return gammastar_small_e(x, result);
  else if (x < GAMMASTAR_X_MID)
    return gammastar_mid1_e(x, result);
  else if (x < GAMMASTAR_X_ASYM)
    return gammastar_mid2_e(x, result);
  else
    return gammastar_asym_e(x, result);
}

double gsl_new_sf_gamma(const double x) { EVAL_RESULT(gsl_new_sf_gamma_e(x, &result)); }
double gsl_new_sf_gammainv(const double x) { EVAL_RESULT(gsl_new_sf_gammainv_e(x, &result)); }
double gsl_new_sf_gammastar(const double x) { EVAL_RESULT(gsl_new_sf_gammastar_e(x, &result)); }

/********************************************
 * INTERNAL ROUTINES FOR GAMMA              *
 ********************************************/

#ifdef HAVE_TGAMMA

/*
 * Evaluate Gamma(x) using the system tgamma() function.
 *
 * The system function supplies only the value, so an absolute error
 * estimate is constructed here.  For normal finite results the estimate
 * is
 *
 *   12 eps |Gamma(x)|,
 *
 * while nonzero subnormal results use 2*GSL_SF_TRUE_MIN.  Exact zero
 * results are reported as underflow and infinite results as overflow.
 * Gamma poles at zero and the negative integers are detected before
 * calling tgamma().
 */

static int
gamma_system_e(const double x, gsl_sf_result * result)
{
  if (x == 0.0 || (x < 0.0 && x == floor(x)))
    DOMAIN_ERROR(result);

  result->val = tgamma(x);

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
      result->err = 2.0 * GSL_SF_TRUE_MIN;
      return GSL_SUCCESS;
    }
  else
    {
      result->err = 12.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
}

#endif /* HAVE_TGAMMA */
