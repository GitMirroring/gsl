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

#include "check.h"
#include "error.h"
#include "eval.h"

#include "chebyshev.h"
#include "cheb_eval.c"

#include "gamma_kernels.c"
#include "gamma_fallback.c"
#include "gammastar_kernels.c"

#define TAYLORCOEFF_POW_NMAX 100

#ifdef HAVE_TGAMMA
static int gamma_system_e(const double x, gsl_sf_result * result);
#endif

static int taylorcoeff_small(const int n, const double x, gsl_sf_result * result);

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

/*
 * Taylor coefficient
 *
 * Computes
 *
 *   x^n / n!
 *
 * for integer n >= 0 and x >= 0.
 *
 * The cases n = 0, n = 1, and x = 0 are returned exactly.  For
 * 2 <= n <= TAYLORCOEFF_POW_NMAX, a scaled power representation is
 * used to reduce accumulated rounding error and avoid intermediate
 * overflow in x^n.  Larger orders use the general product formulation,
 * with logarithmic tests for underflow and overflow.
 */

int
gsl_new_sf_taylorcoeff_e(const int n, const double x, gsl_sf_result * result)
{
  if (x < 0.0 || n < 0)
    DOMAIN_ERROR(result);
  else if (n == 0)
    {
      result->val = 1.0;
      result->err = 0.0;
      return GSL_SUCCESS;
    }
  else if (n == 1)
    {
      result->val = x;
      result->err = 0.0;
      return GSL_SUCCESS;
    }
  else if (x == 0.0)
    {
      result->val = 0.0;
      result->err = 0.0;
      return GSL_SUCCESS;
    }
  else if (n <= TAYLORCOEFF_POW_NMAX)
    {
      return taylorcoeff_small(n, x, result);
    }
  else
    {
      const double log2pi = M_LNPI + M_LN2;
      const double ln_test = n*(log(x) + 1.0) +
                             1.0 -
                             (n + 0.5)*log(n + 1.0) +
                             0.5*log2pi;

      if (ln_test < GSL_LOG_DBL_MIN + 1.0)
        UNDERFLOW_ERROR(result);
      else if (ln_test > GSL_LOG_DBL_MAX - 1.0)
        OVERFLOW_ERROR(result);
      else
        {
          double product = 1.0;
          int k;

          for (k = 1; k <= n; ++k)
            product *= (x / (double) k);

          result->val = product;
          result->err = n * GSL_DBL_EPSILON * product;
          CHECK_UNDERFLOW(result);
          return GSL_SUCCESS;
        }
    }
}

double gsl_new_sf_gamma(const double x) { EVAL_RESULT(gsl_new_sf_gamma_e(x, &result)); }
double gsl_new_sf_gammainv(const double x) { EVAL_RESULT(gsl_new_sf_gammainv_e(x, &result)); }
double gsl_new_sf_gammastar(const double x) { EVAL_RESULT(gsl_new_sf_gammastar_e(x, &result)); }
double gsl_new_sf_taylorcoeff(const int n, const double x) { EVAL_RESULT(gsl_new_sf_taylorcoeff_e(n, x, &result)); }

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

/********************************************
 * INTERNAL ROUTINES FOR TAYLORCOEFF        *
 ********************************************/

/* table of 1/n! */
static const double taylorcoeff_inv_fact[TAYLORCOEFF_POW_NMAX + 1] = {
  1.00000000000000000000e+00,
  1.00000000000000000000e+00,
  5.00000000000000000000e-01,
  1.66666666666666657415e-01,
  4.16666666666666643537e-02,
  8.33333333333333321769e-03,
  1.38888888888888894189e-03,
  1.98412698412698412526e-04,
  2.48015873015873015658e-05,
  2.75573192239858925110e-06,
  2.75573192239858882758e-07,
  2.50521083854417202239e-08,
  2.08767569878681001866e-09,
  1.60590438368216133409e-10,
  1.14707455977297245073e-11,
  7.64716373181981640551e-13,
  4.77947733238738525345e-14,
  2.81145725434552059811e-15,
  1.56192069685862252711e-16,
  8.22063524662432949554e-18,
  4.11031762331216484407e-19,
  1.95729410633912625952e-20,
  8.89679139245057407789e-22,
  3.86817017063068412615e-23,
  1.61173757109611838590e-24,
  6.44695028438447358950e-26,
  2.47959626322479758961e-27,
  9.18368986379554600539e-29,
  3.27988923706983775898e-30,
  1.13099628864477158818e-31,
  3.76998762881590538515e-33,
  1.21612504155351789377e-34,
  3.80039075485474341802e-36,
  1.15163356207719508905e-37,
  3.38715753552116179633e-39,
  9.67759295863189067186e-41,
  2.68822026628663633314e-42,
  7.26546017915307135902e-44,
  1.91196320504028195296e-45,
  4.90246975651354351900e-47,
  1.22561743912838584936e-48,
  2.98931082714240461486e-50,
  7.11740673129143899398e-52,
  1.65521086774219514723e-53,
  3.76184288123226156652e-55,
  8.35965084718280448735e-57,
  1.81731540156147899030e-58,
  3.86662851396059404326e-60,
  8.05547607075123643865e-62,
  1.64397470831657907432e-63,
  3.28794941663315801371e-65,
  6.44695964045717239503e-67,
  1.23979993085714861670e-68,
  2.33924515256065763988e-70,
  4.33193546770492181599e-72,
  7.87624630491803920845e-74,
  1.40647255444964978852e-75,
  2.46749570956078931125e-77,
  4.25430294751860221945e-79,
  7.21068296189593648888e-81,
  1.20178049364932260184e-82,
  1.97013195680216822973e-84,
  3.17763218839059420696e-86,
  5.04386061649300674859e-88,
  7.88103221327032304467e-90,
  1.21246649434928039085e-91,
  1.83707044598375812925e-93,
  2.74189618803546009796e-95,
  4.03220027652273531588e-97,
  5.84376851669961609380e-99,
  8.34824073814230950307e-101,
  1.17580855466793077909e-102,
  1.63306743703879288642e-104,
  2.23707868087505865569e-106,
  3.02307929847980900440e-108,
  4.03077239797307883798e-110,
  5.30364789206984013894e-112,
  6.88785440528550633671e-114,
  8.83058257087885543380e-116,
  1.11779526213656394330e-117,
  1.39724407767070487409e-119,
  1.72499268848235168313e-121,
  2.10364962010042899656e-123,
  2.53451761457883002149e-125,
  3.01728287449860730512e-127,
  3.54974455823365539682e-129,
  4.12760995143448317087e-131,
  4.74437925452239458412e-133,
  5.39134006195726657286e-135,
  6.05768546287333344247e-137,
  6.73076162541481493607e-139,
  7.39644134660968651243e-141,
  8.03961015935835482564e-143,
  8.64474210683694059625e-145,
  9.19653415620951124896e-147,
  9.68056226969422298882e-149,
  1.00839190309314822800e-150,
  1.03957928153932807448e-152,
  1.06079518524421217556e-154,
  1.07151028812546698776e-156,
  1.07151028812546693116e-158,
};

/*
 * Evaluate x^n/n! for 2 <= n <= TAYLORCOEFF_POW_NMAX.
 *
 * Write
 *
 *   x = m 2^e,        0.5 <= m < 1,
 *
 * and evaluate
 *
 *   x^n/n! = [m^n / n!] 2^(n e).
 *
 * Since m < 1, pow(m,n) cannot overflow.  The reciprocal factorial is
 * taken from taylorcoeff_inv_fact[], and ldexp() applies the final
 * power-of-two scaling.  This avoids both an intermediate overflow in
 * x^n and the rounding error accumulated by n successive products.
 *
 * Infinite and zero results are reported as overflow and underflow,
 * respectively.  Nonzero subnormal results are retained.
 */
static int
taylorcoeff_small(const int n, const double x, gsl_sf_result * result)
{
  int e;
  double m, p, y;

  m = frexp(x, &e);
  p = pow(m, (double) n);
  p *= taylorcoeff_inv_fact[n];
  y = ldexp(p, n * e);

  if (gsl_isinf(y))
    {
      result->val = y;
      result->err = GSL_POSINF;
      return GSL_EOVRFLW;
    }
  else if (y == 0.0)
    {
      result->val = 0.0;
      result->err = GSL_SF_TRUE_MIN;
      return GSL_EUNDRFLW;
    }
  else
    {
      result->val = y;

      if (fabs(y) < GSL_DBL_MIN)
        result->err = 2.0 * GSL_SF_TRUE_MIN;
      else
        result->err = 8.0 * GSL_DBL_EPSILON * fabs(y);

      return GSL_SUCCESS;
    }
}
