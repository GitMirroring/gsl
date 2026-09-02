/* specfunc/gamma/lngamma.c
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
#include <gsl/gsl_sf_sincos_pi.h>

#include "error.h"
#include "eval.h"

#ifdef HAVE_LGAMMA
static int lngamma_system_e(const double x, gsl_sf_result * result);
#endif

#ifdef HAVE_LGAMMA_R
static int lngamma_system_sgn_e(const double x, gsl_sf_result * result, double * sgn);
#endif

static int lngamma_near1(const double x, gsl_sf_result * result);
static int lngamma_near2(const double x, gsl_sf_result * result);
static int lngamma_middle(const double x, gsl_sf_result * result);
static int lngamma_reduce(const double x, const double lo,
                          const double hi, const int hi_closed, double * y);
static int lngamma_recurrence(const double x, gsl_sf_result * result);
static int lngamma_stirling(const double x, gsl_sf_result * result);

static void lngamma_neg_mul_split(double *hi, double *lo, const double x, const double y);
static double lgamma_product(const double t, const double x, const double x_eps, const int n);
static double lngamma_neg_sinpi(const double x);
static double lngamma_neg_cospi(const double x);
static double lngamma_neg_cotpi(const double x);
static int lngamma_neg_large(const double x, gsl_sf_result * result);
static int lngamma_neg(const double x, gsl_sf_result * result);
static int lngamma_neg_minus3_minus2(const double x, gsl_sf_result * result);
static int lngamma_neg_minus2_minus1(const double x, gsl_sf_result * result);

static int lngamma_1p_near0(const double x, gsl_sf_result * result);
static int lngamma_1p_near1_y(const double y, gsl_sf_result * result);
static int lngamma_1p_near1(const double x, gsl_sf_result * result);
static int lngamma_1p_recurrence(const double x, gsl_sf_result * result);

/* mathematical constants */
#define LNGAMMA_M_PI2_12                   0.8224670334241132182362075833230125946094749506033992188677791146850037 /* pi^2 / 12 */
#define LNGAMMA_M_PI2_12_MINUS_HALF        0.3224670334241132182362075833230125946094749506033992188677791146850037 /* pi^2 / 12 - 1/2 */
#define LNGAMMA_M_1_MINUS_EULER            0.4227843350984671393934879099175975689578406640600764011942327651151323 /* 1 - M_EULER */
#define LNGAMMA_M_LNGAMMA_3_2             -0.1207822376352452223455184457816472122518527279025994683638684737573247 /* log(Gamma(3/2)) */
#define LNGAMMA_M_ABS_LNGAMMA_3_2          0.1207822376352452223455184457816472122518527279025994683638684737573247 /* |log(Gamma(3/2))| */
#define LNGAMMA_M_LN_SQRT_2PI              0.9189385332046727417803297364056176398613974736377834128171515404827656 /* log(sqrt(2pi)) */

#define LNGAMMA_STIRLING_SCALE             (2.0 * LNGAMMA_X_RECUR_STIRLING * LNGAMMA_X_RECUR_STIRLING)

/* lngamma_1p region boundaries */
#define LNGAMMA_1P_X_RECUR0_RECUR1         (-0.6)
#define LNGAMMA_1P_X_RECUR1_NEAR0          (-0.4)
#define LNGAMMA_1P_X_NEAR0_NEAR1           ( 0.4)
#define LNGAMMA_1P_X_NEAR1_LNGAMMA         ( 2.0)

/* lngamma region boundaries */
#define LNGAMMA_X_NEG0_MIN                 (-28.0)
#define LNGAMMA_X_NEG0_MAX                 (-3.0)
#define LNGAMMA_X_NEG1_MAX                 (-2.0)
#define LNGAMMA_X_NEG2_MAX                 (-1.0)
#define LNGAMMA_X_SMALL_RECUR0_RECUR1      (LNGAMMA_1P_X_NEAR0_NEAR1)
#define LNGAMMA_X_SMALL_RECUR1_NEAR1       (0.8)
#define LNGAMMA_X_NEAR1_MIDDLE             (1.23395)
#define LNGAMMA_X_MIDDLE_NEAR2             (1.67895)
#define LNGAMMA_X_NEAR2_RECUR              (2.3)
#define LNGAMMA_X_RECUR_STIRLING           (8.0)

int
gsl_sf_lngamma_e(const double x, gsl_sf_result * result)
{
#ifdef HAVE_LGAMMA

  /* use system lgamma() */
  return lngamma_system_e(x, result);

#else

  /* use local fallback function */
  return gsl_sf_lngamma_local_e(x, result);

#endif
}

int
gsl_sf_lngamma_sgn_e(const double x, gsl_sf_result * result, double * sgn)
{
#ifdef HAVE_LGAMMA_R

  /* use system lgamma_r() */
  return lngamma_system_sgn_e(x, result, sgn);

#else

  /* use local fallback function */
  return gsl_sf_lngamma_sgn_local_e(x, result);

#endif
}

/*
 * Local fallback function when lgamma() is not available.
 *
 * Compute log|Gamma(x)|.
 *
 * The positive domain is divided among three direct rational
 * approximations, a recurrence region, and a large-x Stirling
 * approximation:
 *
 *   x <= LNGAMMA_X_NEAR1_MIDDLE
 *       near x=1 rational approximation
 *
 *   LNGAMMA_X_NEAR1_MIDDLE < x <= LNGAMMA_X_MIDDLE_NEAR2
 *       rational approximation centered at x=3/2
 *
 *   LNGAMMA_X_MIDDLE_NEAR2 < x < LNGAMMA_X_NEAR2_RECUR
 *       near x=2 rational approximation
 *
 *   LNGAMMA_X_NEAR2_RECUR <= x < LNGAMMA_X_RECUR_STIRLING
 *       recurrence to one of the direct approximation intervals
 *
 *   x >= LNGAMMA_X_RECUR_STIRLING
 *       Stirling approximation
 *
 * The treatment of x below the direct near x=1 interval is handled
 * separately as part of the complete lngamma implementation.
 */
int
gsl_sf_lngamma_local_e(const double x, gsl_sf_result * result)
{
  int status;
  gsl_sf_result r;

  if (x == 0.0)
    {
      DOMAIN_ERROR(result);
    }

  if (x <= LNGAMMA_X_NEG0_MIN)
    {
      return lngamma_neg_large(x, result);
    }
  else if (x > LNGAMMA_X_NEG0_MIN && x < LNGAMMA_X_NEG0_MAX)
    {
      return lngamma_neg(x, result);
    }
  else if (x >= LNGAMMA_X_NEG0_MAX && x < LNGAMMA_X_NEG1_MAX)
    {
      return lngamma_neg_minus3_minus2(x, result);
    }
  else if (x >= LNGAMMA_X_NEG1_MAX && x <= LNGAMMA_X_NEG2_MAX)
    {
      return lngamma_neg_minus2_minus1(x, result);
    }
  else if (x > LNGAMMA_X_NEG2_MAX && x < -LNGAMMA_X_SMALL_RECUR0_RECUR1)
    {
      const double lx = log(-x);

      status = lngamma_1p_recurrence(x, &r);
      if (status)
        return status;

      result->val = r.val - lx;
      result->err = r.err + GSL_DBL_EPSILON * (fabs(r.val) +
                                               2.0 * (fabs(lx) + fabs(result->val)));

      return GSL_SUCCESS;
    }
  else if (x < 0.0)
    {
      const double lx = log(-x);

      status = lngamma_1p_near0(x, &r);
      if (status)
        return status;

      result->val = r.val - lx;
      result->err = r.err + GSL_DBL_EPSILON * (fabs(r.val) +
                                               2.0 * (fabs(lx) + fabs(result->val)));

      return GSL_SUCCESS;
    }
  else if (x <= LNGAMMA_X_SMALL_RECUR0_RECUR1)
    {
      const double lx = log(x);

      status = lngamma_1p_near0(x, &r);
      if (status)
        return status;

      result->val = r.val - lx;
      result->err = r.err + GSL_DBL_EPSILON * (fabs(r.val) +
                                               2.0 * (fabs(lx) + fabs(result->val)));

      return GSL_SUCCESS;
    }
  else if (x <= LNGAMMA_X_SMALL_RECUR1_NEAR1)
    {
      const double lx = log(x);

      status = lngamma_1p_near1(x, &r);
      if (status)
        return status;

      result->val = r.val - lx;
      result->err = r.err + GSL_DBL_EPSILON * (fabs(r.val) +
                                               2.0 * (fabs(lx) + fabs(result->val)));

      return GSL_SUCCESS;
    }
  else if (x <= LNGAMMA_X_NEAR1_MIDDLE)
    return lngamma_near1(x, result);
  else if (x <= LNGAMMA_X_MIDDLE_NEAR2)
    return lngamma_middle(x, result);
  else if (x < LNGAMMA_X_NEAR2_RECUR)
    return lngamma_near2(x, result);
  else if (x < LNGAMMA_X_RECUR_STIRLING)
    return lngamma_recurrence(x, result);
  else
    return lngamma_stirling(x, result);
}

int
gsl_sf_lngamma_sgn_local_e(const double x, gsl_sf_result * result, double * sgn)
{
  if (x > 0.0)
    {
      *sgn = 1.0;
      return gsl_sf_lngamma_local_e(x, result);
    }
  else
    {
      gsl_sf_result rs;
      int status;

      status = gsl_sf_sin_pi_e(x, &rs);
      if (status)
        {
          result->val = GSL_NAN;
          result->err = GSL_NAN;
          *sgn = GSL_NAN;
          return status;
        }

      if (rs.val == 0.0)
        {
          *sgn = 0.0;
          DOMAIN_ERROR(result);
        }

      *sgn = (rs.val < 0.0) ? -1.0 : 1.0;

      return gsl_sf_lngamma_local_e(x, result);
    }
}

/*
 * Compute log(Gamma(1+x)) for x > -1.
 *
 * Dedicated zero-preserving rational approximations are used around
 * the zeros at x=0 and x=1.  For negative arguments away from x=0,
 * recurrence is used:
 *
 *   log(Gamma(1+x)) = log(Gamma(2+x)) - log(1+x).
 *
 * The regions are:
 *
 *   -1 < x < LNGAMMA_1P_X_RECUR0_RECUR1
 *       recurrence through the near x=0 kernel, with t=1+x
 *
 *   LNGAMMA_1P_X_RECUR0_RECUR1
 *       <= x < LNGAMMA_1P_X_RECUR1_NEAR0
 *       recurrence through the near x=1 kernel, supplying its local
 *       coordinate directly to avoid forming (1+x)-1
 *
 *   LNGAMMA_1P_X_RECUR1_NEAR0
 *       <= x < LNGAMMA_1P_X_NEAR0_NEAR1
 *       near x=0 rational approximation
 *
 *   LNGAMMA_1P_X_NEAR0_NEAR1
 *       <= x <= LNGAMMA_1P_X_NEAR1_LNGAMMA
 *       near x=1 rational approximation
 *
 *   x > LNGAMMA_1P_X_NEAR1_LNGAMMA
 *       evaluate log(Gamma(1+x)) with the general lngamma routine
 *
 * The function has a domain error for x <= -1.
 */
int
gsl_sf_lngamma_1p_e(const double x, gsl_sf_result * result)
{
  if (x <= -1.0)
    {
      DOMAIN_ERROR(result);
    }
  else if (x < LNGAMMA_1P_X_RECUR1_NEAR0)
    {
      return lngamma_1p_recurrence(x, result);
    }
  else if (x < LNGAMMA_1P_X_NEAR0_NEAR1)
    {
      return lngamma_1p_near0(x, result);
    }
  else if (x <= LNGAMMA_1P_X_NEAR1_LNGAMMA)
    {
      return lngamma_1p_near1(x, result);
    }
  else
    {
      return gsl_sf_lngamma_e(1.0 + x, result);
    }
}

double gsl_sf_lngamma(const double x) { EVAL_RESULT(gsl_sf_lngamma_e(x, &result)); }
double gsl_sf_lngamma_1p(const double x) { EVAL_RESULT(gsl_sf_lngamma_1p_e(x, &result)); }

/********************************************
 * INTERNAL ROUTINES FOR LNGAMMA            *
 ********************************************/

#ifdef HAVE_LGAMMA

static int
lngamma_system_e(const double x, gsl_sf_result * result)
{
  if (x <= 0.0 && x == floor(x))
    {
      DOMAIN_ERROR(result);
    }

  result->val = lgamma(x);
  result->err = 8.0 * GSL_DBL_EPSILON * (fabs(result->val) + 1.0);

  return GSL_SUCCESS;
}

#endif /* HAVE_LGAMMA */

#ifdef HAVE_LGAMMA_R

static int
lngamma_system_sgn_e(const double x, gsl_sf_result * result, double * sgn)
{
  int signp;

  if (x <= 0.0 && x == floor(x))
    {
      *sgn = 0.0;
      DOMAIN_ERROR(result);
    }

  result->val = lgamma_r(x, &signp);
  result->err = 8.0 * GSL_DBL_EPSILON * (fabs(result->val) + 1.0);

  *sgn = signp;

  return GSL_SUCCESS;
}

#endif /* HAVE_LGAMMA_R */

/*
 * Approximate log(Gamma(x)) near x = 1.
 *
 * With z = x - 1, write
 *
 *   log(Gamma(x)) = z [-gamma + z {c2 + z P(z)/Q(z)}],
 *
 * where
 *
 *   c2 = pi^2/12.
 *
 * The outer form imposes the exact zero at x = 1 and the first two
 * Taylor coefficients there; P/Q is a [7/7] rational approximation
 * to the remaining smooth part.
 *
 * Domain: LNGAMMA_X_SMALL_RECUR1_NEAR1 <= x <= LNGAMMA_X_NEAR1_MIDDLE.
 */
static int
lngamma_near1(const double x, gsl_sf_result * result)
{
  static const double P[8] = { -0.4006856343865314284665793871704833302549954307801662939307571851139461,
                               -1.306579436585641548712271782113446677915081517680219793781534716006438,
                               -1.662972285449665136059155025049677882052643163559739240911916418853459,
                               -1.040461675178213481420397620528836702514774194218551694314378983049904,
                               -0.3296450906559580044090939264629092254215878806975628484120634080197471,
                               -0.04834535578151998623851727447506127039926628170216805759508945202180571,
                               -0.002415967872144939376276310532721770905379280201192752514704711728624949,
                               -0.000001120053248848791893771478545347576393272981023720311321954018480738749 };
  static const double Q[8] = { 1.0,
                               3.93615370670359242696260921708589899866029453565655118501003130739336,
                               6.290802989071669273822419819203842432428094145524342671886642525321006,
                               5.230754050938657980597360934640112327614339705296477127506494147431803,
                               2.405173126777027686859123109248843641215723399353764394224946370879006,
                               0.5977607965221533239754934298461035048188750053944783128477785554124827,
                               0.0718115776797838899005431363233059202747797313176710166776869124604167,
                               0.003037457442975794171971336970639034188704699232607971321950957570313994 };
  const double z = x - 1.0;
  const double p = P[0] + z * (P[1] + z * (P[2] + z * (P[3] + z * (P[4] + z * (P[5] + z * (P[6] + z * P[7]))))));
  const double q = Q[0] + z * (Q[1] + z * (Q[2] + z * (Q[3] + z * (Q[4] + z * (Q[5] + z * (Q[6] + z * Q[7]))))));

  result->val = z * (-M_EULER + z*(LNGAMMA_M_PI2_12 + z*(p/q)));
  result->err = 16.0 * GSL_DBL_EPSILON * fabs(result->val);

  return GSL_SUCCESS;
}

/*
 * Approximate log(Gamma(x)) near x = 2.
 *
 * With z = x - 2, write
 *
 *   log(Gamma(x)) = z [c1 + z {c2 + z P(z)/Q(z)}],
 *
 * where
 *
 *   c1 = 1 - gamma,
 *   c2 = pi^2/12 - 1/2.
 *
 * The outer form imposes the exact zero at x = 2 and the first two
 * Taylor coefficients there; P/Q is a [6/6] rational approximation
 * to the remaining smooth part.
 *
 * This approximation is selected on
 *
 * Domain: LNGAMMA_X_MIDDLE_NEAR2 < x < LNGAMMA_X_NEAR2_RECUR
 */
static int
lngamma_near2(const double x, gsl_sf_result * result)
{
  static const double P[7] = { -0.06735230105319809513324605383714999692166209744683296059742385178061274,
                               -0.09297580366286692610049446044951716899433649598803395565332213918566099,
                               -0.04761589529161025463580567694575034560801296559582730116744556362317244,
                               -0.01105523416343014105197270673755015054690756224904527125037500902343221,
                               -0.001123507034936733143845412470127212431509922791261865953934448841884832,
                               -0.00003777034022023585056876623709205318049903406127231183649632543847810373,
                               -0.00000001755173672055158448181836390840519332592419385385258079740124503033794 };
  static const double Q[7] = { 1.0,
                               1.686009391141053755446712212557175220918127814540515772853914804979866,
                               1.112505131640161913812979704307211899093343756105896551739047540500805,
                               0.36212409537463234264005349962663296263488597766431691979378094943082,
                               0.05999088134290343932014480085264525137135826941805478511415374721240256,
                               0.004637228288666860787633503749797005200804329168986553598899570366110877,
                               0.0001235470951556842591565531689170904136774810451737489336054392223257045 };
  const double z = x - 2.0;
  const double p = P[0] + z * (P[1] + z * (P[2] + z * (P[3] + z * (P[4] + z * (P[5] + z * P[6])))));
  const double q = Q[0] + z * (Q[1] + z * (Q[2] + z * (Q[3] + z * (Q[4] + z * (Q[5] + z * Q[6])))));

  result->val = z * (LNGAMMA_M_1_MINUS_EULER + z*(LNGAMMA_M_PI2_12_MINUS_HALF + z*(p/q)));
  result->err = 16.0 * GSL_DBL_EPSILON * fabs(result->val);

  return GSL_SUCCESS;
}

/*
 * Approximate log(Gamma(x)) on the interval between the expansions
 * centered at x = 1 and x = 2.
 *
 * With t = x - 3/2, write
 *
 *   log(Gamma(x)) = log(Gamma(3/2)) + P(t)/Q(t),
 *
 * where P/Q is a [7/7] rational approximation and P(0) = 0.
 *
 * Domain: LNGAMMA_X_NEAR1_MIDDLE < x <= LNGAMMA_X_MIDDLE_NEAR2
 */
static int
lngamma_middle(const double x, gsl_sf_result * result)
{
  static const double P[8] = { 0.0,
                               0.03648997397857652055902366700124443280684039533956589295287274612834503,
                               0.5447031128143503468569140973921498697590166676909875805778942884868158,
                               0.9150711948302330765088136899219005895539583170088797951807208713875458,
                               0.5983086299313562727894681773071126545645078229675719525770881895635743,
                               0.1784012663527851630066469003794348385367236316776177121335390577620073,
                               0.02305557153131806432348794290686215595157738598253799278167507468952008,
                               0.0009526937144704316742101012945081432361249392587505456153460928607597461 };
  static const double Q[8] = { 1.0,
                               2.118445263550891011838143552137623343233496178468925350429694199903467,
                               1.727606653155870656864627058011710093751952197176559662985728958340636,
                               0.6779389294829769558442237205489904278651350057703522195813128128001802,
                               0.1300253054389699879329483215640318638492551879188947883108606838514124,
                               0.01078767618619780920246915047253657293246158614509083491559310706339744,
                               0.000252861047640508135666699794760537494685483085481773802684733723139857,
                              -0.0000009363176663120461007519319682075362461889989363777612601791741072761643 };
  const double t = x - 1.5;
  const double p = P[0] + t * (P[1] + t * (P[2] + t * (P[3] + t * (P[4] + t * (P[5] + t * (P[6] + t * P[7]))))));
  const double q = Q[0] + t * (Q[1] + t * (Q[2] + t * (Q[3] + t * (Q[4] + t * (Q[5] + t * (Q[6] + t * Q[7]))))));
  const double ratio = p / q;

  result->val = LNGAMMA_M_LNGAMMA_3_2 + ratio;
  result->err = 8.0 * GSL_DBL_EPSILON * (LNGAMMA_M_ABS_LNGAMMA_3_2 + fabs(ratio));

  return GSL_SUCCESS;
}

/*
 * Reduce x by an integer so that y = x - n lies in a target interval.
 *
 * The target interval is [lo,hi] when hi_closed is nonzero and
 * [lo,hi) otherwise.  The initial shift
 *
 *   n = floor(x - lo)
 *
 * places y relative to the lower endpoint.  If y lies beyond the
 * permitted upper endpoint, one additional recurrence step is taken.
 *
 * Return the integer shift n and store the reduced argument in *y.
 */
static int
lngamma_reduce(const double x, const double lo,
               const double hi, const int hi_closed,
               double * y)
{
  int n = (int) floor(x - lo);

  *y = x - (double) n;

  if ((hi_closed && *y > hi) ||
      (!hi_closed && *y >= hi))
    {
      ++n;
      *y -= 1.0;
    }

  return n;
}

/*
 * Evaluate log(Gamma(x)) by recurrence on the interval between the
 * direct rational approximations and the Stirling approximation.
 *
 * For an integer n and y = x - n,
 *
 *   log(Gamma(x)) = log(Gamma(y)) + log(product_{j=0}^{n-1} (y + j)).
 *
 * The reduced argument is chosen in either the middle or near x=2
 * approximation interval.  Integer shifts preserve the fractional
 * part of x, so the target interval can be selected directly from
 * frac(x).  The near x=2 interval is preferred whenever it is
 * reachable; otherwise the middle interval is used.
 *
 * The recurrence factors are accumulated as a product before taking
 * a single logarithm.  This form was found to give lower binary64
 * error than summing the individual logarithms.
 *
 * Domain: LNGAMMA_X_NEAR2_RECUR <= x < LNGAMMA_X_RECUR_STIRLING
 */
static int
lngamma_recurrence(const double x, gsl_sf_result * result)
{
  int status = GSL_SUCCESS;
  const double f = x - floor(x);
  const double f0 = LNGAMMA_X_NEAR2_RECUR - 2.0;
  const double f1 = LNGAMMA_X_MIDDLE_NEAR2 - 1.0;
  gsl_sf_result r;
  double y;
  double p = 1.0;
  double lnp;
  int n;
  int j;

  if (f >= f0 && f < f1)
    {
      n = lngamma_reduce(x, LNGAMMA_X_NEAR1_MIDDLE, LNGAMMA_X_MIDDLE_NEAR2, 1, &y);
      status = lngamma_middle(y, &r);
    }
  else
    {
      n = lngamma_reduce(x, LNGAMMA_X_MIDDLE_NEAR2, LNGAMMA_X_NEAR2_RECUR, 0, &y);
      status = lngamma_near2(y, &r);
    }

  if (status)
    {
      result->val = GSL_NAN;
      result->err = GSL_NAN;
      return status;
    }

  for (j = 0; j < n; ++j)
    p *= y + (double)j;

  lnp = log(p);
  result->val = r.val + lnp;
  result->err = r.err + GSL_DBL_EPSILON * ((double)n + 2.0 * (fabs(lnp) + fabs(result->val)));

  return GSL_SUCCESS;
}

/*
 * Evaluate log(Gamma(x)) for large positive x using Stirling's
 * expansion with an approximated remainder.
 *
 *   log(Gamma(x)) = (x - 1/2)(log(x) - 1) - 1/2 + log(sqrt(2*pi)) + R(x).
 *
 * With
 *
 *   t = 1/x,
 *   u = t^2,
 *   v = LNGAMMA_STIRLING_SCALE*u - 1,
 *
 * approximate the remainder by
 *
 *   R(x) = t P(v),
 *
 * where P is a degree-4 polynomial.  The transformation maps
 *
 *   x in [LNGAMMA_X_RECUR_STIRLING, infinity)
 *
 * onto v in [1,-1], so the approximation covers the complete
 * semi-infinite large-x domain.
 */
static int
lngamma_stirling(const double x, gsl_sf_result * result)
{
  static const double P[5] = { 0.0833116801042943673193300712526394485175055312775091535458318332065380896085724915926523506,
                              -0.0000216053468869444527756868895079526119339815771522477656302614326015892438456960547648136461,
                               0.0000000476073820463398096823885412011280987886273825150989355526167045268157021593489042443611532,
                              -0.000000000271891336223976711247020197256663724316025619734256434929266139301491651564665700065352107,
                               0.00000000000287863899744118155562453381913130753526379718867925114822298836372133477234622468542424781 };
  const double t = 1.0/x;
  const double u = t*t;
  const double v = LNGAMMA_STIRLING_SCALE * u - 1.0;
  const double l = log(x);
  const double p = P[0] + v * (P[1] + v * (P[2] + v * (P[3] + v * P[4])));
  const double a = (x - 0.5) * (l - 1.0);
  const double r = t * p;

  result->val = a - 0.5 + LNGAMMA_M_LN_SQRT_2PI + r;
  result->err = 2.0 * GSL_DBL_EPSILON * (fabs(a) +
                                         0.5 +
                                         LNGAMMA_M_LN_SQRT_2PI +
                                         fabs(r) +
                                         fabs(result->val));

  return GSL_SUCCESS;
}

/***************************************************
 * INTERNAL ROUTINES FOR NEGATIVE-ARGUMENT LNGAMMA *
 * (based on glibc lgamma)                         *
 ***************************************************/

/* The complete glibc negative-axis zero table. */
static const double lgamma_zeros[53][2] =
{
  { -2.457024738220800585963843e+00, -3.707561081551326566914472e-17 },
  { -2.747682646727412691944892e+00,  9.055340329338315233540028e-17 },
  { -3.143580888349979840512560e+00, -2.181817985233171446891983e-16 },
  { -3.955294284858597908538513e+00, -1.999428391746348038798085e-17 },
  { -4.039361839740537085674532e+00,  2.114399550398060184003189e-16 },
  { -4.991544640560047874089378e+00,  1.517441176057172217708900e-16 },
  { -5.008218168322593477626015e+00, -4.392635349101581520792310e-17 },
  { -5.998607480080875298256160e+00, -3.311862478893795229223795e-16 },
  { -6.001385294453155161420455e+00,  6.415847287933041753288166e-17 },
  { -6.999801507890637708442227e+00,  1.055013003740002320402941e-17 },
  { -7.000198333407325002042398e+00,  2.504354173632409234702452e-16 },
  { -7.999975197095820611536965e+00, -5.261737128572353905573783e-17 },
  { -8.000024800270681524239080e+00, -4.354586297860107115062812e-16 },
  { -8.999997244250977246338152e+00, -2.218562050972713201432003e-16 },
  { -9.000002755714822555432875e+00, -9.491348611623207875239968e-17 },
  { -9.999999724426629654772114e+00,  4.883037618642442560127648e-16 },
  { -1.000000027557301329750317e+01, -3.490970833264205701540837e-16 },
  { -1.099999997494789027996376e+01,  1.984399830698540662000643e-16 },
  { -1.100000002505210616732256e+01, -6.850849812286175029300759e-16 },
  { -1.199999999791232418999698e+01, -1.002069392010303614934054e-16 },
  { -1.200000000208767581000302e+01,  1.222254811204818497275399e-16 },
  { -1.299999999983941023629086e+01,  6.747262033096336975224083e-16 },
  { -1.300000000016058976370914e+01, -6.745919484964341975148073e-16 },
  { -1.399999999998853006388799e+01,  8.094860741926607006434527e-16 },
  { -1.400000000001146993611201e+01, -8.094853704222662091495788e-16 },
  { -1.499999999999923616655906e+01,  8.829322414768680367962527e-16 },
  { -1.500000000000076383344094e+01, -8.829322382710274453364822e-16 },
  { -1.599999999999995203836534e+01, -1.668613399265053980612921e-16 },
  { -1.600000000000004618527782e+01, -1.609495499460936677225091e-15 },
  { -1.699999999999999644728632e+01, -7.412564244549575528864648e-16 },
  { -1.700000000000000355271368e+01,  7.412564244550028137809019e-16 },
  { -1.800000000000000000000000e+01,  1.561920696858623266671938e-16 },
  { -1.800000000000000000000000e+01, -1.561920696858622034076774e-16 },
  { -1.900000000000000000000000e+01,  8.220635246624329495537040e-18 },
  { -1.900000000000000000000000e+01, -8.220635246624329495537040e-18 },
  { -2.000000000000000000000000e+01,  4.110317623312164844065017e-19 },
  { -2.000000000000000000000000e+01, -4.110317623312164844065017e-19 },
  { -2.100000000000000000000000e+01,  1.957294106339126259519796e-20 },
  { -2.100000000000000000000000e+01, -1.957294106339126259519796e-20 },
  { -2.200000000000000000000000e+01,  8.896791392450574077889159e-22 },
  { -2.200000000000000000000000e+01, -8.896791392450574077889159e-22 },
  { -2.300000000000000000000000e+01,  3.868170170630684126148688e-23 },
  { -2.300000000000000000000000e+01, -3.868170170630684126148688e-23 },
  { -2.400000000000000000000000e+01,  1.611737571096118385895287e-24 },
  { -2.400000000000000000000000e+01, -1.611737571096118385895287e-24 },
  { -2.500000000000000000000000e+01,  6.446950284384473589498896e-26 },
  { -2.500000000000000000000000e+01, -6.446950284384473589498896e-26 },
  { -2.600000000000000000000000e+01,  2.479596263224797589612253e-27 },
  { -2.600000000000000000000000e+01, -2.479596263224797589612253e-27 },
  { -2.700000000000000000000000e+01,  9.183689863795546005394213e-29 },
  { -2.700000000000000000000000e+01, -9.183689863795546005394213e-29 },
  { -2.800000000000000000000000e+01,  3.279889237069837758976614e-30 }
};

static const double e_hi = 2.718281828459045090795598e+00;
static const double e_lo = 1.445646891729250157831396e-16;

static const double lgamma_coeff[] =
{
  8.333333333333332870740406e-02,
 -2.777777777777777883788657e-03,
  7.936507936507936501052685e-04,
 -5.952380952380952917890600e-04,
  8.417508417508417139715760e-04,
 -1.917526917526917633674555e-03,
  6.410256410256410034009811e-03,
 -2.955065359477124231624146e-02,
  1.796443723688305738050985e-01,
 -1.392432216905901132264489e+00,
  1.340286404416839260989036e+01,
 -1.568482846260020266981883e+02
};

#define NCOEFF (sizeof(lgamma_coeff) / sizeof(lgamma_coeff[0]))

/*
 * glibc special polynomial branch for -3 < x < -2.
 *
 * The polynomials approximate
 *
 *   (|Gamma(x)| - 1) * (x-n) / (x-x0)
 *
 * on eight subintervals.  The coefficients below are decimal
 * representations of glibc's binary64 hexadecimal constants.
 */
static const double lgamma_neg_poly_coeff[] =
{
  /* interval 0 */
  -1.044704781216988953573832e+00,
  -7.782305330904056850371830e-01,
  -1.202431471366234433384435e-01,
  -8.886202355427311738722551e-01,
  -6.342843101814840678365215e-02,
  -9.176282204872978631726710e-01,
  3.785861085708514678360714e-02,
  -9.341733618747365408552241e-01,
  1.520007117912068050902263e-01,
  -9.697252578163961800328252e-01,
  2.788205613765887269650534e-01,
  /* interval 1 */
  -9.475560525532616962252064e-01,
  -7.904814408001735737130389e-01,
  2.257366771626493906488520e-01,
  -1.010340524227663516754205e+00,
  5.895033531088182954249532e-01,
  -1.302783186460660491690078e+00,
  1.109566968329737335352547e+00,
  -1.787705247185346291161068e+00,
  1.843167252209331419265936e+00,
  -2.570250797415364552023220e+00,
  2.962724653240909944429404e+00,
  -3.881766977228747972361589e+00,
  /* interval 2 */
  -8.430565186694192192362607e-01,
  -9.007307074535096935719025e-01,
  6.907850597415765570019630e-01,
  -1.571457453663721803849285e+00,
  1.830491229184272983232518e+00,
  -3.033580501343014823589783e+00,
  4.126868004210012763621762e+00,
  -6.226302470384837306482950e+00,
  8.886294054752328364088498e+00,
  -1.305300533156573727922023e+01,
  1.888301757578075878996060e+01,
  -2.830747086445322580061656e+01,
  4.141446030082260421067986e+01,
  /* interval 3 */
  -7.160435758996516186769554e-01,
  -1.166057428589025501963761e+00,
  1.531281442042658369473429e+00,
  -3.195917070855539332541184e+00,
  5.342476029602517328953581e+00,
  -9.760480704063542845005941e+00,
  1.717216724690779017237219e+01,
  -3.066134609934890775662097e+01,
  5.441945910419376275513059e+01,
  -9.681174791840648197194241e+01,
  1.719325261861630167459225e+02,
  -3.055345717488169157149969e+02,
  5.679666269469928465696285e+02,
  -1.021940037539013587775116e+03,
  /* interval 4 */
  -2.385263769319141413660645e-01,
  1.802085679017021879744220e+00,
  3.506275674054836866133655e+00,
  6.519619683964499756712030e+00,
  1.175663604365609060664610e+01,
  2.104033874311103247123356e+01,
  3.747998021215248343196436e+01,
  6.669270764646357463334425e+01,
  1.185969214540074005981296e+02,
  2.108703290278693316395220e+02,
  3.746001914032927402331552e+02,
  6.649256300730046405078610e+02,
  1.237141804029709646783886e+03,
  2.281385733818075550516369e+03,
  /* interval 5 */
  -4.193900326125993283099547e-01,
  1.159515590372600879476295e+00,
  1.857011060616613473683856e+00,
  2.883164498875658399157373e+00,
  4.269547627846991488809181e+00,
  6.292870578714174989443109e+00,
  9.178552045830782546431692e+00,
  1.338574542288015223334696e+01,
  1.947644230073518301082913e+01,
  2.833631449884830999508267e+01,
  4.120381229805598621851459e+01,
  6.161278288501627997675314e+01,
  9.072645743800207185358886e+01,
  /* interval 6 */
  -5.400649246179795026634451e-01,
  8.033521909763831070350193e-01,
  1.079905352604727974963339e+00,
  1.461407397124978579938670e+00,
  1.822170255131376759294426e+00,
  2.306106582214320965817933e+00,
  2.833821733808250531438944e+00,
  3.519854250618171587206007e+00,
  4.320615359067300431661351e+00,
  5.335365175493834932751724e+00,
  6.675528681738811087598151e+00,
  8.282177097586464142864315e+00,
  /* interval 7 */
  -6.260808948176509458605210e-01,
  5.900197543087655205340525e-01,
  6.661716257481942182749890e-01,
  8.240481145327041012649261e-01,
  8.619638740995989589777082e-01,
  9.880464187192510294721615e-01,
  1.018329045410495403700679e+00,
  1.135579046345332798750860e+00,
  1.173137803102460852144873e+00,
  1.291562533161418890870209e+00,
  1.360338660728580562064849e+00,
  1.492312024525090263793459e+00
};

static const unsigned int lgamma_neg_poly_deg[] =
{
  10, 11, 12, 13, 13, 12, 11, 11
};

static const unsigned int lgamma_neg_poly_end[] =
{
  10, 22, 35, 49, 63, 76, 88, 100
};

/*
 * Evaluate log|Gamma(x)| for -3 < x < -2 using the special
 * negative-argument polynomial approximation from glibc.
 *
 * This interval contains two zeros x0 of log|Gamma(x)|.  For each
 * half-integer region choose the corresponding zero
 *
 *   x0 = x0_hi + x0_lo
 *
 * and write
 *
 *   |Gamma(x)| - 1 = g(x) (x - x0) / (x - xn),
 *
 * where xn is the adjacent negative integer pole.  The smooth factor
 * g(x) is approximated by one of eight local polynomials on subintervals
 * of width 1/8.  The logarithm is then evaluated as
 *
 *   log|Gamma(x)| = log1p(g(x) (x - x0) / (x - xn)),
 *
 * which preserves accuracy near the zeros of log|Gamma(x)|.
 *
 * The zero locations are stored as high/low pairs so that x-x0 can be
 * formed accurately.  Negative integer arguments are poles of Gamma
 * and produce a domain error.
 *
 * This algorithm is based on the negative-argument lgamma algorithm
 * from glibc.
 *
 * Domain: -3 < x < -2, with x=-3 handled as a pole at the dispatcher
 * boundary.
 */
static int
lngamma_neg_minus3_minus2(const double x, gsl_sf_result * result)
{
  int i, j, k;
  double xn;
  double x0_hi;
  double x0_lo;
  double xdiff;
  double xm;
  double x_adj;
  double g;
  double den, q;

  i = (int) floor(-2.0 * x);
  if ((i & 1) == 0 && (double)i == -2.0 * x)
    {
      DOMAIN_ERROR(result);
    }

  xn = ((i & 1) == 0) ? -(double)i / 2.0 : (-(double)i - 1.0) / 2.0;

  /*
   * For -3 < x < -2, i is 4,5,... and subtracting 4 maps the
   * eight half-integer subintervals onto polynomial indices 0..7.
   */
  k = i - 4;

  x0_hi = lgamma_zeros[k][0];
  x0_lo = lgamma_zeros[k][1];
  xdiff = x - x0_hi - x0_lo;

  /*
   * This is exactly glibc's subinterval selection:
   *
   *   j = floor(-8*x) - 16
   *   xm = (-33 - 2*j) / 16.
   */
  j = (int) floor(-8.0 * x) - 16;

  if (j < 0) j = 0;
  if (j > 7) j = 7;

  xm = (-33.0 - 2.0 * (double)j) * 0.0625;
  x_adj = x - xm;

  g = lgamma_neg_poly_coeff[lgamma_neg_poly_end[j]];

  for (k = 1; k <= (int)lgamma_neg_poly_deg[j]; ++k)
    g = g * x_adj + lgamma_neg_poly_coeff[lgamma_neg_poly_end[j] - k];

  den = x - xn;
  q = g * xdiff / den;

  result->val = log1p(q);
  result->err = 16.0 * GSL_DBL_EPSILON * (fabs(result->val) + fabs(q / (1.0 + q)) + 1.0);

  return GSL_SUCCESS;
}

/*
 * Evaluate log|Gamma(x)| for -2 < x < -1 by ordinary reflection:
 *
 *   log|Gamma(x)| = log(pi / |x sin(pi*x)|) - log(Gamma(-x)).
 *
 * The reflected argument -x lies in (1,2), so it is evaluated directly
 * with the appropriate near-1, middle, or near-2 positive-argument
 * kernel.  This keeps the local fallback independent of the system
 * lgamma() implementation.
 */
static int
lngamma_neg_minus2_minus1(const double x, gsl_sf_result * result)
{
  const double ax = -x;
  gsl_sf_result rs, rg;
  double nadj, err_nadj;
  int status;

  status = gsl_sf_sin_pi_e(x, &rs);
  if (status)
    {
      result->val = GSL_NAN;
      result->err = GSL_NAN;
      return status;
    }

  if (rs.val == 0.0)
    {
      DOMAIN_ERROR(result);
    }

  if (ax <= LNGAMMA_X_NEAR1_MIDDLE)
    status = lngamma_near1(ax, &rg);
  else if (ax <= LNGAMMA_X_MIDDLE_NEAR2)
    status = lngamma_middle(ax, &rg);
  else
    status = lngamma_near2(ax, &rg);

  if (status)
    {
      *result = rg;
      return status;
    }

  nadj = log(M_PI / fabs(x * rs.val));
  err_nadj = rs.err / fabs(rs.val) + 4.0 * GSL_DBL_EPSILON * (fabs(nadj) + 1.0);

  result->val = nadj - rg.val;
  result->err = err_nadj + rg.err + GSL_DBL_EPSILON * (fabs(nadj) + fabs(rg.val) + fabs(result->val));

  return GSL_SUCCESS;
}

/*
 * Compute the product x*y as the unevaluated sum
 *
 *   x*y = hi + lo,
 *
 * where hi is the rounded binary64 product and lo contains its rounding
 * residual.  The operands are split into high and low parts using the
 * Dekker constant 2^27+1.
 */
static void
lngamma_neg_mul_split(double *hi, double *lo, const double x, const double y)
{
  const double split = 134217729.0;  /* 2^27 + 1 */
  double tx, ty;
  double xh, xl;
  double yh, yl;

  *hi = x * y;

  tx = split * x;
  xh = tx - (tx - x);
  xl = x - xh;

  ty = split * y;
  yh = ty - (ty - y);
  yl = y - yh;

  *lo = ((xh*yh - *hi) + xh*yl + xl*yh) + xl*yl;
}

/*
 * Compute, with compensated arithmetic,
 *
 *   product_{i=0}^{n-1} (1 + t/(x+i)) - 1.
 *
 * The correction x_eps represents the low-order part of x.  Product
 * and quotient rounding errors are accumulated explicitly using
 * lngamma_neg_mul_split().  Returning the product minus one directly
 * avoids loss of relative accuracy when the product is close to unity.
 *
 * This quantity is used by lngamma_neg() when recurrence is required
 * before evaluating the difference of two Stirling expansions.
 */
static double
lgamma_product(const double t, const double x, const double x_eps,
               const int n)
{
  double ret = 0.0;
  double ret_eps = 0.0;
  int i;

  for (i = 0; i < n; ++i)
    {
      const double xi = x + (double)i;
      const double quot = t / xi;
      double mhi, mlo;
      double quot_lo;
      double rhi, rlo;
      double rpq, rpq_eps, nret, nret_eps;

      lngamma_neg_mul_split(&mhi, &mlo, quot, xi);
      quot_lo = (t - mhi - mlo) / xi
        - t * x_eps / (xi * xi);

      lngamma_neg_mul_split(&rhi, &rlo, ret, quot);
      rpq = ret + quot;
      rpq_eps = (ret - rpq) + quot;
      nret = rpq + rhi;
      nret_eps = (rpq - nret) + rhi;

      ret_eps += rpq_eps + nret_eps + rlo
        + ret_eps * quot + quot_lo
        + quot_lo * (ret + ret_eps);

      ret = nret;
    }

  return ret + ret_eps;
}

/*
 * Compute sin(pi*x) for 0 <= x <= 1/2.
 *
 * For x > 1/4, use
 *
 *   sin(pi*x) = cos(pi*(1/2-x))
 *
 * so that the trigonometric argument remains in [0,pi/4].
 */
static double
lngamma_neg_sinpi(const double x)
{
  if (x <= 0.25)
    return sin(M_PI * x);
  else
    return cos(M_PI * (0.5 - x));
}

/*
 * Compute cos(pi*x) for 0 <= x <= 1/2.
 *
 * For x > 1/4, use
 *
 *   cos(pi*x) = sin(pi*(1/2-x))
 *
 * so that the trigonometric argument remains in [0,pi/4].
 */
static double
lngamma_neg_cospi(const double x)
{
  if (x <= 0.25)
    return cos(M_PI * x);
  else
    return sin(M_PI * (0.5 - x));
}


static double
lngamma_neg_cotpi(const double x)
{
  return lngamma_neg_cospi(x) / lngamma_neg_sinpi(x);
}

/*
 * Evaluate log|Gamma(x)| on the negative interval
 *
 *   LNGAMMA_X_NEG0_MIN < x < LNGAMMA_X_NEG0_MAX.
 *
 * Direct use of the reflection formula loses accuracy near the zeros of
 * log|Gamma(x)| because it subtracts terms of nearly equal magnitude.
 * Instead, choose a nearby tabulated zero
 *
 *   x0 = x0_hi + x0_lo
 *
 * and evaluate the reflection formula relative to x0.  Since
 * log|Gamma(x0)| = 0, write
 *
 *   log|Gamma(x)| = log|sin(pi*x0)/sin(pi*x)| + log(Gamma(1-x0)/Gamma(1-x)).
 *
 * The sine ratio is evaluated either directly or in a log1p form when
 * it is close to unity.  The Gamma ratio is evaluated as the direct
 * difference of two Stirling expansions rather than by subtracting two
 * independently computed log-Gamma values.  The zero locations are
 * stored as high/low pairs so that x-x0 can be formed accurately.
 *
 * Negative integer arguments are poles of Gamma and produce a domain
 * error.
 *
 * This algorithm is based on the negative-argument lgamma algorithm
 * from glibc.
 *
 * Domain: LNGAMMA_X_NEG0_MIN < x < LNGAMMA_X_NEG0_MAX.
 */
static int
lngamma_neg(const double x, gsl_sf_result * result)
{
  int i;
  double xn;
  double x0_hi, x0_lo;
  double xdiff;
  double x_idiff, x0_idiff;
  double log_sinpi_ratio;
  double x0diff2, sx0d2, cx0d2;
  double log_gamma_ratio;
  double y0, y0_eps;
  double y, y_eps;
  double log_gamma_adj = 0.0;
  double log_gamma_high;
  double y0r, yr, y0r2, yr2;
  double rdiff;
  double bterm[NCOEFF];
  double dlast, elast;
  double log_gamma_low = 0.0;
  int j;

  /*
   * Determine the half-integer region exactly as in glibc.
   */
  i = (int)floor(-2.0 * x);

  if ((i & 1) == 0 && (double)i == -2.0 * x)
    {
      /* the input point is a Gamma pole */
      DOMAIN_ERROR(result);
    }

  xn = ((i & 1) == 0)
       ? -(double)i / 2.0
       : (-(double)i - 1.0) / 2.0;

  i -= 4;

  x0_hi = lgamma_zeros[i][0];
  x0_lo = lgamma_zeros[i][1];

  xdiff = x - x0_hi - x0_lo;

  x_idiff = fabs(xn - x);
  x0_idiff = fabs(xn - x0_hi - x0_lo);

  if (x0_idiff < x_idiff * 0.5)
    {
      log_sinpi_ratio =
        log(lngamma_neg_sinpi(x0_idiff) / lngamma_neg_sinpi(x_idiff));
    }
  else
    {
      x0diff2 =
        ((i & 1) == 0 ? xdiff : -xdiff) * 0.5;
      sx0d2 = lngamma_neg_sinpi(x0diff2);
      cx0d2 = lngamma_neg_cospi(x0diff2);

      log_sinpi_ratio =
        log1p(2.0 * sx0d2 *
              (-sx0d2 + cx0d2 * lngamma_neg_cotpi(x_idiff)));
    }

  /*
   * Construct y=1-x and y0=1-x0 together with the low-order
   * rounding residuals, following glibc.
   */
  y0 = 1.0 - x0_hi;
  y0_eps = -x0_hi + (1.0 - y0) - x0_lo;

  y = 1.0 - x;
  y_eps = -x + (1.0 - y);

  /*
   * Move y and y0 upward when required so that the Stirling remainder
   * is sufficiently small.
   */
  if (i < 6)
    {
      const int n_up = (7 - i) / 2;
      double ny0, ny0_eps, ny, ny_eps;
      double prodm1;

      ny0 = y0 + (double)n_up;
      ny0_eps = y0 - (ny0 - (double)n_up) + y0_eps;
      y0 = ny0;
      y0_eps = ny0_eps;

      ny = y + (double)n_up;
      ny_eps = y - (ny - (double)n_up) + y_eps;
      y = ny;
      y_eps = ny_eps;

      prodm1 = lgamma_product(xdiff, y - (double)n_up, y_eps, n_up);

      log_gamma_adj = -log1p(prodm1);
    }

  /*
   * Direct difference of the leading Stirling terms.
   */
  log_gamma_high =
    xdiff * log1p((y0 - e_hi - e_lo + y0_eps) / e_hi)
    + (y - 0.5 + y_eps) * log1p(xdiff / y)
    + log_gamma_adj;

  /*
   * Direct difference of the Stirling correction terms.
   */
  y0r = 1.0 / y0;
  yr = 1.0 / y;
  y0r2 = y0r * y0r;
  yr2 = yr * yr;

  rdiff = -xdiff / (y * y0);
  dlast = rdiff;
  elast = rdiff * yr * (yr + y0r);

  bterm[0] = dlast * lgamma_coeff[0];

  for (j = 1; j < (int)NCOEFF; ++j)
    {
      const double dnext = dlast * y0r2 + elast;
      const double enext = elast * yr2;

      bterm[j] = dnext * lgamma_coeff[j];
      dlast = dnext;
      elast = enext;
    }

  for (j = 0; j < (int)NCOEFF; ++j)
    log_gamma_low += bterm[NCOEFF - 1 - j];

  log_gamma_ratio = log_gamma_high + log_gamma_low;

  result->val = log_sinpi_ratio + log_gamma_ratio;
  result->err = 8.0 * GSL_DBL_EPSILON * (fabs(log_sinpi_ratio) +
                                         fabs(log_gamma_high) +
                                         fabs(log_gamma_low) +
                                         fabs(result->val));

  return GSL_SUCCESS;
}

/*
 * Evaluate log|Gamma(x)| for x <= LNGAMMA_X_NEG0_MIN by ordinary
 * reflection.
 *
 * For the positive Gamma factor, use the compact large-x lgamma
 * approximation from glibc. This helper is restricted to x <= LNGAMMA_X_NEG0_MIN,
 * so the positive argument -x is always at least |LNGAMMA_X_NEG0_MIN|.
 *
 * The error estimate includes the error in sin(pi*x), the rounding
 * in the reflection term, the glibc asymptotic approximation error,
 * the rounding in its evaluation, and the final subtraction.
 */
static int
lngamma_neg_large(const double x, gsl_sf_result * result)
{
  static const double W[7] = { 4.18938533204672725052e-01,
                               8.33333333333329678849e-02,
                              -2.77777777728775536470e-03,
                               7.93650558643019558500e-04,
                              -5.95187557450339963135e-04,
                               8.36339918996282139126e-04,
                              -1.63092934096575273989e-03 };
  const double ax = -x;
  const double z = 1.0 / ax;
  const double z2 = z * z;
  const double w = W[0] + z * (W[1] + z2 * (W[2] + z2 * (W[3] + z2 * (W[4] + z2 * (W[5] + z2 * W[6])))));
  const double l = log(ax);
  const double a = (ax - 0.5) * (l - 1.0);
  const double lg = a + w;
  const double glibc_bound = 2.077294420884752159959086386551806343744833799956e-18;

  gsl_sf_result rs;
  double nadj;
  double err_nadj;
  double err_lg;
  int status;

  status = gsl_sf_sin_pi_e(x, &rs);
  if (status)
    {
      result->val = GSL_NAN;
      result->err = GSL_NAN;
      return status;
    }

  if (rs.val == 0.0)
    {
      DOMAIN_ERROR(result);
    }

  /*
   * Reflection correction:
   *
   *   nadj = log(pi / |x sin(pi*x)|).
   *
   * The sin(pi*x) error propagates approximately as
   *
   *   |d log(s)| <= |ds| / |s|.
   *
   * The remaining terms account conservatively for the arithmetic
   * in x*s, pi/(x*s), log(), and the final representation of nadj.
   */
  nadj = log(M_PI / fabs(x * rs.val));

  err_nadj = rs.err / fabs(rs.val) + 4.0 * GSL_DBL_EPSILON * (fabs(nadj) + 1.0);

  /*
   * glibc documents an approximation error < 2^-58.74 for W(z).
   * Include that approximation error plus conservative binary64
   * rounding for the evaluation of a, w, and lg.
   */
  err_lg = glibc_bound + 4.0 * GSL_DBL_EPSILON * (fabs(a) + fabs(w) + fabs(lg) + 1.0);

  result->val = nadj - lg;
  result->err = err_nadj + err_lg + GSL_DBL_EPSILON * (fabs(nadj) + fabs(lg) + fabs(result->val));

  return GSL_SUCCESS;
}

/********************************************
 * INTERNAL ROUTINES FOR LNGAMMA_1P         *
 ********************************************/

/*
 * Approximate log(Gamma(1+x)) near the zero x = 0.
 *
 * Write
 *
 *   log(Gamma(1+x)) = x [-gamma + x {pi^2/12 + x P(x)/Q(x)}],
 *
 * where P/Q is a [7/7] rational approximation to the remaining
 * smooth part.  The outer factor x preserves the exact zero at
 * x = 0, while the first two Taylor coefficients are imposed
 * explicitly.
 *
 * The approximation is fitted on
 *
 *   -0.4 <= x <= 0.4.
 *
 * Domain: -0.4 <= x < 0.4.
 *
 * This kernel is also used after recurrence for -1 < x < -0.6,
 * where t = 1+x lies in (0,0.4).
 */
static int
lngamma_1p_near0(const double x, gsl_sf_result * result)
{
  static const double P[8] = { -0.40068563438653142842279739082137715531971723438413634748,
                               -1.3317354749778089010529758510057213259262188165396800741,
                               -1.7317839983677253569005843378093380590143862347246751998,
                               -1.1101413603861377226418795688803457412950901847974067339,
                               -0.36154512422449969262768262262466868830005593565402849976,
                               -0.05471576009973619416822545297762117065584763108720527791,
                               -0.0028346322753687156297184743958066785398027885674436732381,
                               -0.0000014216517515929311015630385424515660620839354763260593875 };
  static const double Q[8] = { 1.0,
                               3.9989361881136944668252386856672006982762229553732551787,
                               6.504934567870272697046333712009729305983204275280299118,
                               5.5167623117100826066825847743968411096567313513884888807,
                               2.5936645483667636326525358306667888299545233003675877476,
                               0.66095827270487214429058431150657937698736844341798616954,
                               0.081687702242817240471539558703447530242990334867981767329,
                               0.0035683051537767774700408217297161398926462406660212714498 };
  const double p = P[0] + x*(P[1] + x*(P[2] + x*(P[3] + x*(P[4] + x*(P[5] + x*(P[6] + x*P[7]))))));
  const double q = Q[0] + x*(Q[1] + x*(Q[2] + x*(Q[3] + x*(Q[4] + x*(Q[5] + x*(Q[6] + x*Q[7]))))));

  result->val = x*(-M_EULER + x*(LNGAMMA_M_PI2_12 + x*(p/q)));
  result->err = 16.0 * GSL_DBL_EPSILON * fabs(result->val);

  return GSL_SUCCESS;
}

/*
 * Approximate log(Gamma(2+y)) near the zero y = 0.
 *
 * Write
 *
 *   log(Gamma(2+y)) = y [(1-gamma) + y {pi^2/12 - 1/2 + y P(y)/Q(y)}],
 *
 * where P/Q is a [6/6] rational approximation to the remaining
 * smooth part.  The outer factor y preserves the exact zero at
 * y = 0, while the first two Taylor coefficients are imposed
 * explicitly.
 *
 * The approximation is fitted for
 *
 *   -0.6 <= y <= 1.0,
 *
 * corresponding to 0.4 <= x <= 2.0 when y = x-1 in
 * log(Gamma(1+x)).
 *
 * Supplying the local coordinate y directly is useful after
 * recurrence for -0.6 <= x < -0.4.  In
 *
 *   log(Gamma(1+x)) = log(Gamma(2+x)) - log(1+x),
 *
 * the local coordinate of log(Gamma(2+x)) is exactly y=x.
 * Passing x directly avoids forming (1+x)-1 and thereby preserves
 * the original binary64 argument.
 */
static int
lngamma_1p_near1_y(const double y, gsl_sf_result * result)
{
  static const double P[7] = { -0.06735230105319809413241873726949532774041131955586528944512140896836719651208825,
                               -0.091922234178324071158990694167050886701067813991455433765613531097275700804072305,
                               -0.046406838296033972453259446523471037375493097993223420098713059851278958564078455,
                               -0.010577524223260432937770849828758955087191357448939417556283469105151431574448165,
                               -0.0010489596561228217388845698751476703462058351806453481156588919887387936640807855,
                               -0.00003409145672888353555714164116463001276599660499922299386689323368263239039311506,
                               -0.000000014321845510335209513665088066924677694319415812500231763328552428777877209932555 };
  static const double Q[7] = { 1.0,
                               1.6703667261085600788305526395143512474545939845512971667987937922513515918372555,
                               1.0897739752584232848477544137497983876409517362817620707777561692490725996237289,
                               0.34980074782445768078297546037915877775282844725329309442245121172302505395943531,
                               0.056939686354838956477884714856953320022726806020854962192017852338825034149391612,
                               0.0043030597015131494697236160308987238267138750774403579057195286149394691984354387,
                               0.00011127415620254066622465957274195658915836539842340579892428368657027772464105129 };
  const double p = P[0] + y*(P[1] + y*(P[2] + y*(P[3] + y*(P[4] + y*(P[5] + y*P[6])))));
  const double q = Q[0] + y*(Q[1] + y*(Q[2] + y*(Q[3] + y*(Q[4] + y*(Q[5] + y*Q[6])))));

  result->val = y*(LNGAMMA_M_1_MINUS_EULER + y*(LNGAMMA_M_PI2_12_MINUS_HALF + y*(p/q)));
  result->err = 16.0 * GSL_DBL_EPSILON * fabs(result->val);

  return GSL_SUCCESS;
}

/*
 * Approximate log(Gamma(1+x)) near the zero x = 1.
 *
 * Set
 *
 *   y = x - 1
 *
 * and evaluate the near-y=0 kernel for log(Gamma(2+y)).
 *
 * The approximation preserves the exact zero at x = 1 and is
 * fitted on
 *
 *   0.4 <= x <= 2.0.
 *
 * For recurrence applications where the local coordinate is already
 * known accurately, lngamma_1p_near1_y() should be called directly
 * rather than reconstructing y from a shifted argument.
 */
static int
lngamma_1p_near1(const double x, gsl_sf_result * result)
{
  return lngamma_1p_near1_y(x - 1.0, result);
}

/* x in (-1,LNGAMMA_1P_X_RECUR1_NEAR0) */
static int
lngamma_1p_recurrence(const double x, gsl_sf_result * result)
{
  gsl_sf_result r;
  const double l1px = log1p(x);
  int status;

  if (x < LNGAMMA_1P_X_RECUR0_RECUR1)
    {
      const double t = 1.0 + x;
      status = lngamma_1p_near0(t, &r);
    }
  else
    {
      status = lngamma_1p_near1_y(x, &r);
    }

  if (status)
    {
      result->val = GSL_NAN;
      result->err = GSL_NAN;
      return status;
    }

  result->val = r.val - l1px;
  result->err = r.err + GSL_DBL_EPSILON * (fabs(r.val) +
                                           2.0 * (fabs(l1px) + fabs(result->val)));

  return GSL_SUCCESS;
}
