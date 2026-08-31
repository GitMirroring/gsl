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

#include "error.h"
#include "eval.h"

static int lngamma_near1(const double x, gsl_sf_result * result);
static int lngamma_near2(const double x, gsl_sf_result * result);
static int lngamma_middle(const double x, gsl_sf_result * result);
static int lngamma_reduce(const double x, const double lo,
                          const double hi, const int hi_closed, double * y);
static int lngamma_recurrence(const double x, gsl_sf_result * result);
static int lngamma_stirling(const double x, gsl_sf_result * result);

static int lngamma_1p_near0(const double x, gsl_sf_result * result);
static int lngamma_1p_near1_y(const double y, gsl_sf_result * result);
static int lngamma_1p_near1(const double x, gsl_sf_result * result);

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
#define LNGAMMA_X_SMALL_RECUR0_RECUR1      (LNGAMMA_1P_X_NEAR0_NEAR1)
#define LNGAMMA_X_SMALL_RECUR1_NEAR1       (0.8)
#define LNGAMMA_X_NEAR1_MIDDLE             (1.23395)
#define LNGAMMA_X_MIDDLE_NEAR2             (1.67895)
#define LNGAMMA_X_NEAR2_RECUR              (2.3)
#define LNGAMMA_X_RECUR_STIRLING           (8.0)

/*
 * Compute log(Gamma(x)).
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
gsl_new_sf_lngamma_e(const double x, gsl_sf_result * result)
{
  int status;
  gsl_sf_result r;

  if (x <= LNGAMMA_X_SMALL_RECUR0_RECUR1)
    {
      status = lngamma_1p_near0(x, &r);
      if (status)
        return status;

      result->val = r.val - log(x);

      return GSL_SUCCESS;
    }
  else if (x <= LNGAMMA_X_SMALL_RECUR1_NEAR1)
    {
      status = lngamma_1p_near1(x, &r);
      if (status)
        return status;

      result->val = r.val - log(x);

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
  int status;
  gsl_sf_result r;

  if (x <= -1.0)
    {
      DOMAIN_ERROR(result);
    }
  else if (x < LNGAMMA_1P_X_RECUR0_RECUR1)
    {
      const double t = 1.0 + x;
      const double l1px = log1p(x);

      status = lngamma_1p_near0(t, &r);
      if (status)
        return status;

      result->val = r.val - l1px;
      result->err = r.err + GSL_DBL_EPSILON * (fabs(r.val) +
                                               2.0 * (fabs(l1px) + fabs(result->val)));

      return GSL_SUCCESS;
    }
  else if (x < LNGAMMA_1P_X_RECUR1_NEAR0)
    {
      const double l1px = log1p(x);

      status = lngamma_1p_near1_y(x, &r);
      if (status)
        return status;

      result->val = r.val - l1px;
      result->err = r.err + GSL_DBL_EPSILON * (fabs(r.val) +
                                               2.0 * (fabs(l1px) + fabs(result->val)));

      return GSL_SUCCESS;
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
      return gsl_new_sf_lngamma_e(1.0 + x, result);
    }
}

double gsl_new_sf_lngamma(const double x) { EVAL_RESULT(gsl_new_sf_lngamma_e(x, &result)); }
double gsl_sf_lngamma_1p(const double x) { EVAL_RESULT(gsl_sf_lngamma_1p_e(x, &result)); }

/********************************************
 * INTERNAL ROUTINES FOR LNGAMMA            *
 ********************************************/

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
    return status;

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
