/* specfunc/gamma/lngamma1p_kernels.c
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
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_sincos_pi.h>

#include "constants.h"
#include "error.h"
#include "lngamma.h"

static int lngamma1p_near0(const double x, gsl_sf_result * result);
static int lngamma1p_near1_y(const double y, gsl_sf_result * result);
static int lngamma1p_near1(const double x, gsl_sf_result * result);
static int lngamma1p_recurrence(const double x, gsl_sf_result * result);

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
lngamma1p_near0(const double x, gsl_sf_result * result)
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

  result->val = x*(-M_EULER + x*(GSL_SF_PI2_12 + x*(p/q)));
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
lngamma1p_near1_y(const double y, gsl_sf_result * result)
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

  result->val = y*(GSL_SF_1_MINUS_EULER + y*(LNGAMMA_M_PI2_12_MINUS_HALF + y*(p/q)));
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
 * known accurately, lngamma1p_near1_y() should be called directly
 * rather than reconstructing y from a shifted argument.
 */
static int
lngamma1p_near1(const double x, gsl_sf_result * result)
{
  return lngamma1p_near1_y(x - 1.0, result);
}

/* x in (-1,LNGAMMA_1P_X_RECUR1_NEAR0) */
static int
lngamma1p_recurrence(const double x, gsl_sf_result * result)
{
  gsl_sf_result r;
  const double l1px = gsl_sf_log1p(x);
  int status;

  if (x < LNGAMMA_1P_X_RECUR0_RECUR1)
    {
      const double t = 1.0 + x;
      status = lngamma1p_near0(t, &r);
    }
  else
    {
      status = lngamma1p_near1_y(x, &r);
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
