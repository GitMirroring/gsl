/* specfunc/gamma_inc.c
 *
 * Copyright (C) 2007 Brian Gough
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
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

/* Author:  G. Jungman */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_expint.h>

#include "error.h"

/* Continued fraction which occurs in evaluation
 * of Q(a,x) or Gamma(a,x).
 *
 *              1   (1-a)/x  1/x  (2-a)/x   2/x  (3-a)/x
 *   F(a,x) =  ---- ------- ----- -------- ----- -------- ...
 *             1 +   1 +     1 +   1 +      1 +   1 +
 *
 * Hans E. Plesser, 2002-01-22 (hans dot plesser at itf dot nlh dot no).
 *
 * Split out from gamma_inc_Q_CF() by GJ [Tue Apr  1 13:16:41 MST 2003].
 * See gamma_inc_Q_CF() below.
 *
 */
static int
gamma_inc_F_CF(const double a, const double x, gsl_sf_result * result)
{
  const int    nmax  =  5000;
  const double small =  gsl_pow_3 (GSL_DBL_EPSILON);

  double hn = 1.0;           /* convergent */
  double Cn = 1.0 / small;
  double Dn = 1.0;
  int n;

  /* n == 1 has a_1, b_1, b_0 independent of a,x,
     so that has been done by hand                */
  for ( n = 2 ; n < nmax ; n++ )
  {
    double an;
    double delta;

    if(GSL_IS_ODD(n))
      an = 0.5*(n-1)/x;
    else
      an = (0.5*n-a)/x;

    Dn = 1.0 + an * Dn;
    if ( fabs(Dn) < small )
      Dn = small;
    Cn = 1.0 + an/Cn;
    if ( fabs(Cn) < small )
      Cn = small;
    Dn = 1.0 / Dn;
    delta = Cn * Dn;
    hn *= delta;
    if(fabs(delta-1.0) < GSL_DBL_EPSILON) break;
  }

  result->val = hn;
  result->err = 2.0*GSL_DBL_EPSILON * fabs(hn);
  result->err += GSL_DBL_EPSILON * (2.0 + 0.5*n) * fabs(result->val);

  if(n == nmax)
    GSL_ERROR ("error in CF for F(a,x)", GSL_EMAXITER);
  else
    return GSL_SUCCESS;
}


/* Useful for small a and x. Handles the subtraction analytically.
 */
static
int
gamma_inc_Q_series(const double a, const double x, gsl_sf_result * result)
{
  double term1;  /* 1 - x^a/Gamma(a+1) */
  double sum;    /* 1 + (a+1)/(a+2)(-x)/2! + (a+1)/(a+3)(-x)^2/3! + ... */
  int stat_sum;
  double term2;  /* a temporary variable used at the end */

  {
    /* Evaluate series for 1 - x^a/Gamma(a+1), small a
     */
    const double pg21 = -2.404113806319188570799476;  /* PolyGamma[2,1] */
    const double lnx  = log(x);
    const double el   = M_EULER+lnx;
    const double c1 = -el;
    const double c2 = M_PI*M_PI/12.0 - 0.5*el*el;
    const double c3 = el*(M_PI*M_PI/12.0 - el*el/6.0) + pg21/6.0;
    const double c4 = -0.04166666666666666667
                       * (-1.758243446661483480 + lnx)
                       * (-0.764428657272716373 + lnx)
                       * ( 0.723980571623507657 + lnx)
                       * ( 4.107554191916823640 + lnx);
    const double c5 = -0.0083333333333333333
                       * (-2.06563396085715900 + lnx)
                       * (-1.28459889470864700 + lnx)
                       * (-0.27583535756454143 + lnx)
                       * ( 1.33677371336239618 + lnx)
                       * ( 5.17537282427561550 + lnx);
    const double c6 = -0.0013888888888888889
                       * (-2.30814336454783200 + lnx)
                       * (-1.65846557706987300 + lnx)
                       * (-0.88768082560020400 + lnx)
                       * ( 0.17043847751371778 + lnx)
                       * ( 1.92135970115863890 + lnx)
                       * ( 6.22578557795474900 + lnx);
    const double c7 = -0.00019841269841269841
                       * (-2.5078657901291800 + lnx)
                       * (-1.9478900888958200 + lnx)
                       * (-1.3194837322612730 + lnx)
                       * (-0.5281322700249279 + lnx)
                       * ( 0.5913834939078759 + lnx)
                       * ( 2.4876819633378140 + lnx)
                       * ( 7.2648160783762400 + lnx);
    const double c8 = -0.00002480158730158730
                       * (-2.677341544966400 + lnx)
                       * (-2.182810448271700 + lnx)
                       * (-1.649350342277400 + lnx)
                       * (-1.014099048290790 + lnx)
                       * (-0.191366955370652 + lnx)
                       * ( 0.995403817918724 + lnx)
                       * ( 3.041323283529310 + lnx)
                       * ( 8.295966556941250 + lnx);
    const double c9 = -2.75573192239859e-6
                       * (-2.8243487670469080 + lnx)
                       * (-2.3798494322701120 + lnx)
                       * (-1.9143674728689960 + lnx)
                       * (-1.3814529102920370 + lnx)
                       * (-0.7294312810261694 + lnx)
                       * ( 0.1299079285269565 + lnx)
                       * ( 1.3873333251885240 + lnx)
                       * ( 3.5857258865210760 + lnx)
                       * ( 9.3214237073814600 + lnx);
    const double c10 = -2.75573192239859e-7
                       * (-2.9540329644556910 + lnx)
                       * (-2.5491366926991850 + lnx)
                       * (-2.1348279229279880 + lnx)
                       * (-1.6741881076349450 + lnx)
                       * (-1.1325949616098420 + lnx)
                       * (-0.4590034650618494 + lnx)
                       * ( 0.4399352987435699 + lnx)
                       * ( 1.7702236517651670 + lnx)
                       * ( 4.1231539047474080 + lnx)
                       * ( 10.342627908148680 + lnx);

    term1 = a*(c1+a*(c2+a*(c3+a*(c4+a*(c5+a*(c6+a*(c7+a*(c8+a*(c9+a*c10)))))))));
  }

  {
    /* Evaluate the sum.
     */
    const int nmax = 5000;
    double t = 1.0;
    int n;
    sum = 1.0;

    for(n=1; n<nmax; n++) {
      t *= -x/(n+1.0);
      sum += (a+1.0)/(a+n+1.0)*t;
      if(fabs(t/sum) < GSL_DBL_EPSILON) break;
    }

    if(n == nmax)
      stat_sum = GSL_EMAXITER;
    else
      stat_sum = GSL_SUCCESS;
  }

  term2 = (1.0 - term1) * a/(a+1.0) * x * sum;
  result->val  = term1 + term2;
  result->err  = GSL_DBL_EPSILON * (fabs(term1) + 2.0*fabs(term2));
  result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
  return stat_sum;
}


/* series for small a and x, but not defined for a == 0 */
static int
gamma_inc_series(double a, double x, gsl_sf_result * result)
{
  gsl_sf_result Q;
  gsl_sf_result G;
  const int stat_Q = gamma_inc_Q_series(a, x, &Q);
  const int stat_G = gsl_sf_gamma_e(a, &G);
  result->val = Q.val * G.val;
  result->err = fabs(Q.val * G.err) + fabs(Q.err * G.val);
  result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);

  return GSL_ERROR_SELECT_2(stat_Q, stat_G);
}


static int
gamma_inc_a_gt_0(double a, double x, gsl_sf_result * result)
{
  /* x > 0 and a > 0; use result for Q */
  gsl_sf_result Q;
  gsl_sf_result G;
  const int stat_Q = gsl_sf_gamma_inc_Q_e(a, x, &Q);
  const int stat_G = gsl_sf_gamma_e(a, &G);

  result->val = G.val * Q.val;
  result->err = fabs(G.val * Q.err) + fabs(G.err * Q.val);
  result->err += 2.0*GSL_DBL_EPSILON * fabs(result->val);

  return GSL_ERROR_SELECT_2(stat_G, stat_Q);
}


static int
gamma_inc_CF(double a, double x, gsl_sf_result * result)
{
  gsl_sf_result F;
  gsl_sf_result pre;
  const double am1lgx = (a-1.0)*log(x);
  const int stat_F = gamma_inc_F_CF(a, x, &F);
  const int stat_E = gsl_sf_exp_err_e(am1lgx - x, GSL_DBL_EPSILON*fabs(am1lgx), &pre);

  result->val = F.val * pre.val;
  result->err = fabs(F.err * pre.val) + fabs(F.val * pre.err);
  result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);

  return GSL_ERROR_SELECT_2(stat_F, stat_E);
}


/* evaluate Gamma(0,x), x > 0 */
#define GAMMA_INC_A_0(x, result) gsl_sf_expint_E1_e(x, result)


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_gamma_inc_e(const double a, const double x, gsl_sf_result * result)
{
  if(x < 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(x == 0.0) {
    return gsl_sf_gamma_e(a, result);
  }
  else if(a == 0.0)
  {
    return GAMMA_INC_A_0(x, result);
  }
  else if(a > 0.0)
  {
    return gamma_inc_a_gt_0(a, x, result);
  }
  else if(x > 0.25)
  {
    /* continued fraction seems to fail for x too small; otherwise
       it is ok, independent of the value of |x/a|, because of the
       non-oscillation in the expansion, i.e. the CF is
       un-conditionally convergent for a < 0 and x > 0
     */
    return gamma_inc_CF(a, x, result);
  }
  else if(fabs(a) < 0.5)
  {
    return gamma_inc_series(a, x, result);
  }
  else
  {
    /* a = fa + da; da >= 0 */
    const double fa = floor(a);
    const double da = a - fa;

    gsl_sf_result g_da;
    const int stat_g_da = ( da > 0.0 ? gamma_inc_a_gt_0(da, x, &g_da)
                                     : GAMMA_INC_A_0(x, &g_da));

    double alpha = da;
    double gax = g_da.val;

    /* Gamma(alpha-1,x) = 1/(alpha-1) (Gamma(a,x) - x^(alpha-1) e^-x) */
    do
    {
      const double shift = exp(-x + (alpha-1.0)*log(x));
      gax = (gax - shift) / (alpha - 1.0);
      alpha -= 1.0;
    } while(alpha > a);

    result->val = gax;
    result->err = 2.0*(1.0 + fabs(a))*GSL_DBL_EPSILON*fabs(gax);
    return stat_g_da;
  }

}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

double gsl_sf_gamma_inc(const double a, const double x)
{
   EVAL_RESULT(gsl_sf_gamma_inc_e(a, x, &result));
}
