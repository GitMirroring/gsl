/* specfunc/gamma/lnbeta.c
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
 * Provenance:
 *
 * The positive-argument double-precision algorithm (lnbeta_pos() and
 * the helpers lnbeta_dlgdiv(), lnbeta_dbcorr(), lnbeta_dgsmln()) is a
 * port of the routines DBETLN, DLGDIV, DBCORR and DGSMLN of the NSWC
 * Library of Mathematics Subroutines (A. H. Morris, Jr., Naval Surface
 * Warfare Center), including the DELTA coefficients of DLGDIV/DBCORR.
 * The NSWC library is a work of the U.S. Government and is in the
 * public domain.  The same algorithm is published as BETALN in
 *
 *   A. R. DiDonato and A. H. Morris, Jr., "Algorithm 708: Significant
 *   digit computation of the incomplete beta function ratios", ACM
 *   Trans. Math. Software 18 (1992) 360-373;
 *
 * no code was taken from the ACM (TOMS) distribution.  The
 * double-double refinement near log B = 0, the local expansion about
 * (1,1) and the reflection formulas for negative arguments are original
 * to GSL; their coefficients (Bernoulli numbers, zeta(n)/n) are
 * mathematical constants computed for GSL.
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

/*
 * For min(a,b) < 1, the double-precision result is recomputed in
 * double-double when the terms it is formed from exceed its magnitude
 * by more than this factor (near the log(B) = 0 contour)
 */
#define LNBETA_REFINE_AMP  2.0

/* radius and degree of the local expansion about (a,b) = (1,1) */
#define LNBETA_CORNER_RADIUS 0.04
#define LNBETA_CORNER_NMAX   30

/* smallest argument of the double-double Stirling series */
#define LNBETA_DD_STIRLING_MIN 30.0

/*
 * absolute error of the double-double evaluations relative to the
 * magnitude of their intermediate terms (observed ~1e-32; with margin)
 */
#define LNBETA_DD_REL_ERR 4.0e-32

static int lnbeta_pos(const double a0, const double b0, gsl_sf_result * result);
static int lnbeta_pos_dd(const gsl_dd x, const gsl_dd y, gsl_sf_result * result);
static double lnbeta_log_sin_pi(const gsl_dd x, double * sgn);
static int lnbeta_isnonposint(const double x);
static double lnbeta_small(const double a, const double b, double * t);
static double lnbeta_dlgdiv(const double a, const double b);
static double lnbeta_dbcorr(const double a0, const double b0);
static double lnbeta_dgamln(const double a);
static double lnbeta_dgsmln(const double a, const double b);
static void lnbeta_set_result(gsl_sf_result * r, const double val, const double factor);
static gsl_dd lnbeta_dd(const gsl_dd x0, const gsl_dd y0, double * mag);
static gsl_dd lnbeta_dd_abs_sin_pi(const gsl_dd x);
static gsl_dd lnbeta_neg_dd(const double a, const double b, const gsl_dd s, double * mag);
static gsl_dd lnbeta_dd_delta(const gsl_dd z);
static gsl_dd lnbeta_dd_stirling(const gsl_dd z);
static gsl_dd lnbeta_corner_11(const gsl_dd p, const gsl_dd q);

/*
 * Compute log|B(a,b)| and an absolute error estimate.
 */
int
gsl_sf_lnbeta_e(const double a, const double b, gsl_sf_result * result)
{
  double sgn;
  int status = gsl_sf_lnbeta_sgn_e(a, b, result, &sgn);

  if (sgn == -1)
    DOMAIN_ERROR(result);

  return status;
}

/*
 * Compute log|B(a,b)| and the sign of B(a,b).
 *
 * B(a,b) = Gamma(a) Gamma(b) / Gamma(a+b) has poles where a or b is a
 * non-positive integer (GSL_EDOM) and vanishes where only a+b is one
 * (log|B| = -inf: GSL_EDOM with sgn = 0).
 *
 * Positive arguments are handled by lnbeta_pos().  Otherwise, with
 * a <= b, the reflection formula Gamma(x) = pi / (sin(pi x) Gamma(1-x))
 * reduces B to a positive-argument Beta function and elementary
 * factors, s(x) = sin(pi x):
 *
 *   a < 0 < b, a+b > 0:  B(a,b) = pi / (b s(a) B(1-a, a+b))
 *   a < 0 < b, a+b < 0:  B(a,b) = s(a+b) / s(a) B(b, 1-a-b)
 *   a, b < 0:            B(a,b) = pi s(a+b) / (s(a) s(b) (1-a-b) B(1-a, 1-b))
 *
 * so no large log Gamma values cancel.  The sums a+b, 1-a, 1-a-b, 1-b
 * are formed exactly as double-double pairs, and their low parts are
 * carried to first order.  Near |B| = 1 the elementary terms cancel
 * against log B; when the cancellation exceeds LNBETA_REFINE_AMP the
 * result is recomputed in double-double by lnbeta_neg_dd().
 */
int
gsl_sf_lnbeta_sgn_e(const double a0, const double b0, gsl_sf_result * result, double * sgn)
{
  const double a = GSL_MIN(a0, b0);
  const double b = GSL_MAX(a0, b0);

  if (a > 0.0)
    {
      *sgn = 1.0;
      return lnbeta_pos(a, b, result);
    }
  else if (lnbeta_isnonposint(a) || lnbeta_isnonposint(b))
    {
      *sgn = 0.0;
      DOMAIN_ERROR(result);
    }
  else
    {
      const gsl_dd s = gsl_dd_two_sum(a, b);
      const gsl_dd one = { 1.0, 0.0 };
      gsl_dd ad, bd, x, y;
      gsl_sf_result r;
      double sa, ss, sb, la, ls, lb, t;
      int status;

      if (s.lo == 0.0 && lnbeta_isnonposint(s.hi))
        {
          *sgn = 0.0;          /* B = 0 */
          DOMAIN_ERROR(result);
        }

      ad.hi = a;
      ad.lo = 0.0;
      bd.hi = b;
      bd.lo = 0.0;
      la = lnbeta_log_sin_pi(ad, &sa);

      if (b > 0.0 && (s.hi > 0.0 || (s.hi == 0.0 && s.lo > 0.0)))
        {
          /* log|B| = log(pi) - log(b) - log|s(a)| - log B(1-a, a+b) */
          x = gsl_dd_sub(one, ad);
          status = lnbeta_pos_dd(x, s, &r);
          lb = log(b);
          result->val = ((M_LNPI - lb) - la) - r.val;
          t = M_LNPI + fabs(lb) + fabs(la);
          *sgn = sa;
        }
      else if (b > 0.0)
        {
          /* log|B| = log|s(a+b)| - log|s(a)| + log B(b, 1-a-b) */
          ls = lnbeta_log_sin_pi(s, &ss);
          y = gsl_dd_sub(one, s);
          status = lnbeta_pos_dd(bd, y, &r);
          result->val = (ls - la) + r.val;
          t = fabs(ls) + fabs(la);
          *sgn = ss * sa;
        }
      else
        {
          /* log|B| = log(pi) + log|s(a+b)| - log|s(a)| - log|s(b)| - log(1-a-b) - log B(1-a, 1-b) */
          double l1;

          ls = lnbeta_log_sin_pi(s, &ss);
          lb = lnbeta_log_sin_pi(bd, &sb);
          y = gsl_dd_sub(one, s);
          l1 = log(y.hi) + y.lo / y.hi;
          x = gsl_dd_sub(one, ad);
          status = lnbeta_pos_dd(x, gsl_dd_sub(one, bd), &r);
          result->val = ((((M_LNPI + ls) - la) - lb) - l1) - r.val;
          t = M_LNPI + fabs(ls) + fabs(la) + fabs(lb) + fabs(l1);
          *sgn = ss * sa * sb;
        }

      t += fabs(r.val);

      if (t > LNBETA_REFINE_AMP * fabs(result->val))
        {
          /* near |B| = 1: the terms cancel, recompute in double-double */
          double mag;
          const gsl_dd v = lnbeta_neg_dd(a, b, s, &mag);

          result->val = v.hi + v.lo;
          result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val) + LNBETA_DD_REL_ERR * mag;
        }
      else
        {
          result->err = r.err + 2.0 * GSL_DBL_EPSILON * t + 2.0 * GSL_DBL_EPSILON * fabs(result->val);
        }

      return status;
    }
}

double gsl_sf_lnbeta(const double a, const double b) { EVAL_RESULT(gsl_sf_lnbeta_e(a, b, &result)); }

/*********************************************************
 * INTERNAL ROUTINES FOR NEGATIVE ARGUMENTS              *
 *********************************************************/

/* return 1 if x is zero or a negative integer */
static int
lnbeta_isnonposint(const double x)
{
  return (x <= 0.0 && x == floor(x));
}

/*
 * Compute log|sin(pi x)| and its sign for x = x.hi + x.lo, using
 * sin(pi (x.hi + x.lo)) = sin(pi x.hi) + pi x.lo cos(pi x.hi) to first
 * order in the (tiny) low part, which matters near the zeros of sin.
 */
static double
lnbeta_log_sin_pi(const gsl_dd x, double * sgn)
{
  double v = gsl_sf_sin_pi(x.hi);

  if (x.lo != 0.0)
    v += M_PI * x.lo * gsl_sf_cos_pi(x.hi);

  *sgn = (v < 0.0) ? -1.0 : 1.0;

  return log(fabs(v));
}

/*
 * log B(x,y) for positive x = x.hi + x.lo, y = y.hi + y.lo: lnbeta_pos()
 * at the leading parts, corrected to first order for the low parts,
 *
 *   d log B / dx = psi(x) - psi(x+y),   d log B / dy = psi(y) - psi(x+y).
 */
static int
lnbeta_pos_dd(const gsl_dd x, const gsl_dd y, gsl_sf_result * result)
{
  int status = lnbeta_pos(x.hi, y.hi, result);

  if (x.lo != 0.0 || y.lo != 0.0)
    {
      const double pxy = gsl_sf_psi(x.hi + y.hi);
      const double c = x.lo * (gsl_sf_psi(x.hi) - pxy) + y.lo * (gsl_sf_psi(y.hi) - pxy);

      result->val += c;
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(c);
    }

  return status;
}

/*********************************************************
 * INTERNAL ROUTINES FOR POSITIVE-POSITIVE LNBETA        *
 *********************************************************/

/*
 * Positive-domain log-Beta dispatcher.
 *
 * The arguments are ordered so that a <= b, and the quadrant is divided
 * into the regions of the NSWC DBETLN algorithm (TOMS 708 BETALN),
 * evaluated in double precision.  For a < 1 the log(B) = 0 contour
 * passes through the region; there the double-precision result is
 * recomputed in double-double by lnbeta_dd() whenever the terms
 * it is formed from cancel by more than a factor LNBETA_REFINE_AMP.
 */
static int
lnbeta_pos(const double a0, const double b0, gsl_sf_result * result)
{
  double a = GSL_MIN(a0, b0);
  double b = GSL_MAX(a0, b0);
  double c, h, sn, u, v, w, z, val;
  int n;
  int i;

  if (a < 1.0)
    {
      double t;

      val = lnbeta_small(a, b, &t);

      if (t > LNBETA_REFINE_AMP * fabs(val))
        {
          gsl_dd ad, bd, r;
          double mag;

          ad.hi = a;
          ad.lo = 0.0;
          bd.hi = b;
          bd.lo = 0.0;
          r = lnbeta_dd(ad, bd, &mag);

          result->val = r.hi + r.lo;
          result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val) + LNBETA_DD_REL_ERR * mag;
        }
      else
        {
          result->val = val;
          result->err = 4.0 * GSL_DBL_EPSILON * t;
        }

      return GSL_SUCCESS;
    }

  /* 1 <= a < 10. */
  if (a < 10.0)
    {
      if (a <= 2.0)
        {
          if (b <= 2.0)
            {
              val = lnbeta_dgamln(a) + lnbeta_dgamln(b) - lnbeta_dgsmln(a, b);
              lnbeta_set_result(result, val, 16.0);
              return GSL_SUCCESS;
            }

          if (b >= 10.0)
            {
              val = lnbeta_dgamln(a) + lnbeta_dlgdiv(a, b);
              lnbeta_set_result(result, val, 8.0);
              return GSL_SUCCESS;
            }

          /* 2 < b < 10: reduce b to the small region. */
          w = 0.0;
        }
      else
        {
          /* 2 < a < 10: reduce a first. */
          n = (int)(a - 1.0);
          w = 1.0;
          for (i = 1; i <= n; ++i)
            {
              a -= 1.0;
              h = a / b;
              if (b > 1.0e3)
                w *= a / (1.0 + h);
              else
                w *= h / (1.0 + h);
            }

          w = log(w);

          if (b >= 10.0)
            {
              if (b <= 1.0e3)
                {
                  val = w + lnbeta_dgamln(a) + lnbeta_dlgdiv(a, b);
                  lnbeta_set_result(result, val, 12.0);
                  return GSL_SUCCESS;
                }

              /* b > 1000 */
              sn = (double) n;
              val = (w - sn * log(b)) + (lnbeta_dgamln(a) + lnbeta_dlgdiv(a, b));
              lnbeta_set_result(result, val, 16.0);
              return GSL_SUCCESS;
            }
        }

      /* reduce b when b < 10 */
      n = (int)(b - 1.0);
      z = 1.0;
      for (i = 1; i <= n; ++i)
        {
          b -= 1.0;
          z *= b / (a + b);
        }

      val = w + log(z) + (lnbeta_dgamln(a) + (lnbeta_dgamln(b) - lnbeta_dgsmln(a, b)));
      lnbeta_set_result(result, val, 16.0);
      return GSL_SUCCESS;
    }

  /* a >= 10: Stirling form with the DBCORR correction */
  w = lnbeta_dbcorr(a, b);
  h = a / b;
  c = h / (1.0 + h);
  u = -(a - 0.5) * log(c);
  v = b * gsl_sf_log1p(h);

  if (u <= v)
    val = (((-0.5 * log(b) + GSL_SF_LOG_SQRT_2PI) + w) - u) - v;
  else
    val = (((-0.5 * log(b) + GSL_SF_LOG_SQRT_2PI) + w) - v) - u;

  lnbeta_set_result(result, val, 12.0);
  return GSL_SUCCESS;
}

/*
 * Compute log B(a,b) in double precision for 0 < a < 1, a <= b, and in
 * *t the sum of the magnitudes of the terms it is formed from, so that
 * t / |log B| bounds the amplification of their rounding errors.
 *
 *   log B(a,b) = log Gamma(a) + d,   d = log Gamma(b) - log Gamma(a+b).
 *
 * For b >= 10, d is given by dlgdiv().  For b < 1, d is formed from
 * lngamma1p() and log1p(a/b) (see below).  For 1 <= b < 10, b is
 * reduced to b0 = b - n in [1,2) (exactly) when b >= 2 (TOMS 708 BETALN):
 *
 *   d = log Gamma(b0) - log Gamma(a+b0) - sum_{k<n} log1p(a / (b0+k)),
 *
 * so that every term is small and accurate relative to its magnitude.
 * a + b0 is formed exactly, and log Gamma(a+b0) is corrected for its
 * rounding to first order, d/dx log Gamma(x) = psi(x).
 */
static double
lnbeta_small(const double a, const double b, double * t)
{
  const double ga = lnbeta_dgamln(a);
  double d;

  if (b >= 10.0)
    {
      d = lnbeta_dlgdiv(a, b);
      *t = fabs(ga) + fabs(d);
    }
  else if (b < 1.0)
    {
      /*
       * log Gamma(x) = log Gamma(1+x) - log x for both terms:
       *   d = lngamma1p(b) - lngamma1p(a+b) + log1p(a/b),
       * avoiding the cancellation of -log b against log(a+b); the
       * rounding of a + b is corrected with d/dx lngamma1p(x) = psi(1+x)
       */
      const gsl_dd ab = gsl_dd_two_sum(a, b);
      const double gb = gsl_sf_lngamma1p(b);
      const double gab = gsl_sf_lngamma1p(ab.hi) + ((ab.lo != 0.0) ? gsl_sf_psi(1.0 + ab.hi) * ab.lo : 0.0);
      const double l = gsl_sf_log1p(a / b);

      d = (gb - gab) + l;
      *t = fabs(ga) + fabs(gb) + fabs(gab) + l;
    }
  else
    {
      const int n = (b >= 2.0) ? (int) b - 1 : 0;
      const double b0 = b - (double) n;     /* exact */
      const gsl_dd ab = gsl_dd_two_sum(a, b0);
      const double gb = lnbeta_dgamln(b0);
      const double gab = lnbeta_dgamln(ab.hi) + ((ab.lo != 0.0) ? gsl_sf_psi(ab.hi) * ab.lo : 0.0);
      double s = 0.0;
      int k;

      for (k = 0; k < n; ++k)
        s += gsl_sf_log1p(a / (b0 + (double) k));

      d = (gb - gab) - s;
      *t = fabs(ga) + fabs(gb) + fabs(gab) + s;
    }

  return ga + d;
}

/*
 * Coefficients of the NSWC asymptotic correction series used by
 * DLGDIV and DBCORR, scaled for the variable (10/x)^2.
 */
static const double DELTA[15] = {
  8.333333333333333333333333333333e-02,
 -2.777777777777777777752282e-05,
  7.93650793650793650791732130419e-08,
 -5.95238095238095232389839236182e-10,
  8.41750841850832853294451669990e-12,
 -1.91752691751854612334149171243e-13,
  6.41025640510325475730918472625e-15,
 -2.95506514125338232839867823991e-16,
  1.79643716359402238723287696452e-17,
 -1.39228964661627791231203060395e-18,
  1.33802855014020915603275339093e-19,
 -1.54246009867966094273710216533e-20,
  1.97701992980957427378133333333e-21,
 -2.34065664793997056856992642667e-22,
  1.71348014966398575409015466667e-23
};

/*
 * Compute the NSWC DLGDIV asymptotic representation of
 *
 *   log Gamma(b) - log Gamma(a+b)
 *
 * for a >= 0 and b >= 10.
 *
 * The logarithmic terms are arranged to reduce cancellation, and the
 * DELTA series supplies Delta(b) - Delta(a+b) without subtracting the
 * two corrections.
 */
static double
lnbeta_dlgdiv(const double a, const double b)
{
  double h, c, d, x, x2, w, t, u, v;
  double s[15];
  int j;

  if (a <= b)
    {
      h = a / b;
      c = h / (1.0 + h);
      x = 1.0 / (1.0 + h);
      d = b + (a - 0.5);
    }
  else
    {
      h = b / a;
      c = 1.0 / (1.0 + h);
      x = h / (1.0 + h);
      d = a + (b - 0.5);
    }

  x2 = x * x;
  s[0] = 1.0;
  for (j = 1; j < 15; ++j)
    s[j] = 1.0 + (x + x2 * s[j - 1]);

  t = (10.0 / b) * (10.0 / b);
  w = DELTA[14] * s[14];
  for (j = 14; j >= 1; --j)
    w = t * w + DELTA[j - 1] * s[j - 1];
  w *= c / b;

  u = d * gsl_sf_log1p(a / b);
  v = a * (log(b) - 1.0);

  return (u <= v) ? ((w - u) - v) : ((w - v) - u);
}

/*
 * Compute the asymptotic correction
 *
 *   Delta(a) + Delta(b) - Delta(a+b)
 *
 * for a,b >= 10, using the NSWC DBCORR representation.
 */
static double
lnbeta_dbcorr(const double a0, const double b0)
{
  const double a = GSL_MIN(a0, b0);
  const double b = GSL_MAX(a0, b0);
  const double h = a / b;
  const double c = h / (1.0 + h);
  const double x = 1.0 / (1.0 + h);
  const double x2 = x * x;
  double s[15];
  double w, t;
  int j;

  s[0] = 1.0;
  for (j = 1; j < 15; ++j)
    s[j] = 1.0 + (x + x2 * s[j - 1]);

  t = (10.0 / b) * (10.0 / b);
  w = DELTA[14] * s[14];
  for (j = 14; j >= 1; --j)
    w = t * w + DELTA[j - 1] * s[j - 1];
  w *= c / b;

  return gamma_stirling_delta(a) + w;
}

/*
 * Compute log Gamma(a) for a > 0.
 *
 * Small and moderate arguments are reduced to the lngamma1p and
 * lngamma kernels.  For a >= 10, use Stirling's formula with the
 * correction Delta(a) (stirling.c).
 */
static double
lnbeta_dgamln(const double a)
{
  if (a < 0.5)
    return gsl_sf_lngamma1p(a) - log(a);

  if (a <= 2.5)
    {
      const double x = (a < 1.0) ? (a - 0.5) - 0.5 : a - 1.0;
      return gsl_sf_lngamma1p(x);
    }

  if (a < 10.0)
    {
      const int n = (int)(a - 1.5);
      double x = a;
      double w = 1.0;
      int i;

      for (i = 0; i < n; ++i)
        {
          x -= 1.0;
          w *= x;
        }

      return gsl_sf_lngamma(x) + log(w);
    }

  return (GSL_SF_LOG_SQRT_2PI - 0.5 + gamma_stirling_delta(a)) + (a - 0.5) * (log(a) - 1.0);
}

/*
 * Compute log Gamma(a+b) for 1 <= a,b <= 2 using the NSWC DGSMLN
 * argument reduction and the lngamma1p kernel.
 */
static double
lnbeta_dgsmln(const double a, const double b)
{
  const double x = (a - 1.0) + (b - 1.0);

  if (x <= 0.5)
    return gsl_sf_lngamma1p(1.0 + x);

  if (x <= 1.5)
    return gsl_sf_lngamma1p(x) + gsl_sf_log1p(x);

  return gsl_sf_lngamma1p(x - 1.0) + log(x * (1.0 + x));
}

/* store a result and its branch-dependent absolute error estimate */
static void
lnbeta_set_result(gsl_sf_result * r, const double val, const double factor)
{
  r->val = val;
  r->err = factor * GSL_DBL_EPSILON * (1.0 + fabs(val));
}

/*********************************************************
 * DOUBLE-DOUBLE REFINEMENT NEAR |B(a,b)| = 1            *
 *********************************************************/

/*
 * Compute log B(x,y) in double-double for positive double-double
 * arguments x, y.
 *
 * With x <= y: near (1,1) use the local expansion lnbeta_corner_11().
 * Otherwise shift x and y to z = x + N >= 30 and y' = y + M >= 30 (N = 0,
 * M = 0 for arguments already >= 30):
 *
 *   log Gamma(x) = S(z) - log P,                P = prod_{k<N} (x+k),
 *   log Gamma(y) - log Gamma(x+y)
 *     = D(y') + log Q,                          Q = prod_{k<M} (x+y+k)/(y+k),
 *
 * with the Stirling form S(z) = log Gamma(z) (lnbeta_dd_stirling()) and
 * the Stirling difference
 *
 *   D(y') = -(y' - 1/2) log1p(x/y') - x log(x+y') + x + Delta(y') - Delta(x+y').
 *
 * The products are combined into one logarithm, log(Q/P).  *mag receives
 * the magnitude of the intermediate terms, which bounds the absolute
 * rounding error as LNBETA_DD_REL_ERR * (*mag).
 */
static gsl_dd
lnbeta_dd(const gsl_dd x0, const gsl_dd y0, double * mag)
{
  const gsl_dd one = { 1.0, 0.0 };
  const int swap = (x0.hi > y0.hi || (x0.hi == y0.hi && x0.lo > y0.lo));
  const gsl_dd x = swap ? y0 : x0;
  const gsl_dd y = swap ? x0 : y0;
  const gsl_dd xy = gsl_dd_add(x, y);
  gsl_dd P = one, R, z, yp, xyp, s, t;
  int k, N, M;

  if (fabs(x.hi - 1.0) < LNBETA_CORNER_RADIUS && fabs(y.hi - 1.0) < LNBETA_CORNER_RADIUS)
    {
      const gsl_dd p = gsl_dd_sub_d(x, 1.0);
      const gsl_dd q = gsl_dd_sub_d(y, 1.0);

      *mag = 2.0 * (fabs(p.hi) + fabs(q.hi));    /* |log1p(p+q)| + series terms */
      return lnbeta_corner_11(p, q);
    }

  /* log Gamma(x): shift x up to z >= LNBETA_DD_STIRLING_MIN */
  N = (x.hi < LNBETA_DD_STIRLING_MIN) ? (int) ceil(LNBETA_DD_STIRLING_MIN - x.hi) : 0;
  for (k = 0; k < N; ++k)
    P = gsl_dd_mul(P, gsl_dd_add_d(x, (double) k));
  z = gsl_dd_add_d(x, (double) N);

  /* log Gamma(y) - log Gamma(x+y): shift y up to y' >= LNBETA_DD_STIRLING_MIN */
  M = (y.hi < LNBETA_DD_STIRLING_MIN) ? (int) ceil(LNBETA_DD_STIRLING_MIN - y.hi) : 0;
  {
    gsl_dd num = one, den = P;

    for (k = 0; k < M; ++k)
      {
        num = gsl_dd_mul(num, gsl_dd_add_d(xy, (double) k));
        den = gsl_dd_mul(den, gsl_dd_add_d(y, (double) k));
      }

    R = gsl_dd_div(num, den);    /* Q / P */
    yp = gsl_dd_add_d(y, (double) M);
    xyp = gsl_dd_add_d(xy, (double) M);
  }

  s = lnbeta_dd_stirling(z);

  /* D(y') */
  {
    gsl_dd d;

    d = gsl_dd_mul(gsl_dd_sub_d(yp, 0.5), gsl_dd_log1p(gsl_dd_div(x, yp)));
    d = gsl_dd_add(d, gsl_dd_mul(gsl_dd_log(xyp), x));
    d = gsl_dd_sub(d, x);
    t = gsl_dd_sub(gsl_dd_sub(lnbeta_dd_delta(yp), lnbeta_dd_delta(xyp)), d);
  }

  {
    const gsl_dd lr = gsl_dd_log(R);

    *mag = fabs(s.hi) + fabs(t.hi) + fabs(lr.hi);
    return gsl_dd_add(gsl_dd_add(s, t), lr);
  }
}

/*
 * log|B(a,b)| in double-double for a <= b, a < 0, a+b = s not a
 * non-positive integer, from the reflection formulas of
 * gsl_sf_lnbeta_sgn_e(), with the elementary factors combined into one
 * logarithm.  *mag as for lnbeta_dd().
 */
static gsl_dd
lnbeta_neg_dd(const double a, const double b, const gsl_dd s, double * mag)
{
  const gsl_dd pi = { M_PI, GSL_SF_PI_LO };
  const gsl_dd one = { 1.0, 0.0 };
  const gsl_dd sa = lnbeta_dd_abs_sin_pi(gsl_dd_two_sum(a, 0.0));
  gsl_dd ad, bd;

  ad.hi = a;
  ad.lo = 0.0;
  bd.hi = b;
  bd.lo = 0.0;

  if (b > 0.0 && (s.hi > 0.0 || (s.hi == 0.0 && s.lo > 0.0)))
    {
      /* log(pi / (b |s(a)|)) - log B(1-a, a+b) */
      const gsl_dd le = gsl_dd_log(gsl_dd_div(pi, gsl_dd_mul_d(sa, b)));
      const gsl_dd r = lnbeta_dd(gsl_dd_sub(one, ad), s, mag);

      *mag += fabs(le.hi);
      return gsl_dd_sub(le, r);
    }
  else if (b > 0.0)
    {
      /* log(|s(a+b)| / |s(a)|) + log B(b, 1-a-b) */
      const gsl_dd le = gsl_dd_log(gsl_dd_div(lnbeta_dd_abs_sin_pi(s), sa));
      const gsl_dd r = lnbeta_dd(bd, gsl_dd_sub(one, s), mag);

      *mag += fabs(le.hi);
      return gsl_dd_add(le, r);
    }
  else
    {
      /* log(pi |s(a+b)| / (|s(a)| |s(b)| (1-a-b))) - log B(1-a, 1-b) */
      const gsl_dd num = gsl_dd_mul(pi, lnbeta_dd_abs_sin_pi(s));
      const gsl_dd den = gsl_dd_mul(gsl_dd_mul(sa, lnbeta_dd_abs_sin_pi(bd)), gsl_dd_sub(one, s));
      const gsl_dd le = gsl_dd_log(gsl_dd_div(num, den));
      const gsl_dd r = lnbeta_dd(gsl_dd_sub(one, ad), gsl_dd_sub(one, bd), mag);

      *mag += fabs(le.hi);
      return gsl_dd_sub(le, r);
    }
}

/*
 * |sin(pi x)| in double-double for x = x.hi + x.lo not an integer:
 * with r = x - n, n the nearest integer (|r| <= 1/2, formed exactly),
 * |sin(pi x)| = |u| sinc(u), u = pi r, from the Taylor series of
 * sinc(u) = sum_k (-u^2)^k / (2k+1)!, truncated after k = 18
 * (remainder below 1e-38 for |u| <= pi/2).
 */
static gsl_dd
lnbeta_dd_abs_sin_pi(const gsl_dd x)
{
  const gsl_dd pi = { M_PI, GSL_SF_PI_LO };
  const double n = floor(x.hi + 0.5);
  const gsl_dd r = gsl_dd_two_sum(x.hi - n, x.lo);    /* x.hi - n is exact */
  const gsl_dd u = gsl_dd_mul(pi, r);
  const gsl_dd u2 = gsl_dd_mul(u, u);
  gsl_dd term = { 1.0, 0.0 }, sum = { 1.0, 0.0 };
  int k;

  for (k = 1; k <= 18; ++k)
    {
      term = gsl_dd_div_d(gsl_dd_mul(term, u2), -(double) ((2 * k) * (2 * k + 1)));
      sum = gsl_dd_add(sum, term);
    }

  sum = gsl_dd_mul(sum, u);
  if (sum.hi < 0.0)
    {
      sum.hi = -sum.hi;
      sum.lo = -sum.lo;
    }

  return sum;
}

/*
 * Stirling correction Delta(z) = log Gamma(z) - [(z-1/2) log z - z + log sqrt(2 pi)]
 * in double-double for z >= LNBETA_DD_STIRLING_MIN, from the asymptotic
 * series sum_{k=1}^{15} B_{2k} / (2k (2k-1) z^{2k-1}); the first omitted
 * term is below 1e-38 for z >= 30.
 */
static gsl_dd
lnbeta_dd_delta(const gsl_dd z)
{
  /* B_{2k} / (2k (2k-1)), k = 1, ..., 15, as double-double pairs */
  static const double c[15][2] = {
    {  8.33333333333333287e-02,  4.62592926927148533e-18 },
    { -2.77777777777777788e-03,  1.06010879087471541e-19 },
    {  7.93650793650793650e-04,  6.88382331736828211e-22 },
    { -5.95238095238095292e-04,  5.36938218754726024e-20 },
    {  8.41750841750841714e-04,  3.68701748892376936e-20 },
    { -1.91752691752691763e-03,  1.06757027768724749e-19 },
    {  6.41025641025641003e-03,  2.22400445638052172e-19 },
    { -2.95506535947712423e-02,  4.86176095750885531e-19 },
    {  1.79644372368830574e-01, -6.40160048271094580e-19 },
    { -1.39243221690590113e+00,  1.58370569892303027e-17 },
    {  1.34028640441683926e+01, -6.15411410199396641e-16 },
    { -1.56848284626002027e+02,  9.39182314171538895e-15 },
    {  2.19310333333333347e+03, -1.33392556260029476e-13 },
    { -3.61087712537249899e+04,  5.89758335351436479e-13 },
    {  6.91472268851313042e+05,  2.55852963051579989e-11 }
  };
  const gsl_dd one = { 1.0, 0.0 };
  const gsl_dd w = gsl_dd_div(one, z);
  const gsl_dd w2 = gsl_dd_mul(w, w);
  gsl_dd p, ck;
  int k;

  p.hi = c[14][0];
  p.lo = c[14][1];
  for (k = 13; k >= 0; --k)
    {
      ck.hi = c[k][0];
      ck.lo = c[k][1];
      p = gsl_dd_add(gsl_dd_mul(p, w2), ck);
    }

  return gsl_dd_mul(p, w);
}

/*
 * log Gamma(z) = (z - 1/2) log z - z + log sqrt(2 pi) + Delta(z)
 * in double-double for z >= LNBETA_DD_STIRLING_MIN.
 */
static gsl_dd
lnbeta_dd_stirling(const gsl_dd z)
{
  const gsl_dd log_sqrt_2pi = { GSL_SF_LOG_SQRT_2PI, GSL_SF_LOG_SQRT_2PI_LO };
  gsl_dd s = gsl_dd_mul(gsl_dd_sub_d(z, 0.5), gsl_dd_log(z));

  s = gsl_dd_sub(s, z);
  s = gsl_dd_add(s, log_sqrt_2pi);

  return gsl_dd_add(s, lnbeta_dd_delta(z));
}

/*
 * Local expansion of log B(1+p, 1+q) about p = q = 0,
 *
 *   log B(1+p,1+q) = -log1p(p+q)
 *                    + sum_{n>=2} (-1)^n zeta(n)/n (p^n + q^n - (p+q)^n),
 *
 * evaluated through n = LNBETA_CORNER_NMAX in double-double arithmetic to
 * preserve relative accuracy near log(B) = 0, for |p|, |q| <
 * LNBETA_CORNER_RADIUS.  The truncation error is below
 * (2 LNBETA_CORNER_RADIUS)^(LNBETA_CORNER_NMAX+1), ~1e-34, so the double-double
 * rounding (~1e-32 relative to the largest term) dominates.
 */
static gsl_dd
lnbeta_corner_11(const gsl_dd p, const gsl_dd q)
{
  /* zeta(n)/n, n = 2, ..., LNBETA_CORNER_NMAX, as double-double pairs */
  static const double zn[LNBETA_CORNER_NMAX - 1][2] = {
    { 8.22467033424113203e-01, 1.52033617519923808e-17 },
    { 4.00685634386531431e-01, -2.25074704248750415e-18 },
    { 2.70580808427784536e-01, 1.18712801071384119e-17 },
    { 2.07385551028673981e-01, 4.09976732862181262e-18 },
    { 1.69557176997408188e-01, 2.23938513301672378e-18 },
    { 1.44049896768846108e-01, 9.62314008523255486e-18 },
    { 1.25509669524743045e-01, -2.52146853846723048e-18 },
    { 1.11334265869564686e-01, 4.64399057258292406e-18 },
    { 1.00099457512781806e-01, 2.61024048595832832e-18 },
    { 9.09540171458290414e-02, 8.30670545769188462e-19 },
    { 8.33538405461090037e-02, 2.96383260365264211e-19 },
    { 7.69325164113521948e-02, -3.29003560191811977e-18 },
    { 7.14329462953613298e-02, 6.27880602419149916e-18 },
    { 6.66687058824204648e-02, 3.22958607599663057e-18 },
    { 6.25009551412130382e-02, 2.55109946401931504e-18 },
    { 5.88239786586845850e-02, -2.69129013419663566e-18 },
    { 5.55557676274036141e-02, -3.02618648498309645e-18 },
    { 5.26316793796166582e-02, 2.52384370247121503e-18 },
    { 5.00000476981016934e-02, 2.78944182644587960e-19 },
    { 4.76190703301422255e-02, 2.47963426842933546e-18 },
    { 4.54545562932046690e-02, 4.38293177455007562e-19 },
    { 4.34782660530402612e-02, -1.84622298803959428e-18 },
    { 4.16666691503412082e-02, 2.30817468724826604e-18 },
    { 4.00000011921401374e-02, 3.14569061393772909e-18 },
    { 3.84615390346751823e-02, 3.39272042239591675e-18 },
    { 3.70370373129893238e-02, 1.77099324149498771e-18 },
    { 3.57142858473333547e-02, 3.37720268655954160e-18 },
    { 3.44827586849193041e-02, -3.25998692705954774e-18 },
    { 3.33333333643775834e-02, -2.29368273689617943e-18 }
  };
  const gsl_dd r = gsl_dd_add(p, q);
  gsl_dd termp = p, termq = q, termr = r;
  gsl_dd sum, t, c, u;
  int n;

  sum = gsl_dd_mul_d(gsl_dd_log1p(r), -1.0);

  for (n = 2; n <= LNBETA_CORNER_NMAX; ++n)
    {
      termp = gsl_dd_mul(termp, p);
      termq = gsl_dd_mul(termq, q);
      termr = gsl_dd_mul(termr, r);
      t = gsl_dd_sub(gsl_dd_add(termp, termq), termr);
      c.hi = zn[n - 2][0];
      c.lo = zn[n - 2][1];
      u = gsl_dd_mul(t, c);
      if (n & 1)
        sum = gsl_dd_sub(sum, u);
      else
        sum = gsl_dd_add(sum, u);
    }

  return sum;
}
