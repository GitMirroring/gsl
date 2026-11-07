/* specfunc/hyperg/hyperg_2F1.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
 * Copyright (C) 2009 Brian Gough
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
 * References:
 *
 *   N. Michel and M. V. Stoitsov, "Fast computation of the Gauss
 *   hypergeometric function with all its parameters complex with
 *   application to the Poschl-Teller-Ginocchio potential wave
 *   functions", Comput. Phys. Commun. 178 (2008) 535-551 -- choice of
 *   the transformation by the modulus of the transformed variable, and
 *   the evaluation of the connection formulas uniformly in the distance
 *   of c-a-b (or b-a) from an integer;
 *
 *   J. W. Pearson, S. Olver and M. A. Porter, "Numerical methods for the
 *   computation of the confluent and Gauss hypergeometric functions",
 *   Numer. Algorithms 74 (2017) 821-866;
 *
 *   NIST Digital Library of Mathematical Functions, 15.2, 15.4, 15.8;
 *   M. Abramowitz and I. A. Stegun, Handbook of Mathematical Functions,
 *   15.1, 15.3.
 *
 * No external code was used.
 */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_dd.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_sincos_pi.h>
#include <gsl/gsl_sf_hyperg.h>

#include "error.h"

/*
 * Accept a result without trying alternative transformations when the
 * sum of the magnitudes of its terms exceeds its value by at most this
 * factor (no significant cancellation: the error is rounding only).
 */
#define HYP2F1_GOOD_AMP     8.0

/* ... and its estimated relative error is below this */
#define HYP2F1_GOOD_REL     (64.0 * GSL_DBL_EPSILON)

/* re-sum a series in double-double when its cancellation factor exceeds this */
#define HYP2F1_DD_AMP       1.0e3

/* ... or when the best estimated relative error exceeds this, provided
   the series variable is at most HYP2F1_DD_VAR (a short series: the
   double-double sum costs about 30 ns per term) */
#define HYP2F1_DD_REL       (64.0 * GSL_DBL_EPSILON)
#define HYP2F1_DD_VAR       0.8

/* rounding error of the k-th term of a hypergeometric series, in units
   of eps times its magnitude (random-walk growth, measured) */
#define HYP2F1_TERM_ERR(k)  (2.0 + 6.0 * sqrt((k) + 1.0))

/* the same summed over the terms, sum_k |t_k| (2 + 8 sqrt(k+1)) (a
   slightly larger constant: measured on series dominated by late
   terms), bounded by Cauchy-Schwarz with sabs = sum |t_k|, skabs =
   sum (k+1) |t_k|; for a series whose weight lies in its first terms
   this is much smaller than HYP2F1_TERM_ERR(last k) * sabs.  Once k
   exceeds the parameters, a+k, b+k, c+k round in the same direction
   over long stretches of k, and the drift grows linearly: skabs */
#define HYP2F1_SUM_ERR(sabs, skabs) \
  (2.0 * (sabs) + GSL_MAX_DBL(8.0 * sqrt((sabs) * (skabs)), (skabs)))

/* maximum number of series terms */
#define HYP2F1_MAXITER      200000

/*
 * bound on the relative error of gsl_sf_gamma and gsl_sf_gammainv
 * (measured: below 4.4 ulp), used where their own estimate is larger
 */
#define HYP2F1_GAMMA_REL    (8.0 * GSL_DBL_EPSILON)

/* argument from which lnGamma(x+e)-lnGamma(x) uses Stirling's series */
#define HYP2F1_STIRLING_MIN 12.0

/*
 * A value represented as v * exp(l), with absolute error e on v.  The
 * scaling keeps large prefactors such as (1-x)^(-a) and Gamma function
 * ratios from overflowing before the final product is formed.
 */
typedef struct
{
  double v;
  double e;
  double l;
  double amp;   /* cancellation: (sum of |terms|) / |value| */
} hyp2f1_scaled;

/*
 * The upper parameters enter the algorithms only through the products
 * (a+t)(b+t) and sums (a+t)+(b+t).  A pair is either real (a, b), or
 * the complex conjugate pair a = x + iy, b = x - iy of the _conj
 * functions, for which these are real:
 *
 *   (a+t)(b+t) = (x+t)^2 + y^2,   (a+t)+(b+t) = 2(x+t).
 *
 * The real products are formed as (a+t)*(b+t), which is exactly zero
 * when a series terminates.
 */
typedef struct
{
  int conj;
  double x;   /* a (real pair) or Re a (conjugate pair) */
  double y;   /* b (real pair) or Im a (conjugate pair) */
  double xlo; /* rounding residuals of transformed parameters: */
  double ylo; /* the exact values are x+xlo, y+ylo (ylo = 0 if conj) */
} hyp2f1_pair;

static int hyp2f1_scaled_eval(const double a, const double b, const double c,
                              const double x, hyp2f1_scaled * r);

/*****************************************
 * small helpers                         *
 *****************************************/

/* x is a non-positive integer (exact test) */
static int
hyp2f1_isnpint(const double x)
{
  return (x <= 0.0 && x == floor(x));
}

/* relative error of a scaled result */
static double
hyp2f1_relerr(const hyp2f1_scaled * r)
{
  if (r->v == 0.0)
    return (r->e == 0.0) ? 0.0 : GSL_POSINF;
  return r->e / fabs(r->v);
}

/* r <- r * p * exp(lp), p with relative error pe */
static void
hyp2f1_scale(hyp2f1_scaled * r, const double p, const double pe, const double lp)
{
  r->e = fabs(p) * r->e + pe * fabs(p * r->v);
  r->v *= p;
  r->l += lp;
}

/*
 * Gamma(num[0]) ... Gamma(num[nn-1]) / (Gamma(den[0]) ... Gamma(den[nd-1]))
 * as m * exp(l) with relative error *rel.  The numerator arguments must
 * not be poles; a pole in the denominator gives an exact zero.  The
 * Gamma functions are used directly when no factor can overflow, since
 * gsl_sf_gamma and gsl_sf_gammainv are accurate to a few ulp; otherwise
 * log|Gamma| is used.
 */
static int
hyp2f1_gamma_ratio(const double * num, const int nn,
                   const double * den, const int nd,
                   double * m, double * l, double * rel)
{
  int i, direct = 1;

  for (i = 0; i < nn; ++i)
    if (fabs(num[i]) > 160.0)
      direct = 0;
  for (i = 0; i < nd; ++i)
    {
      if (hyp2f1_isnpint(den[i]))
        {
          *m = 0.0;
          *l = 0.0;
          *rel = 0.0;
          return GSL_SUCCESS;
        }
      if (fabs(den[i]) > 160.0)
        direct = 0;
    }

  if (direct)
    {
      double v = 1.0, re = 0.0;
      gsl_sf_result g;

      for (i = 0; i < nn; ++i)
        {
          int stat = gsl_sf_gamma_e(num[i], &g);
          if (stat != GSL_SUCCESS)
            return stat;
          v *= g.val;
          re += GSL_MIN_DBL(g.err / fabs(g.val), HYP2F1_GAMMA_REL);
        }
      for (i = 0; i < nd; ++i)
        {
          gsl_sf_gammainv_e(den[i], &g);
          v *= g.val;
          re += (g.val != 0.0) ? GSL_MIN_DBL(g.err / fabs(g.val), HYP2F1_GAMMA_REL) : 0.0;
        }

      if (v != 0.0 && gsl_finite(v) && fabs(v) > 1.0e-280 && fabs(v) < 1.0e280)
        {
          *m = v;
          *l = 0.0;
          *rel = re + (nn + nd) * GSL_DBL_EPSILON;
          return GSL_SUCCESS;
        }
    }

  {
    double sgn = 1.0, L = 0.0, Le = 0.0, s;
    gsl_sf_result g;

    for (i = 0; i < nn; ++i)
      {
        int stat = gsl_sf_lngamma_sgn_e(num[i], &g, &s);
        if (stat != GSL_SUCCESS)
          return stat;
        L += g.val;
        Le += g.err;
        sgn *= s;
      }
    for (i = 0; i < nd; ++i)
      {
        gsl_sf_lngamma_sgn_e(den[i], &g, &s);
        L -= g.val;
        Le += g.err;
        sgn *= s;
      }

    *m = sgn;
    *l = L;
    *rel = Le + 2.0 * GSL_DBL_EPSILON * fabs(L);
    return GSL_SUCCESS;
  }
}

/* (a+t)(b+t) */
static double
hyp2f1_pair_prod(const hyp2f1_pair * p, const double t)
{
  if (p->conj)
    return (p->x + t + p->xlo) * (p->x + t + p->xlo) + p->y * p->y;
  return (p->x + t + p->xlo) * (p->y + t + p->ylo);
}

/* (a+t)+(b+t) */
static double
hyp2f1_pair_sum(const hyp2f1_pair * p, const double t)
{
  if (p->conj)
    return 2.0 * (p->x + t + p->xlo);
  return (p->x + t + p->xlo) + (p->y + t + p->ylo);
}

/* min(Re a, Re b), for the sign changes of a series */
static double
hyp2f1_pair_minre(const hyp2f1_pair * p)
{
  return p->conj ? p->x : GSL_MIN_DBL(p->x, p->y);
}

/* the pair (c-a, c-b) of Euler's transformation */
static hyp2f1_pair
hyp2f1_pair_euler(const hyp2f1_pair * p, const double c)
{
  hyp2f1_pair q;
  const gsl_dd cx = gsl_dd_two_sum(c, -p->x);
  q.conj = p->conj;
  q.x = cx.hi;
  q.xlo = cx.lo - p->xlo;
  if (p->conj)
    {
      q.y = -p->y;
      q.ylo = 0.0;
    }
  else
    {
      const gsl_dd cy = gsl_dd_two_sum(c, -p->y);
      q.y = cy.hi;
      q.ylo = cy.lo - p->ylo;
    }
  return q;
}

/*
 * First-order correction for the rounding of Gamma function arguments:
 * Gamma(x + dx) = Gamma(x) (1 + psi(x) dx + ...), where x is a computed
 * parameter combination such as c-a and dx its rounding residual
 * (exact, from TwoSum).  The residuals are of order ulp(x), so the
 * first-order term is exact to working precision; it matters when the
 * Gamma factors are sensitive (large arguments).  Returns the factor
 * multiplying a ratio with x in the numerator (sign +1) or in the
 * denominator (sign -1).
 */
/* psi(x) to a few digits (all the correction needs) */
static double
hyp2f1_psi_approx(double x)
{
  double r = 0.0, x2;

  if (x < 0.5)
    {
      /* reflection: psi(x) = psi(1-x) - pi cot(pi x) */
      gsl_sf_result sx, cx;
      gsl_sf_sin_pi_e(x, &sx);
      gsl_sf_cos_pi_e(x, &cx);
      r = -M_PI * cx.val / sx.val;
      x = 1.0 - x;
    }

  while (x < 6.0)
    {
      r -= 1.0 / x;
      x += 1.0;
    }

  x2 = 1.0 / (x * x);
  return r + log(x) - 0.5 / x - x2 * (1.0 / 12.0 - x2 * (1.0 / 120.0 - x2 / 252.0));
}

static double
hyp2f1_gamma_arg_corr(const double x, const double dx, const double sign)
{
  if (dx == 0.0 || hyp2f1_isnpint(x))
    return 1.0;
  return 1.0 + sign * hyp2f1_psi_approx(x) * dx;
}

/*****************************************
 * Gauss series                          *
 *****************************************/

/*
 * Sum the Gauss series
 *
 *   F(a,b;c;z) = sum_k (a)_k (b)_k / ((c)_k k!) z^k.
 *
 * Either |z| < 1, or the series terminates (a or b a non-positive
 * integer); c must not be a non-positive integer reached before the
 * termination.  The term ratio tends to z, so summation stops once a
 * term is negligible and the ratio has settled below one in magnitude;
 * the remainder is bounded by the geometric tail.  The error estimate
 * is proportional to the sum of the magnitudes of the terms, so that
 * cancellation between terms is reflected in it.
 */
static int
hyp2f1_series_pair(const hyp2f1_pair * ab, const double c,
                   const double z, const double dpab, const double dpc,
                   const double dzr, hyp2f1_scaled * r)
{
  /*
   * Summation may stop only past kcross, where none of a+k, b+k, c+k
   * can change sign any more: before it a term can be momentarily small
   * (a factor near zero) and a small denominator c+k can make the next
   * terms large again.  Past kcross the term ratio approaches z
   * monotonically, and the remainder is bounded by a geometric tail.
   */
  const double kcross = GSL_MAX_DBL(-hyp2f1_pair_minre(ab), GSL_MAX_DBL(-c, 0.0)) + 2.0;
  const int terminating = !ab->conj && ((hyp2f1_isnpint(ab->x) && ab->xlo == 0.0)
                                        || (hyp2f1_isnpint(ab->y) && ab->ylo == 0.0));
  double t = 1.0, s = 1.0, sabs = 1.0;
  const int track = (dpab > 0.0 || dpc > 0.0);
  double sens = 0.0, serr = 0.0, skabs = 0.0;   /* see below */
  double k;

  r->l = 0.0;

  for (k = 0.0; k < HYP2F1_MAXITER; k += 1.0)
    {
      const double num = hyp2f1_pair_prod(ab, k);
      double ratio;

      if (num == 0.0)
        {
          /* exact termination (checked before dividing: c+k may vanish too) */
          r->v = s;
          r->e = GSL_DBL_EPSILON * HYP2F1_SUM_ERR(sabs, skabs + 1.0) + serr + dzr * skabs;
          r->amp = sabs / fabs(s);
          return GSL_SUCCESS;
        }

      ratio = num / ((c + k) * (k + 1.0)) * z;
      t *= ratio;

      /* parameters and variable obtained by a transformation carry
         rounding errors (absolute dpab for a and b, dpc for c, relative
         dzr for z), which perturb term k by the relative amount sens
         + k dzr */
      if (track)
        {
          if (dpab > 0.0)
            {
              const double dnum = ab->conj ? 2.0 * fabs(ab->x + k + ab->xlo)
                                           : fabs(ab->x + k + ab->xlo) + fabs(ab->y + k + ab->ylo);
              sens += dpab * dnum / fabs(num);
            }
          if (dpc > 0.0)
            sens += dpc / fabs(c + k);
          serr += fabs(t) * sens;
        }

      s += t;
      sabs += fabs(t);
      skabs += fabs(t) * (k + 1.0);

      /* the exceptional cases below all fail this test */
      if (!(fabs(t) > 0.25 * GSL_DBL_EPSILON * fabs(s) && fabs(t) <= GSL_DBL_MAX))
        {
          const double rr = GSL_MAX_DBL(fabs(ratio), fabs(z));

          if (!(fabs(t) <= GSL_DBL_MAX))   /* overflow or NaN */
            {
              r->v = GSL_POSINF;
              r->e = GSL_POSINF;
              r->amp = GSL_POSINF;
              return GSL_EOVRFLW;
            }

          if (t == 0.0)
            {
              /* underflow of the terms */
              r->v = s;
              r->e = GSL_DBL_EPSILON * HYP2F1_SUM_ERR(sabs, skabs + 1.0) + serr + dzr * skabs;
              r->amp = sabs / fabs(s);
              return GSL_SUCCESS;
            }

          /* a terminating series is always summed to the end: its terms
             can grow again after a small one when |z| > 1 */
          if (!terminating && rr < 0.999 && k > kcross)
            {
              r->v = s;
              r->e = GSL_DBL_EPSILON * HYP2F1_SUM_ERR(sabs, skabs + 1.0) + serr + dzr * skabs
                     + fabs(t) * (1.0 + sens + (k + 1.0) * dzr) * rr / (1.0 - rr);
              r->amp = sabs / fabs(s);
              return GSL_SUCCESS;
            }
        }
    }

  r->v = s;
  r->e = GSL_DBL_EPSILON * sabs + serr + dzr * skabs + fabs(t);
  r->amp = sabs / fabs(s);
  return GSL_EMAXITER;
}

static int
hyp2f1_series(const double a, const double b, const double c,
              const double z, const double dpab, const double dpc,
              const double dzr, hyp2f1_scaled * r)
{
  hyp2f1_pair ab;
  ab.conj = 0;
  ab.x = a;
  ab.y = b;
  ab.xlo = 0.0;
  ab.ylo = 0.0;
  return hyp2f1_series_pair(&ab, c, z, dpab, dpc, dzr, r);
}

/*
 * The Gauss series summed in double-double arithmetic, for parameters
 * and variable given as double-double values (so that combinations such
 * as c-b and z/(z-1) are exact).  Used only when the double sum shows
 * strong cancellation (see HYP2F1_DD_AMP): the dd rounding errors,
 * about 1e-32 of the terms, then leave the result correct to working
 * precision for a cancellation factor up to about 1e15.
 */
static int
hyp2f1_series_dd(const gsl_dd a, const gsl_dd b, const gsl_dd c, const gsl_dd z,
                 hyp2f1_scaled * r)
{
  const double kcross = GSL_MAX_DBL(GSL_MAX_DBL(-a.hi, -b.hi), GSL_MAX_DBL(-c.hi, 0.0)) + 2.0;
  const double az = fabs(z.hi);
  gsl_dd t, s;
  double sabs = 1.0, tail = 0.0, k;

  t.hi = 1.0;
  t.lo = 0.0;
  s = t;
  r->l = 0.0;

  for (k = 0.0; k < HYP2F1_MAXITER; k += 1.0)
    {
      const gsl_dd num = gsl_dd_mul(gsl_dd_add_d(a, k), gsl_dd_add_d(b, k));
      gsl_dd den, ratio;

      if (num.hi == 0.0)
        break;                  /* exact termination */

      den = gsl_dd_mul_d(gsl_dd_add_d(c, k), k + 1.0);
      ratio = gsl_dd_mul(gsl_dd_div(num, den), z);
      t = gsl_dd_mul(t, ratio);
      s = gsl_dd_add(s, t);
      sabs += fabs(t.hi);

      if (!gsl_finite(t.hi))
        return GSL_EOVRFLW;

      if (t.hi == 0.0)
        break;

      if (az < 1.0 && fabs(t.hi) < 1.0e-20 * fabs(s.hi) && k > kcross
          && GSL_MAX_DBL(fabs(ratio.hi), az) < 0.999)
        {
          const double rr = GSL_MAX_DBL(fabs(ratio.hi), az);
          tail = fabs(t.hi) * rr / (1.0 - rr);
          break;
        }
    }

  if (k >= HYP2F1_MAXITER)
    return GSL_EMAXITER;

  r->v = s.hi + s.lo;
  r->e = 1.0e-31 * HYP2F1_TERM_ERR(k) * sabs + tail + 0.5 * GSL_DBL_EPSILON * fabs(r->v);
  r->amp = 1.0;
  return GSL_SUCCESS;
}

/* the Gauss series for a conjugate pair a = x+iy, b = x-iy, in double-double */
static int
hyp2f1_series_conj_dd(const gsl_dd x, const double y, const gsl_dd c, const double z,
                      hyp2f1_scaled * r)
{
  const double kcross = GSL_MAX_DBL(-x.hi, GSL_MAX_DBL(-c.hi, 0.0)) + 2.0;
  const double az = fabs(z);
  gsl_dd t, s, y2;
  double sabs = 1.0, tail = 0.0, k;

  t.hi = 1.0;
  t.lo = 0.0;
  s = t;
  y2 = gsl_dd_two_prod(y, y);
  r->l = 0.0;

  for (k = 0.0; k < HYP2F1_MAXITER; k += 1.0)
    {
      const gsl_dd xk = gsl_dd_add_d(x, k);
      const gsl_dd num = gsl_dd_add(gsl_dd_mul(xk, xk), y2);
      const gsl_dd den = gsl_dd_mul_d(gsl_dd_add_d(c, k), k + 1.0);
      const gsl_dd ratio = gsl_dd_mul_d(gsl_dd_div(num, den), z);

      t = gsl_dd_mul(t, ratio);
      s = gsl_dd_add(s, t);
      sabs += fabs(t.hi);

      if (!gsl_finite(t.hi))
        return GSL_EOVRFLW;

      if (az < 1.0 && fabs(t.hi) < 1.0e-20 * fabs(s.hi) && k > kcross
          && GSL_MAX_DBL(fabs(ratio.hi), az) < 0.999)
        {
          const double rr = GSL_MAX_DBL(fabs(ratio.hi), az);
          tail = fabs(t.hi) * rr / (1.0 - rr);
          break;
        }
    }

  if (k >= HYP2F1_MAXITER)
    return GSL_EMAXITER;

  r->v = s.hi + s.lo;
  r->e = 1.0e-31 * HYP2F1_TERM_ERR(k) * sabs + tail + 0.5 * GSL_DBL_EPSILON * fabs(r->v);
  r->amp = 1.0;
  return GSL_SUCCESS;
}

static gsl_dd
hyp2f1_dd(const double x)
{
  gsl_dd r;
  r.hi = x;
  r.lo = 0.0;
  return r;
}

/*****************************************
 * lnGamma(x+e) - lnGamma(x)             *
 *****************************************/

/* B_2k / (2k (2k-1)), k = 1..9 */
static const double hyp2f1_stirling_c[9] = {
  8.33333333333333333333e-02,
 -2.77777777777777777778e-03,
  7.93650793650793650794e-04,
 -5.95238095238095238095e-04,
  8.41750841750841750842e-04,
 -1.91752691752691752692e-03,
  6.41025641025641025641e-03,
 -2.95506535947712418301e-02,
  1.79644372368830573165e-01
};

/*
 * lnGamma(x+e) - lnGamma(x) for x >= 1/2, x+e >= 1/2, to full relative
 * accuracy also for tiny e: shift both arguments past
 * HYP2F1_STIRLING_MIN with the exact terms log1p(e/(x+k)), then use the
 * difference of Stirling's series, in which every term is formed as
 * a difference of order e (log1p, expm1).
 */
static double
hyp2f1_dlngamma_pos(double x, const double e, double * err)
{
  double s = 0.0, se = 0.0, y, q, lq, t, yk, w, w2, pm1 = 0.0;
  int k;

  /*
   * sum_k log(1 + e/(x+k)) = log1p(P - 1), P = prod_k (1 + e/(x+k)),
   * with P - 1 accumulated as (P-1)(1+t) + t, which is free of
   * cancellation (all terms have the sign of e): one logarithm
   * instead of one per shift.
   */
  if (x < HYP2F1_STIRLING_MIN)
    {
      while (x < HYP2F1_STIRLING_MIN)
        {
          t = e / x;
          pm1 = pm1 * (1.0 + t) + t;
          x += 1.0;
        }
      t = gsl_sf_log1p(pm1);
      s -= t;
      se += fabs(t) * (HYP2F1_STIRLING_MIN + 1.0);
    }

  y = x;
  q = e / y;
  lq = gsl_sf_log1p(q);

  /* (y-1/2) log(1+e/y) + e (log(y+e) - 1) */
  t = (y - 0.5) * lq + e * (log(y + e) - 1.0);
  s += t;
  se += fabs((y - 0.5) * lq) + fabs(e * log(y + e)) + fabs(e);

  /*
   * sum_k c_k ((y+e)^(1-2k) - y^(1-2k)) = sum_k c_k y^(1-2k) w_k with
   * w_k = q^(1-2k) - 1, q = 1+e/y: w_1 = -e/(y+e), and
   * w_{k+1} = w_k + w2 (1 + w_k), w2 = q^(-2) - 1 = w_1 (2 + w_1);
   * all terms of the recurrence have one sign (that of -e).
   */
  w = -e / (y + e);
  w2 = w * (2.0 + w);
  yk = 1.0 / y;
  for (k = 1; k <= 9; ++k)
    {
      t = hyp2f1_stirling_c[k - 1] * yk * w;
      s += t;
      se += fabs(t) * k;
      if (fabs(t) < 0.1 * GSL_DBL_EPSILON * fabs(s))
        break;
      w += w2 * (1.0 + w);
      yk /= (y * y);
    }

  *err = 4.0 * GSL_DBL_EPSILON * se;
  return s;
}

/*
 * D(x,e) = log|Gamma(x+e)/Gamma(x)| and its sign, for |e| <= 1/2 and
 * x, x+e not non-positive integers.  Arguments below 1/2 are moved up
 * with the reflection formula, whose sine ratio
 *
 *   sin(pi (x+e)) / sin(pi x) = cos(pi e) + sin(pi e) cot(pi x)
 *
 * is formed as 1 - 2 sin^2(pi e/2) + sin(pi e) cot(pi x), again a
 * quantity of order e.
 */
/* sin(pi (hi+lo)) for |lo| <= ulp(hi), to full relative accuracy also near an integer */
static double
hyp2f1_sin_pi_dd(const double hi, const double lo)
{
  gsl_sf_result s, c;
  gsl_sf_sin_pi_e(hi, &s);
  if (lo == 0.0)
    return s.val;
  gsl_sf_cos_pi_e(hi, &c);
  return s.val + M_PI * lo * c.val;
}

/*
 * As hyp2f1_dlngamma(), with x = xhi + xlo and the endpoint x + e =
 * yhi + ylo given exactly (e itself is rounded).  Close to a pole of
 * Gamma(x) or Gamma(x+e) the result is set by the distance to the
 * pole, which x+e formed in double does not resolve (e.g. x = a,
 * x + e = c - b with c - b tiny).
 */
static double
hyp2f1_dlngamma_dd(const double xhi, const double xlo, const double e,
                   const double yhi, const double ylo, double * sgn, double * err)
{
  *sgn = 1.0;

  if (e == 0.0)
    {
      *err = 0.0;
      return 0.0;
    }
  else if (xhi >= 0.5 && yhi >= 0.5)
    {
      return hyp2f1_dlngamma_pos(xhi, e, err);
    }
  else if (xhi < 0.5 && yhi < 0.5)
    {
      gsl_sf_result cx, se2;
      double sx, cot, q, d, de;

      sx = hyp2f1_sin_pi_dd(xhi, xlo);
      gsl_sf_cos_pi_e(xhi, &cx);
      gsl_sf_sin_pi_e(0.5 * e, &se2);
      cot = cx.val / sx;
      q = -2.0 * se2.val * se2.val + 2.0 * se2.val * sqrt(1.0 - se2.val * se2.val) * cot;

      d = hyp2f1_dlngamma_pos(1.0 - xhi, -e, &de);
      *err = de + GSL_DBL_EPSILON * (4.0 * fabs(q) + 4.0 * fabs(se2.val * cot));

      if (fabs(1.0 + q) < 0.5)
        {
          /* x+e near a pole: the sine ratio from the exact endpoint */
          const double ratio = hyp2f1_sin_pi_dd(yhi, ylo) / sx;
          *err = de + 4.0 * GSL_DBL_EPSILON;
          if (ratio < 0.0)
            *sgn = -1.0;
          return -log(fabs(ratio)) - d;
        }

      if (q > -1.0)
        return -gsl_sf_log1p(q) - d;

      *sgn = -1.0;
      return -log(fabs(1.0 + q)) - d;
    }
  else
    {
      /* straddling 1/2: 0 < x, x+e < 3/2; use Gamma(y) = Gamma(y+1)/y */
      double de, t;
      const double d = hyp2f1_dlngamma_pos(xhi + 1.0, e, &de);
      if (yhi < 0.25 * xhi || xhi < 0.25 * yhi)
        t = log(yhi / xhi) + ylo / yhi - xlo / xhi;   /* x or x+e near the pole at 0 */
      else
        t = gsl_sf_log1p(e / xhi);
      *err = de + 2.0 * GSL_DBL_EPSILON * fabs(t);
      return d - t;
    }
}

static double
hyp2f1_dlngamma(const double x, const double e, double * sgn, double * err)
{
  const gsl_dd y = gsl_dd_two_sum(x, e);
  return hyp2f1_dlngamma_dd(x, 0.0, e, y.hi, y.lo, sgn, err);
}

/*
 * Re[lnGamma(z+e) - lnGamma(z)] for complex z = x + iy, real |e| <= 1/2,
 * y != 0, by the same method as hyp2f1_dlngamma(): shift until |z| is
 * past HYP2F1_STIRLING_MIN with the exact terms log|1 + e/(z+k)|, then
 * the difference of Stirling's series; reflection for x < 1/2.
 */
static double
hyp2f1_dlngamma_cre_pos(double x, const double y, const double e, double * err)
{
  double s = 0.0, se = 0.0, t;
  int k;

  while (x * x + y * y < HYP2F1_STIRLING_MIN * HYP2F1_STIRLING_MIN)
    {
      /* log|1 + e w|, w = 1/(x+iy) */
      const double den = x * x + y * y;
      const double wr = x / den, w2 = 1.0 / den;
      t = 0.5 * gsl_sf_log1p(2.0 * e * wr + e * e * w2);
      s -= t;
      se += fabs(t);
      x += 1.0;
    }

  {
    const double den = x * x + y * y;
    const double qr = e * x / den, qi = -e * y / den;           /* q = e/z */
    const double l1r = 0.5 * gsl_sf_log1p(2.0 * qr + qr * qr + qi * qi);
    const double l1i = atan2(qi, 1.0 + qr);                      /* log(1+q) */
    const double lzr = 0.5 * log((x + e) * (x + e) + y * y);     /* Re log(z+e) */
    double zir = x / den, zii = -y / den;                        /* 1/z */
    double pr, pi, z2r, z2i;

    /* Re[(z-1/2) log(1+q) + e (log(z+e) - 1)] */
    t = (x - 0.5) * l1r - y * l1i + e * (lzr - 1.0);
    s += t;
    se += fabs((x - 0.5) * l1r) + fabs(y * l1i) + fabs(e * lzr) + fabs(e);

    /* Re sum_k c_k z^(1-2k) expm1((1-2k) log(1+q)) */
    z2r = zir * zir - zii * zii;
    z2i = 2.0 * zir * zii;
    pr = zir;
    pi = zii;
    for (k = 1; k <= 9; ++k)
      {
        const double wr = (1.0 - 2.0 * k) * l1r, wi = (1.0 - 2.0 * k) * l1i;
        const double sh = sin(0.5 * wi);
        const double emr = gsl_sf_expm1(wr) * cos(wi) - 2.0 * sh * sh;
        const double emi = exp(wr) * sin(wi);
        double tr;

        t = hyp2f1_stirling_c[k - 1] * (pr * emr - pi * emi);
        s += t;
        se += fabs(t);

        tr = pr * z2r - pi * z2i;
        pi = pr * z2i + pi * z2r;
        pr = tr;
      }
  }

  *err = 4.0 * GSL_DBL_EPSILON * se;
  return s;
}

static double
hyp2f1_dlngamma_cre(const double x, const double y, const double e, double * err)
{
  if (e == 0.0)
    {
      *err = 0.0;
      return 0.0;
    }
  else if (x >= 0.5 && x + e >= 0.5)
    {
      return hyp2f1_dlngamma_cre_pos(x, y, e, err);
    }
  else
    {
      /*
       * Re D(z,e) = -log|sin(pi(z+e))/sin(pi z)| - Re D(1-z,-e), with
       * sin(pi(z+e))/sin(pi z) = 1 + q, q = cos(pi e) - 1 + sin(pi e) cot(pi z),
       * cot(pi z) = (sin(2 pi x) - i sinh(2 pi y)) / (cosh(2 pi y) - cos(2 pi x)).
       */
      gsl_sf_result s2x, c2x, se2;
      double den, cr, ci, sinpe, qr, qi, d, de;

      gsl_sf_sin_pi_e(2.0 * x, &s2x);
      gsl_sf_cos_pi_e(2.0 * x, &c2x);
      gsl_sf_sin_pi_e(0.5 * e, &se2);
      den = cosh(2.0 * M_PI * y) - c2x.val;
      cr = s2x.val / den;
      ci = -sinh(2.0 * M_PI * y) / den;
      sinpe = 2.0 * se2.val * sqrt(1.0 - se2.val * se2.val);
      qr = -2.0 * se2.val * se2.val + sinpe * cr;
      qi = sinpe * ci;

      d = hyp2f1_dlngamma_cre_pos(1.0 - x, -y, -e, &de);
      *err = de + 4.0 * GSL_DBL_EPSILON * (fabs(qr) + fabs(qi) + fabs(sinpe * cr));
      return -0.5 * gsl_sf_log1p(2.0 * qr + qr * qr + qi * qi) - d;
    }
}

/*****************************************
 * connection formula around z = 1       *
 *****************************************/

/*
 * F(a,b;c;z) for z = 1-u near 1 (0 < u <= 1/2) from the connection
 * formula (DLMF 15.8.4)
 *
 *   F = Gamma(c)Gamma(d)/(Gamma(c-a)Gamma(c-b)) F(a,b;1-d;u)
 *     + Gamma(c)Gamma(-d)/(Gamma(a)Gamma(b)) u^d F(c-a,c-b;1+d;u),
 *
 * d = c-a-b, rewritten so that it is evaluated uniformly in the distance
 * of d from an integer (after Michel and Stoitsov).  With d = m + e,
 * m >= 0 the nearest integer and |e| <= 1/2 (Euler's transformation is
 * applied first if d < -1/2), the first series splits into the terms
 * k < m, which stay regular, and a tail whose 1/e pole cancels that of
 * Gamma(-d) in the second term:
 *
 *   F = Gamma(c)Gamma(d)/(Gamma(c-a)Gamma(c-b)) [ S1 + C (pi e/sin(pi e)) sum_j R_j u^j ]
 *
 * with
 *
 *   S1 = sum_{k<m} (-1)^k (a)_k (b)_k Gamma(m-k+e)/(Gamma(m+e) k!) u^k,
 *   C  = (-1)^m (a)_m (b)_m u^m / (m! Gamma(m+e)),
 *
 * (for m = 0, S1 = 0 and C = 1/Gamma(d) is absorbed in the prefactor),
 * and R_j = (X_j - Y_j)/e the scaled difference of the two tails,
 *
 *   X_j = Gamma(a+m+j) Gamma(b+m+j) Gamma(1+m) / (Gamma(a+m) Gamma(b+m) Gamma(1-e+j) (m+j)!),
 *   Y_j = u^e Gamma(a+m+e+j) Gamma(b+m+e+j) Gamma(1+m) / (Gamma(a+m) Gamma(b+m) Gamma(1+m+e+j) j!).
 *
 * R_0 is formed from differences lnGamma(x+e) - lnGamma(x) computed to
 * full relative accuracy (expm1 of a sum of such differences), and
 * reduces to the classical psi-function form (A&S 15.3.10-12) at e = 0.
 * The R_j then follow from the exact recurrence
 *
 *   R_{j+1} = rx R_j + ((rx - ry)/e) Y_j,   Y_{j+1} = ry Y_j,
 *
 * rx = X_{j+1}/X_j, ry = Y_{j+1}/Y_j, where (rx - ry)/e is a polynomial
 * expression free of cancellation.
 *
 * Returns GSL_EDOM if the formula does not apply (m too large; the
 * caller then uses a series).
 */
static int
hyp2f1_conn_pair(const hyp2f1_pair * ab, const double c, const double u,
                 hyp2f1_scaled * r)
{
  const double lnu = log(u);
  double d, dhi, dlo, m_d, e, fe, K, Kl, Krel, S1 = 0.0, S1abs = 0.0, Cm = 1.0;
  double R, Re, Y, Yrel = 0.0, sum, sumabs, sumerr, Rerr;
  double jcross, j;
  int m, stat;

  /* d = c-a-b, with its rounding error kept in dlo */
  if (ab->conj)
    {
      const gsl_dd s1 = gsl_dd_two_sum(c, -2.0 * ab->x);
      dhi = s1.hi;
      dlo = s1.lo - 2.0 * ab->xlo;
    }
  else
    {
      const gsl_dd s1 = gsl_dd_two_sum(c, -ab->x);
      const gsl_dd s2 = gsl_dd_two_sum(s1.hi, -ab->y);
      dhi = s2.hi;
      dlo = s1.lo + s2.lo - ab->xlo - ab->ylo;
    }
  d = dhi + dlo;

  if (d < -0.5)
    {
      /* Euler: F(a,b;c;z) = u^d F(c-a,c-b;c;z), and c-(c-a)-(c-b) = -d */
      const hyp2f1_pair eu = hyp2f1_pair_euler(ab, c);
      stat = hyp2f1_conn_pair(&eu, c, u, r);
      r->l += d * lnu;
      r->e += fabs(r->v) * GSL_DBL_EPSILON * fabs(d * lnu);
      return stat;
    }

  if (d > 400.0)
    return GSL_EDOM;

  m_d = floor(d + 0.5);
  m = (int) m_d;
  e = (dhi - m_d) + dlo;

  /* prefactor Gamma(c) Gamma(m+e) / (Gamma(c-a) Gamma(c-b)) (Gamma(m+e) only for m >= 1) */
  if (ab->conj)
    {
      double num[2];
      gsl_sf_result lr, arg;
      num[0] = c;
      num[1] = d;
      stat = hyp2f1_gamma_ratio(num, (m >= 1) ? 2 : 1, num, 0, &K, &Kl, &Krel);
      if (stat != GSL_SUCCESS)
        return stat;
      /* |Gamma(c-a)|^2 */
      stat = gsl_sf_lngamma_complex_e(c - ab->x, -ab->y, &lr, &arg);
      if (stat != GSL_SUCCESS)
        return stat;
      Kl -= 2.0 * lr.val;
      Krel += 2.0 * lr.err + 8.0 * GSL_DBL_EPSILON * (1.0 + fabs(lr.val));
      {
        const gsl_dd cx = gsl_dd_two_sum(c, -ab->x);
        const double cxlo = cx.lo - ab->xlo;
        if (cxlo != 0.0)
          {
            gsl_sf_result pr, pi;
            gsl_sf_complex_psi_e(cx.hi, -ab->y, &pr, &pi);
            K *= 1.0 - 2.0 * pr.val * cxlo;
          }
      }
    }
  else
    {
      double num[2], den[2];
      const gsl_dd ca = gsl_dd_two_sum(c, -ab->x), cb = gsl_dd_two_sum(c, -ab->y);
      num[0] = c;
      num[1] = d;
      den[0] = ca.hi;
      den[1] = cb.hi;
      stat = hyp2f1_gamma_ratio(num, (m >= 1) ? 2 : 1, den, 2, &K, &Kl, &Krel);
      if (stat != GSL_SUCCESS)
        return stat;
      K *= hyp2f1_gamma_arg_corr(ca.hi, ca.lo - ab->xlo, -1.0)
           * hyp2f1_gamma_arg_corr(cb.hi, cb.lo - ab->ylo, -1.0);
    }

  if (m >= 1)
    K *= hyp2f1_gamma_arg_corr(d, (dhi - d) + dlo, 1.0);

  if (m >= 1)
    {
      /* S1 = sum_{k<m} t_k, t_0 = 1, t_{k+1}/t_k = -(a+k)(b+k)u/((k+1)(m-1-k+e)) */
      double t = 1.0;
      gsl_sf_result g;
      int k;

      S1 = 1.0;
      S1abs = 1.0;
      for (k = 0; k < m - 1; ++k)
        {
          t *= -hyp2f1_pair_prod(ab, (double) k) * u / ((k + 1.0) * (m - 1.0 - k + e));
          S1 += t;
          S1abs += fabs(t) * (k + 2.0);
        }

      gsl_sf_gamma_e(1.0 + e, &g);
      Cm = t * (-hyp2f1_pair_prod(ab, m - 1.0) * u / m) / g.val;
    }

  fe = (e == 0.0) ? 1.0 : (M_PI * e) / sin(M_PI * e);

  /* R_0 and Y_0 */
  if (e == 0.0)
    {
      /* R_0 = psi(1) - log u - psi(a+m) - psi(b+m) + psi(1+m) */
      double pab, pabe;
      gsl_sf_result p1;

      if (ab->conj)
        {
          gsl_sf_result pr, pi;
          gsl_sf_complex_psi_e(ab->x + m + ab->xlo, ab->y, &pr, &pi);
          pab = 2.0 * pr.val;
          pabe = 2.0 * pr.err;
        }
      else
        {
          /* a+m = c-b and b+m = c-a exactly here (e = 0) */
          gsl_sf_result pa, pb;
          gsl_dd xm = gsl_dd_two_sum(c, -ab->y), ym = gsl_dd_two_sum(c, -ab->x);
          xm.lo -= ab->ylo;
          ym.lo -= ab->xlo;
          gsl_sf_psi_e(xm.hi, &pa);
          gsl_sf_psi_e(ym.hi, &pb);
          if (xm.lo != 0.0)
            {
              gsl_sf_result p1x;
              gsl_sf_psi_1_e(xm.hi, &p1x);
              pa.val += p1x.val * xm.lo;
            }
          if (ym.lo != 0.0)
            {
              gsl_sf_result p1y;
              gsl_sf_psi_1_e(ym.hi, &p1y);
              pb.val += p1y.val * ym.lo;
            }
          pab = pa.val + pb.val;
          pabe = pa.err + pb.err + GSL_DBL_EPSILON * (fabs(pa.val) + fabs(pb.val));
        }
      gsl_sf_psi_int_e(1 + m, &p1);
      R = -M_EULER - lnu - pab + p1.val;
      Rerr = pabe + p1.err + 2.0 * GSL_DBL_EPSILON * (M_EULER + fabs(lnu) + fabs(pab) + fabs(p1.val));
      Y = 1.0;
    }
  else
    {
      double sab, Dab, eab, s1d, e1, eh, Dh, Dg, Gg, mag;
      const double D1 = hyp2f1_dlngamma(1.0 + m, e, &s1d, &e1);

      if (ab->conj)
        {
          Dab = 2.0 * hyp2f1_dlngamma_cre(ab->x + m + ab->xlo, ab->y, e, &eab);
          eab *= 2.0;
          sab = 1.0;
        }
      else
        {
          /* a+m and b+m, with the exact endpoints a+m+e = c-b, b+m+e = c-a */
          double sa, sb, ea, eb;
          gsl_dd xm = gsl_dd_two_sum(ab->x, m_d), ym = gsl_dd_two_sum(ab->y, m_d);
          gsl_dd xe = gsl_dd_two_sum(c, -ab->y), ye = gsl_dd_two_sum(c, -ab->x);
          xm.lo += ab->xlo;
          ym.lo += ab->ylo;
          xe.lo -= ab->ylo;
          ye.lo -= ab->xlo;
          Dab = hyp2f1_dlngamma_dd(xm.hi, xm.lo, e, xe.hi, xe.lo, &sa, &ea)
                + hyp2f1_dlngamma_dd(ym.hi, ym.lo, e, ye.hi, ye.lo, &sb, &eb);
          eab = ea + eb;
          sab = sa * sb;
        }

      /* Dh = -lnGamma(1-e) */
      Dh = -hyp2f1_dlngamma(1.0, -e, &s1d, &eh);
      Dg = e * lnu + Dab - D1;
      Gg = exp(Dg);
      mag = fabs(Dh) + fabs(e * lnu) + fabs(Dab) + fabs(D1);

      if (sab > 0.0)
        {
          R = Gg * gsl_sf_expm1(Dh - Dg) / e;
          Y = Gg;
        }
      else
        {
          R = (exp(Dh) + Gg) / e;
          Y = -Gg;
        }

      Rerr = (fabs(Gg) + fabs(exp(Dh))) * (eab + e1 + eh + 2.0 * GSL_DBL_EPSILON * mag) / fabs(e)
             + 2.0 * GSL_DBL_EPSILON * fabs(R);
      Yrel = eab + e1 + 2.0 * GSL_DBL_EPSILON * mag;
    }

  /* sum_j R_j u^j, carrying the propagated error of R_j alongside */
  sum = R;
  sumabs = fabs(R);
  sumerr = Rerr;
  Re = Rerr;
  /* stop only past the sign changes of a+m+j, b+m+j (cf. hyp2f1_series_pair) */
  jcross = GSL_MAX_DBL(-(hyp2f1_pair_minre(ab) + m), 0.0) + 2.0;

  for (j = 0.0; j < HYP2F1_MAXITER; j += 1.0)
    {
      const double AB = hyp2f1_pair_prod(ab, m + j), ApB = hyp2f1_pair_sum(ab, m + j);
      const double ABe = hyp2f1_pair_prod(ab, m + j + e);
      const double P = 1.0 + m + j, Q = 1.0 + j;
      const double rx = AB / ((Q - e) * P) * u;
      const double ry = ABe / ((P + e) * Q) * u;
      const double dN = AB * (P + Q) - ApB * P * Q + e * (ApB * P - P * Q) + e * e * P;
      const double dNabs = fabs(AB * (P + Q)) + fabs(ApB * P * Q)
                           + fabs(e) * (fabs(ApB * P) + P * Q) + e * e * P;
      const double dsc = u / ((Q - e) * P * (P + e) * Q);
      const double dr = dN * dsc;
      const double t1 = rx * R, t2 = dr * Y;

      /* rounding of rx, of dN (relative to dNabs, as dN may cancel), and
         the accumulated relative error Yrel of Y */
      Re = fabs(rx) * Re + 4.0 * GSL_DBL_EPSILON * fabs(t1)
           + (6.0 * GSL_DBL_EPSILON * dNabs + Yrel * fabs(dN)) * fabs(dsc * Y);
      Yrel += 4.0 * GSL_DBL_EPSILON;
      R = t1 + t2;
      Y *= ry;

      sum += R;
      sumabs += fabs(R);
      sumerr += Re + GSL_DBL_EPSILON * j * fabs(R);   /* rounding of u */

      if (!gsl_finite(R))
        return GSL_EOVRFLW;

      if (fabs(R) < 0.25 * GSL_DBL_EPSILON * fabs(sum) && fabs(rx) < 0.95
          && fabs(Y * dr) < 0.25 * GSL_DBL_EPSILON * fabs(sum)
          && j > jcross)
        {
          sumerr += (fabs(R) + fabs(Y)) * fabs(rx) / (1.0 - fabs(rx));
          break;
        }
    }

  if (j >= HYP2F1_MAXITER)
    return GSL_EMAXITER;

  {
    const double T = Cm * fe * sum;
    const double Terr = fabs(Cm * fe) * (sumerr + 2.0 * GSL_DBL_EPSILON * sumabs);
    const double v = S1 + T;

    r->v = K * v;
    r->e = fabs(K) * (Terr + 2.0 * GSL_DBL_EPSILON * (S1abs + fabs(T)))
           + Krel * fabs(r->v);
    r->l = Kl;
    r->amp = (S1abs + fabs(Cm * fe) * sumabs) / fabs(v);
  }

  return GSL_SUCCESS;
}

static int
hyp2f1_conn(const double a, const double alo, const double b, const double blo,
            const double c, const double u, hyp2f1_scaled * r)
{
  hyp2f1_pair ab;
  ab.conj = 0;
  ab.x = a;
  ab.y = b;
  ab.xlo = alo;
  ab.ylo = blo;
  return hyp2f1_conn_pair(&ab, c, u, r);
}

/*****************************************
 * terminating series                    *
 *****************************************/

/* (x)_n / (y)_n, the factors y+k nonzero; x has absolute uncertainty dx */
static double
hyp2f1_poch_ratio(const double x, const double y, const int n, const double dx,
                  double * rel)
{
  double p = 1.0, s = 0.0;
  int k;
  for (k = 0; k < n; ++k)
    {
      p *= (x + k) / (y + k);
      if (dx > 0.0 && x + k != 0.0)
        s += dx / fabs(x + k);
    }
  *rel = 2.0 * GSL_DBL_EPSILON * n + s;
  return p;
}

/* c' a non-positive integer reached before the series in which it is
   the lower parameter terminates after n+1 terms */
static int
hyp2f1_pole_before(const double cp, const int n)
{
  return (hyp2f1_isnpint(cp) && -cp < n);
}

/*
 * F(-n,b;c;z), a polynomial of degree n in z, for any real z.  Among the
 * equivalent terminating forms (DLMF 15.8.1, 15.8.6, 15.8.7)
 *
 *   F(-n,b;c;z)
 *     = (1-z)^n F(-n, c-b; c; z/(z-1))
 *     = (c-b)_n/(c)_n F(-n, b; b-c-n+1; 1-z)
 *     = (b)_n/(c)_n (-z)^n F(-n, 1-c-n; 1-b-n; 1/z),
 *
 * the sums are tried in increasing order of the modulus of their
 * variable, until one has a small relative error estimate.
 */
static int
hyp2f1_poly_candidate(const int k, const double n_d, const double b, const double c,
                      const double z, const int dd, hyp2f1_scaled * t)
{
  const int n = (int) n_d;
  const double a = -n_d;
  double rel, p;
  int stat;

  switch (k)
    {
    case 0:
      return dd ? hyp2f1_series_dd(hyp2f1_dd(a), hyp2f1_dd(b), hyp2f1_dd(c), hyp2f1_dd(z), t)
                : hyp2f1_series(a, b, c, z, 0.0, 0.0, 0.0, t);

    case 1:
      {
        const double omz = 1.0 - z;
        const double lo = n_d * log(fabs(omz));
        if (dd)
          {
            const gsl_dd w = gsl_dd_div(hyp2f1_dd(z), gsl_dd_two_sum(z, -1.0));
            stat = hyp2f1_series_dd(hyp2f1_dd(a), gsl_dd_two_sum(c, -b), hyp2f1_dd(c), w, t);
          }
        else
          stat = hyp2f1_series(a, c - b, c, z / (z - 1.0),
                               GSL_DBL_EPSILON * GSL_MAX_DBL(fabs(b), fabs(c)), 0.0,
                               2.0 * GSL_DBL_EPSILON, t);
        t->l += lo;
        if (omz < 0.0 && (n % 2))
          t->v = -t->v;
        t->e += fabs(t->v) * GSL_DBL_EPSILON * fabs(lo);
        return stat;
      }

    case 2:
      if (dd)
        {
          const gsl_dd cp = gsl_dd_add_d(gsl_dd_two_sum(b, -c), 1.0 - n_d);
          stat = hyp2f1_series_dd(hyp2f1_dd(a), hyp2f1_dd(b), cp, gsl_dd_two_sum(1.0, -z), t);
        }
      else
        stat = hyp2f1_series(a, b, b - c - n + 1.0, 1.0 - z, 0.0,
                             2.0 * GSL_DBL_EPSILON * GSL_MAX_DBL(GSL_MAX_DBL(fabs(b), fabs(c)), n_d),
                             GSL_DBL_EPSILON * fabs(z) / fabs(1.0 - z), t);
      p = hyp2f1_poch_ratio(c - b, c, n, GSL_DBL_EPSILON * GSL_MAX_DBL(fabs(b), fabs(c)), &rel);
      hyp2f1_scale(t, p, rel, 0.0);
      return stat;

    default:
      {
        const double lz = n_d * log(fabs(z));
        if (dd)
          {
            const gsl_dd bp = gsl_dd_add_d(hyp2f1_dd(-c), 1.0 - n_d);
            const gsl_dd cp = gsl_dd_add_d(hyp2f1_dd(-b), 1.0 - n_d);
            const gsl_dd iz = gsl_dd_div(hyp2f1_dd(1.0), hyp2f1_dd(z));
            stat = hyp2f1_series_dd(hyp2f1_dd(a), bp, cp, iz, t);
          }
        else
          {
            const double dp = 2.0 * GSL_DBL_EPSILON * GSL_MAX_DBL(GSL_MAX_DBL(fabs(b), fabs(c)), n_d);
            stat = hyp2f1_series(a, 1.0 - c - n, 1.0 - b - n, 1.0 / z, dp, dp,
                                 GSL_DBL_EPSILON, t);
          }
        p = hyp2f1_poch_ratio(b, c, n, 0.0, &rel);
        if (z > 0.0 && (n % 2))
          p = -p;
        hyp2f1_scale(t, p, rel, lz);
        t->e += fabs(t->v) * GSL_DBL_EPSILON * fabs(lz);
        return stat;
      }
    }
}

static int
hyp2f1_poly(const double n_d, const double b, const double c, const double z,
            hyp2f1_scaled * r)
{
  const int n = (int) n_d;
  double vars[4], amps[4];
  int order[4], i, j, have = 0, kbest = -1;
  hyp2f1_scaled best, t;

  best.v = 0.0;
  best.e = GSL_POSINF;
  best.l = 0.0;
  best.amp = GSL_POSINF;

  /* candidate variables (+Inf marks an invalid form) */
  vars[0] = fabs(z);
  vars[1] = (z != 1.0) ? fabs(z / (z - 1.0)) : GSL_POSINF;
  vars[2] = hyp2f1_pole_before(b - c - n + 1.0, n) ? GSL_POSINF : fabs(1.0 - z);
  vars[3] = (z != 0.0 && !hyp2f1_pole_before(1.0 - b - n, n)) ? fabs(1.0 / z) : GSL_POSINF;

  for (i = 0; i < 4; ++i)
    {
      order[i] = i;
      amps[i] = GSL_POSINF;
    }
  for (i = 0; i < 4; ++i)
    for (j = i + 1; j < 4; ++j)
      if (vars[order[j]] < vars[order[i]])
        {
          int tmp = order[i];
          order[i] = order[j];
          order[j] = tmp;
        }

  for (i = 0; i < 4; ++i)
    {
      const int k = order[i];
      int stat;

      if (!gsl_finite(vars[k]))
        continue;

      stat = hyp2f1_poly_candidate(k, n_d, b, c, z, 0, &t);
      if (stat == GSL_SUCCESS && gsl_finite(t.v))
        {
          amps[k] = t.amp;
          if (kbest < 0 || t.amp < amps[kbest])
            kbest = k;
          if (!have || hyp2f1_relerr(&t) < hyp2f1_relerr(&best))
            {
              best = t;
              have = 1;
            }
        }

      if (have && best.amp <= HYP2F1_GOOD_AMP)
        break;
    }

  /* strong cancellation in every form: re-sum the best one in double-double */
  if (kbest >= 0 && best.amp > HYP2F1_DD_AMP
      && hyp2f1_poly_candidate(kbest, n_d, b, c, z, 1, &t) == GSL_SUCCESS
      && gsl_finite(t.v) && hyp2f1_relerr(&t) < hyp2f1_relerr(&best))
    best = t;

  *r = best;
  return have ? GSL_SUCCESS : GSL_EMAXITER;
}

/*****************************************
 * general case                          *
 *****************************************/

enum
{
  HYP2F1_DIRECT,      /* F(a,b;c;z) */
  HYP2F1_EULER,       /* (1-z)^(c-a-b) F(c-a,c-b;c;z) */
  HYP2F1_PFAFF_A,     /* (1-z)^(-a) F(a,c-b;c;z/(z-1)) */
  HYP2F1_PFAFF_B,     /* (1-z)^(-b) F(c-a,b;c;z/(z-1)) */
  HYP2F1_CONN,        /* connection formula in 1-z */
  HYP2F1_PFAFF_A_CONN,/* Pfaff (a), then connection formula in 1/(1-z) */
  HYP2F1_PFAFF_B_CONN /* Pfaff (b), then connection formula in 1/(1-z) */
};

/*
 * F(a,b;c;z) by the series or the terminating series, whichever
 * applies; a, b are transformed parameters with rounding residuals
 * alo, blo (exact values a+alo, b+blo), and z has relative rounding
 * error dzr.
 */
static int
hyp2f1_series_or_poly(const double a, const double alo, const double b, const double blo,
                      const double c, const double z, const double dzr,
                      hyp2f1_scaled * r)
{
  const int ta = hyp2f1_isnpint(a) && alo == 0.0, tb = hyp2f1_isnpint(b) && blo == 0.0;
  hyp2f1_pair ab;

  if (ta || tb)
    {
      const double n = ta && (!tb || a >= b) ? -a : -b;
      const double bb = (n == -a) ? b + blo : a + alo;
      return hyp2f1_poly(n, bb, c, z, r);
    }

  ab.conj = 0;
  ab.x = a;
  ab.y = b;
  ab.xlo = alo;
  ab.ylo = blo;
  return hyp2f1_series_pair(&ab, c, z, 0.0, 0.0, dzr, r);
}

/* c-a-b correctly rounded from c-a = ca (with its residual): c-a-b can
   be much smaller than the parameters, and its absolute error enters
   the Euler factor (1-z)^(c-a-b) as a relative error */
static double
hyp2f1_cab(const gsl_dd ca, const double b)
{
  const gsl_dd s = gsl_dd_two_sum(ca.hi, -b);
  return s.hi + (s.lo + ca.lo);
}

static int
hyp2f1_candidate(const int kind, const double a, const double b, const double c,
                 const double z, hyp2f1_scaled * r)
{
  /* c-a, c-b with their rounding residuals */
  const gsl_dd ca = gsl_dd_two_sum(c, -a), cb = gsl_dd_two_sum(c, -b);
  int stat;

  switch (kind)
    {
    case HYP2F1_DIRECT:
      return hyp2f1_series(a, b, c, z, 0.0, 0.0, 0.0, r);

    case HYP2F1_EULER:
      {
        const double d = hyp2f1_cab(ca, b);
        const double lp = d * gsl_sf_log1p(-z);
        stat = hyp2f1_series_or_poly(ca.hi, ca.lo, cb.hi, cb.lo, c, z, 0.0, r);
        r->l += lp;
        r->e += fabs(r->v) * GSL_DBL_EPSILON * (fabs(lp) + 1.0);
        return stat;
      }

    case HYP2F1_PFAFF_A:
    case HYP2F1_PFAFF_B:
      {
        const double p = (kind == HYP2F1_PFAFF_A) ? a : b;
        const double lp = -p * gsl_sf_log1p(-z);
        if (kind == HYP2F1_PFAFF_A)
          stat = hyp2f1_series_or_poly(a, 0.0, cb.hi, cb.lo, c, z / (z - 1.0),
                                       2.0 * GSL_DBL_EPSILON, r);
        else
          stat = hyp2f1_series_or_poly(ca.hi, ca.lo, b, 0.0, c, z / (z - 1.0),
                                       2.0 * GSL_DBL_EPSILON, r);
        r->l += lp;
        r->e += fabs(r->v) * GSL_DBL_EPSILON * (fabs(lp) + 1.0);
        return stat;
      }

    case HYP2F1_CONN:
      return hyp2f1_conn(a, 0.0, b, 0.0, c, 1.0 - z, r);

    default:
      {
        /* z < 0: w = z/(z-1) in (0,1), 1-w = 1/(1-z) */
        const int pa = (kind == HYP2F1_PFAFF_A_CONN);
        const double u = 1.0 / (1.0 - z);
        const double lp = -(pa ? a : b) * gsl_sf_log1p(-z);
        const double a2 = pa ? a : ca.hi, a2lo = pa ? 0.0 : ca.lo;
        const double b2 = pa ? cb.hi : b, b2lo = pa ? cb.lo : 0.0;

        if ((hyp2f1_isnpint(a2) && a2lo == 0.0) || (hyp2f1_isnpint(b2) && b2lo == 0.0))
          stat = hyp2f1_series_or_poly(a2, a2lo, b2, b2lo, c, z / (z - 1.0),
                                       2.0 * GSL_DBL_EPSILON, r);
        else
          stat = hyp2f1_conn(a2, a2lo, b2, b2lo, c, u, r);
        r->l += lp;
        r->e += fabs(r->v) * GSL_DBL_EPSILON * (fabs(lp) + 1.0);
        return stat;
      }
    }
}

/* series-type candidates re-summed in double-double (exact parameter combinations) */
static int
hyp2f1_candidate_dd(const int kind, const double a, const double b, const double c,
                    const double z, hyp2f1_scaled * r)
{
  const gsl_dd ca = gsl_dd_two_sum(c, -a), cb = gsl_dd_two_sum(c, -b);
  int stat;

  switch (kind)
    {
    case HYP2F1_DIRECT:
      return hyp2f1_series_dd(hyp2f1_dd(a), hyp2f1_dd(b), hyp2f1_dd(c), hyp2f1_dd(z), r);

    case HYP2F1_EULER:
      {
        const double lp = hyp2f1_cab(ca, b) * gsl_sf_log1p(-z);
        stat = hyp2f1_series_dd(ca, cb, hyp2f1_dd(c), hyp2f1_dd(z), r);
        r->l += lp;
        r->e += fabs(r->v) * GSL_DBL_EPSILON * (fabs(lp) + 1.0);
        return stat;
      }

    case HYP2F1_PFAFF_A:
    case HYP2F1_PFAFF_B:
      {
        const double p = (kind == HYP2F1_PFAFF_A) ? a : b;
        const double lp = -p * gsl_sf_log1p(-z);
        const gsl_dd w = gsl_dd_div(hyp2f1_dd(z), gsl_dd_two_sum(z, -1.0));
        if (kind == HYP2F1_PFAFF_A)
          stat = hyp2f1_series_dd(hyp2f1_dd(a), cb, hyp2f1_dd(c), w, r);
        else
          stat = hyp2f1_series_dd(ca, hyp2f1_dd(b), hyp2f1_dd(c), w, r);
        r->l += lp;
        r->e += fabs(r->v) * GSL_DBL_EPSILON * (fabs(lp) + 1.0);
        return stat;
      }

    default:
      return GSL_EDOM;
    }
}

/*
 * F(a,b;c;z) for z < 1, z != 0, with a, b, c-a, c-b not non-positive
 * integers and c not a non-positive integer.  The candidate
 * transformations for the region of z are tried in turn (the first is
 * the one whose variable is smallest); the result with the smallest
 * relative error estimate is kept, stopping as soon as it is small.
 */
static int
hyp2f1_general(const double a, const double b, const double c, const double z,
               hyp2f1_scaled * r)
{
  int kinds[6], nk = 0, i, have = 0, last = GSL_SUCCESS, series = 0;
  /* beyond z = 0.85 the series fallbacks are slow and, being long, no
     more accurate than the connection formula: stop on cancellation alone */
  const double good_rel = (z > 0.85) ? GSL_POSINF : HYP2F1_GOOD_REL;
  double amps[6];
  hyp2f1_scaled best, t;

  /*
   * The connection formulas replace c by 1 -+ (c-a-b) or 1 -+ (b-a) as
   * the lower parameter of the series, so they are kept as fallbacks
   * also where their variable is not the smallest one: they avoid the
   * cancellation of the c-series for large negative c.
   */
  if (z > 0.0 && z <= 0.5)
    {
      kinds[nk++] = HYP2F1_DIRECT;
      kinds[nk++] = HYP2F1_EULER;
      if (z >= 0.05)
        kinds[nk++] = HYP2F1_CONN;
    }
  else if (z > 0.5)
    {
      /*
       * Up to z = 3/4 the series still converges fast, and without
       * cancellation it is more accurate than the connection formula,
       * whose Gamma function prefactors cost a few ulp.
       */
      if (z <= 0.75)
        kinds[nk++] = HYP2F1_DIRECT;
      kinds[nk++] = HYP2F1_CONN;
      if (z > 0.75)
        kinds[nk++] = HYP2F1_DIRECT;
      kinds[nk++] = HYP2F1_EULER;
    }
  else if (z >= -1.0)
    {
      kinds[nk++] = HYP2F1_PFAFF_A;
      kinds[nk++] = HYP2F1_PFAFF_B;
      kinds[nk++] = HYP2F1_DIRECT;
      if (z >= -0.5)
        kinds[nk++] = HYP2F1_EULER;
      if (z <= -0.05)
        {
          kinds[nk++] = HYP2F1_PFAFF_A_CONN;
          kinds[nk++] = HYP2F1_PFAFF_B_CONN;
        }
    }
  else
    {
      /*
       * For -3 <= z < -1 the Pfaff series in w = z/(z-1) <= 3/4 still
       * converges well, and as a single series it avoids the Gamma
       * function prefactors and the cancellation between the two terms
       * of the 1/(1-z) connection formula.
       */
      if (z >= -3.0)
        {
          kinds[nk++] = HYP2F1_PFAFF_A;
          kinds[nk++] = HYP2F1_PFAFF_B;
        }
      kinds[nk++] = HYP2F1_PFAFF_A_CONN;
      kinds[nk++] = HYP2F1_PFAFF_B_CONN;
    }

  for (i = 0; i < nk; ++i)
    amps[i] = GSL_POSINF;   /* also those not reached */

  best.v = 0.0;
  best.e = GSL_POSINF;
  best.l = 0.0;
  best.amp = GSL_POSINF;

  for (i = 0; i < nk; ++i)
    {
      int stat = hyp2f1_candidate(kinds[i], a, b, c, z, &t);

      amps[i] = (stat == GSL_SUCCESS && gsl_finite(t.v)) ? t.amp : GSL_POSINF;
      if (stat == GSL_SUCCESS && gsl_finite(t.v) && gsl_finite(t.e)
          && (!have || hyp2f1_relerr(&t) < hyp2f1_relerr(&best)))
        {
          best = t;
          have = 1;
          series = (kinds[i] <= HYP2F1_PFAFF_B);
        }
      else if (stat != GSL_SUCCESS)
        {
          last = stat;
        }

      /* the error estimates of the connection formulas are pessimistic
         (Gamma function prefactors): for them cancellation decides */
      if (have && best.amp <= HYP2F1_GOOD_AMP
          && (!series || hyp2f1_relerr(&best) <= good_rel))
        break;
    }

  /*
   * Strong cancellation or a large error estimate in every candidate:
   * re-sum the series-type candidate with the least cancellation in
   * double-double.
   */
  {
    const int by_amp = (i == nk && (!have || best.amp > HYP2F1_DD_AMP));
    const int by_rel = (have && hyp2f1_relerr(&best) > HYP2F1_DD_REL);
    const double wmax = by_amp ? GSL_POSINF : HYP2F1_DD_VAR;
    const double zw = fabs(z / (z - 1.0));

    if (by_amp || by_rel)
    {
      int kbest = -1;
      for (i = 0; i < nk; ++i)
        if (kinds[i] <= HYP2F1_PFAFF_B && amps[i] < GSL_POSINF
            && ((kinds[i] <= HYP2F1_EULER) ? fabs(z) : zw) <= wmax
            && (kbest < 0 || amps[i] < amps[kbest]))
          kbest = i;

      if (kbest >= 0 && hyp2f1_candidate_dd(kinds[kbest], a, b, c, z, &t) == GSL_SUCCESS
          && gsl_finite(t.v) && (!have || hyp2f1_relerr(&t) < hyp2f1_relerr(&best)))
        {
          best = t;
          have = 1;
        }
    }
  }

  *r = best;
  return have ? GSL_SUCCESS : (last != GSL_SUCCESS ? last : GSL_EMAXITER);
}

/*
 * F(a,b;c;x) as a scaled value.  Handles the special and degenerate
 * parameter cases, then dispatches to hyp2f1_general().
 */
static int
hyp2f1_scaled_eval(const double a, const double b, const double c,
                   const double x, hyp2f1_scaled * r)
{
  const int a_np = hyp2f1_isnpint(a), b_np = hyp2f1_isnpint(b);

  r->v = 0.0;
  r->e = 0.0;
  r->l = 0.0;
  r->amp = 1.0;

  if (gsl_isnan(a) || gsl_isnan(b) || gsl_isnan(c) || gsl_isnan(x)
      || !gsl_finite(a) || !gsl_finite(b) || !gsl_finite(c) || !gsl_finite(x))
    return GSL_EDOM;

  if (a == 0.0 || b == 0.0)
    {
      r->v = 1.0;
      return GSL_SUCCESS;
    }

  /* terminating series: a polynomial in x, defined for all real x */
  if (a_np || b_np)
    {
      const double n = (a_np && (!b_np || a >= b)) ? -a : -b;
      const double bb = (n == -a) ? b : a;

      if (n > 1.0e7)
        return GSL_EDOM;

      if (hyp2f1_isnpint(c) && -c < n)
        return GSL_EDOM;        /* c pole reached before termination */

      return hyp2f1_poly(n, bb, c, x, r);
    }

  /* a pole of the series (undefined also at x = 0, bug #66850) */
  if (hyp2f1_isnpint(c))
    return GSL_EDOM;

  if (x == 0.0)
    {
      r->v = 1.0;
      return GSL_SUCCESS;
    }

  /* c-a or c-b a non-positive integer: (1-x)^(c-a-b) times a polynomial */
  if (hyp2f1_isnpint(c - a) || hyp2f1_isnpint(c - b))
    {
      const double ca = c - a, cb = c - b, d = c - a - b;
      const double n = (hyp2f1_isnpint(ca) && (!hyp2f1_isnpint(cb) || ca >= cb)) ? -ca : -cb;
      const double bb = (n == -ca) ? cb : ca;
      int stat;

      if (n > 1.0e7)
        return GSL_EDOM;

      if (x == 1.0)
        {
          if (d > 0.0)
            {
              r->v = 0.0;
              return GSL_SUCCESS;
            }
          else if (d == 0.0)
            return hyp2f1_poly(n, bb, c, x, r);
          else
            return GSL_EDOM;
        }

      if (x > 1.0 && d != floor(d))
        return GSL_EDOM;        /* complex */

      stat = hyp2f1_poly(n, bb, c, x, r);
      r->l += d * log(fabs(1.0 - x));
      r->e += fabs(r->v) * GSL_DBL_EPSILON * fabs(d * log(fabs(1.0 - x)));
      if (x > 1.0 && fmod(fabs(d), 2.0) == 1.0)
        r->v = -r->v;
      return stat;
    }

  if (x > 1.0)
    return GSL_EDOM;            /* branch cut */

  if (x == 1.0)
    {
      /* Gauss: F(a,b;c;1) = Gamma(c)Gamma(c-a-b)/(Gamma(c-a)Gamma(c-b)), c-a-b > 0 */
      const gsl_dd s1 = gsl_dd_two_sum(c, -a);
      const gsl_dd s2 = gsl_dd_two_sum(s1.hi, -b);
      const double d = s2.hi + (s1.lo + s2.lo);
      const gsl_dd ca = gsl_dd_two_sum(c, -a), cb = gsl_dd_two_sum(c, -b);
      double num[2], den[2], rel;
      int stat;

      if (d <= 0.0)
        return GSL_EDOM;        /* divergent */

      num[0] = c;
      num[1] = d;
      den[0] = ca.hi;
      den[1] = cb.hi;
      stat = hyp2f1_gamma_ratio(num, 2, den, 2, &r->v, &r->l, &rel);
      r->v *= hyp2f1_gamma_arg_corr(ca.hi, ca.lo, -1.0)
              * hyp2f1_gamma_arg_corr(cb.hi, cb.lo, -1.0)
              * hyp2f1_gamma_arg_corr(d, (s2.hi - d) + (s1.lo + s2.lo), 1.0);
      r->e = rel * fabs(r->v);
      return stat;
    }

  return hyp2f1_general(a, b, c, x, r);
}

/* assemble a gsl_sf_result from a scaled value */
static int
hyp2f1_result(const int stat, const hyp2f1_scaled * r, gsl_sf_result * result)
{
  int stat_e;

  if (stat == GSL_EDOM)
    {
      DOMAIN_ERROR(result);
    }

  if (r->l == 0.0)
    {
      /* no scaling */
      result->val = r->v;
      result->err = r->e + GSL_DBL_EPSILON * fabs(r->v);
      stat_e = (fabs(r->v) <= GSL_DBL_MAX) ? GSL_SUCCESS : GSL_EOVRFLW;
    }
  else
    stat_e = gsl_sf_exp_mult_err_e(r->l, GSL_DBL_EPSILON * fabs(r->l),
                                   r->v, r->e, result);

  if (stat != GSL_SUCCESS)
    {
      if (stat == GSL_EMAXITER || stat == GSL_EOVRFLW)
        GSL_ERROR("error", stat);
      return stat;
    }

  /*
   * No significant digit left (cancellation beyond what the available
   * representations and double-double summation can recover, which can
   * happen for very large parameters): report the loss rather than a
   * meaningless value with GSL_SUCCESS.
   */
  if (stat_e == GSL_SUCCESS && result->err > 0.0 && !(result->err < fabs(result->val)))
    GSL_ERROR("error", GSL_ELOSS);

  return stat_e;
}

/*****************************************
 * conjugate parameters a = x+iy, b = x-iy *
 *****************************************/

/*
 * F(a, c-b; a-b+1; u) term of the 1/(1-z) connection formula (DLMF
 * 15.8.2 with Pfaff's transformation), for the conjugate pair: the two
 * terms of the formula are complex conjugates, so
 *
 *   F(a,b;c;z) = 2 Re[ Gamma(c)Gamma(b-a)/(Gamma(b)Gamma(c-a))
 *                      (1-z)^(-a) F(a, c-b; a-b+1; 1/(1-z)) ],
 *
 * with a-b = 2iy.  The series has complex parameters and is summed in
 * complex arithmetic.  Used for z < -1/2; the two terms cancel as
 * y -> 0 (the error estimate reflects this).
 */
static int
hyp2f1_conj_cplx(const double x, const double y, const double c, const double z,
                 hyp2f1_scaled * r)
{
  const double u = 1.0 / (1.0 - z);
  const double L = gsl_sf_log1p(-z);
  const double kcross = GSL_MAX_DBL(GSL_MAX_DBL(-x, -(c - x)), 0.0) + 2.0;
  double sr = 1.0, si = 0.0, tr = 1.0, ti = 0.0, sabs = 1.0, k;
  double lnB, phi, Berr, th, cth, sth;
  int stat;

  /* B = Gamma(c) Gamma(-2iy) / (Gamma(x-iy) Gamma(c-x-iy)) */
  {
    gsl_sf_result lc, l1, a1, l2, a2, l3, a3;
    double sc;
    stat = gsl_sf_lngamma_sgn_e(c, &lc, &sc);
    stat += gsl_sf_lngamma_complex_e(0.0, -2.0 * y, &l1, &a1);
    stat += gsl_sf_lngamma_complex_e(x, -y, &l2, &a2);
    stat += gsl_sf_lngamma_complex_e(c - x, -y, &l3, &a3);
    if (stat != GSL_SUCCESS)
      return GSL_EDOM;
    lnB = lc.val + l1.val - l2.val - l3.val;
    phi = a1.val - a2.val - a3.val + ((sc < 0.0) ? M_PI : 0.0);
    Berr = lc.err + l1.err + l2.err + l3.err + a1.err + a2.err + a3.err
           + 2.0 * GSL_DBL_EPSILON * (fabs(lnB) + fabs(phi));
  }

  /* S = F(a, c-b; 1+2iy; u), a = x+iy, c-b = c-x+iy */
  for (k = 0.0; k < HYP2F1_MAXITER; k += 1.0)
    {
      /* ratio = (a+k)(c-b+k) / ((1+2iy+k)(k+1)) u */
      const double pr = (x + k) * (c - x + k) - y * y;
      const double pi = y * (c - x + k) + y * (x + k);
      const double qr = 1.0 + k, qi = 2.0 * y;
      const double den = (qr * qr + qi * qi) * (k + 1.0);
      const double rr = (pr * qr + pi * qi) / den * u;
      const double ri = (pi * qr - pr * qi) / den * u;
      const double t2 = tr * rr - ti * ri;

      ti = tr * ri + ti * rr;
      tr = t2;
      sr += tr;
      si += ti;
      sabs += gsl_hypot(tr, ti);

      if (!gsl_finite(sabs))
        return GSL_EOVRFLW;

      if (gsl_hypot(tr, ti) < 0.25 * GSL_DBL_EPSILON * gsl_hypot(sr, si) && k > kcross
          && gsl_hypot(rr, ri) < 0.999)
        {
          const double rm = gsl_hypot(rr, ri);
          sabs += gsl_hypot(tr, ti) * rm / (1.0 - rm);
          break;
        }
    }

  if (k >= HYP2F1_MAXITER)
    return GSL_EMAXITER;

  th = phi - y * L;
  cth = cos(th);
  sth = sin(th);

  r->v = 2.0 * (cth * sr - sth * si);
  r->l = lnB - x * L;
  r->e = 2.0 * (GSL_DBL_EPSILON * HYP2F1_TERM_ERR(k) * sabs
                + (Berr + GSL_DBL_EPSILON * fabs(y * L)) * gsl_hypot(sr, si))
         + GSL_DBL_EPSILON * fabs(x * L) * fabs(r->v);
  r->amp = 2.0 * sabs / fabs(r->v);
  return GSL_SUCCESS;
}

enum
{
  HYP2F1_C_DIRECT,
  HYP2F1_C_EULER,
  HYP2F1_C_CONN,
  HYP2F1_C_CPLX
};

static int
hyp2f1_conj_candidate(const int kind, const hyp2f1_pair * ab, const double c,
                      const double z, hyp2f1_scaled * r)
{
  int stat;

  switch (kind)
    {
    case HYP2F1_C_DIRECT:
      return hyp2f1_series_pair(ab, c, z, 0.0, 0.0, 0.0, r);

    case HYP2F1_C_EULER:
      {
        const hyp2f1_pair eu = hyp2f1_pair_euler(ab, c);
        const double lp = (c - 2.0 * ab->x) * gsl_sf_log1p(-z);
        stat = hyp2f1_series_pair(&eu, c, z, 0.0, 0.0, 0.0, r);
        r->l += lp;
        r->e += fabs(r->v) * GSL_DBL_EPSILON * (fabs(lp) + 1.0);
        return stat;
      }

    case HYP2F1_C_CONN:
      return hyp2f1_conn_pair(ab, c, 1.0 - z, r);

    default:
      return hyp2f1_conj_cplx(ab->x, ab->y, c, z, r);
    }
}

/* F(x+iy, x-iy; c; z) as a scaled value */
static int
hyp2f1_conj_scaled_eval(const double x, const double y, const double c,
                        const double z, hyp2f1_scaled * r)
{
  hyp2f1_pair ab;
  hyp2f1_scaled best, t;
  int kinds[4], nk = 0, i, have = 0, last = GSL_SUCCESS, series = 0;
  /* beyond z = 0.85 the series fallbacks are slow and, being long, no
     more accurate than the connection formula: stop on cancellation alone */
  const double good_rel = (z > 0.85) ? GSL_POSINF : HYP2F1_GOOD_REL;
  double amps[4];

  r->v = 0.0;
  r->e = 0.0;
  r->l = 0.0;
  r->amp = 1.0;

  if (!gsl_finite(x) || !gsl_finite(y) || !gsl_finite(c) || !gsl_finite(z))
    return GSL_EDOM;

  if (y == 0.0)
    return hyp2f1_scaled_eval(x, x, c, z, r);

  /* (a)_k (b)_k = |(a)_k|^2 never vanishes: c a pole is never cancelled */
  if (hyp2f1_isnpint(c) || z > 1.0)
    return GSL_EDOM;

  if (z == 0.0)
    {
      r->v = 1.0;
      return GSL_SUCCESS;
    }

  ab.conj = 1;
  ab.x = x;
  ab.y = y;
  ab.xlo = 0.0;
  ab.ylo = 0.0;

  if (z == 1.0)
    {
      /* Gauss: Gamma(c) Gamma(c-2x) / |Gamma(c-a)|^2, c-2x > 0 */
      const double d = c - 2.0 * x;
      double num[2], rel;
      gsl_sf_result lr, arg;
      int stat;

      if (d <= 0.0)
        return GSL_EDOM;
      num[0] = c;
      num[1] = d;
      stat = hyp2f1_gamma_ratio(num, 2, num, 0, &r->v, &r->l, &rel);
      stat += gsl_sf_lngamma_complex_e(c - x, -y, &lr, &arg);
      if (stat != GSL_SUCCESS)
        return GSL_EDOM;
      r->l -= 2.0 * lr.val;
      r->e = (rel + 2.0 * lr.err + 8.0 * GSL_DBL_EPSILON * (1.0 + fabs(lr.val))) * fabs(r->v);
      return GSL_SUCCESS;
    }

  if (z > 0.5)
    {
      /* as in hyp2f1_general(): the series first while it converges fast */
      if (z <= 0.85)
        kinds[nk++] = HYP2F1_C_DIRECT;
      kinds[nk++] = HYP2F1_C_CONN;
      if (z > 0.85)
        kinds[nk++] = HYP2F1_C_DIRECT;
      kinds[nk++] = HYP2F1_C_EULER;
    }
  else if (z >= -0.5)
    {
      kinds[nk++] = HYP2F1_C_DIRECT;
      kinds[nk++] = HYP2F1_C_EULER;
      if (z >= 0.05)
        kinds[nk++] = HYP2F1_C_CONN;
      if (z <= -0.05)
        kinds[nk++] = HYP2F1_C_CPLX;
    }
  else
    {
      kinds[nk++] = HYP2F1_C_CPLX;
      if (z >= -1.0)
        {
          kinds[nk++] = HYP2F1_C_DIRECT;
          kinds[nk++] = HYP2F1_C_EULER;
        }
    }

  for (i = 0; i < nk; ++i)
    amps[i] = GSL_POSINF;   /* also those not reached */

  best.v = 0.0;
  best.e = GSL_POSINF;
  best.l = 0.0;
  best.amp = GSL_POSINF;

  for (i = 0; i < nk; ++i)
    {
      int stat = hyp2f1_conj_candidate(kinds[i], &ab, c, z, &t);

      amps[i] = (stat == GSL_SUCCESS && gsl_finite(t.v)) ? t.amp : GSL_POSINF;

      if (stat == GSL_SUCCESS && gsl_finite(t.v) && gsl_finite(t.e)
          && (!have || hyp2f1_relerr(&t) < hyp2f1_relerr(&best)))
        {
          best = t;
          have = 1;
          series = (kinds[i] == HYP2F1_C_DIRECT || kinds[i] == HYP2F1_C_EULER);
        }
      else if (stat != GSL_SUCCESS)
        {
          last = stat;
        }

      /* as in hyp2f1_general() */
      if (have && best.amp <= HYP2F1_GOOD_AMP
          && (!series || hyp2f1_relerr(&best) <= good_rel))
        break;
    }

  /* strong cancellation or a large error estimate: re-sum the better of
     the two series in double-double */
  if ((i == nk && (!have || best.amp > HYP2F1_DD_AMP))
      || (have && hyp2f1_relerr(&best) > HYP2F1_DD_REL && fabs(z) <= HYP2F1_DD_VAR))
    {
      int kbest = -1;
      for (i = 0; i < nk; ++i)
        if ((kinds[i] == HYP2F1_C_DIRECT || kinds[i] == HYP2F1_C_EULER) && amps[i] < GSL_POSINF
            && (kbest < 0 || amps[i] < amps[kbest]))
          kbest = i;

      if (kbest >= 0)
        {
          int stat;
          if (kinds[kbest] == HYP2F1_C_DIRECT)
            stat = hyp2f1_series_conj_dd(hyp2f1_dd(x), y, hyp2f1_dd(c), z, &t);
          else
            {
              const double lp = (c - 2.0 * x) * gsl_sf_log1p(-z);
              stat = hyp2f1_series_conj_dd(gsl_dd_two_sum(c, -x), -y, hyp2f1_dd(c), z, &t);
              t.l += lp;
              t.e += fabs(t.v) * GSL_DBL_EPSILON * (fabs(lp) + 1.0);
            }
          if (stat == GSL_SUCCESS && gsl_finite(t.v)
              && (!have || hyp2f1_relerr(&t) < hyp2f1_relerr(&best)))
            {
              best = t;
              have = 1;
            }
        }
    }

  *r = best;
  return have ? GSL_SUCCESS : (last != GSL_SUCCESS ? last : GSL_EMAXITER);
}

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_hyperg_2F1_e(const double a, const double b, const double c,
                    const double x,
                    gsl_sf_result * result)
{
  hyp2f1_scaled r;
  const int stat = hyp2f1_scaled_eval(a, b, c, x, &r);
  return hyp2f1_result(stat, &r, result);
}

/*
 * F(a,b;c;x)/Gamma(c).  For c = -n a non-positive integer the limit is
 * (DLMF 15.2.3_5)
 *
 *   (a)_{n+1} (b)_{n+1} / (n+1)! x^{n+1} F(a+n+1, b+n+1; n+2; x).
 */
int
gsl_sf_hyperg_2F1_renorm_e(const double a, const double b, const double c,
                           const double x,
                           gsl_sf_result * result)
{
  hyp2f1_scaled r;
  int stat;

  if (gsl_isnan(a) || gsl_isnan(b) || gsl_isnan(c) || gsl_isnan(x))
    {
      DOMAIN_ERROR(result);
    }

  if (hyp2f1_isnpint(c))
    {
      const double n1 = 1.0 - c;           /* n+1 */
      double p = 1.0, k;

      if (x == 0.0)
        {
          result->val = 0.0;
          result->err = 0.0;
          return GSL_SUCCESS;
        }

      for (k = 0.0; k < n1; k += 1.0)
        p *= (a + k) * (b + k) / (k + 1.0);

      if (p == 0.0)
        {
          /* a or b a non-positive integer > c: the series terminates first */
          result->val = 0.0;
          result->err = 0.0;
          return GSL_SUCCESS;
        }

      stat = hyp2f1_scaled_eval(a + n1, b + n1, n1 + 1.0, x, &r);
      if (stat == GSL_SUCCESS)
        {
          const double lx = n1 * log(fabs(x));
          if (x < 0.0 && fmod(n1, 2.0) == 1.0)
            p = -p;
          hyp2f1_scale(&r, p, 4.0 * GSL_DBL_EPSILON * n1, lx);
          r.e += fabs(r.v) * GSL_DBL_EPSILON * fabs(lx);
        }
      return hyp2f1_result(stat, &r, result);
    }

  stat = hyp2f1_scaled_eval(a, b, c, x, &r);
  if (stat == GSL_SUCCESS)
    {
      if (fabs(c) < 160.0)
        {
          gsl_sf_result g;
          gsl_sf_gammainv_e(c, &g);
          hyp2f1_scale(&r, g.val, (g.val != 0.0) ? GSL_MIN_DBL(g.err / fabs(g.val), HYP2F1_GAMMA_REL) : 0.0, 0.0);
        }
      else
        {
          gsl_sf_result lg;
          double sgn;
          gsl_sf_lngamma_sgn_e(c, &lg, &sgn);
          hyp2f1_scale(&r, sgn, lg.err + GSL_DBL_EPSILON * fabs(lg.val), -lg.val);
        }
    }

  return hyp2f1_result(stat, &r, result);
}

int
gsl_sf_hyperg_2F1_conj_e(const double aR, const double aI, const double c,
                         const double x,
                         gsl_sf_result * result)
{
  hyp2f1_scaled r;
  const int stat = hyp2f1_conj_scaled_eval(aR, aI, c, x, &r);
  return hyp2f1_result(stat, &r, result);
}

/*
 * F(a,b;c;x)/Gamma(c) for a = aR+i aI, b = aR-i aI.  For c = -n the
 * limit (DLMF 15.2.3_5) has (a)_{n+1} (b)_{n+1} = |(a)_{n+1}|^2 > 0.
 */
int
gsl_sf_hyperg_2F1_conj_renorm_e(const double aR, const double aI, const double c,
                                const double x,
                                gsl_sf_result * result)
{
  hyp2f1_scaled r;
  int stat;

  if (gsl_isnan(aR) || gsl_isnan(aI) || gsl_isnan(c) || gsl_isnan(x))
    {
      DOMAIN_ERROR(result);
    }

  if (hyp2f1_isnpint(c))
    {
      const double n1 = 1.0 - c;
      double p = 1.0, k;

      if (x == 0.0)
        {
          result->val = 0.0;
          result->err = 0.0;
          return GSL_SUCCESS;
        }

      for (k = 0.0; k < n1; k += 1.0)
        p *= ((aR + k) * (aR + k) + aI * aI) / (k + 1.0);

      if (p == 0.0)
        {
          result->val = 0.0;
          result->err = 0.0;
          return GSL_SUCCESS;
        }

      stat = hyp2f1_conj_scaled_eval(aR + n1, aI, n1 + 1.0, x, &r);
      if (stat == GSL_SUCCESS)
        {
          const double lx = n1 * log(fabs(x));
          if (x < 0.0 && fmod(n1, 2.0) == 1.0)
            p = -p;
          hyp2f1_scale(&r, p, 4.0 * GSL_DBL_EPSILON * n1, lx);
          r.e += fabs(r.v) * GSL_DBL_EPSILON * fabs(lx);
        }
      return hyp2f1_result(stat, &r, result);
    }

  stat = hyp2f1_conj_scaled_eval(aR, aI, c, x, &r);
  if (stat == GSL_SUCCESS)
    {
      if (fabs(c) < 160.0)
        {
          gsl_sf_result g;
          gsl_sf_gammainv_e(c, &g);
          hyp2f1_scale(&r, g.val, (g.val != 0.0) ? GSL_MIN_DBL(g.err / fabs(g.val), HYP2F1_GAMMA_REL) : 0.0, 0.0);
        }
      else
        {
          gsl_sf_result lg;
          double sgn;
          gsl_sf_lngamma_sgn_e(c, &lg, &sgn);
          hyp2f1_scale(&r, sgn, lg.err + GSL_DBL_EPSILON * fabs(lg.val), -lg.val);
        }
    }

  return hyp2f1_result(stat, &r, result);
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

double gsl_sf_hyperg_2F1(const double a, const double b, const double c, const double x)
{
  EVAL_RESULT(gsl_sf_hyperg_2F1_e(a, b, c, x, &result));
}

double gsl_sf_hyperg_2F1_conj(const double aR, const double aI, const double c, const double x)
{
  EVAL_RESULT(gsl_sf_hyperg_2F1_conj_e(aR, aI, c, x, &result));
}

double gsl_sf_hyperg_2F1_renorm(const double a, const double b, const double c, const double x)
{
  EVAL_RESULT(gsl_sf_hyperg_2F1_renorm_e(a, b, c, x, &result));
}

double gsl_sf_hyperg_2F1_conj_renorm(const double aR, const double aI, const double c, const double x)
{
  EVAL_RESULT(gsl_sf_hyperg_2F1_conj_renorm_e(aR, aI, c, x, &result));
}
