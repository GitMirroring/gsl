/* specfunc/hyperg/hyperg_1F1.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 * Copyright (C) 2010 Brian Gough
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
 * Kummer's confluent hypergeometric function M(a,b,x) = 1F1(a;b;x) and
 * M(a,b,x)/Gamma(b), for real a, b, x.
 *
 * Methods (the candidates for a case are evaluated and the one with the
 * smallest error estimate is returned):
 *
 *   - the Maclaurin series, with the sum of the magnitudes of its terms
 *     tracked to measure cancellation, re-summed in double-double when
 *     the cancellation is moderate;
 *
 *   - Kummer's transformation M(a,b,x) = e^x M(b-a,b,-x), which for
 *     x < 0 gives a series of positive terms whenever b > 0 and b > a;
 *
 *   - the asymptotic expansion for large x (DLMF 13.7.2), where the
 *     exponentially smaller second term is accounted for in the error;
 *
 *   - three-term recurrences in a (DLMF 13.3.1) and in b (13.3.2), run
 *     as compensated recurrences (accurate to about eps^2 per step,
 *     with the amplification of errors measured by probe solutions) in
 *     their stable directions: for a < 0 < b, x > 0 (the oscillatory
 *     region) M is computed by backward recurrence in a from values with
 *     a in (0, 1], and for b < 0 (not an integer) by backward recurrence
 *     in b from b > 0, since M(a,b+n,x) is the minimal solution as n
 *     increases (Pearson, Olver and Porter 2017, eq. (3.18)).
 *
 * Values are carried as v e^l to avoid premature overflow.
 *
 * References:
 *
 *   J. W. Pearson, S. Olver and M. A. Porter, "Numerical methods for the
 *   computation of the confluent and Gauss hypergeometric functions",
 *   Numer. Algorithms 74 (2017) 821-866;
 *
 *   A. Deano, J. Segura and N. M. Temme, "Computational properties of
 *   three-term recurrence relations for Kummer functions", J. Comput.
 *   Appl. Math. 233 (2010) 1505-1510;
 *
 *   NIST Digital Library of Mathematical Functions, 13.2, 13.3, 13.7.
 */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_dd.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_psi.h>

#include "error.h"

/* maximum number of series terms and recurrence steps */
#define HYP1F1_MAXITER      1000000

/* accept a result without trying further methods when its estimated
   relative error is below this */
#define HYP1F1_GOOD_REL     (64.0 * GSL_DBL_EPSILON)

/* re-sum a series in double-double when its cancellation factor
   (sum of the magnitudes of the terms over the magnitude of the sum) is
   between these */
#define HYP1F1_DD_AMP_MIN   1.0e2
#define HYP1F1_DD_AMP_MAX   1.0e15

/* rounding error of a sum of terms computed as running products: see
   specfunc/hyperg/hyperg_2F1.c (HYP2F1_SUM_ERR) */
#define HYP1F1_SUM_ERR(sabs, skabs) \
  (2.0 * (sabs) + GSL_MAX_DBL(8.0 * sqrt((sabs) * (skabs)), (skabs)))

/* rescale running values beyond this magnitude */
#define HYP1F1_BIG          1.0e150

/*
 * A value v 2^e2 exp(l + lx x): the scaling by a power of two is exact,
 * and the factor e^x (lx = 0 or 1, x the argument of the function that
 * produced the value) is kept apart, so that the e^x of Kummer's
 * transformation and the e^-x of a result for -x cancel exactly instead
 * of leaving the rounding error of x in a logarithm (about eps |x|
 * relative).
 */
typedef struct
{
  double v;
  double e2;     /* integer-valued */
  double l;
  double lx;
  double e;      /* absolute error of v */
  double amp;    /* cancellation factor (1 if none) */
} hyp1f1_scaled;

/* r = 1 exactly */
static void
hyp1f1_set_unit (hyp1f1_scaled * r)
{
  r->v = 1.0;
  r->e2 = 0.0;
  r->l = 0.0;
  r->lx = 0.0;
  r->e = 0.0;
  r->amp = 1.0;
}

static int hyp1f1_pos (const double a, const double alo, const double b, const double blo,
                       const double x, hyp1f1_scaled * r);

/* Throughout, a parameter obtained by a transformation (b - a, a + n,
   b + m) is passed as a + alo or b + blo, alo and blo the rounding
   residuals, so that every method evaluates M at the exact parameter. */

static int
hyp1f1_isnpint (const double x)
{
  return x <= 0.0 && x == floor (x);
}

static double
hyp1f1_relerr (const hyp1f1_scaled * r)
{
  if (r->v == 0.0)
    return (r->e == 0.0) ? 0.0 : GSL_POSINF;
  return r->e / fabs (r->v);
}

/* keep the better of two results (smaller relative error) */
static void
hyp1f1_keep_better (hyp1f1_scaled * best, int * have, const int stat, const hyp1f1_scaled * t)
{
  if (stat == GSL_SUCCESS && gsl_finite (t->v) && gsl_finite (t->e)
      && (!*have || hyp1f1_relerr (t) < hyp1f1_relerr (best)))
    {
      *best = *t;
      *have = 1;
    }
}

/*
 * The Maclaurin series sum_k (a)_k / (b)_k x^k / k!, a = a + alo,
 * b = b + blo, in double precision.  b is not a non-positive integer
 * unless the series terminates first (a a non-positive integer with
 * a >= b).
 */
static int
hyp1f1_series (const double a, const double alo, const double b, const double blo,
               const double x, hyp1f1_scaled * r)
{
  /* the term ratio can change sign or grow again until k passes -a and -b */
  const double kcross = GSL_MAX_DBL (GSL_MAX_DBL (-a, -b), 0.0) + 2.0;
  const int terminating = hyp1f1_isnpint (a) && alo == 0.0;
  double t = 1.0, s = 1.0, sabs = 1.0, skabs = 0.0, k;

  hyp1f1_set_unit (r);

  for (k = 0.0; k < HYP1F1_MAXITER; k += 1.0)
    {
      const double ak = (a + k) + alo;
      double ratio;

      if (ak == 0.0 && terminating)
        break;

      ratio = ak / (((b + k) + blo) * (k + 1.0)) * x;
      t *= ratio;
      s += t;
      sabs += fabs (t);
      skabs += fabs (t) * (k + 1.0);

      if (!(fabs (t) <= HYP1F1_BIG))
        {
          if (!gsl_finite (t))
            {
              r->v = GSL_NAN;
              return GSL_EOVRFLW;
            }
          /* rescale by 2^-500 */
          t = ldexp (t, -500);
          s = ldexp (s, -500);
          sabs = ldexp (sabs, -500);
          skabs = ldexp (skabs, -500);
          r->e2 += 500.0;
        }

      if (!terminating && k > kcross && fabs (t) < 0.25 * GSL_DBL_EPSILON * fabs (s))
        {
          /* geometric tail bound, valid while the ratios decrease */
          const double rr = fabs ((a + k + 1.0) / ((b + k + 1.0) * (k + 2.0)) * x);
          const double rr2 = fabs ((a + k + 2.0) / ((b + k + 2.0) * (k + 3.0)) * x);
          if (rr < 1.0 && rr2 <= rr)
            {
              const double tail = fabs (t) * rr / (1.0 - rr);
              if (tail < 0.25 * GSL_DBL_EPSILON * fabs (s))
                {
                  r->v = s;
                  r->e = GSL_DBL_EPSILON * HYP1F1_SUM_ERR (sabs, skabs + 1.0) + 2.0 * tail;
                  r->amp = sabs / fabs (s);
                  return GSL_SUCCESS;
                }
            }
        }
    }

  if (k >= HYP1F1_MAXITER)
    {
      r->v = s;
      r->e = GSL_POSINF;
      r->amp = GSL_POSINF;
      return GSL_EMAXITER;
    }

  r->v = s;
  r->e = GSL_DBL_EPSILON * HYP1F1_SUM_ERR (sabs, skabs + 1.0);
  r->amp = (s != 0.0) ? sabs / fabs (s) : GSL_POSINF;
  return GSL_SUCCESS;
}

/*
 * The same series summed in double-double, with a, b, x given exactly
 * as double-double values: used when the double sum shows moderate
 * cancellation, which the dd rounding errors (about 1e-32 of the terms)
 * leave harmless up to a cancellation factor of about 1e15.
 */
static int
hyp1f1_series_dd (const gsl_dd a, const gsl_dd b, const gsl_dd x, hyp1f1_scaled * r)
{
  const double kcross = GSL_MAX_DBL (GSL_MAX_DBL (-a.hi, -b.hi), 0.0) + 2.0;
  const int terminating = hyp1f1_isnpint (a.hi) && a.lo == 0.0;
  gsl_dd t, s;
  double sabs = 1.0, tail = 0.0, k;

  t.hi = 1.0;
  t.lo = 0.0;
  s = t;
  hyp1f1_set_unit (r);

  for (k = 0.0; k < HYP1F1_MAXITER; k += 1.0)
    {
      const gsl_dd ak = gsl_dd_add_d (a, k);
      gsl_dd ratio;

      if (ak.hi == 0.0 && terminating)
        break;

      ratio = gsl_dd_div (gsl_dd_mul (ak, x), gsl_dd_mul_d (gsl_dd_add_d (b, k), k + 1.0));
      t = gsl_dd_mul (t, ratio);
      s = gsl_dd_add (s, t);
      sabs += fabs (t.hi);

      if (!(fabs (t.hi) <= HYP1F1_BIG))
        {
          if (!gsl_finite (t.hi))
            return GSL_EOVRFLW;
          t.hi = ldexp (t.hi, -500);
          t.lo = ldexp (t.lo, -500);
          s.hi = ldexp (s.hi, -500);
          s.lo = ldexp (s.lo, -500);
          sabs = ldexp (sabs, -500);
          r->e2 += 500.0;
        }

      if (!terminating && k > kcross && fabs (t.hi) < 1.0e-20 * fabs (s.hi))
        {
          const double rr = fabs ((a.hi + k + 1.0) / ((b.hi + k + 1.0) * (k + 2.0)) * x.hi);
          const double rr2 = fabs ((a.hi + k + 2.0) / ((b.hi + k + 2.0) * (k + 3.0)) * x.hi);
          if (rr < 1.0 && rr2 <= rr && fabs (t.hi) * rr / (1.0 - rr) < 1.0e-20 * fabs (s.hi))
            {
              tail = 2.0 * fabs (t.hi) * rr / (1.0 - rr);
              break;
            }
        }
    }

  if (k >= HYP1F1_MAXITER)
    return GSL_EMAXITER;

  r->v = s.hi + s.lo;
  r->e = 1.0e-31 * (2.0 + 6.0 * sqrt (k + 1.0)) * sabs + tail
         + 0.5 * GSL_DBL_EPSILON * fabs (r->v);
  r->amp = 1.0;
  return GSL_SUCCESS;
}

/*
 * Large x > 0 (DLMF 13.7.2):
 *
 *   M(a,b,x) ~ Gamma(b)/Gamma(a) e^x x^(a-b) sum_s (1-a)_s (b-a)_s / s! x^-s
 *             + Gamma(b)/Gamma(b-a) e^(+-i pi a) x^-a sum_s ...
 *
 * The first term only is summed; the second, smaller by about
 * |Gamma(a)/Gamma(b-a)| x^(b-2a) e^-x, is added to the error, so that
 * the expansion is used only where it is negligible.  The series is
 * summed up to its smallest term, which bounds the truncation error.
 */
static int
hyp1f1_asymp (const double a, const double alo, const double b, const double blo,
              const double x, hyp1f1_scaled * r)
{
  gsl_sf_result lgb, lga, lgba;
  double sgb, sga, sgba;
  double t = 1.0, s = 1.0, sabs = 1.0, k, prev = GSL_POSINF;
  int stat;

  /* the expansion is for x large compared with the parameters */
  if ((hyp1f1_isnpint (a) && alo == 0.0) || x < 2.0 * (fabs (a) + fabs (b)) + 10.0)
    return GSL_EDOM;

  stat = gsl_sf_lngamma_sgn_e (b, &lgb, &sgb);
  stat += gsl_sf_lngamma_sgn_e (a, &lga, &sga);
  if (stat != GSL_SUCCESS)
    return GSL_EDOM;

  for (k = 0.0; k < 1000.0; k += 1.0)
    {
      const double next = t * (1.0 - a + k) * (b - a + k) / ((k + 1.0) * x);
      if (fabs (next) >= prev || fabs (next) >= fabs (t))
        break;                  /* terms start growing: stop at the smallest */
      prev = fabs (t);
      t = next;
      s += t;
      sabs += fabs (t);
      if (fabs (t) < 0.25 * GSL_DBL_EPSILON * fabs (s))
        break;
    }

  /* the smallest term must be small: otherwise the expansion is useless
     here (for a << 0 its terms grow from the start) */
  if (!(fabs (t) < 1.0e-6 * fabs (s)))
    return GSL_EDOM;

  r->v = sgb * sga * s;
  r->e2 = 0.0;
  r->l = lgb.val - lga.val + (a - b) * log (x);
  r->lx = 1.0;
  r->e = GSL_DBL_EPSILON * (4.0 + 2.0 * k) * sabs + 2.0 * fabs (t);

  /* error of the prefactor: lngamma errors and rounding of the exponent */
  r->e += fabs (s) * (lgb.err + lga.err
                      + GSL_DBL_EPSILON * (fabs (lgb.val) + fabs (lga.val)
                                           + 2.0 * fabs ((a - b) * log (x))));

  /* the residual alo of a: d/da ln(x^a / Gamma(a)) = log(x) - psi(a) */
  if (alo != 0.0)
    {
      gsl_sf_result psi;
      if (gsl_sf_psi_e (a, &psi) != GSL_SUCCESS)
        return GSL_EDOM;
      r->e += fabs (s) * fabs (alo) * (fabs (log (x)) + fabs (psi.val));
    }

  /* and of b: d/db ln(Gamma(b) x^-b) = psi(b) - log(x) */
  if (blo != 0.0)
    {
      gsl_sf_result psi;
      if (gsl_sf_psi_e (b, &psi) != GSL_SUCCESS)
        return GSL_EDOM;
      r->e += fabs (s) * fabs (blo) * (fabs (log (x)) + fabs (psi.val));
    }

  /* the omitted second term, relative to the first */
  if (!hyp1f1_isnpint (b - a))
    {
      stat = gsl_sf_lngamma_sgn_e (b - a, &lgba, &sgba);
      if (stat != GSL_SUCCESS)
        return GSL_EDOM;
      {
        const double lrel = lga.val - lgba.val + (b - 2.0 * a) * log (x) - x;
        if (lrel > -2.0)
          return GSL_EDOM;      /* not negligible: other methods */
        r->e += 2.0 * fabs (s) * exp (lrel);
      }
    }

  r->amp = 1.0;
  return GSL_SUCCESS;
}

/*
 * State of a three-term recurrence f_{k+1} = (p f_k + q f_{k-1}) / d,
 * run as a compensated recurrence (f = f.hi + f.lo, the rounding error
 * of each step computed by error-free transformations) with
 * the coefficients p, q, d formed to a few eps^2, so that the result
 * is about as accurate as in double-double arithmetic, together with
 * two probe solutions g, h (in double) of
 * the same recurrence, started from unit relative perturbations of the
 * two starting values (g from f1, h from f0).  Since the recurrence is
 * linear, the errors of the starting values propagate as a combination
 * of g and h, which measure their amplification exactly (about |f|
 * where the computed solution is dominant, more where the recurrence is
 * unstable); the rounding errors of the compensated steps, which are
 * also amplified by these factors, are negligible against them.  (In
 * double precision the rounding errors made on the way can be amplified
 * more than errors of the starting values, near turning points.)
 */
typedef struct
{
  gsl_dd f0, f1;     /* f0 = f_k, f1 = f_{k-1} */
  double g0, g1;
  double h0, h1;
  double e2, l, lx;  /* common scale (as in hyp1f1_scaled) */
  double rerr;       /* accumulated local rounding errors, relative to the
                        amplification scale max(|f|, |g|, |h|) */
} hyp1f1_rec;

static void
hyp1f1_rec_init (hyp1f1_rec * w, const hyp1f1_scaled * m1, const hyp1f1_scaled * m0,
                 const double x)
{
  /* m0 = f_k (the value the recurrence starts from), m1 = f_{k-1} */
  const double dl = (m1->l - m0->l) + (m1->lx - m0->lx) * x;
  w->e2 = m0->e2;
  w->l = m0->l;
  w->lx = m0->lx;
  w->f0.hi = m0->v;
  w->f0.lo = 0.0;
  w->f1.hi = ldexp ((dl == 0.0) ? m1->v : m1->v * exp (dl), (int) (m1->e2 - m0->e2));
  w->f1.lo = 0.0;
  w->g0 = 0.0;
  w->g1 = fabs (w->f1.hi);
  w->h0 = fabs (w->f0.hi);
  w->h1 = 0.0;
  w->rerr = 0.0;
}

static void
hyp1f1_rec_step (hyp1f1_rec * w, const gsl_dd p, const gsl_dd q, const gsl_dd d)
{
  const double pd = p.hi, qd = q.hi, dd = d.hi;
  gsl_dd f;
  const double g = (pd * w->g0 + qd * w->g1) / dd;
  const double h = (pd * w->h0 + qd * w->h1) / dd;
  double scale, m;

  {
    /* f = (p f0 + q f1) / d, compensated: the products and the sum of
       the leading parts are formed exactly (error-free transformations),
       their rounding errors and the low parts are collected in double,
       and the remainder of the division is computed exactly, so that
       f.hi + f.lo has a relative error of a few eps^2 (of the terms) */
    const gsl_dd tp = gsl_dd_two_prod (pd, w->f0.hi);
    const gsl_dd tq = gsl_dd_two_prod (qd, w->f1.hi);
    const gsl_dd sn = gsl_dd_two_sum (tp.hi, tq.hi);
    const double nlo = sn.lo + tp.lo + tq.lo + pd * w->f0.lo + qd * w->f1.lo
                       + p.lo * w->f0.hi + q.lo * w->f1.hi;
    const double fh = sn.hi / dd;
    const gsl_dd fd = gsl_dd_two_prod (fh, dd);
    const double fl = (((sn.hi - fd.hi) - fd.lo) + nlo - fh * d.lo) / dd;
    f = gsl_dd_two_sum (fh, fl);
  }

  scale = GSL_MAX_DBL (fabs (f.hi), GSL_MAX_DBL (fabs (g), fabs (h)));

  /* local rounding of the compensated step (relative to the terms,
     which may cancel), propagating like the probes */
  if (scale > 0.0)
    w->rerr += 1.0e-30 * (fabs (pd * w->f0.hi) + fabs (qd * w->f1.hi)) / (fabs (dd) * scale);

  w->f1 = w->f0;
  w->f0 = f;
  w->g1 = w->g0;
  w->g0 = g;
  w->h1 = w->h0;
  w->h0 = h;

  m = GSL_MAX_DBL (GSL_MAX_DBL (fabs (w->f0.hi), fabs (w->f1.hi)),
                   GSL_MAX_DBL (GSL_MAX_DBL (fabs (w->g0), fabs (w->g1)),
                                GSL_MAX_DBL (fabs (w->h0), fabs (w->h1))));
  if (m > HYP1F1_BIG || (m < 1.0 / HYP1F1_BIG && m > 0.0))
    {
      int e;
      frexp (m, &e);
      w->f0.hi = ldexp (w->f0.hi, -e);
      w->f0.lo = ldexp (w->f0.lo, -e);
      w->f1.hi = ldexp (w->f1.hi, -e);
      w->f1.lo = ldexp (w->f1.lo, -e);
      w->g0 = ldexp (w->g0, -e);
      w->g1 = ldexp (w->g1, -e);
      w->h0 = ldexp (w->h0, -e);
      w->h1 = ldexp (w->h1, -e);
      w->e2 += e;
    }
}

/*
 * Recurrence coefficients, as hi + lo with a relative error of a few
 * eps^2.  They are formed from c + j for an integer-valued j:
 * computing c + j in double would round it with an absolute error of
 * about eps |c + j|, of the same sign over long stretches of j, which
 * acts as a systematic perturbation of the parameter.
 */

/* a + b + lo (lo small), normalized */
static gsl_dd
hyp1f1_shift (const double a, const double b, const double lo)
{
  const gsl_dd t = gsl_dd_two_sum (a, b);
  return gsl_dd_two_sum (t.hi, t.lo + lo);
}

/* (a.hi + a.lo) (b.hi + b.lo) */
static gsl_dd
hyp1f1_prod (const gsl_dd a, const gsl_dd b)
{
  const gsl_dd t = gsl_dd_two_prod (a.hi, b.hi);
  return gsl_dd_two_sum (t.hi, t.lo + (a.hi * b.lo + a.lo * b.hi));
}


/* result of a recurrence of n steps whose starting values have relative
   errors rel0 (sum of both) */
static int
hyp1f1_rec_result (const hyp1f1_rec * w, const double rel0, const double n, hyp1f1_scaled * r)
{
  const double f = w->f0.hi + w->f0.lo;
  const double amp = GSL_MAX_DBL (fabs (f), GSL_MAX_DBL (fabs (w->g0), fabs (w->h0)));

  if (!gsl_finite (f) || !gsl_finite (amp))
    return GSL_EOVRFLW;

  r->v = f;
  r->e2 = w->e2;
  r->l = w->l;
  r->lx = w->lx;
  r->e = (rel0 + w->rerr + 2.0 * GSL_DBL_EPSILON) * amp;
  (void) n;
  r->amp = (f != 0.0) ? amp / fabs (f) : GSL_POSINF;
  return GSL_SUCCESS;
}

/*
 * Backward recurrence in b (DLMF 13.3.2),
 *
 *   M(c-1) = -[c (1 - c - x) M(c) + x (c - a) M(c+1)] / (c (c - 1)),
 *
 * from M(a, b+m) = m0 and M(a, b+m+1) = m1 down to M(a, b).  M is the
 * minimal solution as b increases (Pearson, Olver and Porter 2017, eq.
 * (3.18)), so this direction is stable: errors of the starting values
 * are not amplified.  The coefficients use c = b + j + blo as hi + lo,
 * so that the last step is exactly at b + blo (m0, m1
 * are at the exact b + blo + m, b + blo + m + 1).  b + j must not be 0
 * or 1 (b not an integer).
 */
static int
hyp1f1_run_b_down (const double a, const double alo, const double b, const double blo,
                   const double x, const double m,
                   const hyp1f1_scaled * m0, const hyp1f1_scaled * m1, hyp1f1_scaled * r)
{
  const double rel0 = hyp1f1_relerr (m0) + hyp1f1_relerr (m1)
                      + GSL_DBL_EPSILON * fabs ((m1->l - m0->l) + (m1->lx - m0->lx) * x);
  const gsl_dd xm1 = gsl_dd_two_sum (x, -1.0);
  gsl_dd mx;
  hyp1f1_rec w;
  double j;

  mx.hi = -x;
  mx.lo = 0.0;
  hyp1f1_rec_init (&w, m1, m0, x);

  for (j = m; j > 0.0; j -= 1.0)
    {
      /* f0 = M(c), f1 = M(c+1), c = b + j */
      const gsl_dd cd = hyp1f1_shift (b, j, blo);
      const gsl_dd omcx = hyp1f1_shift (cd.hi, xm1.hi, cd.lo + xm1.lo);   /* c + x - 1 */
      const gsl_dd cma = hyp1f1_shift (cd.hi, -a, cd.lo - alo);           /* c - a */
      const gsl_dd cm1 = hyp1f1_shift (cd.hi, -1.0, cd.lo);               /* c - 1 */
      hyp1f1_rec_step (&w, hyp1f1_prod (cd, omcx), hyp1f1_prod (cma, mx),
                       hyp1f1_prod (cd, cm1));
    }

  return hyp1f1_rec_result (&w, rel0, m, r);
}

/*
 * a < 0 < b, x < b: backward recurrence in a (DLMF 13.3.1),
 *
 *   M(c-1) = [c M(c+1) - (2c - b + x) M(c)] / (b - c),
 *
 * from M(a0) and M(a0+1), a0 = a + n in (0, 1] (positive cases), down to
 * a.  For x >= b this direction is unstable (see hyp1f1_rec_a()); the
 * error estimate follows the amplification of errors (hyp1f1_rec).
 */
static int
hyp1f1_rec_a_direct (const double a, const double alo, const double b, const double blo,
                     const double x, hyp1f1_scaled * r)
{
  /* a0 + alo in [0, 1); a0 = a + n is exact (Sterbenz) */
  const double n = (hyp1f1_isnpint (a) && alo < 0.0) ? 1.0 - a : ceil (-a);
  const double a0 = a + n;
  hyp1f1_scaled m0, m1;
  hyp1f1_rec w;
  gsl_dd bmx;
  double rel0, k;
  int stat;

  if (n > HYP1F1_MAXITER)
    return GSL_EMAXITER;

  if (a0 == 0.0 && alo == 0.0)
    {
      /* a a non-positive integer: M(0) = 1 */
      hyp1f1_set_unit (&m0);
      stat = GSL_SUCCESS;
    }
  else
    {
      stat = hyp1f1_pos (a0, alo, b, blo, x, &m0);
    }
  stat += hyp1f1_pos (a0 + 1.0, alo, b, blo, x, &m1);
  if (stat != GSL_SUCCESS)
    return GSL_EDOM;

  rel0 = hyp1f1_relerr (&m0) + hyp1f1_relerr (&m1)
         + GSL_DBL_EPSILON * fabs ((m1.l - m0.l) + (m1.lx - m0.lx) * x);
  bmx = hyp1f1_shift (b, -x, blo);
  hyp1f1_rec_init (&w, &m1, &m0, x);

  for (k = 0.0; k < n; k += 1.0)
    {
      /* f0 = M(c), f1 = M(c+1), c = a0 - k */
      const gsl_dd cd = hyp1f1_shift (a0, -k, alo);
      const gsl_dd p = hyp1f1_shift (-2.0 * cd.hi, bmx.hi, bmx.lo - 2.0 * cd.lo);  /* -(2c - b + x) */
      const gsl_dd d = hyp1f1_shift (b, -cd.hi, blo - cd.lo);                       /* b - c */
      hyp1f1_rec_step (&w, p, cd, d);
    }

  return hyp1f1_rec_result (&w, rel0, n, r);
}

/*
 * The terminating series for a = -m, b = -n a non-positive integer with
 * b <= a <= -2 (the polynomial, in which the series stops before the
 * zero of (b)_k): backward recurrence in a from M(0) = 1 and
 * M(-1) = 1 - x/b, all of whose values are such polynomials, so that
 * the contiguous relation holds between them.
 */
static int
hyp1f1_rec_a_poly (const double a, const double b, const double x, hyp1f1_scaled * r)
{
  hyp1f1_scaled m0, m1;
  hyp1f1_rec w;
  const gsl_dd bmx = gsl_dd_two_sum (b, -x);
  gsl_dd mm1;
  double k;

  hyp1f1_set_unit (&m1);        /* M(0) */
  hyp1f1_set_unit (&m0);        /* M(-1) = 1 - x/b */
  mm1 = gsl_dd_sub (gsl_dd_two_sum (1.0, 0.0), gsl_dd_div_d (gsl_dd_two_sum (x, 0.0), b));
  m0.v = mm1.hi;
  m0.e = 0.0;

  hyp1f1_rec_init (&w, &m1, &m0, x);
  w.f0 = mm1;                   /* keep M(-1) to double-double accuracy */

  for (k = -1.0; k > a; k -= 1.0)
    {
      /* f0 = M(c), f1 = M(c+1), c = k (exact, as is b - c) */
      const gsl_dd cd = gsl_dd_two_sum (k, 0.0);
      const gsl_dd p = hyp1f1_shift (-2.0 * k, bmx.hi, bmx.lo);   /* -(2c - b + x) */
      const gsl_dd d = gsl_dd_two_sum (b - k, 0.0);                        /* b - c */
      hyp1f1_rec_step (&w, p, cd, d);
    }

  return hyp1f1_rec_result (&w, 1.0e-30, -a, r);
}

/*
 * a < 0 < b, x > 0.  For x >= b the backward recurrence in a is
 * unstable; b is then first raised to b' = b + m > x, M(a, b') and
 * M(a, b'+1) are computed by the recurrence in a, and the stable
 * backward recurrence in b brings them down to b.
 */
static int
hyp1f1_rec_a (const double a, const double alo, const double b, const double blo,
              const double x, hyp1f1_scaled * r)
{
  if (x < b)
    {
      return hyp1f1_rec_a_direct (a, alo, b, blo, x, r);
    }
  else
    {
      const double m = ceil (x - b) + 1.0;
      const gsl_dd bp = hyp1f1_shift (b, m, blo);
      const gsl_dd bp1 = hyp1f1_shift (bp.hi, 1.0, bp.lo);
      hyp1f1_scaled m0, m1;
      int stat;

      if (m > HYP1F1_MAXITER)
        return GSL_EMAXITER;

      stat = hyp1f1_rec_a_direct (a, alo, bp.hi, bp.lo, x, &m0);
      stat += hyp1f1_rec_a_direct (a, alo, bp1.hi, bp1.lo, x, &m1);
      if (stat != GSL_SUCCESS)
        return GSL_EDOM;

      return hyp1f1_run_b_down (a, alo, b, blo, x, m, &m0, &m1, r);
    }
}

/*
 * b < 0, not an integer, x > 0: backward recurrence in b from b + n in
 * (1, 2], n = ceil(-b) + 1.
 */
static int
hyp1f1_rec_b (const double a, const double alo, const double b, const double blo,
              const double x, hyp1f1_scaled * r)
{
  const double n = ceil (-b) + 1.0;
  const gsl_dd bn = hyp1f1_shift (b, n, blo);
  const gsl_dd bn1 = hyp1f1_shift (bn.hi, 1.0, bn.lo);
  hyp1f1_scaled m0, m1;
  int stat;

  if (n > HYP1F1_MAXITER)
    return GSL_EMAXITER;

  stat = hyp1f1_pos (a, alo, bn.hi, bn.lo, x, &m0);
  stat += hyp1f1_pos (a, alo, bn1.hi, bn1.lo, x, &m1);
  if (stat != GSL_SUCCESS)
    return GSL_EDOM;

  return hyp1f1_run_b_down (a, alo, b, blo, x, n, &m0, &m1, r);
}

/* the series, re-summed in double-double if it shows moderate cancellation */
static int
hyp1f1_series_best (const double a, const double alo, const double b, const double blo,
                    const double x, hyp1f1_scaled * r)
{
  int stat = hyp1f1_series (a, alo, b, blo, x, r);

  if (stat == GSL_SUCCESS && r->amp > HYP1F1_DD_AMP_MIN && r->amp < HYP1F1_DD_AMP_MAX)
    {
      hyp1f1_scaled t;
      gsl_dd ad, bd, xd;
      ad.hi = a;
      ad.lo = alo;
      bd.hi = b;
      bd.lo = blo;
      xd.hi = x;
      xd.lo = 0.0;
      if (hyp1f1_series_dd (ad, bd, xd, &t) == GSL_SUCCESS && hyp1f1_relerr (&t) < hyp1f1_relerr (r))
        *r = t;
    }

  return stat;
}

/*
 * M(a + alo, b + blo, x) for x > 0 (x >= 0), b > 0 or b not an integer.
 */
static int
hyp1f1_pos (const double a, const double alo, const double b, const double blo,
            const double x, hyp1f1_scaled * r)
{
  hyp1f1_scaled best, t;
  int have = 0, stat;

  if (x == 0.0 || (a == 0.0 && alo == 0.0))
    {
      hyp1f1_set_unit (r);
      return GSL_SUCCESS;
    }

  hyp1f1_set_unit (&best);
  best.v = 0.0;
  best.e = GSL_POSINF;
  best.amp = GSL_POSINF;

  /* the series: no cancellation for a >= 0, b > 0 */
  stat = hyp1f1_series (a, alo, b, blo, x, &t);
  hyp1f1_keep_better (&best, &have, stat, &t);
  if (have && hyp1f1_relerr (&best) <= HYP1F1_GOOD_REL && x < 30.0)
    {
      *r = best;
      return GSL_SUCCESS;
    }

  if (x >= 10.0)
    {
      stat = hyp1f1_asymp (a, alo, b, blo, x, &t);
      hyp1f1_keep_better (&best, &have, stat, &t);
    }

  if (have && hyp1f1_relerr (&best) <= HYP1F1_GOOD_REL)
    {
      *r = best;
      return GSL_SUCCESS;
    }

  if (b > 0.0 && a + alo < 0.0)
    {
      stat = hyp1f1_rec_a (a, alo, b, blo, x, &t);
      hyp1f1_keep_better (&best, &have, stat, &t);
    }
  else if (b < 0.0)
    {
      stat = hyp1f1_rec_b (a, alo, b, blo, x, &t);
      hyp1f1_keep_better (&best, &have, stat, &t);
    }

  if (!have || hyp1f1_relerr (&best) > HYP1F1_GOOD_REL)
    {
      stat = hyp1f1_series_best (a, alo, b, blo, x, &t);
      hyp1f1_keep_better (&best, &have, stat, &t);
    }

  *r = best;
  return have ? GSL_SUCCESS : GSL_EMAXITER;
}

/*
 * M(a+alo,b,x) for b not a non-positive integer (or a terminating
 * series, with alo = 0).
 */
static int
hyp1f1_eval (const double a, const double alo, const double b, const double x, hyp1f1_scaled * r)
{
  hyp1f1_scaled best, t;
  int have = 0, stat;

  hyp1f1_set_unit (&best);
  best.v = 0.0;
  best.e = GSL_POSINF;
  best.amp = GSL_POSINF;

  if (hyp1f1_isnpint (b))
    {
      /* allowed only when the series terminates first (a = -m >= b):
         the polynomial, directly or by recurrence in a */
      stat = hyp1f1_series_best (a, 0.0, b, 0.0, x, &t);
      hyp1f1_keep_better (&best, &have, stat, &t);
      if ((!have || hyp1f1_relerr (&best) > HYP1F1_GOOD_REL) && a <= -2.0)
        {
          stat = hyp1f1_rec_a_poly (a, b, x, &t);
          hyp1f1_keep_better (&best, &have, stat, &t);
        }
      *r = best;
      return have ? GSL_SUCCESS : GSL_ELOSS;
    }

  if (x >= 0.0)
    return hyp1f1_pos (a, alo, b, 0.0, x, r);

  /* x < 0: the series directly (alternating for a, b > 0) ... */
  stat = hyp1f1_series (a, alo, b, 0.0, x, &t);
  hyp1f1_keep_better (&best, &have, stat, &t);
  if (have && hyp1f1_relerr (&best) <= HYP1F1_GOOD_REL)
    {
      *r = best;
      return GSL_SUCCESS;
    }

  /* ... or Kummer's transformation e^x M(b-a, b, -x) */
  {
    gsl_dd ba = gsl_dd_two_sum (b, -a);
    ba = gsl_dd_two_sum (ba.hi, ba.lo - alo);
    stat = hyp1f1_pos (ba.hi, ba.lo, b, 0.0, -x, &t);
    if (stat == GSL_SUCCESS)
      {
        /* t is v 2^e2 exp(l + lx (-x)); times e^x: exp(l + (1 - lx) x) */
        t.lx = 1.0 - t.lx;
      }
    hyp1f1_keep_better (&best, &have, stat, &t);
  }

  if (!have || hyp1f1_relerr (&best) > HYP1F1_GOOD_REL)
    {
      stat = hyp1f1_series_best (a, alo, b, 0.0, x, &t);
      hyp1f1_keep_better (&best, &have, stat, &t);
    }

  *r = best;
  return have ? GSL_SUCCESS : GSL_EMAXITER;
}

/* assemble a gsl_sf_result from a scaled value */
static int
hyp1f1_result (const int stat, const hyp1f1_scaled * r, const double x, gsl_sf_result * result)
{
  int stat_e;

  if (stat != GSL_SUCCESS)
    {
      result->val = GSL_NAN;
      result->err = GSL_NAN;
      GSL_ERROR ("error", stat);
    }

  /* an exact zero (e.g. a terminating polynomial at a root) carries an
     absolute error bound */
  if (r->v == 0.0 && gsl_finite (r->e))
    {
      double e = ldexp (r->e, (int) GSL_MIN_DBL (r->e2, 1000.0));
      double L = r->l + ((r->lx != 0.0) ? r->lx * x : 0.0);
      result->val = 0.0;
      result->err = (L < GSL_LOG_DBL_MAX) ? e * exp (L) : GSL_POSINF;
      return GSL_SUCCESS;
    }

  /* no significant digit (checked before scaling, which could
     overflow for a meaningless value) */
  if (!(r->e < fabs (r->v)))
    {
      result->val = GSL_NAN;
      result->err = GSL_NAN;
      GSL_ERROR ("error", GSL_ELOSS);
    }

  {
    /* v 2^e2 exp(L): the power of 2 goes into v while it stays a normal
       number, otherwise into L */
    double L = r->l + ((r->lx != 0.0) ? r->lx * x : 0.0);
    double dL = GSL_DBL_EPSILON * (fabs (r->l) + ((r->l != 0.0 && r->lx != 0.0) ? fabs (L) : 0.0));
    double v = ldexp (r->v, (int) r->e2), e = ldexp (r->e, (int) r->e2);

    if (r->e2 != 0.0 && (!gsl_finite (v) || fabs (v) < GSL_DBL_MIN))
      {
        v = r->v;
        e = r->e;
        L += r->e2 * M_LN2;
        dL += GSL_DBL_EPSILON * fabs (r->e2 * M_LN2);
      }

    if (L == 0.0)
      {
        result->val = v;
        result->err = e + GSL_DBL_EPSILON * fabs (v);
        stat_e = gsl_finite (v) ? GSL_SUCCESS : GSL_EOVRFLW;
        if (stat_e != GSL_SUCCESS)
          {
            OVERFLOW_ERROR (result);
          }
      }
    else
      {
        stat_e = gsl_sf_exp_mult_err_e (L, dL, v, e, result);
      }
  }

  /* no significant digit: report the loss rather than a meaningless value */
  if (stat_e == GSL_SUCCESS && result->err > 0.0 && !(result->err < fabs (result->val)))
    GSL_ERROR ("error", GSL_ELOSS);

  return stat_e;
}

/* domain checks shared by the public functions; 1 if r was set */
static int
hyp1f1_special (const double a, const double b, const double x, hyp1f1_scaled * r, int * stat)
{
  hyp1f1_set_unit (r);

  if (gsl_isnan (a) || gsl_isnan (b) || !gsl_finite (a) || !gsl_finite (b))
    {
      *stat = GSL_EDOM;
      return 1;
    }

  /* b a non-positive integer: defined only if the series terminates
     before the zero of (b)_k, i.e. a = -n with b <= a <= 0 */
  if (hyp1f1_isnpint (b) && !(hyp1f1_isnpint (a) && a >= b))
    {
      *stat = GSL_EDOM;
      return 1;
    }

  /* M(0,b,x) = 1 for all x, including x = NaN (as pow(NaN, 0) = 1) */
  if (a == 0.0)
    {
      *stat = GSL_SUCCESS;
      return 1;
    }

  if (gsl_isnan (x) || !gsl_finite (x))
    {
      *stat = GSL_EDOM;
      return 1;
    }

  if (x == 0.0)
    {
      *stat = GSL_SUCCESS;
      return 1;
    }

  if (a == b && !hyp1f1_isnpint (a))
    {
      /* e^x (for a = b = -n the series terminates: sum_{k<=n} x^k/k!) */
      r->lx = 1.0;
      *stat = GSL_SUCCESS;
      return 1;
    }

  return 0;
}

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_hyperg_1F1_e (const double a, const double b, const double x, gsl_sf_result * result)
{
  hyp1f1_scaled r;
  int stat;

  if (!hyp1f1_special (a, b, x, &r, &stat))
    stat = hyp1f1_eval (a, 0.0, b, x, &r);

  if (stat == GSL_EDOM)
    {
      DOMAIN_ERROR (result);
    }

  return hyp1f1_result (stat, &r, x, result);
}

int
gsl_sf_hyperg_1F1_int_e (const int m, const int n, const double x, gsl_sf_result * result)
{
  return gsl_sf_hyperg_1F1_e ((double) m, (double) n, x, result);
}

/*
 * M(a,b,x)/Gamma(b).  For b = -n a non-positive integer the limit is
 * (DLMF 13.2.5)
 *
 *   (a)_{n+1} x^{n+1} / (n+1)! M(a+n+1, n+2, x),
 *
 * which is 0 when a is an integer with b <= a <= 0 (terminating series).
 */
int
gsl_sf_hyperg_1F1_renorm_e (const double a, const double b, const double x,
                            gsl_sf_result * result)
{
  hyp1f1_scaled r;
  int stat;

  if (gsl_isnan (a) || gsl_isnan (b) || !gsl_finite (a) || !gsl_finite (b))
    {
      DOMAIN_ERROR (result);
    }

  if (hyp1f1_isnpint (b))
    {
      const double n = -b;
      double p = 1.0, k;
      gsl_dd ap;

      if (hyp1f1_isnpint (a) && a >= b)
        {
          /* terminating series over Gamma(b) = infinity: 0 for all x */
          result->val = 0.0;
          result->err = 0.0;
          return GSL_SUCCESS;
        }

      if (gsl_isnan (x) || !gsl_finite (x))
        {
          DOMAIN_ERROR (result);
        }

      if (x == 0.0)
        {
          result->val = 0.0;
          result->err = 0.0;
          return GSL_SUCCESS;
        }

      /* a + n + 1 as hi + lo: its rounding is amplified by the
         sensitivity of M to a */
      ap = gsl_dd_two_sum (a, n + 1.0);
      if (ap.lo != 0.0)
        stat = hyp1f1_eval (ap.hi, ap.lo, n + 2.0, x, &r);
      else if (!hyp1f1_special (ap.hi, n + 2.0, x, &r, &stat))
        stat = hyp1f1_eval (ap.hi, 0.0, n + 2.0, x, &r);
      if (stat == GSL_EDOM)
        {
          DOMAIN_ERROR (result);
        }

      /* (a)_{n+1} x^{n+1} / (n+1)!, as a product to keep it scaled */
      for (k = 0.0; k <= n; k += 1.0)
        {
          p *= (a + k) * x / (k + 1.0);
          if (!(fabs (p) <= HYP1F1_BIG) || (fabs (p) < 1.0 / HYP1F1_BIG && p != 0.0))
            {
              int e;
              frexp (p, &e);
              p = ldexp (p, -e);
              r.e2 += e;
            }
        }
      r.v *= p;
      r.e = r.e * fabs (p) + 4.0 * GSL_DBL_EPSILON * (n + 1.0) * fabs (r.v);
      return hyp1f1_result (stat, &r, x, result);
    }
  else
    {
      gsl_sf_result lg;
      double sg;

      if (!hyp1f1_special (a, b, x, &r, &stat))
        stat = hyp1f1_eval (a, 0.0, b, x, &r);
      if (stat == GSL_EDOM)
        {
          DOMAIN_ERROR (result);
        }

      if (gsl_sf_lngamma_sgn_e (b, &lg, &sg) != GSL_SUCCESS)
        {
          DOMAIN_ERROR (result);
        }

      r.v *= sg;
      r.l -= lg.val;
      r.e += fabs (r.v) * (lg.err + GSL_DBL_EPSILON * fabs (lg.val));
      return hyp1f1_result (stat, &r, x, result);
    }
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

double
gsl_sf_hyperg_1F1_int (const int m, const int n, const double x)
{
  EVAL_RESULT (gsl_sf_hyperg_1F1_int_e (m, n, x, &result));
}

double
gsl_sf_hyperg_1F1 (const double a, const double b, const double x)
{
  EVAL_RESULT (gsl_sf_hyperg_1F1_e (a, b, x, &result));
}

double
gsl_sf_hyperg_1F1_renorm (const double a, const double b, const double x)
{
  EVAL_RESULT (gsl_sf_hyperg_1F1_renorm_e (a, b, x, &result));
}
