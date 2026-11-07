/* dd/internal_log_source.c
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
 * Add
 *
 *   w = log(1+r) - r = -r^2/2 + r^3/3 - ...
 *
 * to the expansion sum, for a small double-double argument |r| <= 1/128.
 *
 * Write w = -r^2/2 + r^3 Q(r) with
 *
 *   Q(r) = sum_{n=3}^{18} (-1)^(n+1) r^(n-3) / n.
 *
 * The leading term is added exactly: r^2 = r.hi^2 + 2 r.hi r.lo + r.lo^2
 * is split by TwoProd, so w is not limited by the rounding error of a
 * double-double square.  The term r^3 Q(r) has relative weight below
 * 2|r|/3 in w and is evaluated in double-double: by Horner for the
 * coefficients n = 3,...,9 and in double precision from r.hi for
 * n >= 10, whose relative weight in Q is below |r|^7 <= 2^-49.  Since
 * |r c_(n+1)| < |c_n| in every Horner step, the steps use an addition
 * without the cancellation handling of dd_add().  The truncation error
 * is below |r|^16/19 relative to w.
 *
 * Returning w rather than log(1+r) lets the callers add r and w
 * separately into their final expansion, and gives log1pmx without
 * cancellation.
 */
static void
DD_WRAP(log1pmx_small_add)(dd_ex4 * sum, const gsl_dd r)
{
  static const gsl_dd coef[7] = {
    {  3.33333333333333314830e-01,  1.85037170770859413132e-17 },
    { -2.50000000000000000000e-01,  0.00000000000000000000e+00 },
    {  2.00000000000000011102e-01, -1.11022302462515660205e-17 },
    { -1.66666666666666657415e-01, -9.25185853854297065662e-18 },
    {  1.42857142857142849213e-01,  7.93016446160826056281e-18 },
    { -1.25000000000000000000e-01,  0.00000000000000000000e+00 },
    {  1.11111111111111104943e-01,  6.16790569236198043774e-18 } };
  const double y = r.hi;
  gsl_dd q, p, r2;
  double s;
  int n;

  /* sum_{n=10}^{18} (-1)^(n+1) y^(n-10) / n */
  q.hi = -1.0 / 18.0;
  q.hi = 1.0 / 17.0 + y * q.hi;
  q.hi = -1.0 / 16.0 + y * q.hi;
  q.hi = 1.0 / 15.0 + y * q.hi;
  q.hi = -1.0 / 14.0 + y * q.hi;
  q.hi = 1.0 / 13.0 + y * q.hi;
  q.hi = -1.0 / 12.0 + y * q.hi;
  q.hi = 1.0 / 11.0 + y * q.hi;
  q.hi = -1.0 / 10.0 + y * q.hi;
  q.lo = 0.0;

  for (n = 6; n >= 0; --n)
    {
      p = DD_WRAP(dd_mul)(r, q);
      s = coef[n].hi + p.hi;
      q = fast_two_sum(s, ((coef[n].hi - s) + p.hi) + (coef[n].lo + p.lo));
    }

  /* -r^2/2, exactly */
  r2 = DD_WRAP(two_prod)(r.hi, r.hi);
  ex4_add_d(sum, -0.5 * r2.hi);
  ex4_add_d(sum, -0.5 * r2.lo);
  p = DD_WRAP(two_prod)(r.hi, r.lo);
  ex4_add_d(sum, -p.hi);
  ex4_add_d(sum, -p.lo - 0.5 * r.lo * r.lo);

  /* r^3 Q(r) */
  r2 = fast_two_sum(r2.hi, r2.lo + 2.0 * p.hi);
  p = DD_WRAP(dd_mul)(DD_WRAP(dd_mul)(r2, r), q);
  ex4_add_d(sum, p.hi);
  ex4_add_d(sum, p.lo);
}

/*
 * Accumulate e*log(2) + log(j/128) into a four-component expansion,
 * omitting the terms which vanish (e = 0 or j = 128).
 */
static void
DD_WRAP(log_reconstruct)(dd_ex4 * sum, const int e, const int j)
{
  ex4_zero(sum);

  if (e != 0)
    {
      gsl_dd p;

      p = DD_WRAP(two_prod)((double) e, DD_LN2_HI);
      ex4_add_d(sum, p.hi);
      ex4_add_d(sum, p.lo);

      p = DD_WRAP(two_prod)((double) e, DD_LN2_LO);
      ex4_add_d(sum, p.hi);
      ex4_add_d(sum, p.lo);

      ex4_add_d(sum, (double) e * DD_LN2_XTRA);
    }

  if (j != 128)
    {
      ex4_add_d(sum, log_table[j - 64][0]);
      ex4_add_d(sum, log_table[j - 64][1]);
      ex4_add_d(sum, log_table[j - 64][2]);
    }
}

/**********************************************
 * dd_log routines                            *
 **********************************************/

/*
 * Compute log(x) for a double-double argument x.
 *
 * Normalize
 *
 *   x = 2^e m,    0.75 <= m < 1.5,
 *
 * select the nearest center c = j/128, and form
 *
 *   r = (m - c) / c,    |r| <= 1/192.
 *
 * The subtraction m - c is exact (Sterbenz).  Then
 *
 *   log(x) = e*log(2) + log(c) + r + w,    w = log(1+r) - r,
 *
 * accumulated in a four-component expansion and rounded to
 * double-double once.  For x near 1, e = 0 and c = 1, so the result is
 * r + w with r = x - 1 exact and the relative accuracy is retained.
 *
 * log_sum() accumulates the expansion for positive finite x; dd_log()
 * handles the special arguments and rounds.
 */
static void
DD_WRAP(log_sum)(dd_ex4 * sum, const gsl_dd x)
{
  gsl_dd m, r;
  int e, j;
  double c;

  m = log_normalize(x, &e);
  j = log_center_index(m.hi);
  c = (double) j / 128.0;

  r = dd_sub_d(m, c);
  if (j != 128)
    r = DD_WRAP(dd_div_d)(r, c);

  DD_WRAP(log_reconstruct)(sum, e, j);
  ex4_add_d(sum, r.hi);
  ex4_add_d(sum, r.lo);
  DD_WRAP(log1pmx_small_add)(sum, r);
}

static gsl_dd
DD_WRAP(dd_log)(const gsl_dd x)
{
  dd_ex4 sum;

  if (!(x.hi > 0.0 && x.hi < GSL_POSINF))
    return log_special(x.hi);

  DD_WRAP(log_sum)(&sum, x);

  return ex4_to_dd(sum);
}

/**********************************************
 * dd_log1p routines                          *
 **********************************************/

/*
 * Accumulate log(1+x), x > -1 finite, into a four-component expansion.
 *
 * The exponent e and center c = j/128 are chosen from u = 1+x as in
 * dd_log(); with p = 2^e c the reduced argument is formed from the
 * original x,
 *
 *   r = (x - (p - 1)) / p,
 *
 * where p - 1 is exact in double-double.  This preserves the
 * information in x that is lost in the rounded value of 1+x (x near
 * -1, or |x| tiny), and for small |x| gives p = 1 and r = x exactly.
 * For x > 1e300, where p could overflow, log(x) is used.
 */
static void
DD_WRAP(log1p_sum)(dd_ex4 * sum, const gsl_dd x)
{
  gsl_dd m, r;
  int e, j;
  double p;

  if (x.hi > 1.0e300)
    {
      /* log(1+x) = log(x) + log(1+1/x), the second term below 1e-300 */
      DD_WRAP(log_sum)(sum, x);
      return;
    }

  m = log_normalize(dd_add_d(x, 1.0), &e);
  j = log_center_index(m.hi);
  p = ldexp((double) j / 128.0, e);

  if (p == 1.0)
    r = x;
  else
    {
      r = dd_sub(x, two_sum(p, -1.0));
      r = DD_WRAP(dd_div_d)(r, p);
    }

  DD_WRAP(log_reconstruct)(sum, e, j);
  ex4_add_d(sum, r.hi);
  ex4_add_d(sum, r.lo);
  DD_WRAP(log1pmx_small_add)(sum, r);
}

/*
 * Compute log(1+x) for a double-double argument x.
 */
static gsl_dd
DD_WRAP(dd_log1p)(const gsl_dd x)
{
  dd_ex4 sum;

  if (x.hi == 0.0)
    return x;                   /* +-0 */

  if (!(x.hi > -1.0 && x.hi < GSL_POSINF))
    {
      if (x.hi == -1.0 && x.lo > 0.0)
        ;                       /* 1 + x is a tiny positive number */
      else
        return log_special(x.hi + 1.0);
    }

  DD_WRAP(log1p_sum)(&sum, x);

  return ex4_to_dd(sum);
}

/**********************************************
 * dd_log1pmx routines                        *
 **********************************************/

/*
 * Compute log(1+x) - x on the central interval |x| <= 1/2.
 *
 * Choose the center p = j/128 nearest to 1+x (64 <= j <= 192), so that
 * 1 + x = p (1+r) with
 *
 *   r = (x - (p-1)) / p,    |r| <= 1/128,
 *
 * and use the exact identity
 *
 *   log1pmx(x) = log1pmx(p-1) + log1pmx(r) - (p-1) r.
 *
 * Here log1pmx(p-1) = log(p) - (p-1) is formed from the triple-double
 * table without loss, and log1pmx(r) comes from log1pmx_small_add().
 * The derivative of the last two terms with respect to r is about
 * -(p-1) - r, which is not small relative to the result, so r is used
 * to triple precision, r + r3, and r3 contributes -(p-1+r) r3.  For |x| < 1/256, p = 1 and
 * the result is the kernel value with r = x.
 */
static gsl_dd
DD_WRAP(log1pmx_central)(const gsl_dd x)
{
  gsl_dd d, R, r, t1, t2;
  dd_ex4 sum;
  double p, pm1, q1, q2, r3;
  int j;

  j = log_center_index(1.0 + x.hi);
  if (j == 128)
    {
      ex4_zero(&sum);
      DD_WRAP(log1pmx_small_add)(&sum, x);
      return ex4_to_dd(sum);
    }

  p = (double) j / 128.0;
  pm1 = (double) (j - 128) / 128.0;

  d = dd_sub_d(x, pm1);

  /*
   * r = d/p to triple precision, q1 + q2 + r3.  The remainder of a
   * correctly rounded quotient is representable, so d.hi - q1 p and
   * R.hi - q2 p are exact.
   */
  q1 = d.hi / p;
  t1 = DD_WRAP(two_prod)(q1, p);
  R = two_sum((d.hi - t1.hi) - t1.lo, d.lo);
  q2 = R.hi / p;
  t1 = DD_WRAP(two_prod)(q2, p);
  r3 = (((R.hi - t1.hi) - t1.lo) + R.lo) / p;
  r = fast_two_sum(q1, q2);

  ex4_zero(&sum);
  ex4_add_d(&sum, log_table[j - 64][0]);
  ex4_add_d(&sum, -pm1);
  ex4_add_d(&sum, log_table[j - 64][1]);
  ex4_add_d(&sum, log_table[j - 64][2]);

  t1 = DD_WRAP(two_prod)(pm1, r.hi);
  ex4_add_d(&sum, -t1.hi);
  ex4_add_d(&sum, -t1.lo);
  t2 = DD_WRAP(two_prod)(pm1, r.lo);
  ex4_add_d(&sum, -t2.hi);
  ex4_add_d(&sum, -t2.lo - (pm1 + r.hi) * r3);

  DD_WRAP(log1pmx_small_add)(&sum, r);

  return ex4_to_dd(sum);
}

/*
 * Compute log(1+x) - x for a double-double argument x.
 *
 * On |x| <= 1/2 use log1pmx_central().  Outside this interval the
 * cancellation between log(1+x) and x is at most a factor of about 4,
 * so -x is accumulated directly into the uncompressed log1p expansion
 * and the sum is rounded to double-double once.
 */
static gsl_dd
DD_WRAP(dd_log1pmx)(const gsl_dd x)
{
  dd_ex4 sum;

  if (!(x.hi > -1.0 && x.hi < GSL_POSINF))
    {
      if (x.hi == -1.0 && x.lo > 0.0)
        ;
      else
        {
          gsl_dd r = log_special(x.hi + 1.0);
          if (x.hi == GSL_POSINF)
            r.hi = GSL_NEGINF;  /* log(1+x) - x -> -Inf */
          return r;
        }
    }

  if (x.hi >= -0.5 && x.hi <= 0.5)
    return DD_WRAP(log1pmx_central)(x);

  DD_WRAP(log1p_sum)(&sum, x);
  ex4_add_d(&sum, -x.hi);
  ex4_add_d(&sum, -x.lo);

  return ex4_to_dd(sum);
}
