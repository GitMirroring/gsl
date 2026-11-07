/* dd/internal_exp_source.c
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

#define DD_EXP_ARG_MIN_GUARD         (-746.0)
#define DD_EXP_ARG_MAX_GUARD         (710.0)
#define DD_EXP_MIN_SCALE             (-1075)
#define DD_EXP_MAX_SCALE             (1024)
#define DD_HALF_ULP_ONE              (0.5 * GSL_DBL_EPSILON)
#define DD_QUARTER_ULP_ONE           (0.25 * GSL_DBL_EPSILON)
#define DD_EXP_64_OVER_LN2           (92.332482616893656768)  /* 64/log(2) */
#define DD_EXPM1_TINY                (-76.0)   /* exp(x) < 2^-109 */

/*
 * Add
 *
 *   expm1(r) = r + r^2/2 + r^3/3! + ...
 *
 * to the expansion sum, for a small double-double argument
 * |r| <= log(2)/128 + tiny.
 *
 * Write expm1(r) = r + r^2/2 + r^3 Q(r) with
 *
 *   Q(r) = sum_{n=3}^{12} r^(n-3) / n!.
 *
 * The terms r and r^2/2 are added exactly (r^2 is split by TwoProd).
 * The term r^3 Q(r) has relative weight below r^2/6 < 2^-17 and is
 * evaluated in double-double, by Horner for n = 3,...,6 and in double
 * precision from r.hi for n >= 7.  Since |r c_(n+1)| < |c_n| in every
 * Horner step, the steps use an addition without the cancellation
 * handling of dd_add().  The truncation error is below 2^-111 relative
 * to expm1(r).
 */
static void
DD_WRAP(expm1_small_add)(dd_ex4 * sum, const gsl_dd r)
{
  static const gsl_dd coef[4] = {
    { 1.66666666666666657415e-01,  9.25185853854297065662e-18 },
    { 4.16666666666666643537e-02,  2.31296463463574266415e-18 },
    { 8.33333333333333321769e-03,  1.15648231731787138023e-19 },
    { 1.38888888888888894189e-03, -5.30054395437357705906e-20 } };
  const double y = r.hi;
  gsl_dd q, p, r2;
  double s;
  int n;

  /* sum_{n=7}^{12} y^(n-7) / n! */
  q.hi = 2.08767569878681001866e-09;
  q.hi = 2.50521083854417202239e-08 + y * q.hi;
  q.hi = 2.75573192239858882758e-07 + y * q.hi;
  q.hi = 2.75573192239858925110e-06 + y * q.hi;
  q.hi = 2.48015873015873015658e-05 + y * q.hi;
  q.hi = 1.98412698412698412526e-04 + y * q.hi;
  q.lo = 0.0;

  for (n = 3; n >= 0; --n)
    {
      p = DD_WRAP(dd_mul)(r, q);
      s = coef[n].hi + p.hi;
      q = fast_two_sum(s, ((coef[n].hi - s) + p.hi) + (coef[n].lo + p.lo));
    }

  ex4_add_d(sum, r.hi);
  ex4_add_d(sum, r.lo);

  /* r^2/2, exactly */
  r2 = DD_WRAP(two_prod)(r.hi, r.hi);
  ex4_add_d(sum, 0.5 * r2.hi);
  ex4_add_d(sum, 0.5 * r2.lo);
  p = DD_WRAP(two_prod)(r.hi, r.lo);
  ex4_add_d(sum, p.hi);
  ex4_add_d(sum, p.lo + 0.5 * r.lo * r.lo);

  /* r^3 Q(r) */
  r2 = fast_two_sum(r2.hi, r2.lo + 2.0 * p.hi);
  p = DD_WRAP(dd_mul)(DD_WRAP(dd_mul)(r2, r), q);
  ex4_add_d(sum, p.hi);
  ex4_add_d(sum, p.lo);
}

/*
 * Range reduction and kernel shared by dd_exp() and dd_expm1().
 *
 * Write
 *
 *   a = (64 m + j) log(2)/64 + r,    0 <= j < 64,  |r| <= log(2)/128,
 *
 * so that exp(a) = 2^m T_j (1 + E) with T_j = 2^(j/64) and
 * E = expm1(r).  The reduction r = (64 a - N log(2))/64, N = 64 m + j,
 * uses the three-limb log(2).  If triple is nonzero, r and E are kept
 * to triple precision (their double-double rounding errors r3 and E3
 * are recovered from the expansions): expm1() needs this because for
 * small m its result is of the order of r, so those errors would not be
 * damped; for exp() they are weighted by |E| < 2^-7 and are dropped.
 *
 * On return, *m holds the binary exponent and sum holds T_j + T_j E
 * (without the factor 2^m), with T_j from the triple-double table and
 * the product T_j E formed from TwoProd terms.  If N = 0, sum holds E
 * itself instead, so that expm1() can use it without cancellation, and
 * the return value is 1; otherwise it is 0.
 */
static int
DD_WRAP(exp_core)(dd_ex4 * sum, const gsl_dd a, double * m, const int triple)
{
  const double N = floor(a.hi * DD_EXP_64_OVER_LN2 + 0.5);
  gsl_dd r, E, p;
  double r3, E3;
  const double * T;
  int j;

  /*
   * r + r3 = (64 a - N log(2)) / 64, with the rounding error r3 of the
   * double-double r recovered from the expansion
   */
  ex4_zero(sum);
  ex4_add_d(sum, 64.0 * a.hi);
  p = DD_WRAP(two_prod)(N, DD_LN2_HI);
  ex4_add_d(sum, -p.hi);
  ex4_add_d(sum, 64.0 * a.lo);
  ex4_add_d(sum, -p.lo);
  p = DD_WRAP(two_prod)(N, DD_LN2_LO);
  ex4_add_d(sum, -p.hi);
  ex4_add_d(sum, -p.lo - N * DD_LN2_XTRA);
  r = ex4_to_dd(*sum);
  r3 = 0.0;
  if (triple)
    {
      ex4_add_d(sum, -r.hi);
      ex4_add_d(sum, -r.lo);
      r3 = ex4_to_dd(*sum).hi / 64.0;
    }
  r = mul_pwr2_extred(r, 1.0 / 64.0);

  *m = floor(N / 64.0);
  j = (int) (N - 64.0 * (*m));

  /* E = expm1(r + r3) = expm1(r) + exp(r) r3 */
  ex4_zero(sum);
  DD_WRAP(expm1_small_add)(sum, r);
  if (triple)
    ex4_add_d(sum, r3 + r.hi * r3);

  if (N == 0.0)
    return 1;

  /* E rounded to double-double, and its rounding error E3 */
  E = ex4_to_dd(*sum);
  E3 = 0.0;
  if (triple)
    {
      ex4_add_d(sum, -E.hi);
      ex4_add_d(sum, -E.lo);
      E3 = ex4_to_dd(*sum).hi;
    }
  T = exp_table[j];

  ex4_zero(sum);
  ex4_add_d(sum, T[0]);
  ex4_add_d(sum, T[1]);

  p = DD_WRAP(two_prod)(T[0], E.hi);
  ex4_add_d(sum, p.hi);
  ex4_add_d(sum, p.lo);
  p = DD_WRAP(two_prod)(T[0], E.lo);
  ex4_add_d(sum, p.hi);
  ex4_add_d(sum, p.lo);
  p = DD_WRAP(two_prod)(T[1], E.hi);
  ex4_add_d(sum, p.hi);
  ex4_add_d(sum, p.lo);
  ex4_add_d(sum, T[2] + (T[0] * E3 + T[1] * E.lo + T[2] * E.hi));

  return 0;
}

/*
 * Compute exp(a) in double-double arithmetic.
 *
 * Using exp_core(), exp(a) = 2^m s with s = T_j (1 + E) in
 * [2^(-1/128), 2^(1+1/128)], rounded to double-double once.  The
 * binary scale is restored exactly by dd_ldexp(), with special handling
 * near the binary64 underflow and overflow boundaries which avoids
 * premature overflow/underflow during the final scaling and returns the
 * appropriate double-double representation at the ends of the
 * representable range.
 */
static gsl_dd
DD_WRAP(dd_exp)(const gsl_dd a)
{
  double m;
  gsl_dd s;
  dd_ex4 sum;

  if (a.hi != a.hi)
    return a;

  if (a.hi < DD_EXP_ARG_MIN_GUARD)
    {
      gsl_dd z;
      z.hi = 0.0;
      z.lo = 0.0;
      return z;
    }

  if (a.hi > DD_EXP_ARG_MAX_GUARD)
    {
      gsl_dd z;
      z.hi = GSL_POSINF;
      z.lo = 0.0;
      return z;
    }

  if (DD_WRAP(exp_core)(&sum, a, &m, 0))
    ex4_add_d(&sum, 1.0);       /* sum holds E */
  s = ex4_to_dd(sum);

  if (m > DD_EXP_MAX_SCALE)
    {
      gsl_dd z;
      z.hi = GSL_POSINF;
      z.lo = 0.0;
      return z;
    }

  if (m < DD_EXP_MIN_SCALE)
    {
      gsl_dd z;
      z.hi = 0.0;
      z.lo = 0.0;
      return z;
    }

  if (m == DD_EXP_MIN_SCALE)
    {
      gsl_dd z;
      if (s.hi > 1.0 || (s.hi == 1.0 && s.lo > 0.0))
        z.hi = ldexp(1.0, -1074);
      else
        z.hi = 0.0;
      z.lo = 0.0;
      return z;
    }

  if (m == DD_EXP_MAX_SCALE)
    {
      if (s.hi < 1.0)
        return dd_ldexp(s, 1024);
      else if (s.hi == 1.0 && s.lo <= -DD_QUARTER_ULP_ONE)
        {
          gsl_dd z;
          double d;
          z.hi = GSL_DBL_MAX;
          d = (s.hi - 1.0) + DD_HALF_ULP_ONE + s.lo;
          z.lo = ldexp(d, 1024);
          return z;
        }
      else
        {
          gsl_dd z;
          z.hi = GSL_POSINF;
          z.lo = 0.0;
          return z;
        }
    }

  return dd_ldexp(s, (int) m);
}

/*
 * Compute expm1(x) = exp(x) - 1 in double-double arithmetic.
 *
 * With exp_core(), expm1(x) = 2^m (T_j (1 + E) - 2^-m), where -2^-m is
 * added exactly to the unrounded expansion before the single rounding
 * to double-double, so the cancellation for x near 0 costs no accuracy.
 * For |x| <= log(2)/128, m = j = 0 and the expansion is E itself.
 *
 * For x < -76, exp(x) < 2^-109 and the result is -1 + exp(x); for
 * x > 700 the subtraction of 1 is harmless and exp(x) - 1 is used, which
 * also handles the overflow boundary.
 */
static gsl_dd
DD_WRAP(dd_expm1)(const gsl_dd x)
{
  double m;
  gsl_dd s;
  dd_ex4 sum;

  if (x.hi == 0.0 || x.hi != x.hi)
    return x;                   /* +-0 or NaN */

  if (x.hi < DD_EXPM1_TINY)
    {
      s.hi = -1.0;
      s.lo = DD_WRAP(dd_exp)(x).hi;
      return s;
    }

  if (x.hi > 700.0)
    {
      s = DD_WRAP(dd_exp)(x);
      return (s.hi < GSL_POSINF) ? dd_sub_d(s, 1.0) : s;
    }

  if (DD_WRAP(exp_core)(&sum, x, &m, 1))
    return ex4_to_dd(sum);      /* sum holds E */

  ex4_add_d(&sum, -ldexp(1.0, (int) -m));
  s = ex4_to_dd(sum);

  return dd_ldexp(s, (int) m);
}
